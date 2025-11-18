//+------------------------------------------------------------------+
//| WaveSpecZZ 1.0.4 Core — GPU-only wave extraction                  |
//| Mantém apenas o pipeline Zero-pad → Resample → Remove-DC → FFT  |
//| → filtros espectrais → IFFT e plota a wave dominante + potência   |
//+------------------------------------------------------------------+
#property copyright "Gen2Alglib"
#property link      "https://github.com/gen2alglib"
#property version   "1.0"
#property strict

#property indicator_separate_window
#property indicator_buffers 2
#property indicator_plots   2

#property indicator_label1  "GPU Wave"
#property indicator_type1   DRAW_LINE
#property indicator_color1  clrDeepSkyBlue
#property indicator_style1  STYLE_SOLID
#property indicator_width1  2

#property indicator_label2  "Wave Power"
#property indicator_type2   DRAW_LINE
#property indicator_color2  clrOrange
#property indicator_style2  STYLE_DOT
#property indicator_width2  1

input int    InpFFTWindow             = 4096;
input int    InpPadLeft               = 0;
input int    InpPadRight              = 0;
input bool   InpEnableZeroPad         = false;
input bool   InpEnableResample        = false;
input double InpResampleFactor        = 1.0;
input double InpResampleCutoff        = 0.45;
input int    InpResampleMethod        = 0;
input bool   InpEnableDcRemoval       = true;
input int    InpDcMode                = 0;
input double InpDcAlpha               = 0.98;
input bool   InpEnableSpectralDenoise = true;
input int    InpDenoiseMethod         = 0;
input double InpDenoiseThreshold      = 0.10;
input double InpDenoiseBeta           = 0.75;
input int    InpDenoiseIterations     = 1;
input bool   InpEnableSpectralUpscale = false;
input double InpSpectralUpscaleFactor = 1.0;
input int    InpSpectralUpscaleMode   = 0;
input int    InpSpectralUpscaleNormalize = 1;
input bool   InpEnableSpectralMask    = false;
input double InpMaskLow              = 0.15;
input double InpMaskHigh             = 0.85;
input bool   InpMaskUseZigZag         = false;
input int    InpMaskZigZagWidth       = 2;
input double InpMaskZigZagBlend       = 0.65;
input bool   InpEnablePhaseUnwrap     = false;
input int    InpPhaseMethod           = 0;
input bool   InpEnableSpectralConvolution = false;
input bool   InpEnableSpectralCorrelation = false;
input int    InpConvolutionPeriod     = 32;
input double InpConvolutionBandwidth  = 0.04;
input double InpConvolutionGain       = 1.0;

#import "mt-bridge.dll"
  int  gpu_init(int device_index, int stream_count);
  void gpu_shutdown(void);
  int  gpu_fft_real_forward(const double &in[], int len, double &out[]);
  int  gpu_fft_real_inverse(const double &in_spec[], int len, double &out[]);
  int  gpu_zero_pad_time_series(const double &series_in[], int len, int pad_left, int pad_right, double &out[], int out_cap, int &out_len);
  int  gpu_resample_time_series(const double &series_in[], int len, double factor, double cutoff, int method, double &out[], int out_cap, int &out_len);
  int  gpu_remove_dc_time_series(const double &series_in[], int len, int mode, double alpha, double &out[], int out_len);
  int  gpu_spectral_denoise(const double &spectrum[], int spectrum_len, int method, double threshold, double beta, int iterations, double &out[], int out_len);
  int  gpu_spectral_upscale(const double &spectrum[], int spectrum_len, double factor, int mode, int normalize, double &out[], int out_cap, int &out_len);
  int  gpu_apply_mask(const double &spectrum[], int spectrum_len, const double &mask[], int mask_len, int mask_is_complex, int mode, double &out[], int out_len);
  int  gpu_spectral_phase_unwrap(const double &spectrum[], int spectrum_len, int method, double &out[], int out_len);
  int  gpu_spectral_convolution(const double &lhs[], int lhs_len, const double &rhs[], int rhs_len, int mode, double &out[], int out_len);
  int  gpu_spectral_correlation(const double &lhs[], int lhs_len, const double &rhs[], int rhs_len, double &out[], int out_len);
#import

#define ALGLIB_STATUS_OK 0
#define ALGLIB_STATUS_NOT_READY -5

double WaveBuffer[];
double PowerBuffer[];

double price_data[];
double g_preproc_pipeline[];
double g_preproc_scratch[];
double g_preproc_fft[];
double g_preproc_mask[];
double g_preproc_phase[];
double g_zig_mask_temp[];
int    g_zig_mask_indices[64];
int    g_zig_mask_count = 0;
bool   g_zig_mask_ready = false;

double g_last_power = 0.0;
int    g_pending_fft_job_id = 0;
int    g_pending_fft_len = 0;
bool   g_fft_job_ready = false;

void EnsureArraySize(double &buffer[], const int len)
  {
   if(len <= 0)
      return;
   if(ArraySize(buffer) < len)
      ArrayResize(buffer, len);
  }

void CopySeriesToPipeline(const double &src[], const int len)
  {
   EnsureArraySize(g_preproc_pipeline, len);
   ArrayCopy(g_preproc_pipeline, src, 0, 0, len);
  }

void CopyPipelineToSeries(double &dest[], const int len)
  {
   ArrayCopy(dest, g_preproc_pipeline, 0, 0, len);
  }

void LogPreprocessFailure(const string &op, const int status)
  {
   PrintFormat("[WaveSpec Core] %s falhou (%d)", op, status);
  }

bool ApplyZeroPadStage(int &seriesLen)
  {
   if(!InpEnableZeroPad)
      return true;
   const int left  = MathMax(0, InpPadLeft);
   const int right = MathMax(0, InpPadRight);
   if(left == 0 && right == 0)
      return true;
   EnsureArraySize(g_preproc_scratch, seriesLen + left + right + 8);
   int outLen = 0;
   const int status = gpu_zero_pad_time_series(g_preproc_pipeline,
                                              seriesLen,
                                              left,
                                              right,
                                              g_preproc_scratch,
                                              ArraySize(g_preproc_scratch),
                                              outLen);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_zero_pad_time_series", status);
      return false;
     }
   CopySeriesToPipeline(g_preproc_scratch, outLen);
   seriesLen = outLen;
   return true;
  }

bool ApplyResampleStage(int &seriesLen, const int targetLen = -1)
  {
   if(!InpEnableResample)
     {
      if(targetLen > 0 && seriesLen != targetLen)
        {
         const double factor = (double)targetLen / MathMax(1, seriesLen);
         const double cutoff = MathMax(0.0, MathMin(0.5, InpResampleCutoff));
         EnsureArraySize(g_preproc_scratch, targetLen + 8);
         int outLen = 0;
         const int status = gpu_resample_time_series(g_preproc_pipeline,
                                                     seriesLen,
                                                     factor,
                                                     cutoff,
                                                     InpResampleMethod,
                                                     g_preproc_scratch,
                                                     ArraySize(g_preproc_scratch),
                                                     outLen);
         if(status != ALGLIB_STATUS_OK || outLen <= 0)
           {
            LogPreprocessFailure("gpu_resample_time_series", status);
            return false;
           }
         CopySeriesToPipeline(g_preproc_scratch, MathMin(outLen, targetLen));
         seriesLen = MathMin(outLen, targetLen);
        }
      return true;
     }

   const double factor = MathMax(1e-6, InpResampleFactor);
   const double cutoff = MathMax(0.0, MathMin(0.5, InpResampleCutoff));
   const int expected = (int)MathMax(4.0, MathCeil((double)seriesLen * factor + 4.0));
   EnsureArraySize(g_preproc_scratch, expected + 16);
   int outLen = 0;
   const int status = gpu_resample_time_series(g_preproc_pipeline,
                                               seriesLen,
                                               factor,
                                               cutoff,
                                               InpResampleMethod,
                                               g_preproc_scratch,
                                               ArraySize(g_preproc_scratch),
                                               outLen);
   if(status != ALGLIB_STATUS_OK || outLen <= 0)
     {
      LogPreprocessFailure("gpu_resample_time_series", status);
      return false;
     }
   CopySeriesToPipeline(g_preproc_scratch, outLen);
   seriesLen = outLen;
   if(targetLen > 0 && seriesLen != targetLen)
      return ApplyResampleStage(seriesLen, targetLen);
   return true;
  }

bool ApplyDcStage(const int seriesLen)
  {
   if(!InpEnableDcRemoval)
      return true;
   EnsureArraySize(g_preproc_scratch, seriesLen + 4);
   const int status = gpu_remove_dc_time_series(g_preproc_pipeline,
                                                seriesLen,
                                                InpDcMode,
                                                InpDcAlpha,
                                                g_preproc_scratch,
                                                seriesLen + 4);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_remove_dc_time_series", status);
      return false;
     }
   CopySeriesToPipeline(g_preproc_scratch, seriesLen);
   return true;
  }

void BuildMaskArray(const int len)
  {
   EnsureArraySize(g_preproc_mask, len);
   double low = MathMax(0.0, MathMin(1.0, InpMaskLow));
   double high = MathMax(0.0, MathMin(1.0, InpMaskHigh));
   if(high < low)
      high = low;
   for(int i = 0; i < len; ++i)
     {
      const double ratio = (len <= 1 ? 0.0 : (double)i / (double)(len - 1));
      double value = 1.0;
      if(ratio < low)
         value = 0.0;
      else if(ratio > high)
         value = 0.0;
      g_preproc_mask[i] = value;
     }
   if(InpMaskUseZigZag && g_zig_mask_ready)
     {
      const double blend = MathMax(0.0, MathMin(1.0, InpMaskZigZagBlend));
      EnsureArraySize(g_zig_mask_temp, len);
      ArrayInitialize(g_zig_mask_temp, 0.0);
      const int width = MathMax(1, InpMaskZigZagWidth);
      for(int idx = 0; idx < g_zig_mask_count - 1; ++idx)
        {
         const int delta = g_zig_mask_indices[idx + 1] - g_zig_mask_indices[idx];
         if(delta <= 0)
            continue;
         const int center = g_zig_mask_indices[idx];
         for(int offset = -width; offset <= width; ++offset)
           {
            int bin = center + offset;
            if(bin < 0 || bin >= len)
               continue;
            g_zig_mask_temp[bin] = 1.0;
           }
        }
      for(int i = 0; i < len; ++i)
         g_preproc_mask[i] = (1.0 - blend) * g_preproc_mask[i] + blend * g_zig_mask_temp[i];
     }
  }

void BuildConvolutionKernel(const int len)
  {
   if(len <= 0)
      return;
   EnsureArraySize(g_preproc_mask, len);
   const double period = MathMax(4.0, (double)MathMax(1, InpConvolutionPeriod));
   const double bandwidth = MathMax(0.0001, MathMin(0.5, InpConvolutionBandwidth));
   const double target = 1.0 / period;
   const double gain = MathMax(0.0, InpConvolutionGain);
   const double sigma = 2.0 * bandwidth * bandwidth;
   for(int i = 0; i < len; ++i)
     {
      const double freq = (double)i / (double)len;
      const double delta = freq - target;
      g_preproc_mask[i] = gain * MathExp(-delta * delta / sigma);
     }
  }

bool ApplySpectralConvolutionStage(const int spectrumLen)
  {
   if(!InpEnableSpectralConvolution || spectrumLen <= 0)
      return true;
   BuildConvolutionKernel(spectrumLen);
   EnsureArraySize(g_preproc_scratch, spectrumLen + 16);
   int status = gpu_spectral_convolution(g_preproc_fft,
                                         spectrumLen,
                                         g_preproc_mask,
                                         spectrumLen,
                                         g_preproc_scratch,
                                         spectrumLen + 16);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_spectral_convolution", status);
      return false;
     }
   ArrayCopy(g_preproc_fft, g_preproc_scratch, 0, 0, spectrumLen);
   return true;
  }

bool ApplySpectralCorrelationStage(const int spectrumLen)
  {
   if(!InpEnableSpectralCorrelation || spectrumLen <= 0)
      return true;
   BuildConvolutionKernel(spectrumLen);
   EnsureArraySize(g_preproc_scratch, spectrumLen + 16);
   int status = gpu_spectral_correlation(g_preproc_fft,
                                         spectrumLen,
                                         g_preproc_mask,
                                         spectrumLen,
                                         g_preproc_scratch,
                                         spectrumLen + 16);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_spectral_correlation", status);
      return false;
     }
   ArrayCopy(g_preproc_fft, g_preproc_scratch, 0, 0, spectrumLen);
   return true;
  }

bool AnalyzeSpectrumPower(const int spectrumLen)
  {
   if(spectrumLen <= 0)
      return false;
   double maxMag = 0.0;
   for(int i = 1; i < spectrumLen; ++i)
     {
      double mag = MathAbs(g_preproc_fft[i]);
      if(mag > maxMag)
         maxMag = mag;
     }
   g_last_power = maxMag;
   return true;
  }

bool ApplySpectralStages(int &seriesLen)
  {
   EnsureArraySize(g_preproc_fft, seriesLen + 16);
   ArrayResize(g_preproc_fft, seriesLen);
   int status = gpu_fft_real_forward(g_preproc_pipeline, seriesLen, g_preproc_fft);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_fft_real_forward", status);
      return false;
     }
   int spectrumLen = seriesLen;
   if(InpEnableSpectralDenoise && (spectrumLen % 2 == 0))
     {
      EnsureArraySize(g_preproc_scratch, spectrumLen + 16);
      status = gpu_spectral_denoise(g_preproc_fft,
                                     spectrumLen,
                                     InpDenoiseMethod,
                                     InpDenoiseThreshold,
                                     InpDenoiseBeta,
                                     InpDenoiseIterations,
                                     g_preproc_scratch,
                                     spectrumLen + 16);
      if(status == ALGLIB_STATUS_OK)
         ArrayCopy(g_preproc_fft, g_preproc_scratch, 0, 0, spectrumLen);
      else
         LogPreprocessFailure("gpu_spectral_denoise", status);
     }
   if(InpEnableSpectralUpscale && InpSpectralUpscaleFactor > 1.0)
     {
      const int projected = (int)MathMax(4.0, MathCeil((double)spectrumLen * InpSpectralUpscaleFactor + 4.0));
      EnsureArraySize(g_preproc_scratch, projected + 32);
      int outLen = 0;
      status = gpu_spectral_upscale(g_preproc_fft,
                                    spectrumLen,
                                    InpSpectralUpscaleFactor,
                                    InpSpectralUpscaleMode,
                                    InpSpectralUpscaleNormalize,
                                    g_preproc_scratch,
                                    ArraySize(g_preproc_scratch),
                                    outLen);
      if(status != ALGLIB_STATUS_OK)
        {
         LogPreprocessFailure("gpu_spectral_upscale", status);
         return false;
        }
      spectrumLen = outLen;
      EnsureArraySize(g_preproc_fft, spectrumLen + 16);
      ArrayCopy(g_preproc_fft, g_preproc_scratch, 0, 0, spectrumLen);
     }
   if(InpEnableSpectralMask)
     {
      BuildMaskArray(spectrumLen);
      EnsureArraySize(g_preproc_scratch, spectrumLen + 8);
      status = gpu_apply_mask(g_preproc_fft,
                               spectrumLen,
                               g_preproc_mask,
                               spectrumLen,
                               0,
                               0,
                               g_preproc_scratch,
                               spectrumLen + 8);
      if(status != ALGLIB_STATUS_OK)
        {
         LogPreprocessFailure("gpu_apply_mask", status);
         return false;
        }
      ArrayCopy(g_preproc_fft, g_preproc_scratch, 0, 0, spectrumLen);
     }
   if(!ApplySpectralConvolutionStage(spectrumLen))
      return false;
   if(!ApplySpectralCorrelationStage(spectrumLen))
      return false;
   if(InpEnablePhaseUnwrap)
     {
      EnsureArraySize(g_preproc_phase, spectrumLen + 8);
      status = gpu_spectral_phase_unwrap(g_preproc_fft,
                                         spectrumLen,
                                         InpPhaseMethod,
                                         g_preproc_phase,
                                         spectrumLen + 8);
      if(status != ALGLIB_STATUS_OK)
         LogPreprocessFailure("gpu_spectral_phase_unwrap", status);
     }
   AnalyzeSpectrumPower(spectrumLen);
   EnsureArraySize(g_preproc_scratch, spectrumLen + 16);
   status = gpu_fft_real_inverse(g_preproc_fft, spectrumLen, g_preproc_scratch);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_fft_real_inverse", status);
      return false;
     }
   CopySeriesToPipeline(g_preproc_scratch, spectrumLen);
   seriesLen = spectrumLen;
   return true;
  }

bool SubmitGpuFftJob(const int len)
  {
   if(len <= 0)
      return false;
   int jobId = 0;
   const int status = gpu_submit_fft_real_forward(g_preproc_pipeline, len, jobId);
   if(status != ALGLIB_STATUS_OK)
     {
      LogPreprocessFailure("gpu_submit_fft_real_forward", status);
      g_pending_fft_job_id = 0;
      g_pending_fft_len = 0;
      return false;
     }
   g_pending_fft_job_id = jobId;
   g_pending_fft_len = len;
   return true;
  }

bool CollectGpuFftJob()
  {
   if(g_pending_fft_job_id <= 0 || g_pending_fft_len <= 0)
      return false;
   int ready = 0;
   const int status = gpu_try_get_result(g_pending_fft_job_id,
                                          g_preproc_fft,
                                          g_pending_fft_len,
                                          ready);
   if(status == ALGLIB_STATUS_NOT_READY)
      return false;
   if(status == ALGLIB_STATUS_OK && ready != 0)
     {
      const int free_status = gpu_free_job(g_pending_fft_job_id);
      if(free_status != ALGLIB_STATUS_OK)
         LogPreprocessFailure("gpu_free_job", free_status);
      g_fft_job_ready = true;
      g_pending_fft_job_id = 0;
      return true;
     }
   if(status != ALGLIB_STATUS_OK)
      LogPreprocessFailure("gpu_try_get_result", status);
   if(g_pending_fft_job_id > 0)
      gpu_free_job(g_pending_fft_job_id);
   g_pending_fft_job_id = 0;
   g_pending_fft_len = 0;
   g_fft_job_ready = false;
   return false;
  }

bool EnsureFftResult(const int len)
  {
   if(len <= 0)
      return false;
   if(g_fft_job_ready && g_pending_fft_len == len)
      return true;
   if(g_pending_fft_job_id == 0)
     {
      if(!SubmitGpuFftJob(len))
         return false;
     }
   return CollectGpuFftJob();
  }

void ResetFftJobState()
  {
   g_fft_job_ready = false;
   g_pending_fft_len = 0;
  }

bool RunGpuPreprocessing(double &series[], const int targetLen)
  {
   CopySeriesToPipeline(series, targetLen);
   int seriesLen = targetLen;
   if(!ApplyZeroPadStage(seriesLen))
      return false;
   if(!ApplyResampleStage(seriesLen))
      return false;
   if(!ApplyDcStage(seriesLen))
      return false;
   if(!ApplySpectralStages(seriesLen))
      return false;
   if(seriesLen != targetLen)
     {
      if(!ApplyResampleStage(seriesLen, targetLen))
         return false;
      if(seriesLen != targetLen)
         ArrayResize(g_preproc_pipeline, targetLen);
     }
   CopyPipelineToSeries(series, targetLen);
   return true;
  }

bool EnsureWaveformGpuConfigured(const int length)
  {
   static bool initialized = false;
   if(!initialized)
     {
      const int status = gpu_init(0, 2);
      if(status != ALGLIB_STATUS_OK)
        {
         PrintFormat("[WaveSpec Core] gpu_init failed: %d", status);
         return false;
        }
      initialized = true;
     }
   return true;
  }

int OnInit()
  {
   IndicatorSetInteger(INDICATOR_DIGITS, _Digits);
   return(INIT_SUCCEEDED);
  }

int OnCalculate(const int rates_total,
                const int prev_calculated,
                const int begining,
                const double &open[],
                const double &high[],
                const double &low[],
                const double &close[],
                const long &tick_volume[],
                const long &volume[],
                const int &spread[])
  {
   if(rates_total < InpFFTWindow)
      return 0;
   if(!EnsureWaveformGpuConfigured(InpFFTWindow))
      return 0;
   const int start = MathMax(InpFFTWindow, prev_calculated - 1);
   for(int bar = start; bar < rates_total; ++bar)
     {
      const int offset = bar - InpFFTWindow + 1;
      ArrayCopy(price_data, close, 0, offset, InpFFTWindow);
      ResetFftJobState();
      if(!RunGpuPreprocessing(price_data, InpFFTWindow))
         continue;
      if(!EnsureFftResult(InpFFTWindow))
         continue;
      WaveBuffer[bar] = g_preproc_pipeline[InpFFTWindow - 1];
      PowerBuffer[bar] = g_last_power;
     }
   return rates_total;
  }

void OnDeinit(const int reason)
  {
   if(g_pending_fft_job_id > 0)
     {
      gpu_free_job(g_pending_fft_job_id);
      g_pending_fft_job_id = 0;
     }
   gpu_shutdown();
  }
