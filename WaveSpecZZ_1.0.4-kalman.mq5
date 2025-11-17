//+------------------------------------------------------------------+
//|                            WaveSpecZZ 1.0.4 - Kalman Only       |
//|   Reconstrução espectral usando GPU FFT + filtro de Kalman      |
//|   Sem o pipeline de 12 waves legacy / FollowFirst / ETA.        |
//+------------------------------------------------------------------+
#property copyright "Gen2Alglib"
#property link      ""
#property version   "1.0"
#property indicator_separate_window
#property indicator_buffers 1
#property indicator_plots   1

#property indicator_label1  "WaveKalman"
#property indicator_type1   DRAW_LINE
#property indicator_color1  clrWhite
#property indicator_style1  STYLE_SOLID
#property indicator_width1  2

#import "mt-bridge.dll"
  int  gpu_init(int device_index, int stream_count);
  void gpu_shutdown(void);
  int  gpu_fft_real_forward(const double &in[], int len, double &out[]);
#import

//--- Inputs
input ENUM_APPLIED_PRICE InpAppliedPrice = PRICE_CLOSE;
input int    InpFFTWindow      = 64;
input int    InpMinPeriod      = 12;
input int    InpMaxPeriod      = 256;
input int    InpTopCycles      = 8;
input bool   InpApplyHann      = true;

input double InpKalmanProcessNoise     = 0.25;  // Q
input double InpKalmanMeasurementNoise = 9.0;   // R
input double InpKalmanInitVariance     = 25.0;  // P0

input bool   InpLogDebug = false;

//--- Buffers
double WaveKalman[];

//--- GPU/FFT helpers
double g_price_window[];
double g_fft_interleaved[];
double g_fft_real[];
double g_fft_imag[];
double g_window[];
bool   g_gpu_ready = false;

//--- Kalman state (max 32 cycles) double-check
#define MAX_CYCLES 32
double g_kalman_weights[MAX_CYCLES];
double g_kalman_cov[MAX_CYCLES];
bool   g_kalman_ready = false;

//--- simple struct for ranking
struct CycleInfo { int index; double power; };

//+------------------------------------------------------------------+
int OnInit()
{
   if(InpFFTWindow < 64 || (InpFFTWindow & (InpFFTWindow - 1)) != 0)
   {
      Print("[KalmanWave] InpFFTWindow precisa ser potência de 2 >= 64");
      return INIT_FAILED;
   }
   ArraySetAsSeries(WaveKalman, true);
   SetIndexBuffer(0, WaveKalman, INDICATOR_DATA);
   PlotIndexSetDouble(0, PLOT_EMPTY_VALUE, EMPTY_VALUE);
   IndicatorSetString(INDICATOR_SHORTNAME, "WaveSpecZZ Kalman");
   ArrayResize(g_price_window, InpFFTWindow);
   ArrayResize(g_fft_interleaved, InpFFTWindow);
   ArrayResize(g_fft_real, InpFFTWindow);
   ArrayResize(g_fft_imag, InpFFTWindow);
   ArrayResize(g_window, InpFFTWindow);
   if(InpApplyHann)
   {
      for(int n=0;n<InpFFTWindow;n++)
         g_window[n] = 0.5*(1.0 - MathCos(2.0*M_PI*n/(InpFFTWindow-1)));
   }
   else
      ArrayFill(g_window, 0, InpFFTWindow, 1.0);
   ResetKalmanState();
   return INIT_SUCCEEDED;
}

void OnDeinit(const int)
{
   if(g_gpu_ready)
   {
      gpu_shutdown();
      g_gpu_ready = false;
   }
}

void ResetKalmanState()
{
   for(int i=0;i<MAX_CYCLES;i++)
   {
      g_kalman_weights[i]=0.0;
      g_kalman_cov[i]=MathMax(1e-6, InpKalmanInitVariance);
   }
   g_kalman_ready=false;
}

bool EnsureGpuReady()
{
   if(g_gpu_ready)
      return true;
   const int status = gpu_init(0, 2);
   if(status!=0)
   {
      PrintFormat("[KalmanWave] gpu_init falhou (%d)", status);
      return false;
   }
   g_gpu_ready = true;
   return true;
}

void FillPriceWindow(int start_pos,
                     const double &open[],
                     const double &high[],
                     const double &low[],
                     const double &close[])
{
   for(int j=0;j<InpFFTWindow;j++)
   {
      const int idx = start_pos + j;
      switch(InpAppliedPrice)
      {
         case PRICE_OPEN:     g_price_window[j] = open[idx];  break;
         case PRICE_HIGH:     g_price_window[j] = high[idx];  break;
         case PRICE_LOW:      g_price_window[j] = low[idx];   break;
         case PRICE_MEDIAN:   g_price_window[j] = 0.5*(high[idx]+low[idx]); break;
         case PRICE_TYPICAL:  g_price_window[j] = (high[idx]+low[idx]+close[idx])/3.0; break;
         case PRICE_WEIGHTED: g_price_window[j] = (high[idx]+low[idx]+2.0*close[idx])*0.25; break;
         default:             g_price_window[j] = close[idx]; break;
      }
      g_price_window[j] *= g_window[j];
   }
}

int CollectTopCycles(CycleInfo &cycles[])
{
   const int bins = InpFFTWindow/2;
   const int min_idx = (int)MathCeil((double)InpFFTWindow/InpMaxPeriod);
   const int max_idx = (int)MathFloor((double)InpFFTWindow/InpMinPeriod);
   int count=0;
   for(int k=MathMax(1,min_idx); k<=max_idx && k<bins; ++k)
   {
      const double re = g_fft_real[k];
      const double im = g_fft_imag[k];
      double power = re*re + im*im;
      ArrayResize(cycles, count+1);
      cycles[count].index = k;
      cycles[count].power = power;
      count++;
   }
   // ordenar por power desc
   for(int i=0; i<count-1; ++i)
   {
      int max_idx = i;
      double max_power = cycles[i].power;
      for(int j=i+1; j<count; ++j)
      {
         if(cycles[j].power > max_power)
         {
            max_idx = j;
            max_power = cycles[j].power;
         }
      }
      if(max_idx != i)
      {
         CycleInfo tmp = cycles[i];
         cycles[i] = cycles[max_idx];
         cycles[max_idx] = tmp;
      }
   }
   return count;
}

double ComputeContribution(const int k)
{
   const double re = g_fft_real[k];
   const double im = g_fft_imag[k];
   const double n0 = InpFFTWindow - 1;
   const double angle = 2.0*M_PI*k*n0/InpFFTWindow;
   const double cos_a = MathCos(angle);
   const double sin_a = MathSin(angle);
   const double scale = 2.0/InpFFTWindow;
   return scale*( re*cos_a - im*sin_a );
}

void UpdateKalman(const double &cycle_vals[], int cycle_count, double measurement, int bar_idx)
{
   const double Q = MathMax(1e-9, InpKalmanProcessNoise);
   const double R = MathMax(1e-9, InpKalmanMeasurementNoise);
   if(cycle_count>MAX_CYCLES) cycle_count = MAX_CYCLES;
   if(!g_kalman_ready)
      g_kalman_ready = true;

   double residual = measurement;
   double innovation = R;
   double cov_tmp[MAX_CYCLES];
   double weight_tmp[MAX_CYCLES];

   for(int i=0;i<cycle_count;i++)
   {
      g_kalman_cov[i] += Q;
      cov_tmp[i] = g_kalman_cov[i];
      weight_tmp[i] = g_kalman_weights[i];
      residual -= cycle_vals[i]*weight_tmp[i];
      innovation += cycle_vals[i]*cycle_vals[i]*cov_tmp[i];
   }
   if(innovation<1e-9)
      innovation = R;

   double blended = 0.0;
   for(int i=0;i<cycle_count;i++)
   {
      const double H = cycle_vals[i];
      const double cov = cov_tmp[i];
      const double K = (cov*H)/innovation;
      const double new_w = weight_tmp[i] + K*residual;
      const double new_cov = (1.0 - K*H)*cov;
      g_kalman_weights[i] = new_w;
      g_kalman_cov[i] = MathMax(new_cov, 1e-9);
      blended += new_w*H;
   }
   WaveKalman[bar_idx] = blended;
}

int OnCalculate(const int rates_total,
                const int prev_calculated,
                const datetime &time[],
                const double &open[],
                const double &high[],
                const double &low[],
                const double &close[],
                const long &tick_volume[],
                const long &volume[],
                const int &spread[])
{
   if(!EnsureGpuReady())
      return prev_calculated;
   const int window = InpFFTWindow;
   if(rates_total <= window)
   {
      for(int i=0;i<rates_total;i++) WaveKalman[i]=EMPTY_VALUE;
      return rates_total;
   }

   int start = MathMax(window-1, prev_calculated-1);
   if(start<window-1) start = window-1;

   CycleInfo cycles[];
   for(int bar=start; bar<rates_total; ++bar)
   {
      const int start_pos = bar - window + 1;
      FillPriceWindow(start_pos, open, high, low, close);
      const int status = gpu_fft_real_forward(g_price_window, window, g_fft_interleaved);
      if(status!=0)
      {
         if(InpLogDebug) PrintFormat("[KalmanWave] gpu_fft_real_forward falhou (%d) barra=%d", status, bar);
         WaveKalman[bar] = EMPTY_VALUE;
         continue;
      }
      for(int k=0;k<window/2;k++)
      {
         const int base = 2*k;
         g_fft_real[k] = g_fft_interleaved[base];
         g_fft_imag[k] = (base+1<window)?g_fft_interleaved[base+1]:0.0;
      }
      const int cycle_total = CollectTopCycles(cycles);
      const int use_cycles = MathMin(MathMin(InpTopCycles, MAX_CYCLES), cycle_total);
      if(use_cycles<=0)
      {
         WaveKalman[bar] = 0.0;
         continue;
      }
      double cycle_vals[MAX_CYCLES];
      for(int c=0;c<use_cycles;c++)
         cycle_vals[c] = ComputeContribution(cycles[c].index);
      UpdateKalman(cycle_vals, use_cycles, close[bar], bar);
   }

   return rates_total;
}
//+------------------------------------------------------------------+
