//+------------------------------------------------------------------+
//|  wave4ea_v1.0.5                                                 |
//|  GPU wave/phase view driven entirely by mt-bridge core outputs  |
//+------------------------------------------------------------------+
#property copyright "Copyright 2025"
#property link      ""
#property version   "1.100"
#property indicator_separate_window
#property indicator_buffers 25   // 12 plots * (data+color) + 1 phase buffer
#property indicator_plots   13   // 12 waves + 1 phase view
#property strict

#include "WavePresetDsl.mqh"

//--- Constants ------------------------------------------------------
#define WAVE_SLOTS 12
#define GPU_COLOR_BULL 0
#define GPU_COLOR_BEAR 1
static const int kGpuCycleAttrStride = 12;

//--- View mode ------------------------------------------------------
enum WAVE_SPEC_VIEW_MODE
  {
   WAVE_SPEC_VIEW_WAVES = 0,
   WAVE_SPEC_VIEW_PHASE = 1
  };

//--- Segment mix mode ------------------------------------------------
//--- Inputs ---------------------------------------------------------
input ENUM_APPLIED_PRICE InpAppliedPrice     = PRICE_CLOSE;
input int                InpFFTWindow        = 32768;
input int                InpAsyncWorkers     = 2;          // Max preset jobs in flight
input bool               InpEnableSegmentedFft = true;
input int                InpSegmentLength      = 16384;
input int                InpSegmentOverlap     = 4096;
input bool               InpSegmentAutoTune    = true;
input double             InpSegmentOverlapPct  = 0.25;
input SEGMENT_MIX_MODE   InpSegmentMixMode     = SEGMENT_MIX_ENERGY;
input int                InpWaveSlotsToShow    = 12;
input WAVE_SPEC_VIEW_MODE InpViewMode          = WAVE_SPEC_VIEW_WAVES;
input string             InpPresetTemplate     = "";       // Optional custom preset DSL
input string             InpPresetStageTime    = "";       // Optional override for time-domain stage
input string             InpPresetStageFreq    = "";       // Optional override for frequency-domain stage
input bool               InpLogDebug           = false;
input bool               InpUseTickSeries      = false;    // true => usa CopyTicksRange + gpu_wave_build_tick_series
input int                InpTickIntervalSeconds = 0;       // 0 = auto (usa PeriodSeconds)
input int                InpTickSmoothingWindow = 1;       // media móvel sobre série de ticks
input int                InpTickZigDepth        = 12;      // profundidade ZigZag
input double             InpTickZigDeviationPts = 20.0;    // desvio ZigZag em pontos
input int                InpTickZigBackstep     = 6;       // backstep ZigZag
input int                InpTickZigMode         = 0;       // 0 = estrito, 1 = seguro
input int                InpHudFontSize         = 12;      // HUD font size

//--- Indicator buffers ---------------------------------------------
double WaveBuffer1[],  WaveColor1[];
double WaveBuffer2[],  WaveColor2[];
double WaveBuffer3[],  WaveColor3[];
double WaveBuffer4[],  WaveColor4[];
double WaveBuffer5[],  WaveColor5[];
double WaveBuffer6[],  WaveColor6[];
double WaveBuffer7[],  WaveColor7[];
double WaveBuffer8[],  WaveColor8[];
double WaveBuffer9[],  WaveColor9[];
double WaveBuffer10[], WaveColor10[];
double WaveBuffer11[], WaveColor11[];
double WaveBuffer12[], WaveColor12[];
double PhaseBuffer[];

//--- GPU output scratch buffers ------------------------------------
double detrended_data[];
double price_window[];
double g_fft_interleaved[];
double fft_phase[];
double fft_unwrapped_phase[];
double fft_group_delay[];
double g_gpu_cycle_attrs[];
double g_wave_slot_values[WAVE_SLOTS];
double g_wave_slot_periods[WAVE_SLOTS];
double g_wave_slot_colors[WAVE_SLOTS];
MqlTick g_tick_cache[];
double  g_tick_prices[];
long    g_tick_times[];

//--- GPU session state ---------------------------------------------
bool   g_gpu_initialized = false;
bool   g_tick_warning_logged = false;
bool   g_tick_mode_active = false;
int    g_fft_window_len = 0;
int    g_last_completed_bar = -1;
string g_status_hud_object = "wave4ea_status_hud";
string g_last_bridge_call = "idle";

//--- Async preset job ----------------------------------------------
struct WaveAsyncJob
  {
   int    job_id;
   int    bar_index;
   double seconds_per_bar;
   int    fft_window_len;
   int    segment_len;
   int    overlap;
   int    mix_mode;
   bool   segmentation_enabled;
  };
WaveAsyncJob g_wave_jobs[];

struct CycleView
  {
   double amplitude;
   double phase;
   double period;
   double eta_seconds;
   double eta_bars;
   double score;
   double snr_db;
   double energy;
  };

CycleView g_cycle_views[];
int       g_cycle_order[];

enum GPU_ATTR_FIELD
  {
   GPU_ATTR_AMPLITUDE = 0,
   GPU_ATTR_FREQUENCY,
   GPU_ATTR_PERIOD,
   GPU_ATTR_PHASE,
   GPU_ATTR_ETA_BARS,
   GPU_ATTR_ETA_SECONDS,
   GPU_ATTR_ENERGY_RATIO,
   GPU_ATTR_COHERENCE,
   GPU_ATTR_SNR_DB,
   GPU_ATTR_RESIDUAL_POWER,
   GPU_ATTR_EIGEN_RATIO,
   GPU_ATTR_SCORE
  };

//--- Utility forward declarations ----------------------------------
void   WavePrint(const string text);
bool   EnsureWaveformGpuConfigured(int window_len);
void   ApplyViewMode();
void   EnsureArraySize(double &buffer[], int len);
int    ActiveFftWindowLen();
int    ResolveSegmentLength(int fft_window_len);
int    ResolveSegmentOverlap(int segment_len);
void   EnsureGpuCycleCapacity(int rows);
string ResolveWavePresetTemplate(int fft_window_len,
                                 int segment_len,
                                 int overlap,
                                 SEGMENT_MIX_MODE mix_mode,
                                 int top_cycles);
void   ApplyWavePresetOutputs(int bar_idx,
                              int fft_window_len,
                              int segment_len,
                              int overlap,
                              int mix_mode,
                              bool segmentation_enabled,
                              int cycle_count,
                              double kalman_value,
                              bool &logged_fft,
                              bool &logged_phase,
                              bool &logged_cycle,
                              bool &logged_kalman);
bool   PrepareWindowForBar(int bar_idx,
                           const datetime &time[],
                           const double &open[],
                           const double &high[],
                           const double &low[],
                           const double &close[],
                           double &seconds_per_bar);
bool   PrepareBarsPriceWindow(int bar_idx,
                              const double &open[],
                              const double &high[],
                              const double &low[],
                              const double &close[],
                              int window_len,
                              double &seconds_per_bar);
bool   PrepareTickPriceWindow(int bar_idx,
                              const datetime &time[],
                              int window_len,
                              double &seconds_per_sample);
void   ApplyWindowTransform(int window_len);
int    ResolveTickIntervalSeconds();
double ResolveTickPrice(const MqlTick &tick);
bool   SubmitWaveTemplateJob(int bar_idx,
                             double seconds_per_bar,
                             int fft_window_len,
                             int segment_len,
                             int overlap,
                             int mix_mode,
                             bool segmentation_enabled,
                             const string &preset_text);
bool   TryCompleteWaveTemplateJob(int index,
                                  bool &logged_fft,
                                  bool &logged_phase,
                                  bool &logged_cycle,
                                  bool &logged_kalman);
bool   HarvestWaveTemplateJobs(bool &logged_fft,
                               bool &logged_phase,
                               bool &logged_cycle,
                               bool &logged_kalman);
void   RemoveWaveJobAt(int index);
string FetchLastGpuError();
void   LogGpuFailure(const string &context, int status);
void   PublishStatusHud();
void   UpdateStatusHud(const string &text);
void   DestroyStatusHud();
void   TraceBridgeFlow(const string &text);
string DescribeMixModeLabel(int mix_mode);

//--- Bridge API -----------------------------------------------------
#import "mt-bridge.dll"
  int  mt_gpu_init(int device_index, int stream_count);
  void mt_gpu_shutdown(void);
  int  mt_gpu_wave_submit_template_job(string preset_text,
                                       const double &series[], int len, double seconds_per_sample,
                                       int fft_out_len, int phase_out_len, int unwrapped_len, int group_delay_len,
                                       int cycle_stride, int cycle_capacity, int wave_slots, int &job_id);
  int  mt_gpu_wave_try_get_template_job(int job_id,
                                        double &fft_out[], int fft_out_len,
                                        double &phase_out[], int phase_out_len,
                                        double &unwrapped_out[], int unwrapped_len,
                                        double &group_delay_out[], int group_delay_len,
                                        double &cycles_out[], int cycle_stride, int cycle_capacity, int &cycle_count,
                                        double &wave_values[], double &wave_periods[], double &wave_colors[], int wave_slots,
                                        double &kalman_value,
                                        int &ready);
  int  mt_gpu_wave_free_template_job(int job_id);
  int  mt_gpu_wave_build_tick_series(const double &tick_prices[], const long &tick_times[], int tick_count,
                                     int window_len, int interval_seconds, int smoothing_window,
                                     int zig_depth, double zig_deviation_points, int zig_backstep, int zig_mode,
                                     double point_value, double &out[], int out_len);
  int  mt_gpu_get_last_error_w(ushort &buf[], int buf_len);
#import

#import "tester.dll"
  int  tester_gpu_init(int device_index, int stream_count);
  void tester_gpu_shutdown(void);
  int  tester_gpu_wave_submit_template_job(string preset_text,
                                           const double &series[], int len, double seconds_per_sample,
                                           int fft_out_len, int phase_out_len, int unwrapped_len, int group_delay_len,
                                           int cycle_stride, int cycle_capacity, int wave_slots, int &job_id);
  int  tester_gpu_wave_try_get_template_job(int job_id,
                                            double &fft_out[], int fft_out_len,
                                            double &phase_out[], int phase_out_len,
                                            double &unwrapped_out[], int unwrapped_len,
                                            double &group_delay_out[], int group_delay_len,
                                            double &cycles_out[], int cycle_stride, int cycle_capacity, int &cycle_count,
                                            double &wave_values[], double &wave_periods[], double &wave_colors[], int wave_slots,
                                            double &kalman_value,
                                            int &ready);
  int  tester_gpu_wave_free_template_job(int job_id);
  int  tester_gpu_wave_build_tick_series(const double &tick_prices[], const long &tick_times[], int tick_count,
                                         int window_len, int interval_seconds, int smoothing_window,
                                         int zig_depth, double zig_deviation_points, int zig_backstep, int zig_mode,
                                         double point_value, double &out[], int out_len);
  int  tester_gpu_get_last_error_w(ushort &buf[], int buf_len);
#import

//--- ALGLIB status codes -------------------------------------------
#define ALGLIB_STATUS_OK               0
#define ALGLIB_STATUS_BAD_ARGS        -1
#define ALGLIB_STATUS_BACKEND_UNAVAILABLE -2
#define ALGLIB_STATUS_TIMEOUT         -3
#define ALGLIB_STATUS_INTERNAL_ERROR  -4
#define ALGLIB_STATUS_NOT_READY       -5
#define ALGLIB_STATUS_NO_MEM          -6

bool UseTesterBridge()
{
    static int flag = -1;
    if(flag == -1)
        flag = (int)MQLInfoInteger(MQL_TESTER);
    return flag == 1;
}

int gpu_init(int device_index, int stream_count)
{
    return UseTesterBridge() ? tester_gpu_init(device_index, stream_count)
                             : mt_gpu_init(device_index, stream_count);
}

void gpu_shutdown()
{
    if(UseTesterBridge())
        tester_gpu_shutdown();
    else
        mt_gpu_shutdown();
}

int gpu_wave_submit_template_job(const string &preset_text,
                                 const double &series[], int len, double seconds_per_sample,
                                 int fft_out_len, int phase_out_len, int unwrapped_len, int group_delay_len,
                                 int cycle_stride, int cycle_capacity, int wave_slots, int &job_id)
{
    return UseTesterBridge() ? tester_gpu_wave_submit_template_job(preset_text,
                                                                   series, len, seconds_per_sample,
                                                                   fft_out_len, phase_out_len, unwrapped_len, group_delay_len,
                                                                   cycle_stride, cycle_capacity, wave_slots, job_id)
                             : mt_gpu_wave_submit_template_job(preset_text,
                                                                series, len, seconds_per_sample,
                                                                fft_out_len, phase_out_len, unwrapped_len, group_delay_len,
                                                                cycle_stride, cycle_capacity, wave_slots, job_id);
}

int gpu_wave_try_get_template_job(int job_id,
                                  double &fft_out[], int fft_out_len,
                                  double &phase_out[], int phase_out_len,
                                  double &unwrapped_out[], int unwrapped_len,
                                  double &group_delay_out[], int group_delay_len,
                                  double &cycles_out[], int cycle_stride, int cycle_capacity, int &cycle_count,
                                  double &wave_values[], double &wave_periods[], double &wave_colors[], int wave_slots,
                                  double &kalman_value,
                                  int &ready)
{
    return UseTesterBridge() ? tester_gpu_wave_try_get_template_job(job_id,
                                                                    fft_out, fft_out_len,
                                                                    phase_out, phase_out_len,
                                                                    unwrapped_out, unwrapped_len,
                                                                    group_delay_out, group_delay_len,
                                                                    cycles_out, cycle_stride, cycle_capacity, cycle_count,
                                                                    wave_values, wave_periods, wave_colors, wave_slots,
                                                                    kalman_value,
                                                                    ready)
                             : mt_gpu_wave_try_get_template_job(job_id,
                                                                fft_out, fft_out_len,
                                                                phase_out, phase_out_len,
                                                                unwrapped_out, unwrapped_len,
                                                                group_delay_out, group_delay_len,
                                                                cycles_out, cycle_stride, cycle_capacity, cycle_count,
                                                                wave_values, wave_periods, wave_colors, wave_slots,
                                                                kalman_value,
                                                                ready);
}

int gpu_wave_free_template_job(int job_id)
{
    return UseTesterBridge() ? tester_gpu_wave_free_template_job(job_id)
                             : mt_gpu_wave_free_template_job(job_id);
}

int gpu_wave_build_tick_series(const double &tick_prices[],
                               const long &tick_times[],
                               int tick_count,
                               int window_len,
                               int interval_seconds,
                               int smoothing_window,
                               int zig_depth,
                               double zig_deviation_points,
                               int zig_backstep,
                               int zig_mode,
                               double point_value,
                               double &out[],
                               int out_len)
{
    return UseTesterBridge() ? tester_gpu_wave_build_tick_series(tick_prices,
                                                                 tick_times,
                                                                 tick_count,
                                                                 window_len,
                                                                 interval_seconds,
                                                                 smoothing_window,
                                                                 zig_depth,
                                                                 zig_deviation_points,
                                                                 zig_backstep,
                                                                 zig_mode,
                                                                 point_value,
                                                                 out,
                                                                 out_len)
                             : mt_gpu_wave_build_tick_series(tick_prices,
                                                             tick_times,
                                                             tick_count,
                                                             window_len,
                                                             interval_seconds,
                                                             smoothing_window,
                                                             zig_depth,
                                                             zig_deviation_points,
                                                             zig_backstep,
                                                             zig_mode,
                                                             point_value,
                                                             out,
                                                             out_len);
}

int gpu_get_last_error_w(ushort &buf[], int buf_len)
{
    return UseTesterBridge() ? tester_gpu_get_last_error_w(buf, buf_len)
                             : mt_gpu_get_last_error_w(buf, buf_len);
}

//--- Logging --------------------------------------------------------
void WavePrint(const string text)
{
    if(InpLogDebug)
        Print(text);
}

string FetchLastGpuError()
{
    ushort buf[];
    ArrayResize(buf, 512);
    if(gpu_get_last_error_w(buf, ArraySize(buf)) != ALGLIB_STATUS_OK)
        return "";
    return TrimPresetInput(ShortArrayToString(buf, 0, -1));
}

void LogGpuFailure(const string &context, int status)
{
    string detail = FetchLastGpuError();
    string line = StringFormat("%s (%d)", context, status);
    if(StringLen(detail) > 0)
        line = StringFormat("%s | %s", line, detail);
    WavePrint(line);
    if(!InpLogDebug)
        Print(line);
}

void TraceBridgeFlow(const string &text)
{
    g_last_bridge_call = (StringLen(text) > 64) ? StringSubstr(text, 0, 64) + "..." : text;
    const string ts = TimeToString(TimeLocal(), TIME_SECONDS);
    Print(StringFormat("[Wave4EA/bridge %s] %s", ts, text));
}

string DescribeMixModeLabel(int mix_mode)
{
    switch(mix_mode)
      {
       case SEGMENT_MIX_LATEST: return "latest";
       case SEGMENT_MIX_ENERGY: return "energy";
       default: return "average";
      }
}

void PublishStatusHud()
{
    string windowInfo = StringFormat("window %d", ActiveFftWindowLen());
    const int max_workers = MathMax(1, InpAsyncWorkers);
    string jobInfo = StringFormat("jobs %d/%d", ArraySize(g_wave_jobs), max_workers);
    string modeInfo = InpUseTickSeries ? (g_tick_mode_active ? "mode tick" : "mode tick(waiting)") : "mode bars";
    string gpuInfo = g_gpu_initialized ? "GPU ready" : "GPU init...";
    string lastInfo = (g_last_completed_bar >= 0) ? StringFormat("last %d", g_last_completed_bar) : "last --";
    string bridgeInfo = StringFormat("bridge %s", g_last_bridge_call);
    string status = StringFormat("Wave4EA | %s | %s | %s | %s | %s | %s",
                                 windowInfo,
                                 jobInfo,
                                 modeInfo,
                                 gpuInfo,
                                 lastInfo,
                                 bridgeInfo);
    UpdateStatusHud(status);
    Comment(status);
}

void UpdateStatusHud(const string &text)
{
    const long chart_id = ChartID();
    if(ObjectFind(chart_id, g_status_hud_object) < 0)
      {
       if(!ObjectCreate(chart_id, g_status_hud_object, OBJ_LABEL, 0, 0, 0))
           return;
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_CORNER, CORNER_LEFT_UPPER);
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_XDISTANCE, 8);
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_YDISTANCE, 14);
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_COLOR, clrWhite);
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_SELECTABLE, false);
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_BACK, false);
       ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_HIDDEN, true);
       ObjectSetString(chart_id, g_status_hud_object, OBJPROP_FONT, "Tahoma");
      }
    const int font_size = MathMax(8, InpHudFontSize);
    ObjectSetInteger(chart_id, g_status_hud_object, OBJPROP_FONTSIZE, font_size);
    ObjectSetString(chart_id, g_status_hud_object, OBJPROP_TEXT, text);
}

void DestroyStatusHud()
{
    const long chart_id = ChartID();
    if(ObjectFind(chart_id, g_status_hud_object) >= 0)
        ObjectDelete(chart_id, g_status_hud_object);
}

//--- Helpers --------------------------------------------------------
void EnsureArraySize(double &buffer[], int len)
{
    if(len <= 0)
        return;
    if(ArraySize(buffer) < len)
        ArrayResize(buffer, len);
}

int ActiveFftWindowLen()
{
    return MathMax(32, InpFFTWindow);
}

int ResolveSegmentLength(int fft_window_len)
{
    if(!InpEnableSegmentedFft)
        return fft_window_len;
    if(!InpSegmentAutoTune)
        return MathMin(fft_window_len, MathMax(256, InpSegmentLength));
    int tuned = (int)MathMax(1024, fft_window_len / 4);
    return MathMin(fft_window_len, tuned);
}

int ResolveSegmentOverlap(int segment_len)
{
    if(segment_len <= 0)
        return 0;
    if(!InpSegmentAutoTune)
        return MathMin(segment_len - 1, MathMax(0, InpSegmentOverlap));
    int overlap = (int)MathRound((double)segment_len * MathMax(0.0, MathMin(0.95, InpSegmentOverlapPct)));
    if(overlap >= segment_len)
        overlap = segment_len - 1;
    return MathMax(0, overlap);
}

void EnsureGpuCycleCapacity(int rows)
{
    const int required = MathMax(1, rows) * kGpuCycleAttrStride;
    if(ArraySize(g_gpu_cycle_attrs) < required)
        ArrayResize(g_gpu_cycle_attrs, required);
}

string TrimPresetInput(const string value)
{
    string trimmed = value;
    StringTrimLeft(trimmed);
    StringTrimRight(trimmed);
    return trimmed;
}

string ResolveWavePresetTemplate(int fft_window_len,
                                 int segment_len,
                                 int overlap,
                                 SEGMENT_MIX_MODE mix_mode,
                                 int top_cycles)
{
    const string custom = TrimPresetInput(InpPresetTemplate);
    if(StringLen(custom) > 0)
        return custom;

    const string stage_time = TrimPresetInput(InpPresetStageTime);
    const string stage_freq = TrimPresetInput(InpPresetStageFreq);

    return BuildWavePresetTemplate(segment_len,
                                   overlap,
                                   mix_mode,
                                   top_cycles,
                                   2.0,
                                   MathMax(2.0, (double)fft_window_len),
                                   WAVE_SLOTS,
                                   stage_time,
                                   stage_freq);
}

int LoadCycleViews(int cycle_count)
{
    const int capacity = ArraySize(g_gpu_cycle_attrs) / kGpuCycleAttrStride;
    const int available = MathMax(0, MathMin(cycle_count, capacity));
    ArrayResize(g_cycle_views, available);
    ArrayResize(g_cycle_order, available);
    for(int i = 0; i < available; ++i)
      {
       const int base = i * kGpuCycleAttrStride;
       g_cycle_views[i].amplitude   = g_gpu_cycle_attrs[base + GPU_ATTR_AMPLITUDE];
       g_cycle_views[i].phase       = g_gpu_cycle_attrs[base + GPU_ATTR_PHASE];
       g_cycle_views[i].period      = g_gpu_cycle_attrs[base + GPU_ATTR_PERIOD];
       g_cycle_views[i].eta_bars    = g_gpu_cycle_attrs[base + GPU_ATTR_ETA_BARS];
       g_cycle_views[i].eta_seconds = g_gpu_cycle_attrs[base + GPU_ATTR_ETA_SECONDS];
       g_cycle_views[i].score       = g_gpu_cycle_attrs[base + GPU_ATTR_SCORE];
       g_cycle_views[i].snr_db      = g_gpu_cycle_attrs[base + GPU_ATTR_SNR_DB];
       g_cycle_views[i].energy      = g_gpu_cycle_attrs[base + GPU_ATTR_ENERGY_RATIO];
       g_cycle_order[i] = i;
      }
    return available;
}

bool IsCycleBetter(const int lhs_index, const int rhs_index)
{
    if(lhs_index == rhs_index)
        return false;
    if(lhs_index < 0 || lhs_index >= ArraySize(g_cycle_views))
        return false;
    if(rhs_index < 0 || rhs_index >= ArraySize(g_cycle_views))
        return true;

    CycleView lhs = g_cycle_views[lhs_index];
    CycleView rhs = g_cycle_views[rhs_index];

    if(lhs.score != rhs.score)
        return lhs.score > rhs.score;
    if(lhs.eta_seconds != rhs.eta_seconds)
        return lhs.eta_seconds < rhs.eta_seconds;
    if(lhs.snr_db != rhs.snr_db)
        return lhs.snr_db > rhs.snr_db;
    return lhs.energy > rhs.energy;
}

void SortCycleViews()
{
    const int count = ArraySize(g_cycle_order);
    for(int i = 0; i < count - 1; ++i)
      {
       int best = i;
       for(int j = i + 1; j < count; ++j)
         {
          if(IsCycleBetter(g_cycle_order[j], g_cycle_order[best]))
             best = j;
         }
       if(best != i)
         {
          const int tmp = g_cycle_order[i];
          g_cycle_order[i] = g_cycle_order[best];
          g_cycle_order[best] = tmp;
         }
      }
}

void SetWaveValue(int slot, int bar_idx, double value)
{
    switch(slot)
      {
       case 0:  WaveBuffer1[bar_idx]  = value; break;
       case 1:  WaveBuffer2[bar_idx]  = value; break;
       case 2:  WaveBuffer3[bar_idx]  = value; break;
       case 3:  WaveBuffer4[bar_idx]  = value; break;
       case 4:  WaveBuffer5[bar_idx]  = value; break;
       case 5:  WaveBuffer6[bar_idx]  = value; break;
       case 6:  WaveBuffer7[bar_idx]  = value; break;
       case 7:  WaveBuffer8[bar_idx]  = value; break;
       case 8:  WaveBuffer9[bar_idx]  = value; break;
       case 9:  WaveBuffer10[bar_idx] = value; break;
       case 10: WaveBuffer11[bar_idx] = value; break;
       default: WaveBuffer12[bar_idx] = value; break;
      }
}

void SetWaveColor(int slot, int bar_idx, double color_index)
{
    switch(slot)
      {
       case 0:  WaveColor1[bar_idx]  = color_index; break;
       case 1:  WaveColor2[bar_idx]  = color_index; break;
       case 2:  WaveColor3[bar_idx]  = color_index; break;
       case 3:  WaveColor4[bar_idx]  = color_index; break;
       case 4:  WaveColor5[bar_idx]  = color_index; break;
       case 5:  WaveColor6[bar_idx]  = color_index; break;
       case 6:  WaveColor7[bar_idx]  = color_index; break;
       case 7:  WaveColor8[bar_idx]  = color_index; break;
       case 8:  WaveColor9[bar_idx]  = color_index; break;
       case 9:  WaveColor10[bar_idx] = color_index; break;
       case 10: WaveColor11[bar_idx] = color_index; break;
       default: WaveColor12[bar_idx] = color_index; break;
      }
}

void WriteWaveSlot(int slot, int bar_idx, double value, double color_flag)
{
    SetWaveValue(slot, bar_idx, value);
    SetWaveColor(slot, bar_idx, (color_flag > 0.5) ? GPU_COLOR_BULL : GPU_COLOR_BEAR);
}

void ApplyWavePresetOutputs(int bar_idx,
                            int fft_window_len,
                            int segment_len,
                            int overlap,
                            int mix_mode,
                            bool segmentation_enabled,
                            int cycle_count,
                            double kalman_value,
                            bool &logged_fft,
                            bool &logged_phase,
                            bool &logged_cycle,
                            bool &logged_kalman)
{
    const int cycle_views = LoadCycleViews(cycle_count);
    if(cycle_views > 1)
        SortCycleViews();

    const int slots_to_show = MathMax(1, MathMin(WAVE_SLOTS, InpWaveSlotsToShow));
    for(int slot = 0; slot < slots_to_show; ++slot)
      {
       double wave_value = (slot < WAVE_SLOTS) ? g_wave_slot_values[slot] : 0.0;
       if(!MathIsValidNumber(wave_value))
           wave_value = 0.0;
       double color_flag = (slot < WAVE_SLOTS) ? g_wave_slot_colors[slot] : 0.0;
       SetWaveValue(slot, bar_idx, wave_value);
       SetWaveColor(slot, bar_idx, color_flag);
      }
    for(int slot = slots_to_show; slot < WAVE_SLOTS; ++slot)
        {
         SetWaveValue(slot, bar_idx, 0.0);
         SetWaveColor(slot, bar_idx, 0.0);
        }

    if(cycle_views > 0)
        PhaseBuffer[bar_idx] = g_cycle_views[g_cycle_order[0]].phase;
    else
        PhaseBuffer[bar_idx] = 0.0;

    g_last_completed_bar = bar_idx;
    PublishStatusHud();

    if(!logged_fft && InpLogDebug)
      {
       WavePrint(StringFormat("[Wave4EA] GPU preset len=%d seg=%d ov=%d mix=%d",
                              fft_window_len,
                              segmentation_enabled ? segment_len : fft_window_len,
                              segmentation_enabled ? overlap : 0,
                              mix_mode));
       logged_fft = true;
      }
    if(!logged_phase && InpLogDebug)
      {
       WavePrint("[Wave4EA] Phase buffers updated (GPU preset)");
       logged_phase = true;
      }
    if(!logged_cycle && InpLogDebug)
      {
        const double best_score = (cycle_views > 0) ? g_cycle_views[g_cycle_order[0]].score : 0.0;
        const double best_snr   = (cycle_views > 0) ? g_cycle_views[g_cycle_order[0]].snr_db : 0.0;
       WavePrint(StringFormat("[Wave4EA] GPU cycles=%d (top score=%.3f snr=%.1f dB)", cycle_count, best_score, best_snr));
       logged_cycle = true;
      }
    if(!logged_kalman && InpLogDebug)
      {
        WavePrint(StringFormat("[Wave4EA] Kalman blend value=%.6f", kalman_value));
        logged_kalman = true;
      }
}

bool PrepareWindowForBar(int bar_idx,
                         const datetime &time[],
                         const double &open[],
                         const double &high[],
                         const double &low[],
                         const double &close[],
                         double &seconds_per_bar)
{
    const int window_len = ActiveFftWindowLen();
    if(window_len <= 0)
        return false;

    EnsureArraySize(price_window, window_len);
    EnsureArraySize(detrended_data, window_len);

    bool ready = false;
    if(InpUseTickSeries)
      {
       ready = PrepareTickPriceWindow(bar_idx, time, window_len, seconds_per_bar);
       if(ready)
         {
          if(!g_tick_mode_active)
              WavePrint("[Wave4EA] Tick resampling ativo (CopyTicksRange -> gpu_wave_build_tick_series)");
          g_tick_mode_active = true;
          g_tick_warning_logged = false;
          PublishStatusHud();
         }
       else
         {
          if(!g_tick_warning_logged)
              WavePrint("[Wave4EA] Tick resampling indisponivel, usando serie de barras");
          g_tick_warning_logged = true;
          g_tick_mode_active = false;
          PublishStatusHud();
         }
      }
    if(!ready)
        ready = PrepareBarsPriceWindow(bar_idx, open, high, low, close, window_len, seconds_per_bar);

    if(!ready)
        return false;

    ApplyWindowTransform(window_len);
    return true;
}

bool PrepareBarsPriceWindow(int bar_idx,
                            const double &open[],
                            const double &high[],
                            const double &low[],
                            const double &close[],
                            int window_len,
                            double &seconds_per_bar)
{
    const int start = bar_idx - window_len + 1;
    if(start < 0)
        return false;

    for(int i = 0; i < window_len; ++i)
      {
       const int idx = start + i;
       switch(InpAppliedPrice)
         {
          case PRICE_OPEN:    price_window[i] = open[idx]; break;
          case PRICE_HIGH:    price_window[i] = high[idx]; break;
          case PRICE_LOW:     price_window[i] = low[idx];  break;
          case PRICE_MEDIAN:  price_window[i] = (high[idx] + low[idx]) * 0.5; break;
          case PRICE_TYPICAL: price_window[i] = (high[idx] + low[idx] + close[idx]) / 3.0; break;
          case PRICE_WEIGHTED:
             price_window[i] = (high[idx] + low[idx] + close[idx] + close[idx]) * 0.25;
             break;
          default:
             price_window[i] = close[idx];
             break;
         }
      }

    seconds_per_bar = PeriodSeconds((ENUM_TIMEFRAMES)_Period);
    if(seconds_per_bar <= 0.0)
        seconds_per_bar = 60.0;
    return true;
}

int ResolveTickIntervalSeconds()
{
    if(InpTickIntervalSeconds > 0)
        return InpTickIntervalSeconds;
    int timeframe_seconds = (int)PeriodSeconds((ENUM_TIMEFRAMES)_Period);
    if(timeframe_seconds <= 0)
        timeframe_seconds = 60;
    return timeframe_seconds;
}

double ResolveTickPrice(const MqlTick &tick)
{
    if((tick.flags & TICK_FLAG_LAST) != 0 && tick.last > 0.0)
        return tick.last;
    if((tick.flags & TICK_FLAG_BID) != 0 && tick.bid > 0.0)
        return tick.bid;
    if((tick.flags & TICK_FLAG_ASK) != 0 && tick.ask > 0.0)
        return tick.ask;
    if(tick.last > 0.0)
        return tick.last;
    if(tick.bid > 0.0)
        return tick.bid;
    if(tick.ask > 0.0)
        return tick.ask;
    return 0.0;
}

bool PrepareTickPriceWindow(int bar_idx,
                            const datetime &time[],
                            int window_len,
                            double &seconds_per_sample)
{
    const int interval_seconds = ResolveTickIntervalSeconds();
    seconds_per_sample = (double)interval_seconds;

    const int timeframe_seconds = MathMax(1, (int)PeriodSeconds((ENUM_TIMEFRAMES)_Period));
    const datetime bar_time = time[bar_idx];
    const datetime bar_end = bar_time + timeframe_seconds;
    ulong to_msc = (bar_end > 0 ? (ulong)bar_end : 0) * 1000ULL;
    const ulong span = (ulong)interval_seconds * (ulong)MathMax(0, window_len - 1) * 1000ULL;
    const ulong from_msc = (to_msc > span) ? (to_msc - span) : 0ULL;

    ResetLastError();
    const int copied = CopyTicksRange(_Symbol, g_tick_cache, COPY_TICKS_ALL, from_msc, to_msc);
    if(copied <= 0)
      {
       const int err = GetLastError();
       PrintFormat("[Wave4EA] CopyTicksRange falhou (err=%d)", err);
       return false;
      }

    if(ArraySize(g_tick_prices) < copied)
        ArrayResize(g_tick_prices, copied);
    if(ArraySize(g_tick_times) < copied)
        ArrayResize(g_tick_times, copied);

    for(int i = 0; i < copied; ++i)
      {
       const double price = ResolveTickPrice(g_tick_cache[i]);
       g_tick_prices[i] = (price > 0.0) ? price : (i > 0 ? g_tick_prices[i - 1] : 0.0);
       g_tick_times[i] = (long)g_tick_cache[i].time_msc;
      }

    double point_value = SymbolInfoDouble(_Symbol, SYMBOL_POINT);
    if(point_value <= 0.0)
        point_value = 0.000001;

    const int status = gpu_wave_build_tick_series(g_tick_prices,
                                                  g_tick_times,
                                                  copied,
                                                  window_len,
                                                  interval_seconds,
                                                  InpTickSmoothingWindow,
                                                  InpTickZigDepth,
                                                  InpTickZigDeviationPts,
                                                  InpTickZigBackstep,
                                                  InpTickZigMode,
                                                  point_value,
                                                  price_window,
                                                  window_len);
    if(status != ALGLIB_STATUS_OK)
      {
       LogGpuFailure("[Wave4EA] gpu_wave_build_tick_series failed", status);
       return false;
      }

    return true;
}

void ApplyWindowTransform(int window_len)
{
    if(window_len <= 0)
        return;

    double mean = 0.0;
    for(int i = 0; i < window_len; ++i)
        mean += price_window[i];
    mean /= (double)window_len;

    if(window_len == 1)
      {
       detrended_data[0] = price_window[0] - mean;
       return;
      }

    const double denom = (double)(window_len - 1);
    for(int i = 0; i < window_len; ++i)
      {
       const double w = 0.5 - 0.5 * MathCos((2.0 * M_PI * i) / denom);
        detrended_data[i] = (price_window[i] - mean) * w;
      }
}

void RemoveWaveJobAt(int index)
{
    int last = ArraySize(g_wave_jobs) - 1;
    if(index < 0 || index > last)
        return;
    if(index != last)
        g_wave_jobs[index] = g_wave_jobs[last];
    ArrayResize(g_wave_jobs, last);
}

bool SubmitWaveTemplateJob(int bar_idx,
                           double seconds_per_bar,
                           int fft_window_len,
                           int segment_len,
                           int overlap,
                           int mix_mode,
                           bool segmentation_enabled,
                           const string &preset_text)
{
    int job_id = 0;
    const int status = gpu_wave_submit_template_job(preset_text,
                                                    detrended_data,
                                                    fft_window_len,
                                                    seconds_per_bar,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    kGpuCycleAttrStride,
                                                    ArraySize(g_gpu_cycle_attrs) / kGpuCycleAttrStride,
                                                    WAVE_SLOTS,
                                                    job_id);
    if(status != ALGLIB_STATUS_OK)
      {
       LogGpuFailure("[Wave4EA] gpu_wave_submit_template_job failed", status);
       return false;
      }

    const int next = ArraySize(g_wave_jobs);
    ArrayResize(g_wave_jobs, next + 1);
    g_wave_jobs[next].job_id = job_id;
    g_wave_jobs[next].bar_index = bar_idx;
    g_wave_jobs[next].seconds_per_bar = seconds_per_bar;
    g_wave_jobs[next].fft_window_len = fft_window_len;
    g_wave_jobs[next].segment_len = segment_len;
    g_wave_jobs[next].overlap = overlap;
    g_wave_jobs[next].mix_mode = mix_mode;
    g_wave_jobs[next].segmentation_enabled = segmentation_enabled;
    string preset_preview = preset_text;
    if(StringLen(preset_preview) > 96)
       preset_preview = StringSubstr(preset_preview, 0, 96) + "...";
    const string segment_info = segmentation_enabled ? StringFormat("%d/%d", segment_len, overlap) : "off";
    const string mode_info = InpUseTickSeries ? (g_tick_mode_active ? "tick" : "tick(wait)") : "bars";
    double sample_tail = (fft_window_len > 0) ? detrended_data[fft_window_len - 1] : 0.0;
    double sample_mid  = (fft_window_len > 1) ? detrended_data[fft_window_len / 2] : sample_tail;
    double sample_head = (fft_window_len > 0) ? detrended_data[0] : 0.0;
    TraceBridgeFlow(StringFormat(
        "submit -> gpu_wave_submit_template_job job=%d bar=%d window=%d seg=%s mix=%s seconds=%.2f mode=%s preset=\"%s\" samples=tail%.5f mid%.5f head%.5f",
        job_id,
        bar_idx,
        fft_window_len,
        segment_info,
        DescribeMixModeLabel(mix_mode),
        seconds_per_bar,
        mode_info,
        preset_preview,
        sample_tail,
        sample_mid,
        sample_head));
    return true;
}

bool TryCompleteWaveTemplateJob(int index,
                                bool &logged_fft,
                                bool &logged_phase,
                                bool &logged_cycle,
                                bool &logged_kalman)
{
    if(index < 0 || index >= ArraySize(g_wave_jobs))
        return false;

    WaveAsyncJob job = g_wave_jobs[index];
    int ready = 0;
    int cycle_count = 0;
    double kalman_value = 0.0;

    int status = gpu_wave_try_get_template_job(job.job_id,
                                                     g_fft_interleaved,
                                                     0,
                                                     fft_phase,
                                                     0,
                                                     fft_unwrapped_phase,
                                                     0,
                                                     fft_group_delay,
                                                     0,
                                                     g_gpu_cycle_attrs,
                                                     kGpuCycleAttrStride,
                                                     ArraySize(g_gpu_cycle_attrs) / kGpuCycleAttrStride,
                                                     cycle_count,
                                                     g_wave_slot_values,
                                                     g_wave_slot_periods,
                                                     g_wave_slot_colors,
                                                     WAVE_SLOTS,
                                                     kalman_value,
                                                     ready);

    if(status == ALGLIB_STATUS_NOT_READY || ready == 0)
      {
       TraceBridgeFlow(StringFormat(
           "poll <- gpu_wave_try_get_template_job job=%d status=%s(%d) ready=%d",
           job.job_id,
           StatusToString(status),
           status,
           ready));
        return false;
      }

    gpu_wave_free_template_job(job.job_id);
    RemoveWaveJobAt(index);

    if(status != ALGLIB_STATUS_OK)
      {
       LogGpuFailure("[Wave4EA] gpu_wave_try_get_template_job failed", status);
       TraceBridgeFlow(StringFormat(
           "result <- gpu_wave_try_get_template_job job=%d FAILED status=%s(%d)",
           job.job_id,
           StatusToString(status),
           status));
       return true;
      }

    ApplyWavePresetOutputs(job.bar_index,
                           job.fft_window_len,
                           job.segment_len,
                           job.overlap,
                           job.mix_mode,
                           job.segmentation_enabled,
                           cycle_count,
                           kalman_value,
                           logged_fft,
                           logged_phase,
                           logged_cycle,
                           logged_kalman);
    const double wave0 = (WAVE_SLOTS > 0) ? g_wave_slot_values[0] : 0.0;
    const double wave1 = (WAVE_SLOTS > 1) ? g_wave_slot_values[1] : 0.0;
    double cycle_amp = 0.0, cycle_period = 0.0;
    if(cycle_count > 0 && ArraySize(g_cycle_views) > 0 && ArraySize(g_cycle_order) > 0)
      {
       const int best_index = g_cycle_order[0];
        cycle_amp = g_cycle_views[best_index].amplitude;
        cycle_period = g_cycle_views[best_index].period;
      }
    TraceBridgeFlow(StringFormat(
        "result <- gpu_wave_try_get_template_job job=%d bar=%d status=%s cycles=%d kalman=%.6f waves=%.5f/%.5f top_cycle(A=%.5f,P=%.2f)",
        job.job_id,
        job.bar_index,
        StatusToString(status),
        cycle_count,
        kalman_value,
        wave0,
        wave1,
        cycle_amp,
        cycle_period));
    return true;
}

bool HarvestWaveTemplateJobs(bool &logged_fft,
                             bool &logged_phase,
                             bool &logged_cycle,
                             bool &logged_kalman)
{
    bool completed = false;
    for(int idx = ArraySize(g_wave_jobs) - 1; idx >= 0; --idx)
      {
       if(TryCompleteWaveTemplateJob(idx, logged_fft, logged_phase, logged_cycle, logged_kalman))
          completed = true;
      }
    return completed;
}

//--- Indicator lifecycle -------------------------------------------
int GetTopCyclesToShow()
{
    return MathMax(1, MathMin(WAVE_SLOTS, InpWaveSlotsToShow));
}

void ConfigureWavePlot(int plot, double &value_buffer[], double &color_buffer[])
{
    const int data_index  = plot * 2;
    const int color_index = data_index + 1;
    SetIndexBuffer(data_index,  value_buffer,  INDICATOR_DATA);
    SetIndexBuffer(color_index, color_buffer, INDICATOR_COLOR_INDEX);
    PlotIndexSetInteger(plot, PLOT_DRAW_TYPE, DRAW_COLOR_LINE);
    PlotIndexSetInteger(plot, PLOT_LINE_WIDTH, 1);
    PlotIndexSetInteger(plot, PLOT_COLOR_INDEXES, 2);
    PlotIndexSetInteger(plot, PLOT_LINE_COLOR, 0, clrDodgerBlue);
    PlotIndexSetInteger(plot, PLOT_LINE_COLOR, 1, clrTomato);
    PlotIndexSetString(plot, PLOT_LABEL, StringFormat("Wave %d", plot + 1));
}

int OnInit()
{
    g_status_hud_object = StringFormat("wave4ea_status_%I64d", ChartID());
    g_last_completed_bar = -1;
    DestroyStatusHud();
    if(!EnsureWaveformGpuConfigured(InpFFTWindow))
        return INIT_FAILED;

    ConfigureWavePlot(0, WaveBuffer1,  WaveColor1);
    ConfigureWavePlot(1, WaveBuffer2,  WaveColor2);
    ConfigureWavePlot(2, WaveBuffer3,  WaveColor3);
    ConfigureWavePlot(3, WaveBuffer4,  WaveColor4);
    ConfigureWavePlot(4, WaveBuffer5,  WaveColor5);
    ConfigureWavePlot(5, WaveBuffer6,  WaveColor6);
    ConfigureWavePlot(6, WaveBuffer7,  WaveColor7);
    ConfigureWavePlot(7, WaveBuffer8,  WaveColor8);
    ConfigureWavePlot(8, WaveBuffer9,  WaveColor9);
    ConfigureWavePlot(9, WaveBuffer10, WaveColor10);
    ConfigureWavePlot(10, WaveBuffer11, WaveColor11);
    ConfigureWavePlot(11, WaveBuffer12, WaveColor12);

    const int phase_plot = WAVE_SLOTS;
    SetIndexBuffer(phase_plot * 2, PhaseBuffer, INDICATOR_DATA);
    PlotIndexSetInteger(phase_plot, PLOT_DRAW_TYPE, DRAW_LINE);
    PlotIndexSetInteger(phase_plot, PLOT_LINE_WIDTH, 1);
    PlotIndexSetInteger(phase_plot, PLOT_LINE_COLOR, 0, clrLimeGreen);
    PlotIndexSetString(phase_plot, PLOT_LABEL, "Phase");

    ApplyViewMode();
    ArrayResize(g_wave_jobs, 0);
    g_fft_window_len = InpFFTWindow;
    PublishStatusHud();
    return INIT_SUCCEEDED;
}

void OnDeinit(const int reason)
{
    for(int idx = 0; idx < ArraySize(g_wave_jobs); ++idx)
      {
       if(g_wave_jobs[idx].job_id > 0)
          gpu_wave_free_template_job(g_wave_jobs[idx].job_id);
      }
    ArrayResize(g_wave_jobs, 0);

    if(g_gpu_initialized)
      {
       gpu_shutdown();
       g_gpu_initialized = false;
      }
    DestroyStatusHud();
}

bool EnsureWaveformGpuConfigured(int window_len)
{
    if(!g_gpu_initialized)
      {
       const int init_status = gpu_init(0, 2);
       if(init_status != ALGLIB_STATUS_OK)
         {
          LogGpuFailure("[Wave4EA] gpu_init failed", init_status);
          return false;
         }
       g_gpu_initialized = true;
       g_tick_warning_logged = false;
       g_tick_mode_active = false;
       WavePrint("[Wave4EA] GPU session ready");
      }
    g_fft_window_len = window_len;
    return true;
}

void ApplyViewMode()
{
    const bool show_phase = (InpViewMode == WAVE_SPEC_VIEW_PHASE);
    for(int plot = 0; plot < WAVE_SLOTS; ++plot)
        PlotIndexSetInteger(plot, PLOT_SHOW_DATA, show_phase ? false : true);
    PlotIndexSetInteger(WAVE_SLOTS, PLOT_SHOW_DATA, show_phase ? true : false);
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
    if(rates_total <= InpFFTWindow)
        return prev_calculated;

    static WAVE_SPEC_VIEW_MODE s_last_view = InpViewMode;
    if(s_last_view != InpViewMode)
      {
       ApplyViewMode();
       s_last_view = InpViewMode;
      }

    bool logged_fft            = false;
    bool logged_phase_analysis = false;
    bool logged_cycle_scan     = false;
    bool logged_kalman         = false;

    int start = (prev_calculated > 0) ? prev_calculated - 1 : InpFFTWindow;
    if(start < InpFFTWindow)
        start = InpFFTWindow;

    bool pending_exit = false;
    const bool segmentation_enabled = InpEnableSegmentedFft;
    const int max_workers = MathMax(1, InpAsyncWorkers);

    for(int bar = start; bar < rates_total && !IsStopped(); ++bar)
      {
       const bool completed = HarvestWaveTemplateJobs(logged_fft,
                                                      logged_phase_analysis,
                                                      logged_cycle_scan,
                                                      logged_kalman);
       if(ArraySize(g_wave_jobs) >= max_workers && !completed)
         {
          pending_exit = true;
          break;
         }
       if(ArraySize(g_wave_jobs) >= max_workers)
           continue;

       double seconds_per_bar = 0.0;
       if(!PrepareWindowForBar(bar, time, open, high, low, close, seconds_per_bar))
           continue;

       const int fft_window_len = ActiveFftWindowLen();
       int segment_len = segmentation_enabled ? ResolveSegmentLength(fft_window_len) : fft_window_len;
       if(segment_len <= 0)
           segment_len = fft_window_len;
       int overlap = segmentation_enabled ? ResolveSegmentOverlap(segment_len) : 0;
       const string preset_text = ResolveWavePresetTemplate(fft_window_len,
                                                            segment_len,
                                                            overlap,
                                                            InpSegmentMixMode,
                                                            GetTopCyclesToShow());

       EnsureGpuCycleCapacity(MathMax(WAVE_SLOTS, InpWaveSlotsToShow));
       if(!SubmitWaveTemplateJob(bar,
                                 seconds_per_bar,
                                 fft_window_len,
                                 segment_len,
                                 overlap,
                                 (int)InpSegmentMixMode,
                                 segmentation_enabled,
                                 preset_text))
         {
          pending_exit = true;
          break;
         }
      }

    const bool final_completed = HarvestWaveTemplateJobs(logged_fft,
                                                         logged_phase_analysis,
                                                         logged_cycle_scan,
                                                         logged_kalman);
    if(ArraySize(g_wave_jobs) > 0)
        pending_exit = true;
    else if(final_completed)
        pending_exit = false;

    PublishStatusHud();
    if(pending_exit)
        return prev_calculated;
    return rates_total;
}
