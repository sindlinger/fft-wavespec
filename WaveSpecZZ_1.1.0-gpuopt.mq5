// WaveSpecZZ minimal GPU+CPU FFT top-8 plotter (nodetrend, no ETA/leak/trackers)
#property copyright ""
#property link      ""
#property version   "1.101"
#property indicator_separate_window
#property indicator_buffers 66  // waves, feed, periods, forecast up/down, eta, phase, zigzag calc, MUSIC attrs
#property indicator_plots   25  // 8 waves + feed + 8 up + 8 down

// Buffers ZigZag mínimo (picos/fundos) para construir feed (calculations only)
double ZigzagPeakBuffer[], ZigzagBottomBuffer[], ColorBuffer[], HighMapBuffer[], LowMapBuffer[];

#include "Include\\imports.mqh"
#include "Include\\FeedCache.mqh"

#define ALGLIB_STATUS_OK 0
#define ALGLIB_STATUS_NOT_READY -5
int g_gpu_init_fail_counter = 0;
datetime g_gpu_init_last_log = 0;

// Inputs
input int  InpFFTWindow   = 4096;
input int  InpMinPeriod   = 18;
input int  InpMaxPeriod   = 200;

enum FEED_DATA_MODE { FEED_PLA = 0, FEED_ZIGZAG = 1, FEED_CLOSE = 2 };
input FEED_DATA_MODE InpFeedData = FEED_PLA;
// Feed sempre usa o timeframe do gráfico (Period()); input removido.

// Desenho das waves
enum DRAW_MODE { DRAW_POINTS = 0, DRAW_SINE_RECON = 1 };
input DRAW_MODE InpDrawMode = DRAW_SINE_RECON;

input group "PLA"
input int    InpPlaMaxSegments = 32;
input double InpPlaMaxError    = 0.0005;

input group "ZigZag"
input int InpZigZagDepth    = 12;
input int InpZigZagDeviation= 5;
input int InpZigZagBackstep = 3;
enum ZIG_MODE { ZIG_STEP = 0, ZIG_INTERP = 1, ZIG_MID = 2 };
input ZIG_MODE InpZigZagMode = ZIG_STEP;

input group "Feed View"
enum VIEW_MODE { VIEW_WAVES = 0, VIEW_FEED = 1 };
input VIEW_MODE InpViewMode = VIEW_WAVES; // Waves OU feed (exclusivo)

input group "GPU Cycle Extractor"
input int    InpGpuTopK        = 4;        // Número de ciclos extraídos da GPU (1–8)
input int    InpGpuMethod      = 1;        // 0=FFT ridge, 1=MUSIC/ESPRIT, -1=auto
input int    InpGpuArOrder     = 10;       // ordem AR para MUSIC/ESPRIT (ajustada p/ 2 ciclos “perfeitos”)
input double InpGpuMinPeriod   = 9;        // períodos em barras
input double InpGpuMaxPeriod   = 200;      // períodos em barras
input int    InpGpuStreams     = 64;       // streams GPU (1–512); mais streams = mais overlap GPU/CPU
input bool   InpEtaCountdown   = true;     // mostra contagem regressiva (barras) até próximo pico/fundo nos buffers de período
input bool   InpForecastMarks  = true;     // plota marcador na barra prevista do próximo pico/fundo (ETA) por wave
input int    InpBackfillWindows= 1;        // quantas janelas FFT recalcular no primeiro run (0 = todo histórico)
input int    InpMaxProcessBars = 20000;    // limite de barras a processar no primeiro run (0 = ilimitado)
input bool   InpEnableFeedCache = true;    // grava/usa cache de feed em arquivo (por símbolo/TF)
input bool   InpDebugCycles    = false;    // loga ciclos retornados pela GPU periodicamente
input bool   InpAsyncCycles    = true;     // modo assíncrono (submit/poll); false = síncrono
input int    InpAsyncDepth     = 64;       // máximo de jobs pendentes quando assíncrono (fila de jobs GPU)
input int    InpPrefetchBars   = 0;        // 0 = apenas o necessário; >0 força mínimo de cache
input bool   InpMusicOnly      = true;     // plotar apenas ciclos MUSIC/ESPRIT (descarta FFT ridge)
input bool   InpCycleCache     = true;     // cacheia waves/ciclos reconstruídos para warmup instantâneo
input bool   InpBatchWarmup    = true;     // usa api batch no primeiro run (hop=1) antes do modo online
input bool   InpForceBatch     = true;     // se true, ignora cache e refaz batch no attach
input int    InpBatchBarsLimit = 20000;    // máximo de barras usadas no warmup batch (0 = usa InpMaxLiveBars)
input int    InpBatchWaitMs    = 120000;   // tempo máximo para aguardar batch (ms); 0 = aguarda indefinidamente
input int    InpMaxLiveBars    = 120000;   // limite de barras processadas no gráfico (0=todas)

input group "MUSIC Weights"
input bool   InpUseMusicWeights = true;    // pondera amplitude por energia*coerência*SNR*score
input double InpMinCoherence    = 0.05;    // descarta (ou zera peso) abaixo deste valor
input double InpMinScore        = 0.01;    // idem para score
input double InpMinSnrDb        = -40.0;   // SNR mínimo considerado (dB)
input double InpMinEtaConf      = 0.0;     // confiança mínima para exibir forecast/ETA

input group "Kalman"
input bool   InpEnableKalman    = false;
input double InpKalmanFollowStrength  = 1.0;

// Buffers
#property indicator_label1  "Wave1"
#property indicator_type1   DRAW_LINE
#property indicator_color1  clrRed
#property indicator_width1  1
#property indicator_label2  "Wave2"
#property indicator_type2   DRAW_LINE
#property indicator_color2  clrOrangeRed
#property indicator_width2  1
#property indicator_label3  "Wave3"
#property indicator_type3   DRAW_LINE
#property indicator_color3  clrOrange
#property indicator_width3  1
#property indicator_label4  "Wave4"
#property indicator_type4   DRAW_LINE
#property indicator_color4  clrGold
#property indicator_width4  1
#property indicator_label5  "Wave5"
#property indicator_type5   DRAW_LINE
#property indicator_color5  clrYellow
#property indicator_width5  1
#property indicator_label6  "Wave6"
#property indicator_type6   DRAW_LINE
#property indicator_color6  clrChartreuse
#property indicator_width6  1
#property indicator_label7  "Wave7"
#property indicator_type7   DRAW_LINE
#property indicator_color7  clrLime
#property indicator_width7  1
#property indicator_label8  "Wave8"
#property indicator_type8   DRAW_LINE
#property indicator_color8  clrSpringGreen
#property indicator_width8  1

#property indicator_label9  "FeedTrace"
#property indicator_type9   DRAW_LINE
#property indicator_color9  clrWhite
#property indicator_style9  STYLE_DOT
#property indicator_width9  1

// Forecast markers UP (plots 10-17) - cores das waves
#property indicator_label10 "FUp1"
#property indicator_type10  DRAW_NONE
#property indicator_color10 clrRed
#property indicator_width10 1
#property indicator_label11 "FUp2"
#property indicator_type11  DRAW_NONE
#property indicator_color11 clrOrangeRed
#property indicator_width11 1
#property indicator_label12 "FUp3"
#property indicator_type12  DRAW_NONE
#property indicator_color12 clrOrange
#property indicator_width12 1
#property indicator_label13 "FUp4"
#property indicator_type13  DRAW_NONE
#property indicator_color13 clrGold
#property indicator_width13 1
#property indicator_label14 "FUp5"
#property indicator_type14  DRAW_NONE
#property indicator_color14 clrYellow
#property indicator_width14 1
#property indicator_label15 "FUp6"
#property indicator_type15  DRAW_NONE
#property indicator_color15 clrChartreuse
#property indicator_width15 1
#property indicator_label16 "FUp7"
#property indicator_type16  DRAW_NONE
#property indicator_color16 clrLime
#property indicator_width16 1
#property indicator_label17 "FUp8"
#property indicator_type17  DRAW_NONE
#property indicator_color17 clrSpringGreen
#property indicator_width17 1

// Forecast markers DOWN (plots 18-25) - cores das waves
#property indicator_label18 "FDn1"
#property indicator_type18  DRAW_NONE
#property indicator_color18 clrRed
#property indicator_width18 1
#property indicator_label19 "FDn2"
#property indicator_type19  DRAW_NONE
#property indicator_color19 clrOrangeRed
#property indicator_width19 1
#property indicator_label20 "FDn3"
#property indicator_type20  DRAW_NONE
#property indicator_color20 clrOrange
#property indicator_width20 1
#property indicator_label21 "FDn4"
#property indicator_type21  DRAW_NONE
#property indicator_color21 clrGold
#property indicator_width21 1
#property indicator_label22 "FDn5"
#property indicator_type22  DRAW_NONE
#property indicator_color22 clrYellow
#property indicator_width22 1
#property indicator_label23 "FDn6"
#property indicator_type23  DRAW_NONE
#property indicator_color23 clrChartreuse
#property indicator_width23 1
#property indicator_label24 "FDn7"
#property indicator_type24  DRAW_NONE
#property indicator_color24 clrLime
#property indicator_width24 1
#property indicator_label25 "FDn8"
#property indicator_type25  DRAW_NONE
#property indicator_color25 clrSpringGreen
#property indicator_width25 1


double WaveBuffer1[],WaveBuffer2[],WaveBuffer3[],WaveBuffer4[];
double WaveBuffer5[],WaveBuffer6[],WaveBuffer7[],WaveBuffer8[];
double WavePeriod1[],WavePeriod2[],WavePeriod3[],WavePeriod4[];
double WavePeriod5[],WavePeriod6[],WavePeriod7[],WavePeriod8[];
double WaveKalman[];
double FeedTrace[];

// ETA countdown por slot (barras até próximo pico/fundo) - corrente e por‑barra
double EtaCountdown[8];
double EtaCount1[],EtaCount2[],EtaCount3[],EtaCount4[],EtaCount5[],EtaCount6[],EtaCount7[],EtaCount8[];
// Fase (radianos) por wave
double PhaseVal1[],PhaseVal2[],PhaseVal3[],PhaseVal4[],PhaseVal5[],PhaseVal6[],PhaseVal7[],PhaseVal8[];

// Atributos avançados do MUSIC para consulta por EA / Data Window
double MusEnergy1[], MusEnergy2[];
double MusCoher1[], MusCoher2[];
double MusSnrDb1[], MusSnrDb2[];
double MusScore1[], MusScore2[];
double MusEigen1[], MusEigen2[];
double MusEtaConf1[], MusEtaConf2[];
// Marcadores de previsão (um buffer por wave, up e down)
double ForecastMark1[],ForecastMark2[],ForecastMark3[],ForecastMark4[];
double ForecastMark5[],ForecastMark6[],ForecastMark7[],ForecastMark8[];
double ForecastDn1[],ForecastDn2[],ForecastDn3[],ForecastDn4[];
double ForecastDn5[],ForecastDn6[],ForecastDn7[],ForecastDn8[];

void ClearForecastBuffersAt(const int i)
{
    ForecastMark1[i]=ForecastMark2[i]=ForecastMark3[i]=ForecastMark4[i]=ForecastMark5[i]=ForecastMark6[i]=ForecastMark7[i]=ForecastMark8[i]=EMPTY_VALUE;
    ForecastDn1[i]=ForecastDn2[i]=ForecastDn3[i]=ForecastDn4[i]=ForecastDn5[i]=ForecastDn6[i]=ForecastDn7[i]=ForecastDn8[i]=EMPTY_VALUE;
}

string CycleCacheName()
{
    // WaveSpecZZ_cycles_<sym>_<tf>_w<win>_m<method>_ar<ar>_k<topk>.bin
    string tf = EnumToString(FeedTF());
    return StringFormat("WaveSpecZZ_cycles_%s_%s_w%d_m%d_ar%d_k%d.bin", _Symbol, tf, InpFFTWindow, InpGpuMethod, InpGpuArOrder, InpGpuTopK);
}

bool LoadCycleCache(const int rates_total)
{
    if(!InpCycleCache) return false;
    g_cycle_cache_file = CycleCacheName();
    int h = FileOpen(g_cycle_cache_file, FILE_READ|FILE_BIN|FILE_COMMON);
    if(h==INVALID_HANDLE) return false;

    int version = FileReadInteger(h, INT_VALUE);
    if(version != 1) { FileClose(h); return false; }
    int bars = FileReadInteger(h, INT_VALUE);
    int topk = FileReadInteger(h, INT_VALUE);
    if(topk < 1 || topk > 2) { FileClose(h); return false; }
    int count = MathMin(bars, rates_total);

    // garante tamanho
    int need = count;
    ArrayResize(WaveBuffer1, need); ArrayResize(WaveBuffer2, need);
    ArrayResize(WavePeriod1, need); ArrayResize(WavePeriod2, need);
    ArrayResize(EtaCount1, need);   ArrayResize(EtaCount2, need);
    ArrayResize(PhaseVal1, need);   ArrayResize(PhaseVal2, need);
    ArrayResize(MusEnergy1,need);   ArrayResize(MusEnergy2,need);
    ArrayResize(MusCoher1, need);   ArrayResize(MusCoher2, need);
    ArrayResize(MusSnrDb1,need);    ArrayResize(MusSnrDb2,need);
    ArrayResize(MusScore1,need);    ArrayResize(MusScore2,need);
    ArrayResize(MusEigen1,need);    ArrayResize(MusEigen2,need);
    ArrayResize(MusEtaConf1,need);  ArrayResize(MusEtaConf2,need);

    for(int i=0;i<count;i++)
    {
        WaveBuffer1[i] = FileReadDouble(h);
        WaveBuffer2[i] = FileReadDouble(h);
        WavePeriod1[i] = FileReadDouble(h);
        WavePeriod2[i] = FileReadDouble(h);
        EtaCount1[i]   = FileReadDouble(h);
        EtaCount2[i]   = FileReadDouble(h);
        PhaseVal1[i]   = FileReadDouble(h);
        PhaseVal2[i]   = FileReadDouble(h);
        MusEnergy1[i]  = FileReadDouble(h);
        MusEnergy2[i]  = FileReadDouble(h);
        MusCoher1[i]   = FileReadDouble(h);
        MusCoher2[i]   = FileReadDouble(h);
        MusSnrDb1[i]   = FileReadDouble(h);
        MusSnrDb2[i]   = FileReadDouble(h);
        MusScore1[i]   = FileReadDouble(h);
        MusScore2[i]   = FileReadDouble(h);
        MusEigen1[i]   = FileReadDouble(h);
        MusEigen2[i]   = FileReadDouble(h);
        MusEtaConf1[i] = FileReadDouble(h);
        MusEtaConf2[i] = FileReadDouble(h);
    }
    FileClose(h);
    g_cycle_cache_loaded = true;
    PrintFormat("[WaveSpecZZ][CACHE] cycles loaded (%d bars) file=%s", count, g_cycle_cache_file);
    return true;
}

void SaveCycleCache(const int bars)
{
    if(!InpCycleCache) return;
    if(g_cycle_cache_file == "")
        g_cycle_cache_file = CycleCacheName();

    int h = FileOpen(g_cycle_cache_file, FILE_WRITE|FILE_BIN|FILE_COMMON);
    if(h==INVALID_HANDLE) return;
    FileWriteInteger(h, 1, INT_VALUE);     // version
    FileWriteInteger(h, bars, INT_VALUE);  // bars
    FileWriteInteger(h, 2, INT_VALUE);     // topK (fixo 2 waves)

    for(int i=0;i<bars;i++)
    {
        FileWriteDouble(h, WaveBuffer1[i]);
        FileWriteDouble(h, WaveBuffer2[i]);
        FileWriteDouble(h, WavePeriod1[i]);
        FileWriteDouble(h, WavePeriod2[i]);
        FileWriteDouble(h, EtaCount1[i]);
        FileWriteDouble(h, EtaCount2[i]);
        FileWriteDouble(h, PhaseVal1[i]);
        FileWriteDouble(h, PhaseVal2[i]);
        FileWriteDouble(h, MusEnergy1[i]);
        FileWriteDouble(h, MusEnergy2[i]);
        FileWriteDouble(h, MusCoher1[i]);
        FileWriteDouble(h, MusCoher2[i]);
        FileWriteDouble(h, MusSnrDb1[i]);
        FileWriteDouble(h, MusSnrDb2[i]);
        FileWriteDouble(h, MusScore1[i]);
        FileWriteDouble(h, MusScore2[i]);
        FileWriteDouble(h, MusEigen1[i]);
        FileWriteDouble(h, MusEigen2[i]);
        FileWriteDouble(h, MusEtaConf1[i]);
        FileWriteDouble(h, MusEtaConf2[i]);
    }
    FileClose(h);
    PrintFormat("[WaveSpecZZ][CACHE] cycles saved (%d bars) file=%s", bars, g_cycle_cache_file);
}

double feed_data[], detrended_data[];
double g_fft_interleaved[], fft_real[], fft_imag[], spectrum[];
// buffers para chamada extract_cycles GPU
double g_cycles_raw[]; // interleaved: amplitude, freq, period, phase, eta_bars, eta_seconds, energy_ratio, coherence, snr_db, residual_power, eigen_ratio, score, kalman_pred, eta_confidence, method

// Cache de ondas/ciclos já reconstruídos (warmup rápido)
bool   g_cycle_cache_loaded = false;
string g_cycle_cache_file = "";
// Feed timeframe helper: sempre usa o timeframe do gráfico
ENUM_TIMEFRAMES FeedTF(){ return (ENUM_TIMEFRAMES)Period(); }
bool g_gpu_session = false;
bool g_gpu_init_logged = false;
int g_zig_handle = INVALID_HANDLE;
const int kFeedLogEvery = 128; // log feed status a cada 128 barras

// Cache de feed (Close) para reduzir chamadas repetidas
static FeedCache g_feed_cache;
bool g_cache_use_logged = false;
long g_cycles_job_id = 0;
long g_cycles_jobs[];            // fila de jobs pendentes (assíncronos)
int  g_cycles_job_count = 0;
int  g_cycles_job_capacity = 0;
double g_cycles_job_buf[];
int g_cycles_job_out = 0;
bool g_cycles_sync_boot = true; // faz uma rodada síncrona inicial para popular buffers
bool g_mode_logged = false;     // log único de modo/params
// Progresso de backfill
int  g_prog_start = -1;
int  g_prog_total = 0;
bool g_prog_done  = false;
bool g_prog_logged_done = false;

//---------------- OO Helper structs ----------------
class ZigZagFeed
{
public:
    bool LoadWindow(const int shift_end_feed,int len,double &main_ch[],double &high_ch[],double &low_ch[])
    {
        if(shift_end_feed<0) return false;
        static double zz_main[], zz_high[], zz_low[], zz_peak[], zz_bottom[];
        ArraySetAsSeries(zz_main, true);  ArraySetAsSeries(zz_high, true);  ArraySetAsSeries(zz_low, true);
        ArraySetAsSeries(zz_peak, true);  ArraySetAsSeries(zz_bottom, true);
        ArrayResize(zz_main, len); ArrayResize(zz_high, len); ArrayResize(zz_low, len);
        ArrayResize(zz_peak, len); ArrayResize(zz_bottom, len);

        int copied_main   = CopyBuffer(g_zig_handle, 0, shift_end_feed, len, zz_main);
        int copied_high   = CopyBuffer(g_zig_handle, 1, shift_end_feed, len, zz_high);
        int copied_low    = CopyBuffer(g_zig_handle, 2, shift_end_feed, len, zz_low);
        int copied_peak   = CopyBuffer(g_zig_handle, 0, shift_end_feed, len, zz_peak);
        int copied_bottom = CopyBuffer(g_zig_handle, 1, shift_end_feed, len, zz_bottom);
        if(copied_main != len || copied_high != len || copied_low != len) return false;

        ArrayResize(main_ch, len); ArrayResize(high_ch, len); ArrayResize(low_ch, len);
        for(int j=0;j<len;j++)
        {
            int src = len-1-j; // cronológico
            double peak_v   = (copied_peak==len   ? zz_peak[src]   : 0.0);
            double bottom_v = (copied_bottom==len ? zz_bottom[src] : 0.0);
            double main_v   = zz_main[src];
            double pick = (peak_v!=0.0) ? peak_v : (bottom_v!=0.0 ? bottom_v : main_v);
            main_ch[j] = pick;
            high_ch[j] = zz_high[src];
            low_ch[j]  = zz_low[src];
        }
        return true;
    }

    bool BuildFeed(const int shift_end_feed,
                   const double &high[],
                   const double &low[],
                   ZIG_MODE mode)
    {
        int len = InpFFTWindow;
        static double main_ch[], high_ch[], low_ch[];
        if(!LoadWindow(shift_end_feed, len, main_ch, high_ch, low_ch))
            return false;

        int last_ext = -1; double last_val = 0.0;
        for(int k=0;k<len;k++){ if(main_ch[k]!=0.0){ last_ext=k; last_val=main_ch[k]; break; } }
        if(last_ext==-1) last_val = (high[0]+low[0])*0.5;

        for(int j=0; j<len; ++j)
        {
            double v=last_val;
            switch(mode)
            {
                case ZIG_STEP:
                    if(main_ch[j]!=0.0){ last_ext=j; last_val=main_ch[j]; }
                    v = last_val;
                    break;
                case ZIG_INTERP:
                    {
                        // interp apenas entre extremos confirmados
                        static int ext_pos[]; static double ext_val[];
                        ArrayResize(ext_pos,0); ArrayResize(ext_val,0);
                        for(int k=0;k<len;k++)
                            if(main_ch[k]!=0.0){
                                int sz=ArraySize(ext_pos);
                                ArrayResize(ext_pos,sz+1); ArrayResize(ext_val,sz+1);
                                ext_pos[sz]=k; ext_val[sz]=main_ch[k];
                            }
                        int n=ArraySize(ext_pos);
                        if(n==0){ v=last_val; }
                        else if(j<=ext_pos[0]) v=ext_val[0];
                        else if(j>=ext_pos[n-1]) v=ext_val[n-1];
                        else{
                            int kseg=-1;
                            for(int kk=0;kk<n-1;kk++) if(j>=ext_pos[kk] && j<ext_pos[kk+1]) { kseg=kk; break; }
                            if(kseg==-1) v=ext_val[n-1];
                            else{
                                int a=ext_pos[kseg], b=ext_pos[kseg+1];
                                double va=ext_val[kseg], vb=ext_val[kseg+1];
                                double t=(double)(j-a)/(double)(b-a);
                                v = va + (vb - va)*t;
                            }
                        }
                    }
                    break;
                case ZIG_MID:
                    v = (high_ch[j]+low_ch[j])*0.5;
                    break;
            }
            feed_data[j]=v;
        }
        return true;
    }
};

class FeedBuilder
{
public:
    bool EnsureCache(const int needed_bars)
    {
        if(InpFeedData!=FEED_PLA && InpFeedData!=FEED_CLOSE)
            return true; // cache só para close/PLA
        int delta=0; bool from_file=false;
        int need_adj = (InpPrefetchBars>0) ? MathMax(needed_bars, InpPrefetchBars) : needed_bars;
        bool ok = EnsureFeedCache(g_feed_cache, _Symbol, FeedTF(), need_adj, InpEnableFeedCache, "WaveSpecZZ", delta, from_file);
        if(!ok) return false;
        if(InpEnableFeedCache && g_feed_cache.loaded && !g_cache_use_logged)
        {
            string src = from_file ? "cache" : "fresh";
            PrintFormat("[WaveSpecZZ][CACHE] using %s cache size=%d need=%d tf=%s (+%d)", src, ArraySize(g_feed_cache.close), need_adj, EnumToString(FeedTF()), delta);
            g_cache_use_logged=true;
        }
        return true;
    }

    bool Build(const int shift_end_feed,
               const int start_pos,
               const datetime &time[],
               const double &close[],
               const double &high[],
               const double &low[])
    {
        if(InpFeedData==FEED_CLOSE)
        {
            int need = shift_end_feed + InpFFTWindow;
            need = MathMax(need, InpPrefetchBars);
            if(!EnsureCache(need)) return false;
            if(InpEnableFeedCache && g_feed_cache.loaded && !g_cache_use_logged)
            {
                string src = g_feed_cache.from_file ? "file" : "fresh";
                PrintFormat("[WaveSpecZZ][CACHE] using %s cache size=%d need=%d tf=%s", src, ArraySize(g_feed_cache.close), need, EnumToString(FeedTF()));
                g_cache_use_logged=true;
            }
            for(int j=0;j<InpFFTWindow;j++)
            {
                int idx = shift_end_feed + (InpFFTWindow-1-j);
                feed_data[j]=g_feed_cache.close[idx];
            }
            return true;
        }
        if(InpFeedData==FEED_PLA)
        {
            // log único já feito na EnsureCache; evitar spam
            return BuildPlaPriceSeries(shift_end_feed);
        }
        // ZigZag
        return zig.BuildFeed(shift_end_feed, high, low, InpZigZagMode);
    }
private:
    ZigZagFeed zig;
};

class FftProcessor
{
public:
    bool Ensure(int len)
    {
        return EnsureGpu(len);
    }
    bool Run(const double &in[], int len)
    {
        int st = gpu_fft_real_forward(in, len, g_fft_interleaved);
        if(st!=ALGLIB_STATUS_OK) return false;
        int bins=len/2;
        for(int k=0;k<bins;k++)
        {
            int base=2*k;
            fft_real[k]=g_fft_interleaved[base];
            fft_imag[k]=(base+1<len)?g_fft_interleaved[base+1]:0.0;
        }
        for(int k=0;k<bins;k++)
            spectrum[k]=fft_real[k]*fft_real[k]+fft_imag[k]*fft_imag[k];
        return true;
    }
};

class ViewRouter
{
public:
    void ShowFeedOnly(int i)
    {
        FeedTrace[i] = feed_data[InpFFTWindow-1];
        WaveBuffer1[i]=WaveBuffer2[i]=WaveBuffer3[i]=WaveBuffer4[i]=EMPTY_VALUE;
        WaveBuffer5[i]=WaveBuffer6[i]=WaveBuffer7[i]=WaveBuffer8[i]=EMPTY_VALUE;
        WavePeriod1[i]=WavePeriod2[i]=WavePeriod3[i]=WavePeriod4[i]=EMPTY_VALUE;
        WavePeriod5[i]=WavePeriod6[i]=WavePeriod7[i]=WavePeriod8[i]=EMPTY_VALUE;
        ClearForecastBuffersAt(i);
    }
    void HideFeed(int i){ FeedTrace[i]=EMPTY_VALUE; }
};

static FeedBuilder   s_feed;
static FftProcessor  s_fft;
static ViewRouter    s_view;

int OnInit()
{
    Print("[WaveSpecZZ] OnInit start");
    ENUM_TIMEFRAMES tf_feed = (ENUM_TIMEFRAMES)Period(); // feed sempre no timeframe do gráfico
    if(InpEnableFeedCache)
    {
        string cache_file = FeedCacheFileName("WaveSpecZZ", _Symbol, tf_feed);
        Print("[WaveSpecZZ][CACHE] feed cache enabled (file=", cache_file, ") prefetch_bars=", InpPrefetchBars);
    }
    PrintFormat("[WaveSpecZZ][MODE] %s (importing tester.dll)", MQLInfoInteger(MQL_TESTER) ? "Strategy Tester" : "Live/visual");
    // Handle ZigZag somente se esse modo de feed for usado
    if(InpFeedData == FEED_ZIGZAG)
    {
        g_zig_handle = iCustom(_Symbol, tf_feed, "ZigZag", InpZigZagDepth, InpZigZagDeviation, InpZigZagBackstep);
        if(g_zig_handle == INVALID_HANDLE)
        {
            Print("[WaveSpecZZ][WARN] ZigZag handle inválido em OnInit; feed ZigZag ficará indisponível");
            // Não cancelar teste: permaneceremos, mas Build falhará e pulará barras até obter handle
        }
    }

    // Waves (plots 0-7)
    SetIndexBuffer(0, WaveBuffer1, INDICATOR_DATA);
    SetIndexBuffer(1, WaveBuffer2, INDICATOR_DATA);
    SetIndexBuffer(2, WaveBuffer3, INDICATOR_DATA);
    SetIndexBuffer(3, WaveBuffer4, INDICATOR_DATA);
    SetIndexBuffer(4, WaveBuffer5, INDICATOR_DATA);
    SetIndexBuffer(5, WaveBuffer6, INDICATOR_DATA);
    SetIndexBuffer(6, WaveBuffer7, INDICATOR_DATA);
    SetIndexBuffer(7, WaveBuffer8, INDICATOR_DATA);

    // Feed trace (plot 8)
    SetIndexBuffer(8, FeedTrace, INDICATOR_DATA);

    // Forecast markers UP (plots 10-17)
    SetIndexBuffer(10, ForecastMark1, INDICATOR_DATA);
    SetIndexBuffer(11, ForecastMark2, INDICATOR_DATA);
    SetIndexBuffer(12, ForecastMark3, INDICATOR_DATA);
    SetIndexBuffer(13, ForecastMark4, INDICATOR_DATA);
    SetIndexBuffer(14, ForecastMark5, INDICATOR_DATA);
    SetIndexBuffer(15, ForecastMark6, INDICATOR_DATA);
    SetIndexBuffer(16, ForecastMark7, INDICATOR_DATA);
    SetIndexBuffer(17, ForecastMark8, INDICATOR_DATA);

    // Forecast markers DOWN (plots 18-25)
    SetIndexBuffer(18, ForecastDn1, INDICATOR_DATA);
    SetIndexBuffer(19, ForecastDn2, INDICATOR_DATA);
    SetIndexBuffer(20, ForecastDn3, INDICATOR_DATA);
    SetIndexBuffer(21, ForecastDn4, INDICATOR_DATA);
    SetIndexBuffer(22, ForecastDn5, INDICATOR_DATA);
    SetIndexBuffer(23, ForecastDn6, INDICATOR_DATA);
    SetIndexBuffer(24, ForecastDn7, INDICATOR_DATA);
    // plot 25 não existe (indicator_plots=25), então DN8 fica sem plot visível

    // Period buffers (slots 25-32) - expostos para EA (período ou countdown). Não têm plot.
    SetIndexBuffer(25, WavePeriod1, INDICATOR_CALCULATIONS);
    SetIndexBuffer(26, WavePeriod2, INDICATOR_CALCULATIONS);
    SetIndexBuffer(27, WavePeriod3, INDICATOR_CALCULATIONS);
    SetIndexBuffer(28, WavePeriod4, INDICATOR_CALCULATIONS);
    SetIndexBuffer(29, WavePeriod5, INDICATOR_CALCULATIONS);
    SetIndexBuffer(30, WavePeriod6, INDICATOR_CALCULATIONS);
    SetIndexBuffer(31, WavePeriod7, INDICATOR_CALCULATIONS);
    SetIndexBuffer(32, WavePeriod8, INDICATOR_CALCULATIONS);

    // ETA countdown (slots 33-40)
    SetIndexBuffer(33, EtaCount1, INDICATOR_DATA);
    SetIndexBuffer(34, EtaCount2, INDICATOR_DATA);
    SetIndexBuffer(35, EtaCount3, INDICATOR_DATA);
    SetIndexBuffer(36, EtaCount4, INDICATOR_DATA);
    SetIndexBuffer(37, EtaCount5, INDICATOR_DATA);
    SetIndexBuffer(38, EtaCount6, INDICATOR_DATA);
    SetIndexBuffer(39, EtaCount7, INDICATOR_DATA);
    SetIndexBuffer(40, EtaCount8, INDICATOR_DATA);

    // Fase (radianos) (slots 41-48)
    SetIndexBuffer(41, PhaseVal1, INDICATOR_DATA);
    SetIndexBuffer(42, PhaseVal2, INDICATOR_DATA);
    SetIndexBuffer(43, PhaseVal3, INDICATOR_DATA);
    SetIndexBuffer(44, PhaseVal4, INDICATOR_DATA);
    SetIndexBuffer(45, PhaseVal5, INDICATOR_DATA);
    SetIndexBuffer(46, PhaseVal6, INDICATOR_DATA);
    SetIndexBuffer(47, PhaseVal7, INDICATOR_DATA);
    SetIndexBuffer(48, PhaseVal8, INDICATOR_DATA);

    // Atributos MUSIC (cálculo / DataWindow via iCustom)
    SetIndexBuffer(54, MusEnergy1,   INDICATOR_CALCULATIONS);
    SetIndexBuffer(55, MusEnergy2,   INDICATOR_CALCULATIONS);
    SetIndexBuffer(56, MusCoher1,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(57, MusCoher2,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(58, MusSnrDb1,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(59, MusSnrDb2,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(60, MusScore1,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(61, MusScore2,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(62, MusEigen1,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(63, MusEigen2,    INDICATOR_CALCULATIONS);
    SetIndexBuffer(64, MusEtaConf1,  INDICATOR_CALCULATIONS);
    SetIndexBuffer(65, MusEtaConf2,  INDICATOR_CALCULATIONS);

    // ZigZag buffers (slots 49-53)
    SetIndexBuffer(49, ZigzagPeakBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(50, ZigzagBottomBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(51, ColorBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(52, HighMapBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(53, LowMapBuffer,INDICATOR_CALCULATIONS);

    // Exibição exclusiva: waves OU feed
    if(InpViewMode == VIEW_FEED)
    {
        for(int p=0; p<8; p++) PlotIndexSetInteger(p, PLOT_DRAW_TYPE, DRAW_NONE);
        PlotIndexSetInteger(8, PLOT_DRAW_TYPE, DRAW_LINE);
    }
    else // VIEW_WAVES
    {
        PlotIndexSetInteger(8, PLOT_DRAW_TYPE, DRAW_NONE);
        for(int p=0; p<8; p++) PlotIndexSetInteger(p, PLOT_DRAW_TYPE, DRAW_LINE);
    }

    // Labels e Data Window habilitados
    const string wave_labels[8] = {"Wave1","Wave2","Wave3","Wave4","Wave5","Wave6","Wave7","Wave8"};
    for(int p=0; p<8; p++)
    {
        PlotIndexSetString(p, PLOT_LABEL, wave_labels[p]);
        PlotIndexSetInteger(p, PLOT_SHOW_DATA, true);
    }
    PlotIndexSetString(8, PLOT_LABEL, "Feed");
    PlotIndexSetInteger(8, PLOT_SHOW_DATA, true);
    // Forecast markers também visíveis na Data Window (ETA em barras)
    const string fmark_labels[8] = {"FMark1","FMark2","FMark3","FMark4","FMark5","FMark6","FMark7","FMark8"};
    for(int p=0; p<8; p++)
    {
        PlotIndexSetString(10+p, PLOT_LABEL, fmark_labels[p]);
        PlotIndexSetInteger(10+p, PLOT_SHOW_DATA, true);
    }
    IndicatorSetString(INDICATOR_SHORTNAME, "WaveSpecZZ 1.1.0");

    // fila de jobs assíncronos
    g_cycles_job_capacity = MathMax(1, InpAsyncDepth);
    ArrayResize(g_cycles_jobs, g_cycles_job_capacity);
    g_cycles_job_count = 0;

    Print("[WaveSpecZZ] OnInit ok");
    return(INIT_SUCCEEDED);
}

void OnDeinit(const int reason)
{
    PrintFormat("[WaveSpecZZ][DEINIT] reason=%d", reason);
    PrintFormat("[WaveSpecZZ] OnDeinit reason=%d", reason);
    if(g_zig_handle != INVALID_HANDLE)
        IndicatorRelease(g_zig_handle);

    // libera todos os jobs pendentes
    for(int j=0; j<g_cycles_job_count; ++j)
        gpu_free_job(g_cycles_jobs[j]);
    g_cycles_job_count = 0;
    g_cycles_job_id = 0;
    if(g_gpu_session)
    {
        Print("[WaveSpecZZ][DEINIT] shutting down GPU session");
        gpu_shutdown();
        g_gpu_session=false;
        g_gpu_init_logged=false;
    }

    ArrayInitialize(EtaCountdown, 0.0);
    g_cycle_cache_loaded=false;
}

bool EnsureGpu(int length)
{
    if(length<=0) return false;
    if(g_gpu_session && length==ArraySize(g_fft_interleaved)) return true;
    if(!g_gpu_session)
    {
        // Core aceita até 512 streams; aplica piso mínimo para evitar valores muito baixos herdados do tester
        int streams = MathMax(16, MathMin(512, InpGpuStreams));
        if(!g_gpu_init_logged)
        {
            PrintFormat("[WaveSpecZZ][GPU] gpu_init start streams=%d", streams);
            g_gpu_init_logged=true;
        }
        int st = gpu_init(0, streams);
        if(st!=ALGLIB_STATUS_OK)
        {
            g_gpu_init_fail_counter++;
            g_gpu_init_logged=false; // permite novo log ao tentar novamente
            if((g_gpu_init_fail_counter % 50)==1 && (TimeCurrent()-g_gpu_init_last_log)>=5)
            {
                ushort wbuf[256];
                int n = gpu_get_last_error_w(wbuf, ArraySize(wbuf));
                string reason = (n>0 ? ShortArrayToString(wbuf, 0, n-1) : "n/a");
                PrintFormat("[WaveSpecZZ][ERR] gpu_init failed st=%d (every 50 bars) reason=%s", st, reason);
                g_gpu_init_last_log = TimeCurrent();
            }
            return false;
        }
        g_gpu_init_fail_counter = 0; // reset after success
        g_gpu_session=true;
    }
    ArrayResize(g_fft_interleaved, length);
    ArrayResize(fft_real, length);
    ArrayResize(fft_imag, length);
    return true;
}

// PLA builder (sem fallback) alinhado ao timeframe de feed
bool BuildPlaPriceSeries(const int shift_end_feed)
{
    if(shift_end_feed<0) return false;
    int need = shift_end_feed + InpFFTWindow;
    if(!s_feed.EnsureCache(need)) return false;
    for(int j=0;j<InpFFTWindow;j++)
    {
        int idx = shift_end_feed + (InpFFTWindow-1-j);
        feed_data[j]=g_feed_cache.close[idx];
    }
    return true;
}

// ZigZag price series builder (3 modos: STEP, INTERP, MID) no feed timeframe
bool LoadZigZagWindow(const int shift_end_feed,int len,double &main_ch[],double &high_ch[],double &low_ch[])
{
    static double zz_main[], zz_high[], zz_low[], zz_peak[], zz_bottom[];
    ArraySetAsSeries(zz_main, true);
    ArraySetAsSeries(zz_high, true);
    ArraySetAsSeries(zz_low,  true);
    ArraySetAsSeries(zz_peak, true);
    ArraySetAsSeries(zz_bottom, true);
    ArrayResize(zz_main, len); ArrayResize(zz_high, len); ArrayResize(zz_low, len);
    ArrayResize(zz_peak, len); ArrayResize(zz_bottom, len);

    int copied_main   = CopyBuffer(g_zig_handle, 0, shift_end_feed, len, zz_main);
    int copied_high   = CopyBuffer(g_zig_handle, 1, shift_end_feed, len, zz_high);
    int copied_low    = CopyBuffer(g_zig_handle, 2, shift_end_feed, len, zz_low);
    int copied_peak   = CopyBuffer(g_zig_handle, 0, shift_end_feed, len, zz_peak);   // picos (separados)
    int copied_bottom = CopyBuffer(g_zig_handle, 1, shift_end_feed, len, zz_bottom); // fundos (separados)
    if(copied_main != len || copied_high != len || copied_low != len)
        return false;

    ArrayResize(main_ch, len); ArrayResize(high_ch, len); ArrayResize(low_ch, len);
    for(int j=0;j<len;j++)
    {
        int src = len-1-j; // reverte para cronológico
        double peak_v   = (copied_peak==len   ? zz_peak[src]   : 0.0);
        double bottom_v = (copied_bottom==len ? zz_bottom[src] : 0.0);
        double main_v   = zz_main[src];
        double pick = (peak_v!=0.0) ? peak_v : (bottom_v!=0.0 ? bottom_v : main_v);
        main_ch[j] = pick;
        high_ch[j] = zz_high[src];
        low_ch[j]  = zz_low[src];
    }
    return true;
}

bool BuildZigZagPriceSeries(const int shift_end_feed,
                            const double &high[],
                            const double &low[],
                            const datetime &time[],
                            ZIG_MODE mode)
{
    if(shift_end_feed < 0) return false;

    int len = InpFFTWindow;
    static double main_ch[], high_ch[], low_ch[];
    if(!LoadZigZagWindow(shift_end_feed, len, main_ch, high_ch, low_ch))
        return false;

    int last_ext = -1;
    double last_val = 0.0;
    // inicializar last_ext com primeiro extremo à frente, se existir
    for(int k=0;k<len;k++){ if(main_ch[k]!=0.0){ last_ext=k; last_val=main_ch[k]; break; } }
    if(last_ext==-1) last_val = (high[0]+low[0])*0.5; // fallback: preço atual do gráfico

    for(int j=0; j<len; ++j)
    {
        double v=last_val;
        switch(mode)
        {
            case ZIG_STEP:
                if(main_ch[j]!=0.0){ last_ext=j; last_val=main_ch[j]; }
                v = last_val;
                break;
            case ZIG_INTERP:
                {
                    // Interpola apenas entre extremos já confirmados (prev->curr), mantendo valor após o último
                    static int ext_pos[]; static double ext_val[];
                    ArrayResize(ext_pos, 0); ArrayResize(ext_val, 0);
                    for(int k=0;k<len;k++)
                    {
                        if(main_ch[k]!=0.0)
                        {
                            int sz = ArraySize(ext_pos);
                            ArrayResize(ext_pos, sz+1);
                            ArrayResize(ext_val, sz+1);
                            ext_pos[sz]=k;
                            ext_val[sz]=main_ch[k];
                        }
                    }
                    int n = ArraySize(ext_pos);
                    if(n==0){ v = last_val; }
                    else if(j <= ext_pos[0]) { v = ext_val[0]; }
                    else if(j >= ext_pos[n-1]) { v = ext_val[n-1]; }
                    else
                    {
                        // encontra o segmento ext_pos[k] <= j < ext_pos[k+1]
                        int kseg = -1;
                        for(int kidx=0;kidx<n-1;kidx++)
                        {
                            if(j >= ext_pos[kidx] && j < ext_pos[kidx+1]) { kseg = kidx; break; }
                        }
                        if(kseg==-1){ v = ext_val[n-1]; }
                        else {
                            int a = ext_pos[kseg], b = ext_pos[kseg+1];
                            double va = ext_val[kseg], vb = ext_val[kseg+1];
                            double t = (double)(j - a) / (double)(b - a);
                            v = va + (vb - va)*t;
                        }
                    }
                }
                break;
            case ZIG_MID:
                v = (high_ch[j]+low_ch[j])*0.5;
                break;
        }
        feed_data[j] = v;
    }
    return true;
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
    ENUM_TIMEFRAMES InpFeedTimeframe_local = FeedTF(); // feed sempre no timeframe atual
    if(rates_total < InpFFTWindow) return prev_calculated;

    ArrayResize(feed_data, InpFFTWindow);
    ArrayResize(detrended_data, InpFFTWindow);
    ArrayResize(spectrum, InpFFTWindow/2);
    ArrayResize(WaveKalman, rates_total);
    ArrayResize(FeedTrace, rates_total);
    // Garantir capacidade dos buffers de plot (defensivo)
    if(ArraySize(WaveBuffer1) < rates_total) { ArrayResize(WaveBuffer1, rates_total); ArrayResize(WaveBuffer2, rates_total); ArrayResize(WaveBuffer3, rates_total); ArrayResize(WaveBuffer4, rates_total); ArrayResize(WaveBuffer5, rates_total); ArrayResize(WaveBuffer6, rates_total); ArrayResize(WaveBuffer7, rates_total); ArrayResize(WaveBuffer8, rates_total); }
    if(ArraySize(WavePeriod1) < rates_total) { ArrayResize(WavePeriod1, rates_total); ArrayResize(WavePeriod2, rates_total); ArrayResize(WavePeriod3, rates_total); ArrayResize(WavePeriod4, rates_total); ArrayResize(WavePeriod5, rates_total); ArrayResize(WavePeriod6, rates_total); ArrayResize(WavePeriod7, rates_total); ArrayResize(WavePeriod8, rates_total); }


    // Garantir capacidade dos buffers de forecast (tamanho >= rates_total)
    if(ArraySize(ForecastMark1) < rates_total)
    {
        ArrayResize(ForecastMark1, rates_total);
        ArrayResize(ForecastMark2, rates_total);
        ArrayResize(ForecastMark3, rates_total);
        ArrayResize(ForecastMark4, rates_total);
        ArrayResize(ForecastMark5, rates_total);
        ArrayResize(ForecastMark6, rates_total);
        ArrayResize(ForecastMark7, rates_total);
        ArrayResize(ForecastMark8, rates_total);
        ArrayInitialize(ForecastMark1, EMPTY_VALUE);
        ArrayInitialize(ForecastMark2, EMPTY_VALUE);
        ArrayInitialize(ForecastMark3, EMPTY_VALUE);
        ArrayInitialize(ForecastMark4, EMPTY_VALUE);
        ArrayInitialize(ForecastMark5, EMPTY_VALUE);
        ArrayInitialize(ForecastMark6, EMPTY_VALUE);
        ArrayInitialize(ForecastMark7, EMPTY_VALUE);
        ArrayInitialize(ForecastMark8, EMPTY_VALUE);
    }
    if(ArraySize(ForecastDn1) < rates_total)
    {
        ArrayResize(ForecastDn1, rates_total);
        ArrayResize(ForecastDn2, rates_total);
        ArrayResize(ForecastDn3, rates_total);
        ArrayResize(ForecastDn4, rates_total);
        ArrayResize(ForecastDn5, rates_total);
        ArrayResize(ForecastDn6, rates_total);
        ArrayResize(ForecastDn7, rates_total);
        ArrayResize(ForecastDn8, rates_total);
        ArrayInitialize(ForecastDn1, EMPTY_VALUE);
        ArrayInitialize(ForecastDn2, EMPTY_VALUE);
        ArrayInitialize(ForecastDn3, EMPTY_VALUE);
        ArrayInitialize(ForecastDn4, EMPTY_VALUE);
        ArrayInitialize(ForecastDn5, EMPTY_VALUE);
        ArrayInitialize(ForecastDn6, EMPTY_VALUE);
        ArrayInitialize(ForecastDn7, EMPTY_VALUE);
        ArrayInitialize(ForecastDn8, EMPTY_VALUE);
    }
    if(ArraySize(EtaCount1) < rates_total)
    {
        ArrayResize(EtaCount1, rates_total);
        ArrayResize(EtaCount2, rates_total);
        ArrayResize(EtaCount3, rates_total);
        ArrayResize(EtaCount4, rates_total);
        ArrayResize(EtaCount5, rates_total);
        ArrayResize(EtaCount6, rates_total);
        ArrayResize(EtaCount7, rates_total);
        ArrayResize(EtaCount8, rates_total);
    }
    if(ArraySize(PhaseVal1) < rates_total)
    {
        ArrayResize(PhaseVal1, rates_total);
        ArrayResize(PhaseVal2, rates_total);
        ArrayResize(PhaseVal3, rates_total);
        ArrayResize(PhaseVal4, rates_total);
        ArrayResize(PhaseVal5, rates_total);
        ArrayResize(PhaseVal6, rates_total);
        ArrayResize(PhaseVal7, rates_total);
        ArrayResize(PhaseVal8, rates_total);
    }

    if(ArraySize(MusEnergy1) < rates_total)
    {
        ArrayResize(MusEnergy1, rates_total);  ArrayResize(MusEnergy2, rates_total);
        ArrayResize(MusCoher1,  rates_total);  ArrayResize(MusCoher2,  rates_total);
        ArrayResize(MusSnrDb1,  rates_total);  ArrayResize(MusSnrDb2,  rates_total);
        ArrayResize(MusScore1,  rates_total);  ArrayResize(MusScore2,  rates_total);
        ArrayResize(MusEigen1,  rates_total);  ArrayResize(MusEigen2,  rates_total);
        ArrayResize(MusEtaConf1,rates_total);  ArrayResize(MusEtaConf2,rates_total);
    }

    // Warmup: garante cache mínimo de 1 janela antes de decidir o start (evita backfill longo se cache pronto)
    if(prev_calculated==0 && InpEnableFeedCache)
    {
        int d0=0; bool from0=false;
        int warm_need = (InpPrefetchBars>0) ? MathMax(InpFFTWindow, InpPrefetchBars) : InpFFTWindow;
        EnsureFeedCache(g_feed_cache, _Symbol, InpFeedTimeframe_local, warm_need, InpEnableFeedCache, "WaveSpecZZ", d0, from0);
    }

    // Tenta carregar cache; se InpForceBatch=true, ignora cache e refaz batch
    if(prev_calculated==0 && !g_cycle_cache_loaded)
    {
        bool cache_loaded = false;
        if(!InpForceBatch && LoadCycleCache(rates_total))
        {
            // já carregou buffers -> considera tudo calculado
            return(rates_total);
        }
        // Batch warmup sempre que habilitado ou se cache não estava disponível
        if(InpBatchWarmup)
        {
            Print("[WaveSpecZZ][BATCH][WARMUP] cache não encontrado; iniciando batch warmup...");
            // Garante backend carregado antes de chamar a API batch (senão retorna BACKEND_UNAVAILABLE)
            EnsureGpu(InpFFTWindow);
            // Prepara série de preços (close) no TF atual
            double prices[]; ArraySetAsSeries(prices, true);
            int total = Bars(_Symbol, FeedTF());
            if(total <= 0) total = rates_total;
            // Limite dedicado para warmup batch (pode ser bem menor para acelerar o attach)
            if(InpBatchBarsLimit > 0)
                total = MathMin(total, InpBatchBarsLimit);
            else if(InpMaxLiveBars > 0)
                total = MathMin(total, InpMaxLiveBars);
            int got = CopyClose(_Symbol, FeedTF(), 0, total, prices); // WHOLE_ARRAY limitado
            if(got >= InpFFTWindow)
            {
                const int hop = 1;
                const int stride = 15;
                int nwin = 1 + (got - InpFFTWindow) / hop;
                int buf_cap = nwin * InpGpuTopK * stride;
                PrintFormat("[WaveSpecZZ][BATCH][WARMUP] len=%d win=%d hop=%d nwin=%d topk=%d stride=%d buf_cap=%d", got, InpFFTWindow, hop, nwin, InpGpuTopK, stride, buf_cap);
                double cycles[]; ArrayResize(cycles, buf_cap);
                long jid=0; int st = gpu_submit_extract_cycles_batch(prices, got, InpFFTWindow, hop, InpGpuTopK,
                                                                     InpGpuMinPeriod, InpGpuMaxPeriod, (double)PeriodSeconds(FeedTF()),
                                                                     InpGpuMethod, InpGpuArOrder, stride, jid);
                if(st==ALGLIB_STATUS_OK && jid!=0)
                {
                    PrintFormat("[WaveSpecZZ][BATCH] submitted job=%I64d len=%d windows=%d", jid, got, nwin);
                    int ready=0,out_len=0;
                    ulong wait_start = GetTickCount64();
                    const ulong wait_limit = (InpBatchWaitMs > 0 ? (ulong)InpBatchWaitMs : 0);
                    while(true)
                    {
                        st = gpu_try_get_cycles_batch(jid, cycles, buf_cap, out_len, ready);
                        if(st==ALGLIB_STATUS_OK && ready==1)
                            break;
                        if(st!=ALGLIB_STATUS_OK && st!=ALGLIB_STATUS_NOT_READY)
                            break;
                        if(wait_limit > 0 && (GetTickCount64() - wait_start) >= wait_limit)
                            break;
                        Sleep(5);
                    }
                    gpu_free_job(jid);
                    if(st==ALGLIB_STATUS_OK && ready==1 && out_len>0)
                    {
                        PrintFormat("[WaveSpecZZ][BATCH][WARMUP] job ready cycles=%d", out_len);
                        // zera buffers
                        ArrayResize(WaveBuffer1, got); ArrayResize(WaveBuffer2, got);
                        ArrayResize(WavePeriod1, got); ArrayResize(WavePeriod2, got);
                        ArrayResize(EtaCount1, got); ArrayResize(EtaCount2, got);
                        ArrayResize(PhaseVal1, got); ArrayResize(PhaseVal2, got);
                        ArrayResize(MusEnergy1, got); ArrayResize(MusEnergy2, got);
                        ArrayResize(MusCoher1, got);  ArrayResize(MusCoher2, got);
                        ArrayResize(MusSnrDb1, got);  ArrayResize(MusSnrDb2, got);
                        ArrayResize(MusScore1, got);  ArrayResize(MusScore2, got);
                        ArrayResize(MusEigen1, got);  ArrayResize(MusEigen2, got);
                        ArrayResize(MusEtaConf1, got);ArrayResize(MusEtaConf2, got);
                        ArrayInitialize(WaveBuffer1, EMPTY_VALUE); ArrayInitialize(WaveBuffer2, EMPTY_VALUE);
                        ArrayInitialize(WavePeriod1, EMPTY_VALUE); ArrayInitialize(WavePeriod2, EMPTY_VALUE);
                        ArrayInitialize(EtaCount1, EMPTY_VALUE);  ArrayInitialize(EtaCount2, EMPTY_VALUE);
                        ArrayInitialize(PhaseVal1, EMPTY_VALUE);  ArrayInitialize(PhaseVal2, EMPTY_VALUE);
                        ArrayInitialize(MusEnergy1, EMPTY_VALUE); ArrayInitialize(MusEnergy2, EMPTY_VALUE);
                        ArrayInitialize(MusCoher1, EMPTY_VALUE);  ArrayInitialize(MusCoher2, EMPTY_VALUE);
                        ArrayInitialize(MusSnrDb1, EMPTY_VALUE);  ArrayInitialize(MusSnrDb2, EMPTY_VALUE);
                        ArrayInitialize(MusScore1, EMPTY_VALUE);  ArrayInitialize(MusScore2, EMPTY_VALUE);
                        ArrayInitialize(MusEigen1, EMPTY_VALUE);  ArrayInitialize(MusEigen2, EMPTY_VALUE);
                        ArrayInitialize(MusEtaConf1, EMPTY_VALUE);ArrayInitialize(MusEtaConf2, EMPTY_VALUE);

                        const double two_pi = 6.28318530717958647692;
                        for(int c=0; c<out_len; ++c)
                        {
                            int base = c*stride;
                            int method_id = (stride>14 ? (int)cycles[base+14] : 0);
                            if(InpMusicOnly && method_id!=1) continue;
                            double amp=cycles[base+0], freq=cycles[base+1], period=cycles[base+2], phase=cycles[base+3];
                            double eta_sec=cycles[base+5];
                            double energy=cycles[base+6], coher=cycles[base+7], snr=cycles[base+8], eigen=cycles[base+10], score=cycles[base+11], etac=cycles[base+13];

                            double w_energy = MathMax(energy,0.0);
                            double w_coher  = MathMax(coher,0.0);
                            double w_score  = MathMax(score,0.0);
                            double snr_eff  = MathMax(snr, InpMinSnrDb);
                            double w_snr    = 1.0 / (1.0 + MathPow(10.0, -snr_eff/10.0));
                            double weight_total = InpUseMusicWeights ? (w_energy*w_coher*w_score*w_snr) : 1.0;
                            if(coher < InpMinCoherence || score < InpMinScore) weight_total = 0.0;

                            int window_idx = c / InpGpuTopK;
                            int start_bar = window_idx * hop;
                            if(start_bar >= got) continue;
                            double omega = two_pi * freq;
                            int recon_span = MathMin(InpFFTWindow-1, got - start_bar - 1);
                            int slot = c % InpGpuTopK; // 0/1

                            for(int k=0; k<=recon_span; ++k)
                            {
                                int idx = start_bar + k;
                                double theta = phase - omega * k;
                                double val = amp * weight_total * MathSin(theta);
                                if(slot==0){ WaveBuffer1[idx]=val; WavePeriod1[idx]=period; EtaCount1[idx]=MathMax(eta_sec - k*PeriodSeconds(FeedTF()),0.0); PhaseVal1[idx]=theta; MusEnergy1[idx]=energy; MusCoher1[idx]=coher; MusSnrDb1[idx]=snr; MusScore1[idx]=score; MusEigen1[idx]=eigen; MusEtaConf1[idx]=etac; }
                                else       { WaveBuffer2[idx]=val; WavePeriod2[idx]=period; EtaCount2[idx]=MathMax(eta_sec - k*PeriodSeconds(FeedTF()),0.0); PhaseVal2[idx]=theta; MusEnergy2[idx]=energy; MusCoher2[idx]=coher; MusSnrDb2[idx]=snr; MusScore2[idx]=score; MusEigen2[idx]=eigen; MusEtaConf2[idx]=etac; }
                            }
                        }

                        g_cycle_cache_loaded = true; // evita recomputar no mesmo attach
                        SaveCycleCache(got);
                        PrintFormat("[WaveSpecZZ][BATCH] cache salvo (%d barras)", got);
                        return got;
                    }
                    else
                    {
                        ulong waited = GetTickCount64() - wait_start;
                        PrintFormat("[WaveSpecZZ][BATCH][ERR] st=%d ready=%d out=%d waited_ms=%I64u", st, ready, out_len, waited);
                    }
                }
                else
                {
                    PrintFormat("[WaveSpecZZ][BATCH][ERR] submit st=%d", st);
                }
            }
            else
            {
                PrintFormat("[WaveSpecZZ][BATCH][ERR] poucas barras=%d (win=%d)", got, InpFFTWindow);
            }
        }
    }

    // Log único de modo/parametrização (síncrono/assíncrono, streams, depth, wait)
if(!g_mode_logged)
{
    int streams_clamped = MathMax(16, MathMin(512, InpGpuStreams));
    PrintFormat("[WaveSpecZZ][MODECFG][ASYNC] submit=%s depth=%d streams=%d wait_loop_ms=5 fallback=sync_after_wait feed_tf=%s fft=%d topk=%d",
                (InpAsyncCycles ? "on" : "off"), InpAsyncDepth, streams_clamped, EnumToString(FeedTF()), InpFFTWindow, InpGpuTopK);
    g_mode_logged = true;
}

    // Limite de barras em tempo real (evita travar carregando 1M+); live_limit <= rates_total
    int live_limit = rates_total;
    if(InpMaxLiveBars > 0)
        live_limit = MathMin(rates_total, InpMaxLiveBars);

    // Processa em ordem decrescente (index 0 = barra atual, índices maiores = barras mais antigas)
    int start, end;

    if(prev_calculated==0)
    {
        int span = MathMax(InpFFTWindow, InpFFTWindow * InpBackfillWindows); // mínimo: 1 janela
        if(InpMaxProcessBars > 0) span = MathMin(span, InpMaxProcessBars);
        span = MathMin(span, live_limit);

        // Processa apenas o necessário: últimas 'span' barras (mais recentes)
        start = MathMin(live_limit - 1, span - 1); // índice mais antigo que será processado
        end   = 0;                                  // até a barra atual

        g_prog_start = start;
        g_prog_total = (start - end) + 1;
        g_prog_done  = false;
        g_prog_logged_done = false;
        PrintFormat("[WaveSpecZZ][PROG][WARMUP] span=%d windows=%d maxbars=%d start=%d end=%d live_limit=%d", span, InpBackfillWindows, InpMaxProcessBars, start, end, live_limit);
    }
    else
    {
        start = MathMax(MathMin(prev_calculated-1, live_limit-1), InpFFTWindow-1);
        end   = 0;

        g_prog_start = start;
        g_prog_total = (start - end) + 1;
        g_prog_done  = false;
        g_prog_logged_done = false;
    }

    // Garante histórico suficiente desde end (mais antigo que vamos usar) até a barra atual (time[0])
    // garante histórico até a barra mais antiga processada (start) menos a janela FFT
    datetime earliest_needed = (datetime)((long)time[start] - (long)(InpFFTWindow-1) * (long)PeriodSeconds(FeedTF()));
    if(earliest_needed < 0) earliest_needed = 0;
    HistorySelect(earliest_needed, time[0]);

    // Itera de barras mais antigas (start) até mais novas (end) em ordem decrescente
    for(int i=start; i>=end; i--)
    {
        int start_pos = i - InpFFTWindow + 1;
        if(start_pos<0) continue;

        // Índice (shift) no timeframe de feed para a barra corrente (fim da janela)
        int shift_end_feed = iBarShift(_Symbol, FeedTF(), time[i], false); // tolera barra próxima
        if(shift_end_feed < 0)
        {
            // tentativa de forçar carregamento da barra alvo
            static double preload_try[];
            int got = CopyClose(_Symbol, FeedTF(), time[i], 1, preload_try);
            shift_end_feed = iBarShift(_Symbol, FeedTF(), time[i], false);
        }
        if(shift_end_feed < 0)
        {
            PrintFormat("[WaveSpecZZ][ERR] i=%d shift(feed)=-1 tf=%s time=%s (sem barra no feed)", i, EnumToString(FeedTF()), TimeToString(time[i], TIME_DATE|TIME_SECONDS));
            continue;
        }

        // feed selection
        bool feed_ok = s_feed.Build(shift_end_feed, start_pos, time, close, high, low);
        if(!feed_ok)
        {
            PrintFormat("[WaveSpecZZ][ERR] Build feed failed i=%d shift=%d tf=%s", i, shift_end_feed, EnumToString(FeedTF()));
            continue;
        }

        // Log periódico do feed/timeframe + progresso de backfill
        if( ((i % kFeedLogEvery)==0 && i <= g_prog_start) )
        {
            string feed_mode = (InpFeedData==FEED_PLA ? "PLA" : (InpFeedData==FEED_ZIGZAG ? "ZIGZAG" : "CLOSE"));
            string view_mode = (InpViewMode==VIEW_FEED ? "FEED" : "WAVES");
            int done = (g_prog_start - i) + 1;
            double pct = (g_prog_total>0) ? (100.0 * done / g_prog_total) : 100.0;
            PrintFormat("[WaveSpecZZ][FEED][LIVE] i=%d tf=%s shift=%d mode=%s view=%s window=%d progress=%d/%d (%.1f%%)",
                        i, EnumToString(FeedTF()), shift_end_feed, feed_mode, view_mode, InpFFTWindow,
                        done, (g_prog_total>0 ? g_prog_total : done), pct);
        }
        // log de conclusão do backfill
        if(!g_prog_logged_done && g_prog_total>0 && (g_prog_start - i + 1)>=g_prog_total)
        {
            g_prog_logged_done = true;
            PrintFormat("[WaveSpecZZ][PROG] backfill done %d bars", g_prog_total);
            if(InpCycleCache && !g_cycle_cache_loaded)
                SaveCycleCache(g_prog_total);
        }

        // Exibição exclusiva: se for só feed, publicar feed e pular cálculo de ondas
        if(InpViewMode == VIEW_FEED){ s_view.ShowFeedOnly(i); continue; }
        s_view.HideFeed(i);

        // Decrementa ETA (contagem regressiva) por slot, se habilitado
        if(InpEtaCountdown)
        {
            for(int s=0; s<8; s++)
                if(EtaCountdown[s] > 0.0) EtaCountdown[s] -= 1.0;
        }

        ArrayCopy(detrended_data, feed_data, 0, 0, InpFFTWindow);

        // windowing: none (can add later)

        if(!s_fft.Ensure(InpFFTWindow))
        {
            if((i % 256)==0) Print("[WaveSpecZZ][GPU] EnsureGpu failed");
            continue; // sem GPU, não processa
        }

        if(!s_fft.Run(detrended_data, InpFFTWindow)) continue;

        // GPU: extrai ciclos já classificados (FFT ou MUSIC no core)
        const int stride = 15; // amplitude,freq,period,phase,eta_bars,eta_seconds,energy,coherence,snr,residual,eigen,score,kalman,eta_conf,method
        int needed = InpGpuTopK * stride;
        if(ArraySize(g_cycles_raw) < needed) ArrayResize(g_cycles_raw, needed);
        if(ArraySize(g_cycles_job_buf) < needed) ArrayResize(g_cycles_job_buf, needed);

        int cycles_out = 0;
        bool have_cycles = false;

        bool allow_async = InpAsyncCycles; // sempre permite async; garantimos preenchimento com fallback síncrono

        // sample_rate em segundos reais do timeframe do gráfico
        double sample_rate_sec = (double)PeriodSeconds(Period());
        if(sample_rate_sec <= 0) sample_rate_sec = 1.0;

        // poll async jobs (fila até InpAsyncDepth)
        if(allow_async && g_cycles_job_count>0)
        {
            for(int idx=0; idx<g_cycles_job_count; ++idx)
            {
                long jid = g_cycles_jobs[idx];
                int ready=0, out_len=0;
                int stp = gpu_try_get_cycles(jid, g_cycles_job_buf, stride, InpGpuTopK, out_len, ready);
                if(stp==ALGLIB_STATUS_NOT_READY)
                    continue;

                gpu_free_job(jid);

                // remove job da fila
                for(int m=idx; m<g_cycles_job_count-1; ++m)
                    g_cycles_jobs[m] = g_cycles_jobs[m+1];
                g_cycles_job_count--;

                if(stp==ALGLIB_STATUS_OK && ready==1 && out_len>0)
                {
                    ArrayCopy(g_cycles_raw, g_cycles_job_buf, 0, 0, out_len*stride);
                    cycles_out = out_len;
                    have_cycles = true;
                    if(InpDebugCycles)
                    {
                        string msg = StringFormat("[WaveSpecZZ][GPU][ASYNC] job ready cycles=%d i=%d", cycles_out, i);
                        for(int s=0; s<MathMin(cycles_out, InpGpuTopK); s++)
                        {
                            int base = s*stride;
                            double amp    = g_cycles_raw[base+0];
                            double period = g_cycles_raw[base+2];
                            double score  = g_cycles_raw[base+11];
                            msg += StringFormat(" |s%d amp=%.5f per=%.2f score=%.3f", s+1, amp, period, score);
                        }
                        Print(msg);
                    }
                }
                else if(InpDebugCycles)
                {
                    PrintFormat("[WaveSpecZZ][GPU][ASYNC][ERR] job=%I64d st=%d ready=%d out=%d", jid, stp, ready, out_len);
                }
                // processa no máximo um resultado por barra
                if(have_cycles) break;
            }
        }

        // submit new job se há capacidade livre
        if(allow_async && !have_cycles && InpAsyncDepth>0 && g_cycles_job_count < InpAsyncDepth)
        {
            long job_new = 0;
            int sts = gpu_submit_extract_cycles(detrended_data,
                                               InpFFTWindow,
                                               InpGpuTopK,
                                               InpGpuMinPeriod,
                                               InpGpuMaxPeriod,
                                               sample_rate_sec,
                                               InpGpuMethod,
                                               InpGpuArOrder,
                                               job_new);
            if(sts==ALGLIB_STATUS_OK && job_new!=0)
            {
                if(g_cycles_job_count < g_cycles_job_capacity)
                {
                    g_cycles_jobs[g_cycles_job_count++] = job_new;
                    if(InpDebugCycles && (i%256)==0)
                        PrintFormat("[WaveSpecZZ][GPU][ASYNC] submitted job=%I64d queue=%d/%d i=%d", job_new, g_cycles_job_count, g_cycles_job_capacity, i);
                }
                else
                {
                    // capacidade foi excedida (por segurança) -> libera o job
                    gpu_free_job(job_new);
                }
            }
        }

        // modo assíncrono: espera até o job entregar (sem fallback)
        if(allow_async && !have_cycles && g_cycles_job_count>0)
        {
            while(!have_cycles && g_cycles_job_count>0)
            {
                long jid = g_cycles_jobs[0];
                int ready=0, out_len=0;
                int stp = gpu_try_get_cycles(jid, g_cycles_job_buf, stride, InpGpuTopK, out_len, ready);
                if(stp==ALGLIB_STATUS_OK && ready==1 && out_len>0)
                {
                    ArrayCopy(g_cycles_raw, g_cycles_job_buf, 0, 0, out_len*stride);
                    cycles_out = out_len;
                    have_cycles = true;
                    gpu_free_job(jid);
                    for(int m=0; m<g_cycles_job_count-1; ++m) g_cycles_jobs[m]=g_cycles_jobs[m+1];
                    g_cycles_job_count--;
                    if(InpDebugCycles && (i%256)==0)
                        PrintFormat("[WaveSpecZZ][GPU][ASYNC] job ready cycles=%d i=%d", cycles_out, i);
                }
                else if(stp==ALGLIB_STATUS_NOT_READY)
                {
                    Sleep(1);
                }
                else
                {
                    if((i%256)==0)
                        PrintFormat("[WaveSpecZZ][GPU][ASYNC][ERR] job=%I64d st=%d ready=%d out=%d", jid, stp, ready, out_len);
                    gpu_free_job(jid);
                    for(int m=0; m<g_cycles_job_count-1; ++m) g_cycles_jobs[m]=g_cycles_jobs[m+1];
                    g_cycles_job_count--;
                    break;
                }
            }
        }

        // fallback apenas quando modo assíncrono está desligado
        if(!have_cycles)
        {
            if(!allow_async)
            {
                int st_sync = gpu_extract_cycles(detrended_data,
                                                 InpFFTWindow,
                                                 InpGpuTopK,
                                                 InpGpuMinPeriod,
                                                 InpGpuMaxPeriod,
                                                 sample_rate_sec,
                                                 InpGpuMethod,
                                                 InpGpuArOrder,
                                                 g_cycles_raw,
                                                 stride,
                                                 InpGpuTopK,
                                                 cycles_out);
                if(st_sync==ALGLIB_STATUS_OK && cycles_out>0)
                {
                    have_cycles = true;
                    if(InpDebugCycles && (i%256)==0)
                        PrintFormat("[WaveSpecZZ][GPU][SYNC] cycles=%d i=%d", cycles_out, i);
                }
                else
                {
                    if((i%256)==0)
                        PrintFormat("[WaveSpecZZ][GPU][SYNC_FAIL] st=%d cycles=%d len=%d method=%d ar=%d", st_sync, cycles_out, InpFFTWindow, InpGpuMethod, InpGpuArOrder);
                    continue;
                }
            }
            else
            {
                // Sem resultado async (job falhou). Prossegue para próxima barra para evitar travar.
                continue;
            }
        }

        if((i % 256)==0) // log periódico obrigatório para depuração
        {
            string msg = StringFormat("[WaveSpecZZ][GPU][DBG] i=%d cycles=%d method=%d ar=%d", i, cycles_out, InpGpuMethod, InpGpuArOrder);
            for(int s=0; s<MathMin(cycles_out, InpGpuTopK); s++)
            {
                int base = s*stride;
                double amp    = g_cycles_raw[base+0];
                double freq   = g_cycles_raw[base+1];
                double period = g_cycles_raw[base+2];
                double phase  = g_cycles_raw[base+3];
                double eta    = g_cycles_raw[base+4];
                double score  = g_cycles_raw[base+11];
                msg += StringFormat(" | s%d amp=%.5f per=%.2f freq=%.6f phase=%.2f eta=%.1f score=%.3f", s+1, amp, period, freq, phase, eta, score);
            }
            Print(msg);
        }

        // Garante capacidade mínima até i (defensivo contra barras recém-adicionadas)
        if(ArraySize(WaveBuffer1) <= i)
        {
            int new_sz = i+1;
            ArrayResize(WaveBuffer1, new_sz); ArrayResize(WaveBuffer2, new_sz);
            ArrayResize(WaveBuffer3, new_sz); ArrayResize(WaveBuffer4, new_sz);
            ArrayResize(WaveBuffer5, new_sz); ArrayResize(WaveBuffer6, new_sz);
            ArrayResize(WaveBuffer7, new_sz); ArrayResize(WaveBuffer8, new_sz);
            ArrayResize(WavePeriod1, new_sz); ArrayResize(WavePeriod2, new_sz);
            ArrayResize(WavePeriod3, new_sz); ArrayResize(WavePeriod4, new_sz);
            ArrayResize(WavePeriod5, new_sz); ArrayResize(WavePeriod6, new_sz);
            ArrayResize(WavePeriod7, new_sz); ArrayResize(WavePeriod8, new_sz);
            ArrayResize(MusEnergy1, new_sz);  ArrayResize(MusEnergy2, new_sz);
            ArrayResize(MusCoher1,  new_sz);  ArrayResize(MusCoher2,  new_sz);
            ArrayResize(MusSnrDb1,  new_sz);  ArrayResize(MusSnrDb2,  new_sz);
            ArrayResize(MusScore1,  new_sz);  ArrayResize(MusScore2,  new_sz);
            ArrayResize(MusEigen1,  new_sz);  ArrayResize(MusEigen2,  new_sz);
            ArrayResize(MusEtaConf1,new_sz);  ArrayResize(MusEtaConf2,new_sz);
        }

        // Limpa buffers da barra antes de plotar (evita linhas horizontais quando top_k<8)
        if(i < ArraySize(WaveBuffer1)) WaveBuffer1[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer2)) WaveBuffer2[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer3)) WaveBuffer3[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer4)) WaveBuffer4[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer5)) WaveBuffer5[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer6)) WaveBuffer6[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer7)) WaveBuffer7[i]=EMPTY_VALUE;
        if(i < ArraySize(WaveBuffer8)) WaveBuffer8[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod1)) WavePeriod1[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod2)) WavePeriod2[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod3)) WavePeriod3[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod4)) WavePeriod4[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod5)) WavePeriod5[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod6)) WavePeriod6[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod7)) WavePeriod7[i]=EMPTY_VALUE;
        if(i < ArraySize(WavePeriod8)) WavePeriod8[i]=EMPTY_VALUE;
        if(i < ArraySize(EtaCount1)) { EtaCount1[i]=EMPTY_VALUE; EtaCount2[i]=EMPTY_VALUE; EtaCount3[i]=EMPTY_VALUE; EtaCount4[i]=EMPTY_VALUE; EtaCount5[i]=EMPTY_VALUE; EtaCount6[i]=EMPTY_VALUE; EtaCount7[i]=EMPTY_VALUE; EtaCount8[i]=EMPTY_VALUE; }
        if(i < ArraySize(PhaseVal1)) { PhaseVal1[i]=EMPTY_VALUE; PhaseVal2[i]=EMPTY_VALUE; PhaseVal3[i]=EMPTY_VALUE; PhaseVal4[i]=EMPTY_VALUE; PhaseVal5[i]=EMPTY_VALUE; PhaseVal6[i]=EMPTY_VALUE; PhaseVal7[i]=EMPTY_VALUE; PhaseVal8[i]=EMPTY_VALUE; }
        if(i < ArraySize(MusEnergy1)) { MusEnergy1[i]=MusEnergy2[i]=EMPTY_VALUE; MusCoher1[i]=MusCoher2[i]=EMPTY_VALUE; MusSnrDb1[i]=MusSnrDb2[i]=EMPTY_VALUE; MusScore1[i]=MusScore2[i]=EMPTY_VALUE; MusEigen1[i]=MusEigen2[i]=EMPTY_VALUE; MusEtaConf1[i]=MusEtaConf2[i]=EMPTY_VALUE; }
        ArrayInitialize(EtaCountdown, 0.0);

        // Reconstrução avançada (MUSIC-only): plota no máx 2 waves, propagando fase
        int plotted = 0;
        for(int s=0; s<MathMin(cycles_out, InpGpuTopK) && plotted < 2; s++)
        {
            int base = s*stride;
            double amp        = g_cycles_raw[base+0];
            double freq       = g_cycles_raw[base+1];      // ciclos por barra
            double period     = g_cycles_raw[base+2];      // barras
            double phase      = g_cycles_raw[base+3];      // radianos
            double eta_bars   = g_cycles_raw[base+4];      // barras até horizonte
            double eta_sec    = g_cycles_raw[base+5];      // segundos até horizonte
            double energy     = g_cycles_raw[base+6];      // energia relativa (0..1)
            double coherence  = g_cycles_raw[base+7];      // coerência (0..1)
            double snr_db     = g_cycles_raw[base+8];      // SNR dB
            double eigen_rat  = g_cycles_raw[base+10];     // razão de autovalores
            double score      = g_cycles_raw[base+11];     // score composto
            double kalman_pr  = g_cycles_raw[base+12];     // predição Kalman 1-pass
            double eta_conf   = g_cycles_raw[base+13];     // confiança ETA
            int    method_id  = (stride>14 ? (int)g_cycles_raw[base+14] : 0);

            // Modo MUSIC-only: ignora FFT ridge
            if(method_id != 1)
                continue;

            int slot = plotted; // 0 ou 1

            // Fase contínua: retropropaga até o tamanho da janela (limitado por desempenho)
            // Peso opcional pela qualidade MUSIC
            double w_energy = MathMax(energy, 0.0);
            double w_coher  = MathMax(coherence, 0.0);
            double w_score  = MathMax(score, 0.0);
            double snr_eff_db = MathMax(snr_db, InpMinSnrDb);
            double w_snr    = 1.0 / (1.0 + MathPow(10.0, -snr_eff_db/10.0)); // mapeia dB -> [0,1)
            double weight_total = InpUseMusicWeights ? (w_energy * w_coher * w_score * w_snr) : 1.0;
            if(weight_total < 0.0) weight_total = 0.0;

            const double two_pi = 6.28318530717958647692;
            double omega = two_pi * freq; // radianos por barra
            int recon_span = (int)MathMin((int)MathRound(MathMax(eta_bars, 1.0)), MathMin(InpFFTWindow-1, 512));

            for(int k=0; k<=recon_span; ++k)
            {
                int idx = i - k;
                if(idx < 0) break;

                double theta = phase - omega * k;               // fase correspondente à barra idx
                double amp_w = amp * ( (InpUseMusicWeights && (coherence < InpMinCoherence || score < InpMinScore)) ? 0.0 : weight_total );
                double wave_value = (InpDrawMode == DRAW_SINE_RECON && period>0.0)
                                    ? amp_w * MathSin(theta)
                                    : amp_w;

                // Preenche buffers do slot alvo
                switch(slot)
                {
                    case 0:
                        if(idx < ArraySize(WaveBuffer1))  WaveBuffer1[idx]  = wave_value;
                        if(idx < ArraySize(WavePeriod1))  WavePeriod1[idx]  = period;
                        if(idx < ArraySize(EtaCount1))    EtaCount1[idx]    = MathMax(eta_sec - k*sample_rate_sec, 0.0);
                        if(idx < ArraySize(PhaseVal1))    PhaseVal1[idx]    = theta;
                        if(k==0 && idx < ArraySize(MusEnergy1)) { MusEnergy1[idx]=energy; MusCoher1[idx]=coherence; MusSnrDb1[idx]=snr_db; MusScore1[idx]=score; MusEigen1[idx]=eigen_rat; MusEtaConf1[idx]=eta_conf; }
                        break;
                    case 1:
                        if(idx < ArraySize(WaveBuffer2))  WaveBuffer2[idx]  = wave_value;
                        if(idx < ArraySize(WavePeriod2))  WavePeriod2[idx]  = period;
                        if(idx < ArraySize(EtaCount2))    EtaCount2[idx]    = MathMax(eta_sec - k*sample_rate_sec, 0.0);
                        if(idx < ArraySize(PhaseVal2))    PhaseVal2[idx]    = theta;
                        if(k==0 && idx < ArraySize(MusEnergy2)) { MusEnergy2[idx]=energy; MusCoher2[idx]=coherence; MusSnrDb2[idx]=snr_db; MusScore2[idx]=score; MusEigen2[idx]=eigen_rat; MusEtaConf2[idx]=eta_conf; }
                        break;
                }
            }

            // Marca de previsão apenas na barra atual usando horizonte em barras (e confiança mínima)
            bool forecast_ok = (eta_conf >= InpMinEtaConf);
            if(InpForecastMarks && eta_bars > 1 && forecast_ok)
            {
                int t_forecast = i + (int)MathRound(eta_bars);
                double wave_now = (InpDrawMode == DRAW_SINE_RECON && period>0.0)
                                  ? (amp * (InpUseMusicWeights ? weight_total : 1.0)) * MathSin(phase)
                                  : (amp * (InpUseMusicWeights ? weight_total : 1.0));
                if(slot==0)
                {
                    if(t_forecast >= ArraySize(ForecastMark1)) ArrayResize(ForecastMark1, t_forecast+1);
                    if(t_forecast <  ArraySize(ForecastMark1)) ForecastMark1[t_forecast] = wave_now;
                }
                else if(slot==1)
                {
                    if(t_forecast >= ArraySize(ForecastMark2)) ArrayResize(ForecastMark2, t_forecast+1);
                    if(t_forecast <  ArraySize(ForecastMark2)) ForecastMark2[t_forecast] = wave_now;
                }
            }

            plotted++;
        }

    }

    return(rates_total);
}
