// WaveSpecZZ minimal GPU+CPU FFT top-8 plotter (nodetrend, no ETA/leak/trackers)
#property copyright ""
#property link      ""
#property version   "1.100"
#property indicator_separate_window
#property indicator_buffers 22
#property indicator_plots   9

// Buffers ZigZag mínimo (picos/fundos) para construir feed (calculations only)
double ZigzagPeakBuffer[], ZigzagBottomBuffer[], ColorBuffer[], HighMapBuffer[], LowMapBuffer[];

#import "mt-bridge.dll"
  int  gpu_init(int device_index, int stream_count);
  void gpu_shutdown(void);
  int  gpu_fft_real_forward(const double &in[], int len, double &out[]);
#import

#define ALGLIB_STATUS_OK 0

// Inputs
input int  InpFFTWindow   = 16384;
input int  InpMinPeriod   = 18;
input int  InpMaxPeriod   = 200;

enum FEED_DATA_MODE { FEED_PLA = 0, FEED_ZIGZAG = 1, FEED_CLOSE = 2 };
input FEED_DATA_MODE InpFeedData = FEED_PLA;
input ENUM_TIMEFRAMES InpFeedTimeframe = PERIOD_M1;

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
input int    InpGpuTopK        = 2;        // ciclos extraídos da GPU
input int    InpGpuMethod      = 1;        // 0=FFT ridge, 1=MUSIC/ESPRIT, -1=auto
input int    InpGpuArOrder     = 8;        // ordem AR para MUSIC/ESPRIT

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
#property indicator_color2  clrOrange
#property indicator_width2  1
#property indicator_label3  "Wave3"
#property indicator_type3   DRAW_LINE
#property indicator_color3  clrYellow
#property indicator_width3  1
#property indicator_label4  "Wave4"
#property indicator_type4   DRAW_LINE
#property indicator_color4  clrGreen
#property indicator_width4  1
#property indicator_label5  "Wave5"
#property indicator_type5   DRAW_LINE
#property indicator_color5  clrAqua
#property indicator_width5  1
#property indicator_label6  "Wave6"
#property indicator_type6   DRAW_LINE
#property indicator_color6  clrBlue
#property indicator_width6  1
#property indicator_label7  "Wave7"
#property indicator_type7   DRAW_LINE
#property indicator_color7  clrFuchsia
#property indicator_width7  1
#property indicator_label8  "Wave8"
#property indicator_type8   DRAW_LINE
#property indicator_color8  clrGray
#property indicator_width8  1

#property indicator_label9  "FeedTrace"
#property indicator_type9   DRAW_LINE
#property indicator_color9  clrWhite
#property indicator_style9  STYLE_DOT
#property indicator_width9  1


double WaveBuffer1[],WaveBuffer2[],WaveBuffer3[],WaveBuffer4[];
double WaveBuffer5[],WaveBuffer6[],WaveBuffer7[],WaveBuffer8[];
double WavePeriod1[],WavePeriod2[],WavePeriod3[],WavePeriod4[];
double WavePeriod5[],WavePeriod6[],WavePeriod7[],WavePeriod8[];
double WaveKalman[];
double FeedTrace[];

double feed_data[], detrended_data[];
double g_fft_interleaved[], fft_real[], fft_imag[], spectrum[];
bool g_gpu_session = false;
int g_zig_handle = INVALID_HANDLE;

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
    bool Build(const int shift_end_feed,
               const int start_pos,
               const datetime &time[],
               const double &close[],
               const double &high[],
               const double &low[])
    {
        if(InpFeedData==FEED_CLOSE)
        {
            static double buf[]; ArrayResize(buf, InpFFTWindow); ArraySetAsSeries(buf,true);
            if(CopyClose(_Symbol, InpFeedTimeframe, shift_end_feed, InpFFTWindow, buf)!=InpFFTWindow) return false;
            for(int j=0;j<InpFFTWindow;j++) feed_data[j]=buf[InpFFTWindow-1-j];
            return true;
        }
        if(InpFeedData==FEED_PLA)
            return BuildPlaPriceSeries(shift_end_feed);
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
    }
    void HideFeed(int i){ FeedTrace[i]=EMPTY_VALUE; }
};

static FeedBuilder   s_feed;
static FftProcessor  s_fft;
static ViewRouter    s_view;

int OnInit()
{
    // Handle ZigZag para feed (usar indicador padrão ZigZag)
    g_zig_handle = iCustom(_Symbol, InpFeedTimeframe, "ZigZag", InpZigZagDepth, InpZigZagDeviation, InpZigZagBackstep);
    if(InpFeedData == FEED_ZIGZAG && g_zig_handle == INVALID_HANDLE)
        return(INIT_FAILED);

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

    // Period buffers (calc slots 9-16)
    SetIndexBuffer(9, WavePeriod1, INDICATOR_CALCULATIONS);
    SetIndexBuffer(10, WavePeriod2, INDICATOR_CALCULATIONS);
    SetIndexBuffer(11, WavePeriod3, INDICATOR_CALCULATIONS);
    SetIndexBuffer(12, WavePeriod4, INDICATOR_CALCULATIONS);
    SetIndexBuffer(13, WavePeriod5, INDICATOR_CALCULATIONS);
    SetIndexBuffer(14, WavePeriod6, INDICATOR_CALCULATIONS);
    SetIndexBuffer(15, WavePeriod7, INDICATOR_CALCULATIONS);
    SetIndexBuffer(16, WavePeriod8, INDICATOR_CALCULATIONS);

    // ZigZag buffers (calc slots 17-21)
    SetIndexBuffer(17, ZigzagPeakBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(18, ZigzagBottomBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(19, ColorBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(20, HighMapBuffer,INDICATOR_CALCULATIONS);
    SetIndexBuffer(21, LowMapBuffer,INDICATOR_CALCULATIONS);

    // Exibição exclusiva: waves OU feed
    if(InpViewMode == VIEW_FEED)
    {
        for(int p=0; p<8; p++) PlotIndexSetInteger(p, PLOT_DRAW_TYPE, DRAW_NONE);
        PlotIndexSetInteger(8, PLOT_DRAW_TYPE, DRAW_LINE);
    }
    else // VIEW_WAVES
    {
        // Feed totalmente oculto
        PlotIndexSetInteger(8, PLOT_DRAW_TYPE, DRAW_NONE);
        for(int p=0; p<8; p++) PlotIndexSetInteger(p, PLOT_DRAW_TYPE, DRAW_LINE);
    }

    // Remover rótulos/valores no Data Window e na barra de título da subjanela
    for(int p=0; p<indicator_plots; p++)
    {
        PlotIndexSetString(p, PLOT_LABEL, "");
        PlotIndexSetInteger(p, PLOT_SHOW_DATA, false);
    }
    IndicatorSetString(INDICATOR_SHORTNAME, "WaveSpecZZ 1.1.0");

    return(INIT_SUCCEEDED);
}

void OnDeinit(const int reason)
{
    if(g_zig_handle != INVALID_HANDLE)
        IndicatorRelease(g_zig_handle);

    if(g_gpu_session)
    {
        gpu_shutdown();
        g_gpu_session=false;
    }
}

bool EnsureGpu(int length)
{
    if(length<=0) return false;
    if(g_gpu_session && length==ArraySize(g_fft_interleaved)) return true;
    if(!g_gpu_session)
    {
        if(gpu_init(0,2)!=ALGLIB_STATUS_OK) return false;
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
    static double feed_close_tf[]; ArrayResize(feed_close_tf, InpFFTWindow); ArraySetAsSeries(feed_close_tf, true);
    if(shift_end_feed<0) return false;
    if(CopyClose(_Symbol, InpFeedTimeframe, shift_end_feed, InpFFTWindow, feed_close_tf)!=InpFFTWindow) return false;
    for(int j=0;j<InpFFTWindow;j++) feed_data[j]=feed_close_tf[InpFFTWindow-1-j];
    // trivial PLA: copy as-is
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
    // Garante histórico mínimo no timeframe de feed para evitar janela vazia
    int need_feed = InpFFTWindow + 256;
    HistorySelect(TimeCurrent() - (long)need_feed * PeriodSeconds(InpFeedTimeframe), TimeCurrent());
    static double preload_buf[];
    int preload = CopyClose(_Symbol, InpFeedTimeframe, 0, need_feed, preload_buf);
    if(preload < InpFFTWindow)
    {
        PrintFormat("[WaveSpecZZ][WARN] histórico insuficiente no feed %s (%d/%d); aguardando download.",
                    EnumToString(InpFeedTimeframe), preload, InpFFTWindow);
        return prev_calculated; // não avança até completar
    }

    if(rates_total < InpFFTWindow) return prev_calculated;

    ArrayResize(feed_data, InpFFTWindow);
    ArrayResize(detrended_data, InpFFTWindow);
    ArrayResize(spectrum, InpFFTWindow/2);
    ArrayResize(WaveKalman, rates_total);
    ArrayResize(FeedTrace, rates_total);

    int start = MathMax(prev_calculated-1, InpFFTWindow-1);
    for(int i=start;i<rates_total;i++)
    {
        int start_pos = i - InpFFTWindow + 1;
        if(start_pos<0) continue;

        // Índice (shift) no timeframe de feed para a barra corrente (fim da janela)
        int shift_end_feed = iBarShift(_Symbol, InpFeedTimeframe, time[i]);
        if(shift_end_feed < 0) continue;

        // feed selection
        if(!s_feed.Build(shift_end_feed, start_pos, time, close, high, low))
            continue;

        // Exibição exclusiva: se for só feed, publicar feed e pular cálculo de ondas
        if(InpViewMode == VIEW_FEED){ s_view.ShowFeedOnly(i); continue; }
        s_view.HideFeed(i);

        ArrayCopy(detrended_data, feed_data, 0, 0, InpFFTWindow);

        // windowing: none (can add later)

        if(!s_fft.Ensure(InpFFTWindow))
            continue; // sem GPU, não processa

        if(!s_fft.Run(detrended_data, InpFFTWindow)) continue;

        int bins=InpFFTWindow/2;
        for(int k=0;k<bins;k++)
        {
            int base=2*k;
            fft_real[k]=g_fft_interleaved[base];
            fft_imag[k]=(base+1<InpFFTWindow)?g_fft_interleaved[base+1]:0.0;
        }

        int bins2 = InpFFTWindow/2;
        for(int k=0;k<bins2;k++)
            spectrum[k]=fft_real[k]*fft_real[k]+fft_imag[k]*fft_imag[k];

        // top-8 bins by power within period range
        double top_pow[8]; int top_bin[8];
        for(int s=0;s<8;s++){top_pow[s]=-1.0; top_bin[s]=-1;}
        int min_index = (int)MathCeil((double)InpFFTWindow / (double)InpMaxPeriod);
        int max_index = (int)MathFloor((double)InpFFTWindow / (double)InpMinPeriod);
        if(max_index>=bins) max_index=bins-1;
        for(int b=min_index;b<=max_index;b++)
        {
            double p=spectrum[b];
            // insert into top list
            for(int s=0;s<8;s++)
            {
                if(p>top_pow[s])
                {
                    for(int t=7;t>s;t--){top_pow[t]=top_pow[t-1]; top_bin[t]=top_bin[t-1];}
                    top_pow[s]=p; top_bin[s]=b; break;
                }
            }
        }

        // fill buffers slot a slot (MQL não suporta ponteiros de array)
        for(int s=0;s<8;s++)
        {
            double amp=0.0, period=0.0;
            if(top_bin[s]>0)
            {
                period = (double)InpFFTWindow / (double)top_bin[s];
                double mag = MathSqrt(top_pow[s]);
                // Reconstrução aproximada no último ponto da janela usando fase da FFT
                double phase = MathArctan2(fft_imag[top_bin[s]], fft_real[top_bin[s]]);
                double n = (double)(InpFFTWindow-1);
                amp = (mag / (double)InpFFTWindow) * MathCos(phase + 2.0*M_PI*(double)top_bin[s]*n/(double)InpFFTWindow);
            }
            switch(s)
            {
                case 0: WaveBuffer1[i]=amp; WavePeriod1[i]=period; break;
                case 1: WaveBuffer2[i]=amp; WavePeriod2[i]=period; break;
                case 2: WaveBuffer3[i]=amp; WavePeriod3[i]=period; break;
                case 3: WaveBuffer4[i]=amp; WavePeriod4[i]=period; break;
                case 4: WaveBuffer5[i]=amp; WavePeriod5[i]=period; break;
                case 5: WaveBuffer6[i]=amp; WavePeriod6[i]=period; break;
                case 6: WaveBuffer7[i]=amp; WavePeriod7[i]=period; break;
                case 7: WaveBuffer8[i]=amp; WavePeriod8[i]=period; break;
            }
        }

        // optional Kalman smoothing on the feed_data last point
        if(InpEnableKalman)
        {
            if(ArraySize(WaveKalman)<rates_total) ArrayResize(WaveKalman,rates_total);
            double meas = feed_data[InpFFTWindow-1];
            WaveKalman[i]=meas; // placeholder (Kalman removed)
        }
        else if(ArraySize(WaveKalman)>=rates_total)
        {
            WaveKalman[i]=0.0;
        }

    }

    return(rates_total);
}
