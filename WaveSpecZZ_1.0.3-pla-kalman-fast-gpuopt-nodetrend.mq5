// WaveSpecZZ minimal GPU+CPU FFT top-8 plotter (nodetrend, no ETA/leak/trackers)
#property copyright ""
#property link      ""
#property version   "1.005"
#property indicator_separate_window
#property indicator_buffers 16
#property indicator_plots   8

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

enum FEED_DATA_MODE { FEED_PLA = 0, FEED_ZIGZAG_CONTINUOUS = 1, FEED_ZIGZAG_ALTERNATING = 2, FEED_CLOSE = 3 };
input FEED_DATA_MODE InpFeedData = FEED_PLA;
input ENUM_TIMEFRAMES InpFeedTimeframe = PERIOD_M1;

input group "PLA"
input int    InpPlaMaxSegments = 32;
input double InpPlaMaxError    = 0.0005;

input group "ZigZag"
input int InpZigZagDepth    = 12;
input int InpZigZagDeviation= 5;
input int InpZigZagBackstep = 3;

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


double WaveBuffer1[],WaveBuffer2[],WaveBuffer3[],WaveBuffer4[];
double WaveBuffer5[],WaveBuffer6[],WaveBuffer7[],WaveBuffer8[];
double WavePeriod1[],WavePeriod2[],WavePeriod3[],WavePeriod4[];
double WavePeriod5[],WavePeriod6[],WavePeriod7[],WavePeriod8[];
double WaveKalman[];

double feed_data[], detrended_data[];
double g_fft_interleaved[], fft_real[], fft_imag[], spectrum[];
bool g_gpu_session = false;

int OnInit()
{
    SetIndexBuffer(0, WaveBuffer1, INDICATOR_DATA);
    SetIndexBuffer(1, WaveBuffer2, INDICATOR_DATA);
    SetIndexBuffer(2, WaveBuffer3, INDICATOR_DATA);
    SetIndexBuffer(3, WaveBuffer4, INDICATOR_DATA);
    SetIndexBuffer(4, WaveBuffer5, INDICATOR_DATA);
    SetIndexBuffer(5, WaveBuffer6, INDICATOR_DATA);
    SetIndexBuffer(6, WaveBuffer7, INDICATOR_DATA);
    SetIndexBuffer(7, WaveBuffer8, INDICATOR_DATA);

    SetIndexBuffer(8, WavePeriod1, INDICATOR_CALCULATIONS);
    SetIndexBuffer(9, WavePeriod2, INDICATOR_CALCULATIONS);
    SetIndexBuffer(10,WavePeriod3, INDICATOR_CALCULATIONS);
    SetIndexBuffer(11,WavePeriod4, INDICATOR_CALCULATIONS);
    SetIndexBuffer(12,WavePeriod5, INDICATOR_CALCULATIONS);
    SetIndexBuffer(13,WavePeriod6, INDICATOR_CALCULATIONS);
    SetIndexBuffer(14,WavePeriod7, INDICATOR_CALCULATIONS);
    SetIndexBuffer(15,WavePeriod8, INDICATOR_CALCULATIONS);

    return(INIT_SUCCEEDED);
}

void OnDeinit(const int reason)
{
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

// Simple CPU FFT fallback (radix-2, real-to-complex)
void CpuFftRealForward(const double &in[], int len, double &out_r[], double &out_i[])
{
    if(len<=0 || (len & (len-1))!=0) return;
    ArrayResize(out_r,len); ArrayResize(out_i,len);
    for(int i=0;i<len;i++){ out_r[i]=in[i]; out_i[i]=0.0; }
    int j=0;
    for(int i=0;i<len;i++){
        if(i<j){ double tr=out_r[i]; out_r[i]=out_r[j]; out_r[j]=tr; double ti=out_i[i]; out_i[i]=out_i[j]; out_i[j]=ti; }
        int m=len>>1; while(m>=1 && j>=m){ j-=m; m>>=1; } j+=m;
    }
    for(int step=1; step<len; step<<=1){
        double theta=-M_PI/step; double wpr=-2.0*sin(0.5*theta)*sin(0.5*theta); double wpi=sin(theta);
        for(int m=0;m<step;m++){
            double wr=1.0, wi=0.0;
            for(int k=m;k<len;k+=(step<<1)){
                int l=k+step; double tr=wr*out_r[l]-wi*out_i[l]; double ti=wr*out_i[l]+wi*out_r[l];
                out_r[l]=out_r[k]-tr; out_i[l]=out_i[k]-ti; out_r[k]+=tr; out_i[k]+=ti;
            }
            double wrn=wr*wpr - wi*wpi + wr; double win=wi*wpr + wr*wpi + wi; wr=wrn; wi=win;
        }
    }
}

// PLA builder (fallback: flat segment)
bool BuildPlaPriceSeries(const int start_pos, const datetime &time[], const double &close[])
{
    static double feed_close_tf[]; ArrayResize(feed_close_tf, InpFFTWindow); ArraySetAsSeries(feed_close_tf, true);
    int shift = iBarShift(_Symbol, InpFeedTimeframe, time[start_pos]);
    if(shift<0) return false;
    if(CopyClose(_Symbol, InpFeedTimeframe, shift, InpFFTWindow, feed_close_tf)!=InpFFTWindow) return false;
    for(int j=0;j<InpFFTWindow;j++) feed_data[j]=feed_close_tf[InpFFTWindow-1-j];
    // trivial PLA: copy as-is
    return true;
}

// ZigZag price series builder (continuous or alternating)
bool BuildZigZagPriceSeries(const int start_pos,
                            const double &high[],
                            const double &low[],
                            const datetime &time[],
                            bool alternating)
{
    if(start_pos < 0) return false;
    int len = InpFFTWindow;
    for(int j=0; j<len; ++j)
    {
        int idx = start_pos + j;
        double v = (alternating ? ((j%2==0)? high[idx]: low[idx]) : (high[idx]+low[idx])*0.5);
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
    if(rates_total < InpFFTWindow) return prev_calculated;

    ArrayResize(feed_data, InpFFTWindow);
    ArrayResize(detrended_data, InpFFTWindow);
    ArrayResize(spectrum, InpFFTWindow/2);
    ArrayResize(WaveKalman, rates_total);

    int start = MathMax(prev_calculated-1, InpFFTWindow-1);
    for(int i=start;i<rates_total;i++)
    {
        int start_pos = i - InpFFTWindow + 1;
        if(start_pos<0) continue;

        // feed selection
        if(InpFeedData==FEED_CLOSE)
        {
            int shift = iBarShift(_Symbol, InpFeedTimeframe, time[start_pos]);
            if(shift<0) continue;
            static double buf[]; ArrayResize(buf, InpFFTWindow); ArraySetAsSeries(buf,true);
            if(CopyClose(_Symbol, InpFeedTimeframe, shift, InpFFTWindow, buf)!=InpFFTWindow) continue;
            for(int j=0;j<InpFFTWindow;j++) feed_data[j]=buf[InpFFTWindow-1-j];
        }
        else if(InpFeedData==FEED_PLA)
        {
            if(!BuildPlaPriceSeries(start_pos, time, close)) continue;
        }
        else // ZigZag continuous/alternating
        {
            bool alternating = (InpFeedData == FEED_ZIGZAG_ALTERNATING);
            if(!BuildZigZagPriceSeries(start_pos, high, low, time, alternating))
                continue;
        }

        ArrayCopy(detrended_data, feed_data, 0, 0, InpFFTWindow);

        // windowing: none (can add later)

        if(!EnsureGpu(InpFFTWindow))
            continue; // sem GPU, não processa

        int st = gpu_fft_real_forward(detrended_data, InpFFTWindow, g_fft_interleaved);
        if(st!=ALGLIB_STATUS_OK)
            continue;

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
