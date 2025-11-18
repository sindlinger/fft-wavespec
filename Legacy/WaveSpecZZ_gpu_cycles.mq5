// WaveSpecZZ_gpu_cycles.mq5
// Protótipo focado em extrair ciclos na GPU via gpu_extract_cycles / wave preset.
// Mantém lógica mínima: coleta janela, detrend/window na GPU, extrai top K ciclos e popula buffers.

#property indicator_separate_window
#property indicator_buffers 12
#property indicator_plots   12

#import "mt-bridge.dll"
  int  gpu_init(int device_index, int stream_count);
  void gpu_shutdown(void);
  int  gpu_fft_real_forward(const double &in[], int len, double &out[]);
  int  gpu_wave_fft_segmented(const double &in[], int len, int segment_len, int overlap, int mix_mode, double &out[], int out_len);
  int  gpu_fft_real_forward_batch(const double &in[], int window_len, int n_windows, double &out[]);
  int  gpu_extract_cycles(const double* series,
                         int           len,
                         int           top_k,
                         double        min_period,
                         double        max_period,
                         double        sample_rate_seconds,
                         double*       out,
                         int           out_stride,
                         int           out_capacity,
                         int*          out_len);
#import

input int  InpFFTWindow   = 8192;
input int  InpTopCycles   = 12;
input double InpMinPeriod = 4.0;
input double InpMaxPeriod = 64.0;
input bool InpUseSegmented = true;
input int  InpSegmentLen = 4096;
input int  InpSegmentOverlap = 1024;
input int  InpSegmentMixMode = 0; // 0=energy

double g_window[]; // preço 
double g_det[];    // detrended (aqui copiamos direto, sem detrend)
double g_fft[];
double g_cycles[12*4]; // stride 4: period, power, re, im (ou reservado)

double Wave1[],Wave2[],Wave3[],Wave4[],Wave5[],Wave6[],Wave7[],Wave8[],Wave9[],Wave10[],Wave11[],Wave12[];

int OnInit()
{
    if(gpu_init(0,2)!=0) return INIT_FAILED;
    SetIndexBuffer(0,Wave1,INDICATOR_DATA);
    SetIndexBuffer(1,Wave2,INDICATOR_DATA);
    SetIndexBuffer(2,Wave3,INDICATOR_DATA);
    SetIndexBuffer(3,Wave4,INDICATOR_DATA);
    SetIndexBuffer(4,Wave5,INDICATOR_DATA);
    SetIndexBuffer(5,Wave6,INDICATOR_DATA);
    SetIndexBuffer(6,Wave7,INDICATOR_DATA);
    SetIndexBuffer(7,Wave8,INDICATOR_DATA);
    SetIndexBuffer(8,Wave9,INDICATOR_DATA);
    SetIndexBuffer(9,Wave10,INDICATOR_DATA);
    SetIndexBuffer(10,Wave11,INDICATOR_DATA);
    SetIndexBuffer(11,Wave12,INDICATOR_DATA);
    IndicatorSetInteger(INDICATOR_DIGITS,_Digits);
    return INIT_SUCCEEDED;
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
    if(rates_total < InpFFTWindow)
        return prev_calculated;

    int start = MathMax(prev_calculated-1,0);
    for(int i=start; i<rates_total && !IsStopped(); ++i)
    {
        int start_pos = i - InpFFTWindow + 1;
        if(start_pos < 0) continue;

        ArrayResize(g_window, InpFFTWindow);
        ArrayCopy(g_window, close, 0, start_pos, InpFFTWindow);
        ArrayResize(g_det, InpFFTWindow);
        ArrayCopy(g_det, g_window, 0, 0, InpFFTWindow); // sem detrend por simplicidade

        // FFT (apenas para manter compatibilidade; ciclos virão do extract)
        ArrayResize(g_fft, InpFFTWindow);

        // EXTRAÇÃO DE CICLOS NA GPU
        ArrayInitialize(g_cycles, 0.0);
        int cycles_len=0;
        int status = gpu_extract_cycles(g_det,
                                        InpFFTWindow,
                                        InpTopCycles,
                                        InpMinPeriod,
                                        InpMaxPeriod,
                                        PeriodSeconds((ENUM_TIMEFRAMES)_Period),
                                        g_cycles,
                                        4,
                                        InpTopCycles,
                                        &cycles_len);
        if(status!=0 || cycles_len<=0)
            continue;

        // Preencher buffers: usa amplitude (power) como valor da wave
        for(int c=0;c<InpTopCycles && c<12; ++c)
        {
            int idx = c*4;
            double period = g_cycles[idx];
            double power  = g_cycles[idx+1];
            double val = power; // simples: usar power
            switch(c)
            {
                case 0: Wave1[i]=val; break; case 1: Wave2[i]=val; break; case 2: Wave3[i]=val; break;
                case 3: Wave4[i]=val; break; case 4: Wave5[i]=val; break; case 5: Wave6[i]=val; break;
                case 6: Wave7[i]=val; break; case 7: Wave8[i]=val; break; case 8: Wave9[i]=val; break;
                case 9: Wave10[i]=val; break; case 10: Wave11[i]=val; break; case 11: Wave12[i]=val; break;
            }
        }
    }
    return rates_total;
}

void OnDeinit(const int)
{
    gpu_shutdown();
}

