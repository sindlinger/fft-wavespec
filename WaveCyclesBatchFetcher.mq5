//+------------------------------------------------------------------+
//| WaveCyclesBatchFetcher                                           |
//| Gera e salva cache de ciclos (FFT/MUSIC) em batch para WaveSpec. |
//| Não plota nada; roda uma vez no attach.                          |
//+------------------------------------------------------------------+
#property indicator_chart_window
#property indicator_buffers 1
#property indicator_plots   0
#property strict

#include "Include/imports.mqh"

// Garantia: códigos de status (caso não venham dos headers por algum motivo)
#ifndef ALGLIB_STATUS_OK
  #define ALGLIB_STATUS_OK           0
  #define ALGLIB_STATUS_BAD_ARGS    -1
  #define ALGLIB_STATUS_BACKEND_UNAVAILABLE -2
  #define ALGLIB_STATUS_TIMEOUT     -3
  #define ALGLIB_STATUS_INTERNAL_ERROR -4
  #define ALGLIB_STATUS_NOT_READY   -5
  #define ALGLIB_STATUS_NO_MEM      -6
#endif

double __dummy[]; // buffer silencioso

input string InpSymbol        = "";   // vazio = _Symbol
input ENUM_TIMEFRAMES InpTF   = PERIOD_M1;
input int    InpFFTWindow     = 4096;
input int    InpTopK          = 4;     // 1..8
input double InpMinPeriod     = 9;
input double InpMaxPeriod     = 200;
input int    InpMethod        = 1;     // 1 = MUSIC/ESPRIT
input int    InpArOrder       = 10;
input int    InpHop           = 1;
input int    InpStride        = 15;    // layout compatível com gpu_extract_cycles
input int    InpBars          = 500000; // número máximo de barras a baixar (0=todas)

int OnInit()
{
   IndicatorSetString(INDICATOR_SHORTNAME, "WaveCyclesBatchFetcher");
   // executa imediatamente no attach
   EventSetTimer(1);
   return(INIT_SUCCEEDED);
}

void OnDeinit(const int)
{
   EventKillTimer();
}

// Nome do cache igual ao usado pelo WaveSpec
string CycleCacheName()
{
   string tf = EnumToString(InpTF);
    return StringFormat("WaveSpecZZ_cycles_%s_%s_w%d_m%d_ar%d_k%d.bin",
                        InpSymbol, tf, InpFFTWindow, InpMethod, InpArOrder, InpTopK);
}

void SaveCycleCache(const double &cycles[], int out_len, int stride, int bars)
{
   string file = CycleCacheName();
   int h = FileOpen(file, FILE_WRITE|FILE_BIN|FILE_COMMON);
   if(h==INVALID_HANDLE)
   {
      PrintFormat("[BatchFetcher][ERR] FileOpen failed %s", file);
      return;
   }
   FileWriteInteger(h, 1, INT_VALUE);    // version
   FileWriteInteger(h, bars, INT_VALUE); // bars
   FileWriteInteger(h, MathMin(InpTopK,2), INT_VALUE); // topk stored (2 waves)
   // Apenas duas waves (slots 0 e 1) com o stride completo (amplitude/period/eta/…)
   for(int i=0;i<out_len;i++)
   {
      int base=i*stride;
      FileWriteDouble(h, cycles[base+0]); // amp
      FileWriteDouble(h, cycles[base+1]); // freq
      FileWriteDouble(h, cycles[base+2]); // period
      FileWriteDouble(h, cycles[base+3]); // phase
      FileWriteDouble(h, cycles[base+5]); // eta_seconds
      FileWriteDouble(h, cycles[base+6]); // energy
      FileWriteDouble(h, cycles[base+7]); // coherence
      FileWriteDouble(h, cycles[base+8]); // snr_db
      FileWriteDouble(h, cycles[base+10]); // eigen_ratio
      FileWriteDouble(h, cycles[base+11]); // score
      FileWriteDouble(h, cycles[base+13]); // eta_conf
   }
   FileClose(h);
   PrintFormat("[BatchFetcher] cache salvo %s (bars=%d cycles=%d)", file, bars, out_len);
}

void OnTimer()
{
   EventKillTimer(); // roda uma vez só
   string sym = (InpSymbol=="" ? _Symbol : InpSymbol);

   // prepara série
   double prices[]; ArraySetAsSeries(prices,true);
   int total = Bars(sym, InpTF);
   if(total<=0){ Print("[BatchFetcher][ERR] sem barras"); return; }
   if(InpBars>0) total = MathMin(total, InpBars);
   int got = CopyClose(sym, InpTF, 0, total, prices);
   if(got < InpFFTWindow){ Print("[BatchFetcher][ERR] poucas barras"); return; }

   // garante backend
   int streams = 64;
   if(gpu_init(0, streams)!=ALGLIB_STATUS_OK)
   {
      Print("[BatchFetcher][ERR] gpu_init falhou");
      return;
   }

   int nwin = 1 + (got - InpFFTWindow)/InpHop;
   int buf_cap = nwin * InpTopK * InpStride;
   double cycles[]; ArrayResize(cycles, buf_cap);

   long jid=0; int st = gpu_submit_extract_cycles_batch(prices, got, InpFFTWindow, InpHop, InpTopK,
                                                        InpMinPeriod, InpMaxPeriod, (double)PeriodSeconds(InpTF),
                                                        InpMethod, InpArOrder, InpStride, jid);
   if(st!=ALGLIB_STATUS_OK || jid==0)
   {
      PrintFormat("[BatchFetcher][ERR] submit st=%d", st);
      return;
   }
   PrintFormat("[BatchFetcher] submitted job=%I64d len=%d win=%d hop=%d", jid, got, InpFFTWindow, InpHop);

   int ready=0,out_len=0;
   for(int tries=0; tries<4000 && ready==0; ++tries)
   {
      st = gpu_try_get_cycles_batch(jid, cycles, buf_cap, out_len, ready);
      if(st==ALGLIB_STATUS_OK && ready==0) Sleep(5);
      else if(st!=ALGLIB_STATUS_OK && st!=ALGLIB_STATUS_NOT_READY) break;
   }
   gpu_free_job(jid);

   if(st==ALGLIB_STATUS_OK && ready==1 && out_len>0)
   {
      SaveCycleCache(cycles, out_len, InpStride, got);
   }
   else
   {
      PrintFormat("[BatchFetcher][ERR] st=%d ready=%d out=%d", st, ready, out_len);
   }
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
   // nada a plotar
   return(rates_total);
}
