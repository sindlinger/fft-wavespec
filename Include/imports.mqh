#ifndef __WAVESPEC_IMPORTS_MQH__
#define __WAVESPEC_IMPORTS_MQH__

// Ajuste o nome da DLL aqui quando necessário (MQL não permite macro em #import)
#import "mt-bridge.dll"
  int  gpu_init(int device_index, int stream_count);
  void gpu_shutdown(void);
  int  gpu_fft_real_forward(const double &in[], int len, double &out[]);
  int  gpu_extract_cycles(const double &series[], int len, int top_k, double min_period, double max_period,
                          double sample_rate_seconds, int method, int ar_order,
                          double &out[], int out_stride, int out_capacity, int &out_len);
  int  gpu_submit_extract_cycles(const double &series[], int len, int top_k, double min_period, double max_period,
                                 double sample_rate_seconds, int method, int ar_order, long &job_id);
  int  gpu_try_get_cycles(long job_id, double &out[], int out_stride, int out_capacity, int &out_len, int &ready);
  int  gpu_submit_extract_cycles_batch(const double &series[], int series_len, int window_len, int hop, int top_k,
                                       double min_period, double max_period, double sample_rate_seconds,
                                       int method, int ar_order, int stride, long &job_id);
  int  gpu_try_get_cycles_batch(long job_id, double &out[], int out_cap, int &out_len, int &ready);
  int  gpu_free_job(long job_id);
  int  gpu_get_last_error_w(ushort &buf[], int buf_len);
#import

#endif // __WAVESPEC_IMPORTS_MQH__
