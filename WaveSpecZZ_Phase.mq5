//+------------------------------------------------------------------+
//| WaveSpecZZ Phase-only variant                                   |
//+------------------------------------------------------------------+
#property copyright "Copyright 2025"
#property link      ""
#property version   "2.000"
#property indicator_separate_window
#property indicator_buffers 1
#property indicator_plots   1
#property strict

#define WAVESPEC_VARIANT_PHASE

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
  int  mt_gpu_extract_cycles(const double &series[], int len, int top_k,
                             double min_period, double max_period, double seconds_per_sample,
                             double &out[], int out_stride, int out_capacity, int &out_len);
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
  int  tester_gpu_extract_cycles(const double &series[], int len, int top_k,
                                 double min_period, double max_period, double seconds_per_sample,
                                 double &out[], int out_stride, int out_capacity, int &out_len);
#import

bool UseTesterBridge()
{
    static int flag = -1;
    if(flag == -1)
        flag = (int)MQLInfoInteger(MQL_TESTER);
    return flag == 1;
}

int gpu_init_impl(int device_index, int stream_count)
{
    return UseTesterBridge() ? tester_gpu_init(device_index, stream_count)
                             : mt_gpu_init(device_index, stream_count);
}

void gpu_shutdown_impl()
{
    if(UseTesterBridge())
        tester_gpu_shutdown();
    else
        mt_gpu_shutdown();
}

int gpu_wave_submit_template_job_impl(const string &preset_text,
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

int gpu_wave_try_get_template_job_impl(int job_id,
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

int gpu_wave_free_template_job_impl(int job_id)
{
    return UseTesterBridge() ? tester_gpu_wave_free_template_job(job_id)
                             : mt_gpu_wave_free_template_job(job_id);
}

int gpu_wave_build_tick_series_impl(const double &tick_prices[], const long &tick_times[], int tick_count,
                                    int window_len, int interval_seconds, int smoothing_window,
                                    int zig_depth, double zig_deviation_points, int zig_backstep, int zig_mode,
                                    double point_value, double &out[], int out_len)
{
    return UseTesterBridge() ? tester_gpu_wave_build_tick_series(tick_prices,
                                                                 tick_times, tick_count,
                                                                 window_len, interval_seconds, smoothing_window,
                                                                 zig_depth, zig_deviation_points, zig_backstep, zig_mode,
                                                                 point_value, out, out_len)
                             : mt_gpu_wave_build_tick_series(tick_prices,
                                                           tick_times, tick_count,
                                                           window_len, interval_seconds, smoothing_window,
                                                           zig_depth, zig_deviation_points, zig_backstep, zig_mode,
                                                           point_value, out, out_len);
}

int gpu_extract_cycles_impl(const double &series[], int len, int top_k,
                            double min_period, double max_period, double seconds_per_sample,
                            double &out[], int out_stride, int out_capacity, int &out_len)
{
    return UseTesterBridge() ? tester_gpu_extract_cycles(series, len, top_k,
                                                         min_period, max_period, seconds_per_sample,
                                                         out, out_stride, out_capacity, out_len)
                             : mt_gpu_extract_cycles(series, len, top_k,
                                                     min_period, max_period, seconds_per_sample,
                                                     out, out_stride, out_capacity, out_len);
}

int gpu_get_last_error_w_impl(ushort &buf[], int buf_len)
{
    return UseTesterBridge() ? tester_gpu_get_last_error_w(buf, buf_len)
                             : mt_gpu_get_last_error_w(buf, buf_len);
}

#define gpu_init                      gpu_init_impl
#define gpu_shutdown                  gpu_shutdown_impl
#define gpu_wave_submit_template_job  gpu_wave_submit_template_job_impl
#define gpu_wave_try_get_template_job gpu_wave_try_get_template_job_impl
#define gpu_wave_free_template_job    gpu_wave_free_template_job_impl
#define gpu_wave_build_tick_series    gpu_wave_build_tick_series_impl
#define gpu_extract_cycles            gpu_extract_cycles_impl
#define gpu_get_last_error_w          gpu_get_last_error_w_impl

#include "WavePresetDsl.mqh"
#include "WaveSpecZZ_gpu_core.mqh"
