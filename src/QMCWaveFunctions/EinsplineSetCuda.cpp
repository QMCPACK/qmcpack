//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Christos Kartsaklis, kartsaklisc@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

// #define SPLIT_SPLINE_DEBUG
// #define USE_SPLIT_SPLINES_MEM_PREFETCH

#include "QMCWaveFunctions/EinsplineSet.h"
#include "einspline/multi_bspline.h"
#include "einspline/multi_bspline_eval_cuda.h"
#include "Configuration.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/AtomicOrbitalCuda.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/PhaseFactors.h"
#include "CPU/math.hpp"
#ifdef HAVE_MKL
#include <mkl_vml.h>
#endif

namespace qmcplusplus
{
inline void create_multi_UBspline_3d_cuda(multi_UBspline_3d_d* in, multi_UBspline_3d_s_cuda*& out)
{
  out = create_multi_UBspline_3d_s_cuda_conv(in);
}

inline void create_multi_UBspline_3d_cuda(multi_UBspline_3d_d* in, multi_UBspline_3d_d_cuda*& out)
{
  out = create_multi_UBspline_3d_d_cuda(in);
}

inline void create_multi_UBspline_3d_cuda(multi_UBspline_3d_z* in, multi_UBspline_3d_c_cuda*& out)
{
  out = create_multi_UBspline_3d_c_cuda_conv(in);
}

inline void create_multi_UBspline_3d_cuda(multi_UBspline_3d_z* in, multi_UBspline_3d_z_cuda*& out)
{
  out = create_multi_UBspline_3d_z_cuda(in);
}

inline void create_multi_UBspline_3d_cuda(multi_UBspline_3d_z* in, multi_UBspline_3d_d_cuda*& out)
{
  app_error() << "Attempted to convert complex CPU spline into a real "
              << " GPU spline.\n";
  abort();
}

inline void create_multi_UBspline_3d_cuda(multi_UBspline_3d_z* in, multi_UBspline_3d_s_cuda*& out)
{
  app_error() << "Attempted to convert complex CPU spline into a real "
              << " GPU spline.\n";
  abort();
}

// Real evaluation functions
inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline, TinyVector<double, 3> r, Vector<double>& psi)
{
  eval_multi_UBspline_3d_d(spline, r[0], r[1], r[2], psi.data());
}

inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline, TinyVector<double, 3> r, std::vector<double>& psi)
{
  eval_multi_UBspline_3d_d(spline, r[0], r[1], r[2], &(psi[0]));
}


inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline,
                               TinyVector<double, 3> r,
                               Vector<double>& psi,
                               Vector<TinyVector<double, 3>>& grad,
                               Vector<Tensor<double, 3>>& hess)
{
  eval_multi_UBspline_3d_d_vgh(spline, r[0], r[1], r[2], psi.data(), (double*)grad.data(), (double*)hess.data());
}

//////////////////////////////////
// Complex evaluation functions //
//////////////////////////////////
inline void EinsplineMultiEval(multi_UBspline_3d_z* restrict spline,
                               TinyVector<double, 3> r,
                               Vector<std::complex<double>>& psi)
{
  eval_multi_UBspline_3d_z(spline, r[0], r[1], r[2], psi.data());
}


inline void EinsplineMultiEval(multi_UBspline_3d_z* restrict spline,
                               TinyVector<double, 3> r,
                               Vector<std::complex<double>>& psi,
                               Vector<TinyVector<std::complex<double>, 3>>& grad,
                               Vector<Tensor<std::complex<double>, 3>>& hess)
{
  eval_multi_UBspline_3d_z_vgh(spline, r[0], r[1], r[2], psi.data(), (std::complex<double>*)grad.data(),
                               (std::complex<double>*)hess.data());
}


inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_s_cuda* spline,
                                              float* pos,
                                              float* sign,
                                              float* phi[],
                                              int N)
{
  eval_multi_multi_UBspline_3d_s_sign_cuda(spline, pos, sign, phi, N);
}

inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_s_cuda* spline,
                                              float* pos,
                                              float* sign,
                                              float* phi[],
                                              int N,
                                              float* spline_coefs[],
                                              cudaEvent_t spline_events[],
                                              cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
#ifdef SPLIT_SPLINE_DEBUG
    std::cerr << "Rank " << OHMMS::Controller->rank() << ", GPU mem before: " << spline->coefs
              << "; after: " << spline_coefs[devicenr] << "\n";
#endif
    eval_multi_multi_UBspline_3d_s_sign_cudasplit(spline, pos, sign, phi, N, spline_coefs[devicenr], devicenr,
                                                  spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}


inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_d_cuda* spline,
                                              double* pos,
                                              double* sign,
                                              double* phi[],
                                              int N)
{
  eval_multi_multi_UBspline_3d_d_sign_cuda(spline, pos, sign, phi, N);
}

inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_d_cuda* spline,
                                              double* pos,
                                              double* sign,
                                              double* phi[],
                                              int N,
                                              double* spline_coefs[],
                                              cudaEvent_t spline_events[],
                                              cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
#ifdef SPLIT_SPLINE_DEBUG
    std::cerr << "Rank " << OHMMS::Controller->rank() << ", GPU mem before: " << spline->coefs
              << "; after: " << spline_coefs[devicenr] << "\n";
#endif
    eval_multi_multi_UBspline_3d_d_sign_cudasplit(spline, pos, sign, phi, N, spline_coefs[devicenr], devicenr,
                                                  spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}


inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_d_cuda* spline, double* pos, double* phi[], int N)
{
  eval_multi_multi_UBspline_3d_d_cuda(spline, pos, phi, N);
}

inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_d_cuda* spline,
                                              double* pos,
                                              double* phi[],
                                              int N,
                                              double* spline_coefs[],
                                              cudaEvent_t spline_events[],
                                              cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
#ifdef SPLIT_SPLINE_DEBUG
    std::cerr << "Rank " << OHMMS::Controller->rank() << ", GPU mem before: " << spline->coefs
              << "; after: " << spline_coefs[devicenr] << "\n";
#endif
    eval_multi_multi_UBspline_3d_d_cudasplit(spline, pos, phi, N, spline_coefs[devicenr], devicenr,
                                             spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_s_cuda* spline,
                                                  float* pos,
                                                  float* sign,
                                                  float Linv[],
                                                  float* phi[],
                                                  float* grad_lapl[],
                                                  int N,
                                                  int row_stride)
{
  eval_multi_multi_UBspline_3d_s_vgl_sign_cuda(spline, pos, sign, Linv, phi, grad_lapl, N, row_stride);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_s_cuda* spline,
                                                  float* pos,
                                                  float* sign,
                                                  float Linv[],
                                                  float* phi[],
                                                  float* grad_lapl[],
                                                  int N,
                                                  int row_stride,
                                                  float* spline_coefs[],
                                                  cudaEvent_t spline_events[],
                                                  cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
    eval_multi_multi_UBspline_3d_s_vgl_sign_cudasplit(spline, pos, sign, Linv, phi, grad_lapl, N, row_stride,
                                                      spline_coefs[devicenr], devicenr, spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_vgl_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}


inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_d_cuda* spline,
                                                  double* pos,
                                                  double* sign,
                                                  double Linv[],
                                                  double* phi[],
                                                  double* grad_lapl[],
                                                  int N,
                                                  int row_stride)
{
  eval_multi_multi_UBspline_3d_d_vgl_sign_cuda(spline, pos, sign, Linv, phi, grad_lapl, N, row_stride);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_d_cuda* spline,
                                                  double* pos,
                                                  double* sign,
                                                  double Linv[],
                                                  double* phi[],
                                                  double* grad_lapl[],
                                                  int N,
                                                  int row_stride,
                                                  double* spline_coefs[],
                                                  cudaEvent_t spline_events[],
                                                  cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
    eval_multi_multi_UBspline_3d_d_vgl_sign_cudasplit(spline, pos, sign, Linv, phi, grad_lapl, N, row_stride,
                                                      spline_coefs[devicenr], devicenr, spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_vgl_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}


inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_d_cuda* spline,
                                                  double* pos,
                                                  double Linv[],
                                                  double* phi[],
                                                  double* grad_lapl[],
                                                  int N,
                                                  int row_stride)
{
  eval_multi_multi_UBspline_3d_d_vgl_cuda(spline, pos, Linv, phi, grad_lapl, N, row_stride);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_d_cuda* spline,
                                                  double* pos,
                                                  double Linv[],
                                                  double* phi[],
                                                  double* grad_lapl[],
                                                  int N,
                                                  int row_stride,
                                                  double* spline_coefs[],
                                                  cudaEvent_t spline_events[],
                                                  cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
    eval_multi_multi_UBspline_3d_d_vgl_cudasplit(spline, pos, Linv, phi, grad_lapl, N, row_stride,
                                                 spline_coefs[devicenr], devicenr, spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_vgl_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}

// Complex CUDA routines
inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_c_cuda* spline,
                                              float* pos,
                                              std::complex<float>* phi[],
                                              int N)
{
  eval_multi_multi_UBspline_3d_c_cuda(spline, pos, phi, N);
}


inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_c_cuda* spline,
                                              float* pos,
                                              std::complex<float>* phi[],
                                              int N,
                                              std::complex<float>* spline_coefs[],
                                              cudaEvent_t spline_events[],
                                              cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
#ifdef SPLIT_SPLINE_DEBUG
    std::cerr << "Rank " << OHMMS::Controller->rank() << ", GPU mem before: " << spline->coefs
              << "; after: " << spline_coefs[devicenr] << "\n";
#endif
    eval_multi_multi_UBspline_3d_c_cudasplit(spline, pos, phi, N, (float*)spline_coefs[devicenr], devicenr,
                                             spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}

inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_z_cuda* spline,
                                              double* pos,
                                              std::complex<double>* phi[],
                                              int N)
{
  eval_multi_multi_UBspline_3d_z_cuda(spline, pos, phi, N);
}

inline void eval_multi_multi_UBspline_3d_cuda(multi_UBspline_3d_z_cuda* spline,
                                              double* pos,
                                              std::complex<double>* phi[],
                                              int N,
                                              std::complex<double>* spline_coefs[],
                                              cudaEvent_t spline_events[],
                                              cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(double), curr_gpu, spline_streams[devicenr]);
#endif
#ifdef SPLIT_SPLINE_DEBUG
    std::cerr << "Rank " << OHMMS::Controller->rank() << ", GPU mem before: " << spline->coefs
              << "; after: " << spline_coefs[devicenr] << "\n";
#endif
    eval_multi_multi_UBspline_3d_z_cudasplit(spline, pos, phi, N, (double*)spline_coefs[devicenr], devicenr,
                                             spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_c_cuda* spline,
                                                  float* pos,
                                                  float Linv[],
                                                  std::complex<float>* phi[],
                                                  std::complex<float>* grad_lapl[],
                                                  int N,
                                                  int row_stride)
{
  eval_multi_multi_UBspline_3d_c_vgl_cuda(spline, pos, Linv, phi, grad_lapl, N, row_stride);
}


inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_c_cuda* spline,
                                                  float* pos,
                                                  float Linv[],
                                                  std::complex<float>* phi[],
                                                  std::complex<float>* grad_lapl[],
                                                  int N,
                                                  int row_stride,
                                                  std::complex<float>* spline_coefs[],
                                                  cudaEvent_t spline_events[],
                                                  cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(float), curr_gpu, spline_streams[devicenr]);
#endif
    eval_multi_multi_UBspline_3d_c_vgl_cudasplit(spline, pos, Linv, phi, grad_lapl, N, row_stride,
                                                 (float*)spline_coefs[devicenr], devicenr, spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_vgl_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_z_cuda* spline,
                                                  double* pos,
                                                  double Linv[],
                                                  std::complex<double>* phi[],
                                                  std::complex<double>* grad_lapl[],
                                                  int N,
                                                  int row_stride)
{
  eval_multi_multi_UBspline_3d_z_vgl_cuda(spline, pos, Linv, phi, grad_lapl, N, row_stride);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda(multi_UBspline_3d_z_cuda* spline,
                                                  double* pos,
                                                  double Linv[],
                                                  std::complex<double>* phi[],
                                                  std::complex<double>* grad_lapl[],
                                                  int N,
                                                  int row_stride,
                                                  std::complex<double>* spline_coefs[],
                                                  cudaEvent_t spline_events[],
                                                  cudaStream_t spline_streams[])
{
  int curr_gpu;
  int devicenr = 0;
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
  {
    devicenr = (i + gpu::relative_rank + 1) % gpu::device_group_size;
    curr_gpu = gpu::device_group_numbers[devicenr];
    cudaSetDevice(curr_gpu);
#ifdef USE_SPLIT_SPLINES_MEM_PREFETCH
    cudaMemPrefetchAsync(pos, 3 * N * sizeof(double), curr_gpu, spline_streams[devicenr]);
#endif
    eval_multi_multi_UBspline_3d_z_vgl_cudasplit(spline, pos, Linv, phi, grad_lapl, N, row_stride,
                                                 (double*)spline_coefs[devicenr], devicenr, spline_streams[devicenr]);
    cudaEventRecord(spline_events[devicenr], spline_streams[devicenr]);
  }
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in eval_multi_multi_UBspline_3d_vgl_cuda (rank %i):\n  %s\n", OHMMS::Controller->rank(),
            cudaGetErrorString(err));
    abort();
  }
  for (unsigned int i = 0; i < gpu::device_group_size; i++)
    cudaStreamWaitEvent(gpu::kernelStream, spline_events[i], 0);
}

//////////////////////////////////////////////
// Vectorized evaluation routines using GPU //
//////////////////////////////////////////////

#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<Walker_t*>& walkers,
                                            int iat,
                                            gpu::device_vector<CTS::RealType*>& phi)
{
  //  app_log() << "Start EinsplineSet CUDA evaluation\n";
  int N                       = walkers.size();
  CTS::RealType plus_minus[2] = {1.0, -1.0};
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    hostSign.resize(N);
    cudaSign.resize(N, 1.0, split_splines);
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = walkers[iw]->R[iat];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i = 0; i < OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int)img;
    }
    int sign = 0;
    for (int i = 0; i < OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    hostSign[iw] = plus_minus[sign & 1];
    hostPos[iw]  = ru;
  }
  cudapos  = hostPos;
  cudaSign = hostSign;
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)(cudapos.data()), cudaSign.data(), phi.data(), N,
                                      spline_rank_pointers.data(), spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)(cudapos.data()), cudaSign.data(), phi.data(),
                                      N);
  }
}

template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                          int iat,
                                                          gpu::device_vector<CTS::RealType*>& phi)
{
  // app_log() << "Eval 1.\n";
  int N = walkers.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    hostPhasePos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    cudaphasepos.resize(N);
    if (split_splines)
    {
      cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetReadMostly, 0);
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
        cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
    }
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = walkers[iw]->R[iat];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor(ru[0]);
    ru[1] -= std::floor(ru[1]);
    ru[2] -= std::floor(ru[2]);
    hostPos[iw]      = ru;
    hostPhasePos[iw] = r;
  }
  cudapos = hostPos;
  cudaphasepos.asyncCopy(hostPhasePos); // for phase factors kernel
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N,
                                      spline_rank_pointers.data(), spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N);
  }
  // Now, add on phases
  apply_phase_factors((CTS::RealType*)CudakPoints.data(), CudaMakeTwoCopies.data(), (CTS::RealType*)cudaphasepos.data(),
                      (CTS::RealType**)CudaValuePointers.data(), phi.data(), CudaMultiSpline->num_splines,
                      walkers.size());
}


#else
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<Walker_t*>& walkers,
                                            int iat,
                                            gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "Code should not arrive at this point: "
              << "EinsplineSetExtended<double>::evaluate at line " << __LINE__ << " in file " << __FILE__ << "\n";
  abort();
}

template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                          int iat,
                                                          gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetExtended<std::complex<double> >::evaluate at line " << __LINE__ << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
}
#endif


template<typename T>
void EinsplineSetExtended<T>::get_split_spline_pointers()
{
  split_splines = (gpu::device_group_size > 1);
  int size      = OHMMS::Controller->size(); // how many MPI ranks
#ifdef HAVE_MPI
  if ((size > 1) && split_splines)
  {
    app_log() << "Gathering einspline GPU memory pointers from all ranks.\n";
    std::vector<cudaIpcMemHandle_t> rank_handle(1);
    cudaIpcGetMemHandle(rank_handle.data(), CudaMultiSpline->coefs);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error setting memory handle:\n  %s\n", cudaGetErrorString(err));
      abort();
    }
    if (spline_rank_handles.size() < size)
      spline_rank_handles.resize(size);

    MPI_Allgather(&rank_handle[0], sizeof(cudaIpcMemHandle_t), MPI_CHAR, &spline_rank_handles[0],
                  sizeof(cudaIpcMemHandle_t), MPI_CHAR, OHMMS::Controller->getMPI());

    if ((spline_rank_pointers.size() < gpu::device_group_size) || (spline_events.size() < gpu::device_group_size) ||
        (spline_streams.size() < gpu::device_group_size))
    {
      spline_rank_pointers.resize(gpu::device_group_size);
      spline_events.resize(gpu::device_group_size);
      spline_streams.resize(gpu::device_group_size);
    }
    for (unsigned int i = 0; i < gpu::device_group_size; i++)
    {
      cudaSetDevice(gpu::device_group_numbers[i]);
      if (i != gpu::relative_rank % gpu::device_group_size)
      {
        cudaIpcOpenMemHandle((void**)&spline_rank_pointers[i], spline_rank_handles[gpu::device_rank_numbers[i]],
                             cudaIpcMemLazyEnablePeerAccess);
        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "CUDA error getting memory handles between GPUs %i and %i:\n  %s\n",
                  gpu::device_group_numbers[gpu::relative_rank % gpu::device_group_size], gpu::device_group_numbers[i],
                  cudaGetErrorString(err));
          abort();
        }
      }
      else
        spline_rank_pointers[i] = CudaMultiSpline->coefs;
      cudaStreamCreate(&spline_streams[i]);
      cudaEventCreateWithFlags(&spline_events[i], cudaEventDisableTiming);
    }
    cudaSetDevice(gpu::device_group_numbers[gpu::relative_rank % gpu::device_group_size]);
    app_log() << "Successful allgather.\n";
#ifdef SPLIT_SPLINE_DEBUG
    std::cerr << "Rank " << OHMMS::Controller->rank() << " pointers: ";
    for (unsigned int i = 0; i < gpu::device_group_size; i++)
      std::cerr << spline_rank_pointers[i] << " ";
    std::cerr << "\n";
#endif
  }
#endif
}

template<typename T>
void EinsplineSetExtended<T>::resize_cuda(int numWalkers)
{
  CudaValuePointers.resize(numWalkers, 1.0, split_splines); // use managed memory with split splines
  CudaGradLaplPointers.resize(numWalkers, 1.0, split_splines);
  int N      = CudaMultiSpline->num_splines;
  int Nsplit = CudaMultiSpline->num_split_splines;
  if (split_splines)
  {
    cudaMemAdvise(CudaValuePointers.data(), numWalkers * sizeof(CudaStorageType*), cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(CudaGradLaplPointers.data(), numWalkers * sizeof(CudaStorageType*), cudaMemAdviseSetReadMostly, 0);
    N = Nsplit * gpu::device_group_size;
  }
  CudaValueVector.resize(N * numWalkers, 1.0, split_splines);
  CudaGradLaplVector.resize(4 * N * numWalkers, 1.0, split_splines);
  if (split_splines)
    for (unsigned int i = 0; i < gpu::device_group_size; i++)
    {
      cudaMemAdvise(CudaValuePointers.data(), numWalkers * sizeof(CudaStorageType*), cudaMemAdviseSetAccessedBy,
                    gpu::device_group_numbers[i]);
      cudaMemAdvise(CudaGradLaplPointers.data(), numWalkers * sizeof(CudaStorageType*), cudaMemAdviseSetAccessedBy,
                    gpu::device_group_numbers[i]);
      if (i == gpu::relative_rank)
      { // all of the output memory is accessed by this rank's GPU
        cudaMemAdvise(CudaValueVector.data(), N * numWalkers * sizeof(CudaStorageType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
        cudaMemAdvise(CudaGradLaplVector.data(), 4 * N * numWalkers * sizeof(CudaStorageType),
                      cudaMemAdviseSetAccessedBy, gpu::device_group_numbers[i]);
      }
      else
      { // only a subsection of the output memory is written to by each other GPU
        cudaMemAdvise(&CudaValueVector.data()[i * numWalkers * Nsplit], Nsplit * numWalkers * sizeof(CudaStorageType),
                      cudaMemAdviseSetAccessedBy, gpu::device_group_numbers[i]);
        cudaMemAdvise(&CudaGradLaplVector.data()[4 * i * numWalkers * Nsplit],
                      4 * Nsplit * numWalkers * sizeof(CudaStorageType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
      }
    }
  gpu::host_vector<CudaStorageType*> hostValuePointers(numWalkers);
  gpu::host_vector<CudaStorageType*> hostGradLaplPointers(numWalkers);
  for (int i = 0; i < numWalkers; i++)
  {
    hostValuePointers[i]    = &(CudaValueVector.data()[i * Nsplit]);
    hostGradLaplPointers[i] = &(CudaGradLaplVector.data()[4 * i * Nsplit]);
  }
  CudaValuePointers.asyncCopy(hostValuePointers);
  CudaGradLaplPointers.asyncCopy(hostGradLaplPointers);
  int M = MakeTwoCopies.size();
  CudaMakeTwoCopies.resize(M);
  CudaTwoCopiesIndex.resize(M);
  gpu::host_vector<int> hostMakeTwoCopies(M);
  gpu::host_vector<int> hostTwoCopiesIndex(M);
  int TwoCopiesIndexCounter = 0;
  for (int i = 0; i < M; i++)
  {
    hostMakeTwoCopies[i]  = MakeTwoCopies[i];
    hostTwoCopiesIndex[i] = TwoCopiesIndexCounter;
    TwoCopiesIndexCounter = MakeTwoCopies[i] ? TwoCopiesIndexCounter + 2 : TwoCopiesIndexCounter + 1;
  }
  CudaMakeTwoCopies.asyncCopy(hostMakeTwoCopies);
  CudaTwoCopiesIndex.asyncCopy(hostTwoCopiesIndex);
  CudakPoints.resize(M);
  CudakPoints_reduced.resize(M);
  gpu::host_vector<TinyVector<CUDA_PRECISION, OHMMS_DIM>> hostkPoints(M), hostkPoints_reduced(M);
  for (int i = 0; i < M; i++)
  {
    //      PosType k_red1 = PrimLattice.toCart(kPoints[i]);
    PosType k_red2(dot(kPoints[i], PrimLattice.a(0)), dot(kPoints[i], PrimLattice.a(1)),
                   dot(kPoints[i], PrimLattice.a(2)));
    //       fprintf (stderr, "kred1 = %8.3f %8.3f %8.3f\n", k_red1[0], k_red1[1], k_red1[2]);
    //       fprintf (stderr, "kred2 = %8.3f %8.3f %8.3f\n", k_red2[0], k_red2[1], k_red2[2]);
    for (int j = 0; j < OHMMS_DIM; j++)
    {
      hostkPoints[i][j]         = kPoints[i][j];
      hostkPoints_reduced[i][j] = k_red2[j];
    }
  }
  CudakPoints.asyncCopy(hostkPoints);
  CudakPoints_reduced =
      hostkPoints_reduced; // make sure the last copy is synchronous so when we return from this function everything's finished
}

#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<Walker_t*>& walkers,
                                            std::vector<PosType>& newpos,
                                            gpu::device_vector<CTS::RealType*>& phi)
{
  //  app_log() << "Start EinsplineSet CUDA evaluation\n";
  int N                       = newpos.size();
  CTS::RealType plus_minus[2] = {1.0, -1.0};
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    hostSign.resize(N);
    cudaSign.resize(N, 1.0, split_splines);
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i = 0; i < OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int)img;
    }
    int sign = 0;
    for (int i = 0; i < OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    hostSign[iw] = plus_minus[sign & 1];
    hostPos[iw]  = ru;
  }
  cudapos  = hostPos;
  cudaSign = hostSign;
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)(cudapos.data()), cudaSign.data(), phi.data(), N,
                                      spline_rank_pointers.data(), spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)(cudapos.data()), cudaSign.data(), phi.data(),
                                      N);
  }
  //app_log() << "End EinsplineSet CUDA evaluation\n";
}


template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                          std::vector<PosType>& newpos,
                                                          gpu::device_vector<CTS::RealType*>& phi)
{
  // app_log() << "Eval 2.\n";
  int N = walkers.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    hostPhasePos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    cudaphasepos.resize(N);
    if (split_splines)
    {
      cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetReadMostly, 0);
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
        cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
    }
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor(ru[0]);
    ru[1] -= std::floor(ru[1]);
    ru[2] -= std::floor(ru[2]);
    hostPos[iw]      = ru;
    hostPhasePos[iw] = r;
  }
  cudapos = hostPos;
  cudaphasepos.asyncCopy(hostPhasePos); // for phase factors kernel
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N,
                                      spline_rank_pointers.data(), spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N);
  }
  // Now, add on phases
  apply_phase_factors((CTS::RealType*)CudakPoints.data(), CudaMakeTwoCopies.data(), (CTS::RealType*)cudaphasepos.data(),
                      (CTS::RealType**)CudaValuePointers.data(), phi.data(), CudaMultiSpline->num_splines,
                      walkers.size());
}

#else
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<Walker_t*>& walkers,
                                            std::vector<PosType>& newpos,
                                            gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "Code should not arrive at this point: "
              << "EinsplineSetExtended<double>::evaluate at line " << __LINE__ << " in file " << __FILE__ << "\n";
  abort();
}

template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                          std::vector<PosType>& newpos,
                                                          gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetExtended<std::complex<double>>::evaluate at line " << __LINE__ << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
}
#endif


#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<Walker_t*>& walkers,
                                            std::vector<PosType>& newpos,
                                            gpu::device_vector<CTS::RealType*>& phi,
                                            gpu::device_vector<CTS::RealType*>& grad_lapl,
                                            int row_stride,
                                            int k,
                                            bool klinear)
{
  int nw     = walkers.size();
  int N      = newpos.size();
  int offset = 0;
  if ((nw != N) && klinear)
  {
    offset = k * nw;
    N      = nw;
  }
  CTS::RealType plus_minus[2] = {1.0, -1.0};
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    hostSign.resize(N);
    cudaSign.resize(N, 1.0, split_splines);
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = newpos[iw + offset];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i = 0; i < OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int)img;
    }
    int sign = 0;
    for (int i = 0; i < OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    hostSign[iw] = plus_minus[sign & 1];
    hostPos[iw]  = ru;
  }
  cudapos  = hostPos;
  cudaSign = hostSign;
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_vgl_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), cudaSign.data(),
                                          Linv_cuda.data(), &(phi.data()[offset]), &(grad_lapl.data()[offset]), N,
                                          row_stride, spline_rank_pointers.data(), spline_events.data(),
                                          spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_vgl_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), cudaSign.data(),
                                          Linv_cuda.data(), &(phi.data()[offset]), &(grad_lapl.data()[offset]), N,
                                          row_stride);
  }
  //gpu::host_vector<CTS::RealType*> pointers;
  //pointers = phi;
  //CTS::RealType data[N];
  //cudaMemcpy (data, pointers[0], N*sizeof(CTS::RealType), cudaMemcpyDeviceToHost);
  //for (int i=0; i<N; i++)
  //  fprintf (stderr, "%1.12e\n", data[i]);
}


template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                          std::vector<PosType>& newpos,
                                                          gpu::device_vector<CTS::RealType*>& phi,
                                                          gpu::device_vector<CTS::RealType*>& grad_lapl,
                                                          int row_stride,
                                                          int k,
                                                          bool klinear)
{
  // app_log() << "Eval 3.\n";
  int nw     = walkers.size();
  int N      = newpos.size(); // should work out-of-the-box for k-delayed positions now
  int offset = 0;
  if ((nw != N) && klinear)
  {
    offset = k * nw;
    N      = nw;
  }
  int M = CudaMultiSpline->num_splines;
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    hostPhasePos.resize(N);
    cudapos.resize(N, 1.0, split_splines); // use managed memory here for split splines
    cudaphasepos.resize(N);
    if (split_splines)
    {
      cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetReadMostly, 0);
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
        cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
    }
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = newpos[iw + offset];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor(ru[0]);
    ru[1] -= std::floor(ru[1]);
    ru[2] -= std::floor(ru[2]);
    hostPos[iw]      = ru;
    hostPhasePos[iw] = r;
  }
  cudapos = hostPos;
  cudaphasepos.asyncCopy(hostPhasePos); // for phase factors kernel
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_vgl_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), Linv_cuda.data(),
                                          CudaValuePointers.data(), CudaGradLaplPointers.data(), N,
                                          CudaMultiSpline->num_split_splines, spline_rank_pointers.data(),
                                          spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_vgl_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), Linv_cuda.data(),
                                          CudaValuePointers.data(), CudaGradLaplPointers.data(), N,
                                          CudaMultiSpline->num_splines);
  }
  // DEBUG
  //  TinyVector<double,OHMMS_DIM> r(hostPos[0][0], hostPos[0][1], hostPos[0][2]);
  //  Vector<std::complex<double > > psi(M);
  //  Vector<TinyVector<std::complex<double>,3> > grad(M);
  //  Vector<Tensor<std::complex<double>,3> > hess(M);
  //  EinsplineMultiEval (MultiSpline, r, psi, grad, hess);
  //  // std::complex<double> cpuSpline[M];
  //  // TinyVector<std::complex<double>,OHMMS_DIM> std::complex<double> cpuGrad[M];
  //  // Tensor cpuHess[M];
  //  // eval_multi_UBspline_3d_z_vgh (MultiSpline, hostPos[0][0], hostPos[0][1], hostPos[0][2],
  //  // 				  cpuSpline);
  //  gpu::host_vector<CudaStorageType*> pointers;
  //  pointers = CudaGradLaplPointers;
  //  std::complex<CTS::RealType> gpuSpline[4*M];
  //  cudaMemcpy(gpuSpline, pointers[10], 4*M * sizeof(std::complex<CTS::RealType>), cudaMemcpyDeviceToHost);
  //  for (int i=0; i<M; i++)
  //    fprintf (stderr, "real: %10.6e %10.6e %10.6e , imag: %10.6e %10.6e %10.6e .\n",
  //             trace(hess[i],GGt).real(), gpuSpline[3*M+i].real(), trace(hess[i],GGt).real() - gpuSpline[3*M+i].real(),
  //             trace(hess[i], GGt).imag(), gpuSpline[3*M+i].imag(), trace(hess[i], GGt).imag() - gpuSpline[3*M+i].imag());
  //  fprintf (stderr, "\n");

// AT debug:
#ifdef SPLIT_SPLINE_DEBUG
  if (gpu::rank == 1)
  {
    if (split_splines)
    {
      int mygpu = gpu::device_group_numbers[gpu::relative_rank % gpu::device_group_size];
      int curr_gpu;
      int devicenr = 0;
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
      {
        devicenr = (i + gpu::relative_rank) % gpu::device_group_size;
        curr_gpu = gpu::device_group_numbers[devicenr];
        cudaSetDevice(curr_gpu);
        cudaDeviceSynchronize();
      }
      cudaSetDevice(mygpu); // set device back to original GPU for this rank
    }
    else
      cudaDeviceSynchronize();
    for (unsigned int iw = 0; iw < N; iw++)
      for (unsigned int g = 0; g < gpu::device_group_size; g++)
      {
        int off = g * CudaMultiSpline->num_split_splines * N;
        for (unsigned int j = 0; j < CudaMultiSpline->num_split_splines; j++)
          fprintf(stderr, "walker %i, orbital %i: %f+%fi, (%f+%fi, %f+%fi, %f+%fi) | %f+%fi\n", iw,
                  j + g * CudaMultiSpline->num_split_splines, CudaValuePointers[iw][j + off].real(),
                  CudaValuePointers[iw][j + off].imag(),
                  CudaGradLaplPointers[iw][0 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][0 * CudaMultiSpline->num_split_splines + j + 4 * off].imag(),
                  CudaGradLaplPointers[iw][1 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][1 * CudaMultiSpline->num_split_splines + j + 4 * off].imag(),
                  CudaGradLaplPointers[iw][2 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][2 * CudaMultiSpline->num_split_splines + j + 4 * off].imag(),
                  CudaGradLaplPointers[iw][3 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][3 * CudaMultiSpline->num_split_splines + j + 4 * off].imag());
      }
  }
#endif

  // Now, add on phases
  /* Original implementation
  apply_phase_factors ((CTS::RealType*) CudakPoints.data(),
                       CudaMakeTwoCopies.data(),
                       (CTS::RealType*)cudapos.data(),
                       (CTS::RealType**)CudaValuePointers.data(), phi.data(),
                       (CTS::RealType**)CudaGradLaplPointers.data(), grad_lapl.data(),
                       CudaMultiSpline->num_splines,  walkers.size(), row_stride);
  */
  // Ye: optimized memory access.
  apply_phase_factors((CTS::RealType*)CudakPoints.data(), CudaMakeTwoCopies.data(), CudaTwoCopiesIndex.data(),
                      (CTS::RealType*)cudaphasepos.data(), (CTS::RealType**)CudaValuePointers.data(),
                      &(phi.data()[offset]), (CTS::RealType**)CudaGradLaplPointers.data(), &(grad_lapl.data()[offset]),
                      CudaMultiSpline->num_splines, N, row_stride);
// AT debug:
#ifdef SPLIT_SPLINE_DEBUG
  if (abort_counter >= 3)
  {
    _exit(0);
  }
  else
    abort_counter++;
#endif
}


#else
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<Walker_t*>& walkers,
                                            std::vector<PosType>& newpos,
                                            gpu::device_vector<CTS::ComplexType*>& phi,
                                            gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                                            int row_stride,
                                            int k,
                                            bool klinear)
{
  app_error() << "Code should not arrive at this point: "
              << "EinsplineSetExtended<double>::evaluate at line " << __LINE__ << " in file " << __FILE__ << "\n";
  abort();
}

template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                          std::vector<PosType>& newpos,
                                                          gpu::device_vector<CTS::ComplexType*>& phi,
                                                          gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                                                          int row_stride,
                                                          int k,
                                                          bool klinear)
{
  int nw     = walkers.size();
  int N      = newpos.size(); // should work out-of-the-box for k-delayed positions now
  int offset = 0;
  if ((nw != N) && klinear)
  {
    offset = k * nw;
    N      = nw;
  }
  int M = CudaMultiSpline->num_splines;
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    hostPhasePos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    cudaphasepos.resize(N);
    if (split_splines)
    {
      cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetReadMostly, 0);
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
        cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
    }
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = newpos[iw + offset];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor(ru[0]);
    ru[1] -= std::floor(ru[1]);
    ru[2] -= std::floor(ru[2]);
    hostPos[iw]      = ru;
    hostPhasePos[iw] = r;
  }
  cudapos = hostPos;
  cudaphasepos.asyncCopy(hostPhasePos); // for phase factors kernel
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_vgl_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), Linv_cuda.data(),
                                          CudaValuePointers.data(), CudaGradLaplPointers.data(), N,
                                          CudaMultiSpline->num_split_splines, spline_rank_pointers.data(),
                                          spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_vgl_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), Linv_cuda.data(),
                                          CudaValuePointers.data(), CudaGradLaplPointers.data(), N,
                                          CudaMultiSpline->num_splines);
  }
// AT debug:
#ifdef SPLIT_SPLINE_DEBUG
  if (gpu::rank == 1)
  {
    if (split_splines)
    {
      int mygpu = gpu::device_group_numbers[gpu::relative_rank % gpu::device_group_size];
      int curr_gpu;
      int devicenr = 0;
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
      {
        devicenr = (i + gpu::relative_rank) % gpu::device_group_size;
        curr_gpu = gpu::device_group_numbers[devicenr];
        cudaSetDevice(curr_gpu);
        cudaDeviceSynchronize();
      }
      cudaSetDevice(mygpu); // set device back to original GPU for this rank
    }
    else
      cudaDeviceSynchronize();
    for (unsigned int iw = 0; iw < N; iw++)
      for (unsigned int g = 0; g < gpu::device_group_size; g++)
      {
        int off = g * CudaMultiSpline->num_split_splines * N;
        for (unsigned int j = 0; j < CudaMultiSpline->num_split_splines; j++)
          fprintf(stderr, "walker %i, orbital %i: %f+%fi, (%f+%fi, %f+%fi, %f+%fi) | %f+%fi\n", iw,
                  j + g * CudaMultiSpline->num_split_splines, CudaValuePointers[iw][j + off].real(),
                  CudaValuePointers[iw][j + off].imag(),
                  CudaGradLaplPointers[iw][0 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][0 * CudaMultiSpline->num_split_splines + j + 4 * off].imag(),
                  CudaGradLaplPointers[iw][1 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][1 * CudaMultiSpline->num_split_splines + j + 4 * off].imag(),
                  CudaGradLaplPointers[iw][2 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][2 * CudaMultiSpline->num_split_splines + j + 4 * off].imag(),
                  CudaGradLaplPointers[iw][3 * CudaMultiSpline->num_split_splines + j + 4 * off].real(),
                  CudaGradLaplPointers[iw][3 * CudaMultiSpline->num_split_splines + j + 4 * off].imag());
      }
  }
#endif
  // Now, add on phases
  apply_phase_factors((CTS::RealType*)CudakPoints.data(), (CTS::RealType*)cudaphasepos.data(),
                      (CTS::ValueType**)CudaValuePointers.data(), (CTS::ValueType**)&(phi.data()[offset]),
                      (CTS::ValueType**)CudaGradLaplPointers.data(), (CTS::ValueType**)&(grad_lapl.data()[offset]),
                      CudaMultiSpline->num_splines, N, row_stride);
// AT debug:
/*  if((gpu::rank==1) && (abort_counter%768==0))
  {
    gpu::host_vector<CTS::ValueType*> pointers;
    pointers = CudaGradLaplPointers;
    CTS::ValueType data[4*N*M], data_new[4*N*M], phi_orig[N*M], phi_new[N*M];
    cudaMemcpy (data, pointers[0], 4*N*M*sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
    pointers = grad_lapl;
    cudaMemcpy (data_new, pointers[0], 4*N*M*sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
    pointers = CudaValuePointers;
    cudaMemcpy (phi_orig, pointers[0], N*M*sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
    pointers = phi;
    cudaMemcpy (phi_new, pointers[0], N*M*sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
    std::cerr << "CudaGradLaplPointers -> grad_lapl (# splines: " << M << "):\n";
    for (int iw=0; iw<N; iw++)
      for(unsigned int g=0; g<gpu::device_group_size; g++)
      {
        int off = g*CudaMultiSpline->num_split_splines*N;
        for(unsigned int j=0; j<CudaMultiSpline->num_split_splines; j++)
          fprintf(stderr,"walker %i, orbital %i: %f+%fi -> %f+%fi (%f+%fi, %f+%fi, %f+%fi) | %f+%fi -> (%f+%fi, %f+%fi, %f+%fi) | %f+%fi\n",iw,j+g*CudaMultiSpline->num_split_splines,
                        phi_orig[iw*CudaMultiSpline->num_split_splines+off+j].real(),phi_orig[iw*CudaMultiSpline->num_split_splines+off+j].imag(),
                        phi_new[g*CudaMultiSpline->num_split_splines+j].real(),phi_new[g*CudaMultiSpline->num_split_splines+j].imag(),
                        data[iw*4*CudaMultiSpline->num_split_splines+0*CudaMultiSpline->num_split_splines+j+4*off].real(),data[iw*4*CudaMultiSpline->num_split_splines+0*CudaMultiSpline->num_split_splines+j+4*off].imag(),
                        data[iw*4*CudaMultiSpline->num_split_splines+1*CudaMultiSpline->num_split_splines+j+4*off].real(),data[iw*4*CudaMultiSpline->num_split_splines+1*CudaMultiSpline->num_split_splines+j+4*off].imag(),
                        data[iw*4*CudaMultiSpline->num_split_splines+2*CudaMultiSpline->num_split_splines+j+4*off].real(),data[iw*4*CudaMultiSpline->num_split_splines+2*CudaMultiSpline->num_split_splines+j+4*off].imag(),
                        data[iw*4*CudaMultiSpline->num_split_splines+3*CudaMultiSpline->num_split_splines+j+4*off].real(),data[iw*4*CudaMultiSpline->num_split_splines+3*CudaMultiSpline->num_split_splines+j+4*off].imag(),
                        data_new[iw*4*M+0*M+j+g*CudaMultiSpline->num_split_splines].real(),data_new[iw*4*M+0*M+j+g*CudaMultiSpline->num_split_splines].imag(),
                        data_new[iw*4*M+1*M+j+g*CudaMultiSpline->num_split_splines].real(),data_new[iw*4*M+1*M+j+g*CudaMultiSpline->num_split_splines].imag(),
                        data_new[iw*4*M+2*M+j+g*CudaMultiSpline->num_split_splines].real(),data_new[iw*4*M+2*M+j+g*CudaMultiSpline->num_split_splines].imag(),
                        data_new[iw*4*M+3*M+j+g*CudaMultiSpline->num_split_splines].real(),data_new[iw*4*M+3*M+j+g*CudaMultiSpline->num_split_splines].imag());
      }
  }
  abort_counter++;*/
#ifdef SPLIT_SPLINE_DEBUG
  if (abort_counter >= 3)
  {
    _exit(0);
  }
  else
    abort_counter++;
#endif
}
#endif


#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetExtended<double>::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi)
{
  int N                       = pos.size();
  CTS::RealType plus_minus[2] = {1.0, -1.0};
  if (NLcudapos.size() < N)
  {
    NLhostPos.resize(N);
    NLcudapos.resize(N, 1.0, split_splines);
    NLhostSign.resize(N);
    NLcudaSign.resize(N, 1.0, split_splines);
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = pos[iw];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i = 0; i < OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int)img;
    }
    int sign = 0;
    for (int i = 0; i < OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    NLhostSign[iw] = plus_minus[sign & 1];
    NLhostPos[iw]  = ru;
  }
  NLcudapos  = NLhostPos;
  NLcudaSign = NLhostSign;
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)(NLcudapos.data()), NLcudaSign.data(),
                                      phi.data(), N, spline_rank_pointers.data(), spline_events.data(),
                                      spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)(NLcudapos.data()), NLcudaSign.data(),
                                      phi.data(), N);
  }
}

template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<PosType>& pos,
                                                          gpu::device_vector<CTS::RealType*>& phi)
{
  // app_log() << "Eval 4.\n";
  int N = pos.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    hostPhasePos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    cudaphasepos.resize(N);
    if (split_splines)
    {
      cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetReadMostly, 0);
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
        cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
    }
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = pos[iw];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor(ru[0]);
    ru[1] -= std::floor(ru[1]);
    ru[2] -= std::floor(ru[2]);
    hostPos[iw]      = ru;
    hostPhasePos[iw] = r;
  }
  cudapos = hostPos;
  cudaphasepos.asyncCopy(hostPhasePos); // for phase factors kernel
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N,
                                      spline_rank_pointers.data(), spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N);
  }
  // Now, add on phases
  apply_phase_factors((CTS::RealType*)CudakPoints.data(), CudaMakeTwoCopies.data(), (CTS::RealType*)cudaphasepos.data(),
                      (CTS::RealType**)CudaValuePointers.data(), phi.data(), CudaMultiSpline->num_splines, N);
}

#else

template<>
void EinsplineSetExtended<double>::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetExtended<std::complex<double> >::evaluate at " << __LINE__ << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
}

template<>
void EinsplineSetExtended<std::complex<double>>::evaluate(std::vector<PosType>& pos,
                                                          gpu::device_vector<CTS::ComplexType*>& phi)
{
  int N = pos.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    hostPhasePos.resize(N);
    cudapos.resize(N, 1.0, split_splines);
    cudaphasepos.resize(N);
    if (split_splines)
    {
      cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetReadMostly, 0);
      for (unsigned int i = 0; i < gpu::device_group_size; i++)
        cudaMemAdvise(cudapos.data(), N * sizeof(CTS::PosType), cudaMemAdviseSetAccessedBy,
                      gpu::device_group_numbers[i]);
    }
  }
  for (int iw = 0; iw < N; iw++)
  {
    PosType r = pos[iw];
    PosType ru(PrimLattice.toUnit(r));
    for (int i = 0; i < OHMMS_DIM; i++)
      ru[i] -= std::floor(ru[i]);
    hostPos[iw]      = ru;
    hostPhasePos[iw] = r;
  }
  cudapos = hostPos;
  cudaphasepos.asyncCopy(hostPhasePos); // for phase factors kernel
  if (split_splines)
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N,
                                      spline_rank_pointers.data(), spline_events.data(), spline_streams.data());
  }
  else
  {
    eval_multi_multi_UBspline_3d_cuda(CudaMultiSpline, (CTS::RealType*)cudapos.data(), CudaValuePointers.data(), N);
  }
  // Now, add on phases
  apply_phase_factors((CTS::RealType*)CudakPoints.data(), (CTS::RealType*)cudaphasepos.data(), CudaValuePointers.data(),
                      phi.data(), CudaMultiSpline->num_splines, N);
  // AT debug:
  /*  if(gpu::rank==1)
  {
    int M = CudaMultiSpline->num_splines;
    CTS::ValueType phi_orig[N*M], phi_new[N*M];
    gpu::host_vector<CTS::ValueType*> pointers;
    pointers = CudaValuePointers;
    cudaMemcpy (phi_orig, pointers[0], N*M*sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
    pointers = phi;
    cudaMemcpy (phi_new, pointers[0], N*M*sizeof(CTS::ValueType), cudaMemcpyDeviceToHost);
    std::cerr << "CudaValues -> phi (# splines: " << M << "):\n";
    for (int iw=0; iw<N; iw++)
      for(unsigned int g=0; g<gpu::device_group_size; g++)
      {
        int off = g*CudaMultiSpline->num_split_splines*N;
        for(unsigned int j=0; j<CudaMultiSpline->num_split_splines; j++)
          fprintf(stderr,"walker %i, orbital %i: %f+%fi -> %f+%fi\n",iw,j+g*CudaMultiSpline->num_split_splines,
                        phi_orig[iw*CudaMultiSpline->num_split_splines+off+j].real(),phi_orig[iw*CudaMultiSpline->num_split_splines+off+j].imag(),
                        phi_new[g*CudaMultiSpline->num_split_splines+j].real(),phi_new[g*CudaMultiSpline->num_split_splines+j].imag());
      }
  }
  abort();*/
}
#endif


////////////////////////////////
// Hybrid evaluation routines //
////////////////////////////////
// template<typename StorageType> void
// EinsplineSetHybrid<StorageType>::sort_electrons(std::vector<PosType> &pos)
// {
//   int nw = pos.size();
//   if (nw > CurrentWalkers)
//     resize_cuda(nw);

//   AtomicPolyJobs_CPU.clear();
//   AtomicSplineJobs_CPU.clear();
//   rhats_CPU.clear();
//   int numAtomic = 0;

//   // First, sort electrons into three categories:
//   // 1) Interstitial region with 3D-Bsplines
//   // 2) Atomic region near origin:      polynomial  radial functions
//   // 3) Atomic region not near origin:  1D B-spline radial functions

//   for (int i=0; i<newpos.size(); i++) {
//     PosType r = newpos[i];
//     // Note: this assumes that the atomic radii are smaller than the simulation cell radius.
//     for (int j=0; j<AtomicOrbitals.size(); j++) {
// 	AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[j];
// 	PosType dr = r - orb.Pos;
// 	PosType u = PrimLattice.toUnit(dr);
// 	for (int k=0; k<OHMMS_DIM; k++)
// 	  u[k] -= round(u[k]);
// 	dr = PrimLattice.toCart(u);
// 	RealType dist2 = dot (dr,dr);
// 	if (dist2 < orb.PolyRadius * orb.PolyRadius) {
// 	  AtomicPolyJob<CTS::RealType> job;
// 	  RealType dist = std::sqrt(dist2);
// 	  job.dist = dist;
// 	  RealType distInv = 1.0/dist;
// 	  for (int k=0; k<OHMMS_DIM; k++) {
// 	    CTS::RealType x = distInv*dr[k];
// 	    job.rhat[k] = distInv * dr[k];
// 	    rhats_CPU.push_back(x);
// 	  }
// 	  job.lMax = orb.lMax;
// 	  job.YlmIndex  = i;
// 	  job.PolyOrder = orb.PolyOrder;
// 	  //job.PolyCoefs = orb.PolyCoefs;
// 	  AtomicPolyJobs_CPU.push_back(job);
// 	  numAtomic++;
// 	}
// 	else if (dist2 < orb.CutoffRadius*orb.CutoffRadius) {
// 	  AtomicSplineJob<CTS::RealType> job;
// 	   RealType dist = std::sqrt(dist2);
// 	  job.dist = dist;
// 	  RealType distInv = 1.0/dist;
// 	  for (int k=0; k<OHMMS_DIM; k++) {
// 	    CTS::RealType x = distInv*dr[k];
// 	    job.rhat[k] = distInv * dr[k];
// 	    rhats_CPU.push_back(x);
// 	  }
// 	  job.lMax      = orb.lMax;
// 	  job.YlmIndex  = i;
// 	  job.phi       = phi[i];
// 	  job.grad_lapl = grad_lapl[i];
// 	  //job.PolyCoefs = orb.PolyCoefs;
// 	  AtomicSplineJobs_CPU.push_back(job);
// 	  numAtomic++;
// 	}
// 	else { // Regular 3D B-spline job
// 	  BsplinePos_CPU.push_back (r);

// 	}
//     }
//   }


// }

template<>
void EinsplineSetHybrid<double>::resize_cuda(int numwalkers)
{
  EinsplineSetExtended<double>::resize_cuda(numwalkers);
  CurrentWalkers = numwalkers;
  // Resize Ylm temporaries
  // Find lMax;
  lMax = -1;
  for (int i = 0; i < AtomicOrbitals.size(); i++)
    lMax = std::max(AtomicOrbitals[i].lMax, lMax);
  numlm  = (lMax + 1) * (lMax + 1);
  Ylm_BS = ((numlm + 15) / 16) * 16;
  Ylm_GPU.resize(numwalkers * Ylm_BS * 3);
  Ylm_ptr_GPU.resize(numwalkers);
  Ylm_ptr_CPU.resize(numwalkers);
  dYlm_dtheta_ptr_GPU.resize(numwalkers);
  dYlm_dtheta_ptr_CPU.resize(numwalkers);
  dYlm_dphi_ptr_GPU.resize(numwalkers);
  dYlm_dphi_ptr_CPU.resize(numwalkers);
  rhats_CPU.resize(OHMMS_DIM * numwalkers);
  rhats_GPU.resize(OHMMS_DIM * numwalkers);
  HybridJobs_GPU.resize(numwalkers);
  HybridData_GPU.resize(numwalkers);
  for (int iw = 0; iw < numwalkers; iw++)
  {
    Ylm_ptr_CPU[iw]         = Ylm_GPU.data() + (3 * iw + 0) * Ylm_BS;
    dYlm_dtheta_ptr_CPU[iw] = Ylm_GPU.data() + (3 * iw + 1) * Ylm_BS;
    dYlm_dphi_ptr_CPU[iw]   = Ylm_GPU.data() + (3 * iw + 2) * Ylm_BS;
  }
  Ylm_ptr_GPU         = Ylm_ptr_CPU;
  dYlm_dtheta_ptr_GPU = dYlm_dtheta_ptr_CPU;
  dYlm_dphi_ptr_GPU   = dYlm_dphi_ptr_CPU;
  // Resize AtomicJob temporaries
  // AtomicPolyJobs_GPU.resize(numwalkers);
  // AtomicSplineJobs_GPU.resize(numwalkers);
}


template<>
void EinsplineSetHybrid<std::complex<double>>::resize_cuda(int numwalkers)
{
  EinsplineSetExtended<std::complex<double>>::resize_cuda(numwalkers);
  CurrentWalkers = numwalkers;
  // Resize Ylm temporaries
  // Find lMax;
  lMax = -1;
  for (int i = 0; i < AtomicOrbitals.size(); i++)
    lMax = std::max(AtomicOrbitals[i].lMax, lMax);
  numlm  = (lMax + 1) * (lMax + 1);
  Ylm_BS = ((2 * numlm + 15) / 16) * 16;
  Ylm_GPU.resize(numwalkers * Ylm_BS * 3);
  Ylm_ptr_GPU.resize(numwalkers);
  Ylm_ptr_CPU.resize(numwalkers);
  dYlm_dtheta_ptr_GPU.resize(numwalkers);
  dYlm_dtheta_ptr_CPU.resize(numwalkers);
  dYlm_dphi_ptr_GPU.resize(numwalkers);
  dYlm_dphi_ptr_CPU.resize(numwalkers);
  rhats_CPU.resize(OHMMS_DIM * numwalkers);
  rhats_GPU.resize(OHMMS_DIM * numwalkers);
  HybridJobs_GPU.resize(numwalkers);
  HybridData_GPU.resize(numwalkers);
  for (int iw = 0; iw < numwalkers; iw++)
  {
    Ylm_ptr_CPU[iw]         = Ylm_GPU.data() + (3 * iw + 0) * Ylm_BS;
    dYlm_dtheta_ptr_CPU[iw] = Ylm_GPU.data() + (3 * iw + 1) * Ylm_BS;
    dYlm_dphi_ptr_CPU[iw]   = Ylm_GPU.data() + (3 * iw + 2) * Ylm_BS;
  }
  Ylm_ptr_GPU         = Ylm_ptr_CPU;
  dYlm_dtheta_ptr_GPU = dYlm_dtheta_ptr_CPU;
  dYlm_dphi_ptr_GPU   = dYlm_dphi_ptr_CPU;
  // Resize AtomicJob temporaries
  // AtomicPolyJobs_GPU.resize(numwalkers);
  // AtomicSplineJobs_GPU.resize(numwalkers);
}


// Vectorized evaluation functions
#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<Walker_t*>& walkers,
                                          int iat,
                                          gpu::device_vector<CTS::RealType*>& phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, int iat,\n"
              << " gpu::device_vector<CTS::RealType*> &phi) not implemented.\n";
  abort();
}

template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<Walker_t*>& walkers,
                                          std::vector<PosType>& newpos,
                                          gpu::device_vector<CTS::RealType*>& phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << " (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos, \n"
              << " gpu::device_vector<CTS::RealType*> &phi) not implemented.\n";
  abort();
}

#else
template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<Walker_t*>& walkers,
                                          int iat,
                                          gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, int iat,\n"
              << " gpu::device_vector<CTS::ComplexType*> &phi) not implemented.\n";
  abort();
}

template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<Walker_t*>& walkers,
                                          std::vector<PosType>& newpos,
                                          gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << "  (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,\n"
              << "   gpu::device_vector<CTS::ComplexType*> &phi) not implemented.\n";
  abort();
}
#endif


#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<Walker_t*>& walkers,
                                          std::vector<PosType>& newpos,
                                          gpu::device_vector<CTS::RealType*>& phi,
                                          gpu::device_vector<CTS::RealType*>& grad_lapl,
                                          int row_stride)
{
  int N = newpos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw = 0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;
  // hostPos = cudapos;
  // for (int i=0; i<newpos.size(); i++)
  //   std::cerr << "newPos[" << i << "] = " << newpos[i] << std::endl;
  // gpu::host_vector<CTS::RealType> IonPos_CPU(IonPos_GPU.size());
  // IonPos_CPU = IonPos_GPU;
  // for (int i=0; i<IonPos_CPU.size()/3; i++)
  //   fprintf (stderr, "ion[%d] = [%10.6f %10.6f %10.6f]\n",
  // 	       i, IonPos_CPU[3*i+0], IonPos_CPU[3*i+2], IonPos_CPU[3*i+2]);
  // std::cerr << "cudapos.size()        = " << cudapos.size() << std::endl;
  // std::cerr << "IonPos.size()         = " << IonPos_GPU.size() << std::endl;
  // std::cerr << "PolyRadii.size()      = " << PolyRadii_GPU.size() << std::endl;
  // std::cerr << "CutoffRadii.size()    = " << CutoffRadii_GPU.size() << std::endl;
  // std::cerr << "AtomicOrbitals.size() = " << AtomicOrbitals.size() << std::endl;
  // std::cerr << "L_cuda.size()         = " << L_cuda.size() << std::endl;
  // std::cerr << "Linv_cuda.size()      = " << Linv_cuda.size() << std::endl;
  // std::cerr << "HybridJobs_GPU.size() = " << HybridJobs_GPU.size() << std::endl;
  // std::cerr << "rhats_GPU.size()      = " << rhats_GPU.size() << std::endl;

  MakeHybridJobList<CTS::RealType>((CTS::RealType*)cudapos.data(), N, IonPos_GPU.data(), PolyRadii_GPU.data(),
                                   CutoffRadii_GPU.data(), AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(), HybridData_GPU.data());

  CalcYlmRealCuda<CTS::RealType>(rhats_GPU.data(), HybridJobs_GPU.data(), Ylm_ptr_GPU.data(),
                                 dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(), lMax, newpos.size());

  evaluate3DSplineReal<CTS::RealType>(HybridJobs_GPU.data(), (CTS::RealType*)cudapos.data(),
                                      (CTS::RealType*)CudakPoints_reduced.data(), CudaMultiSpline, Linv_cuda.data(),
                                      phi.data(), grad_lapl.data(), row_stride, NumOrbitals, newpos.size());

  evaluateHybridSplineReal<CTS::RealType>(HybridJobs_GPU.data(), rhats_GPU.data(), Ylm_ptr_GPU.data(),
                                          dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(),
                                          AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
                                          (CTS::RealType*)CudakPoints_reduced.data(), phi.data(), grad_lapl.data(),
                                          row_stride, NumOrbitals, newpos.size(), lMax);
#ifdef HYBRID_DEBUG
  gpu::host_vector<CTS::RealType*> phi_CPU(phi.size()), grad_lapl_CPU(phi.size());
  phi_CPU       = phi;
  grad_lapl_CPU = grad_lapl;
  gpu::host_vector<CTS::RealType> vals_CPU(NumOrbitals), GL_CPU(4 * row_stride);
  gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  HybridJobs_CPU = HybridJobs_GPU;
  gpu::host_vector<HybridData<CTS::RealType>> HybridData_CPU(HybridData_GPU.size());
  HybridData_CPU = HybridData_GPU;
  rhats_CPU      = rhats_GPU;
  for (int iw = 0; iw < newpos.size(); iw++)
    if (false && HybridJobs_CPU[iw] == ATOMIC_POLY_JOB)
    {
      ValueVector CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector CPUgrad(NumOrbitals);
      HybridData<CTS::RealType>& d = HybridData_CPU[iw];
      AtomicOrbital<double>& atom  = AtomicOrbitals[d.ion];
      atom.evaluate(newpos[iw], CPUvals, CPUgrad, CPUlapl);
      cudaMemcpy(&vals_CPU[0], phi_CPU[iw], NumOrbitals * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(&GL_CPU[0], grad_lapl_CPU[iw], 4 * row_stride * sizeof(float), cudaMemcpyDeviceToHost);
      // fprintf (stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n",
      // 	 iw, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
      double mindist = 1.0e5;
      for (int ion = 0; ion < AtomicOrbitals.size(); ion++)
      {
        PosType disp = newpos[iw] - AtomicOrbitals[ion].Pos;
        PosType u    = PrimLattice.toUnit(disp);
        PosType img;
        for (int i = 0; i < OHMMS_DIM; i++)
          u[i] -= round(u[i]);
        disp        = PrimLattice.toCart(u);
        double dist = std::sqrt(dot(disp, disp));
        if (dist < AtomicOrbitals[ion].CutoffRadius)
          mindist = dist;
      }
      if (std::fabs(mindist - d.dist) > 1.0e-3)
        fprintf(stderr, "CPU dist = %1.8f  GPU dist = %1.8f\n", mindist, d.dist);
      for (int j = 0; j < NumOrbitals; j++)
      {
        //	  if (isnan(vals_CPU[j])) {
        if (true || isnan(GL_CPU[0 * row_stride + j]))
        {
          std::cerr << "iw = " << iw << std::endl;
          fprintf(stderr, "rhat[%d] = [%10.6f %10.6f %10.6f] dist = %10.6e\n", iw, rhats_CPU[3 * iw + 0],
                  rhats_CPU[3 * iw + 1], rhats_CPU[3 * iw + 2], d.dist);
          fprintf(stderr, "val[%2d]  = %10.5e %10.5e\n", j, vals_CPU[j], CPUvals[j]);
          fprintf(stderr, "grad[%2d] = %10.5e %10.5e  %10.5e %10.5e  %10.5e %10.5e\n", j, GL_CPU[0 * row_stride + j],
                  CPUgrad[j][0], GL_CPU[1 * row_stride + j], CPUgrad[j][1], GL_CPU[2 * row_stride + j], CPUgrad[j][2]);
          fprintf(stderr, "lapl[%2d] = %10.5e %10.5e\n", j, GL_CPU[3 * row_stride + j], CPUlapl[j]);
        }
      }
    }
    else if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
    {
      std::cerr << "HalfG = " << HalfG << std::endl;
      ValueVector CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector CPUgrad(NumOrbitals);
      PosType ru(PrimLattice.toUnit(newpos[iw]));
      PosType img;
      int sign = 0;
      for (int i = 0; i < 3; i++)
      {
        img[i] = std::floor(ru[i]);
        ru[i] -= img[i];
        sign += HalfG[i] * (int)img[i];
      }
      EinsplineMultiEval(MultiSpline, ru, CPUvals, CPUgrad, storage_hess_vector_);
      cudaMemcpy(&vals_CPU[0], phi_CPU[iw], NumOrbitals * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(&GL_CPU[0], grad_lapl_CPU[iw], 4 * row_stride * sizeof(float), cudaMemcpyDeviceToHost);
      for (int j = 0; j < NumOrbitals; j++)
      {
        CPUgrad[j] = dot(PrimLattice.G, CPUgrad[j]);
        CPUlapl[j] = trace(storage_hess_vector_[j], GGt);
        if (sign & 1)
        {
          CPUvals[j] *= -1.0;
          CPUgrad[j] *= -1.0;
          CPUlapl[j] *= -1.0;
        }
        fprintf(stderr, "\nGPU=%10.6f  %10.6f %10.6f %10.6f  %10.6f\n", vals_CPU[j], GL_CPU[0 * row_stride + j],
                GL_CPU[1 * row_stride + j], GL_CPU[2 * row_stride + j], GL_CPU[3 * row_stride + j]);
        fprintf(stderr, "CPU=%10.6f  %10.6f %10.6f %10.6f  %10.6f sign = %d\n", CPUvals[j], CPUgrad[j][0],
                CPUgrad[j][1], CPUgrad[j][2], CPUlapl[j], sign);
        if (std::isnan(GL_CPU[0 * row_stride + j]))
        {
          std::cerr << "r[" << iw << "] = " << newpos[iw] << std::endl;
          std::cerr << "iw = " << iw << std::endl;
          fprintf(stderr, "3D val[%2d]  = %10.5e %10.5e\n", j, vals_CPU[j], CPUvals[j]);
          fprintf(stderr, "3D grad[%2d] = %10.5e %10.5e  %10.5e %10.5e  %10.5e %10.5e\n", j, GL_CPU[0 * row_stride + j],
                  CPUgrad[j][0], GL_CPU[1 * row_stride + j], CPUgrad[j][1], GL_CPU[2 * row_stride + j], CPUgrad[j][2]);
          fprintf(stderr, "3D lapl[%2d] = %10.5e %10.5e\n", j, GL_CPU[3 * row_stride + j], CPUlapl[j]);
        }
      }
    }
  gpu::host_vector<float> Ylm_CPU(Ylm_GPU.size());
  Ylm_CPU   = Ylm_GPU;
  rhats_CPU = rhats_GPU;
  for (int i = 0; i < rhats_CPU.size() / 3; i++)
    fprintf(stderr, "rhat[%d] = [%10.6f %10.6f %10.6f]\n", i, rhats_CPU[3 * i + 0], rhats_CPU[3 * i + 1],
            rhats_CPU[3 * i + 2]);
  //    gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  HybridJobs_CPU = HybridJobs_GPU;
  //    gpu::host_vector<HybridDataFloat> HybridData_CPU(HybridData_GPU.size());
  HybridData_CPU = HybridData_GPU;
  std::cerr << "Before loop.\n";
  for (int i = 0; i < newpos.size(); i++)
    if (HybridJobs_CPU[i] != BSPLINE_3D_JOB)
    {
      std::cerr << "Inside if.\n";
      PosType rhat(rhats_CPU[3 * i + 0], rhats_CPU[3 * i + 1], rhats_CPU[3 * i + 2]);
      AtomicOrbital<double>& atom = AtomicOrbitals[HybridData_CPU[i].ion];
      int numlm                   = (atom.lMax + 1) * (atom.lMax + 1);
      std::vector<double> Ylm(numlm), dYlm_dtheta(numlm), dYlm_dphi(numlm);
      atom.CalcYlm(rhat, Ylm, dYlm_dtheta, dYlm_dphi);
      for (int lm = 0; lm < numlm; lm++)
      {
        fprintf(stderr, "lm=%3d  Ylm_CPU=%8.5f  Ylm_GPU=%8.5f\n", lm, Ylm[lm], Ylm_CPU[3 * i * Ylm_BS + lm]);
      }
    }
  fprintf(stderr, " N  img      dist    ion    lMax\n");
  for (int i = 0; i < HybridData_CPU.size(); i++)
  {
    HybridData<CTS::RealType>& d = HybridData_CPU[i];
    fprintf(stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n", i, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
  }
#endif
  // int N = newpos.size();
  // if (N > CurrentWalkers)
  //   resize_cuda(N);
  // AtomicPolyJobs_CPU.clear();
  // AtomicSplineJobs_CPU.clear();
  // rhats_CPU.clear();
  // BsplinePos_CPU.clear();
  // BsplineVals_CPU.clear();
  // BsplineGradLapl_CPU.clear();
  // int numAtomic = 0;
  // // First, sort electrons into three categories:
  // // 1) Interstitial region with 3D B-splines
  // // 2) Atomic region near origin:      polynomial  radial functions
  // // 3) Atomic region not near origin:  1D B-spline radial functions
  // for (int i=0; i<newpos.size(); i++) {
  //   PosType r = newpos[i];
  //   // Note: this assumes that the atomic radii are smaller than the simulation cell radius.
  //   for (int j=0; j<AtomicOrbitals.size(); j++) {
  // 	AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[j];
  // 	PosType dr = r - orb.Pos;
  // 	PosType u = PrimLattice.toUnit(dr);
  // 	for (int k=0; k<OHMMS_DIM; k++)
  // 	  u[k] -= round(u[k]);
  // 	dr = PrimLattice.toCart(u);
  // 	RealType dist2 = dot (dr,dr);
  // 	if (dist2 < orb.PolyRadius * orb.PolyRadius) {
  // 	  AtomicPolyJob<CTS::RealType> job;
  // 	  RealType dist = std::sqrt(dist2);
  // 	  job.dist = dist;
  // 	  RealType distInv = 1.0/dist;
  // 	  for (int k=0; k<OHMMS_DIM; k++) {
  // 	    CTS::RealType x = distInv*dr[k];
  // 	    job.rhat[k] = distInv * dr[k];
  // 	    rhats_CPU.push_back(x);
  // 	  }
  // 	  job.lMax = orb.lMax;
  // 	  job.YlmIndex  = i;
  // 	  job.PolyOrder = orb.PolyOrder;
  // 	  //job.PolyCoefs = orb.PolyCoefs;
  // 	  AtomicPolyJobs_CPU.push_back(job);
  // 	  numAtomic++;
  // 	}
  // 	else if (dist2 < orb.CutoffRadius*orb.CutoffRadius) {
  // 	  AtomicSplineJob<CTS::RealType> job;
  // 	   RealType dist = std::sqrt(dist2);
  // 	  job.dist = dist;
  // 	  RealType distInv = 1.0/dist;
  // 	  for (int k=0; k<OHMMS_DIM; k++) {
  // 	    CTS::RealType x = distInv*dr[k];
  // 	    job.rhat[k] = distInv * dr[k];
  // 	    rhats_CPU.push_back(x);
  // 	  }
  // 	  job.lMax      = orb.lMax;
  // 	  job.YlmIndex  = i;
  // 	  job.phi       = phi[i];
  // 	  job.grad_lapl = grad_lapl[i];
  // 	  //job.PolyCoefs = orb.PolyCoefs;
  // 	  AtomicSplineJobs_CPU.push_back(job);
  // 	  numAtomic++;
  // 	}
  // 	else { // Regular 3D B-spline job
  // 	  BsplinePos_CPU
  // 	}
  //   }
  // }
  // //////////////////////////////////
  // // First, evaluate 3D B-splines //
  // //////////////////////////////////
  // int N = newpos.size();
  // CTS::RealType plus_minus[2] = {1.0, -1.0};
  // if (cudapos.size() < N) {
  //   hostPos.resize(N);
  //   cudapos.resize(N);
  //   hostSign.resize(N);
  //   cudaSign.resize(N);
  // }
  // for (int iw=0; iw < N; iw++) {
  //   PosType r = newpos[iw];
  //   PosType ru(PrimLattice.toUnit(r));
  //   int image[OHMMS_DIM];
  //   for (int i=0; i<OHMMS_DIM; i++) {
  // 	RealType img = std::floor(ru[i]);
  // 	ru[i] -= img;
  // 	image[i] = (int) img;
  //   }
  //   int sign = 0;
  //   for (int i=0; i<OHMMS_DIM; i++)
  // 	sign += HalfG[i] * image[i];
  //   hostSign[iw] = plus_minus[sign&1];
  //   hostPos[iw] = ru;
  // }
  // cudapos = hostPos;
  // cudaSign = hostSign;
  // eval_multi_multi_UBspline_3d_cuda
  //   (CudaMultiSpline, (CTS::RealType*)(cudapos.data()), cudaSign.data(),
  //    phi.data(), N);
  // ////////////////////////////////////////////////////////////
  // // Next, evaluate spherical harmonics for atomic orbitals //
  // ////////////////////////////////////////////////////////////
  // // Evaluate Ylms
  // if (rhats_CPU.size()) {
  //   rhats_GPU = rhats_CPU;
  //   CalcYlmComplexCuda(rhats_GPU.data(), Ylm_ptr_GPU.data(), dYlm_dtheta_ptr_GPU.data(),
  // 			 dYlm_dphi_ptr_GPU.data(), lMax, numAtomic);
  //   std::cerr << "Calculated Ylms.\n";
  // }
  // ///////////////////////////////////////
  // // Next, evaluate 1D spline orbitals //
  // ///////////////////////////////////////
  // if (AtomicSplineJobs_CPU.size()) {
  //   AtomicSplineJobs_GPU = AtomicSplineJobs_CPU;
  // }
  // ///////////////////////////////////////////
  // // Next, evaluate 1D polynomial orbitals //
  // ///////////////////////////////////////////
  // if (AtomicSplineJobs_CPU.size()) {
  //   AtomicPolyJobs_GPU = AtomicPolyJobs_CPU;
  // }
}

#else
template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<Walker_t*>& walkers,
                                          std::vector<PosType>& newpos,
                                          gpu::device_vector<CTS::ComplexType*>& phi,
                                          gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                                          int row_stride)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << "(std::vector<std::unique_ptr<Walker_t<>> &walkers, std::vector<PosType> &newpos, \n"
              << " gpu::device_vector<CTS::ComplexType*> &phi,\n"
              << " gpu::device_vector<CTS::ComplexType*> &grad_lapl, int row_stride)\n"
              << "     is not yet implemented.\n";
  abort();
}
#endif

#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi)
{
  int N = pos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw = 0; iw < N; iw++)
    hostPos[iw] = pos[iw];
  cudapos = hostPos;

  MakeHybridJobList<CTS::RealType>((CTS::RealType*)cudapos.data(), N, IonPos_GPU.data(), PolyRadii_GPU.data(),
                                   CutoffRadii_GPU.data(), AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(), HybridData_GPU.data());

  CalcYlmRealCuda<CTS::RealType>(rhats_GPU.data(), HybridJobs_GPU.data(), Ylm_ptr_GPU.data(),
                                 dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(), lMax, pos.size());

  evaluateHybridSplineReal<CTS::RealType>(HybridJobs_GPU.data(), Ylm_ptr_GPU.data(), AtomicOrbitals_GPU.data(),
                                          HybridData_GPU.data(), (CTS::RealType*)CudakPoints_reduced.data(), phi.data(),
                                          NumOrbitals, pos.size(), lMax);

  evaluate3DSplineReal<CTS::RealType>(HybridJobs_GPU.data(), (CTS::RealType*)cudapos.data(),
                                      (CTS::RealType*)CudakPoints_reduced.data(), CudaMultiSpline, Linv_cuda.data(),
                                      phi.data(), NumOrbitals, pos.size());
}
#else
template<>
void EinsplineSetHybrid<double>::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << "(std::vector<PosType> &pos, gpu::device_vector<CTS::ComplexType*> &phi)\n"
              << "     is not yet implemented.\n";
  abort();
}
#endif

#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                        int iat,
                                                        gpu::device_vector<CTS::RealType*>& phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> "
                 "&walkers, int iat,\n"
              << "			                            gpu::device_vector<CTS::RealType*> &phi)\n"
              << "not yet implemented.\n";
}

template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                        std::vector<PosType>& newpos,
                                                        gpu::device_vector<CTS::RealType*>& phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, "
                 "std::vector<PosType> &newpos,\n"
              << "			                            gpu::device_vector<CTS::RealType*> &phi)\n"
              << "not yet implemented.\n";
}

#else
template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                        int iat,
                                                        gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> "
                 "&walkers, int iat,\n"
              << "			                            gpu::device_vector<CTS::ComplexType*> &phi)\n"
              << "not yet implemented.\n";
}

template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                        std::vector<PosType>& newpos,
                                                        gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, "
                 "std::vector<PosType> ,\n"
              << "			                            gpu::device_vector<CTS::ComplexType*> &phi)\n"
              << "not yet implemented.\n";
}
#endif


#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                        std::vector<PosType>& newpos,
                                                        gpu::device_vector<CTS::RealType*>& phi,
                                                        gpu::device_vector<CTS::RealType*>& grad_lapl,
                                                        int row_stride)
{
  static int numAtomic = 0;
  static int num3D     = 0;
  int N                = newpos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw = 0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;

  MakeHybridJobList<CTS::RealType>((CTS::RealType*)cudapos.data(), N, IonPos_GPU.data(), PolyRadii_GPU.data(),
                                   CutoffRadii_GPU.data(), AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(), HybridData_GPU.data());

  CalcYlmComplexCuda<CTS::RealType>(rhats_GPU.data(), HybridJobs_GPU.data(), Ylm_ptr_GPU.data(),
                                    dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(), lMax, newpos.size());

  evaluate3DSplineComplexToReal<CTS::RealType>(HybridJobs_GPU.data(), (CTS::RealType*)cudapos.data(),
                                               (CTS::RealType*)CudakPoints.data(), CudaMakeTwoCopies.data(),
                                               CudaMultiSpline, Linv_cuda.data(), phi.data(), grad_lapl.data(),
                                               row_stride, CudaMakeTwoCopies.size(), newpos.size());

  evaluateHybridSplineComplexToReal<CTS::RealType>(HybridJobs_GPU.data(), rhats_GPU.data(), Ylm_ptr_GPU.data(),
                                                   dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(),
                                                   AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
                                                   (CTS::RealType*)CudakPoints_reduced.data(), CudaMakeTwoCopies.data(),
                                                   (CTS::RealType**)phi.data(), grad_lapl.data(), row_stride,
                                                   CudaMakeTwoCopies.size(), newpos.size(), lMax);

#ifdef HYBRID_DEBUG
  // gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  // HybridJobs_CPU = HybridJobs_GPU;
  // int M = MakeTwoCopies.size();
  // // ComplexValueVector CPUzvals(M), CPUzlapl(M);
  // // ComplexGradVector CPUzgrad(M);
  // for (int iw=0; iw<newpos.size(); iw++) {
  //   if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  // 	num3D++;
  //   else
  // 	numAtomic++;
  //   // bool atomic=false;
  //   // for (int ion=0; ion<AtomicOrbitals.size(); ion++)
  //   // 	if (AtomicOrbitals[ion].evaluate
  //   // 	    (newpos[iw], CPUzvals, CPUzgrad, CPUzlapl)) {
  //   // 	  atomic = true;
  //   // 	  break;
  //   // 	}
  //   // if (atomic)
  //   // 	numAtomic++;
  //   // else
  //   // 	num3D++;
  //   // if (atomic && HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  //   // 	cerr << "Error!  Used BSPLINE_3D when I should have used atomic.\n";
  //   // else if (!atomic && HybridJobs_CPU[iw] != BSPLINE_3D_JOB)
  //   // 	cerr << "Error!  Used atomic when I should have used 3D.\n";
  // }
  // // fprintf (stderr, "Num atomic = %d  Num 3D = %d\n",
  // // 	     numAtomic, num3D);
  // if (numAtomic + num3D > 100000) {
  //   fprintf (stderr, "Num atomic = %d  Num 3D = %d\n",
  // 	       numAtomic, num3D);
  //   fprintf (stderr, "Percent atomic = %1.5f\%\n",
  // 	       100.0*(double)numAtomic / (double)(numAtomic+num3D));
  //   numAtomic = num3D = 0;
  // }
  gpu::host_vector<CTS::RealType*> phi_CPU(phi.size()), grad_lapl_CPU(phi.size());
  phi_CPU       = phi;
  grad_lapl_CPU = grad_lapl;
  gpu::host_vector<CTS::RealType> vals_CPU(NumOrbitals), GL_CPU(4 * row_stride);
  gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  HybridJobs_CPU = HybridJobs_GPU;
  gpu::host_vector<HybridData<CTS::RealType>> HybridData_CPU(HybridData_GPU.size());
  HybridData_CPU = HybridData_GPU;
  rhats_CPU      = rhats_GPU;
  // for (int iw=0; iw<newpos.size(); iw++)
  //   fprintf (stderr, "rhat[%d] = [%10.6f %10.6f %10.6f]\n",
  // 	       iw, rhats_CPU[3*iw+0], rhats_CPU[3*iw+1], rhats_CPU[3*iw+2]);
  for (int iw = 0; iw < newpos.size(); iw++)
    if (false && HybridJobs_CPU[iw] == ATOMIC_POLY_JOB)
    {
      //if (HybridJobs_CPU[iw] != BSPLINE_3D_JOB && std::abs(rhats_CPU[3*iw+2]) < 1.0e-6) {
      int M = MakeTwoCopies.size();
      ComplexValueVector CPUzvals(M), CPUzlapl(M);
      ComplexGradVector CPUzgrad(M);
      ValueVector CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector CPUgrad(NumOrbitals);
      HybridData<CTS::RealType>& d              = HybridData_CPU[iw];
      AtomicOrbital<std::complex<double>>& atom = AtomicOrbitals[d.ion];
      atom.evaluate(newpos[iw], CPUzvals, CPUzgrad, CPUzlapl);
      int index = 0;
      for (int i = 0; i < storage_value_vector_.size(); i++)
      {
        CPUvals[index] = CPUzvals[i].real();
        CPUlapl[index] = CPUzlapl[i].real();
        for (int j = 0; j < OHMMS_DIM; j++)
          CPUgrad[index][j] = CPUzgrad[i][j].real();
        index++;
        if (MakeTwoCopies[i])
        {
          CPUvals[index] = CPUzvals[i].imag();
          CPUlapl[index] = CPUzlapl[i].imag();
          for (int j = 0; j < OHMMS_DIM; j++)
            CPUgrad[index][j] = CPUzgrad[i][j].imag();
          index++;
        }
      }
      cudaMemcpy(&vals_CPU[0], phi_CPU[iw], NumOrbitals * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(&GL_CPU[0], grad_lapl_CPU[iw], 4 * row_stride * sizeof(float), cudaMemcpyDeviceToHost);
      // fprintf (stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n",
      // 	 iw, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
      double mindist = 1.0e5;
      for (int ion = 0; ion < AtomicOrbitals.size(); ion++)
      {
        PosType disp = newpos[iw] - AtomicOrbitals[ion].Pos;
        PosType u    = PrimLattice.toUnit(disp);
        PosType img;
        for (int i = 0; i < OHMMS_DIM; i++)
          u[i] -= round(u[i]);
        disp        = PrimLattice.toCart(u);
        double dist = std::sqrt(dot(disp, disp));
        if (dist < AtomicOrbitals[ion].CutoffRadius)
          mindist = dist;
      }
      if (std::fabs(mindist - d.dist) > 1.0e-3)
        fprintf(stderr, "CPU dist = %1.8f  GPU dist = %1.8f\n", mindist, d.dist);
      fprintf(stderr, "rhat[%d] = [%10.6f %10.6f %10.6f] dist = %10.6e\n", iw, rhats_CPU[3 * iw + 0],
              rhats_CPU[3 * iw + 1], rhats_CPU[3 * iw + 2], d.dist);
      for (int j = 0; j < NumOrbitals; j++)
      {
        //	  if (isnan(vals_CPU[j])) {
        if (true || isnan(GL_CPU[0 * row_stride + j]))
        {
          std::cerr << "iw = " << iw << std::endl;
          fprintf(stderr, "val[%2d]  = %10.5e %10.5e\n", j, vals_CPU[j], CPUvals[j]);
          fprintf(stderr, "grad[%2d] = %10.5e %10.5e  %10.5e %10.5e  %10.5e %10.5e\n", j, GL_CPU[0 * row_stride + j],
                  CPUgrad[j][0], GL_CPU[1 * row_stride + j], CPUgrad[j][1], GL_CPU[2 * row_stride + j], CPUgrad[j][2]);
          fprintf(stderr, "lapl[%2d] = %10.5e %10.5e\n", j, GL_CPU[3 * row_stride + j], CPUlapl[j]);
        }
      }
    }
    else if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
    {
      ComplexValueVector CPUzvals(NumOrbitals), CPUzlapl(NumOrbitals);
      ComplexGradVector CPUzgrad(NumOrbitals);
      ValueVector CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector CPUgrad(NumOrbitals);
      PosType ru(PrimLattice.toUnit(newpos[iw]));
      for (int i = 0; i < 3; i++)
        ru[i] -= std::floor(ru[i]);
      EinsplineMultiEval(MultiSpline, ru, CPUzvals, CPUzgrad, storage_hess_vector_);
      for (int j = 0; j < MakeTwoCopies.size(); j++)
      {
        CPUzgrad[j] = dot(PrimLattice.G, CPUzgrad[j]);
        CPUzlapl[j] = trace(storage_hess_vector_[j], GGt);
      }
      // Add e^-ikr phase to B-spline orbitals
      std::complex<double> eye(0.0, 1.0);
      for (int j = 0; j < MakeTwoCopies.size(); j++)
      {
        std::complex<double> u                            = CPUzvals[j];
        TinyVector<std::complex<double>, OHMMS_DIM> gradu = CPUzgrad[j];
        std::complex<double> laplu                        = CPUzlapl[j];
        PosType k                                         = kPoints[j];
        TinyVector<std::complex<double>, OHMMS_DIM> ck;
        for (int n = 0; n < OHMMS_DIM; n++)
          ck[n] = k[n];
        double s, c;
        double phase = -dot(newpos[iw], k);
        qmcplusplus::sincos(phase, &s, &c);
        std::complex<double> e_mikr(c, s);
        CPUzvals[j] = e_mikr * u;
        CPUzgrad[j] = e_mikr * (-eye * u * ck + gradu);
        CPUzlapl[j] = e_mikr * (-dot(k, k) * u - 2.0 * eye * dot(ck, gradu) + laplu);
      }
      int index = 0;
      for (int i = 0; i < MakeTwoCopies.size(); i++)
      {
        CPUvals[index] = CPUzvals[i].real();
        CPUlapl[index] = CPUzlapl[i].real();
        for (int j = 0; j < OHMMS_DIM; j++)
          CPUgrad[index][j] = CPUzgrad[i][j].real();
        index++;
        if (MakeTwoCopies[i])
        {
          CPUvals[index] = CPUzvals[i].imag();
          CPUlapl[index] = CPUzlapl[i].imag();
          for (int j = 0; j < OHMMS_DIM; j++)
            CPUgrad[index][j] = CPUzgrad[i][j].imag();
          index++;
        }
      }
      cudaMemcpy(&vals_CPU[0], phi_CPU[iw], NumOrbitals * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(&GL_CPU[0], grad_lapl_CPU[iw], 4 * row_stride * sizeof(float), cudaMemcpyDeviceToHost);
      // for (int i=0; i<4*row_stride; i++)
      //   fprintf (stderr, "%d %10.5e\n", i, GL_CPU[i]);
      static long int numgood = 0, numbad = 0;
      for (int j = 0; j < NumOrbitals; j++)
      {
        double lap_ratio = GL_CPU[3 * row_stride + j] / CPUlapl[j];
        if (std::abs(GL_CPU[3 * row_stride + j] - CPUlapl[j]) > 1.0e-4)
        {
          fprintf(stderr, "Error:  CPU laplacian = %1.8e  GPU = %1.8e\n", CPUlapl[j], GL_CPU[3 * row_stride + j]);
          fprintf(stderr, "        CPU value     = %1.8e  GPU = %1.8e\n", CPUvals[j], vals_CPU[j]);
          fprintf(stderr, "u = %1.8f %1.8f %1.8f \n", ru[0], ru[1], ru[2]);
          std::cerr << "iw = " << iw << std::endl;
          numbad++;
        }
        else
          numgood++;
        if (numbad + numgood >= 100000)
        {
          double percent_bad = 100.0 * (double)numbad / (double)(numbad + numgood);
          fprintf(stderr, "Percent bad = %1.8f\n", percent_bad);
          numbad = numgood = 0;
        }
        // if (true || isnan(GL_CPU[0*row_stride+j])) {
        //   // std::cerr << "r[" << iw << "] = " << newpos[iw] << std::endl;
        //   // std::cerr << "iw = " << iw << std::endl;
        //   fprintf (stderr, "\n3D      val[%2d]  = %10.5e %10.5e\n",
        // 	     j, vals_CPU[j], CPUvals[j]);
        //   fprintf (stderr, "3D GPU grad[%2d] = %10.5e %10.5e %10.5e\n", j,
        // 	     GL_CPU[0*row_stride+j],
        // 	     GL_CPU[1*row_stride+j],
        // 	     GL_CPU[2*row_stride+j]);
        //   fprintf (stderr, "3D CPU grad[%2d] = %10.5e %10.5e %10.5e\n", j,
        // 	     CPUgrad[j][0],
        // 	     CPUgrad[j][1],
        // 	     CPUgrad[j][2]);
        //   fprintf (stderr, "3D     lapl[%2d] = %10.5e %10.5e\n",
        // 	     j, GL_CPU[3*row_stride+j], CPUlapl[j]);
        // }
      }
    }
  gpu::host_vector<float> Ylm_CPU(Ylm_GPU.size());
  Ylm_CPU   = Ylm_GPU;
  rhats_CPU = rhats_GPU;
  for (int i = 0; i < rhats_CPU.size() / 3; i++)
    fprintf(stderr, "rhat[%d] = [%10.6f %10.6f %10.6f]\n", i, rhats_CPU[3 * i + 0], rhats_CPU[3 * i + 1],
            rhats_CPU[3 * i + 2]);
  //    gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  //    HybridJobs_CPU = HybridJobs_GPU;
  //    gpu::host_vector<HybridDataFloat> HybridData_CPU(HybridData_GPU.size());
  //    HybridData_CPU = HybridData_GPU;
  std::cerr << "Before loop.\n";
  for (int i = 0; i < newpos.size(); i++)
    if (HybridJobs_CPU[i] != BSPLINE_3D_JOB)
    {
      std::cerr << "Inside if.\n";
      PosType rhat(rhats_CPU[3 * i + 0], rhats_CPU[3 * i + 1], rhats_CPU[3 * i + 2]);
      AtomicOrbital<std::complex<double>>& atom = AtomicOrbitals[HybridData_CPU[i].ion];
      int numlm                                 = (atom.lMax + 1) * (atom.lMax + 1);
      std::vector<double> Ylm(numlm), dYlm_dtheta(numlm), dYlm_dphi(numlm);
      atom.CalcYlm(rhat, Ylm, dYlm_dtheta, dYlm_dphi);
      for (int lm = 0; lm < numlm; lm++)
      {
        fprintf(stderr, "lm=%3d  Ylm_CPU=%8.5f  Ylm_GPU=%8.5f\n", lm, Ylm[lm], Ylm_CPU[3 * i * Ylm_BS + lm]);
      }
    }
  fprintf(stderr, " N  img      dist    ion    lMax\n");
  for (int i = 0; i < HybridData_CPU.size(); i++)
  {
    HybridData<CTS::RealType>& d = HybridData_CPU[i];
    fprintf(stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n", i, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
  }
#endif
}

#else
template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<Walker_t*>& walkers,
                                                        std::vector<PosType>& newpos,
                                                        gpu::device_vector<CTS::ComplexType*>& phi,
                                                        gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                                                        int row_stride)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate \n"
              << "(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos, \n"
              << " gpu::device_vector<CTS::ComplexType*> &phi,\n"
              << " gpu::device_vector<CTS::ComplexType*> &grad_lapl,\n"
              << " int row_stride)\n"
              << "not yet implemented.\n";
  abort();
}
#endif

#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<PosType>& pos,
                                                        gpu::device_vector<CTS::RealType*>& phi)
{
  int N = pos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw = 0; iw < N; iw++)
    hostPos[iw] = pos[iw];
  cudapos = hostPos;

  MakeHybridJobList<CTS::RealType>((CTS::RealType*)cudapos.data(), N, IonPos_GPU.data(), PolyRadii_GPU.data(),
                                   CutoffRadii_GPU.data(), AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(), HybridData_GPU.data());

  CalcYlmComplexCuda<CTS::RealType>(rhats_GPU.data(), HybridJobs_GPU.data(), Ylm_ptr_GPU.data(),
                                    dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(), lMax, pos.size());

  evaluate3DSplineComplexToReal<CTS::RealType>(HybridJobs_GPU.data(), (CTS::RealType*)cudapos.data(),
                                               (CTS::RealType*)CudakPoints.data(), CudaMakeTwoCopies.data(),
                                               CudaMultiSpline, Linv_cuda.data(), phi.data(), CudaMakeTwoCopies.size(),
                                               pos.size());

  evaluateHybridSplineComplexToReal<CTS::RealType> //NLPP
      (HybridJobs_GPU.data(), Ylm_ptr_GPU.data(), AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
       (CTS::RealType*)CudakPoints_reduced.data(), CudaMakeTwoCopies.data(), (CTS::RealType**)phi.data(),
       CudaMakeTwoCopies.size(), pos.size(), lMax);

  // gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  // HybridJobs_CPU = HybridJobs_GPU;
  // int M = CudaMakeTwoCopies.size();
  // ComplexValueVector CPUzvals(M);
  // for (int iw=0; iw<pos.size(); iw++)
  //   if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  // 	cerr << "Error:  used BSPLINE_3D for PP eval.  Walker = "
  // 	     << iw << "\n";
  // std::cerr << "pos.size() = " << pos.size() << std::endl;
  // for (int iw=0; iw<pos.size(); iw++) {
  //   // if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  //   // 	num3D++;
  //   // else
  //   // 	numAtomic++;
  //   bool atomic=false;
  //   for (int ion=0; ion<AtomicOrbitals.size(); ion++)
  //   	if (AtomicOrbitals[ion].evaluate (pos[iw], CPUzvals)) {
  //   	  atomic = true;
  //   	  break;
  //   	}
  //   // if (atomic)
  //   // 	numAtomic++;
  //   // else
  //   // 	num3D++;
  //   if (atomic && HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  //   	cerr << "Error!  Used BSPLINE_3D when I should have used atomic.\n";
  //   else if (!atomic && HybridJobs_CPU[iw] != BSPLINE_3D_JOB)
  //   	cerr << "Error!  Used atomic when I should have used 3D.\n";
  // }
  //////////////////////////
  // This must be tested! //
  //////////////////////////
}

#else
template<>
void EinsplineSetHybrid<std::complex<double>>::evaluate(std::vector<PosType>& pos,
                                                        gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate \n"
              << "(std::vector<PosType> &pos, gpu::device_vector<CTS::ComplexType*> &phi)\n"
              << "not yet implemented.\n";
}
#endif

template<>
void EinsplineSetExtended<double>::finalizeConstruction()
{
  app_log() << "Copying einspline orbitals to GPU.\n";
  create_multi_UBspline_3d_cuda(MultiSpline, CudaMultiSpline);
  app_log() << "Successful copy.\n";
  get_split_spline_pointers();
  // Destroy original CPU spline
  // HACK HACK HACK
  //destroy_Bspline (MultiSpline);
  L_host.resize(9);
  Linv_host.resize(9);
  Linv_cuda.resize(9, 1.0, split_splines);
  L_cuda.resize(9, 1.0, split_splines);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      L_host[i * 3 + j]    = PrimLattice.R(i, j);
      Linv_host[i * 3 + j] = PrimLattice.G(i, j);
    }
  L_cuda    = L_host;
  Linv_cuda = Linv_host;
  if (split_splines)
  {
    cudaMemAdvise(L_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(Linv_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetReadMostly, 0);
    for (unsigned int i = 0; i < gpu::device_group_size; i++)
    {
      cudaMemAdvise(L_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetAccessedBy, gpu::device_group_numbers[i]);
      cudaMemAdvise(Linv_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetAccessedBy,
                    gpu::device_group_numbers[i]);
      cudaMemPrefetchAsync(Linv_cuda.data(), 9 * sizeof(CTS::RealType), gpu::device_group_numbers[i]);
    }
  }
}

template<>
void EinsplineSetExtended<std::complex<double>>::finalizeConstruction()
{
  app_log() << "Copying einspline orbitals to GPU.\n";
  create_multi_UBspline_3d_cuda(MultiSpline, CudaMultiSpline);
  app_log() << "Successful copy.\n";
  get_split_spline_pointers();
  // Destroy original CPU spline
  // HACK HACK HACK
  //destroy_Bspline (MultiSpline);
  L_host.resize(9);
  Linv_host.resize(9);
  Linv_cuda.resize(9, 1.0, split_splines);
  L_cuda.resize(9, 1.0, split_splines);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      L_host[i * 3 + j]    = PrimLattice.R(i, j);
      Linv_host[i * 3 + j] = PrimLattice.G(i, j);
    }
  L_cuda    = L_host;
  Linv_cuda = Linv_host;
  if (split_splines)
  {
    cudaMemAdvise(L_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(Linv_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetReadMostly, 0);
    for (unsigned int i = 0; i < gpu::device_group_size; i++)
    {
      cudaMemAdvise(L_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetAccessedBy, gpu::device_group_numbers[i]);
      cudaMemAdvise(Linv_cuda.data(), 9 * sizeof(CTS::RealType), cudaMemAdviseSetAccessedBy,
                    gpu::device_group_numbers[i]);
      cudaMemPrefetchAsync(Linv_cuda.data(), 9 * sizeof(CTS::RealType), gpu::device_group_numbers[i]);
    }
  }
}


template<>
void EinsplineSetHybrid<double>::finalizeConstruction()
{
  EinsplineSetExtended<double>::finalizeConstruction();
  // Setup B-spline Acuda matrix in constant memory
  init_atomic_cuda();
  gpu::host_vector<AtomicOrbitalCuda<CTS::RealType>> AtomicOrbitals_CPU;
  const int BS = 16;
  NumOrbitals  = getOrbitalSetSize();
  // Bump up the stride to be a multiple of 512-bit bus width
  int lm_stride = ((NumOrbitals + BS - 1) / BS) * BS;
  AtomicSplineCoefs_GPU.resize(AtomicOrbitals.size());
  AtomicPolyCoefs_GPU.resize(AtomicOrbitals.size());
  std::vector<CTS::RealType> IonPos_CPU, CutoffRadii_CPU, PolyRadii_CPU;
  for (int iat = 0; iat < AtomicOrbitals.size(); iat++)
  {
    app_log() << "Copying real atomic orbitals for ion " << iat << " to GPU memory.\n";
    AtomicOrbital<double>& atom = AtomicOrbitals[iat];
    for (int i = 0; i < OHMMS_DIM; i++)
      IonPos_CPU.push_back(atom.Pos[i]);
    CutoffRadii_CPU.push_back(atom.CutoffRadius);
    PolyRadii_CPU.push_back(atom.PolyRadius);
    AtomicOrbitalCuda<CTS::RealType> atom_cuda;
    atom_cuda.lMax                                = atom.lMax;
    int numlm                                     = (atom.lMax + 1) * (atom.lMax + 1);
    atom_cuda.lm_stride                           = lm_stride;
    atom_cuda.spline_stride                       = numlm * lm_stride;
    AtomicOrbital<double>::SplineType& cpu_spline = *atom.get_radial_spline();
    atom_cuda.spline_dr_inv                       = cpu_spline.x_grid.delta_inv;
    int Ngrid                                     = cpu_spline.x_grid.num;
    int spline_size                               = 2 * atom_cuda.spline_stride * (Ngrid + 2);
    gpu::host_vector<CTS::RealType> spline_coefs(spline_size);
    AtomicSplineCoefs_GPU[iat].resize(spline_size);
    atom_cuda.spline_coefs = AtomicSplineCoefs_GPU[iat].data();
    // Reorder and copy splines to GPU memory
    for (int igrid = 0; igrid < Ngrid; igrid++)
      for (int lm = 0; lm < numlm; lm++)
        for (int orb = 0; orb < NumOrbitals; orb++)
        {
          // Convert splines to Hermite spline form, because
          // B-spline form in single precision with dense grids
          // leads to large truncation error
          double u = 1.0 / 6.0 *
              (1.0 * cpu_spline.coefs[(igrid + 0) * cpu_spline.x_stride + orb * numlm + lm] +
               4.0 * cpu_spline.coefs[(igrid + 1) * cpu_spline.x_stride + orb * numlm + lm] +
               1.0 * cpu_spline.coefs[(igrid + 2) * cpu_spline.x_stride + orb * numlm + lm]);
          double d2u = cpu_spline.x_grid.delta_inv * cpu_spline.x_grid.delta_inv *
              (1.0 * cpu_spline.coefs[(igrid + 0) * cpu_spline.x_stride + orb * numlm + lm] +
               -2.0 * cpu_spline.coefs[(igrid + 1) * cpu_spline.x_stride + orb * numlm + lm] +
               1.0 * cpu_spline.coefs[(igrid + 2) * cpu_spline.x_stride + orb * numlm + lm]);
          spline_coefs[(2 * igrid + 0) * atom_cuda.spline_stride + lm * atom_cuda.lm_stride + orb] = u;
          spline_coefs[(2 * igrid + 1) * atom_cuda.spline_stride + lm * atom_cuda.lm_stride + orb] = d2u;
        }
    AtomicSplineCoefs_GPU[iat] = spline_coefs;
    atom_cuda.poly_stride      = numlm * atom_cuda.lm_stride;
    atom_cuda.poly_order       = atom.PolyOrder;
    int poly_size              = (atom.PolyOrder + 1) * atom_cuda.poly_stride;
    gpu::host_vector<CTS::RealType> poly_coefs(poly_size);
    AtomicPolyCoefs_GPU[iat].resize(poly_size);
    atom_cuda.poly_coefs = AtomicPolyCoefs_GPU[iat].data();
    for (int lm = 0; lm < numlm; lm++)
      for (int n = 0; n < atom.PolyOrder; n++)
        for (int orb = 0; orb < NumOrbitals; orb++)
          poly_coefs[n * atom_cuda.poly_stride + lm * atom_cuda.lm_stride + orb] = atom.get_poly_coefs()(n, orb, lm);
    AtomicPolyCoefs_GPU[iat] = poly_coefs;
    AtomicOrbitals_CPU.push_back(atom_cuda);
  }
  AtomicOrbitals_GPU = AtomicOrbitals_CPU;
  IonPos_GPU         = IonPos_CPU;
  CutoffRadii_GPU    = CutoffRadii_CPU;
  PolyRadii_GPU      = PolyRadii_CPU;
}

template<>
void EinsplineSetHybrid<std::complex<double>>::finalizeConstruction()
{
  EinsplineSetExtended<std::complex<double>>::finalizeConstruction();
  // Setup B-spline Acuda matrix in constant memory
  init_atomic_cuda();
  gpu::host_vector<AtomicOrbitalCuda<CTS::RealType>> AtomicOrbitals_CPU;
  const int BS = 16;
  NumOrbitals  = getOrbitalSetSize();
  // Bump up the stride to be a multiple of 512-bit bus width
  int lm_stride = ((2 * NumOrbitals + BS - 1) / BS) * BS;
  AtomicSplineCoefs_GPU.resize(AtomicOrbitals.size());
  AtomicPolyCoefs_GPU.resize(AtomicOrbitals.size());
  std::vector<CTS::RealType> IonPos_CPU, CutoffRadii_CPU, PolyRadii_CPU;
  for (int iat = 0; iat < AtomicOrbitals.size(); iat++)
  {
    app_log() << "Copying atomic orbitals for ion " << iat << " to GPU memory.\n";
    AtomicOrbital<std::complex<double>>& atom = AtomicOrbitals[iat];
    for (int i = 0; i < OHMMS_DIM; i++)
      IonPos_CPU.push_back(atom.Pos[i]);
    CutoffRadii_CPU.push_back(atom.CutoffRadius);
    PolyRadii_CPU.push_back(atom.PolyRadius);
    AtomicOrbitalCuda<CTS::RealType> atom_cuda;
    atom_cuda.lMax                                              = atom.lMax;
    int numlm                                                   = (atom.lMax + 1) * (atom.lMax + 1);
    atom_cuda.lm_stride                                         = lm_stride;
    atom_cuda.spline_stride                                     = numlm * lm_stride;
    AtomicOrbital<std::complex<double>>::SplineType& cpu_spline = *atom.get_radial_spline();
    atom_cuda.spline_dr_inv                                     = cpu_spline.x_grid.delta_inv;
    int Ngrid                                                   = cpu_spline.x_grid.num;
    int spline_size                                             = 2 * atom_cuda.spline_stride * (Ngrid + 2);
    gpu::host_vector<CTS::RealType> spline_coefs(spline_size);
    AtomicSplineCoefs_GPU[iat].resize(spline_size);
    atom_cuda.spline_coefs = AtomicSplineCoefs_GPU[iat].data();
    // Reorder and copy splines to GPU memory
    for (int igrid = 0; igrid < Ngrid; igrid++)
      for (int lm = 0; lm < numlm; lm++)
        for (int orb = 0; orb < NumOrbitals; orb++)
        {
          // Convert splines to Hermite spline form, because
          // B-spline form in single precision with dense grids
          // leads to large truncation error
          std::complex<double> u = 1.0 / 6.0 *
              (1.0 * cpu_spline.coefs[(igrid + 0) * cpu_spline.x_stride + orb * numlm + lm] +
               4.0 * cpu_spline.coefs[(igrid + 1) * cpu_spline.x_stride + orb * numlm + lm] +
               1.0 * cpu_spline.coefs[(igrid + 2) * cpu_spline.x_stride + orb * numlm + lm]);
          std::complex<double> d2u = cpu_spline.x_grid.delta_inv * cpu_spline.x_grid.delta_inv *
              (1.0 * cpu_spline.coefs[(igrid + 0) * cpu_spline.x_stride + orb * numlm + lm] +
               -2.0 * cpu_spline.coefs[(igrid + 1) * cpu_spline.x_stride + orb * numlm + lm] +
               1.0 * cpu_spline.coefs[(igrid + 2) * cpu_spline.x_stride + orb * numlm + lm]);
          spline_coefs[((2 * igrid + 0) * atom_cuda.spline_stride + lm * atom_cuda.lm_stride + 2 * orb) + 0] = u.real();
          spline_coefs[((2 * igrid + 0) * atom_cuda.spline_stride + lm * atom_cuda.lm_stride + 2 * orb) + 1] = u.imag();
          spline_coefs[((2 * igrid + 1) * atom_cuda.spline_stride + lm * atom_cuda.lm_stride + 2 * orb) + 0] =
              d2u.real();
          spline_coefs[((2 * igrid + 1) * atom_cuda.spline_stride + lm * atom_cuda.lm_stride + 2 * orb) + 1] =
              d2u.imag();
        }
    AtomicSplineCoefs_GPU[iat] = spline_coefs;
    atom_cuda.poly_stride      = numlm * atom_cuda.lm_stride;
    atom_cuda.poly_order       = atom.PolyOrder;
    int poly_size              = (atom.PolyOrder + 1) * atom_cuda.poly_stride;
    gpu::host_vector<CTS::RealType> poly_coefs(poly_size);
    AtomicPolyCoefs_GPU[iat].resize(poly_size);
    atom_cuda.poly_coefs = &AtomicPolyCoefs_GPU[iat].data()[0];
    for (int lm = 0; lm < numlm; lm++)
      for (int n = 0; n < atom.PolyOrder; n++)
        for (int orb = 0; orb < NumOrbitals; orb++)
        {
          poly_coefs[n * atom_cuda.poly_stride + lm * atom_cuda.lm_stride + 2 * orb + 0] =
              atom.get_poly_coefs()(n, orb, lm).real();
          poly_coefs[n * atom_cuda.poly_stride + lm * atom_cuda.lm_stride + 2 * orb + 1] =
              atom.get_poly_coefs()(n, orb, lm).imag();
        }
    AtomicPolyCoefs_GPU[iat] = poly_coefs;
    AtomicOrbitals_CPU.push_back(atom_cuda);
  }
  AtomicOrbitals_GPU = AtomicOrbitals_CPU;
  IonPos_GPU         = IonPos_CPU;
  CutoffRadii_GPU    = CutoffRadii_CPU;
  PolyRadii_GPU      = PolyRadii_CPU;
}


} // namespace qmcplusplus
