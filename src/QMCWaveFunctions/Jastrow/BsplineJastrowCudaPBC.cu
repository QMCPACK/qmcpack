//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#define MAX_SPLINES 100
#include <stdio.h>
#include <config.h>
#include "BsplineJastrowCudaPBC.h"
#include "../../CUDA/gpu_misc.h"


bool AisInitializedPBC = false;

static bool CMC_profile = false;

void CMC_profileSample(const char *function, float msec)
{
  if (strcmp(function, "two_body_NLratios_PBC"))
  {
    return;
  }
  printf("%s: %1.3e msec\n", function, msec);
}

#define CMC_PROFILING_BEGIN() \
  cudaMemcpyToSymbolAsync(CMC_L,    lattice,    sizeof(CMC_L),    0, cudaMemcpyDeviceToDevice, gpu::kernelStream); \
  cudaMemcpyToSymbolAsync(CMC_Linv, latticeInv, sizeof(CMC_Linv), 0, cudaMemcpyDeviceToDevice, gpu::kernelStream); \
  cudaEvent_t start; \
  cudaEvent_t stop; \
  if (CMC_profile) { \
    cudaEventCreate(&start); \
    cudaEventCreate(&stop); \
    cudaGetLastError(); \
    cudaEventRecord(start); \
  }

#define CMC_PROFILING_END(lineno) \
  if (CMC_profile) { \
    cudaEventRecord(stop); \
    cudaEventSynchronize(stop); \
    float time = 0.0f; \
    cudaEventElapsedTime(&time, start, stop); \
    cudaEventDestroy(start); \
    cudaEventDestroy(stop); \
    CMC_profileSample(__FUNCTION__, time); \
  } \
  cudaError_t error = cudaGetLastError(); \
  if (error) { printf("%s\nCUDA ERROR!!! Detected at end of CMC_PROFILING_END in BsplineJastrowCudaPBC line %d!!!\n", cudaGetErrorString(error), lineno); exit(1); }

#define COPY_LATTICE_DP_TO_SP() \
  cudaMemcpyToSymbolAsync(CMC_L,    lattice,    sizeof(CMC_L),    0, cudaMemcpyDeviceToDevice, gpu::kernelStream); \
  cudaMemcpyToSymbolAsync(CMC_Linv, latticeInv, sizeof(CMC_Linv), 0, cudaMemcpyDeviceToDevice, gpu::kernelStream); \
  cudaError_t error = cudaGetLastError(); \
  if (error) { printf("%s\nCUDA ERROR!!! Detected at end of COPY_LATTICE_DP_TO_SP in BsplineJastrowCudaPBC!!!\n", cudaGetErrorString(error)); exit(1); }

static __constant__ CUDA_PRECISION CMC_L[3][3];
static __constant__ CUDA_PRECISION CMC_Linv[3][3];

template<typename T>
__device__ __forceinline__
T CMC_min_dist_fast(T& __restrict__ x, T& __restrict__ y, T& __restrict__ z)
{
  T u0 = CMC_Linv[0][0]*x + CMC_Linv[1][0]*y + CMC_Linv[2][0]*z;
  T u1 = CMC_Linv[0][1]*x + CMC_Linv[1][1]*y + CMC_Linv[2][1]*z;
  T u2 = CMC_Linv[0][2]*x + CMC_Linv[1][2]*y + CMC_Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  x = CMC_L[0][0]*u0 + CMC_L[1][0]*u1 + CMC_L[2][0]*u2;
  y = CMC_L[0][1]*u0 + CMC_L[1][1]*u1 + CMC_L[2][1]*u2;
  z = CMC_L[0][2]*u0 + CMC_L[1][2]*u1 + CMC_L[2][2]*u2;
  return sqrt (x*x + y*y + z*z);
}

template<typename T>
__device__ __forceinline__
T CMC_min_dist_only(T x, T y, T z)
{
  T u0 = CMC_Linv[0][0]*x + CMC_Linv[1][0]*y + CMC_Linv[2][0]*z;
  T u1 = CMC_Linv[0][1]*x + CMC_Linv[1][1]*y + CMC_Linv[2][1]*z;
  T u2 = CMC_Linv[0][2]*x + CMC_Linv[1][2]*y + CMC_Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  x = CMC_L[0][0]*u0 + CMC_L[1][0]*u1 + CMC_L[2][0]*u2;
  y = CMC_L[0][1]*u0 + CMC_L[1][1]*u1 + CMC_L[2][1]*u2;
  z = CMC_L[0][2]*u0 + CMC_L[1][2]*u1 + CMC_L[2][2]*u2;
  T d2min = x*x + y*y + z*z;
#pragma unroll
  for (int i = -1; i <= 1; i++)
  {
#pragma unroll
    for (int j = -1; j <= 1; j++)
    {
#pragma unroll
      for (int k = -1; k <= 1; k++)
      {
        T xnew = CMC_L[0][0]*(u0+i) + CMC_L[1][0]*(u1+j) + CMC_L[2][0]*(u2+k);
        T ynew = CMC_L[0][1]*(u0+i) + CMC_L[1][1]*(u1+j) + CMC_L[2][1]*(u2+k);
        T znew = CMC_L[0][2]*(u0+i) + CMC_L[1][2]*(u1+j) + CMC_L[2][2]*(u2+k);
        T d2 = xnew*xnew + ynew*ynew + znew*znew;
        d2min = min(d2, d2min);
      }
    }
  }
  return sqrt(d2min);
}


template<typename T>
__device__ __forceinline__
T CMC_min_dist(T& __restrict__ x, T& __restrict__ y, T& __restrict__ z)
{
  T u0 = CMC_Linv[0][0]*x + CMC_Linv[1][0]*y + CMC_Linv[2][0]*z;
  T u1 = CMC_Linv[0][1]*x + CMC_Linv[1][1]*y + CMC_Linv[2][1]*z;
  T u2 = CMC_Linv[0][2]*x + CMC_Linv[1][2]*y + CMC_Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  x = CMC_L[0][0]*u0 + CMC_L[1][0]*u1 + CMC_L[2][0]*u2;
  y = CMC_L[0][1]*u0 + CMC_L[1][1]*u1 + CMC_L[2][1]*u2;
  z = CMC_L[0][2]*u0 + CMC_L[1][2]*u1 + CMC_L[2][2]*u2;
  T d2min = x*x + y*y + z*z;
#pragma unroll
  for (int i = -1; i <= 1; i++)
  {
#pragma unroll
    for (int j = -1; j <= 1; j++)
    {
#pragma unroll
      for (int k = -1; k <= 1; k++)
      {
        T xnew = CMC_L[0][0]*(u0+i) + CMC_L[1][0]*(u1+j) + CMC_L[2][0]*(u2+k);
        T ynew = CMC_L[0][1]*(u0+i) + CMC_L[1][1]*(u1+j) + CMC_L[2][1]*(u2+k);
        T znew = CMC_L[0][2]*(u0+i) + CMC_L[1][2]*(u1+j) + CMC_L[2][2]*(u2+k);
        T d2new = xnew*xnew + ynew*ynew + znew*znew;
        if (d2new < d2min)
        {
          x = xnew;
          y = ynew;
          z = znew;
          d2min = d2new;
        }
      }
    }
  }
  return sqrt(d2min);
}


// void
// createCudaSplines (float rmax, int N,
// 		   float f[], float df[], float d2f[],
// 		   int &fSpline, int &dfSpline, int &d2fSpline)
// {
//   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
//   cudaArray *fArray, *dfArray, *d2fArray;
//   cudaMallocArray(  &fArray, &channelDesc, N);
//   cudaMallocArray( &dfArray, &channelDesc, N);
//   cudaMallocArray(&d2fArray, &channelDesc, N);

//   cudaMemcpyToArray(fArray,  N,1,  f,N*sizeof(float),cudaMemcpyHostToDevice);
//   cudaMemcpyToArray(dfArray, N,1, df,N*sizeof(float),cudaMemcpyHostToDevice);
//   cudaMemcpyToArray(d2fArray,N,1,d2f,N*sizeof(float),cudaMemcpyHostToDevice);


//   cudaBindTextureToArray(texSplines[fSpline=curTex++], fArray);
//   cudaBindTextureToArray(texSplines[dfSpline=curTex++], dfArray);
//   cudaBindTextureToArray(texSplines[d2fSpline=curTex++], d2fArray);
// }


template<typename T>
__device__ __forceinline__
T min_dist (T& __restrict__ x, T& __restrict__ y, T& __restrict__ z,
            T const L[3][3], T const Linv[3][3])
{
  T u0 = Linv[0][0]*x + Linv[1][0]*y + Linv[2][0]*z;
  T u1 = Linv[0][1]*x + Linv[1][1]*y + Linv[2][1]*z;
  T u2 = Linv[0][2]*x + Linv[1][2]*y + Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  x = L[0][0]*u0 + L[1][0]*u1 + L[2][0]*u2;
  y = L[0][1]*u0 + L[1][1]*u1 + L[2][1]*u2;
  z = L[0][2]*u0 + L[1][2]*u1 + L[2][2]*u2;
  T d2min = x*x + y*y + z*z;
#pragma unroll
  for (int i = -1; i <= 1; i++)
  {
#pragma unroll
    for (int j = -1; j <= 1; j++)
    {
#pragma unroll
      for (int k = -1; k <= 1; k++)
      {
        T xnew = L[0][0]*(u0+i) + L[1][0]*(u1+j) + L[2][0]*(u2+k);
        T ynew = L[0][1]*(u0+i) + L[1][1]*(u1+j) + L[2][1]*(u2+k);
        T znew = L[0][2]*(u0+i) + L[1][2]*(u1+j) + L[2][2]*(u2+k);
        T d2 = xnew*xnew + ynew*ynew + znew*znew;
        if (d2 < d2min)
        {
          d2min = d2;
          x = xnew;
          y = ynew;
          z = znew;
        }
      }
    }
  }
  return sqrt(d2min);
}

template<typename T>
__device__ __forceinline__
T min_dist_fast (T& __restrict__ x, T& __restrict__ y, T& __restrict__ z,
                 T const L[3][3], T const Linv[3][3])
{
  T u0 = Linv[0][0]*x + Linv[1][0]*y + Linv[2][0]*z;
  T u1 = Linv[0][1]*x + Linv[1][1]*y + Linv[2][1]*z;
  T u2 = Linv[0][2]*x + Linv[1][2]*y + Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  x = L[0][0]*u0 + L[1][0]*u1 + L[2][0]*u2;
  y = L[0][1]*u0 + L[1][1]*u1 + L[2][1]*u2;
  z = L[0][2]*u0 + L[1][2]*u1 + L[2][2]*u2;
  return sqrt(x*x + y*y + z*z);
}

template<typename T>
__device__ __forceinline__
T min_dist (T& __restrict__ x, T& __restrict__ y, T& __restrict__ z,
            T const L[3][3], T const Linv[3][3],
            T const images[27][3])
{
  T u0 = Linv[0][0]*x + Linv[1][0]*y + Linv[2][0]*z;
  T u1 = Linv[0][1]*x + Linv[1][1]*y + Linv[2][1]*z;
  T u2 = Linv[0][2]*x + Linv[1][2]*y + Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  T xtmp = L[0][0]*u0 + L[1][0]*u1 + L[2][0]*u2;
  T ytmp = L[0][1]*u0 + L[1][1]*u1 + L[2][1]*u2;
  T ztmp = L[0][2]*u0 + L[1][2]*u1 + L[2][2]*u2;
  x = xtmp;
  y = ytmp;
  z = ztmp;
  T d2min = xtmp*xtmp + ytmp*ytmp + ztmp*ztmp;
  for (int i=0; i<27; i++)
  {
    T xnew = xtmp + images[i][0];
    T ynew = ytmp + images[i][1];
    T znew = ztmp + images[i][2];
    T d2 = xnew*xnew + ynew*ynew + znew*znew;
    if (d2 < d2min)
    {
      x = xnew;
      y = ynew;
      z = znew;
      d2min = d2;
    }
    // __syncthreads(); // XXXJCW: this doesn't appear to be needed
  }
  return sqrt(d2min);
}


template<typename T>
__device__ __forceinline__
T min_dist_only (T x, T y, T z,
                 T const L[3][3], T const Linv[3][3],
                 T const images[27][3])
{
  T u0 = Linv[0][0]*x + Linv[1][0]*y + Linv[2][0]*z;
  T u1 = Linv[0][1]*x + Linv[1][1]*y + Linv[2][1]*z;
  T u2 = Linv[0][2]*x + Linv[1][2]*y + Linv[2][2]*z;
  u0 -= rint(u0);
  u1 -= rint(u1);
  u2 -= rint(u2);
  x = L[0][0]*u0 + L[1][0]*u1 + L[2][0]*u2;
  y = L[0][1]*u0 + L[1][1]*u1 + L[2][1]*u2;
  z = L[0][2]*u0 + L[1][2]*u1 + L[2][2]*u2;
  T d2min = x*x + y*y + z*z;
  for (int i=0; i<27; i++)
  {
    T xnew = x + images[i][0];
    T ynew = y + images[i][1];
    T znew = z + images[i][2];
    T d2 = xnew*xnew + ynew*ynew + znew*znew;
    d2min = min (d2min, d2);
    // __syncthreads(); // XXXJCW: this doesn't appear to be needed
  }
  return sqrt(d2min);
}



__constant__ float AcudaSpline[48];
__constant__ double AcudaSpline_double[48];

void
cuda_spline_init_PBC()
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
                    3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
                    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
                    1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
                    0.0,     -0.5,      1.0,    -0.5,
                    0.0,      1.5,     -2.0,     0.0,
                    0.0,     -1.5,      1.0,     0.5,
                    0.0,      0.5,      0.0,     0.0,
                    0.0,      0.0,     -1.0,     1.0,
                    0.0,      0.0,      3.0,    -2.0,
                    0.0,      0.0,     -3.0,     1.0,
                    0.0,      0.0,      1.0,     0.0
                  };
  cudaMemcpyToSymbol(AcudaSpline, A_h, 48*sizeof(float), 0,
                     cudaMemcpyHostToDevice);
  double A_d[48] = {-1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
                    3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
                    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
                    1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
                    0.0,     -0.5,      1.0,    -0.5,
                    0.0,      1.5,     -2.0,     0.0,
                    0.0,     -1.5,      1.0,     0.5,
                    0.0,      0.5,      0.0,     0.0,
                    0.0,      0.0,     -1.0,     1.0,
                    0.0,      0.0,      3.0,    -2.0,
                    0.0,      0.0,     -3.0,     1.0,
                    0.0,      0.0,      1.0,     0.0
                   };
  cudaMemcpyToSymbol(AcudaSpline_double, A_d, 48*sizeof(double), 0,
                     cudaMemcpyHostToDevice);
  AisInitializedPBC = true;
}


template<typename T>
__device__ __forceinline__
T eval_1d_spline (T dist, T rmax, T drInv, T const A[4][4],
                  T const * __restrict__ coefs)
{
  T res;
  if (dist >= rmax)
  {
    res = (T)0;
  }
  else
  {
    T s = dist * drInv;
    T sf = floor(s);
    int index = (int)sf;
    T t = s - sf;
    T t2 = t*t;
    T t3 = t*t2;
    res = (coefs[index+0]*(A[0][0]*t3 + A[0][1]*t2 + A[0][2]*t + A[0][3]) +
           coefs[index+1]*(A[1][0]*t3 + A[1][1]*t2 + A[1][2]*t + A[1][3]) +
           coefs[index+2]*(A[2][0]*t3 + A[2][1]*t2 + A[2][2]*t + A[2][3]) +
           coefs[index+3]*(A[3][0]*t3 + A[3][1]*t2 + A[3][2]*t + A[3][3]));
  }
  return res;
}


template<typename T>
__device__ __forceinline__
T CMC_eval_1d_spline (T dist, T rmax, T drInv, T const * __restrict__ coefs)
{
  T res;
  if (dist >= rmax)
  {
    res = (T)0;
  }
  else
  {
    T s = dist * drInv;
    T sf = floor(s);
    int index = (int)sf;
    T t = s - sf;
    res = (coefs[index+0] * (((AcudaSpline[ 0] * t + AcudaSpline[ 1]) * t + AcudaSpline[ 2]) * t + AcudaSpline[ 3]) +
           coefs[index+1] * (((AcudaSpline[ 4] * t + AcudaSpline[ 5]) * t + AcudaSpline[ 6]) * t + AcudaSpline[ 7]) +
           coefs[index+2] * (((AcudaSpline[ 8] * t + AcudaSpline[ 9]) * t + AcudaSpline[10]) * t + AcudaSpline[11]) +
           coefs[index+3] * (((AcudaSpline[12] * t + AcudaSpline[13]) * t + AcudaSpline[14]) * t + AcudaSpline[15]) );
  }
  return res;
}

template<typename T>
__device__ __forceinline__
void eval_1d_spline_vgl (T dist, T rmax, T drInv, T const A[12][4],
                         T const * __restrict__ coefs,
                         T& __restrict__ u, T& __restrict__ du,
                         T& __restrict__ d2u)
{
  if (dist >= rmax)
  {
    u = du = d2u = (T)0;
  }
  else
  {
    T s = dist * drInv;
    T sf = floor (s);
    int index = (int)sf;
    T t = s - sf;
    T t2 = t*t;
    T t3 = t*t2;
    T c0 = coefs[index+0];
    T c1 = coefs[index+1];
    T c2 = coefs[index+2];
    T c3 = coefs[index+3];
    u = (c0 * (A[0][0]*t3 + A[0][1]*t2 + A[0][2]*t + A[0][3]) +
         c1 * (A[1][0]*t3 + A[1][1]*t2 + A[1][2]*t + A[1][3]) +
         c2 * (A[2][0]*t3 + A[2][1]*t2 + A[2][2]*t + A[2][3]) +
         c3 * (A[3][0]*t3 + A[3][1]*t2 + A[3][2]*t + A[3][3]));
    du = drInv *
         (c0 * (A[4][0]*t3 + A[4][1]*t2 + A[4][2]*t + A[4][3]) +
          c1 * (A[5][0]*t3 + A[5][1]*t2 + A[5][2]*t + A[5][3]) +
          c2 * (A[6][0]*t3 + A[6][1]*t2 + A[6][2]*t + A[6][3]) +
          c3 * (A[7][0]*t3 + A[7][1]*t2 + A[7][2]*t + A[7][3]));
    d2u = drInv*drInv *
          (c0 * (A[ 8][0]*t3 + A[ 8][1]*t2 + A[ 8][2]*t + A[ 8][3]) +
           c1 * (A[ 9][0]*t3 + A[ 9][1]*t2 + A[ 9][2]*t + A[ 9][3]) +
           c2 * (A[10][0]*t3 + A[10][1]*t2 + A[10][2]*t + A[10][3]) +
           c3 * (A[11][0]*t3 + A[11][1]*t2 + A[11][2]*t + A[11][3]));
  }
}

#define NAIVE_SCHEME  0
#define HORNER_SCHEME 1
#define ESTRIN_SCHEME 2

#define SCHEME2 HORNER_SCHEME
template<typename T>
__device__ __forceinline__
void CMC_eval_1d_spline_vgl (T dist, T rmax, T drInv,
                             T const * __restrict__ coefs,
                             T& __restrict__ u,
                             T& __restrict__ du,
                             T& __restrict__ d2u)
{
  if (dist >= rmax)
  {
    u = du = d2u = (T)0;
  }
  else
  {
    T s = dist * drInv;
    T sf = floor (s);
    int index = (int)sf;
    T t = s - sf;
    T c0 = coefs[index+0];
    T c1 = coefs[index+1];
    T c2 = coefs[index+2];
    T c3 = coefs[index+3];
#if (SCHEME2 == HORNER_SCHEME)
    u = (c0 * (((AcudaSpline[ 0] * t + AcudaSpline[ 1]) * t + AcudaSpline[ 2]) * t + AcudaSpline[ 3]) +
         c1 * (((AcudaSpline[ 4] * t + AcudaSpline[ 5]) * t + AcudaSpline[ 6]) * t + AcudaSpline[ 7]) +
         c2 * (((AcudaSpline[ 8] * t + AcudaSpline[ 9]) * t + AcudaSpline[10]) * t + AcudaSpline[11]) +
         c3 * (((AcudaSpline[12] * t + AcudaSpline[13]) * t + AcudaSpline[14]) * t + AcudaSpline[15]));
    du = drInv *
         (c0 * (((AcudaSpline[16] * t + AcudaSpline[17]) * t + AcudaSpline[18]) * t + AcudaSpline[19]) +
          c1 * (((AcudaSpline[20] * t + AcudaSpline[21]) * t + AcudaSpline[22]) * t + AcudaSpline[23]) +
          c2 * (((AcudaSpline[24] * t + AcudaSpline[25]) * t + AcudaSpline[26]) * t + AcudaSpline[27]) +
          c3 * (((AcudaSpline[28] * t + AcudaSpline[29]) * t + AcudaSpline[30]) * t + AcudaSpline[31]));
    d2u = drInv * drInv *
          (c0 * (((AcudaSpline[32] * t + AcudaSpline[33]) * t + AcudaSpline[34]) * t + AcudaSpline[35]) +
           c1 * (((AcudaSpline[36] * t + AcudaSpline[37]) * t + AcudaSpline[38]) * t + AcudaSpline[39]) +
           c2 * (((AcudaSpline[40] * t + AcudaSpline[41]) * t + AcudaSpline[42]) * t + AcudaSpline[43]) +
           c3 * (((AcudaSpline[44] * t + AcudaSpline[45]) * t + AcudaSpline[46]) * t + AcudaSpline[47]));
#elif (SCHEME2 == ESTRIN_SCHEME)
    T t2 = t*t;
    u = (c0 * ((AcudaSpline[ 0*4 + 0] * t + AcudaSpline[ 0*4 + 1]) * t2 + (AcudaSpline[ 0*4 + 2] * t + AcudaSpline[ 0*4 + 3])) +
         c1 * ((AcudaSpline[ 1*4 + 0] * t + AcudaSpline[ 1*4 + 1]) * t2 + (AcudaSpline[ 1*4 + 2] * t + AcudaSpline[ 1*4 + 3])) +
         c2 * ((AcudaSpline[ 2*4 + 0] * t + AcudaSpline[ 2*4 + 1]) * t2 + (AcudaSpline[ 2*4 + 2] * t + AcudaSpline[ 2*4 + 3])) +
         c3 * ((AcudaSpline[ 3*4 + 0] * t + AcudaSpline[ 3*4 + 1]) * t2 + (AcudaSpline[ 3*4 + 2] * t + AcudaSpline[ 3*4 + 3])) );
    du = drInv *
         (c0 * ((AcudaSpline[ 4*4 + 0] * t + AcudaSpline[ 4*4 + 1]) * t2 + (AcudaSpline[ 4*4 + 2] * t + AcudaSpline[ 4*4 + 3])) +
          c1 * ((AcudaSpline[ 5*4 + 0] * t + AcudaSpline[ 5*4 + 1]) * t2 + (AcudaSpline[ 5*4 + 2] * t + AcudaSpline[ 5*4 + 3])) +
          c2 * ((AcudaSpline[ 6*4 + 0] * t + AcudaSpline[ 6*4 + 1]) * t2 + (AcudaSpline[ 6*4 + 2] * t + AcudaSpline[ 6*4 + 3])) +
          c3 * ((AcudaSpline[ 7*4 + 0] * t + AcudaSpline[ 7*4 + 1]) * t2 + (AcudaSpline[ 7*4 + 2] * t + AcudaSpline[ 7*4 + 3])) );
    d2u = drInv * drInv *
          (c0 * ((AcudaSpline[ 8*4 + 0] * t + AcudaSpline[ 8*4 + 1]) * t2 + (AcudaSpline[ 8*4 + 2] * t + AcudaSpline[ 8*4 + 3])) +
           c1 * ((AcudaSpline[ 9*4 + 0] * t + AcudaSpline[ 9*4 + 1]) * t2 + (AcudaSpline[ 9*4 + 2] * t + AcudaSpline[ 9*4 + 3])) +
           c2 * ((AcudaSpline[10*4 + 0] * t + AcudaSpline[10*4 + 1]) * t2 + (AcudaSpline[10*4 + 2] * t + AcudaSpline[10*4 + 3])) +
           c3 * ((AcudaSpline[11*4 + 0] * t + AcudaSpline[11*4 + 1]) * t2 + (AcudaSpline[11*4 + 2] * t + AcudaSpline[11*4 + 3])) );
#else
    T t2 = t*t;
    T t3 = t*t2;
    u = (c0 * (AcudaSpline[ 0*4 + 0] * t3 + AcudaSpline[ 0*4 + 1] * t2 + AcudaSpline[ 0*4 + 2] * t  + AcudaSpline[ 0*4 + 3]) +
         c1 * (AcudaSpline[ 1*4 + 0] * t3 + AcudaSpline[ 1*4 + 1] * t2 + AcudaSpline[ 1*4 + 2] * t  + AcudaSpline[ 1*4 + 3]) +
         c2 * (AcudaSpline[ 2*4 + 0] * t3 + AcudaSpline[ 2*4 + 1] * t2 + AcudaSpline[ 2*4 + 2] * t  + AcudaSpline[ 2*4 + 3]) +
         c3 * (AcudaSpline[ 3*4 + 0] * t3 + AcudaSpline[ 3*4 + 1] * t2 + AcudaSpline[ 3*4 + 2] * t  + AcudaSpline[ 3*4 + 3]));
    du = drInv *
         (c0 * (AcudaSpline[ 4*4 + 0] * t3 + AcudaSpline[ 4*4 + 1] * t2 + AcudaSpline[ 4*4 + 2] * t  + AcudaSpline[ 4*4 + 3]) +
          c1 * (AcudaSpline[ 5*4 + 0] * t3 + AcudaSpline[ 5*4 + 1] * t2 + AcudaSpline[ 5*4 + 2] * t  + AcudaSpline[ 5*4 + 3]) +
          c2 * (AcudaSpline[ 6*4 + 0] * t3 + AcudaSpline[ 6*4 + 1] * t2 + AcudaSpline[ 6*4 + 2] * t  + AcudaSpline[ 6*4 + 3]) +
          c3 * (AcudaSpline[ 7*4 + 0] * t3 + AcudaSpline[ 7*4 + 1] * t2 + AcudaSpline[ 7*4 + 2] * t  + AcudaSpline[ 7*4 + 3]));
    d2u = drInv * drInv *
          (c0 * (AcudaSpline[ 8*4 + 0] * t3 + AcudaSpline[ 8*4 + 1] * t2 + AcudaSpline[ 8*4 + 2] * t  + AcudaSpline[ 8*4 + 3]) +
           c1 * (AcudaSpline[ 9*4 + 0] * t3 + AcudaSpline[ 9*4 + 1] * t2 + AcudaSpline[ 9*4 + 2] * t  + AcudaSpline[ 9*4 + 3]) +
           c2 * (AcudaSpline[10*4 + 0] * t3 + AcudaSpline[10*4 + 1] * t2 + AcudaSpline[10*4 + 2] * t  + AcudaSpline[10*4 + 3]) +
           c3 * (AcudaSpline[11*4 + 0] * t3 + AcudaSpline[11*4 + 1] * t2 + AcudaSpline[11*4 + 2] * t  + AcudaSpline[11*4 + 3]));
#endif
  }
}

#define MAX_COEFS 32
template<typename T, int BS >
__global__ void
two_body_sum_PBC_kernel(T **R, int e1_first, int e1_last,
                        int e2_first, int e2_last,
                        T *spline_coefs, int numCoefs, T rMax,
                        T *lattice, T* latticeInv, T* sum)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  if (tid == 0)
    myR = R[blockIdx.x];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  __shared__ T r1[BS][3], r2[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[4][4];
  if (tid < 16)
    A[tid>>2][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int N1 = e1_last - e1_first + 1;
  int N2 = e2_last - e2_first + 1;
  int NB1 = N1/BS + ((N1 % BS) ? 1 : 0);
  int NB2 = N2/BS + ((N2 % BS) ? 1 : 0);
  T mysum = (T)0.0;
  for (int b1=0; b1 < NB1; b1++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*b1+i)*BS + tid < 3*N1)
        r1[0][i*BS + tid] = myR[3*e1_first + (3*b1+i)*BS + tid];
    __syncthreads();
    int ptcl1 = e1_first+b1*BS + tid;
    for (int b2=0; b2 < NB2; b2++)
    {
      // Load block of positions from global memory
      for (int i=0; i<3; i++)
        if ((3*b2+i)*BS + tid < 3*N2)
          r2[0][i*BS + tid] = myR[3*e2_first + (3*b2+i)*BS + tid];
      __syncthreads();
      // Now, loop over particles
      int end = (b2+1)*BS < N2 ? BS : N2-b2*BS;
      for (int j=0; j<end; j++)
      {
        int ptcl2 = e2_first + b2*BS+j;
        T dx, dy, dz;
        dx = r2[j][0] - r1[tid][0];
        dy = r2[j][1] - r1[tid][1];
        dz = r2[j][2] - r1[tid][2];
        T dist = min_dist(dx, dy, dz, L, Linv);
        if (ptcl1 != ptcl2 && (ptcl1 < (N1+e1_first) ) && (ptcl2 < (N2+e2_first)))
          mysum += eval_1d_spline (dist, rMax, drInv, A, coefs);
      }
      __syncthreads();
    }
  }
  __shared__ T shared_sum[BS];
  shared_sum[tid] = mysum;
  __syncthreads();
  for (int s=BS>>1; s>0; s >>=1)
  {
    if (tid < s)
      shared_sum[tid] += shared_sum[tid+s];
    __syncthreads();
  }
  T factor = (e1_first == e2_first) ? 0.5 : 1.0;
  if (tid==0)
    sum[blockIdx.x] += factor*shared_sum[0];
}

void
two_body_sum_PBC (float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                  float spline_coefs[], int numCoefs, float rMax,
                  float lattice[], float latticeInv[], float sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  two_body_sum_PBC_kernel<float,BS><<<dimGrid,dimBlock>>>
  (R, e1_first, e1_last, e2_first, e2_last,
   spline_coefs, numCoefs, rMax, lattice, latticeInv, sum);
}


void
two_body_sum_PBC (double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                  double spline_coefs[], int numCoefs, double rMax,
                  double lattice[], double latticeInv[], double sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  two_body_sum_PBC_kernel<double,BS><<<dimGrid,dimBlock>>>
  (R, e1_first, e1_last, e2_first, e2_last,
   spline_coefs, numCoefs, rMax, lattice, latticeInv, sum);
}




template<typename T, int BS>
__global__ void
two_body_ratio_PBC_kernel(T **R, int first, int last,
                          T *Rnew, int inew,
                          T *spline_coefs, int numCoefs, T rMax,
                          T *lattice, T* latticeInv, T* sum)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = ((T)1)/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  __shared__ T myRnew[3], myRold[3];
  if (tid == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
  {
    myRnew[tid] = Rnew[3*blockIdx.x+tid];
    myRold[tid] = myR[3*inew+tid];
  }
  __syncthreads();
  __shared__ T coefs[MAX_COEFS];
  __shared__ T r1[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[4][4];
  if (tid < 16)
    A[(tid>>2)][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ T shared_sum[BS];
  shared_sum[tid] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*N)
        r1[0][n] = myR[3*first + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz;
    dx = myRnew[0] - r1[tid][0];
    dy = myRnew[1] - r1[tid][1];
    dz = myRnew[2] - r1[tid][2];
    T dist = min_dist(dx, dy, dz, L, Linv);
    T delta = eval_1d_spline (dist, rMax, drInv, A, coefs);
    dx = myRold[0] - r1[tid][0];
    dy = myRold[1] - r1[tid][1];
    dz = myRold[2] - r1[tid][2];
    dist = min_dist(dx, dy, dz, L, Linv);
    delta -= eval_1d_spline (dist, rMax, drInv, A, coefs);
    if (ptcl1 != inew && (ptcl1 < (N+first) ))
      shared_sum[tid] += delta;
    __syncthreads();
  }
  __syncthreads();
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
      shared_sum[tid] += shared_sum[tid+s];
    __syncthreads();
  }
  if (tid==0)
    sum[blockIdx.x] += shared_sum[0];
}




void
two_body_ratio_PBC (float *R[], int first, int last,
                    float Rnew[], int inew,
                    float spline_coefs[], int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  two_body_ratio_PBC_kernel<float,BS><<<dimGrid,dimBlock>>>
  (R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, sum);
}



void
two_body_ratio_PBC (double *R[], int first, int last,
                    double Rnew[], int inew,
                    double spline_coefs[], int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  dim3 dimBlock(128);
  dim3 dimGrid(numWalkers);
  two_body_ratio_PBC_kernel<double,128><<<dimGrid,dimBlock>>>
  (R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, sum);
}

template<typename T, int BS>
__global__ void
two_body_ratio_grad_PBC_kernel(T const * const * __restrict__ R,
                               int first, int last,
                               T const * __restrict__ Rnew, int inew,
                               T const * __restrict__ spline_coefs,
                               int numCoefs, T rMax,
                               T const * __restrict__ lattice,
                               T const * __restrict__ latticeInv,
                               bool zero, T *__restrict__ ratio_grad)
{
  __shared__ T shared_grad[BS][4];
  __shared__ T coefs[MAX_COEFS];
  int tid = threadIdx.x;
  T dr = rMax /(T)(numCoefs-3);
  T drInv = ((T)1)/dr;
  // Safety for rounding error
  rMax *= (T)0.999999;
  if (tid < numCoefs)
  {
    coefs[tid] = spline_coefs[tid];
  }
  shared_grad[tid][0] = (T)0;
  shared_grad[tid][1] = (T)0;
  shared_grad[tid][2] = (T)0;
  shared_grad[tid][3] = (T)0;
  __syncthreads();
  T const * __restrict__ myR = R[blockIdx.x];
  T rnew_x = Rnew[3*blockIdx.x+0];
  T rnew_y = Rnew[3*blockIdx.x+1];
  T rnew_z = Rnew[3*blockIdx.x+2];
  T rold_x = myR[3*inew+0];
  T rold_y = myR[3*inew+1];
  T rold_z = myR[3*inew+2];
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    T r1[3];
    int n = b*BS + tid;
    if ( n < N )
    {
      r1[0] = myR[3*(first+n)  ];
      r1[1] = myR[3*(first+n)+1];
      r1[2] = myR[3*(first+n)+2];
    }
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz, u, du, d2u, delta, dist;
    dx = rold_x - r1[0];
    dy = rold_y - r1[1];
    dz = rold_z - r1[2];
    dist = CMC_min_dist_only(dx, dy, dz/*, L, Linv, images*/);
    delta = -CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs);
    dx = rnew_x - r1[0];
    dy = rnew_y - r1[1];
    dz = rnew_z - r1[2];
    dist = CMC_min_dist(dx, dy, dz/*, L, Linv, images*/);
    CMC_eval_1d_spline_vgl (dist, rMax, drInv/*, A*/, coefs, u, du, d2u);
    delta += u;
    if ((ptcl1 != inew) && (ptcl1 < (N + first) ))
    {
      du /= dist;
      shared_grad[tid][0] += delta;
      shared_grad[tid][1] += du * dx;
      shared_grad[tid][2] += du * dy;
      shared_grad[tid][3] += du * dz;
    }
  }
  __syncthreads();
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
    {
      shared_grad[tid][0] += shared_grad[tid+s][0];
      shared_grad[tid][1] += shared_grad[tid+s][1];
      shared_grad[tid][2] += shared_grad[tid+s][2];
      shared_grad[tid][3] += shared_grad[tid+s][3];
    }
    __syncthreads();
  }
  if (zero)
  {
    if ( tid < 4 ) ratio_grad[4*blockIdx.x+tid] = shared_grad[0][tid];
  }
  else
  {
    if ( tid < 4 ) ratio_grad[4*blockIdx.x+tid] += shared_grad[0][tid];
  }
}



template<typename T, int BS>
__global__ void
two_body_ratio_grad_PBC_kernel_fast (T **R, int first, int last,
                                     T *Rnew, int inew,
                                     T *spline_coefs, int numCoefs, T rMax,
                                     T *lattice, T* latticeInv,
                                     bool zero, T* ratio_grad)
{
  int tid = threadIdx.x;
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  __shared__ T *myR;
  __shared__ T myRnew[3], myRold[3];
  if (tid == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
  {
    myRnew[tid] = Rnew[3*blockIdx.x+tid];
    myRold[tid] = myR[3*inew+tid];
  }
  __syncthreads();
  __shared__ T coefs[MAX_COEFS];
  __shared__ T r1[BS][3];
  /*
  __shared__ T L[3][3], Linv[3][3];
  */
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  /*
  if (tid < 9) {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  */
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  int N = last - first + 1;
  int NB = (N+BS-1)/BS;
  __shared__ T shared_sum[BS];
  __shared__ T shared_grad[BS][3];
  shared_sum[tid] = (T)0.0;
  shared_grad[tid][0] = shared_grad[tid][1] = shared_grad[tid][2] = 0.0f;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*N)
        r1[0][n] = myR[3*first + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz, u, du, d2u, delta, dist;
    dx = myRold[0] - r1[tid][0];
    dy = myRold[1] - r1[tid][1];
    dz = myRold[2] - r1[tid][2];
    dist = CMC_min_dist_fast(dx, dy, dz/*, L, Linv*/);
    delta = -eval_1d_spline (dist, rMax, drInv, A, coefs);
    dx = myRnew[0] - r1[tid][0];
    dy = myRnew[1] - r1[tid][1];
    dz = myRnew[2] - r1[tid][2];
    dist = CMC_min_dist_fast(dx, dy, dz/*, L, Linv*/);
    eval_1d_spline_vgl (dist, rMax, drInv, A, coefs,
                        u, du, d2u);
    delta += u;
    if (ptcl1 != inew && (ptcl1 < (N+first) ))
    {
      du /= dist;
      shared_sum[tid] += delta;
      shared_grad[tid][0] += du * dx;
      shared_grad[tid][1] += du * dy;
      shared_grad[tid][2] += du * dz;
    }
    __syncthreads();
  }
  __syncthreads();
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
    {
      shared_sum[tid] += shared_sum[tid+s];
      shared_grad[tid][0] += shared_grad[tid+s][0];
      shared_grad[tid][1] += shared_grad[tid+s][1];
      shared_grad[tid][2] += shared_grad[tid+s][2];
    }
    __syncthreads();
  }
  if (tid==0)
  {
    if (zero)
    {
      ratio_grad[4*blockIdx.x+0] = shared_sum[0];
      ratio_grad[4*blockIdx.x+1] = shared_grad[0][0];
      ratio_grad[4*blockIdx.x+2] = shared_grad[0][1];
      ratio_grad[4*blockIdx.x+3] = shared_grad[0][2];
    }
    else
    {
      ratio_grad[4*blockIdx.x+0] += shared_sum[0];
      ratio_grad[4*blockIdx.x+1] += shared_grad[0][0];
      ratio_grad[4*blockIdx.x+2] += shared_grad[0][1];
      ratio_grad[4*blockIdx.x+3] += shared_grad[0][2];
    }
  }
}




// use_fast_image indicates that Rmax < simulation cell radius.  In
// this case, we don't have to search over 27 images.
void
two_body_ratio_grad_PBC(float *R[], int first, int last,
                        float  Rnew[], int inew,
                        float spline_coefs[], int numCoefs, float rMax,
                        float lattice[], float latticeInv[], bool zero,
                        float ratio_grad[], int numWalkers,
                        bool use_fast_image)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  CMC_PROFILING_BEGIN();
  // fprintf(stderr, "first = %d\n", first);
  // fprintf(stderr, "last  = %d\n", last);
  // fprintf(stderr, "inew  = %d\n", inew);
  // fprintf(stderr, "rMax = %1.3f\n", rMax);
  if (use_fast_image)
  {
    two_body_ratio_grad_PBC_kernel_fast<float,BS><<<dimGrid,dimBlock>>>
    (R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
     lattice, latticeInv, zero, ratio_grad);
  }
  else
  {
    two_body_ratio_grad_PBC_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
    (R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
     lattice, latticeInv, zero, ratio_grad);
  }
  CMC_PROFILING_END(__LINE__);
}


void
two_body_ratio_grad_PBC(double *R[], int first, int last,
                        double  Rnew[], int inew,
                        double spline_coefs[], int numCoefs, double rMax,
                        double lattice[], double latticeInv[], bool zero,
                        double ratio_grad[], int numWalkers,
                        bool use_fast_image)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  if (use_fast_image)
  {
    two_body_ratio_grad_PBC_kernel_fast<double,BS><<<dimGrid,dimBlock>>>
    (R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
     lattice, latticeInv, zero, ratio_grad);
  }
  else
  {
    two_body_ratio_grad_PBC_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
    (R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
     lattice, latticeInv, zero, ratio_grad);
  }
  //CMC_PROFILING_END(__LINE__);
}




template<typename T, int BS>
__global__ void
two_body_NLratio_PBC_kernel(NLjobGPU<T> const * __restrict__ jobs,
                            int first, int last,
                            T const * const * __restrict__ spline_coefs,
                            int const * __restrict__ numCoefs,
                            T const * __restrict__ rMaxList,
                            T const * __restrict__ lattice,
                            T const * __restrict__ latticeInv,
                            T sim_cell_radius)
{
  const int MAX_RATIOS = 18;
  __shared__ T shared_sum[MAX_RATIOS][BS+1];
  __shared__ T myRnew[MAX_RATIOS][3];
  __shared__ T coefs[MAX_COEFS];
  __shared__ T r1[BS][3];
  T const * __restrict__ myCoefs = spline_coefs[blockIdx.x];
  NLjobGPU<T> myJob = jobs[blockIdx.x];
  const int myNumCoefs = numCoefs[blockIdx.x];
  const int tid = threadIdx.x;
  if (tid < myNumCoefs)
  {
    coefs[tid] = myCoefs[tid];
  }
  for (int i = 0; i < 3; i++)
  {
    if (i*BS + tid < 3*myJob.NumQuadPoints)
    {
      myRnew[0][i*BS+tid] = myJob.QuadPoints[i*BS+tid];
    }
  }
  for (int i = 0; i < myJob.NumQuadPoints; i++)
  {
    shared_sum[i][tid] = (T)0;
  }
  __syncthreads();
  const T rMax = rMaxList[blockIdx.x];
  const T dr = rMax / (myNumCoefs - 3);
  const T drInv = (T)1.0 / dr;
  const int use_fast = sim_cell_radius >= rMax;
  const T rold_x = myJob.R[3*myJob.Elec+0];
  const T rold_y = myJob.R[3*myJob.Elec+1];
  const T rold_z = myJob.R[3*myJob.Elec+2];
  const int N = last - first + 1;
  const int NB = N / BS + ((N % BS) ? 1 : 0);
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if (((3*b+i)*BS + tid) < (3*N))
      {
        r1[0][n] = myJob.R[3*first + (3*b+i)*BS + tid];
      }
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx = rold_x - r1[tid][0];
    T dy = rold_y - r1[tid][1];
    T dz = rold_z - r1[tid][2];
    T dist;
    if (use_fast)
    {
      dist = CMC_min_dist_fast(dx, dy, dz/*, L, Linv*/);
    }
    else
    {
      dist = CMC_min_dist_only(dx, dy, dz/*, L, Linv, images*/);
    }
    T uOld = CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs);
    if (use_fast)
    {
      for (int iq=0; iq<myJob.NumQuadPoints; iq++)
      {
        dx = myRnew[iq][0] - r1[tid][0];
        dy = myRnew[iq][1] - r1[tid][1];
        dz = myRnew[iq][2] - r1[tid][2];
        dist = CMC_min_dist_fast(dx, dy, dz/*, L, Linv*/);
        if ((ptcl1 != myJob.Elec) && (ptcl1 < (N + first)))
        {
          shared_sum[iq][tid] += CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs) - uOld;
        }
      }
    }
    else
    {
      for (int iq=0; iq<myJob.NumQuadPoints; iq++)
      {
        dx = myRnew[iq][0] - r1[tid][0];
        dy = myRnew[iq][1] - r1[tid][1];
        dz = myRnew[iq][2] - r1[tid][2];
        dist = CMC_min_dist_only(dx, dy, dz/*, L, Linv, images*/);
        if ((ptcl1 != myJob.Elec) && (ptcl1 < (N + first)))
        {
          shared_sum[iq][tid] += CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs) - uOld;
        }
      }
    }
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
    {
      for (int iq=0; iq < myJob.NumQuadPoints; iq++)
      {
        shared_sum[iq][tid] += shared_sum[iq][tid+s];
      }
    }
    __syncthreads();
  }
  if (tid < myJob.NumQuadPoints)
  {
    myJob.Ratios[tid] *= exp(-shared_sum[tid][0]);
  }
}


/* no longer needed, use the template version instead
template<int BS>
__global__ void
two_body_NLratio_PBC_kernel(NLjobGPU<double> *jobs, int first, int last,
                            double **spline_coefs, int *numCoefs,
                            double *rMaxList,
                            double *lattice, double *latticeInv,
                            double sim_cell_radius)
{
  const int MAX_RATIOS = 18;
  int tid = threadIdx.x;
  __shared__ NLjobGPU<double> myJob;
  __shared__ double myRnew[MAX_RATIOS][3], myRold[3];
  __shared__ double* myCoefs;
  __shared__ int myNumCoefs;
  __shared__ double rMax;
  if (tid == 0)
  {
    myJob = jobs[blockIdx.x];
    myCoefs = spline_coefs[blockIdx.x];
    myNumCoefs = numCoefs[blockIdx.x];
    rMax = rMaxList[blockIdx.x];
  }
  __syncthreads();
  if (tid < 3 )
    myRold[tid] = myJob.R[3*myJob.Elec+tid];
  for (int i=0; i<3; i++)
    if (i*BS + tid < 3*myJob.NumQuadPoints)
      myRnew[0][i*BS+tid] = myJob.QuadPoints[i*BS+tid];
  __syncthreads();
  double dr = rMax/(double)(myNumCoefs-3);
  double drInv = 1.0/dr;
  __shared__ double coefs[MAX_COEFS];
  __shared__ double r1[BS][3];
  __shared__ double L[3][3], Linv[3][3];
  if (tid < myNumCoefs)
    coefs[tid] = myCoefs[tid];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ double A[4][4];
  if (tid < 16)
    A[(tid>>2)][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ double shared_sum[MAX_RATIOS][BS+1];
  for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    shared_sum[iq][tid] = (double)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*N)
        r1[0][n] = myJob.R[3*first + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    double dx, dy, dz;
    dx = myRold[0] - r1[tid][0];
    dy = myRold[1] - r1[tid][1];
    dz = myRold[2] - r1[tid][2];
    double dist = min_dist(dx, dy, dz, L, Linv);
    double uOld = eval_1d_spline (dist, rMax, drInv, A, coefs);
    for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    {
      dx = myRnew[iq][0] - r1[tid][0];
      dy = myRnew[iq][1] - r1[tid][1];
      dz = myRnew[iq][2] - r1[tid][2];
      dist = min_dist(dx, dy, dz, L, Linv);
      if (ptcl1 != myJob.Elec && (ptcl1 < (N+first)))
        shared_sum[iq][tid] += eval_1d_spline (dist, rMax, drInv, A, coefs) - uOld;
    }
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
      for (int iq=0; iq < myJob.NumQuadPoints; iq++)
        shared_sum[iq][tid] += shared_sum[iq][tid+s];
    __syncthreads();
  }
  if (tid < myJob.NumQuadPoints)
    myJob.Ratios[tid] *= exp(-shared_sum[tid][0]);
}

*/


void
two_body_NLratios_PBC(NLjobGPU<float> jobs[], int first, int last,
                      float* spline_coefs[], int numCoefs[], float rMax[],
                      float lattice[], float latticeInv[], float sim_cell_radius,
                      int numjobs)
{
  if (numjobs==0) return;
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=32;
  dim3 dimBlock(BS);
  CMC_PROFILING_BEGIN();
  while (numjobs > 65535)
  {
    dim3 dimGrid(65535);
    two_body_NLratio_PBC_kernel<float, BS><<<dimGrid,dimBlock>>>
    (jobs, first, last, spline_coefs, numCoefs, rMax,
     lattice, latticeInv, sim_cell_radius);
    jobs += 65535;
    numjobs -= 65535;
  }
  dim3 dimGrid(numjobs);
  //fprintf(stdout, " numjobs %d\n", numjobs);
  two_body_NLratio_PBC_kernel<float, BS><<<dimGrid,dimBlock>>>
  (jobs, first, last, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, sim_cell_radius);
  CMC_PROFILING_END(__LINE__);
}

void
two_body_NLratios_PBC(NLjobGPU<double> jobs[], int first, int last,
                      double* spline_coefs[], int numCoefs[], double rMax[],
                      double lattice[], double latticeInv[],
                      double sim_cell_radius, int numjobs)
{
  if (numjobs==0) return;
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=32;
  dim3 dimBlock(BS);
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  while (numjobs > 65535)
  {
    dim3 dimGrid(65535);
    two_body_NLratio_PBC_kernel<double, BS><<<dimGrid,dimBlock>>>
    (jobs, first, last, spline_coefs, numCoefs, rMax,
     lattice, latticeInv, sim_cell_radius);
    jobs += 65535;
    numjobs -= 65535;
  }
  dim3 dimGrid(numjobs);
  two_body_NLratio_PBC_kernel<double, BS><<<dimGrid,dimBlock>>>
  (jobs, first, last, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, sim_cell_radius);
  //CMC_PROFILING_END(__LINE__);
}



template<typename T>
__global__ void
two_body_update_PBC_kernel (T **R, int N, int iat)
{
  __shared__ T* myR;
  if (threadIdx.x == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (threadIdx.x < 3)
    myR[3*iat + threadIdx.x] = myR[3*N + threadIdx.x];
}


void
two_body_update_PBC(float *R[], int N, int iat, int numWalkers)
{
  dim3 dimBlock(32);
  dim3 dimGrid(numWalkers);
  two_body_update_PBC_kernel<float><<<dimGrid, dimBlock>>> (R, N, iat);
}

void
two_body_update_PBC(double *R[], int N, int iat, int numWalkers)
{
  dim3 dimBlock(3);
  dim3 dimGrid(numWalkers);
  two_body_update_PBC_kernel<double><<<dimGrid, dimBlock>>> (R, N, iat);
}

#define MAX_COEFS 32

template<typename T>
__global__ void
two_body_grad_lapl_PBC_kernel(T **R, int e1_first, int e1_last,
                              int e2_first, int e2_last,
                              T *spline_coefs, int numCoefs, T rMax,
                              T *lattice, T *latticeInv,
                              T *gradLapl, int row_stride)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  if (tid == 0)
    myR = R[blockIdx.x];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  /*
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9) {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }


  __shared__ T A[12][4];
  if (tid < 16) {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  */
  __syncthreads();
  int N1 = e1_last - e1_first + 1;
  int N2 = e2_last - e2_first + 1;
  T sGradLapl[4];
  T r1[3], r2[3];
  int b1 = blockIdx.y;
  int BS1 = blockDim.x;
  // Load positions from global memory
  int ptcl1 = e1_first + b1*BS1 + tid;
  if (b1*BS1 + tid < N1)
  {
    r1[0] = myR[3*ptcl1    ];
    r1[1] = myR[3*ptcl1 + 1];
    r1[2] = myR[3*ptcl1 + 2];
  }
  int offset = blockIdx.x * row_stride + 4 * ( b1*BS1 + e1_first );
  sGradLapl[0] = sGradLapl[1] = sGradLapl[2] = sGradLapl[3] = (T)0.0;
  // Now, loop over particles
  for (int j=0; j<N2; j++)
  {
    int ptcl2 = e2_first + j;
    r2[0] = myR[3*ptcl2    ];
    r2[1] = myR[3*ptcl2 + 1];
    r2[2] = myR[3*ptcl2 + 2];
    T dx, dy, dz, u, du, d2u;
    dx = r2[0] - r1[0];
    dy = r2[1] - r1[1];
    dz = r2[2] - r1[2];
    T dist = CMC_min_dist(dx, dy, dz/*, L, Linv*/);
    CMC_eval_1d_spline_vgl (dist, rMax, drInv/*, A*/, coefs, u, du, d2u);
    //if (ptcl1 != ptcl2 && ptcl1 < (N1+e1_first))
    if (ptcl1 != ptcl2 && ptcl1 < (N1+e1_first))
    {
      du /= dist;
      sGradLapl[0] += du * dx;
      sGradLapl[1] += du * dy;
      sGradLapl[2] += du * dz;
      sGradLapl[3] -= d2u + 2.0*du;
    }
  }
  if (b1*BS1 + tid < N1)
  {
    gradLapl[offset + tid*4    ] += sGradLapl[0];
    gradLapl[offset + tid*4 + 1] += sGradLapl[1];
    gradLapl[offset + tid*4 + 2] += sGradLapl[2];
    gradLapl[offset + tid*4 + 3] += sGradLapl[3];
  }
}


template<typename T, int BS>
__global__ void
two_body_grad_lapl_PBC_kernel_fast(T **R, int e1_first, int e1_last,
                                   int e2_first, int e2_last,
                                   T *spline_coefs, int numCoefs, T rMax,
                                   T *lattice, T *latticeInv,
                                   T *gradLapl, int row_stride)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  if (tid == 0)
    myR = R[blockIdx.x];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  __shared__ T r1[BS][3], r2[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  int N1 = e1_last - e1_first + 1;
  int N2 = e2_last - e2_first + 1;
  int NB1 = N1/BS + ((N1 % BS) ? 1 : 0);
  int NB2 = N2/BS + ((N2 % BS) ? 1 : 0);
  __shared__ T sGradLapl[BS][4];
  for (int b1=0; b1 < NB1; b1++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*b1+i)*BS + tid < 3*N1)
        r1[0][i*BS + tid] = myR[3*e1_first + (3*b1+i)*BS + tid];
    __syncthreads();
    int ptcl1 = e1_first+b1*BS + tid;
    int offset = blockIdx.x * row_stride + 4*b1*BS + 4*e1_first;
    sGradLapl[tid][0] = sGradLapl[tid][1] =
                          sGradLapl[tid][2] = sGradLapl[tid][3] = (T)0.0;
    for (int b2=0; b2 < NB2; b2++)
    {
      // Load block of positions from global memory
      for (int i=0; i<3; i++)
        if ((3*b2+i)*BS + tid < 3*N2)
          r2[0][i*BS + tid] = myR[3*e2_first + (3*b2+i)*BS + tid];
      __syncthreads();
      // Now, loop over particles
      int end = (b2+1)*BS < N2 ? BS : N2-b2*BS;
      for (int j=0; j<end; j++)
      {
        int ptcl2 = e2_first + b2*BS+j;
        T dx, dy, dz, u, du, d2u;
        dx = r2[j][0] - r1[tid][0];
        dy = r2[j][1] - r1[tid][1];
        dz = r2[j][2] - r1[tid][2];
        T dist = min_dist(dx, dy, dz, L, Linv);
        eval_1d_spline_vgl (dist, rMax, drInv, A, coefs, u, du, d2u);
        if (ptcl1 != ptcl2 && (ptcl1 < (N1+e1_first) ) && (ptcl2 < (N2+e2_first)))
        {
          du /= dist;
          sGradLapl[tid][0] += du * dx;
          sGradLapl[tid][1] += du * dy;
          sGradLapl[tid][2] += du * dz;
          sGradLapl[tid][3] -= d2u + 2.0*du;
        }
      }
      __syncthreads();
    }
    for (int i=0; i<4; i++)
      if ((4*b1+i)*BS + tid < 4*N1)
        gradLapl[offset + i*BS +tid] += sGradLapl[0][i*BS+tid];
    __syncthreads();
  }
}



void
two_body_grad_lapl_PBC(float *R[], int e1_first, int e1_last,
                       int e2_first, int e2_last,
                       float spline_coefs[], int numCoefs, float rMax,
                       float lattice[], float latticeInv[], float sim_cell_radius,
                       float gradLapl[], int row_stride, int numWalkers)
{
  CMC_PROFILING_BEGIN();
  if (sim_cell_radius >= rMax)
  {
    const int BS=32;
    dim3 dimBlock(BS);
    dim3 dimGrid(numWalkers);
    two_body_grad_lapl_PBC_kernel_fast<float,BS><<<dimGrid,dimBlock>>>
    (R, e1_first, e1_last, e2_first, e2_last, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  gradLapl, row_stride);
  }
  else
  {
    const int BS1=128;
    dim3 dimBlock(BS1);
    dim3 dimGrid(numWalkers, (e1_last - e1_first + BS1)/BS1);
    two_body_grad_lapl_PBC_kernel<float><<<dimGrid,dimBlock>>>
    (R, e1_first, e1_last, e2_first, e2_last, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  gradLapl, row_stride);
  }
  CMC_PROFILING_END(__LINE__);
}


void
two_body_grad_lapl_PBC(double *R[], int e1_first, int e1_last,
                       int e2_first, int e2_last,
                       double spline_coefs[], int numCoefs, double rMax,
                       double lattice[], double latticeInv[], double sim_cell_radius,
                       double gradLapl[], int row_stride, int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  if (sim_cell_radius >= rMax)
  {
    const int BS=32;
    dim3 dimBlock(BS);
    dim3 dimGrid(numWalkers);
    two_body_grad_lapl_PBC_kernel_fast<double,BS><<<dimGrid,dimBlock>>>
    (R, e1_first, e1_last, e2_first, e2_last, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  gradLapl, row_stride);
  }
  else
  {
    const int BS1=128;
    dim3 dimBlock(BS1);
    dim3 dimGrid(numWalkers, (e1_last - e1_first + BS1)/BS1);
    two_body_grad_lapl_PBC_kernel<double><<<dimGrid,dimBlock>>>
    (R, e1_first, e1_last, e2_first, e2_last, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  gradLapl, row_stride);
  }
  //CMC_PROFILING_END(__LINE__);
}



template<typename T, int BS>
__global__ void
two_body_grad_PBC_kernel (T const * const * __restrict__ R,
                          int first, int last, int iat,
                          T const * __restrict__ spline_coefs,
                          int numCoefs, T rMax,
                          T const * __restrict__ lattice,
                          T const * __restrict__ latticeInv,
                          bool zeroOut, T * __restrict__ grad)
{
  __shared__ T sGrad[BS][3];
  __shared__ T coefs[MAX_COEFS];
  T dr = rMax/(T)(numCoefs-3);
  T drInv = ((T)1)/dr;
  int tid = threadIdx.x;
  // Safety for rounding error
  rMax *= (T)0.999999;
  if (tid < numCoefs)
  {
    coefs[tid] = spline_coefs[tid];
  }
  sGrad[tid][0] = (T)0;
  sGrad[tid][1] = (T)0;
  sGrad[tid][2] = (T)0;
  __syncthreads();
  T const * __restrict__ myR = R[blockIdx.x];
  T r2_x = myR[3*iat+0];
  T r2_y = myR[3*iat+1];
  T r2_z = myR[3*iat+2];
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  for (int b = 0; b < NB; b++)
  {
    // Load block of positions from global memory
    int n = b*BS + tid;
    T r1[3];
    if ( n < N )
    {
      r1[0] = myR[3*(first+n)  ];
      r1[1] = myR[3*(first+n)+1];
      r1[2] = myR[3*(first+n)+2];
    }
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz, u, du, d2u;
    dx = r2_x - r1[0];
    dy = r2_y - r1[1];
    dz = r2_z - r1[2];
    T dist = CMC_min_dist(dx, dy, dz/*, L, Linv, images*/);
    CMC_eval_1d_spline_vgl (dist, rMax, drInv/*, A*/, coefs, u, du, d2u);
    if (ptcl1 != iat && ptcl1 < (N+first))
    {
      du /= dist;
      sGrad[tid][0] += du * dx;
      sGrad[tid][1] += du * dy;
      sGrad[tid][2] += du * dz;
    }
  }
  __syncthreads();
  // Do reduction across threads in block
  for (int s=BS>>1; s>0; s>>=1)
  {
    if (tid < s)
    {
      sGrad[tid][0] += sGrad[tid+s][0];
      sGrad[tid][1] += sGrad[tid+s][1];
      sGrad[tid][2] += sGrad[tid+s][2];
    }
    __syncthreads();
  }
  if (tid < 3)
  {
    if (zeroOut)
    {
      grad[3*blockIdx.x + tid]  = sGrad[0][tid];
    }
    else
    {
      grad[3*blockIdx.x + tid] += sGrad[0][tid];
    }
  }
}


template<typename T, int BS>
__global__ void
two_body_grad_PBC_kernel_fast(T **R, int first, int last, int iat,
                              T *spline_coefs, int numCoefs, T rMax,
                              T *lattice, T *latticeInv, bool zeroOut, T *grad)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR, r2[3];
  if (tid == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (tid < 3)
    r2[tid] = myR[3*iat+tid];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  __shared__ T r1[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ T sGrad[BS][3];
  sGrad[tid][0]   = sGrad[tid][1] = sGrad[tid][2] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*b+i)*BS + tid < 3*N)
        r1[0][i*BS + tid] = myR[3*first + (3*b+i)*BS + tid];
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz, u, du, d2u;
    dx = r2[0] - r1[tid][0];
    dy = r2[1] - r1[tid][1];
    dz = r2[2] - r1[tid][2];
    T dist = min_dist_fast(dx, dy, dz, L, Linv);
    eval_1d_spline_vgl (dist, rMax, drInv, A, coefs, u, du, d2u);
    if (ptcl1 != iat && ptcl1 < (N+first))
    {
      du /= dist;
      sGrad[tid][0] += du * dx;
      sGrad[tid][1] += du * dy;
      sGrad[tid][2] += du * dz;
    }
    __syncthreads();
  }
  // Do reduction across threads in block
  for (int s=BS>>1; s>0; s>>=1)
  {
    if (tid < s)
    {
      sGrad[tid][0] += sGrad[tid+s][0];
      sGrad[tid][1] += sGrad[tid+s][1];
      sGrad[tid][2] += sGrad[tid+s][2];
    }
    __syncthreads();
  }
  if (tid < 3)
  {
    if (zeroOut)
      grad[3*blockIdx.x + tid] = sGrad[0][tid];
    else
      grad[3*blockIdx.x + tid] += sGrad[0][tid];
  }
}





void
two_body_gradient_PBC (float *R[], int first, int last, int iat,
                       float spline_coefs[], int numCoefs, float rMax,
                       float lattice[], float latticeInv[], float sim_cell_radius,
                       bool zeroOut,
                       float grad[], int numWalkers)
{
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  CMC_PROFILING_BEGIN();
  if (sim_cell_radius >= rMax)
    two_body_grad_PBC_kernel_fast<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
    (R, first, last, iat, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  zeroOut, grad);
  else
    two_body_grad_PBC_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
    (R, first, last, iat, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  zeroOut, grad);
  CMC_PROFILING_END(__LINE__);
}


void
two_body_gradient_PBC (double *R[], int first, int last, int iat,
                       double spline_coefs[], int numCoefs, double rMax,
                       double lattice[], double latticeInv[],
                       double sim_cell_radius, bool zeroOut,
                       double grad[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  if (sim_cell_radius >= rMax)
    two_body_grad_PBC_kernel_fast<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
    (R, first, last, iat, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  zeroOut, grad);
  else
    two_body_grad_PBC_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
    (R, first, last, iat, spline_coefs, numCoefs,
     rMax, lattice, latticeInv,  zeroOut, grad);
  //CMC_PROFILING_END(__LINE__);
}




template<typename T, int BS, unsigned COMPLEX>
__global__ void
two_body_derivs_PBC_kernel(T **R, T **gradLogPsi,
                           int e1_first, int e1_last,
                           int e2_first, int e2_last,
                           int numCoefs, T rMax,
                           T *lattice, T *latticeInv,
                           T **derivs)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0f/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR, *myGrad, *myDerivs;
  if (tid == 0)
  {
    myR      =          R[blockIdx.x];
    myGrad   = gradLogPsi[blockIdx.x];
    myDerivs =     derivs[blockIdx.x];
  }
  __shared__ T sderivs[MAX_COEFS][2];
  // __shared__ T coefs[MAX_COEFS];
  // if (tid < numCoefs)
  //   coefs[tid] = spline_coefs[tid];
  __shared__ T r1[BS][3], r2[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  sderivs[tid][0] = T();
  sderivs[tid][1] = T();
  int N1 = e1_last - e1_first + 1;
  int N2 = e2_last - e2_first + 1;
  int NB1 = N1/BS + ((N1 % BS) ? 1 : 0);
  int NB2 = N2/BS + ((N2 % BS) ? 1 : 0);
  __shared__ T sGrad[BS][3];
  for (int b1=0; b1 < NB1; b1++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*b1+i)*BS + tid < 3*N1)
      {
        int outoff = i*BS+tid;
        int inoff  = outoff + 3*e1_first + 3*b1*BS;
        r1[0][outoff]    =    myR[inoff];//[3*e1_first + (3*b1+i)*BS + tid];
        sGrad[0][outoff] = myGrad[inoff*COMPLEX];
      }
    __syncthreads();
    int ptcl1 = e1_first+b1*BS + tid;
    for (int b2=0; b2 < NB2; b2++)
    {
      // Load block of positions from global memory
      for (int i=0; i<3; i++)
        if ((3*b2+i)*BS + tid < 3*N2)
          r2[0][i*BS + tid] = myR[3*e2_first + (3*b2+i)*BS + tid];
      __syncthreads();
      // Now, loop over particles
      int end = (b2+1)*BS < N2 ? BS : N2-b2*BS;
      for (int j=0; j<end; j++)
      {
        int ptcl2 = e2_first + b2*BS+j;
        T dx, dy, dz;
        dx = r2[j][0] - r1[tid][0];
        dy = r2[j][1] - r1[tid][1];
        dz = r2[j][2] - r1[tid][2];
        T dist = min_dist(dx, dy, dz, L, Linv);
        T distInv = 1.0f/dist;
        T s = dist * drInv;
        T sf = floorf (s);
        int index = (int)sf;
        T t = s - sf;
        T t2 = t*t;
        T t3 = t*t2;
        T v0, v1, v2, v3;
        // sderivs[index+0][0] += (A[0][0]*t3 + A[0][1]*t2 + A[0][2]*t + A[0][3]);
        // sderivs[index+1][0] += (A[1][0]*t3 + A[1][1]*t2 + A[1][2]*t + A[1][3]);
        // sderivs[index+2][0] += (A[2][0]*t3 + A[2][1]*t2 + A[2][2]*t + A[2][3]);
        // sderivs[index+3][0] += (A[3][0]*t3 + A[3][1]*t2 + A[3][2]*t + A[3][3]);
        v0 = (A[0][0]*t3 + A[0][1]*t2 + A[0][2]*t + A[0][3]);
        v1 = (A[1][0]*t3 + A[1][1]*t2 + A[1][2]*t + A[1][3]);
        v2 = (A[2][0]*t3 + A[2][1]*t2 + A[2][2]*t + A[2][3]);
        v3 = (A[3][0]*t3 + A[3][1]*t2 + A[3][2]*t + A[3][3]);
        for (int id=0; id<BS; id++)
          if (tid == id && ptcl1 != ptcl2 && ptcl1 <= e1_last && (dist < rMax))
          {
            sderivs[index+0][0] += v0;
            sderivs[index+1][0] += v1;
            sderivs[index+2][0] += v2;
            sderivs[index+3][0] += v3;
          }
        T prefact = (dx*sGrad[tid][0] + dy*sGrad[tid][1] + dz*sGrad[tid][2])*distInv;
        T du0 = drInv * (A[4][0]*t3 + A[4][1]*t2 + A[4][2]*t + A[4][3]);
        T du1 = drInv * (A[5][0]*t3 + A[5][1]*t2 + A[5][2]*t + A[5][3]);
        T du2 = drInv * (A[6][0]*t3 + A[6][1]*t2 + A[6][2]*t + A[6][3]);
        T du3 = drInv * (A[7][0]*t3 + A[7][1]*t2 + A[7][2]*t + A[7][3]);
        // This is the dot (gradu, grad_log_psi) term.
        v0 = 2.0f* prefact * du0;
        v1 = 2.0f* prefact * du1;
        v2 = 2.0f* prefact * du2;
        v3 = 2.0f* prefact * du3;
        // This is the lapl u term
        v0 -= drInv*drInv*(A[ 8][0]*t3 + A[ 8][1]*t2 + A[ 8][2]*t + A[ 8][3]) + 2.0f*du0*distInv;
        v1 -= drInv*drInv*(A[ 9][0]*t3 + A[ 9][1]*t2 + A[ 9][2]*t + A[ 9][3]) + 2.0f*du1*distInv;
        v2 -= drInv*drInv*(A[10][0]*t3 + A[10][1]*t2 + A[10][2]*t + A[10][3]) + 2.0f*du2*distInv;
        v3 -= drInv*drInv*(A[11][0]*t3 + A[11][1]*t2 + A[11][2]*t + A[11][3]) + 2.0f*du3*distInv;
        for (int id=0; id<BS; id++)
          if (tid == id && ptcl1 != ptcl2 && ptcl1 <= e1_last && (dist < rMax))
          {
            sderivs[index+0][1] += v0;
            sderivs[index+1][1] += v1;
            sderivs[index+2][1] += v2;
            sderivs[index+3][1] += v3;
          }
      }
      __syncthreads();
    }
  }
  //  if (e1_first == e2_first)
  sderivs[tid][0] *= 0.5f;
  sderivs[tid][1] *= 0.5f;
  if (tid < 2*numCoefs)
    myDerivs[tid] = -sderivs[0][tid];
  if (tid+BS < 2*numCoefs)
    myDerivs[tid+BS] = sderivs[0][tid+BS];
}


void
two_body_derivs_PBC(float *R[], float *gradLogPsi[], int e1_first, int e1_last,
                    int e2_first, int e2_last,
                    int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sim_cell_radius,
                    float *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  if (sim_cell_radius >= rMax)
    two_body_derivs_PBC_kernel<float,BS,1><<<dimGrid,dimBlock>>>
    (R, gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    two_body_derivs_PBC_kernel<float,BS,1><<<dimGrid,dimBlock>>>
    (R, gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
}

void
two_body_derivs_PBC(double *R[], double *gradLogPsi[], int e1_first, int e1_last,
                    int e2_first, int e2_last,
                    int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sim_cell_radius,
                    double *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  if (sim_cell_radius >= rMax)
    two_body_derivs_PBC_kernel<double,BS,1><<<dimGrid,dimBlock>>>
    (R, gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    two_body_derivs_PBC_kernel<double,BS,1><<<dimGrid,dimBlock>>>
    (R, gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
}


// Ye: use offset to recycle the old routines
// block size can be further optimized.
#ifdef QMC_COMPLEX
void
two_body_derivs_PBC(float *R[], std::complex<float> *gradLogPsi[], int e1_first, int e1_last,
                    int e2_first, int e2_last,
                    int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sim_cell_radius,
                    float *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);

  if (sim_cell_radius >= rMax)
    two_body_derivs_PBC_kernel<float,BS,2><<<dimGrid,dimBlock>>>
    (R, (float**)gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    two_body_derivs_PBC_kernel<float,BS,2><<<dimGrid,dimBlock>>>
    (R, (float**)gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
}

void
two_body_derivs_PBC(double *R[], std::complex<double> *gradLogPsi[], int e1_first, int e1_last,
                    int e2_first, int e2_last,
                    int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sim_cell_radius,
                    double *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);

  if (sim_cell_radius >= rMax)
    two_body_derivs_PBC_kernel<double,BS,2><<<dimGrid,dimBlock>>>
    (R, (double**)gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    two_body_derivs_PBC_kernel<double,BS,2><<<dimGrid,dimBlock>>>
    (R, (double**)gradLogPsi, e1_first, e1_last, e2_first, e2_last, numCoefs,
     rMax, lattice, latticeInv, derivs);
}

#endif

////////////////////////////////////////////////////////////////
//                      One-body routines                     //
////////////////////////////////////////////////////////////////

template<typename T, int BS >
__global__ void
one_body_sum_PBC_kernel(T *C, T **R, int cfirst, int clast,
                        int efirst, int elast,
                        T *spline_coefs, int numCoefs, T rMax,
                        T *lattice, T *latticeInv, T *sum)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  if (tid == 0)
    myR = R[blockIdx.x];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  __shared__ T rc[BS][3], re[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[4][4];
  if (tid < 16)
    A[tid>>2][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int Nc = clast - cfirst + 1;
  int Ne = elast - efirst + 1;
  int NBc = Nc/BS + ((Nc % BS) ? 1 : 0);
  int NBe = Ne/BS + ((Ne % BS) ? 1 : 0);
  T mysum = (T)0.0;
  for (int bc=0; bc < NBc; bc++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*bc+i)*BS + tid < 3*Nc)
        rc[0][i*BS + tid] = C[3*cfirst + (3*bc+i)*BS + tid];
    __syncthreads();
    int ptcl1 = cfirst+bc*BS + tid;
    for (int be=0; be < NBe; be++)
    {
      // Load block of positions from global memory
      for (int i=0; i<3; i++)
        if ((3*be+i)*BS + tid < 3*Ne)
          re[0][i*BS + tid] = myR[3*efirst + (3*be+i)*BS + tid];
      __syncthreads();
      // Now, loop over particles
      int end = (be+1)*BS < Ne ? BS : Ne-be*BS;
      for (int j=0; j<end; j++)
      {
        int ptcl2 = efirst + be*BS+j;
        T dx, dy, dz;
        dx = re[j][0] - rc[tid][0];
        dy = re[j][1] - rc[tid][1];
        dz = re[j][2] - rc[tid][2];
        T dist = min_dist(dx, dy, dz, L, Linv);
        if ((ptcl1 < (Nc+cfirst) ) && (ptcl2 < (Ne+efirst)))
          mysum += eval_1d_spline (dist, rMax, drInv, A, coefs);
      }
    }
    __syncthreads();
  }
  __shared__ T shared_sum[BS];
  shared_sum[tid] = mysum;
  __syncthreads();
  for (int s=BS>>1; s>0; s >>=1)
  {
    if (tid < s)
      shared_sum[tid] += shared_sum[tid+s];
    __syncthreads();
  }
  if (tid==0)
    sum[blockIdx.x] += shared_sum[0];
}

void
one_body_sum_PBC (float C[], float *R[], int cfirst, int clast, int efirst, int elast,
                  float spline_coefs[], int numCoefs, float rMax,
                  float lattice[], float latticeInv[], float sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  one_body_sum_PBC_kernel<float,BS><<<dimGrid,dimBlock>>>
  (C, R, cfirst, clast, efirst, elast,
   spline_coefs, numCoefs, rMax, lattice, latticeInv, sum);
}


void
one_body_sum_PBC (double C[], double *R[], int cfirst, int clast, int efirst, int elast,
                  double spline_coefs[], int numCoefs, double rMax,
                  double lattice[], double latticeInv[], double sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  one_body_sum_PBC_kernel<double,BS><<<dimGrid,dimBlock>>>
  (C, R, cfirst, clast, efirst, elast,
   spline_coefs, numCoefs, rMax, lattice, latticeInv, sum);
}



template<typename T, int BS>
__global__ void
one_body_ratio_PBC_kernel(T *C, T **R, int cfirst, int clast,
                          T *Rnew, int inew,
                          T *spline_coefs, int numCoefs, T rMax,
                          T *lattice, T *latticeInv, T *sum)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  __shared__ T myRnew[3], myRold[3];
  if (tid == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
  {
    myRnew[tid] = Rnew[3*blockIdx.x+tid];
    myRold[tid] = myR[3*inew+tid];
  }
  __syncthreads();
  __shared__ T coefs[MAX_COEFS];
  __shared__ T c[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[4][4];
  if (tid < 16)
    A[(tid>>2)][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int Nc = clast - cfirst + 1;
  int NB = Nc/BS + ((Nc % BS) ? 1 : 0);
  __shared__ T shared_sum[BS];
  shared_sum[tid] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*Nc)
        c[0][n] = C[3*cfirst + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = cfirst+b*BS + tid;
    T dx, dy, dz;
    dx = myRnew[0] - c[tid][0];
    dy = myRnew[1] - c[tid][1];
    dz = myRnew[2] - c[tid][2];
    T dist = min_dist(dx, dy, dz, L, Linv);
    T delta = eval_1d_spline (dist, rMax, drInv, A, coefs);
    dx = myRold[0] - c[tid][0];
    dy = myRold[1] - c[tid][1];
    dz = myRold[2] - c[tid][2];
    dist = min_dist(dx, dy, dz, L, Linv);
    delta -= eval_1d_spline (dist, rMax, drInv, A, coefs);
    if (ptcl1 < (Nc+cfirst) )
      shared_sum[tid] += delta;
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
      shared_sum[tid] += shared_sum[tid+s];
    __syncthreads();
  }
  if (tid==0)
    sum[blockIdx.x] += shared_sum[0];
}




void
one_body_ratio_PBC (float C[], float *R[], int first, int last,
                    float Rnew[], int inew,
                    float spline_coefs[], int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  one_body_ratio_PBC_kernel<float,BS><<<dimGrid,dimBlock>>>
  (C, R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, sum);
}



void
one_body_ratio_PBC (double C[], double *R[], int first, int last,
                    double Rnew[], int inew,
                    double spline_coefs[], int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sum[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  dim3 dimBlock(128);
  dim3 dimGrid(numWalkers);
  one_body_ratio_PBC_kernel<double,128><<<dimGrid,dimBlock>>>
  (C, R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, sum);
}


template<typename T, int BS>
__global__ void
one_body_ratio_grad_PBC_kernel(T const * __restrict__ C,
                               T const * const * __restrict__ R,
                               int cfirst, int clast,
                               T const * __restrict__ Rnew, int inew,
                               T const * __restrict__ spline_coefs,
                               int numCoefs, T rMax,
                               T const * __restrict__ lattice,
                               T const * __restrict__ latticeInv,
                               bool zero, T * __restrict__ ratio_grad)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  T const * __restrict__ myR = R[blockIdx.x];
  T myRnew[3], myRold[3];
  myRnew[0] = Rnew[3*blockIdx.x  ];
  myRnew[1] = Rnew[3*blockIdx.x+1];
  myRnew[2] = Rnew[3*blockIdx.x+2];
  myRold[0] = myR[3*inew  ];
  myRold[1] = myR[3*inew+1];
  myRold[2] = myR[3*inew+2];
  __shared__ T coefs[MAX_COEFS];
  /*
  __shared__ T L[3][3], Linv[3][3];
  */
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  /*
  if (tid < 9) {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }

  int index=0;
  __shared__ T images[27][3];
  if (tid < 3)
    for (T i=-1.0; i<=1.001; i+=1.0)
      for (T j=-1.0; j<=1.001; j+=1.0)
  for (T k=-1.0; k<=1.001; k+=1.0) {
    images[index][tid] =
      i*L[0][tid] + j*L[1][tid] + k*L[2][tid];
      index++;
  }

  __shared__ T A[12][4];
  if (tid < 16) {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  */
  __syncthreads();
  int Nc = clast - cfirst + 1;
  int NB = Nc/BS + ((Nc % BS) ? 1 : 0);
  __shared__ T shared_grad[BS][4];
  shared_grad[tid][0] = shared_grad[tid][1]
    = shared_grad[tid][2] = shared_grad[tid][3] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    T c[3];
    // Load block of positions from global memory
    int n = b*BS + tid;
    int ptcl1 = cfirst+b*BS + tid;
    if (n<Nc)
    {
        c[0] = C[3*ptcl1  ];
        c[1] = C[3*ptcl1+1];
        c[2] = C[3*ptcl1+2];
    }
    T dx, dy, dz, dist, delta, u, du, d2u;
    dx = myRold[0] - c[0];
    dy = myRold[1] - c[1];
    dz = myRold[2] - c[2];
    dist = CMC_min_dist(dx, dy, dz/*, L, Linv, images*/);
    delta =- CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs);
    dx = myRnew[0] - c[0];
    dy = myRnew[1] - c[1];
    dz = myRnew[2] - c[2];
    dist = CMC_min_dist(dx, dy, dz/*, L, Linv, images*/);
    CMC_eval_1d_spline_vgl (dist, rMax, drInv/*, A*/, coefs, u, du, d2u);
    delta += u;
    if (n < Nc)
    {
      du /= dist;
      shared_grad[tid][0] += delta;
      shared_grad[tid][1] += du * dx;
      shared_grad[tid][2] += du * dy;
      shared_grad[tid][3] += du * dz;
    }
  }
  __syncthreads();
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
    {
      shared_grad[tid][0] += shared_grad[tid+s][0];
      shared_grad[tid][1] += shared_grad[tid+s][1];
      shared_grad[tid][2] += shared_grad[tid+s][2];
      shared_grad[tid][3] += shared_grad[tid+s][3];
    }
    __syncthreads();
  }
  if (zero)
  {
    if ( tid < 4 ) ratio_grad[4*blockIdx.x+tid] = shared_grad[0][tid];
  }
  else
  {
    if ( tid < 4 ) ratio_grad[4*blockIdx.x+tid] += shared_grad[0][tid];
  }
}




template<typename T, int BS>
__global__ void
one_body_ratio_grad_PBC_kernel_fast(T *C, T **R, int cfirst, int clast,
                                    T *Rnew, int inew,
                                    T *spline_coefs, int numCoefs, T rMax,
                                    T *lattice, T *latticeInv, bool zero,
                                    T *ratio_grad)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  __shared__ T myRnew[3], myRold[3];
  if (tid == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
  {
    myRnew[tid] = Rnew[3*blockIdx.x+tid];
    myRold[tid] = myR[3*inew+tid];
  }
  __syncthreads();
  __shared__ T coefs[MAX_COEFS];
  __shared__ T c[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  int Nc = clast - cfirst + 1;
  int NB = Nc/BS + ((Nc % BS) ? 1 : 0);
  __shared__ T shared_sum[BS];
  __shared__ T shared_grad[BS][3];
  shared_sum[tid] = (T)0.0;
  shared_grad[tid][0] = shared_grad[tid][1] = shared_grad[tid][2] = 0.0f;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*Nc)
        c[0][n] = C[3*cfirst + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = cfirst+b*BS + tid;
    T dx, dy, dz, dist, delta, u, du, d2u;
    dx = myRold[0] - c[tid][0];
    dy = myRold[1] - c[tid][1];
    dz = myRold[2] - c[tid][2];
    dist = min_dist_fast(dx, dy, dz, L, Linv);
    delta =- eval_1d_spline (dist, rMax, drInv, A, coefs);
    dx = myRnew[0] - c[tid][0];
    dy = myRnew[1] - c[tid][1];
    dz = myRnew[2] - c[tid][2];
    dist = min_dist_fast(dx, dy, dz, L, Linv);
    eval_1d_spline_vgl (dist, rMax, drInv, A, coefs, u, du, d2u);
    delta += u;
    if (ptcl1 < (Nc+cfirst) )
    {
      du /= dist;
      shared_sum[tid] += delta;
      shared_grad[tid][0] += du * dx;
      shared_grad[tid][1] += du * dy;
      shared_grad[tid][2] += du * dz;
    }
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
    {
      shared_sum[tid] += shared_sum[tid+s];
      shared_grad[tid][0] += shared_grad[tid+s][0];
      shared_grad[tid][1] += shared_grad[tid+s][1];
      shared_grad[tid][2] += shared_grad[tid+s][2];
    }
    __syncthreads();
  }
  if (tid==0)
  {
    if (zero)
    {
      ratio_grad[4*blockIdx.x+0] = shared_sum[0];
      ratio_grad[4*blockIdx.x+1] = shared_grad[0][0];
      ratio_grad[4*blockIdx.x+2] = shared_grad[0][1];
      ratio_grad[4*blockIdx.x+3] = shared_grad[0][2];
    }
    else
    {
      ratio_grad[4*blockIdx.x+0] += shared_sum[0];
      ratio_grad[4*blockIdx.x+1] += shared_grad[0][0];
      ratio_grad[4*blockIdx.x+2] += shared_grad[0][1];
      ratio_grad[4*blockIdx.x+3] += shared_grad[0][2];
    }
  }
}




void
one_body_ratio_grad_PBC (float C[], float *R[], int first, int last,
                         float Rnew[], int inew,
                         float spline_coefs[], int numCoefs, float rMax,
                         float lattice[], float latticeInv[], bool zero,
                         float ratio_grad[], int numWalkers,
                         bool use_fast_image)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  CMC_PROFILING_BEGIN();
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  // if (use_fast_image)
  //   one_body_ratio_grad_kernel_fast<float,BS><<<dimGrid,dimBlock>>>
  //     (C, R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
  //      lattice, latticeInv, zero, ratio_grad);
  // else
  one_body_ratio_grad_PBC_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (C, R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, zero, ratio_grad);
  CMC_PROFILING_END(__LINE__);
}


void
one_body_ratio_grad_PBC (double C[], double *R[], int first, int last,
                         double Rnew[], int inew,
                         double spline_coefs[], int numCoefs, double rMax,
                         double lattice[], double latticeInv[], bool zero,
                         double ratio_grad[], int numWalkers, bool use_fast_image)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS = 128;
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  // if (use_fast_image)
  //   one_body_ratio_grad_kernel_fast<double,BS><<<dimGrid,dimBlock>>>
  //     (C, R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
  //      lattice, latticeInv, zero, ratio_grad);
  // else
  one_body_ratio_grad_PBC_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (C, R, first, last, Rnew, inew, spline_coefs, numCoefs, rMax,
   lattice, latticeInv, zero, ratio_grad);
  //CMC_PROFILING_END(__LINE__);
}




template<typename T>
__global__ void
one_body_update_kernel (T **R, int N, int iat)
{
  __shared__ T* myR;
  if (threadIdx.x == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (threadIdx.x < 3)
    myR[3*iat + threadIdx.x] = myR[3*N + threadIdx.x];
}

void
one_body_update(float *R[], int N, int iat, int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  dim3 dimBlock(32);
  dim3 dimGrid(numWalkers);
  one_body_update_kernel<float><<<dimGrid, dimBlock>>> (R, N, iat);
}

void
one_body_update(double *R[], int N, int iat, int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  dim3 dimBlock(3);
  dim3 dimGrid(numWalkers);
  one_body_update_kernel<double><<<dimGrid, dimBlock>>> (R, N, iat);
}





template<typename T, int BS>
__global__ void
one_body_grad_lapl_PBC_kernel(T *C, T **R, int cfirst, int clast,
                              int efirst, int elast,
                              T *spline_coefs, int numCoefs, T rMax,
                              T *lattice, T* latticeInv,
                              T *gradLapl, int row_stride)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR;
  if (tid == 0)
    myR = R[blockIdx.x];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  __shared__ T r[BS][3], c[BS][3];
  /*
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9) {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }

  __syncthreads();
  // if (tid == 31)
  //   printf ("1) coefs[] = %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f\n",
  // 	    coefs[0], coefs[1], coefs[2], coefs[3],
  // 	    coefs[4], coefs[5], coefs[6], coefs[7]);

  int index=0;
  __shared__ T images[27][3];
  if (tid < 3)
    for (T i=-1.0; i<=1.001; i+=1.0)
      for (T j=-1.0; j<=1.001; j+=1.0)
  for (T k=-1.0; k<=1.001; k+=1.0) {
    images[index][tid] =
      i*L[0][tid] + j*L[1][tid] + k*L[2][tid];
      index++;
  }
  __syncthreads();


  __shared__ T A[12][4];
  if (tid < 16) {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  */
  __syncthreads();
  int Nc = clast - cfirst + 1;
  int Ne = elast - efirst + 1;
  int NBc = (Nc+BS-1)/BS;
  int NBe = (Ne+BS-1)/BS;
  __shared__ T sGradLapl[BS][4];
  for (int be=0; be < NBe; be++)
  {
    // if (tid == 31)
    //   printf ("2) coefs[] = %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f\n",
    // 	    coefs[0], coefs[1], coefs[2], coefs[3],
    // 	    coefs[4], coefs[5], coefs[6], coefs[7]);
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*be+i)*BS + tid < 3*Ne)
        r[0][i*BS + tid] = myR[3*efirst + (3*be+i)*BS + tid];
    __syncthreads();
    int eptcl = efirst+be*BS + tid;
    int offset = blockIdx.x * row_stride + 4*be*BS + 4*efirst;
    sGradLapl[tid][0] = sGradLapl[tid][1] =
                          sGradLapl[tid][2] = sGradLapl[tid][3] = (T)0.0;
    for (int bc=0; bc < NBc; bc++)
    {
      // Load block of positions from global memory
      for (int i=0; i<3; i++)
        if ((3*bc+i)*BS + tid < 3*Nc)
          c[0][i*BS + tid] = C[3*cfirst + (3*bc+i)*BS + tid];
      __syncthreads();
      // Now, loop over particles
      int end = ((bc+1)*BS < Nc) ? BS : Nc-bc*BS;
      for (int j=0; j<end; j++)
      {
        int cptcl = cfirst + bc*BS+j;
        T dx, dy, dz, u, du, d2u;
        dx = r[tid][0] - c[j][0];
        dy = r[tid][1] - c[j][1];
        dz = r[tid][2] - c[j][2];
        T dist = CMC_min_dist(dx, dy, dz/*, L, Linv, images*/);
        // if (isinf(coefs[0]))
        //   printf ("3) c0=%1.5f coefs[] = %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f tid=%d\n", c0,
        // 	  coefs[0], spline_coefs[1], coefs[2], coefs[3],
        // 	  coefs[4], coefs[5], coefs[6], coefs[7], tid);
        CMC_eval_1d_spline_vgl (dist, rMax, drInv/*, A*/, coefs, u, du, d2u);
        // if (isinf(coefs[0]))
        //   printf ("4) coefs[] = %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f tid=%d\n",
        // 	  coefs[0], coefs[1], coefs[2], coefs[3],
        // 	  coefs[4], coefs[5], coefs[6], coefs[7], tid);
        // printf("drInv=%1.5f  dist=%1.5f coefs[1]=%1.5f A[0]=%1.5f\n",
        // 	 drInv, dist, coefs[1], A[0]);
        if (cptcl < (Nc+cfirst)  && (eptcl < (Ne+efirst)))
        {
          du /= dist;
          sGradLapl[tid][0] -= du * dx;
          sGradLapl[tid][1] -= du * dy;
          sGradLapl[tid][2] -= du * dz;
          sGradLapl[tid][3] -= d2u + 2.0*du;
        }
      }
      __syncthreads();
    }
    __syncthreads();
    for (int i=0; i<4; i++)
      if ((4*be+i)*BS + tid < 4*Ne)
        gradLapl[offset + i*BS +tid] += sGradLapl[0][i*BS+tid];
    __syncthreads();
  }
}


void
one_body_grad_lapl_PBC(float C[], float *R[], int e1_first, int e1_last,
                       int e2_first, int e2_last,
                       float spline_coefs[], int numCoefs, float rMax,
                       float lattice[], float latticeInv[],
                       float gradLapl[], int row_stride, int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  one_body_grad_lapl_PBC_kernel<float,BS><<<dimGrid,dimBlock>>>
  (C, R, e1_first, e1_last, e2_first, e2_last, spline_coefs, numCoefs,
   rMax, lattice, latticeInv,  gradLapl, row_stride);
}


void
one_body_grad_lapl_PBC(double C[], double *R[], int e1_first, int e1_last,
                       int e2_first, int e2_last,
                       double spline_coefs[], int numCoefs, double rMax,
                       double lattice[], double latticeInv[],
                       double gradLapl[], int row_stride, int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  one_body_grad_lapl_PBC_kernel<double,BS><<<dimGrid,dimBlock>>>
  (C, R, e1_first, e1_last, e2_first, e2_last, spline_coefs, numCoefs,
   rMax, lattice, latticeInv,  gradLapl, row_stride);
}


template<typename T, int BS>
__global__ void
one_body_NLratio_PBC_kernel(NLjobGPU<T> *jobs, T *C, int first, int last,
                            T *spline_coefs, int numCoefs, T rMax,
                            T *lattice, T *latticeInv)
{
  const int MAX_RATIOS = 18;
  int tid = threadIdx.x;
  __shared__ NLjobGPU<T> myJob;
  __shared__ T myRnew[MAX_RATIOS][3], myRold[3];
  if (tid == 0)
    myJob = jobs[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
    myRold[tid] = myJob.R[3*myJob.Elec+tid];
  for (int i=0; i<3; i++)
    if (i*BS + tid < 3*myJob.NumQuadPoints)
      myRnew[0][i*BS+tid] = myJob.QuadPoints[i*BS+tid];
  __syncthreads();
  T dr = rMax/(T)(numCoefs-3);
  T drInv = (T)1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  __shared__ T coefs[MAX_COEFS];
  __shared__ T c[BS][3];
  /*
  __shared__ T L[3][3], Linv[3][3];
  */
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  /*
  if (tid < 9) {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }

  __syncthreads();

  int index=0;
  __shared__ T images[27][3];
  if (tid < 3)
    for (T i=-1.0; i<=1.001; i+=1.0)
      for (T j=-1.0; j<=1.001; j+=1.0)
  for (T k=-1.0; k<=1.001; k+=1.0) {
    images[index][tid] =
      i*L[0][tid] + j*L[1][tid] + k*L[2][tid];
      index++;
  }
  __syncthreads();


  __shared__ T A[4][4];
  if (tid < 16)
    A[(tid>>2)][tid&3] = AcudaSpline[tid];
  */
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ T shared_sum[MAX_RATIOS][BS+1];
  for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    shared_sum[iq][tid] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*N)
        c[0][n] = C[3*first + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz;
    dx = myRold[0] - c[tid][0];
    dy = myRold[1] - c[tid][1];
    dz = myRold[2] - c[tid][2];
    T dist = CMC_min_dist_only(dx, dy, dz/*, L, Linv, images*/);
    T uOld = CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs);
    for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    {
      dx = myRnew[iq][0] - c[tid][0];
      dy = myRnew[iq][1] - c[tid][1];
      dz = myRnew[iq][2] - c[tid][2];
      dist = CMC_min_dist_only(dx, dy, dz/*, L, Linv, images*/);
      if (ptcl1 < (N+first))
        shared_sum[iq][tid] += CMC_eval_1d_spline (dist, rMax, drInv/*, A*/, coefs) - uOld;
    }
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
      for (int iq=0; iq < myJob.NumQuadPoints; iq++)
        shared_sum[iq][tid] += shared_sum[iq][tid+s];
    __syncthreads();
  }
  if (tid < myJob.NumQuadPoints)
    myJob.Ratios[tid] *= exp(-shared_sum[tid][0]);
}



template<typename T, int BS>
__global__ void
one_body_NLratio_PBC_kernel_fast(NLjobGPU<T> *jobs, T *C, int first, int last,
                                 T *spline_coefs, int numCoefs, T rMax,
                                 T *lattice, T *latticeInv)
{
  const int MAX_RATIOS = 18;
  int tid = threadIdx.x;
  __shared__ NLjobGPU<T> myJob;
  __shared__ T myRnew[MAX_RATIOS][3], myRold[3];
  if (tid == 0)
    myJob = jobs[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
    myRold[tid] = myJob.R[3*myJob.Elec+tid];
  for (int i=0; i<3; i++)
    if (i*BS + tid < 3*myJob.NumQuadPoints)
      myRnew[0][i*BS+tid] = myJob.QuadPoints[i*BS+tid];
  __syncthreads();
  T dr = rMax/(T)(numCoefs-3);
  T drInv = (T)1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  __shared__ T coefs[MAX_COEFS];
  __shared__ T c[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[4][4];
  if (tid < 16)
    A[(tid>>2)][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ T shared_sum[MAX_RATIOS][BS+1];
  for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    shared_sum[iq][tid] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*N)
        c[0][n] = C[3*first + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz;
    dx = myRold[0] - c[tid][0];
    dy = myRold[1] - c[tid][1];
    dz = myRold[2] - c[tid][2];
    T dist = min_dist_fast(dx, dy, dz, L, Linv);
    T uOld = eval_1d_spline (dist, rMax, drInv, A, coefs);
    for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    {
      dx = myRnew[iq][0] - c[tid][0];
      dy = myRnew[iq][1] - c[tid][1];
      dz = myRnew[iq][2] - c[tid][2];
      dist = min_dist_fast(dx, dy, dz, L, Linv);
      if (ptcl1 < (N+first))
        shared_sum[iq][tid] += eval_1d_spline (dist, rMax, drInv, A, coefs) - uOld;
    }
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
      for (int iq=0; iq < myJob.NumQuadPoints; iq++)
        shared_sum[iq][tid] += shared_sum[iq][tid+s];
    __syncthreads();
  }
  if (tid < myJob.NumQuadPoints)
    myJob.Ratios[tid] *= exp(-shared_sum[tid][0]);
}

/* no longer needed, use the template version instead
template<int BS>
__global__ void
one_body_NLratio_PBC_kernel(NLjobGPU<double> *jobs, double *C, int first, int last,
                            double *spline_coefs, int numCoefs, double rMax,
                            double *lattice, double *latticeInv)
{
  const int MAX_RATIOS = 18;
  int tid = threadIdx.x;
  __shared__ NLjobGPU<double> myJob;
  __shared__ double myRnew[MAX_RATIOS][3], myRold[3];
  if (tid == 0)
    myJob = jobs[blockIdx.x];
  __syncthreads();
  if (tid < 3 )
    myRold[tid] = myJob.R[3*myJob.Elec+tid];
  for (int i=0; i<3; i++)
    if (i*BS + tid < 3*myJob.NumQuadPoints)
      myRnew[0][i*BS+tid] = myJob.QuadPoints[i*BS+tid];
  __syncthreads();
  double dr = rMax/(double)(numCoefs-3);
  double drInv = 1.0/dr;
  __shared__ double coefs[MAX_COEFS];
  __shared__ double c[BS][3];
  __shared__ double L[3][3], Linv[3][3];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ double images[27][3];
  int index=0;
  if (tid < 3)
    for (float i=-1.0; i<=1.001; i+=1.0)
      for (float j=-1.0; j<=1.001; j+=1.0)
        for (float k=-1.0; k<=1.001; k+=1.0)
        {
          images[index][tid] =
            i*L[0][tid] + j*L[1][tid] + k*L[2][tid];
          index++;
        }
  __syncthreads();
  __shared__ double A[4][4];
  if (tid < 16)
    A[(tid>>2)][tid&3] = AcudaSpline[tid];
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ double shared_sum[MAX_RATIOS][BS+1];
  for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    shared_sum[iq][tid] = (double)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
    {
      int n = i*BS + tid;
      if ((3*b+i)*BS + tid < 3*N)
        c[0][n] = C[3*first + (3*b+i)*BS + tid];
    }
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    double dx, dy, dz;
    dx = myRold[0] - c[tid][0];
    dy = myRold[1] - c[tid][1];
    dz = myRold[2] - c[tid][2];
    double dist = min_dist(dx, dy, dz, L, Linv, images);
    double uOld = eval_1d_spline (dist, rMax, drInv, A, coefs);
    for (int iq=0; iq<myJob.NumQuadPoints; iq++)
    {
      dx = myRnew[iq][0] - c[tid][0];
      dy = myRnew[iq][1] - c[tid][1];
      dz = myRnew[iq][2] - c[tid][2];
      dist = min_dist(dx, dy, dz, L, Linv, images);
      if (ptcl1 != myJob.Elec && (ptcl1 < (N+first)))
        shared_sum[iq][tid] += eval_1d_spline (dist, rMax, drInv, A, coefs) - uOld;
    }
    __syncthreads();
  }
  for (int s=(BS>>1); s>0; s>>=1)
  {
    if (tid < s)
      for (int iq=0; iq < myJob.NumQuadPoints; iq++)
        shared_sum[iq][tid] += shared_sum[iq][tid+s];
    __syncthreads();
  }
  if (tid < myJob.NumQuadPoints)
    myJob.Ratios[tid] *= exp(-shared_sum[tid][0]);
}

*/



void
one_body_NLratios_PBC(NLjobGPU<float> jobs[], float C[], int first, int last,
                      float spline_coefs[], int numCoefs, float rMax,
                      float lattice[], float latticeInv[], float sim_cell_radius,
                      int numjobs)
{
  if (numjobs==0) return;
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=32;
  dim3 dimBlock(BS);
  CMC_PROFILING_BEGIN();
  while (numjobs > 65535)
  {
    dim3 dimGrid(65535);
    if (rMax <= sim_cell_radius)
    {
      // fprintf (stderr, "Using fast J1 NL kernel.\n");
      one_body_NLratio_PBC_kernel_fast<float, BS><<<dimGrid,dimBlock>>>
      (jobs, C, first, last, spline_coefs, numCoefs, rMax,
       lattice, latticeInv);
    }
    else
    {
      // fprintf (stderr, "Using slow J1 NL kernel.\n");
      one_body_NLratio_PBC_kernel<float, BS><<<dimGrid,dimBlock>>>
      (jobs, C, first, last, spline_coefs, numCoefs, rMax,
       lattice, latticeInv);
    }
    numjobs -= 65535;
    jobs += 65535;
  }
  dim3 dimGrid(numjobs);
  if (rMax <= sim_cell_radius)
  {
    // fprintf (stderr, "Using fast J1 NL kernel.\n");
    one_body_NLratio_PBC_kernel_fast<float, BS><<<dimGrid,dimBlock>>>
    (jobs, C, first, last, spline_coefs, numCoefs, rMax,
     lattice, latticeInv);
  }
  else
  {
    // fprintf (stderr, "Using slow J1 NL kernel.\n");
    one_body_NLratio_PBC_kernel<float, BS><<<dimGrid,dimBlock>>>
    (jobs, C, first, last, spline_coefs, numCoefs, rMax,
     lattice, latticeInv);
  }
  CMC_PROFILING_END(__LINE__);
}


void
one_body_NLratios_PBC(NLjobGPU<double> jobs[], double C[], int first, int last,
                      double spline_coefs[], int numCoefs, double rMax,
                      double lattice[], double latticeInv[], double sim_cell_radius,
                      int numjobs)
{
  if (numjobs==0) return;
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  const int BS=32;
  dim3 dimBlock(BS);
  /*
  int blockx = numjobs % 65535;
  int blocky = numjobs / 65535 + 1;
  dim3 dimGrid(blockx, blocky);
  one_body_NLratio_PBC_kernel<double, BS><<<dimGrid,dimBlock>>>
  (jobs, C, first, last, spline_coefs, numCoefs, rMax,
   lattice, latticeInv);
  */
  while (numjobs > 65535)
  {
    dim3 dimGrid(65535);
    if (rMax <= sim_cell_radius)
    {
      // fprintf (stderr, "Using fast J1 NL kernel.\n");
      one_body_NLratio_PBC_kernel_fast<double, BS><<<dimGrid,dimBlock>>>
      (jobs, C, first, last, spline_coefs, numCoefs, rMax,
       lattice, latticeInv);
    }
    else
    {
      // fprintf (stderr, "Using slow J1 NL kernel.\n");
      one_body_NLratio_PBC_kernel<double, BS><<<dimGrid,dimBlock>>>
      (jobs, C, first, last, spline_coefs, numCoefs, rMax,
       lattice, latticeInv);
    }
    numjobs -= 65535;
    jobs += 65535;
  }
  dim3 dimGrid(numjobs);
  if (rMax <= sim_cell_radius)
  {
    // fprintf (stderr, "Using fast J1 NL kernel.\n");
    one_body_NLratio_PBC_kernel_fast<double, BS><<<dimGrid,dimBlock>>>
    (jobs, C, first, last, spline_coefs, numCoefs, rMax,
     lattice, latticeInv);
  }
  else
  {
    // fprintf (stderr, "Using slow J1 NL kernel.\n");
    one_body_NLratio_PBC_kernel<double, BS><<<dimGrid,dimBlock>>>
    (jobs, C, first, last, spline_coefs, numCoefs, rMax,
     lattice, latticeInv);
  }
  //CMC_PROFILING_END(__LINE__);
}



template<typename T, int BS>
__global__ void
one_body_grad_PBC_kernel(T const * const * __restrict__ R, int iat,
                         T const * __restrict__ C, int first, int last,
                         T const * __restrict__ spline_coefs, int numCoefs, T rMax,
                         T const * __restrict__ lattice,
                         T const * __restrict__ latticeInv,
                         bool zeroOut, T * __restrict__ grad)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  T const * __restrict__ myR = R[blockIdx.x];
  T r[3];
  r[0] = myR[3*iat  ];
  r[1] = myR[3*iat+1];
  r[2] = myR[3*iat+2];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  /*
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9) {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }

  __shared__ T A[12][4];
  if (tid < 16) {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();

  int index=0;
  __shared__ T images[27][3];
  if (tid < 3)
    for (T i=-1.0; i<=1.001; i+=1.0)
      for (T j=-1.0; j<=1.001; j+=1.0)
  for (T k=-1.0; k<=1.001; k+=1.0) {
    images[index][tid] =
      i*L[0][tid] + j*L[1][tid] + k*L[2][tid];
      index++;
  }

  */
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ T sGrad[BS][3];
  sGrad[tid][0]   = sGrad[tid][1] = sGrad[tid][2] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    T c[3];
    int n = b*BS + tid;
    int ptcl1 = first + n;
    // Load block of positions from global memory
    if (n < N)
    {
      c[0] = C[3*ptcl1  ];
      c[1] = C[3*ptcl1+1];
      c[2] = C[3*ptcl1+2];
    }
    T dx, dy, dz, u, du, d2u;
    dx = r[0] - c[0];
    dy = r[1] - c[1];
    dz = r[2] - c[2];
    T dist = CMC_min_dist(dx, dy, dz/*, L, Linv, images*/);
    CMC_eval_1d_spline_vgl (dist, rMax, drInv/*, A*/, coefs, u, du, d2u);
    if (n < N)
    {
      du /= dist;
      sGrad[tid][0] += du * dx;
      sGrad[tid][1] += du * dy;
      sGrad[tid][2] += du * dz;
    }
  }
  __syncthreads();
  // Do reduction across threads in block
  for (int s=BS>>1; s>0; s>>=1)
  {
    if (tid < s)
    {
      sGrad[tid][0] += sGrad[tid+s][0];
      sGrad[tid][1] += sGrad[tid+s][1];
      sGrad[tid][2] += sGrad[tid+s][2];
    }
    __syncthreads();
  }
  if (tid < 3)
  {
    if (zeroOut)
      grad[3*blockIdx.x + tid] = sGrad[0][tid];
    else
      grad[3*blockIdx.x + tid] += sGrad[0][tid];
  }
}


template<typename T, int BS>
__global__ void
one_body_grad_PBC_kernel_fast(T **R, int iat, T *C, int first, int last,
                              T *spline_coefs, int numCoefs, T rMax,
                              T *lattice, T *latticeInv, bool zeroOut, T *grad)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR, r[3];
  if (tid == 0)
    myR = R[blockIdx.x];
  __syncthreads();
  if (tid < 3)
    r[tid] = myR[3*iat+tid];
  __shared__ T coefs[MAX_COEFS];
  if (tid < numCoefs)
    coefs[tid] = spline_coefs[tid];
  __shared__ T c[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  int N = last - first + 1;
  int NB = N/BS + ((N % BS) ? 1 : 0);
  __shared__ T sGrad[BS][3];
  sGrad[tid][0]   = sGrad[tid][1] = sGrad[tid][2] = (T)0.0;
  for (int b=0; b < NB; b++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*b+i)*BS + tid < 3*N)
        c[0][i*BS + tid] = C[3*first + (3*b+i)*BS + tid];
    __syncthreads();
    int ptcl1 = first+b*BS + tid;
    T dx, dy, dz, u, du, d2u;
    dx = r[0] - c[tid][0];
    dy = r[1] - c[tid][1];
    dz = r[2] - c[tid][2];
    T dist = min_dist_fast(dx, dy, dz, L, Linv);
    eval_1d_spline_vgl (dist, rMax, drInv, A, coefs, u, du, d2u);
    if (ptcl1 < (N+first))
    {
      du /= dist;
      sGrad[tid][0] += du * dx;
      sGrad[tid][1] += du * dy;
      sGrad[tid][2] += du * dz;
    }
    __syncthreads();
  }
  // Do reduction across threads in block
  for (int s=BS>>1; s>0; s>>=1)
  {
    if (tid < s)
    {
      sGrad[tid][0] += sGrad[tid+s][0];
      sGrad[tid][1] += sGrad[tid+s][1];
      sGrad[tid][2] += sGrad[tid+s][2];
    }
    __syncthreads();
  }
  if (tid < 3)
  {
    if (zeroOut)
      grad[3*blockIdx.x + tid] = sGrad[0][tid];
    else
      grad[3*blockIdx.x + tid] += sGrad[0][tid];
  }
}




void
one_body_gradient_PBC (float *Rlist[], int iat, float C[], int first, int last,
                       float spline_coefs[], int num_coefs, float rMax,
                       float lattice[], float latticeInv[], bool zeroSum,
                       float grad[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  CMC_PROFILING_BEGIN();
  // if (sim_cell_radius >= rMax)
  //   one_body_grad_kernel_fast<float,BS><<<dimGrid,dimBlock>>>
  //     (Rlist, iat, C, first, last, spline_coefs, num_coefs, rMax,
  //      L, Linv, zeroSum, grad);
  // else
  one_body_grad_PBC_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (Rlist, iat, C, first, last, spline_coefs, num_coefs, rMax,
   lattice, latticeInv, zeroSum, grad);
  CMC_PROFILING_END(__LINE__);
}


void
one_body_gradient_PBC (double *Rlist[], int iat, double C[], int first, int last,
                       double spline_coefs[], int num_coefs, double rMax,
                       double lattice[], double latticeInv[], bool zeroSum,
                       double grad[], int numWalkers)
{
  if (!AisInitializedPBC)
    cuda_spline_init_PBC();
  const int BS=128;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  //CMC_PROFILING_BEGIN();
  COPY_LATTICE_DP_TO_SP();
  one_body_grad_PBC_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (Rlist, iat, C, first, last, spline_coefs, num_coefs, rMax,
   lattice, latticeInv, zeroSum, grad);
  //CMC_PROFILING_END(__LINE__);
}



template<typename T, int BS, unsigned COMPLEX>
__global__ void
one_body_derivs_PBC_kernel(T* C, T **R, T **gradLogPsi,
                           int cfirst, int clast,
                           int efirst, int elast,
                           int numCoefs, T rMax,
                           T *lattice, T *latticeInv,
                           T **derivs)
{
  T dr = rMax/(T)(numCoefs-3);
  T drInv = 1.0/dr;
  __syncthreads();
  // Safety for rounding error
  rMax *= (T)0.999999;
  int tid = threadIdx.x;
  __shared__ T *myR, *myGrad, *myDerivs;
  if (tid == 0)
  {
    myR      =          R[blockIdx.x];
    myGrad   = gradLogPsi[blockIdx.x];
    myDerivs =     derivs[blockIdx.x];
  }
  __shared__ T sderivs[MAX_COEFS][2];
  __shared__ T r[BS][3], c[BS][3];
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9)
  {
    L[0][tid] = lattice[tid];
    Linv[0][tid] = latticeInv[tid];
  }
  __shared__ T A[12][4];
  if (tid < 16)
  {
    A[0+(tid>>2)][tid&3] = AcudaSpline[tid+0];
    A[4+(tid>>2)][tid&3] = AcudaSpline[tid+16];
    A[8+(tid>>2)][tid&3] = AcudaSpline[tid+32];
  }
  __syncthreads();
  sderivs[tid][0] = T();
  sderivs[tid][1] = T();
  int Nc = clast - cfirst + 1;
  int Ne = elast - efirst + 1;
  int NBc = (Nc+BS-1)/BS;
  int NBe = (Ne+BS-1)/BS;
  __shared__ T sGrad[BS][3];
  for (int be=0; be < NBe; be++)
  {
    // Load block of positions from global memory
    for (int i=0; i<3; i++)
      if ((3*be+i)*BS + tid < 3*Ne)
      {
        int outoff = i*BS+tid;
        int inoff  = outoff + 3*efirst + 3*be*BS;
        r[0][outoff]    =     myR[inoff];
        sGrad[0][outoff] = myGrad[inoff*COMPLEX];
      }
    __syncthreads();
    int eptcl = efirst+be*BS + tid;
    for (int bc=0; bc < NBc; bc++)
    {
      // Load block of positions from global memory
      for (int i=0; i<3; i++)
        if ((3*bc+i)*BS + tid < 3*Nc)
          c[0][i*BS + tid] = C[3*cfirst + (3*bc+i)*BS + tid];
      __syncthreads();
      // Now, loop over particles
      int end = min(BS, Nc-bc*BS);
      for (int j=0; j<end; j++)
      {
        T dx, dy, dz;
        dx = c[j][0] - r[tid][0];
        dy = c[j][1] - r[tid][1];
        dz = c[j][2] - r[tid][2];
        T dist = min_dist(dx, dy, dz, L, Linv);
        T distInv = 1.0f/dist;
        T s = dist * drInv;
        T sf = floorf (s);
        int index = (int)sf;
        T t = s - sf;
        T t2 = t*t;
        T t3 = t*t2;
        T v0 = (A[0][0]*t3 + A[0][1]*t2 + A[0][2]*t + A[0][3]);
        T v1 = (A[1][0]*t3 + A[1][1]*t2 + A[1][2]*t + A[1][3]);
        T v2 = (A[2][0]*t3 + A[2][1]*t2 + A[2][2]*t + A[2][3]);
        T v3 = (A[3][0]*t3 + A[3][1]*t2 + A[3][2]*t + A[3][3]);
        for (int id=0; id<BS; id++)
          if (tid == id && eptcl <= elast && (dist < rMax))
          {
            sderivs[index+0][0] += v0;
            sderivs[index+1][0] += v1;
            sderivs[index+2][0] += v2;
            sderivs[index+3][0] += v3;
          }
        T prefact = (dx*sGrad[tid][0] + dy*sGrad[tid][1] + dz*sGrad[tid][2])*distInv;
        T du0 = drInv * (A[4][0]*t3 + A[4][1]*t2 + A[4][2]*t + A[4][3]);
        T du1 = drInv * (A[5][0]*t3 + A[5][1]*t2 + A[5][2]*t + A[5][3]);
        T du2 = drInv * (A[6][0]*t3 + A[6][1]*t2 + A[6][2]*t + A[6][3]);
        T du3 = drInv * (A[7][0]*t3 + A[7][1]*t2 + A[7][2]*t + A[7][3]);
        // This is the dot (gradu, grad_log_psi) term.
        v0 = 2.0f* prefact * du0;
        v1 = 2.0f* prefact * du1;
        v2 = 2.0f* prefact * du2;
        v3 = 2.0f* prefact * du3;
        // This is the lapl u term
        v0 -= drInv*drInv*(A[ 8][0]*t3 + A[ 8][1]*t2 + A[ 8][2]*t + A[ 8][3]) + 2.0f*du0*distInv;
        v1 -= drInv*drInv*(A[ 9][0]*t3 + A[ 9][1]*t2 + A[ 9][2]*t + A[ 9][3]) + 2.0f*du1*distInv;
        v2 -= drInv*drInv*(A[10][0]*t3 + A[10][1]*t2 + A[10][2]*t + A[10][3]) + 2.0f*du2*distInv;
        v3 -= drInv*drInv*(A[11][0]*t3 + A[11][1]*t2 + A[11][2]*t + A[11][3]) + 2.0f*du3*distInv;
        for (int id=0; id<BS; id++)
          if (tid == id && eptcl <= elast && (dist < rMax))
          {
            sderivs[index+0][1] += v0;
            sderivs[index+1][1] += v1;
            sderivs[index+2][1] += v2;
            sderivs[index+3][1] += v3;
          }
      }
      __syncthreads();
    }
  }
  sderivs[tid][1] *= 0.5f;
  if (tid < 2*numCoefs)
    myDerivs[tid] = -sderivs[0][tid];
  if (tid+BS < 2*numCoefs)
    myDerivs[tid+BS] = -sderivs[0][tid+BS];
}


void
one_body_derivs_PBC(float C[], float *R[], float *gradLogPsi[],
                    int cfirst, int clast,
                    int efirst, int elast,
                    int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sim_cell_radius,
                    float *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  if (sim_cell_radius >= rMax)
    one_body_derivs_PBC_kernel<float,BS,1><<<dimGrid,dimBlock>>>
    (C, R, gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    one_body_derivs_PBC_kernel<float,BS,1><<<dimGrid,dimBlock>>>
    (C, R, gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
}



void
one_body_derivs_PBC(double C[], double *R[], double *gradLogPsi[],
                    int cfirst, int clast,
                    int efirst, int elast,
                    int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sim_cell_radius,
                    double *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  if (sim_cell_radius >= rMax)
    one_body_derivs_PBC_kernel<double,BS,1><<<dimGrid,dimBlock>>>
    (C, R, gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    one_body_derivs_PBC_kernel<double,BS,1><<<dimGrid,dimBlock>>>
    (C, R, gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
}


// Ye: use offset to recycle the old routines
// block size can be further optimized.
#ifdef QMC_COMPLEX
void
one_body_derivs_PBC(float C[], float *R[], std::complex<float> *gradLogPsi[],
                    int cfirst, int clast,
                    int efirst, int elast,
                    int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sim_cell_radius,
                    float *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);

  if (sim_cell_radius >= rMax)
    one_body_derivs_PBC_kernel<float,BS,2><<<dimGrid,dimBlock>>>
    (C, R, (float**)gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    one_body_derivs_PBC_kernel<float,BS,2><<<dimGrid,dimBlock>>>
    (C, R, (float**)gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
}



void
one_body_derivs_PBC(double C[], double *R[], std::complex<double> *gradLogPsi[],
                    int cfirst, int clast,
                    int efirst, int elast,
                    int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sim_cell_radius,
                    double *derivs[], int numWalkers)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);

  if (sim_cell_radius >= rMax)
    one_body_derivs_PBC_kernel<double,BS,2><<<dimGrid,dimBlock>>>
    (C, R, (double**)gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
  else
    one_body_derivs_PBC_kernel<double,BS,2><<<dimGrid,dimBlock>>>
    (C, R, (double**)gradLogPsi, cfirst, clast, efirst, elast, numCoefs,
     rMax, lattice, latticeInv, derivs);
}
#endif


void testPBC()
{
  dim3 dimBlock(32);
  dim3 dimGrid(1000);
  float *R[1000];
  float L[9], Linv[9];
  float spline_coefs[10];
  float dr = 0.1;
  float sum[1000];
  two_body_sum_PBC_kernel<float,32><<<dimGrid,dimBlock>>>(R, 0, 100, 0, 100, spline_coefs, 10, dr,
      L, Linv, sum);
}
