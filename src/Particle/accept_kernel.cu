//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



template<typename T, int BS>
__global__ void accept_kernel (T** Rlist, T* Rnew,
                               int* toAccept, int iat, int N)
{
  int tid = threadIdx.x;
  __shared__ T* myR[BS];
  __shared__ T Rnew_shared[BS*3];
  __shared__ int accept_shared[BS];
  int block = blockIdx.x;
  for (int i=0; i<3; i++)
    if ((3*block+i)*BS + tid < 3*N)
      Rnew_shared[i*BS + tid] = Rnew[(3*block+i)*BS + tid];
  __syncthreads();
  if (block*BS + tid < N)
  {
    myR[tid] = Rlist[block*BS+tid] + 3*iat;
    accept_shared[tid] = toAccept[block*BS+tid];
  }
  __syncthreads();
  // if (block*BS + tid < N && accept_shared[tid]) {
  //   myR[tid][0] = Rnew_shared[3*tid+0];
  //   myR[tid][1] = Rnew_shared[3*tid+1];
  //   myR[tid][2] = Rnew_shared[3*tid+2];
  // }
  // return;
  for (int i=0; i<3; i++)
  {
    int index = i*BS + tid;
    int iw = index / 3;
    int dim = index % 3;
    if (iw+block*BS < N && accept_shared[iw])
      myR[iw%BS][dim] = Rnew_shared[index];
  }
}

#include <cstdio>
#include "../CUDA/gpu_misc.h"

void
accept_move_GPU_cuda (float** Rlist, float new_pos[],
                      int toAccept[], int iat, int N)
{
  const int BS=32;
  int NB = N / BS + ((N % BS) ? 1 : 0);
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  accept_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (Rlist, new_pos, toAccept, iat, N);
}


void
accept_move_GPU_cuda (double** Rlist, double* new_pos,
                      int* toAccept, int iat, int N)
{
  const int BS=32;
  int NB = N / BS + ((N % BS) ? 1 : 0);
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  accept_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (Rlist, new_pos, toAccept, iat, N);
}


template<typename T, int BS>
__global__ void NL_move_kernel (T** Rlist, T* Rnew, int N)
{
  __shared__ T Rnew_shared[BS][3];
  __shared__ T* Rlist_shared[BS];
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x+i)*BS+threadIdx.x;
    if (off < 3*N)
      Rnew_shared[0][i*BS+threadIdx.x] = Rnew[off];
  }
  int off = blockIdx.x * BS + threadIdx.x;
  if (off < N)
    Rlist_shared[threadIdx.x] = Rlist[off];
  __syncthreads();
  int end = min (BS, N-blockIdx.x*BS);
  for (int i=0; i<end; i++)
    if (threadIdx.x < 3)
      Rlist_shared[i][threadIdx.x] = Rnew_shared[i][threadIdx.x];
}

void
NL_move_cuda (float **Rlist, float* new_pos, int N)
{
  const int BS=32;
  int NB = (N+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  NL_move_kernel<float,BS><<<dimGrid,dimBlock>>>
  (Rlist, new_pos, N);
}

void
NL_move_cuda (double **Rlist, double new_pos[], int N)
{
  const int BS=32;
  int NB = (N+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  NL_move_kernel<double,BS><<<dimGrid,dimBlock>>>
  (Rlist, new_pos, N);
}
