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
    
    


#include <cuda.h>
#include <stdio.h>

const int BS=32;

__global__ void kernel(float *a, size_t N)
{
  int tid = threadIdx.x;
  __shared__ float s[BS];
  int blocks = (N+BS-1)/BS;
  float sum = 0.0f;
  for (int ib=0; ib<blocks; ib++)
  {
    int off = ib*BS+tid;
    s[tid] = a[off];
    for (int skip=16; skip>0; skip>>=1)
      if (tid+skip < N && tid < skip)
        s[tid] += s[tid+skip];
    sum += s[0];
  }
  a[0] = sum;
}


main()
{
  CUdevice cuDevice[4];
  CUcontext context[4];
  float *buffer[4];
  int dev_ids[4];
  int j=0;
  for (int i=0; i<5; i++)
  {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    fprintf (stderr, "Device %d  Compute ability = %d.%d\n", i, prop.major, prop.minor);
    if (prop.minor > 1)
    {
      dev_ids[j] = i;
      j++;
    }
  }
  size_t N = 20000000;
  float *host = (float*)malloc(N*sizeof(float));
  float sum = 0.0;
  for (int i=0; i<N; i++)
  {
    host[i] = drand48();
    sum += host[i];
  }
  for (int i=0; i<4; i++)
  {
    CUcontext ctx;
    cuDeviceGet(&(cuDevice[i]), dev_ids[i]);
    cuCtxCreate (&(context[i]), CU_CTX_SCHED_SPIN, cuDevice[i]);
    // = cudaMalloc(&(buffer[i]), (size_t)1000000*i);
    cudaError_t result = cudaMalloc(&(buffer[i]), (size_t)N*sizeof(float));
    cudaMemcpyAsync(buffer[i], host, N*sizeof(float), cudaMemcpyHostToDevice);
    if (result != cudaSuccess)
    {
      fprintf (stderr, "Error allocating memory on device %d\n", i);
      abort();
    }
    else
      fprintf (stderr, "Device pointer = %lp\n", buffer[i]);
    cuCtxPopCurrent(&ctx);
  }
  #pragma omp parallel for
  for (int i=0; i<4; i++)
  {
    CUcontext ctx;
    fprintf (stderr, "Before kernel, i=%d\n", i);
    cuCtxPushCurrent (context[i]);
    dim3 dimBlock(BS);
    dim3 dimGrid(1);
    kernel<<<dimGrid,dimBlock>>>(buffer[i], N);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf (stderr, "CUDA error in kernel, i=%d\n", i);
    }
    float devsum;
    cudaMemcpy(&devsum, buffer[i], 4, cudaMemcpyDeviceToHost);
    fprintf (stderr, "Exact sum = %f  device sum = %f\n", sum, devsum);
    cuCtxPopCurrent(&ctx);
    fprintf (stderr, "After kernel, i=%d\n", i);
  }
}
