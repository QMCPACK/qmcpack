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
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <assert.h>
#include <thrust/complex.h>

/*                   !!!WARNING!!!
   Kernels in this file strongly depend on warp-synchronous behavior
   of current generations of GPUs. Any change to that behavior
   as well as any change to the warp size will break the code!
   In such case extra synchronizations are necessary.
*/

template<typename T, int BS> __global__
void phase_factor_kernel (T *kPoints, int *makeTwoCopies,
                          T *pos, T **phi_in, T **phi_out,
                          int num_splines, int num_walkers)
{
  __shared__ T in_shared[2*BS+1], kPoints_s[BS][3],
             pos_s[BS][3];
  volatile __shared__ T out_shared[2*BS+1];
  __shared__ T *phi_in_ptr[BS], *phi_out_ptr[BS];
  int tid = threadIdx.x;
  assert(warpSize == 32);
#pragma unroll
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x+i)*BS + tid;
    if (off < 3*num_walkers)
      pos_s[0][i*BS + tid] =  pos[off];
  }
  if (blockIdx.x*BS+tid < num_walkers)
  {
    phi_in_ptr[tid]  = phi_in[blockIdx.x*BS+tid];
    phi_out_ptr[tid] = phi_out[blockIdx.x*BS+tid];
  }
  //__syncthreads();
  int nb = (num_splines + BS-1)/BS;
  int outIndex=0;
  int outBlock=0;
  int m2c;
  volatile __shared__ int m2c_ps[BS];
  int numWrite = min(BS, num_walkers-blockIdx.x*BS);
  for (int block=0; block<nb; block++)
  {
    // Load kpoints into shared memory
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*num_splines)
        kPoints_s[0][i*BS+tid] = kPoints[off];
    }
    // Load makeTwoCopies with coallesced reads
    if (block*BS+tid < num_splines)
    {
      if(makeTwoCopies[block*BS + tid])
        m2c = 1;
      else
        m2c = 0;
    }
    else
      m2c = 0;
    //prefix sum of m2c array
    m2c_ps[tid] = m2c+1;
    if(tid >= 1)
      m2c_ps[tid] += m2c_ps[tid-1];
    if(tid >= 2)
      m2c_ps[tid] += m2c_ps[tid-2];
    if(tid >= 4)
      m2c_ps[tid] += m2c_ps[tid-4];
    if(tid >= 8)
      m2c_ps[tid] += m2c_ps[tid-8];
    if(tid >= 16)
      m2c_ps[tid] += m2c_ps[tid-16];
    if(tid > 0)
      outIndex  = m2c_ps[tid-1];
    T s, c;
    int end = min (BS, num_splines-block*BS);
    if(tid < end)
      for (int i=0; i<numWrite; i++)
      {
        if ((2*block)*BS+tid < 2*num_splines)
          in_shared[tid   ] = phi_in_ptr[i][(2*block+0)*BS+tid];
        if ((2*block)*BS+tid + end < 2*num_splines)
          in_shared[tid+end] = phi_in_ptr[i][(2*block)*BS+tid + end];
        // Compute e^{-ikr}
        T phase = -(pos_s[i][0]*kPoints_s[tid][0] +
                    pos_s[i][1]*kPoints_s[tid][1] +
                    pos_s[i][2]*kPoints_s[tid][2]);
        sincos(phase, &s, &c);
        T phi_real = in_shared[2*tid]*c - in_shared[2*tid+1]*s;
        T phi_imag = in_shared[2*tid]*s + in_shared[2*tid+1]*c;
        out_shared[outIndex] = phi_real;
        if(m2c)
        {
          out_shared[outIndex + 1] = phi_imag;
        }
        phi_out_ptr[i][outBlock+tid]= out_shared[tid];
        if(tid + end < m2c_ps[end-1])
        {
          phi_out_ptr[i][outBlock + tid + end] = out_shared[tid+end];
        }
      }
    outBlock+= m2c_ps[end-1];
  }
}


template<typename T, int BS> __global__
void phase_factor_kernel_new (T *kPoints, int *makeTwoCopies,
                              T *pos, T **phi_in, T **phi_out,
                              int num_splines)
{
  __shared__ T in_shared[2*BS], out_shared[2*BS], kPoints_s[BS][3],
             pos_s[3];
  __shared__ int m2c[BS];
  __shared__ T *phi_in_ptr, *phi_out_ptr;
  int tid = threadIdx.x;
  if (tid < 3)
    pos_s[tid] = pos[3*blockIdx.x+tid];
  if (tid == 0)
  {
    phi_in_ptr = phi_in[blockIdx.x];
    phi_out_ptr = phi_out[blockIdx.x];
  }
  int NB = (num_splines+BS-1)/BS;
  int outIndex=0, outBlock=0;
  for (int ib=0; ib<NB; ib++)
  {
    for (int i=0; i<3; i++)
      kPoints_s[0][i*BS+tid] = kPoints[(3*ib+i)*BS+tid];
    T phase = -(kPoints_s[tid][0]*pos_s[0] +
                kPoints_s[tid][1]*pos_s[1] +
                kPoints_s[tid][2]*pos_s[2]);
    T s, c;
    sincosf (phase, &s, &c);
    int off = 2*ib*BS + tid;
    in_shared[tid]    = phi_in_ptr[off];
    in_shared[tid+BS] = phi_in_ptr[off+BS];
    T phi_real = in_shared[2*tid]*c - in_shared[2*tid+1]*s;
    T phi_imag = in_shared[2*tid]*s + in_shared[2*tid+1]*c;
    m2c[tid] = makeTwoCopies[ib*BS + tid];
    int iend = min (BS, num_splines - ib*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == i)
        out_shared[outIndex] = phi_real;
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        phi_out_ptr[outBlock*BS+tid] = out_shared[tid];
        outIndex = 0;
        outBlock++;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == i)
          out_shared[outIndex] = phi_imag;
        outIndex++;
      }
      __syncthreads();
      if (outIndex == BS)
      {
        phi_out_ptr[outBlock*BS+tid] = out_shared[tid];
        outIndex = 0;
        outBlock++;
      }
      __syncthreads();
    }
  }
  if (tid < outIndex)
    phi_out_ptr[outBlock*BS+tid] = out_shared[tid];
}


// Original implementation

template<typename T, int BS> __global__
void phase_factor_kernel (T *kPoints, int *makeTwoCopies,
                          T *pos, T **phi_in, T **phi_out,
                          T **grad_lapl_in, T **grad_lapl_out,
                          int num_splines, int num_walkers,
                          int row_stride)
{
  volatile __shared__ T in_shared[5][2*BS+1], out_shared[5][BS+1], kPoints_s[BS][3];
  __shared__ T  pos_s[3];
  __shared__ T *my_phi_in, *my_phi_out, *my_GL_in, *my_GL_out;
  int tid = threadIdx.x;
  if (tid == 0)
  {
    my_phi_in  = phi_in[blockIdx.x];
    my_phi_out = phi_out[blockIdx.x];
    my_GL_in   = grad_lapl_in[blockIdx.x];
    my_GL_out  = grad_lapl_out[blockIdx.x];
  }
  if (tid < 3)
    pos_s[tid] = pos[3*blockIdx.x+tid];
  //__syncthreads();
  int nb = (num_splines + BS-1)/BS;
  int outIndex=0;
  int outBlock=0;
  __shared__ int m2c[BS];
  for (int block=0; block<nb; block++)
  {
    // Load kpoints into shared memory
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*num_splines)
        kPoints_s[0][i*BS+tid] = kPoints[off];
    }
    // Load phi_in with coallesced reads
    if ((2*block+0)*BS+tid < 2*num_splines)
    {
      in_shared[0][tid+ 0] = my_phi_in[(2*block+0)*BS+tid];
      for (int j=0; j<4; j++)
        in_shared[j+1][tid+ 0] = my_GL_in[2*j*num_splines+(2*block+0)*BS+tid];
    }
    if ((2*block+1)*BS+tid < 2*num_splines)
    {
      in_shared[0][tid+BS] = my_phi_in[(2*block+1)*BS+tid];
      for (int j=0; j<4; j++)
        in_shared[j+1][tid+BS] = my_GL_in[2*j*num_splines+(2*block+1)*BS+tid];
    }
    //__syncthreads();
    // Now add on phase factors
    T phase = -(pos_s[0]*kPoints_s[tid][0] +
                pos_s[1]*kPoints_s[tid][1] +
                pos_s[2]*kPoints_s[tid][2]);
    T s, c;
    sincos (phase, &s, &c);
    T u_re, u_im, gradu_re[3], gradu_im[3], laplu_re, laplu_im;
    u_re        = in_shared[0][2*tid+0];
    u_im        = in_shared[0][2*tid+1];
    gradu_re[0] = in_shared[1][2*tid+0];
    gradu_im[0] = in_shared[1][2*tid+1];
    gradu_re[1] = in_shared[2][2*tid+0];
    gradu_im[1] = in_shared[2][2*tid+1];
    gradu_re[2] = in_shared[3][2*tid+0];
    gradu_im[2] = in_shared[3][2*tid+1];
    laplu_re    = in_shared[4][2*tid+0];
    laplu_im    = in_shared[4][2*tid+1];
    in_shared[0][2*tid+0] = u_re*c - u_im*s;
    in_shared[0][2*tid+1] = u_re*s + u_im*c;
    // Gradient = e^(-ikr)*(-i*u*k + gradu)
    for (int dim=0; dim<3; dim++)
    {
      T gre, gim;
      gre = gradu_re[dim] + kPoints_s[tid][dim]*u_im;
      gim = gradu_im[dim] - kPoints_s[tid][dim]*u_re;
      in_shared[dim+1][2*tid+0] = gre*c - gim*s;
      in_shared[dim+1][2*tid+1] = gre*s + gim*c;
    }
    // Add phase contribution to laplacian
    T k2 = (kPoints_s[tid][0]*kPoints_s[tid][0] +
            kPoints_s[tid][1]*kPoints_s[tid][1] +
            kPoints_s[tid][2]*kPoints_s[tid][2]);
    T lre = laplu_re - k2*u_re + 2.0*(kPoints_s[tid][0]*gradu_im[0]+
                                      kPoints_s[tid][1]*gradu_im[1]+
                                      kPoints_s[tid][2]*gradu_im[2]);
    T lim = laplu_im - k2*u_im - 2.0*(kPoints_s[tid][0]*gradu_re[0]+
                                      kPoints_s[tid][1]*gradu_re[1]+
                                      kPoints_s[tid][2]*gradu_re[2]);
    in_shared[4][2*tid+0] = lre*c - lim*s;
    in_shared[4][2*tid+1] = lre*s + lim*c;
    // Load makeTwoCopies with coallesced reads
    if (block*BS+tid < num_splines)
      m2c[tid] = makeTwoCopies[block*BS + tid];
    //__syncthreads();
    // Now, serialize to output buffer
    int end = min (BS, num_splines - block*BS);
    for (int i=0; i<end; i++)
    {
      if (tid < 5)
        out_shared[tid][outIndex] = in_shared[tid][2*i+0];
      outIndex++;
      //__syncthreads();
      if (outIndex == BS)
      {
        // Write back to global memory
        my_phi_out[             outBlock*BS+tid] = out_shared[0][tid];
        my_GL_out[0*row_stride +outBlock*BS+tid] = out_shared[1][tid];
        my_GL_out[1*row_stride +outBlock*BS+tid] = out_shared[2][tid];
        my_GL_out[2*row_stride +outBlock*BS+tid] = out_shared[3][tid];
        my_GL_out[3*row_stride +outBlock*BS+tid] = out_shared[4][tid];
        outIndex = 0;
        outBlock++;
      }
      if (m2c[i])
      {
        if (tid < 5)
          out_shared[tid][outIndex] = in_shared[tid][2*i+1];
        outIndex++;
        //__syncthreads();
        if (outIndex == BS)
        {
          // Write back to global memory
          my_phi_out[             outBlock*BS+tid] = out_shared[0][tid];
          my_GL_out[0*row_stride +outBlock*BS+tid] = out_shared[1][tid];
          my_GL_out[1*row_stride +outBlock*BS+tid] = out_shared[2][tid];
          my_GL_out[2*row_stride +outBlock*BS+tid] = out_shared[3][tid];
          my_GL_out[3*row_stride +outBlock*BS+tid] = out_shared[4][tid];
          outIndex = 0;
          outBlock++;
          //__syncthreads();
        }
      }
    }
    //__syncthreads();
  }
  if (tid < outIndex)
  {
    my_phi_out[             outBlock*BS+tid] = out_shared[0][tid];
    my_GL_out[0*row_stride +outBlock*BS+tid] = out_shared[1][tid];
    my_GL_out[1*row_stride +outBlock*BS+tid] = out_shared[2][tid];
    my_GL_out[2*row_stride +outBlock*BS+tid] = out_shared[3][tid];
    my_GL_out[3*row_stride +outBlock*BS+tid] = out_shared[4][tid];
  }
}

// Ye: optimized memory access.

template<typename T, int BS> __global__
void phase_factor_kernel (T *kPoints, int *makeTwoCopies,
                          int *TwoCopiesIndex,
                          T *pos, T **phi_in, T **phi_out,
                          T **grad_lapl_in, T **grad_lapl_out,
                          int num_splines, int num_walkers,
                          int row_stride)
{
  T in_shared[5][2], kPoints_s[3];
  __shared__ T  pos_s[3];
  __shared__ T *my_phi_in, *my_phi_out, *my_GL_in, *my_GL_out;
  int tid = threadIdx.x;
  if (tid == 0)
  {
    my_phi_in  = phi_in[blockIdx.x];
    my_phi_out = phi_out[blockIdx.x];
    my_GL_in   = grad_lapl_in[blockIdx.x];
    my_GL_out  = grad_lapl_out[blockIdx.x];
  }
  if (tid < 3)
    pos_s[tid] = pos[3*blockIdx.x+tid];
  __syncthreads();
  int nb = (num_splines + BS-1)/BS;
  __shared__ int m2c[BS], m2cIndex[BS];
  for (int block=0; block<nb; block++)
  {
    int off = block*BS + tid;
    if (off < num_splines)
    {
      // Load kpoints
      kPoints_s[0] = kPoints[off*3  ];
      kPoints_s[1] = kPoints[off*3+1];
      kPoints_s[2] = kPoints[off*3+2];
      // Load phi_in
      in_shared[0][0] = my_phi_in[off*2];
      in_shared[0][1] = my_phi_in[off*2+1];
      for (int j=0; j<4; j++)
      {
        in_shared[j+1][0] = my_GL_in[2*j*num_splines+2*off  ];
        in_shared[j+1][1] = my_GL_in[2*j*num_splines+2*off+1];
      }
      // Load makeTwoCopies
      m2c[tid] = makeTwoCopies[off];
      m2cIndex[tid] = TwoCopiesIndex[off];
    }
    // Now add on phase factors
    T phase = -(pos_s[0]*kPoints_s[0] +
                pos_s[1]*kPoints_s[1] +
                pos_s[2]*kPoints_s[2]);
    T s, c;
    sincos (phase, &s, &c);
    T u_re, u_im, gradu_re[3], gradu_im[3], laplu_re, laplu_im;
    u_re        = in_shared[0][0];
    u_im        = in_shared[0][1];
    gradu_re[0] = in_shared[1][0];
    gradu_im[0] = in_shared[1][1];
    gradu_re[1] = in_shared[2][0];
    gradu_im[1] = in_shared[2][1];
    gradu_re[2] = in_shared[3][0];
    gradu_im[2] = in_shared[3][1];
    laplu_re    = in_shared[4][0];
    laplu_im    = in_shared[4][1];
    in_shared[0][0] = u_re*c - u_im*s;
    in_shared[0][1] = u_re*s + u_im*c;
    // Gradient = e^(-ikr)*(-i*u*k + gradu)
    for (int dim=0; dim<3; dim++)
    {
      T gre, gim;
      gre = gradu_re[dim] + kPoints_s[dim]*u_im;
      gim = gradu_im[dim] - kPoints_s[dim]*u_re;
      in_shared[dim+1][0] = gre*c - gim*s;
      in_shared[dim+1][1] = gre*s + gim*c;
    }
    // Add phase contribution to laplacian
    T k2 = (kPoints_s[0]*kPoints_s[0] +
            kPoints_s[1]*kPoints_s[1] +
            kPoints_s[2]*kPoints_s[2]);
    T lre = laplu_re - k2*u_re + 2.0*(kPoints_s[0]*gradu_im[0]+
                                      kPoints_s[1]*gradu_im[1]+
                                      kPoints_s[2]*gradu_im[2]);
    T lim = laplu_im - k2*u_im - 2.0*(kPoints_s[0]*gradu_re[0]+
                                      kPoints_s[1]*gradu_re[1]+
                                      kPoints_s[2]*gradu_re[2]);
    in_shared[4][0] = lre*c - lim*s;
    in_shared[4][1] = lre*s + lim*c;
    // Now prepare the output buffer
    int end = min (BS, num_splines - block*BS);
    if ( tid < end )
    {
      my_phi_out[              m2cIndex[tid]] = in_shared[0][0];
      my_GL_out[0*row_stride + m2cIndex[tid]] = in_shared[1][0];
      my_GL_out[1*row_stride + m2cIndex[tid]] = in_shared[2][0];
      my_GL_out[2*row_stride + m2cIndex[tid]] = in_shared[3][0];
      my_GL_out[3*row_stride + m2cIndex[tid]] = in_shared[4][0];
      if (m2c[tid])
      {
        my_phi_out[              m2cIndex[tid]+1] = in_shared[0][1];
        my_GL_out[0*row_stride + m2cIndex[tid]+1] = in_shared[1][1];
        my_GL_out[1*row_stride + m2cIndex[tid]+1] = in_shared[2][1];
        my_GL_out[2*row_stride + m2cIndex[tid]+1] = in_shared[3][1];
        my_GL_out[3*row_stride + m2cIndex[tid]+1] = in_shared[4][1];
      }
    }
  }
}

#include <cstdio>
#include <complex>
#include <iostream>
#include "../CUDA/gpu_misc.h"

void apply_phase_factors(float kPoints[], int makeTwoCopies[],
                         float pos[], float *phi_in[], float *phi_out[],
                         int num_splines, int num_walkers)
{
//   float kPoints_h[3*num_splines];
//   int makeTwoCopies_h[num_splines];
//   float pos_h[3*num_walkers];
//   float *phi_in_ptr[num_walkers];
//   float *phi_out_ptr[num_walkers];
//   cudaMemcpy (kPoints_h, kPoints, 3*num_splines*sizeof(float), cudaMemcpyDeviceToHost);
//   cudaMemcpy (makeTwoCopies_h, makeTwoCopies, num_splines*sizeof(int), cudaMemcpyDeviceToHost);
//   cudaMemcpy (pos_h, pos, 3*num_walkers*sizeof(float), cudaMemcpyDeviceToHost);
//   cudaMemcpy (phi_in_ptr,  phi_in,  num_walkers*sizeof(float*), cudaMemcpyDeviceToHost);
//   cudaMemcpy (phi_out_ptr, phi_out, num_walkers*sizeof(float*), cudaMemcpyDeviceToHost);
//   for (int iw=0; iw<num_walkers; iw++) {
//     cudaMemcpy (kPoints_h, kPoints, 3*num_splines*sizeof(float), cudaMemcpyDeviceToHost);
//     std::complex<float> phi_in_h[num_splines];
//     float phi_out_h[num_splines*2];
//     cudaMemcpy (phi_in_h, phi_in_ptr[iw], num_splines*2*sizeof(float), cudaMemcpyDeviceToHost);
//     int iout = 0;
//     for (int isp=0; isp < num_splines; isp++) {
//       float phase = -(kPoints_h[3*isp+0] * pos_h[3*iw+0] +
// 		      kPoints_h[3*isp+1] * pos_h[3*iw+1] +
// 		      kPoints_h[3*isp+2] * pos_h[3*iw+2]);
//       float s,c;
//       sincosf(phase, &s, &c);
//       std::complex<float> z(c,s);
//       std::complex<float> out = z*phi_in_h[isp];
//       phi_out_h[iout++] = out.real();
//       if (makeTwoCopies_h[isp])
// 	phi_out_h[iout++] = out.imag();
//     }
//     cudaMemcpyAsync (phi_out_ptr[iw], phi_out_h, iout*sizeof(float), cudaMemcpyHostToDevice);
//   }
//   return;
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid ((num_walkers+BS-1)/BS);
  phase_factor_kernel<float,BS><<<dimGrid,dimBlock>>>
  (kPoints, makeTwoCopies, pos, phi_in, phi_out, num_splines, num_walkers);
  // dim3 dimGrid (num_walkers);
  // phase_factor_kernel_new<float,BS><<<dimGrid,dimBlock>>>
  //   (kPoints, makeTwoCopies, pos, phi_in, phi_out, num_splines);
}


void apply_phase_factors(float kPoints[], int makeTwoCopies[],
                         float pos[], float *phi_in[], float *phi_out[],
                         float *GL_in[], float *GL_out[],
                         int num_splines, int num_walkers, int row_stride)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);
  phase_factor_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, makeTwoCopies, pos, phi_in, phi_out,
   GL_in, GL_out, num_splines, num_walkers, row_stride);
}

void apply_phase_factors(float kPoints[], int makeTwoCopies[], int TwoCopiesIndex[],
                         float pos[], float *phi_in[], float *phi_out[],
                         float *GL_in[], float *GL_out[],
                         int num_splines, int num_walkers, int row_stride)
{
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);
  phase_factor_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, makeTwoCopies, TwoCopiesIndex, pos, phi_in, phi_out,
   GL_in, GL_out, num_splines, num_walkers, row_stride);
  /*
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);
  phase_factor_kernel<float,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, makeTwoCopies, pos, phi_in, phi_out,
   GL_in, GL_out, num_splines, num_walkers, row_stride);
  */
}

void apply_phase_factors(double kPoints[], int makeTwoCopies[],
                         double pos[], double *phi_in[], double *phi_out[],
                         int num_splines, int num_walkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid ((num_walkers+BS-1)/BS);
  phase_factor_kernel<double,BS><<<dimGrid,dimBlock>>>
  (kPoints, makeTwoCopies, pos, phi_in, phi_out, num_splines, num_walkers);
}


void apply_phase_factors(double kPoints[], int makeTwoCopies[],
                         double pos[], double *phi_in[], double *phi_out[],
                         double *GL_in[], double *GL_out[],
                         int num_splines, int num_walkers, int row_stride)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);
  phase_factor_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, makeTwoCopies, pos, phi_in, phi_out,
   GL_in, GL_out, num_splines, num_walkers, row_stride);
}


void apply_phase_factors(double kPoints[], int makeTwoCopies[], int TwoCopiesIndex[],
                         double pos[], double *phi_in[], double *phi_out[],
                         double *GL_in[], double *GL_out[],
                         int num_splines, int num_walkers, int row_stride)
{
  const int BS = 128;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);
  phase_factor_kernel<double,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, makeTwoCopies, TwoCopiesIndex, pos, phi_in, phi_out,
   GL_in, GL_out, num_splines, num_walkers, row_stride);
}


// YingWai's implementation for complex wavefunctions
// * The use of shared memory can be further optimized.
// * Thread count can be further optimized.
// Ye's optimization eliminates the shared memory buffer and potential race conditions.
// * BS>32 is safe but the optimal choice still needs further attempts.

#ifdef QMC_COMPLEX
// blockIdx.x = counter for walker
// BS threads per block
template<typename T, typename T2, int BS> __global__
void phase_factor_kernel (T *kPoints,
                          T *pos, T2 **phi_in, T2 **phi_out,
                          int num_splines)
{
  __shared__ T pos_s[3];
  __shared__ T2 *phi_in_ptr, *phi_out_ptr;

  int tid = threadIdx.x;

  // Copy electron position into shared memory
  if (tid < 3)
    pos_s[tid] = pos[3*blockIdx.x+tid];

  // Set pointers to point to the address (of the 1st element) where phi for that walker is stored
  if (tid == 0)
  {
    phi_in_ptr = phi_in[blockIdx.x];
    phi_out_ptr = phi_out[blockIdx.x];
  }

  // BS>warpSize sync needed 
  //__syncthreads();

  // "ib" is NOT the thread block! This is just a "block" of BS splines that the BS threads
  // process together at one time, so it has to repeat NB times to take care of all the splines
  int NB = (num_splines+BS-1)/BS;
  for (int ib=0; ib<NB; ib++)
  {
    int off = ib*BS + tid;
    if (off < num_splines)
    {
      // Now compute phase factors
      T phase = -(kPoints[off*3  ]*pos_s[0] +
                  kPoints[off*3+1]*pos_s[1] +
                  kPoints[off*3+2]*pos_s[2]);
      T s, c;
      sincos (phase, &s, &c);
      T2 e_mikr = T2(c,s);

      // Write back to global memory directly
      phi_out_ptr[ib*BS+tid] = phi_in_ptr[ib*BS+tid]*e_mikr;
    }
  }
}


template<typename T, typename T2, int BS> __global__
void phase_factor_kernel (T *kPoints,
                          T *pos, T2 **phi_in, T2 **phi_out,
                          T2 **grad_lapl_in, T2 **grad_lapl_out,
                          int num_splines, int num_walkers,
                          int row_stride)
{
  // Dummy quantities to make use of shared memory
  __shared__ T  pos_s[3];
  __shared__ T2 *my_phi_in, *my_phi_out, *my_GL_in, *my_GL_out;

  int tid = threadIdx.x;

  // Set pointers to point to the address (of the 1st element) where phi for that walker is stored
  if (tid == 0)
  {
    my_phi_in  = phi_in[blockIdx.x];
    my_phi_out = phi_out[blockIdx.x];
    my_GL_in   = grad_lapl_in[blockIdx.x];
    my_GL_out  = grad_lapl_out[blockIdx.x];
  }
  // Copy electron position into shared memory
  if (tid < 3)
    pos_s[tid] = pos[3*blockIdx.x+tid];

  __syncthreads();

  // "block" is NOT the thread block! This is just a "block" of BS splines that the BS threads
  // process together at one time, so it has to repeat nb times to take care of all the splines
  int nb = (num_splines + BS-1)/BS;
  for (int block=0; block<nb; block++)
  {
    T kPoints_s[3];
    int off = block*BS + tid;
    if (off < num_splines)
    {
      // Load kpoints
      kPoints_s[0] = kPoints[off*3  ];
      kPoints_s[1] = kPoints[off*3+1];
      kPoints_s[2] = kPoints[off*3+2];
      // Now compute phase factors
      T phase = -(kPoints_s[0]*pos_s[0] +
                  kPoints_s[1]*pos_s[1] +
                  kPoints_s[2]*pos_s[2]);
      T s, c;
      sincos (phase, &s, &c);
      T2 e_mikr = T2(c,s);

      // Temp. storage for phi, grad, lapl in the registers
      T2 u, gradlaplu[4];
      u = my_phi_in[off];
      for (int j=0; j<4; j++)
        gradlaplu[j] = my_GL_in[j*num_splines+off];

      // Apply e^(-ikr) to the wavefunction
      my_phi_out[off] = u * e_mikr;

      // Laplacian = e^(-ikr)*(laplu - u*k^2 - 2ik*gradu)
      T k2 = (kPoints_s[0]*kPoints_s[0] +
              kPoints_s[1]*kPoints_s[1] +
              kPoints_s[2]*kPoints_s[2]);
      gradlaplu[3] = e_mikr * ( gradlaplu[3] - k2*u
                               - T2(0.0,2.0)*(kPoints_s[0]*gradlaplu[0]+
                                              kPoints_s[1]*gradlaplu[1]+
                                              kPoints_s[2]*gradlaplu[2]) );

      // Gradient = e^(-ikr)*(-i*u*k + gradu)
      for (int dim=0; dim<3; dim++)
      {
        T2 g;
        g = gradlaplu[dim] - (T2(0.0,1.0) * u * kPoints_s[dim]);
        gradlaplu[dim] = g * e_mikr;
      }

      // Write back to global memory directly
      for (int j=0; j<4; j++)
        my_GL_out[j*row_stride + off] = gradlaplu[j];
    }
  }
}


// YingWai's implementation for complex wavefunction

void apply_phase_factors(float kPoints[], float pos[],
                         std::complex<float>* phi_in[], std::complex<float>* phi_out[],
                         int num_splines, int num_walkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);

  phase_factor_kernel<float,thrust::complex<float>,BS><<<dimGrid,dimBlock>>>
  (kPoints, pos, (thrust::complex<float>**)phi_in, (thrust::complex<float>**)phi_out, num_splines);
}

void apply_phase_factors(double kPoints[], double pos[],
                         std::complex<double>* phi_in[], std::complex<double>* phi_out[],
                         int num_splines, int num_walkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);

  phase_factor_kernel<double,thrust::complex<double>,BS><<<dimGrid,dimBlock>>>
  (kPoints, pos, (thrust::complex<double>**)phi_in, (thrust::complex<double>**)phi_out, num_splines);
}

void apply_phase_factors(float kPoints[], float pos[],
                         std::complex<float>* phi_in[], std::complex<float>* phi_out[],
                         std::complex<float>* GL_in[], std::complex<float>* GL_out[],
                         int num_splines, int num_walkers, int row_stride)
{

  const int BS = 128; 
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);

  phase_factor_kernel<float,thrust::complex<float>,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, pos, (thrust::complex<float>**)phi_in, (thrust::complex<float>**)phi_out, (thrust::complex<float>**)GL_in, (thrust::complex<float>**)GL_out, 
  num_splines, num_walkers, row_stride);
}


void apply_phase_factors(double kPoints[], double pos[],
                         std::complex<double>* phi_in[], std::complex<double>* phi_out[],
                         std::complex<double>* GL_in[], std::complex<double>* GL_out[],
                         int num_splines, int num_walkers, int row_stride)
{
  const int BS = 128; 
  dim3 dimBlock(BS);
  dim3 dimGrid (num_walkers);

  phase_factor_kernel<double,thrust::complex<double>,BS><<<dimGrid,dimBlock, 0, gpu::kernelStream>>>
  (kPoints, pos, (thrust::complex<double>**)phi_in, (thrust::complex<double>**)phi_out, (thrust::complex<double>**)GL_in, (thrust::complex<double>**)GL_out, 
  num_splines, num_walkers, row_stride);
}

#endif

