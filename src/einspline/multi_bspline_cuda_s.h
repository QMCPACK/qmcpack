//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef MULTI_BSPLINE_CUDA_S_H
#define MULTI_BSPLINE_CUDA_S_H

#include "multi_bspline_structs_cuda.h"


__global__ static void
eval_multi_multi_UBspline_3d_s_cuda (float *pos, float3 drInv,
                                     const float *coefs, float *vals[], uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+thr;
  __shared__ float *myval;
  __shared__ float abc[64];
  __shared__ float3 r;
  if (thr == 0)
  {
    r.x = pos[4*ir+0];
    r.y = pos[4*ir+1];
    r.z = pos[4*ir+2];
    myval = vals[ir];
  }
  __syncthreads();
  int3 index;
  float3 t;
  float s, sf;
  float4 tp[3];
  s = r.x * drInv.x;
  sf = floor(s);
  index.x = (int)sf;
  t.x = s - sf;
  s = r.y * drInv.y;
  sf = floor(s);
  index.y = (int)sf;
  t.y = s - sf;
  s = r.z * drInv.z;
  sf = floor(s);
  index.z = (int)sf;
  t.z = s - sf;
  tp[0] = make_float4(t.x*t.x*t.x, t.x*t.x, t.x, 1.0);
  tp[1] = make_float4(t.y*t.y*t.y, t.y*t.y, t.y, 1.0);
  tp[2] = make_float4(t.z*t.z*t.z, t.z*t.z, t.z, 1.0);
  __shared__ float a[4], b[4], c[4];
  if (thr < 4)
  {
    a[thr] = A[4*thr+0]*tp[0].x + A[4*thr+1]*tp[0].y + A[4*thr+2]*tp[0].z + A[4*thr+3]*tp[0].w;
    b[thr] = A[4*thr+0]*tp[1].x + A[4*thr+1]*tp[1].y + A[4*thr+2]*tp[1].z + A[4*thr+3]*tp[1].w;
    c[thr] = A[4*thr+0]*tp[2].x + A[4*thr+1]*tp[2].y + A[4*thr+2]*tp[2].z + A[4*thr+3]*tp[2].w;
  }
  __syncthreads();
  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);
  if (thr < 64)
    abc[thr] = a[i]*b[j]*c[k];
  __syncthreads();
  float val = 0.0;
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    {
      float *base = coefs + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      for (int k=0; k<4; k++)
        val += abc[16*i+4*j+k] * base[off+k*strides.z];
    }
  }
  myval[off] = val;
}



__global__ static void
eval_multi_multi_UBspline_3d_s_vgh_cuda (float *pos, float3 drInv,  const float *coefs,
    float *vals[], float *grads[], float *hess[],
    uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+threadIdx.x;
  __shared__ float *myval, *mygrad, *myhess;
  __shared__ float3 r;
  if (thr == 0)
  {
    r.x = pos[4*ir+0];
    r.y = pos[4*ir+1];
    r.z = pos[4*ir+2];
    myval  = vals[ir];
    mygrad = grads[ir];
    myhess = hess[ir];
  }
  __syncthreads();
  int3 index;
  float3 t;
  float s, sf;
  float4 tp[3];
  s = r.x * drInv.x;
  sf = floor(s);
  index.x = (int)sf;
  t.x = s - sf;
  s = r.y * drInv.y;
  sf = floor(s);
  index.y = (int)sf;
  t.y = s - sf;
  s = r.z * drInv.z;
  sf = floor(s);
  index.z = (int)sf;
  t.z = s - sf;
  tp[0] = make_float4(t.x*t.x*t.x, t.x*t.x, t.x, 1.0);
  tp[1] = make_float4(t.y*t.y*t.y, t.y*t.y, t.y, 1.0);
  tp[2] = make_float4(t.z*t.z*t.z, t.z*t.z, t.z, 1.0);
  // First 4 of a are value, second 4 are derivative, last four are
  // second derivative.
  __shared__ float a[12], b[12], c[12];
  if (thr < 12)
  {
    a[thr] = A[4*thr+0]*tp[0].x + A[4*thr+1]*tp[0].y + A[4*thr+2]*tp[0].z + A[4*thr+3]*tp[0].z;
    b[thr] = A[4*thr+0]*tp[1].x + A[4*thr+1]*tp[1].y + A[4*thr+2]*tp[1].z + A[4*thr+3]*tp[1].z;
    c[thr] = A[4*thr+0]*tp[2].x + A[4*thr+1]*tp[2].y + A[4*thr+2]*tp[2].z + A[4*thr+3]*tp[2].z;
  }
  __syncthreads();
  __shared__ float abc[640];
  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);
  abc[(16*i+4*j+k)+0]   = a[i+0]*b[j+0]*c[k+0]; // val
  abc[(16*i+4*j+k)+64]  = a[i+4]*b[j+0]*c[k+0]; // d/dx
  abc[(16*i+4*j+k)+128] = a[i+0]*b[j+4]*c[k+0]; // d/dy
  abc[(16*i+4*j+k)+192] = a[i+0]*b[j+0]*c[k+4]; // d/dz
  abc[(16*i+4*j+k)+256] = a[i+8]*b[j+0]*c[k+0]; // d2/dx2
  abc[(16*i+4*j+k)+320] = a[i+4]*b[j+4]*c[k+0]; // d2/dxdy
  abc[(16*i+4*j+k)+384] = a[i+4]*b[j+0]*c[k+4]; // d2/dxdz
  abc[(16*i+4*j+k)+448] = a[i+0]*b[j+8]*c[k+0]; // d2/dy2
  abc[(16*i+4*j+k)+512] = a[i+0]*b[j+4]*c[k+4]; // d2/dydz
  abc[(16*i+4*j+k)+576] = a[i+0]*b[j+0]*c[k+8]; // d2/dz2
  __syncthreads();
  float v = 0.0, g0=0.0,  g1=0.0, g2=0.0,
        h00=0.0, h01=0.0, h02=0.0, h11=0.0, h12=0.0, h22=0.0;
  int n = 0;
  float *b0 = coefs + index.x*strides.x + index.y*strides.y + index.z*strides.z + off;
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    {
      float *base = b0 + i*strides.x + j*strides.y;
      for (int k=0; k<4; k++)
      {
        float c  = base[k*strides.z];
        v   += abc[n+0] * c;
        g0  += abc[n+1] * c;
        g1  += abc[n+2] * c;
        g2  += abc[n+3] * c;
        h00 += abc[n+4] * c;
        h01 += abc[n+5] * c;
        h02 += abc[n+6] * c;
        h11 += abc[n+7] * c;
        h12 += abc[n+8] * c;
        h22 += abc[n+9] * c;
        n += 10;
      }
    }
  }
  g0 *= drInv.x;
  g1 *= drInv.y;
  g2 *= drInv.z;
  h00 *= drInv.x * drInv.x;
  h01 *= drInv.x * drInv.y;
  h02 *= drInv.x * drInv.z;
  h11 *= drInv.y * drInv.y;
  h12 *= drInv.y * drInv.z;
  h22 *= drInv.z * drInv.z;
  //  __shared__ float buff[6*SPLINE_BLOCK_SIZE];
  // Note, we can reuse abc, by replacing buff with abc.
  myval[off] = v;
  abc[3*thr+0] = g0;
  abc[3*thr+1] = g1;
  abc[3*thr+2] = g2;
  __syncthreads();
  for (int i=0; i<3; i++)
    mygrad[(3*block+i)*SPLINE_BLOCK_SIZE+thr] = abc[i*SPLINE_BLOCK_SIZE+thr];
  __syncthreads();
  // Write first half of Hessians
  abc[6*thr+0]  = h00;
  abc[6*thr+1]  = h01;
  abc[6*thr+2]  = h02;
  abc[6*thr+3]  = h11;
  abc[6*thr+4]  = h12;
  abc[6*thr+5]  = h22;
  __syncthreads();
  for (int i=0; i<6; i++)
    myhess[(6*block+i)*SPLINE_BLOCK_SIZE+thr] = abc[i*SPLINE_BLOCK_SIZE+thr];
}


#endif
