//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#define BLOCK_SIZE 64

#include "multi_bspline.h"
#include "multi_bspline_create_cuda.h"

//__constant__ float A[48];

// typedef struct
// {
//   float *coefs_real, *coefs_imag;
//   uint3 stride;
//   float3 gridInv;
//   int num_splines;
// } multi_UBspline_3d_c_cuda;

#ifndef NO_CUDA_MAIN
extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda (multi_UBspline_3d_c* spline)
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
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(A, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  multi_UBspline_3d_c_cuda *cuda_spline =
    (multi_UBspline_3d_c_cuda*) malloc (sizeof (multi_UBspline_3d_c_cuda*));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%BLOCK_SIZE) != 0)
    N += 64 - (N%BLOCK_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  size_t size = Nx*Ny*Nz+N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs_real), size);
  cudaMalloc((void**)&(cuda_spline->coefs_imag), size);
  
  float *spline_buff = (float*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp].real();
	}
  cudaMemcpy(cuda_spline->coefs_real, spline_buff, size, cudaMemcpyHostToDevice);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp].imag();
	}
  cudaMemcpy(cuda_spline->coefs_imag, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}
#endif


__global__ static void
eval_multi_multi_UBspline_3d_c_cuda (float *pos, float3 drInv, 
				     const float *coefs_real, const float *coefs_imag,
				     float *vals[], uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*BLOCK_SIZE+thr;

  __shared__ float *myval;
  __shared__ float abc[64], coefs[2*BLOCK_SIZE];

  // __shared__ float pos_s[BLOCK_SIZE];
  // int ir1 = (ir >> 4)*64;
  // int ir2 = (ir & 15)*4;
  // pos_s[thr] = pos[ir1+thr];
  // __syncthreads();
  // float3 r;
  // r.x = pos_s[ir2+0];
  // r.y = pos_s[ir2+1];
  // r.z = pos_s[ir2+2];
  __shared__ float3 r;
  if (thr == 0) {
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
  if (thr < 4) {
    a[thr] = Acuda[4*thr+0]*tp[0].x + Acuda[4*thr+1]*tp[0].y + Acuda[4*thr+2]*tp[0].z + Acuda[4*thr+3]*tp[0].w;
    b[thr] = Acuda[4*thr+0]*tp[1].x + Acuda[4*thr+1]*tp[1].y + Acuda[4*thr+2]*tp[1].z + Acuda[4*thr+3]*tp[1].w;
    c[thr] = Acuda[4*thr+0]*tp[2].x + Acuda[4*thr+1]*tp[2].y + Acuda[4*thr+2]*tp[2].z + Acuda[4*thr+3]*tp[2].w;
  }
  __syncthreads();

  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);
  
  abc[thr] = a[i]*b[j]*c[k];
  __syncthreads();


  float val_real = 0.0;
  float val_imag = 0.0;
  val_real = val_imag = 0.0;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      for (int k=0; k<4; k++) {
  	val_real += abc[16*i+4*j+k] * base_real[off+k*strides.z];
  	val_imag += abc[16*i+4*j+k] * base_imag[off+k*strides.z];
      }
    }
  }
  // for (int i=0; i<4; i++) {
  //   for (int j=0; j<4; j++) {
  //     float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
  //     float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
  //     for (int k=0; k<4; k++) {
  // 	coefs[thr]            = base_real[(2*block+0)*BLOCK_SIZE+thr];
  // 	coefs[thr+BLOCK_SIZE] = base_real[(2*block+1)*BLOCK_SIZE+thr];
  // 	__syncthreads();
  // 	val_real += abc[16*i+4*j+k] * coefs[2*thr+0];
  // 	val_imag += abc[16*i+4*j+k] * coefs[2*thr+1];
  //     }
  //   }
  // }
  __shared__ float buff[2*BLOCK_SIZE];
  buff[2*thr+0] = val_real;
  buff[2*thr+1] = val_imag;
  __syncthreads();
  myval[off] = buff[thr];
  myval[off+BLOCK_SIZE] = buff[thr+BLOCK_SIZE];

//   myval[2*off+0] = val_real;
//   myval[2*off+1] = val_imag;
  //myval[off+BLOCK_SIZE] = val_imag;
  //vals_real[ir][offset] = val_real;
  //vals_imag[ir][offset] = val_imag;
}



__global__ static void
eval_multi_multi_UBspline_3d_c_vgh_cuda (float *pos, float3 drInv, 
					 const float *coefs_real, const float *coefs_imag,
					 float *vals[], float *grads[], float *hess[],
					 uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*BLOCK_SIZE+thr;

  __shared__ float *myval, *mygrad, *myhess;
  __shared__ float3 r;
  if (thr == 0) {
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
  if (thr < 12) {
    a[thr] = Acuda[4*thr+0]*tp[0].x + Acuda[4*thr+1]*tp[0].y + Acuda[4*thr+2]*tp[0].z + Acuda[4*thr+3]*tp[0].w;
    b[thr] = Acuda[4*thr+0]*tp[1].x + Acuda[4*thr+1]*tp[1].y + Acuda[4*thr+2]*tp[1].z + Acuda[4*thr+3]*tp[1].w;
    c[thr] = Acuda[4*thr+0]*tp[2].x + Acuda[4*thr+1]*tp[2].y + Acuda[4*thr+2]*tp[2].z + Acuda[4*thr+3]*tp[2].w;
  }
  __syncthreads();

  __shared__ float abc[640];
  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);

  abc[10*(16*i+4*j+k)+0] = a[i+0]*b[j+0]*c[k+0]; // val
  abc[10*(16*i+4*j+k)+1] = a[i+4]*b[j+0]*c[k+0]; // d/dx
  abc[10*(16*i+4*j+k)+2] = a[i+0]*b[j+4]*c[k+0]; // d/dy
  abc[10*(16*i+4*j+k)+3] = a[i+0]*b[j+0]*c[k+4]; // d/dz
  abc[10*(16*i+4*j+k)+4] = a[i+8]*b[j+0]*c[k+0]; // d2/dx2
  abc[10*(16*i+4*j+k)+5] = a[i+4]*b[j+4]*c[k+0]; // d2/dxdy
  abc[10*(16*i+4*j+k)+6] = a[i+4]*b[j+0]*c[k+4]; // d2/dxdz
  abc[10*(16*i+4*j+k)+7] = a[i+0]*b[j+8]*c[k+0]; // d2/dy2
  abc[10*(16*i+4*j+k)+8] = a[i+0]*b[j+4]*c[k+4]; // d2/dydz
  abc[10*(16*i+4*j+k)+9] = a[i+0]*b[j+0]*c[k+8]; // d2/dz2

  __syncthreads();

  float v_r = 0.0;
  float v_i = 0.0;
  float g0_r=0.0, g0_i=0.0, g1_r=0.0, g1_i=0.0, g2_r=0.0, g2_i=0.0, 
    h00_r=0.0, h00_i=0.0, h01_r=0.0, h01_i=0.0, h02_r=0.0, h02_i=0.0, 
    h11_r=0.0, h11_i=0.0, h12_r=0.0, h12_i=0.0, h22_r=0.0, h22_i=0.0;
  int n = 0;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
//       float c0_r, c0_i, c1_r, c1_i, c2_r, c2_i, c3_r, c3_i;
//       c0_r = base_real[off+0*strides.z];  c0_i = base_imag[off+0*strides.z];
//       c1_r = base_real[off+1*strides.z];  c1_i = base_imag[off+1*strides.z];
//       c2_r = base_real[off+2*strides.z];  c2_i = base_imag[off+2*strides.z];
//       c3_r = base_real[off+3*strides.z];  c3_i = base_imag[off+3*strides.z];

//       v_r   += abc[n+0] * c0_r;  v_i   += abc[n+0] * c0_i;
//       g0_r  += abc[n+1] * c0_r;  g0_i  += abc[n+1] * c0_i;
//       g1_r  += abc[n+2] * c0_r;  g1_i  += abc[n+2] * c0_i;
//       g2_r  += abc[n+3] * c0_r;  g2_i  += abc[n+3] * c0_i;
//       h00_r += abc[n+4] * c0_r;  h00_i += abc[n+4] * c0_i;
//       h01_r += abc[n+5] * c0_r;  h01_i += abc[n+5] * c0_i;
//       h02_r += abc[n+6] * c0_r;  h02_i += abc[n+6] * c0_i;
//       h11_r += abc[n+7] * c0_r;  h11_i += abc[n+7] * c0_i;
//       h12_r += abc[n+8] * c0_r;  h12_i += abc[n+8] * c0_i;
//       h22_r += abc[n+9] * c0_r;  h22_i += abc[n+9] * c0_i;       

//       v_r   += abc[n+10] * c1_r;  v_i   += abc[n+10] * c1_i;
//       g0_r  += abc[n+11] * c1_r;  g0_i  += abc[n+11] * c1_i;
//       g1_r  += abc[n+12] * c1_r;  g1_i  += abc[n+12] * c1_i;
//       g2_r  += abc[n+13] * c1_r;  g2_i  += abc[n+13] * c1_i;
//       h00_r += abc[n+14] * c1_r;  h00_i += abc[n+14] * c1_i;
//       h01_r += abc[n+15] * c1_r;  h01_i += abc[n+15] * c1_i;
//       h02_r += abc[n+16] * c1_r;  h02_i += abc[n+16] * c1_i;
//       h11_r += abc[n+17] * c1_r;  h11_i += abc[n+17] * c1_i;
//       h12_r += abc[n+18] * c1_r;  h12_i += abc[n+18] * c1_i;
//       h22_r += abc[n+19] * c1_r;  h22_i += abc[n+19] * c1_i;       

//       v_r   += abc[n+20] * c2_r;  v_i   += abc[n+20] * c2_i;
//       g0_r  += abc[n+21] * c2_r;  g0_i  += abc[n+21] * c2_i;
//       g1_r  += abc[n+22] * c2_r;  g1_i  += abc[n+22] * c2_i;
//       g2_r  += abc[n+23] * c2_r;  g2_i  += abc[n+23] * c2_i;
//       h00_r += abc[n+24] * c2_r;  h00_i += abc[n+24] * c2_i;
//       h01_r += abc[n+25] * c2_r;  h01_i += abc[n+25] * c2_i;
//       h02_r += abc[n+26] * c2_r;  h02_i += abc[n+26] * c2_i;
//       h11_r += abc[n+27] * c2_r;  h11_i += abc[n+27] * c2_i;
//       h12_r += abc[n+28] * c2_r;  h12_i += abc[n+28] * c2_i;
//       h22_r += abc[n+29] * c2_r;  h22_i += abc[n+29] * c2_i;       

//       v_r   += abc[n+30] * c3_r;  v_i   += abc[n+30] * c3_i;
//       g0_r  += abc[n+31] * c3_r;  g0_i  += abc[n+31] * c3_i;
//       g1_r  += abc[n+32] * c3_r;  g1_i  += abc[n+32] * c3_i;
//       g2_r  += abc[n+33] * c3_r;  g2_i  += abc[n+33] * c3_i;
//       h00_r += abc[n+34] * c3_r;  h00_i += abc[n+34] * c3_i;
//       h01_r += abc[n+35] * c3_r;  h01_i += abc[n+35] * c3_i;
//       h02_r += abc[n+36] * c3_r;  h02_i += abc[n+36] * c3_i;
//       h11_r += abc[n+37] * c3_r;  h11_i += abc[n+37] * c3_i;
//       h12_r += abc[n+38] * c3_r;  h12_i += abc[n+38] * c3_i;
//       h22_r += abc[n+39] * c3_r;  h22_i += abc[n+39] * c3_i;       
//       n += 40;

      for (int k=0; k<4; k++) {
	float cr = base_real[off+k*strides.z];
	float ci = base_imag[off+k*strides.z];
	v_r   += abc[n+0] * cr;  v_i   += abc[n+0] * ci;
	g0_r  += abc[n+1] * cr;  g0_i  += abc[n+1] * ci;
	g1_r  += abc[n+2] * cr;  g1_i  += abc[n+2] * ci;
	g2_r  += abc[n+3] * cr;  g2_i  += abc[n+3] * ci;
	h00_r += abc[n+4] * cr;  h00_i += abc[n+4] * ci;
	h01_r += abc[n+5] * cr;  h01_i += abc[n+5] * ci;
	h02_r += abc[n+6] * cr;  h02_i += abc[n+6] * ci;
	h11_r += abc[n+7] * cr;  h11_i += abc[n+7] * ci;
	h12_r += abc[n+8] * cr;  h12_i += abc[n+8] * ci;
	h22_r += abc[n+9] * cr;  h22_i += abc[n+9] * ci; 
	n += 10;
      }
    }
  }
  g0_r *= drInv.x; g0_i *= drInv.x;
  g1_r *= drInv.y; g1_i *= drInv.y;
  g2_r *= drInv.z; g2_i *= drInv.z;

  h00_r *= drInv.x * drInv.x;  h00_i *= drInv.x * drInv.x;
  h01_r *= drInv.x * drInv.y;  h01_i *= drInv.x * drInv.y;
  h02_r *= drInv.x * drInv.z;  h02_i *= drInv.x * drInv.z;
  h11_r *= drInv.y * drInv.y;  h11_i *= drInv.y * drInv.y;
  h12_r *= drInv.y * drInv.z;  h12_i *= drInv.y * drInv.z;
  h22_r *= drInv.z * drInv.z;  h22_i *= drInv.z * drInv.z;

  
  __shared__ float buff[6*BLOCK_SIZE];
  // Note, we can reuse abc, by replacing buff with abc.
  
  buff[2*thr+0] = v_r;  buff[2*thr+1] = v_i;
  __syncthreads();
  myval[off] = buff[thr];    
  myval[off+BLOCK_SIZE] = buff[thr+BLOCK_SIZE];

  buff[6*thr+0] = g0_r;  buff[6*thr+1] = g0_i;
  buff[6*thr+2] = g1_r;  buff[6*thr+3] = g1_i;
  buff[6*thr+4] = g2_r;  buff[6*thr+5] = g2_i;
  __syncthreads();
  for (int i=0; i<6; i++) 
    mygrad[(6*block+i)*BLOCK_SIZE+thr] = buff[i*BLOCK_SIZE+thr]; 
  __syncthreads();

  // Write first half of Hessians
  if (thr < 32) {
    buff[12*thr+0]  = h00_r;    buff[12*thr+1]  = h00_i;
    buff[12*thr+2]  = h01_r;    buff[12*thr+3]  = h01_i;
    buff[12*thr+4]  = h02_r;    buff[12*thr+5]  = h02_i;
    buff[12*thr+6]  = h11_r;    buff[12*thr+7]  = h11_i;
    buff[12*thr+8]  = h12_r;    buff[12*thr+9]  = h12_i;
    buff[12*thr+10] = h22_r;    buff[12*thr+11] = h22_i;
  }
  __syncthreads();
  if (thr < 32) 
    for (int i=0; i<6; i++) 
      myhess[(12*block+i)*BLOCK_SIZE+thr] = buff[i*BLOCK_SIZE+thr];

  __syncthreads();
  int th2 = thr-32;
  if (thr >= 32) {
    buff[12*th2+0]  = h00_r;    buff[12*th2+1]  = h00_i;
    buff[12*th2+2]  = h01_r;    buff[12*th2+3]  = h01_i;
    buff[12*th2+4]  = h02_r;    buff[12*th2+5]  = h02_i;
    buff[12*th2+6]  = h11_r;    buff[12*th2+7]  = h11_i;
    buff[12*th2+8]  = h12_r;    buff[12*th2+9]  = h12_i;
    buff[12*th2+10] = h22_r;    buff[12*th2+11] = h22_i;
  }
  __syncthreads();
  if (thr >= 32) {
    for (int i=0; i<6; i++) 
      myhess[(12*block+i+6)*BLOCK_SIZE+th2] = buff[i*BLOCK_SIZE+th2];
  }

}

				    

#ifndef NO_CUDA_MAIN
static void *
test_multi_cuda(void *thread)
{
//   CUcontext ctx;
//   CUdevice dev;
//   cuDeviceGet (&dev, (int)(size_t)thread);
//   cuCtxCreate(&ctx, CU_CTX_SCHED_YIELD, dev);

//   int deviceCount;
//   cudaGetDeviceCount(&deviceCount);

  cudaSetDevice((int)(size_t)thread);
  fprintf (stderr, "In thread %p\n", thread);

  int numWalkers = 200;
  float *coefs  ,  __device__ *vals[numWalkers], *grads[numWalkers], *hess[numWalkers];
  float *coefs_real_d, *coefs_imag_d, __device__ **vals_d, **grads_d, **hess_d;
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
		         0.0,      0.0,      1.0,     0.0 };

  // Copy A to host
  cudaMemcpy(Acuda, A_h, 48*sizeof(float), cudaMemcpyHostToDevice); 

  float *r_d, *r_h;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 128;
  Nx = Ny = Nz = 16;
  xs = Ny*Nz*N;
  ys = Nz*N;
  zs = N;

  float3 drInv;
  drInv.x = 1.0/float(Nx);
  drInv.y = 1.0/float(Ny);
  drInv.z = 1.0/float(Nz);

  // Setup Bspline coefficients
  int size = Nx*Ny*Nz*N*sizeof(float);
  posix_memalign((void**)&coefs, 16, size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int n=0; n<N; n++)
	  coefs[ix*xs + iy*ys + iz*zs + n] = drand48();


  fprintf (stderr, "Filled in coefs.\n");

  // Setup values
  //posix_memalign((void**)&vals, 16, N*sizeof(float));

  // cudaMemcpy(r_d, r, numWalkers*sizeof(float3), cudaMemcpyHostToDevice);

  
  fprintf (stderr, "size = %d\n", size);
  
  // Setup CUDA coefficients
  fprintf (stderr, "Before first CUDA mallocs.\n");
  cudaMalloc((void**)&coefs_real_d, 2*size);
  cudaMalloc((void**)&coefs_imag_d, 2*size);
  fprintf (stderr, "Before Memcpy.\n");
  cudaMemcpy(coefs_real_d, coefs, size, cudaMemcpyHostToDevice);
  cudaMemcpy(coefs_imag_d, coefs, size, cudaMemcpyHostToDevice);
  fprintf (stderr, "After Memcpy.\n");  

  // Setup device value storage
  int numVals = 2*N*numWalkers*10;
  float *valBlock_d, *valBlock_h;
  cudaMalloc((void**)&(valBlock_d),     numVals*sizeof(float));
  cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(float));
  cudaMalloc((void**)&(vals_d), 2*numWalkers*sizeof(float*));
  cudaMalloc((void**)&(grads_d), 2*numWalkers*sizeof(float*));
  cudaMalloc((void**)&(hess_d), 2*numWalkers*sizeof(float*));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals[i]  = valBlock_d + 2*i*N;
    grads[i] = valBlock_d + 2*N*numWalkers + 6*i*N;
    hess[i]  = valBlock_d + 8*N*numWalkers + 12*i*N;
  }
  cudaMemcpy(vals_d,  vals,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(grads_d, grads, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(hess_d,  hess,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  
  fprintf (stderr, "Finished cuda allocations.\n");


  // Setup walker positions
  cudaMalloc((void**)&(r_d),     4*numWalkers*sizeof(float));
  cudaMallocHost((void**)&(r_h), 4*numWalkers*sizeof(float));

  for (int ir=0; ir<numWalkers; ir++) {
    r_h[4*ir+0] = 0.5*drand48();
    r_h[4*ir+1] = 0.5*drand48();
    r_h[4*ir+2] = 0.5*drand48();
  }

  uint3 strides;
  strides.x = xs;
  strides.y = ys;
  strides.z = zs;

  dim3 dimBlock(BLOCK_SIZE);
  dim3 dimGrid(N/BLOCK_SIZE,numWalkers);
  
  clock_t start, end;

  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 4*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_c_cuda<<<dimGrid,dimBlock>>> 
       (r_d, drInv, coefs_real_d, coefs_imag_d, 
        vals_d, strides);
    // eval_multi_multi_UBspline_3d_cuda_c<<<dimGrid,dimBlock>>> 
    //   (r_d, drInv, coefs_real_d, coefs_imag_d, 
    //    valBlock_d, valBlock_d+numVals/2, strides);
    //cudaMemcpy(valBlock_h, valBlock_d, numVals*sizeof(float), cudaMemcpyDeviceToHost);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "VGH evals per second = %1.8e\n", 1.0/time);


  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 4*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_c_vgh_cuda<<<dimGrid,dimBlock>>> 
       (r_d, drInv, coefs_real_d, coefs_imag_d, 
        vals_d, grads_d, hess_d, strides);
  }
  end = clock();
  time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (valBlock_d);
  cudaFree (vals_d);
  cudaFree (coefs_real_d);
  cudaFree (coefs_imag_d);
  cudaFree (r_d);

  return NULL;

  // cudaMemcpy (vals, vals_d, N*sizeof(float), cudaMemcpyDeviceToHost);

  // float vals2[N];
  
  // for (int n=0; n<N; n++) {
  //   vals2[n] = 0.0;
  //   int index=0;
  //   for(int i=0; i<4; i++)
  //     for (int j=0; j<4; j++)
  // 	for (int k=0; k<4; k++)  {
  // 	  vals2[n] += abc[index] * coefs[(ix+i)*xs+(iy+j)*ys+(iz+k)*zs+n];
  // 	  index++;
  // 	}
  // }


  // for (int i=0; i<N/256; i++)	
  //   fprintf (stderr, "%1.9f %1.9f\n", vals[i], vals2[i]); 


  // cudaFree(abc_d);
  // cudaFree(coefs_d);
  // cudaFree(vals_d);
}
#endif

#ifndef NO_CUDA_MAIN

main()
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  fprintf (stderr, "Detected %d CUDA devices.\n", deviceCount);

  // test_cuda();

  for (int device = 0; device < deviceCount; ++device) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    fprintf (stderr, "Device %d:\n", device);
    fprintf (stderr, "  Global memory:     %10d\n",
	     deviceProp.totalGlobalMem);
    fprintf (stderr, "  MultiProcessors:   %10d\n",
	     deviceProp.multiProcessorCount);
    fprintf (stderr, "  Registers:         %10d\n", 
	     deviceProp.regsPerBlock);
    fprintf (stderr, "  Constant memory:   %10d\n", 
	     deviceProp.totalConstMem);
    fprintf (stderr, "  Shared memory:     %10d\n", 
	     deviceProp.sharedMemPerBlock);
  }

  //  pthread_t threads[deviceCount];

  // for (int device = 0; device < deviceCount; device++) 
  //   pthread_create (&(threads[device]), NULL, test_multi_cuda, (void*)device);
  // cutStartThread((CUT_THREADROUTINE)test_multi_cuda,(void*)device);
  test_multi_cuda((void*)0);

  //  pthread_exit(NULL);
  //test_multi_cuda();
}

#endif
