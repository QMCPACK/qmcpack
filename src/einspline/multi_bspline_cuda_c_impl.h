#ifndef MULTI_BSPLINE_CUDA_C_IMPL_H
#define MULTI_BSPLINE_CUDA_C_IMPL_H

#include "multi_bspline.h"
#include "multi_bspline_create_cuda.h"


__global__ static void
eval_multi_multi_UBspline_1d_c_kernel
(float *pos, float drInv, const float *coefs, float **vals, 
 uint dim, uint stride, int N)
{
  int tid   = threadIdx.x;
  int ir    = blockIdx.x;

  __shared__ float *ourval;
  __shared__ float r;
  if (tid == 0) {
    r = pos[ir];
    ourval = vals[ir];
  }
  __syncthreads();
  
  int index;
  float t;
  float s, sf;
  float4 tp;

  s = r * drInv;
  sf = floor(s);
  index = min(max(0,(int)sf), dim-1);
  t = s - sf;
  tp = make_float4(t*t*t, t*t, t, 1.0);

  __shared__ float a[4];
  if (tid < 4) 
    a[tid] = Acuda[4*tid+0]*tp.x + Acuda[4*tid+1]*tp.y + Acuda[4*tid+2]*tp.z + Acuda[4*tid+3]*tp.w;
  __syncthreads();

  int numBlocks = 2*N / SPLINE_BLOCK_SIZE;
  const float *c = coefs + index*stride + tid;
  float *myval = ourval + tid;
  int stride2 = 2*stride;
  int stride3 = 3*stride;
  for (int block=0; block < numBlocks; block++) {
     *myval = (a[0] * c[0] +
	       a[1] * c[stride] +
	       a[2] * c[stride2] +
	       a[3] * c[stride3]);
     myval += SPLINE_BLOCK_SIZE;    c += SPLINE_BLOCK_SIZE;
  }
      
  int remainder = 2*N - numBlocks*SPLINE_BLOCK_SIZE;
  if (tid < remainder) {
    *myval = (a[0] * c[0] +
	      a[1] * c[stride] +
	      a[2] * c[stride2] +
	      a[3] * c[stride3]);
  }
}

extern "C" void
eval_multi_multi_UBspline_1d_c_cuda (const multi_UBspline_1d_c_cuda *spline,
				     float *pos_d, float *vals_d[], int num)
{
  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(num);

  eval_multi_multi_UBspline_1d_c_kernel<<<dimGrid,dimBlock>>>
    (pos_d, spline->gridInv, (float*)spline->coefs, vals_d, spline->dim, spline->stride, spline->num_splines);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in eval_multi_multi_UBspline_1d_c_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}





__global__ static void
eval_multi_multi_UBspline_1d_c_vgl_kernel
(float *pos, float drInv, const float *coefs, float **vals, 
float **grads, float **lapl, 
 uint dim, uint stride, int N)
{
  int tid   = threadIdx.x;
  int ir    = blockIdx.x;

  __shared__ float *ourval, *ourgrad, *ourlapl;
  __shared__ float r;
  if (tid == 0) {
    r = pos[ir];
    ourval = vals[ir];
    ourgrad = grads[ir];
    ourlapl = lapl[ir];
  }
  __syncthreads();
  
  int index;
  float t;
  float s, sf;
  float4 tp;

  s = r * drInv;
  sf = floor(s);
  index = min(max(0,(int)sf), dim-1);
  t = s - sf;
  tp = make_float4(t*t*t, t*t, t, 1.0);

  __shared__ float a[12];
  if (tid < 12) 
    a[tid] = Acuda[4*tid+0]*tp.x + Acuda[4*tid+1]*tp.y + Acuda[4*tid+2]*tp.z + Acuda[4*tid+3]*tp.w;
  __syncthreads();

  int numBlocks = 2*N / SPLINE_BLOCK_SIZE;
  const float *c = coefs + index*stride + tid;
  float *myval  = ourval + tid;
  float *mygrad = ourgrad + tid;
  float *mylapl = ourlapl + tid;
  int stride2 = 2*stride;
  int stride3 = 3*stride;
  __shared__ float coef[SPLINE_BLOCK_SIZE][5];
  for (int block=0; block < numBlocks; block++) {
    coef[tid][0] = c[0];
    coef[tid][1] = c[stride];
    coef[tid][2] = c[stride2];
    coef[tid][3] = c[stride3];
    *myval = (a[0] * coef[tid][0]   + a[1] * coef[tid][1] +
	      a[2] * coef[tid][2]   + a[3] * coef[tid][3]);
    *mygrad = (a[4] * coef[tid][0]  + a[5] * coef[tid][1] +
	       a[6] * coef[tid][2]  + a[7] * coef[tid][3]);
    *mylapl = (a[8] * coef[tid][0]  + a[9] * coef[tid][1] +
	       a[10] * coef[tid][2] + a[11]* coef[tid][3]);
    myval  += SPLINE_BLOCK_SIZE;    
    mygrad += SPLINE_BLOCK_SIZE;    
    mylapl += SPLINE_BLOCK_SIZE;    
    c += SPLINE_BLOCK_SIZE;
  }
      
  int remainder = 2*N - numBlocks*SPLINE_BLOCK_SIZE;
  if (tid < remainder) {
    *myval = (a[0] * c[0] +
	      a[1] * c[stride] +
	      a[2] * c[stride2] +
	      a[3] * c[stride3]);
  }
}

extern "C" void
eval_multi_multi_UBspline_1d_c_vgl_cuda (const multi_UBspline_1d_c_cuda *spline,
					 float *pos_d, float *vals_d[], 
					 float *grads_d[], float *lapl_d[], int num)
{
  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(num);

  eval_multi_multi_UBspline_1d_c_vgl_kernel<<<dimGrid,dimBlock>>>
    (pos_d, spline->gridInv, (float*)spline->coefs, vals_d, grads_d, lapl_d,
     spline->dim, spline->stride, spline->num_splines);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in eval_multi_multi_UBspline_1d_c_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}



__global__ static void
eval_multi_multi_UBspline_3d_c_kernel 
(float *pos, float3 drInv, const float *coefs, float *vals[], 
 uint3 dim, uint3 strides, int N)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+thr;

  __shared__ float *myval;
  __shared__ float abc[64];

  __shared__ float3 r;
  if (thr == 0) {
    r.x = pos[3*ir+0];
    r.y = pos[3*ir+1];
    r.z = pos[3*ir+2];
    myval = vals[ir];
  }
  __syncthreads();
  
  int3 index;
  float3 t;
  float s, sf;
  float4 tp[3];

  s = r.x * drInv.x;
  sf = floor(s);
  index.x = min(max(0,(int)sf), dim.x-1);
  //index.x = (int)sf;
  t.x = s - sf;

  s = r.y * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);
  //index.y = (int)sf;
  t.y = s - sf;

  s = r.z * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
  //index.z = (int)sf;
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
  
  if (thr < 64)
    abc[thr] = a[i]*b[j]*c[k];
  __syncthreads();

  if (off < 2*N) {
    float val = 0.0;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	const float *base = coefs + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
	for (int k=0; k<4; k++) 
	  val += abc[16*i+4*j+k] * base[off+k*strides.z];
      }
    }
    myval[off] = val;
  }
}



__global__ static void
eval_multi_multi_UBspline_3d_c_vgh_kernel 
(float *pos, float3 drInv,  const float *coefs, 
 float *vals[], float *grads[], float *hess[], 
 uint3 dim, uint3 strides, int N)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+threadIdx.x;

  __shared__ float *myval, *mygrad, *myhess;
  __shared__ float3 r;
  if (thr == 0) {
    r.x = pos[3*ir+0];
    r.y = pos[3*ir+1];
    r.z = pos[3*ir+2];
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
  index.x = min(max(0,(int)sf), dim.x-1);
  t.x = s - sf;

  s = r.y * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);

  t.y = s - sf;

  s = r.z * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
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
  const float *b0 = coefs + index.x*strides.x + index.y*strides.y + index.z*strides.z + off;
  if (off < 2*N) {
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	const float *base = b0 + i*strides.x + j*strides.y;
	float c0  = base[0*strides.z];
	float c1  = base[1*strides.z];
	float c2  = base[2*strides.z];
	float c3  = base[3*strides.z];
	v   += abc[n+  0]*c0 + abc[n+  1]*c1  + abc[n+  2]*c2  + abc[n+  3]*c3;
	g0  += abc[n+ 64]*c0 + abc[n+ 65]*c1  + abc[n+ 66]*c2  + abc[n+ 67]*c3;
	g1  += abc[n+128]*c0 + abc[n+129]*c1  + abc[n+130]*c2  + abc[n+131]*c3;
	g2  += abc[n+192]*c0 + abc[n+193]*c1  + abc[n+194]*c2  + abc[n+195]*c3;
	h00 += abc[n+256]*c0 + abc[n+257]*c1  + abc[n+258]*c2  + abc[n+259]*c3;
	h01 += abc[n+320]*c0 + abc[n+321]*c1  + abc[n+322]*c2  + abc[n+323]*c3;
	h02 += abc[n+384]*c0 + abc[n+385]*c1  + abc[n+386]*c2  + abc[n+387]*c3;
	h11 += abc[n+448]*c0 + abc[n+449]*c1  + abc[n+450]*c2  + abc[n+451]*c3;
	h12 += abc[n+512]*c0 + abc[n+513]*c1  + abc[n+514]*c2  + abc[n+515]*c3;
	h22 += abc[n+576]*c0 + abc[n+577]*c1  + abc[n+578]*c2  + abc[n+579]*c3;
	n += 4;
	// for (int k=0; k<4; k++) {
	//   float c  = base[k*strides.z];
	//   v   += abc[n+0] * c;
	//   g0  += abc[n+64] * c;
	//   g1  += abc[n+128] * c;
	//   g2  += abc[n+192] * c;
	//   h00 += abc[n+256] * c;
	//   h01 += abc[n+320] * c;
	//   h02 += abc[n+384] * c;
	//   h11 += abc[n+448] * c;
	//   h12 += abc[n+512] * c;
	//   h22 += abc[n+576] * c;
	//   n += 1;
	// }
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
  }
  abc[3*thr+0] = g0; 
  abc[3*thr+1] = g1; 
  abc[3*thr+2] = g2; 
  __syncthreads();
  for (int i=0; i<3; i++) {
    int myoff = (3*block+i)*SPLINE_BLOCK_SIZE + thr;
    if (myoff < 6*N)
      mygrad[myoff] = abc[i*SPLINE_BLOCK_SIZE+thr]; 
  }
  __syncthreads();

  // Write Hessians
  abc[6*thr+0]  = h00;
  abc[6*thr+1]  = h01;
  abc[6*thr+2]  = h02;
  abc[6*thr+3]  = h11;
  abc[6*thr+4]  = h12;
  abc[6*thr+5]  = h22;
  __syncthreads();
  for (int i=0; i<6; i++) {
    int myoff = (6*block+i)*SPLINE_BLOCK_SIZE + thr;
    if (myoff < 12*N)
      myhess[myoff] = abc[i*SPLINE_BLOCK_SIZE+thr];
  }
}


extern "C" void
eval_multi_multi_UBspline_3d_c_cuda (const multi_UBspline_3d_c_cuda *spline,
				     float *pos_d, complex_float *vals_d[], int num)
{
  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(2*spline->num_splines/SPLINE_BLOCK_SIZE, num);

  if (2*spline->num_splines % SPLINE_BLOCK_SIZE)
    dimGrid.x++;
  

  eval_multi_multi_UBspline_3d_c_kernel<<<dimGrid,dimBlock>>>
    (pos_d, spline->gridInv, (float*)spline->coefs, (float**)vals_d, spline->dim, spline->stride, spline->num_splines);


  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in eval_multi_multi_UBspline_3d_c_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

}

extern "C" void
eval_multi_multi_UBspline_3d_c_vgh_cuda (const multi_UBspline_3d_c_cuda *spline,
					 float *pos_d, complex_float *vals_d[], complex_float *grads_d[],
					 complex_float *hess_d[], int num)
{
  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(2*spline->num_splines/SPLINE_BLOCK_SIZE, num);

  if ((2*spline->num_splines) % SPLINE_BLOCK_SIZE)
    dimGrid.x++;

  eval_multi_multi_UBspline_3d_c_vgh_kernel<<<dimGrid,dimBlock>>>
    (pos_d, spline->gridInv, (float*)spline->coefs, 
     (float**)vals_d, (float**)grads_d, (float**)hess_d,
     spline->dim, spline->stride, spline->num_splines);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in eval_multi_multi_UBspline_3d_c_vgh_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}


__global__ static void
eval_multi_multi_UBspline_3d_c_vgl_kernel 
(float *pos, float3 drInv,  const float *coefs,  float Linv[],
 float *vals[], float *grad_lapl[], uint3 dim, uint3 strides,
 int N, int row_stride)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+threadIdx.x;

  __shared__ float *myval, *mygrad_lapl;
  __shared__ float3 r;
  if (thr == 0) {
    r.x = pos[3*ir+0];
    r.y = pos[3*ir+1];
    r.z = pos[3*ir+2];
    myval  = vals[ir];
    mygrad_lapl = grad_lapl[ir];
  }
  __syncthreads();
  
  int3 index;
  float3 t;
  float s, sf;
  float4 tp[3];

  s = r.x * drInv.x;
  sf = floor(s);
  index.x = min(max(0,(int)sf), dim.x-1);
  t.x = s - sf;

  s = r.y * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);
  t.y = s - sf;

  s = r.z * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
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
  const float *b0 = coefs + index.x*strides.x + index.y*strides.y + index.z*strides.z + off;
  if (off < 2*N) {
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	const float *base = b0 + i*strides.x + j*strides.y;
	float c0  = base[0*strides.z];
	float c1  = base[1*strides.z];
	float c2  = base[2*strides.z];
	float c3  = base[3*strides.z];
	v   += abc[n+  0]*c0 + abc[n+  1]*c1  + abc[n+  2]*c2  + abc[n+  3]*c3;
	g0  += abc[n+ 64]*c0 + abc[n+ 65]*c1  + abc[n+ 66]*c2  + abc[n+ 67]*c3;
	g1  += abc[n+128]*c0 + abc[n+129]*c1  + abc[n+130]*c2  + abc[n+131]*c3;
	g2  += abc[n+192]*c0 + abc[n+193]*c1  + abc[n+194]*c2  + abc[n+195]*c3;
	h00 += abc[n+256]*c0 + abc[n+257]*c1  + abc[n+258]*c2  + abc[n+259]*c3;
	h01 += abc[n+320]*c0 + abc[n+321]*c1  + abc[n+322]*c2  + abc[n+323]*c3;
	h02 += abc[n+384]*c0 + abc[n+385]*c1  + abc[n+386]*c2  + abc[n+387]*c3;
	h11 += abc[n+448]*c0 + abc[n+449]*c1  + abc[n+450]*c2  + abc[n+451]*c3;
	h12 += abc[n+512]*c0 + abc[n+513]*c1  + abc[n+514]*c2  + abc[n+515]*c3;
	h22 += abc[n+576]*c0 + abc[n+577]*c1  + abc[n+578]*c2  + abc[n+579]*c3;
	n += 4;
	// for (int k=0; k<4; k++) {
	//   float c  = base[k*strides.z];
	//   v   += abc[n+  0] * c;
	//   g0  += abc[n+ 64] * c;
	//   g1  += abc[n+128] * c;
	//   g2  += abc[n+192] * c;
	//   h00 += abc[n+256] * c;
	//   h01 += abc[n+320] * c;
	//   h02 += abc[n+384] * c;
	//   h11 += abc[n+448] * c;
	//   h12 += abc[n+512] * c;
	//   h22 += abc[n+576] * c;
	//   n += 1;
	// }
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
  }

  __shared__ float G[3][3], GGt[3][3];
  int i0 = threadIdx.x/3;
  int i1 = threadIdx.x - 3*i0;
  if (threadIdx.x < 9) 
    G[i0][i1] = Linv[threadIdx.x];
  __syncthreads();
  if (threadIdx.x < 9)   
    GGt[i0][i1] = (G[0][i0]*G[0][i1] + 
		   G[1][i0]*G[1][i1] + 
		   G[2][i0]*G[2][i1]);

  __syncthreads();
  if (off < 2*N) {
    // Store gradients back to global memory
    mygrad_lapl[off+0*row_stride] = G[0][0]*g0 + G[0][1]*g1 + G[0][2]*g2;
    mygrad_lapl[off+2*row_stride] = G[1][0]*g0 + G[1][1]*g1 + G[1][2]*g2;
    mygrad_lapl[off+4*row_stride] = G[2][0]*g0 + G[2][1]*g1 + G[2][2]*g2;
    
    // Store laplacians back to global memory
    // Hessian = H00 H01 H02 H11 H12 H22
    // Matrix = [0 1 2]
    //          [1 3 4]
    //          [2 4 5]
    // laplacian = Trace(GGt*Hessian)
    mygrad_lapl[off+6*row_stride] =
      (GGt[0][0]*h00 + GGt[1][0]*h01 + GGt[2][0]*h02 +
       GGt[0][1]*h01 + GGt[1][1]*h11 + GGt[2][1]*h12 +
       GGt[0][2]*h02 + GGt[1][2]*h12 + GGt[2][2]*h22);
  }
}


extern "C" void
eval_multi_multi_UBspline_3d_c_vgl_cuda 
(const multi_UBspline_3d_c_cuda *spline, float *pos_d, float *Linv_d, 
 float *vals_d[], float *grad_lapl_d[], int num, int row_stride)
{
  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(2*spline->num_splines/SPLINE_BLOCK_SIZE, num);

  if ((2*spline->num_splines) % SPLINE_BLOCK_SIZE)
    dimGrid.x++;

  eval_multi_multi_UBspline_3d_c_vgl_kernel<<<dimGrid,dimBlock>>>
    (pos_d, spline->gridInv, (float*)spline->coefs, Linv_d, (float**)vals_d, 
     (float**)grad_lapl_d, spline->dim, spline->stride, spline->num_splines, row_stride);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in eval_multi_multi_UBspline_3d_c_vgl_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}












/*


__global__ static void
eval_multi_multi_UBspline_3d_c_cuda (float *pos, float3 drInv, 
				     float *coefs_real, float *coefs_imag,
				     float *vals[], uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+thr;

  __shared__ float *myval;
  __shared__ float abc[64];

  // __shared__ float pos_s[SPLINE_BLOCK_SIZE];
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
  __shared__ float buff[2*SPLINE_BLOCK_SIZE];
  buff[2*thr+0] = val_real;
  buff[2*thr+1] = val_imag;
  __syncthreads();
  myval[off] = buff[thr];
  myval[off+SPLINE_BLOCK_SIZE] = buff[thr+SPLINE_BLOCK_SIZE];
}



__global__ static void
eval_multi_multi_UBspline_3d_c_vgh_cuda (float *pos, float3 drInv, 
					 float *coefs_real, float *coefs_imag,
					 float *vals[], float *grads[], 
					 float *hess[], uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+thr;

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

  
  __shared__ float buff[6*SPLINE_BLOCK_SIZE];
  // Note, we can reuse abc, by replacing buff with abc.
  
  buff[2*thr+0] = v_r;  buff[2*thr+1] = v_i;
  __syncthreads();
  myval[off] = buff[thr];    
  myval[off+SPLINE_BLOCK_SIZE] = buff[thr+SPLINE_BLOCK_SIZE];

  buff[6*thr+0] = g0_r;  buff[6*thr+1] = g0_i;
  buff[6*thr+2] = g1_r;  buff[6*thr+3] = g1_i;
  buff[6*thr+4] = g2_r;  buff[6*thr+5] = g2_i;
  __syncthreads();
  for (int i=0; i<6; i++) 
    mygrad[(6*block+i)*SPLINE_BLOCK_SIZE+thr] = buff[i*SPLINE_BLOCK_SIZE+thr]; 
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
      myhess[(12*block+i)*SPLINE_BLOCK_SIZE+thr] = buff[i*SPLINE_BLOCK_SIZE+thr];

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
      myhess[(12*block+i+6)*SPLINE_BLOCK_SIZE+th2] = buff[i*SPLINE_BLOCK_SIZE+th2];
  }
}
*/
#endif
