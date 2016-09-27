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

#include <stdio.h>
#include <pthread.h>
#include <cuda.h>
#include <cutil.h>
#include <multithreading.h>

__global__ void 
eval_multi_UBspline_3d_cuda_c (const float *coefs, float *abc, float *vals,
			       int ix, int iy, int iz,
			       int xs, int ys, int zs, int N)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int offset = block*BLOCK_SIZE+thr;
  __shared__ float abcs[64];
  abcs[thr] = abc[thr];
  
  __syncthreads();

  float val= 0.0;
  //int index=0;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      for (int k=0; k<4; k++) {	
	float *base_addr = coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs;
	//val += abc[(16*i+4*j+k)*BLOCK_SIZE + thr] * base_addr[offset];
	val += abcs[16*i+4*j+k] * base_addr[offset];	
	//index++;
      }
    }
  vals[offset] = val;
}


__constant__ float A[16], dA[16], d2A[16];

__global__ static void
eval_multi_multi_UBspline_3d_cuda_c (float *pos, float3 drInv, 
				     const float *coefs_real,const  float *coefs_imag,
				     float *vals_real, float *vals_imag, 
				     int3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int offset = block*BLOCK_SIZE+thr;

  __shared__ float abc[64];

  __shared__ float pos_s[BLOCK_SIZE];
  int ir1 = (ir >> 4)*64;
  int ir2 = (ir & 15)*4;
  pos_s[thr] = pos[ir1+thr];
  __syncthreads();
  float3 r;
  r.x = pos_s[ir2+0];
  r.y = pos_s[ir2+1];
  r.z = pos_s[ir2+2];
  
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
  
  tp[0] = make_float4(1.0, t.x, t.x*t.x, t.x*t.x*t.x);
  tp[1] = make_float4(1.0, t.y, t.y*t.y, t.y*t.y*t.y);
  tp[2] = make_float4(1.0, t.z, t.z*t.z, t.z*t.z*t.z);

  __shared__ float a[4], b[4], c[4];
  if (thr == 0) {
    a[0] = A[ 0]*tp[0].x + A[ 1]*tp[0].y + A[ 2]*tp[0].z + A[ 3]*tp[0].w;
    a[1] = A[ 4]*tp[0].x + A[ 5]*tp[0].y + A[ 6]*tp[0].z + A[ 7]*tp[0].w;
    a[2] = A[ 8]*tp[0].x + A[ 9]*tp[0].y + A[10]*tp[0].z + A[11]*tp[0].w;
    a[3] = A[12]*tp[0].x + A[13]*tp[0].y + A[14]*tp[0].z + A[15]*tp[0].w;
    
    b[0] = A[ 0]*tp[1].x + A[ 1]*tp[1].y + A[ 2]*tp[1].z + A[ 3]*tp[1].w;
    b[1] = A[ 4]*tp[1].x + A[ 5]*tp[1].y + A[ 6]*tp[1].z + A[ 7]*tp[1].w;
    b[2] = A[ 8]*tp[1].x + A[ 9]*tp[1].y + A[10]*tp[1].z + A[11]*tp[1].w;
    b[3] = A[12]*tp[1].x + A[13]*tp[1].y + A[14]*tp[1].z + A[15]*tp[1].w;
    
    c[0] = A[ 0]*tp[2].x + A[ 1]*tp[2].y + A[ 2]*tp[2].z + A[ 3]*tp[2].w;
    c[1] = A[ 4]*tp[2].x + A[ 5]*tp[2].y + A[ 6]*tp[2].z + A[ 7]*tp[2].w;
    c[2] = A[ 8]*tp[2].x + A[ 9]*tp[2].y + A[10]*tp[2].z + A[11]*tp[2].w;
    c[3] = A[12]*tp[2].x + A[13]*tp[2].y + A[14]*tp[2].z + A[15]*tp[2].w;
  }

  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);
  
  abc[thr] = a[i]*b[j]*c[k];
  __syncthreads();

  float val_real = 0.0;
  float val_imag = 0.0;
  //int index=0;
  val_real = val_imag = 0.0;
//   int di = strides.x - 4*strides.y;
//   int dj = strides.y - 4*strides.z;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      for (int k=0; k<4; k++) {
	// 	float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + (index.z+k)*strides.z;
	// 	float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + (index.z+k)*strides.z;
	val_real += abc[16*i+4*j+k] * base_real[offset+k*strides.z];
	val_imag += abc[16*i+4*j+k] * base_imag[offset+k*strides.z];
 	// base_real += strides.z;
 	// base_imag += strides.z;
      }
//       base_real += dj;
//       base_imag += dj;
    }
//     base_real += di;
//     base_imag += di;
  }
  vals_real[offset+ir*128] = val_real;
  vals_imag[offset+ir*128] = val_imag;
  //vals_real[ir][offset] = val_real;
  // vals_imag[ir][offset] = val_imag;
}
				    


// __global__ void 
// eval_multi_UBspline_3d_cuda_c2 (float3 r,
// 				float *coefs, float *vals,
// 				int xs, int ys, int zs, int N)
// {
//   int block = blockIdx.x;
//   int thr   = threadIdx.x;

//   __shared__ float abcs[64];
//   abcs[thr] = abc[thr];

//   float dxInv = 0.0625f;
//   float v, dv;

//   v = floor(dxInv*r.x);
//   dv = dxInv*r.x - v;
//   int ix = (int) v;

//   v = floor(dxInv*r.x);
//   dv = dxInv*r.x - v;
//   int iy = (int) v;

//   v = floor(dxInv*r.y);
//   dv = dxInv*r.y - v;
//   int iz = (int) v;

//   int offset = block*BLOCK_SIZE+thr;
//   __shared__ float abcs[64];
//   abcs[thr] = abc[thr];
  

//   float val= 0.0;
//   //int index=0;
//   val = 0.0;
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++)
//       for (int k=0; k<4; k++) {
// 	float *base_addr = coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs;
// 	//val += abc[(16*i+4*j+k)*BLOCK_SIZE + thr] * base_addr[offset];
// 	val += abcs[16*i+4*j+k] * base_addr[offset];	
// 	//index++;
//       }
//   vals[offset] = val;
// }


void
test_cuda()
{
  float *coefs  , *abc  , *abc2, *vals;
  float *coefs_d, *abc_d, *vals_d;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 4096;
  Nx = Ny = Nz = 16;
  xs = Nx*Ny*Nz;
  ys = Ny*Nz;
  zs = Nz;
  
  int size = Nx*Ny*Nz*N*sizeof(float);
  posix_memalign((void**)&coefs, 16, size);
  cudaMalloc((void**)&coefs_d, size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int n=0; n<N; n++)
	  coefs[ix*xs + iy*ys + iz*zs + n] = drand48();
  cudaMemcpy(coefs_d, coefs, size, cudaMemcpyHostToDevice);

  posix_memalign ((void**)&abc, 16, 64*sizeof(float));
  posix_memalign ((void**)&abc2, 16, 64*BLOCK_SIZE*sizeof(float));
  cudaMalloc((void**)&abc_d, 64*BLOCK_SIZE*sizeof(float));
  for (int i=0; i<64; i++) {
    abc[i] = drand48();
    for (int j=0; j<BLOCK_SIZE; j++)
      abc2[i*BLOCK_SIZE+j] = abc[i];
  }
  //  cudaMemcpy(abc_d, abc2, 64*BLOCK_SIZE*sizeof(float), 
  //     cudaMemcpyHostToDevice);
  cudaMemcpy(abc_d, abc, 64*sizeof(float), 
	     cudaMemcpyHostToDevice);

  posix_memalign((void**)&vals, 16, N*sizeof(float));
  cudaMalloc((void**)&vals_d, N*sizeof(float));

  dim3 dimBlock(BLOCK_SIZE);
  dim3 dimGrid(N/BLOCK_SIZE);

  int ix=1; 
  int iy=2;
  int iz=3;
  
  clock_t start, end;
  start = clock();
  for (int i=0; i<100000; i++) {
    eval_multi_UBspline_3d_cuda_c<<<dimGrid,dimBlock>>> 
      (coefs_d, abc_d, vals_d, ix, iy, iz, xs, ys, zs, N);
  }
  end = clock();
  double time = (double)(end-start)/(double)(CLOCKS_PER_SEC*100000*N);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);

  cudaMemcpy (vals, vals_d, N*sizeof(float), cudaMemcpyDeviceToHost);

  float vals2[N];
  
  for (int n=0; n<N; n++) {
    vals2[n] = 0.0;
    int index=0;
    for(int i=0; i<4; i++)
      for (int j=0; j<4; j++)
	for (int k=0; k<4; k++)  {
	  vals2[n] += abc[index] * coefs[(ix+i)*xs+(iy+j)*ys+(iz+k)*zs+n];
	  index++;
	}
  }


  for (int i=0; i<N/256; i++)	
    fprintf (stderr, "%1.9f %1.9f\n", vals[i], vals2[i]); 


  cudaFree(abc_d);
  cudaFree(coefs_d);
  cudaFree(vals_d);
}


static void *
test_multi_cuda(void *thread)
{
//   CUcontext ctx;
//   CUdevice dev;
//   cuDeviceGet (&dev, (int)(size_t)thread);
//   cuCtxCreate(&ctx, CU_CTX_SCHED_YIELD, dev);

//   int deviceCount;
//   cudaGetDeviceCount(&deviceCount);

  CUDA_SAFE_CALL(cudaSetDevice((int)(size_t)thread));
  fprintf (stderr, "In thread %p\n", thread);


  int numWalkers = 2000;
  float *coefs  ,  __device__ *vals_real[numWalkers],   __device__ *vals_imag[numWalkers];
  float *coefs_real_d, *coefs_imag_d, __device__ *vals_real_d[numWalkers], __device__ *vals_imag_d[numWalkers];
  float *r_d, *r_h;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 128;
  Nx = Ny = Nz = 64;
  xs = Ny*Nz*N;
  ys = Nz*N;
  zs = N;

  float3 drInv;
  drInv.x = 1.0/float(Nx);
  drInv.y = 1.0/float(Ny);
  drInv.z = 1.0/float(Nz);

  // Setup Bspline coefficients
  int size = Nx*Ny*Nz*N*sizeof(float);
  CUT_SAFE_MALLOC(posix_memalign((void**)&coefs, 16, size));
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
  CUDA_SAFE_CALL(cudaMalloc((void**)&coefs_real_d, size));
  CUDA_SAFE_CALL(cudaMalloc((void**)&coefs_imag_d, size));
  fprintf (stderr, "Before Memcpy.\n");
  CUDA_SAFE_CALL(cudaMemcpy(coefs_real_d, coefs, size, cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(coefs_imag_d, coefs, size, cudaMemcpyHostToDevice));
  fprintf (stderr, "After Memcpy.\n");  

  // Setup device value storage
  int numVals = 2*N*numWalkers;
  float *valBlock_d, *valBlock_h;
  CUDA_SAFE_CALL(cudaMalloc((void**)&(valBlock_d), numVals*sizeof(float)));
  CUDA_SAFE_CALL(cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(float)));
  CUDA_SAFE_CALL(cudaMalloc((void**)&(vals_real_d), numWalkers*sizeof(float*)));
  CUDA_SAFE_CALL(cudaMalloc((void**)&(vals_imag_d), numWalkers*sizeof(float*)));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals_real[i] = valBlock_d + 2*i*N;
    vals_imag[i] = valBlock_d + (2*i+1)*N;
  }
  CUDA_SAFE_CALL(cudaMemcpy(vals_real_d, vals_real, numWalkers*sizeof(float*), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(vals_imag_d, vals_imag, numWalkers*sizeof(float*), cudaMemcpyHostToDevice));
  
  fprintf (stderr, "Finished cuda allocations.\n");


  // Setup walker positions
  CUDA_SAFE_CALL(cudaMalloc((void**)&(r_d),     4*numWalkers*sizeof(float)));
  CUDA_SAFE_CALL(cudaMallocHost((void**)&(r_h), 4*numWalkers*sizeof(float)));

  for (int ir=0; ir<numWalkers; ir++) {
    r_h[4*ir+0] = 0.75*drand48();
    r_h[4*ir+1] = 0.75*drand48();
    r_h[4*ir+2] = 0.75*drand48();
  }

  
  int3 strides;
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
    CUDA_SAFE_CALL(cudaMemcpy(r_d, r_h, 4*numWalkers*sizeof(float), cudaMemcpyHostToDevice));
    // eval_multi_multi_UBspline_3d_cuda_c<<<dimGrid,dimBlock>>> 
    //   (r_d, drInv, coefs_real_d, coefs_imag_d, 
    //    vals_real_d, vals_imag_d, strides);
    eval_multi_multi_UBspline_3d_cuda_c<<<dimGrid,dimBlock>>> 
      (r_d, drInv, coefs_real_d, coefs_imag_d, 
       valBlock_d, valBlock_d+numVals/2, strides);
    //cudaMemcpy(valBlock_h, valBlock_d, numVals*sizeof(float), cudaMemcpyDeviceToHost);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (valBlock_d);
  cudaFree (vals_real_d);
  cudaFree (vals_imag_d);
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
  }

  //  pthread_t threads[deviceCount];

  // for (int device = 0; device < deviceCount; device++) 
  //   pthread_create (&(threads[device]), NULL, test_multi_cuda, (void*)device);
  // cutStartThread((CUT_THREADROUTINE)test_multi_cuda,(void*)device);
  test_multi_cuda((void*)0);

  //  pthread_exit(NULL);
  //test_multi_cuda();
}
