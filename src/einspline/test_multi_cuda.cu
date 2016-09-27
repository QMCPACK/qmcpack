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


#include "multi_bspline.h"
#include "multi_bspline_create_cuda.h"
#include "multi_bspline_structs_cuda.h"
#include "multi_bspline_eval_cuda.h"


void
test_float_1d()
{
  int numWalkers = 1000;
  float *vals[numWalkers], *grads[numWalkers], *hess[numWalkers];
  float *coefs,  **vals_d, **grads_d, **hess_d;
  float *r_d, *r_h;
  int xs, N;
  int Nx;

  N = 128*36;
  Nx = 100;
  xs = N;
  // Setup Bspline coefficients
  int size = Nx*N*sizeof(float);
  posix_memalign((void**)&coefs, 16, size);
  for (int ix=0; ix<Nx; ix++)
    for (int n=0; n<N; n++)
      coefs[ix*xs+ n] = drand48();

  Ugrid x_grid;
  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nx;
  BCtype_s xBC;
  xBC.lCode = xBC.rCode = PERIODIC;

  multi_UBspline_1d_s *spline = 
    create_multi_UBspline_1d_s (x_grid, xBC, N);
  for (int i=0; i<N; i++) 
    set_multi_UBspline_1d_s (spline, i, coefs);

  multi_UBspline_1d_s_cuda *cudaspline = 
    create_multi_UBspline_1d_s_cuda (spline);

  // Setup device value storage
  int numVals = N*numWalkers*3;
  float *valBlock_d, *valBlock_h;
  cudaMalloc((void**)&(valBlock_d),     numVals*sizeof(float));
  cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(float));
  cudaMalloc((void**)&(vals_d),  numWalkers*sizeof(float*));
  cudaMalloc((void**)&(grads_d), numWalkers*sizeof(float*));
  cudaMalloc((void**)&(hess_d),  numWalkers*sizeof(float*));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals[i]  = valBlock_d + i*N;
    grads[i] = valBlock_d + N*numWalkers   + i*N;
    hess[i]  = valBlock_d + 2*N*numWalkers + i*N;
  }
  cudaMemcpy(vals_d,  vals,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(grads_d, grads, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(hess_d,  hess,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  fprintf (stderr, "Finished cuda allocations.\n");

  // Setup walker positions
  cudaMalloc((void**)&(r_d),     numWalkers*sizeof(float));
  cudaMallocHost((void**)&(r_h), numWalkers*sizeof(float));
  fprintf (stderr, "r_h = %p\n", r_h);

  for (int ir=0; ir<numWalkers; ir++) 
    r_h[ir] = 0.5*drand48();

  float vals_host[N], vals_cuda[N];

  // Check value
  for (int w=0; w<numWalkers; w++) {
    eval_multi_UBspline_1d_s (spline, r_h[w], vals_host);
    cudaMemcpy(r_d, r_h, numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_1d_s_cuda (cudaspline, r_d, vals_d, numWalkers);
    cudaMemcpy(vals_cuda, valBlock_d+(N*w), N*sizeof(float), cudaMemcpyDeviceToHost);
    //for (int i=0; i<N; i++)
    if (w < 10)
      fprintf (stderr, "%3i  %15.8e %15.8e\n", w, vals_host[0], vals_cuda[0]);
  }


  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_1d_s_cuda (cudaspline, r_d, vals_d, numWalkers);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);

  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_1d_s_vgl_cuda (cudaspline, r_d, vals_d, grads_d, hess_d, numWalkers);
  }
  end = clock();
  time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "VGL Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (cudaspline->coefs);
  cudaFree (valBlock_d);
  cudaFree (vals_d);
  cudaFree (grads_d);
  cudaFree (hess_d);
  cudaFree (r_d);
}



void
test_float()
{
  int numWalkers = 1024;
  float *vals[numWalkers], *grads[numWalkers], *hess[numWalkers];
  float *coefs,  **vals_d, **grads_d, **hess_d;
  float *r_d, *r_h;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 256;
  Nx = Ny = Nz = 32;
  xs = Ny*Nz*N;
  ys = Nz*N;
  zs = N;

  // Setup Bspline coefficients
  int size = Nx*Ny*Nz*N*sizeof(float);
  posix_memalign((void**)&coefs, 16, size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int n=0; n<N; n++)
	  coefs[ix*xs + iy*ys + iz*zs + n] = drand48();

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nx;
  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num = Ny;
  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num = Nz;
  BCtype_s xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;
  

  multi_UBspline_3d_s *spline = 
    create_multi_UBspline_3d_s (x_grid, y_grid, z_grid, xBC, yBC, zBC, N);
  for (int i=0; i<N; i++) 
    set_multi_UBspline_3d_s (spline, i, coefs);

  multi_UBspline_3d_s_cuda *cudaspline = 
    create_multi_UBspline_3d_s_cuda (spline);

  // Setup device value storage
  int numVals = N*numWalkers*10;
  float *valBlock_d, *valBlock_h;
  cudaMalloc((void**)&(valBlock_d),     numVals*sizeof(float));
  cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(float));
  cudaMalloc((void**)&(vals_d),  numWalkers*sizeof(float*));
  cudaMalloc((void**)&(grads_d), numWalkers*sizeof(float*));
  cudaMalloc((void**)&(hess_d),  numWalkers*sizeof(float*));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals[i]  = valBlock_d + i*N;
    grads[i] = valBlock_d + N*numWalkers + 3*i*N;
    hess[i]  = valBlock_d + 4*N*numWalkers + 6*i*N;
  }
  cudaMemcpy(vals_d,  vals,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(grads_d, grads, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(hess_d,  hess,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  fprintf (stderr, "Finished cuda allocations.\n");

  // Setup walker positions
  cudaMalloc((void**)&(r_d),     3*numWalkers*sizeof(float));
  cudaMallocHost((void**)&(r_h), 3*numWalkers*sizeof(float));

  for (int ir=0; ir<numWalkers; ir++) {
    r_h[3*ir+0] = 0.5*drand48();
    r_h[3*ir+1] = 0.5*drand48();
    r_h[3*ir+2] = 0.5*drand48();
  }

  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(N/SPLINE_BLOCK_SIZE,numWalkers);
  
  float vals_host[N], vals_cuda[N];

  // Check value
  for (int w=0; w<numWalkers; w++) {
    eval_multi_UBspline_3d_s (spline, r_h[3*w+0], r_h[3*w+1], r_h[3*w+2], vals_host);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_s_cuda (cudaspline, r_d, vals_d, numWalkers);
    cudaMemcpy(vals_cuda, valBlock_d+(N*w), N*sizeof(float), cudaMemcpyDeviceToHost);
    //for (int i=0; i<N; i++)
    if (w < 10)
      fprintf (stderr, "%3i  %15.8e %15.8e\n", w, vals_host[0], vals_cuda[0]);
  }


  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_s_cuda (cudaspline, r_d, vals_d, numWalkers);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);

  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_s_vgh_cuda (cudaspline, r_d, vals_d, grads_d, hess_d, numWalkers);
  }
  end = clock();
  time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "VGH Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (cudaspline->coefs);
  cudaFree (valBlock_d);
  cudaFree (vals_d);
  cudaFree (grads_d);
  cudaFree (hess_d);
  cudaFree (r_d);
}



void
test_complex_float()
{
  int numWalkers = 1000;
  complex_float *vals[numWalkers], *grads[numWalkers], *hess[numWalkers];
  complex_float *coefs,  **vals_d, **grads_d, **hess_d;
  float *Linv_d;
  float *r_d, *r_h;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 128;
  Nx = Ny = Nz = 32;
  xs = Ny*Nz*N;
  ys = Nz*N;
  zs = N;

  // Setup Bspline coefficients
  int size = Nx*Ny*Nz*N*sizeof(complex_float);
  posix_memalign((void**)&coefs, 16, size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int n=0; n<N; n++)
	  coefs[ix*xs + iy*ys + iz*zs + n] = std::complex<float>(drand48(), drand48());

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nx;
  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num = Ny;
  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num = Nz;
  BCtype_c xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;
  

  multi_UBspline_3d_c *spline = 
    create_multi_UBspline_3d_c (x_grid, y_grid, z_grid, xBC, yBC, zBC, N);
  for (int i=0; i<N; i++) 
    set_multi_UBspline_3d_c (spline, i, coefs);

  multi_UBspline_3d_c_cuda *cudaspline = 
    create_multi_UBspline_3d_c_cuda (spline);

  // Setup device value storage
  int numVals = N*numWalkers*10;
  complex_float *valBlock_d, *valBlock_h;
  cudaMalloc((void**)&(valBlock_d),     numVals*sizeof(complex_float));
  cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(complex_float));
  cudaMalloc((void**)&(vals_d),  numWalkers*sizeof(complex_float*));
  cudaMalloc((void**)&(grads_d), numWalkers*sizeof(complex_float*));
  cudaMalloc((void**)&(hess_d),  numWalkers*sizeof(complex_float*));
  cudaMalloc((void**)&(Linv_d), 9*sizeof(float));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals[i]  = valBlock_d + i*N;
    grads[i] = valBlock_d + N*numWalkers + 3*i*N;
    hess[i]  = valBlock_d + 4*N*numWalkers + 6*i*N;
  }
  float Linv[9] = { 1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0 };
  cudaMemcpy(vals_d,  vals,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(grads_d, grads, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(hess_d,  hess,  numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(Linv_d,  Linv,  9*sizeof(float), cudaMemcpyHostToDevice);
  fprintf (stderr, "Finished cuda allocations.\n");

  // Setup walker positions
  cudaMalloc((void**)&(r_d),     3*numWalkers*sizeof(float));
  cudaMallocHost((void**)&(r_h), 3*numWalkers*sizeof(float));

  for (int ir=0; ir<numWalkers; ir++) {
    r_h[3*ir+0] = 0.5*drand48();
    r_h[3*ir+1] = 0.5*drand48();
    r_h[3*ir+2] = 0.5*drand48();
  }

  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(N/SPLINE_BLOCK_SIZE,numWalkers);
  
  complex_float vals_host[N], vals_cuda[N];

  

  // Check value
  for (int w=0; w<numWalkers; w++) {
    eval_multi_UBspline_3d_c (spline, r_h[3*w+0], r_h[3*w+1], r_h[3*w+2], vals_host);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    //eval_multi_multi_UBspline_3d_c_cuda (cudaspline, r_d, vals_d, numWalkers);
    //eval_multi_multi_UBspline_3d_c_vgh_cuda (cudaspline, r_d, vals_d, grads_d, hess_d, numWalkers);
    eval_multi_multi_UBspline_3d_c_vgl_cuda (cudaspline, r_d, Linv_d, vals_d, grads_d, numWalkers, N);
    cudaMemcpy(vals_cuda, valBlock_d+(N*w), N*sizeof(float), cudaMemcpyDeviceToHost);
    //for (int i=0; i<N; i++)
    if (w < 10)
      fprintf (stderr, "%3i  %15.8e %15.8e  %15.8e %15.8e\n", w, 
	       vals_host[0].real(), vals_cuda[0].real(),
	       vals_host[0].imag(), vals_cuda[0].imag());
  }


  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_c_cuda (cudaspline, r_d, vals_d, numWalkers);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);

  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(float), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_c_vgh_cuda (cudaspline, r_d, vals_d, grads_d, hess_d, numWalkers);
  }
  end = clock();
  time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "VGH Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (cudaspline->coefs);
  cudaFree (valBlock_d);
  cudaFree (vals_d);
  cudaFree (grads_d);
  cudaFree (hess_d);
  cudaFree (r_d);
}



void
test_double()
{
  int numWalkers = 1000;
  double *vals[numWalkers], *grads[numWalkers], *hess[numWalkers];
  double *coefs,  **vals_d, **grads_d, **hess_d;
  double *r_d, *r_h;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 128;
  Nx = Ny = Nz = 32;
  xs = Ny*Nz*N;
  ys = Nz*N;
  zs = N;

  // Setup Bspline coefficients
  int size = Nx*Ny*Nz*N*sizeof(double);
  posix_memalign((void**)&coefs, 16, size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int n=0; n<N; n++)
	  coefs[ix*xs + iy*ys + iz*zs + n] = drand48();

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nx;
  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num = Ny;
  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num = Nz;
  BCtype_d xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;
  

  multi_UBspline_3d_d *spline = 
    create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, N);
  for (int i=0; i<N; i++) 
    set_multi_UBspline_3d_d (spline, i, coefs);

  multi_UBspline_3d_d_cuda *cudaspline = 
    create_multi_UBspline_3d_d_cuda (spline);

  // Setup device value storage
  int numVals = N*numWalkers*10;
  double *valBlock_d, *valBlock_h;
  cudaMalloc((void**)&(valBlock_d),     numVals*sizeof(double));
  cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(double));
  cudaMalloc((void**)&(vals_d),  numWalkers*sizeof(double*));
  cudaMalloc((void**)&(grads_d), numWalkers*sizeof(double*));
  cudaMalloc((void**)&(hess_d),  numWalkers*sizeof(double*));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals[i]  = valBlock_d + i*N;
    grads[i] = valBlock_d + N*numWalkers + 3*i*N;
    hess[i]  = valBlock_d + 4*N*numWalkers + 6*i*N;
  }
  cudaMemcpy(vals_d,  vals,  numWalkers*sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(grads_d, grads, numWalkers*sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(hess_d,  hess,  numWalkers*sizeof(double*), cudaMemcpyHostToDevice);
  fprintf (stderr, "Finished cuda allocations.\n");

  // Setup walker positions
  cudaMalloc((void**)&(r_d),     3*numWalkers*sizeof(double));
  cudaMallocHost((void**)&(r_h), 3*numWalkers*sizeof(double));

  for (int ir=0; ir<numWalkers; ir++) {
    r_h[3*ir+0] = 0.5*drand48();
    r_h[3*ir+1] = 0.5*drand48();
    r_h[3*ir+2] = 0.5*drand48();
  }

  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(N/SPLINE_BLOCK_SIZE,numWalkers);
  
  double vals_host[N], vals_cuda[N];

  // Check value
  for (int w=0; w<numWalkers; w++) {
    eval_multi_UBspline_3d_d (spline, r_h[3*w+0], r_h[3*w+1], r_h[3*w+2], vals_host);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(double), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_d_cuda (cudaspline, r_d, vals_d, numWalkers);
    cudaMemcpy(vals_cuda, valBlock_d+(N*w), N*sizeof(double), cudaMemcpyDeviceToHost);
    //for (int i=0; i<N; i++)
    if (w < 10)
      fprintf (stderr, "%3i  %15.8e %15.8e\n", w, vals_host[0], vals_cuda[0]);
  }


  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(double), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_d_cuda (cudaspline, r_d, vals_d, numWalkers);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);

  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(double), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_d_vgh_cuda (cudaspline, r_d, vals_d, grads_d, hess_d, numWalkers);
  }
  end = clock();
  time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "VGH Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (cudaspline->coefs);
  cudaFree (valBlock_d);
  cudaFree (vals_d);
  cudaFree (grads_d);
  cudaFree (hess_d);
  cudaFree (r_d);
}



void
test_complex_double()
{
  int numWalkers = 1000;
  complex_double *vals[numWalkers], *grads[numWalkers], *hess[numWalkers];
  complex_double *coefs, **vals_d, **grads_d, **hess_d;
  double *r_d, *r_h;
  int xs, ys, zs, N;
  int Nx, Ny, Nz;

  N = 128;
  Nx = Ny = Nz = 32;
  xs = Ny*Nz*N;
  ys = Nz*N;
  zs = N;

  // Setup Bspline coefficients
  int size = Nx*Ny*Nz*N*sizeof(complex_double);
  posix_memalign((void**)&coefs, 16, size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int n=0; n<N; n++)
	  coefs[ix*xs + iy*ys + iz*zs + n] = std::complex<double>(drand48(), drand48());

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nx;
  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num = Ny;
  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num = Nz;
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;
  

  multi_UBspline_3d_z *spline = 
    create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, N);
  for (int i=0; i<N; i++) 
    set_multi_UBspline_3d_z (spline, i, coefs);

  multi_UBspline_3d_z_cuda *cudaspline = 
    create_multi_UBspline_3d_z_cuda (spline);

  // Setup device value storage
  int numVals = N*numWalkers*10;
  complex_double *valBlock_d, *valBlock_h;
  cudaMalloc((void**)&(valBlock_d),     numVals*sizeof(complex_double));
  cudaMallocHost((void**)&(valBlock_h), numVals*sizeof(complex_double));
  cudaMalloc((void**)&(vals_d),  numWalkers*sizeof(complex_double*));
  cudaMalloc((void**)&(grads_d), numWalkers*sizeof(complex_double*));
  cudaMalloc((void**)&(hess_d),  numWalkers*sizeof(complex_double*));
  fprintf (stderr, "valBlock_d = %p\n", valBlock_d);
  for (int i=0; i<numWalkers; i++) {
    vals[i]  = valBlock_d + i*N;
    grads[i] = valBlock_d + N*numWalkers + 3*i*N;
    hess[i]  = valBlock_d + 4*N*numWalkers + 6*i*N;
  }
  cudaMemcpy(vals_d,  vals,  numWalkers*sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(grads_d, grads, numWalkers*sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(hess_d,  hess,  numWalkers*sizeof(double*), cudaMemcpyHostToDevice);
  fprintf (stderr, "Finished cuda allocations.\n");

  // Setup walker positions
  cudaMalloc((void**)&(r_d),     3*numWalkers*sizeof(double));
  cudaMallocHost((void**)&(r_h), 3*numWalkers*sizeof(double));

  for (int ir=0; ir<numWalkers; ir++) {
    r_h[3*ir+0] = 0.5*drand48();
    r_h[3*ir+1] = 0.5*drand48();
    r_h[3*ir+2] = 0.5*drand48();
  }

  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid(N/SPLINE_BLOCK_SIZE,numWalkers);
  
  complex_double vals_host[N], vals_cuda[N];

  // Check value
  for (int w=0; w<numWalkers; w++) {
    eval_multi_UBspline_3d_z (spline, r_h[3*w+0], r_h[3*w+1], r_h[3*w+2], vals_host);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(double), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_z_cuda (cudaspline, r_d, vals_d, numWalkers);
    cudaMemcpy(vals_cuda, valBlock_d+(N*w), N*sizeof(double), cudaMemcpyDeviceToHost);
    //for (int i=0; i<N; i++)
    if (w < 10)
      fprintf (stderr, "%3i  %15.8e %15.8e  %15.8e %15.8e\n", w, 
	       vals_host[0].real(), vals_cuda[0].real(),
	       vals_host[0].imag(), vals_cuda[0].imag());
  }


  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(double), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_z_cuda (cudaspline, r_d, vals_d, numWalkers);
  }
  end = clock();
  double time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "Evals per second = %1.8e\n", 1.0/time);

  start = clock();
  for (int i=0; i<10000; i++) {
    if ((i%1000) == 0) 
      fprintf (stderr, "i = %d\n", i);
    cudaMemcpy(r_d, r_h, 3*numWalkers*sizeof(double), cudaMemcpyHostToDevice);
    eval_multi_multi_UBspline_3d_z_vgh_cuda (cudaspline, r_d, vals_d, grads_d, hess_d, numWalkers);
  }
  end = clock();
  time = (double)(end-start)/(double)((double)CLOCKS_PER_SEC*(double)10000*N*numWalkers);
  fprintf (stderr, "VGH Evals per second = %1.8e\n", 1.0/time);
  
  cudaFree (cudaspline->coefs);
  cudaFree (valBlock_d);
  cudaFree (vals_d);
  cudaFree (grads_d);
  cudaFree (hess_d);
  cudaFree (r_d);
}



main() 
{
//   int deviceCount;
//   cudaGetDeviceCount(&deviceCount);
//   int num_appropriate=0;
//   for (int device=0; device < deviceCount; ++device) {
//     cudaDeviceProp deviceProp;
//     cudaGetDeviceProperties(&deviceProp, device);
//     fprintf (stderr, "Device %d has architecture %d.%d\n",
// 	     device, deviceProp.major, deviceProp.minor);
//   }
//   cudaSetDevice(0);	
  // fprintf(stderr, "Testing 1D single-precision real routines:\n");
  // test_float_1d();
  fprintf(stderr, "Testing 3D single-precision real routines:\n");
  test_float();
  // fprintf(stderr, "Testing 3D single-precision complex routines:\n");
  // test_complex_float();
  // fprintf(stderr, "Testing 3D double-precision real routines:\n");
  // test_double();
  // fprintf(stderr, "Testing 3D double-precision complex routines:\n");
  // test_complex_double();
}
