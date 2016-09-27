//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <cstdio>
#include <vector>
#include <complex>
#include "AtomicOrbitalCuda.h"


__constant__ float  Acuda[48];

const int MaxQuad = 12;

bool atomic_cuda_initialized = false;

// type traits for cuda variable make types
template<typename T> __device__ typename cudaTypeTraits<T>::realType4 cudaMakeType4(T a, T b, T c, T d); 

template<> __device__ typename cudaTypeTraits<float>::realType4 cudaMakeType4(float a, float b, float c, float d)
{
  return make_float4(a, b, c, d); 
};

template<> __device__ typename cudaTypeTraits<double>::realType4 cudaMakeType4(double a, double b, double c, double d)
{
  return make_double4(a, b, c, d); 
};


void
init_atomic_cuda()
{
  fprintf (stderr, "Initializing B-spline matrix.\n");
  if (atomic_cuda_initialized)
    return;
  atomic_cuda_initialized = true;
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
  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);
}


template<typename T, int BS> __global__ void
MakeHybridJobList_kernel (T* elec_list, int num_elecs, T* ion_list,
                          T* poly_radii, T* spline_radii,
                          int num_ions, T *L, T *Linv,
                          HybridJobType *job_list, T *rhat_list,
                          HybridData<T> *data_list)
{
  __shared__ T epos[BS][3], ipos[BS][3], L_s[3][3], Linv_s[3][3];
  __shared__ T rhat[BS][3];
  __shared__ T r_spline[BS], r_poly[BS];
  __shared__ HybridJobType jobs_shared[BS];
  int tid = threadIdx.x;
  if (tid < 9)
  {
    L_s[0][tid]    = L[tid];
    Linv_s[0][tid] = Linv[tid];
  }
  // Load electron positions
  for (int i=0; i<3; i++)
    if ((3*blockIdx.x+i)*BS + tid < 3*num_elecs)
      epos[0][i*BS+tid] = elec_list[(3*blockIdx.x+i)*BS + tid];
  int iBlocks = (num_ions+BS-1)/BS;
  jobs_shared[tid] = BSPLINE_3D_JOB;
  __shared__ HybridData<T> data[BS];
  __syncthreads();
  for (int ib=0; ib<iBlocks; ib++)
  {
    // Fetch ion positions into shared memory
    for (int j=0; j<3; j++)
    {
      int off = (3*ib+j)*BS+tid;
      if (off < 3*num_ions)
        ipos[0][j*BS+tid] = ion_list[off];
    }
    // Fetch radii into shared memory
    if (ib*BS+tid < num_ions)
    {
      r_spline[tid] = spline_radii[ib*BS+tid];
      r_poly[tid]   = poly_radii  [ib*BS+tid];
    }
    __syncthreads();
    int iend = min (BS, num_ions - ib*BS);
    T dr0, dr1, dr2, u0, u1, u2, img0, img1, img2;
    rhat[tid][0] = 1.0f;
    rhat[tid][1] = 0.0f;
    rhat[tid][2] = 0.0f;
    for (int ion=0; ion<iend; ion++)
    {
      // Find mininum image displacement
      dr0 = epos[tid][0] - ipos[ion][0];
      dr1 = epos[tid][1] - ipos[ion][1];
      dr2 = epos[tid][2] - ipos[ion][2];
      u0 = Linv_s[0][0]*dr0 + Linv_s[1][0]*dr1 + Linv_s[2][0]*dr2;
      u1 = Linv_s[0][1]*dr0 + Linv_s[1][1]*dr1 + Linv_s[2][1]*dr2;
      u2 = Linv_s[0][2]*dr0 + Linv_s[1][2]*dr1 + Linv_s[2][2]*dr2;
      img0 = rintf(u0);
      img1 = rintf(u1);
      img2 = rintf(u2);
      u0  -= img0;
      u1  -= img1;
      u2  -= img2;
      dr0 = L_s[0][0]*u0 + L_s[1][0]*u1 + L_s[2][0]*u2;
      dr1 = L_s[0][1]*u0 + L_s[1][1]*u1 + L_s[2][1]*u2;
      dr2 = L_s[0][2]*u0 + L_s[1][2]*u1 + L_s[2][2]*u2;
      T dist2 = dr0*dr0 + dr1*dr1 + dr2*dr2;
      T dist = sqrtf(dist2);
      // Compare with radii
      if (dist < r_poly[ion])
        jobs_shared[tid] =  ATOMIC_POLY_JOB;
      else if (dist < r_spline[ion])
        jobs_shared[tid] =  ATOMIC_SPLINE_JOB;
      // Compute rhat
      if (dist < r_spline[ion])
      {
        data[tid].dist = dist;
        data[tid].img[0] = img0;
        data[tid].img[1] = img1;
        data[tid].img[2] = img2;
        data[tid].ion = ion;
        dist = 1.0f/dist;
        rhat[tid][0] =  dr0 * dist;
        rhat[tid][1] =  dr1 * dist;
        rhat[tid][2] =  dr2 * dist;
      }
    }
  }
  __syncthreads();
  // Now write rhats and job types to global memory
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x+i)*BS + tid;
    if (off < 3*num_elecs)
      rhat_list[off] = rhat[0][i*BS+tid];
  }
  if (blockIdx.x*BS+tid < num_elecs)
    job_list[blockIdx.x*BS+tid] = jobs_shared[tid];
  const int m = (sizeof(HybridData<T>)+sizeof(T)-1)/sizeof(T);
  T *data_f = (T*)data_list;
  for (int i=0; i<m; i++)
  {
    int off = (blockIdx.x*m+i)*BS + tid;
    if (off < m*num_elecs)
      data_f[off] = ((T*)data)[i*BS+tid];
  }
}


template<typename T> void
MakeHybridJobList (T* elec_list, int num_elecs, T* ion_list,
                   T* poly_radii, T* spline_radii,
                   int num_ions, T *L, T *Linv,
                   HybridJobType *job_list, T *rhat_list,
                   HybridData<T> *data_list)
{
  const int BS=32;
  int numBlocks = (num_elecs+BS-1)/BS;
  dim3 dimGrid(numBlocks);
  dim3 dimBlock(BS);
  MakeHybridJobList_kernel<T,BS><<<dimGrid,dimBlock>>>
  (elec_list, num_elecs, ion_list, poly_radii, spline_radii,
   num_ions, L, Linv, job_list, rhat_list, data_list);
}


// The spline coefficients should be reordered so that the
// orbital is the fastest index.
template<typename T, int BS, int LMAX> __global__ void
evaluateHybridSplineReal_kernel (HybridJobType *job_types,
                                 T **YlmReal, AtomicOrbitalCuda<T> *orbitals,
                                 HybridData<T> *data, T *k_reduced,
                                 T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_SPLINE_JOB)
    return;
  __shared__ T *myYlm, *myCoefs, *myVal;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm   = YlmReal[blockIdx.x];
    myVal   = vals[blockIdx.x];
    myCoefs = myOrbital.spline_coefs;
  }
  __syncthreads();
  // Compute spline basis functions
  T unit = myData.dist * myOrbital.spline_dr_inv;
  T sf = floor(unit);
  T t  = unit - sf;
  T v  = 1.0f - t;
  T dr = 1.0f / myOrbital.spline_dr_inv;
  int index= (int) sf;
  // float4 tp;
  // tp = make_float4(t*t*t, t*t, t, 1.0f);
  __shared__ T a[4];
  if (tid == 0)
  {
    a[0]  = v;
    a[1]  = 0.166666666666666666f*(v*v*v-v)*dr*dr;
    a[2]  = t;
    a[3]  = 0.166666666666666666f*(t*t*t-t)*dr*dr;
  }
  // if (tid < 4)
  //   a[tid] = (Acuda[4*tid+0]*tp.x + Acuda[4*tid+1]*tp.y +
  // 	      Acuda[4*tid+2]*tp.z + Acuda[4*tid+3]*tp.w);
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < numlm)
      Ylm[ib*BS+tid] = myYlm[ib*BS + tid];
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int ustride  = 1*myOrbital.spline_stride;
  int ustride2 = 2*myOrbital.spline_stride;
  int ustride3 = 3*myOrbital.spline_stride;
  for (int block=0; block<numBlocks; block++)
  {
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __syncthreads();
    T sign = cos(-(k_red[tid][0]*myData.img[0]+
                   k_red[tid][1]*myData.img[1]+
                   k_red[tid][2]*myData.img[2]));
    T *c0 =  myCoefs + 2*index*ustride + block*BS + tid;
    T val = T();
    for (int lm=0; lm<numlm; lm++)
    {
      T *c = c0 + lm*myOrbital.lm_stride;
      T u = (a[0] * c[0] +
             a[1] * c[ustride] +
             a[2] * c[ustride2] +
             a[3] * c[ustride3]);
      val +=  u * Ylm[lm];
    }
    int off = block*BS + tid;
    if (off < N)
      myVal[off] = sign*val;
  }
  __syncthreads();
}


// The spline coefficients should be reordered so that the
// orbital is the fastest index.
template<typename T, int BS, int LMAX> __global__ void
evaluateHybridPolyReal_kernel (HybridJobType *job_types,
                               T **YlmReal, AtomicOrbitalCuda<T> *orbitals,
                               HybridData<T> *data, T *k_reduced,
                               T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_POLY_JOB)
    return;
  __shared__ T *myYlm, *myCoefs, *myVal;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm   = YlmReal[blockIdx.x];
    myVal   = vals[blockIdx.x];
    myCoefs = myOrbital.poly_coefs;
  }
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < numlm)
      Ylm[ib*BS+tid] = myYlm[ib*BS + tid];
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  __shared__ T r2n[16];
  if (tid < 16)
    r2n[tid] = pow(myData.dist, (T)tid);
  for (int block=0; block<numBlocks; block++)
  {
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __syncthreads();
    T sign = cos(-(k_red[tid][0]*myData.img[0]+
                   k_red[tid][1]*myData.img[1]+
                   k_red[tid][2]*myData.img[2]));
    T *c0 =  myCoefs + block*BS + tid;
    T val = T();
    for (int lm=0; lm<numlm; lm++)
    {
      T *c = c0 + lm*myOrbital.lm_stride;
      T u = 0.0f;
      for (int n=0; n<=myOrbital.poly_order; n++)
        u += r2n[n] * c[n*myOrbital.poly_stride];
      val +=  u * Ylm[lm];
    }
    int off = block*BS + tid;
    if (off < N)
      myVal[off] = sign*val;
  }
  __syncthreads();
}


template<typename T> void
evaluateHybridSplineReal (HybridJobType *job_types,
                          T **Ylm_real, AtomicOrbitalCuda<T> *orbitals,
                          HybridData<T> *data, T *k_reduced,
                          T **vals, int N, int numWalkers, int lMax)
{
  const int BS=32;
  dim3 dimGrid(numWalkers);
  dim3 dimBlock(BS);
  if (lMax == 0)
  {
    evaluateHybridSplineReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 1)
  {
    evaluateHybridSplineReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 2)
  {
    evaluateHybridSplineReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 3)
  {
    evaluateHybridSplineReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 4)
  {
    evaluateHybridSplineReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 5)
  {
    evaluateHybridSplineReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 6)
  {
    evaluateHybridSplineReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 7)
  {
    evaluateHybridSplineReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 8)
  {
    evaluateHybridSplineReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
  else if (lMax == 9)
  {
    evaluateHybridSplineReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
    evaluateHybridPolyReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, vals, N);
  }
}


template<typename T,int BS,int LMAX> __global__ void
evaluateHybridSplineReal_kernel (HybridJobType *job_types, T* rhats,
                                 T **YlmReal, T **dYlm_dTheta, T **dYlm_dphi,
                                 AtomicOrbitalCuda<T> *orbitals, HybridData<T> *data,
                                 T *k_reduced, T **vals, T **grad_lapl,
                                 int row_stride, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_SPLINE_JOB)
    return;
  __shared__ T *myYlm, *mydTheta, *mydPhi, *myCoefs, *myVal, *myGradLapl;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm      = YlmReal[blockIdx.x];
    mydTheta   = dYlm_dTheta[blockIdx.x];
    mydPhi     = dYlm_dphi[blockIdx.x];
    myVal      = vals[blockIdx.x];
    myGradLapl = grad_lapl[blockIdx.x];
    myCoefs    = myOrbital.spline_coefs;
  }
  __shared__ T rhat[3], thetahat[3], phihat[3];
  __shared__ T sintheta, cosphi, sinphi, rInv, sinthetaInv;
  if (tid < 3)
    rhat[tid] = rhats[3*blockIdx.x+tid];
  if (tid ==0)
  {
    rInv = 1.0f/myData.dist;
    sintheta = sqrtf(1.0-rhat[2]*rhat[2]);
    sinthetaInv = 1.0/sintheta;
    cosphi = rhat[0]*sinthetaInv;
    sinphi = rhat[1]*sinthetaInv;
    thetahat[0] = rhat[2]*cosphi;
    thetahat[1] = rhat[2]*sinphi;
    thetahat[2] = -sintheta;
    phihat[0]   = -sinphi;
    phihat[1]   = cosphi;
    phihat[2]   = 0.0f;
  }
  __syncthreads();
  // Compute spline basis functions
  T unit = myData.dist * myOrbital.spline_dr_inv;
  T sf = floor(unit);
  T t  = unit - sf;
  T v  = 1.0f - t;
  int index= (int) sf;
  T dr = 1.0f / myOrbital.spline_dr_inv;
  // float4 tp;
  // tp = make_float4(t*t*t, t*t, t, 1.0f);
  __shared__ T a[12];
  if (tid == 0)
  {
    a[0]  = v;
    a[1]  = 0.166666666666666666f*(v*v*v-v)*dr*dr;
    a[2]  = t;
    a[3]  = 0.166666666666666666f*(t*t*t-t)*dr*dr;
    a[4]  = -myOrbital.spline_dr_inv;
    a[5]  = -0.16666666666666666f * dr * (3.0*v*v-1.0);
    a[6]  =  myOrbital.spline_dr_inv;
    a[7]  =  0.16666666666666666f * dr * (3.0*t*t-1.0);
    a[8]  = 0.0f;
    a[9]  = v;
    a[10] = 0.0f;
    a[11] = t;
  }
  // if (tid < 12)
  //   a[tid] = (Acuda[4*tid+0]*tp.x + Acuda[4*tid+1]*tp.y +
  // 	      Acuda[4*tid+2]*tp.z + Acuda[4*tid+3]*tp.w);
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)], dTheta[(LMAX+1)*(LMAX+1)], dPhi[(LMAX+1)*(LMAX+1)],
             lpref[(LMAX+1)*(LMAX+1)];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < numlm)
    {
      Ylm[ib*BS+tid]    = myYlm[ib*BS + tid];
      dTheta[ib*BS+tid] = mydTheta[ib*BS + tid];
      dPhi[ib*BS+tid]   = mydPhi[ib*BS + tid];
    }
  for (int l=0; l<=myOrbital.lMax; l++)
    if (tid < 2*l+1)
    {
      int lm = l*l + tid;
      lpref[lm] = -rInv*rInv*(T)(l*(l+1));
    }
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int ustride  = 1*myOrbital.spline_stride;
  int ustride2 = 2*myOrbital.spline_stride;
  int ustride3 = 3*myOrbital.spline_stride;
  for (int block=0; block<numBlocks; block++)
  {
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __syncthreads();
    T sign = cos(-(k_red[tid][0]*myData.img[0]+
                      k_red[tid][1]*myData.img[1]+
                      k_red[tid][2]*myData.img[2]));
    T *c0 =  myCoefs + 2*index*ustride + block*BS + tid;
    T val = T();
    T g_rhat=T(), g_thetahat=T(), g_phihat=T(), lap=T();
    for (int lm=0; lm<numlm; lm++)
    {
      T *c = c0 + lm*myOrbital.lm_stride;
      T u, du, d2u, coef;
      coef = c[0];
      // Evaluate spline value and derivatives
      u  = a[0]*coef;
      du  = a[4]*coef;
      d2u  = a[8]*coef;
      coef = c[ustride];
      u += a[1]*coef;
      du += a[5]*coef;
      d2u += a[9]*coef;
      coef = c[ustride2];
      u += a[2]*coef;
      du += a[6]*coef;
      d2u += a[10]*coef;
      coef = c[ustride3];
      u += a[3]*coef;
      du += a[7]*coef;
      d2u += a[11]*coef;
      // du  *= myOrbital.spline_dr_inv;
      // d2u *= myOrbital.spline_dr_inv * myOrbital.spline_dr_inv;
      // Now accumulate values
      val    +=   u * Ylm[lm];
      g_rhat +=  du*Ylm[lm];
      g_thetahat += u*rInv*dTheta[lm];
      g_phihat   += u*rInv*sinthetaInv*dPhi[lm];
      lap += Ylm[lm] * (lpref[lm] * u + d2u + 2.0f*rInv*du);
    }
    int off = block*BS + tid;
    if (off < N)
    {
      myVal[off] = sign * val;
      myGradLapl[0*row_stride+off] = sign *
                                     (g_rhat*rhat[0]+g_thetahat*thetahat[0]+g_phihat*phihat[0]);
      myGradLapl[1*row_stride+off] = sign *
                                     (g_rhat*rhat[1]+g_thetahat*thetahat[1]+g_phihat*phihat[1]);
      myGradLapl[2*row_stride+off] = sign *
                                     (g_rhat*rhat[2]+g_thetahat*thetahat[2]+g_phihat*phihat[2]);
      myGradLapl[3*row_stride+off] = sign * lap;
    }
  }
  __syncthreads();
}


template<typename T,int BS,int LMAX> __global__ void
evaluateHybridPolyReal_kernel (HybridJobType *job_types, T* rhats,
                               T **YlmReal, T **dYlm_dTheta, T **dYlm_dphi,
                               AtomicOrbitalCuda<T> *orbitals, HybridData<T> *data,
                               T *k_reduced, T **vals, T **grad_lapl,
                               int row_stride, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_POLY_JOB)
    return;
  __shared__ T *myYlm, *mydTheta, *mydPhi, *myCoefs, *myVal, *myGradLapl;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm      = YlmReal[blockIdx.x];
    mydTheta   = dYlm_dTheta[blockIdx.x];
    mydPhi     = dYlm_dphi[blockIdx.x];
    myVal      = vals[blockIdx.x];
    myGradLapl = grad_lapl[blockIdx.x];
    myCoefs    = myOrbital.poly_coefs;
  }
  __shared__ T rhat[3], thetahat[3], phihat[3];
  __shared__ T sintheta, cosphi, sinphi, rInv, sinthetaInv;
  if (tid < 3)
    rhat[tid] = rhats[3*blockIdx.x+tid];
  if (tid ==0)
  {
    rInv = 1.0f/myData.dist;
    sintheta = sqrtf(1.0-rhat[2]*rhat[2]);
    sinthetaInv = 1.0/sintheta;
    cosphi = rhat[0]*sinthetaInv;
    sinphi = rhat[1]*sinthetaInv;
    thetahat[0] = rhat[2]*cosphi;
    thetahat[1] = rhat[2]*sinphi;
    thetahat[2] = -sintheta;
    phihat[0]   = -sinphi;
    phihat[1]   = cosphi;
    phihat[2]   = 0.0f;
  }
  __syncthreads();
  // Compute polynomial basis functions
  // Note maximum polynomial order of 16
  __shared__ T polyfuncs[16][3];
  // if (tid == 0) {
  //   polyfuncs[0][0] = 1.0f; polyfuncs[1][0] = myData.dist;
  //   polyfuncs[0][1] = 0.0f; polyfuncs[1][1] = 1.0f;
  //   polyfuncs[0][2] = 0.0f; polyfuncs[1][2] = 0.0f;
  //   float dn = 2.0f;
  //   for (int n=2; n<=myOrbital.poly_order; n++) {
  //     polyfuncs[n][0] = myData.dist*polyfuncs[n-1][0];
  //     polyfuncs[n][1] = dn*polyfuncs[n-1][0];
  //     polyfuncs[n][2] = dn*(dn-1.0f)*polyfuncs[n-2][0];
  //     dn += 1.0f;
  //   }
  // }
  if (tid < 16)
  {
    polyfuncs[tid][0] = pow(myData.dist,(T)tid);
    polyfuncs[tid][1] = (T)tid    *polyfuncs[tid][0] / myData.dist;
    polyfuncs[tid][2] = (T)(tid-1)*polyfuncs[tid][1] / myData.dist;
  }
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)], dTheta[(LMAX+1)*(LMAX+1)], dPhi[(LMAX+1)*(LMAX+1)],
             lpref[(LMAX+1)*(LMAX+1)];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < numlm)
    {
      Ylm[ib*BS+tid]    = myYlm[ib*BS + tid];
      dTheta[ib*BS+tid] = mydTheta[ib*BS + tid];
      dPhi[ib*BS+tid]   = mydPhi[ib*BS + tid];
    }
  for (int l=0; l<=myOrbital.lMax; l++)
    if (tid < 2*l+1)
    {
      int lm = l*l + tid;
      lpref[lm] = -rInv*rInv*(T)(l*(l+1));
    }
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  for (int block=0; block<numBlocks; block++)
  {
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __syncthreads();
    T sign = cos(-(k_red[tid][0]*myData.img[0]+
                   k_red[tid][1]*myData.img[1]+
                   k_red[tid][2]*myData.img[2]));
    T *c0 =  myCoefs + block*BS + tid;
    T val = T();
    T g_rhat=T(), g_thetahat=T(), g_phihat=T(), lap=T();
    for (int lm=0; lm<numlm; lm++)
    {
      T *c = c0 + lm*myOrbital.lm_stride;
      T u=0.0f, du=0.0f, d2u=0.0f, coef;
      for (int n=0; n<=myOrbital.poly_order; n++)
      {
        coef = c[n*myOrbital.poly_stride];
        u   += coef * polyfuncs[n][0];
        du  += coef * polyfuncs[n][1];
        d2u += coef * polyfuncs[n][2];
      }
      // Now accumulate values
      val    +=   u * Ylm[lm];
      g_rhat +=  du*Ylm[lm];
      g_thetahat += u*rInv*dTheta[lm];
      g_phihat   += u*rInv*sinthetaInv*dPhi[lm];
      lap += Ylm[lm] * (lpref[lm] * u + d2u + 2.0f*rInv*du);
    }
    int off = block*BS + tid;
    if (off < N)
    {
      myVal[off] = sign * val;
      myGradLapl[0*row_stride+off] = sign *
                                     (g_rhat*rhat[0]+g_thetahat*thetahat[0]+g_phihat*phihat[0]);
      myGradLapl[1*row_stride+off] = sign *
                                     (g_rhat*rhat[1]+g_thetahat*thetahat[1]+g_phihat*phihat[1]);
      myGradLapl[2*row_stride+off] = sign *
                                     (g_rhat*rhat[2]+g_thetahat*thetahat[2]+g_phihat*phihat[2]);
      myGradLapl[3*row_stride+off] = sign * lap;
    }
  }
  __syncthreads();
}


template<typename T> void
evaluateHybridSplineReal (HybridJobType *job_types, T *rhats,
                          T **Ylm_real, T **dYlm_dTheta, T **dYlm_dphi,
                          AtomicOrbitalCuda<T> *orbitals,
                          HybridData<T> *data, T *k_reduced,
                          T **vals, T **grad_lapl,
                          int row_stride, int N, int numWalkers, int lMax)
{
  const int BS=32;
  dim3 dimGrid(numWalkers);
  dim3 dimBlock(BS);
  if (lMax == 0)
  {
    evaluateHybridSplineReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 1)
  {
    evaluateHybridSplineReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 2)
  {
    evaluateHybridSplineReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 3)
  {
    evaluateHybridSplineReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 4)
  {
    evaluateHybridSplineReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 5)
  {
    evaluateHybridSplineReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 6)
  {
    evaluateHybridSplineReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 7)
  {
    evaluateHybridSplineReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 8)
  {
    evaluateHybridSplineReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 9)
  {
    evaluateHybridSplineReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm_real, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, vals, grad_lapl, row_stride, N);
  }
}


/************************************************************
*************************************************************
//// Complex splines to real orbitals evaluation routines ////
************************************************************
************************************************************/



// The spline coefficients should be reordered so that the
// orbital is the fastest index.
template<typename T, int BS, int LMAX> __global__ void
evaluateHybridSplineComplexToReal_kernel
(HybridJobType *job_types,
 T **YlmComplex, AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int *make2copies,
 T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_SPLINE_JOB)
    return;
  __shared__ T *myYlm, *myCoefs, *myVal;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm   = YlmComplex[blockIdx.x];
    myVal   = vals[blockIdx.x];
    myCoefs = myOrbital.spline_coefs;
  }
  __syncthreads();
  // Compute spline basis functions
  T unit = myData.dist * myOrbital.spline_dr_inv;
  T sf = floor(unit);
  T t  = unit - sf;
  T v  = 1.0f - t;
  int index= (int) sf;
  // float4 tp;
  // tp = make_float4(t*t*t, t*t, t, 1.0f);
  __shared__ T a[4];
  T dr = 1.0f / myOrbital.spline_dr_inv;
  if (tid == 0)
  {
    a[0]  = v;
    a[1]  = 0.166666666666666666f*(v*v*v-v)*dr*dr;
    a[2]  = t;
    a[3]  = 0.166666666666666666f*(t*t*t-t)*dr*dr;
  }
  // if (tid < 4)
  //   a[tid] = (Acuda[4*tid+0]*tp.x + Acuda[4*tid+1]*tp.y +
  // 	      Acuda[4*tid+2]*tp.z + Acuda[4*tid+3]*tp.w);
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)][2];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (2*numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < 2*numlm)
      Ylm[0][ib*BS+tid] = myYlm[ib*BS + tid];
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int ustride  = 1*myOrbital.spline_stride;
  int ustride2 = 2*myOrbital.spline_stride;
  int ustride3 = 3*myOrbital.spline_stride;
  int outIndex=0, outBlock=0;
  __shared__ T outbuff[BS];
  for (int block=0; block<numBlocks; block++)
  {
    T phase_re, phase_im;
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __shared__ int m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    sincos(-(k_red[tid][0]*myData.img[0]+
             k_red[tid][1]*myData.img[1]+
             k_red[tid][2]*myData.img[2]),
             &phase_im, &phase_re);
    T *c0 =  myCoefs + 2*index*ustride + 2*block*BS;
    __shared__ T c[BS][2], val[BS][2];
    val[0][tid]    = T();
    val[0][tid+BS] = T();
    T v_re=T(), v_im=T();
    for (int lm=0; lm<numlm; lm++)
    {
      T u_re, u_im;
      T *coef = c0 + lm*myOrbital.lm_stride;
      c[0][tid]    = coef[tid];
      c[0][BS+tid] = coef[BS+tid];
      u_re = a[0]*c[tid][0];
      u_im = a[0]*c[tid][1];
      c[0][tid]    = coef[ustride      + tid];
      c[0][BS+tid] = coef[ustride + BS + tid];
      u_re += a[1]*c[tid][0];
      u_im += a[1]*c[tid][1];
      c[0][tid]    = coef[ustride2      + tid];
      c[0][BS+tid] = coef[ustride2 + BS + tid];
      u_re += a[2]*c[tid][0];
      u_im += a[2]*c[tid][1];
      c[0][tid]    = coef[ustride3      + tid];
      c[0][BS+tid] = coef[ustride3 + BS + tid];
      u_re += a[3]*c[tid][0];
      u_im += a[3]*c[tid][1];
      v_re +=  u_re*Ylm[lm][0] - u_im*Ylm[lm][1];
      v_im +=  u_im*Ylm[lm][0] + u_re*Ylm[lm][1];
    }
    val[tid][0] = phase_re * v_re - phase_im * v_im;
    val[tid][1] = phase_re * v_im + phase_im * v_re;
    __syncthreads();
    int iend = min (BS, N-block*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == 0)
        outbuff[outIndex] = val[i][0];
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        myVal[outBlock++*BS+tid] = outbuff[tid];
        outIndex = 0;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == 0)
          outbuff[outIndex] = val[i][1];
        outIndex++;
        __syncthreads();
      }
      if (outIndex == BS)
      {
        myVal[outBlock++*BS+tid] = outbuff[tid];
        outIndex = 0;
      }
      __syncthreads();
    }
    // int off = block*BS + tid;
    // if (off < N)
    //   myVal[off] = phase_re*val;
  }
  if (tid < outIndex)
    myVal[outBlock*BS+tid] = outbuff[tid];
  __syncthreads();
}

// This version is used for pseudopotential quadrature.  If there are
// 12 quadrature points, the radial parts only need to be evaluated
// once.
template<typename T, int BS, int LMAX> __global__ void
evaluateHybridSplineComplexToReal_NLPP_kernel
(HybridJobType *job_types, T **YlmComplex, int numQuad,
 AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int *make2copies,
 T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[numQuad*blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_SPLINE_JOB)
    return;
  __shared__ T *myYlm[MaxQuad], *myCoefs, *myVal[MaxQuad];
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[numQuad*blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myCoefs = myOrbital.spline_coefs;
  }
  if (tid < numQuad)
  {
    myYlm[tid] = YlmComplex[numQuad*blockIdx.x+tid];
    myVal[tid] = vals[numQuad*blockIdx.x+tid];
  }
  __syncthreads();
  // Compute spline basis functions
  T unit = myData.dist * myOrbital.spline_dr_inv;
  T sf = floor(unit);
  T t  = unit - sf;
  T v  = 1.0f - t;
  int index= (int) sf;
  __shared__ T a[4];
  T dr = 1.0f / myOrbital.spline_dr_inv;
  if (tid == 0)
  {
    a[0]  = v;
    a[1]  = 0.166666666666666666f*(v*v*v-v)*dr*dr;
    a[2]  = t;
    a[3]  = 0.166666666666666666f*(t*t*t-t)*dr*dr;
  }
  __syncthreads();
  __shared__ T Ylm[MaxQuad][(LMAX+1)*(LMAX+1)][2];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (2*numlm+BS-1)/BS;
  for (int iq=0; iq<numQuad; iq++)
  {
    for (int ib=0; ib<Yblocks; ib++)
      if (ib*BS + tid < 2*numlm)
        Ylm[iq][0][ib*BS+tid] = myYlm[iq][ib*BS + tid];
  }
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int ustride  = 1*myOrbital.spline_stride;
  int ustride2 = 2*myOrbital.spline_stride;
  int ustride3 = 3*myOrbital.spline_stride;
  int outIndex=0, outBlock=0;
  __shared__ T outbuff[MaxQuad][BS];
  for (int block=0; block<numBlocks; block++)
  {
    T phase_re, phase_im;
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __shared__ int m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    sincos(-(k_red[tid][0]*myData.img[0]+
             k_red[tid][1]*myData.img[1]+
             k_red[tid][2]*myData.img[2]),
           &phase_im, &phase_re);
    T *c0 =  myCoefs + 2*index*ustride + 2*block*BS;
    __shared__ T c[BS][2], val[BS][2];
    val[0][tid]    = T();
    val[0][tid+BS] = T();
    __shared__ T v[MaxQuad][2][BS];
    for (int iq=0; iq<MaxQuad; iq++)
    {
      v[iq][0][tid] = T();
      v[iq][1][tid] = T();
    }
    for (int lm=0; lm<numlm; lm++)
    {
      T u_re, u_im;
      T *coef = c0 + lm*myOrbital.lm_stride;
      c[0][tid]    = coef[tid];
      c[0][BS+tid] = coef[BS+tid];
      u_re = a[0]*c[tid][0];
      u_im = a[0]*c[tid][1];
      c[0][tid]    = coef[ustride      + tid];
      c[0][BS+tid] = coef[ustride + BS + tid];
      u_re += a[1]*c[tid][0];
      u_im += a[1]*c[tid][1];
      c[0][tid]    = coef[ustride2      + tid];
      c[0][BS+tid] = coef[ustride2 + BS + tid];
      u_re += a[2]*c[tid][0];
      u_im += a[2]*c[tid][1];
      c[0][tid]    = coef[ustride3      + tid];
      c[0][BS+tid] = coef[ustride3 + BS + tid];
      u_re += a[3]*c[tid][0];
      u_im += a[3]*c[tid][1];
      for (int iq=0; iq<numQuad; iq++)
      {
        v[iq][0][tid] +=  u_re*Ylm[iq][lm][0] - u_im*Ylm[iq][lm][1];
        v[iq][1][tid] +=  u_im*Ylm[iq][lm][0] + u_re*Ylm[iq][lm][1];
      }
    }
    for (int iq=0; iq<numQuad; iq++)
    {
      int saveIndex = outIndex;
      int saveBlock = outBlock;
      val[tid][0] = phase_re * v[iq][0][tid] - phase_im * v[iq][1][tid];
      val[tid][1] = phase_re * v[iq][1][tid] + phase_im * v[iq][0][tid];
      __syncthreads();
      int iend = min (BS, N-block*BS);
      for (int i=0; i<iend; i++)
      {
        if (tid == 0)
          outbuff[iq][outIndex] = val[i][0];
        outIndex++;
        __syncthreads();
        if (outIndex == BS)
        {
          myVal[iq][outBlock++*BS+tid] = outbuff[iq][tid];
          outIndex = 0;
        }
        __syncthreads();
        if (m2c[i])
        {
          if (tid == 0)
            outbuff[iq][outIndex] = val[i][1];
          outIndex++;
          __syncthreads();
        }
        if (outIndex == BS)
        {
          myVal[iq][outBlock++*BS+tid] = outbuff[iq][tid];
          outIndex = 0;
        }
        __syncthreads();
      }
      // If it's not the last quad point, restore value from
      // the beginning of the iq loop iteration.
      if (iq < (numQuad-1))
      {
        outIndex = saveIndex;
        outBlock = saveBlock;
      }
    }
    // int off = block*BS + tid;
    // if (off < N)
    //   myVal[off] = phase_re*val;
  }
  for (int iq=0; iq<numQuad; iq++)
    if (tid < outIndex)
      myVal[iq][outBlock*BS+tid] = outbuff[iq][tid];
  __syncthreads();
}


// The spline coefficients should be reordered so that the
// orbital is the fastest index.
template<typename T, int BS, int LMAX> __global__ void
evaluateHybridPolyComplexToReal_kernel
(HybridJobType *job_types,
 T **YlmComplex, AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int *make2copies,
 T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_POLY_JOB)
    return;
  __shared__ T *myYlm, *myCoefs, *myVal;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm   = YlmComplex[blockIdx.x];
    myVal   = vals[blockIdx.x];
    myCoefs = myOrbital.poly_coefs;
  }
  __syncthreads();
  __shared__ T r2n[16];
  if (tid < 16)
    r2n[tid] = pow(myData.dist, (T)tid);
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)][2];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (2*numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < 2*numlm)
      Ylm[0][ib*BS+tid] = myYlm[ib*BS + tid];
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int outIndex=0, outBlock=0;
  __shared__ T outbuff[BS];
  for (int block=0; block<numBlocks; block++)
  {
    T phase_re, phase_im;
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __shared__ int m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    sincos(-(k_red[tid][0]*myData.img[0]+
             k_red[tid][1]*myData.img[1]+
             k_red[tid][2]*myData.img[2]),
           &phase_im, &phase_re);
    T *c0 =  myCoefs + 2*block*BS + tid;
    __shared__ T c[BS][2], val[BS][2];
    val[0][tid]    = T();
    val[0][tid+BS] = T();
    T v_re=T(), v_im=T();
    for (int lm=0; lm<numlm; lm++)
    {
      T u_re=0.0f, u_im=0.0f;
      T *coef = c0 + lm*myOrbital.lm_stride;
      for (int n=0; n<=myOrbital.poly_order; n++)
      {
        c[0][tid]    = coef[n*myOrbital.poly_stride];
        c[0][BS+tid] = coef[n*myOrbital.poly_stride+BS];
        __syncthreads();
        u_re += r2n[n]*c[tid][0];
        u_im += r2n[n]*c[tid][1];
      }
      v_re +=  u_re*Ylm[lm][0] - u_im*Ylm[lm][1];
      v_im +=  u_im*Ylm[lm][0] + u_re*Ylm[lm][1];
    }
    val[tid][0] = phase_re * v_re - phase_im * v_im;
    val[tid][1] = phase_re * v_im + phase_im * v_re;
    __syncthreads();
    int iend = min (BS, N-block*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == 0)
        outbuff[outIndex] = val[i][0];
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        myVal[outBlock++*BS+tid] = outbuff[tid];
        outIndex = 0;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == 0)
          outbuff[outIndex] = val[i][1];
        outIndex++;
        __syncthreads();
      }
      if (outIndex == BS)
      {
        myVal[outBlock++*BS+tid] = outbuff[tid];
        outIndex = 0;
      }
      __syncthreads();
    }
    // int off = block*BS + tid;
    // if (off < N)
    //   myVal[off] = phase_re*val;
  }
  if (tid < outIndex)
    myVal[outBlock*BS+tid] = outbuff[tid];
  __syncthreads();
}


template<typename T> void
evaluateHybridSplineComplexToReal
(HybridJobType *job_types, T **Ylm_real,
 AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int* make2copies,
 T **vals, int N, int numWalkers, int lMax)
{
  const int BS=32;
  dim3 dimGrid(numWalkers);
  dim3 dimBlock(BS);
  if (lMax == 0)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 1)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 2)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 3)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 4)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 5)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 6)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 7)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 8)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 9)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
}

// Ying Wai: this seems to be unused.  (Oct 1, 15)
/*
template<typename T> void
evaluateHybridSplineComplexToRealNLPP
(HybridJobType *job_types,
 T **Ylm_real, AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int* make2copies,
 T **vals, int N, int numWalkers, int lMax)
{
  const int BS=32;
  const int numQuad = 12;
  dim3 polyGrid(numWalkers), NLPPGrid((numWalkers+numQuad-1)/numQuad);
  dim3 dimBlock(BS);
  if (lMax == 0)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,0><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,0><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 1)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,1><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,1><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 2)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,2><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,2><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 3)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,3><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,3><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 4)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,4><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,4><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 5)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,5><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,5><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 6)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,6><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,6><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 7)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,7><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,7><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 8)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,8><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,8><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
  else if (lMax == 9)
  {
    evaluateHybridSplineComplexToReal_NLPP_kernel<T,BS,9><<<NLPPGrid,dimBlock>>>
    (job_types, Ylm_real, numQuad, orbitals, data, k_reduced, make2copies, vals, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,9><<<polyGrid,dimBlock>>>
    (job_types, Ylm_real, orbitals, data, k_reduced, make2copies, vals, N);
  }
}
*/

template<typename T,int BS,int LMAX> __global__ void
evaluateHybridSplineComplexToReal_kernel
(HybridJobType *job_types, T* rhats,
 T **Ylm_complex, T **dYlm_dTheta, T **dYlm_dphi,
 AtomicOrbitalCuda<T> *orbitals, HybridData<T> *data,
 T *k_reduced, int *make2copies, T **vals, T **grad_lapl,
 int row_stride, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_SPLINE_JOB)
    return;
  __shared__ T *myYlm, *mydTheta, *mydPhi, *myCoefs, *myVal, *myGradLapl;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm      = Ylm_complex[blockIdx.x];
    mydTheta   = dYlm_dTheta[blockIdx.x];
    mydPhi     = dYlm_dphi[blockIdx.x];
    myVal      = vals[blockIdx.x];
    myGradLapl = grad_lapl[blockIdx.x];
    myCoefs    = myOrbital.spline_coefs;
  }
  __shared__ T rhat[3], thetahat[3], phihat[3];
  __shared__ T sintheta, cosphi, sinphi, rInv, sinthetaInv;
  if (tid < 3)
    rhat[tid] = rhats[3*blockIdx.x+tid];
  if (tid ==0)
  {
    rInv = 1.0f/myData.dist;
    sintheta = sqrtf(1.0-rhat[2]*rhat[2]);
    sinthetaInv = 1.0/sintheta;
    cosphi = rhat[0]*sinthetaInv;
    sinphi = rhat[1]*sinthetaInv;
    thetahat[0] = rhat[2]*cosphi;
    thetahat[1] = rhat[2]*sinphi;
    thetahat[2] = -sintheta;
    phihat[0]   = -sinphi;
    phihat[1]   = cosphi;
    phihat[2]   = 0.0f;
  }
  __syncthreads();
  // Compute spline basis functions
  T unit = myData.dist * myOrbital.spline_dr_inv;
  T sf = floor(unit);
  T t  = unit - sf;
  T v = 1.0f - t;
  int index= (int) sf;
  // float4 tp;
  // tp = make_float4(t*t*t, t*t, t, 1.0f);
  __shared__ T a[12];
  T dr = 1.0f / myOrbital.spline_dr_inv;
  if (tid == 0)
  {
    a[0]  = v;
    a[1]  = 0.166666666666666666f*(v*v*v-v)*dr*dr;
    a[2]  = t;
    a[3]  = 0.166666666666666666f*(t*t*t-t)*dr*dr;
    a[4]  = -myOrbital.spline_dr_inv;
    a[5]  = -0.16666666666666666f * dr * (3.0*v*v-1.0);
    a[6]  =  myOrbital.spline_dr_inv;
    a[7]  =  0.16666666666666666f * dr * (3.0*t*t-1.0);
    a[8]  = 0.0f;
    a[9]  = v;
    a[10] = 0.0f;
    a[11] = t;
  }
  // if (tid < 12)
  //   a[tid] = (Acuda[4*tid+0]*tp.x + Acuda[4*tid+1]*tp.y +
  // 	      Acuda[4*tid+2]*tp.z + Acuda[4*tid+3]*tp.w);
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)][2],
             dTheta[(LMAX+1)*(LMAX+1)][2], dPhi[(LMAX+1)*(LMAX+1)][2],
             lpref[(LMAX+1)*(LMAX+1)];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (2*numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < 2*numlm)
    {
      Ylm[0][ib*BS+tid]    = myYlm[ib*BS + tid];
      dTheta[0][ib*BS+tid] = mydTheta[ib*BS + tid];
      dPhi[0][ib*BS+tid]   = mydPhi[ib*BS + tid];
    }
  for (int l=0; l<=myOrbital.lMax; l++)
    if (tid < 2*l+1)
    {
      int lm = l*l + tid;
      lpref[lm] = -rInv*rInv*(T)(l*(l+1));
    }
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int ustride  = 1*myOrbital.spline_stride;
  int ustride2 = 2*myOrbital.spline_stride;
  int ustride3 = 3*myOrbital.spline_stride;
  int outIndex=0, outBlock=0;
  __shared__ T outval[BS], outgrad[BS][3], outlap[BS];
  for (int block=0; block<numBlocks; block++)
  {
    T phase_re, phase_im;
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __shared__ T m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    sincos(-(k_red[tid][0]*myData.img[0]+
             k_red[tid][1]*myData.img[1]+
             k_red[tid][2]*myData.img[2]),
           &phase_im, &phase_re);
    T *c0 =  myCoefs + 2*index*ustride + 2*block*BS;
    T v_re         =T(), v_im         =T();
    T g_rhat_re    =T(), g_rhat_im    =T(),
      g_thetahat_re=T(), g_thetahat_im=T(),
      g_phihat_re  =T(), g_phihat_im  =T(),
      lap_re       =T(), lap_im       =T();
    __shared__ T c[BS][2];
    for (int lm=0; lm<numlm; lm++)
    {
      T u_re, u_im, du_re, du_im, d2u_re, d2u_im;
      T *coef = c0 + lm*myOrbital.lm_stride;
      c[0][tid]    = coef[tid];
      c[0][BS+tid] = coef[BS+tid];
      u_re   = a[0]*c[tid][0];
      u_im   = a[0]*c[tid][1];
      du_re  = a[4]*c[tid][0];
      du_im  = a[4]*c[tid][1];
      d2u_re = a[8]*c[tid][0];
      d2u_im = a[8]*c[tid][1];
      c[0][tid]    = coef[ustride      + tid];
      c[0][BS+tid] = coef[ustride + BS + tid];
      u_re   += a[1]*c[tid][0];
      u_im   += a[1]*c[tid][1];
      du_re  += a[5]*c[tid][0];
      du_im  += a[5]*c[tid][1];
      d2u_re += a[9]*c[tid][0];
      d2u_im += a[9]*c[tid][1];
      c[0][tid]    = coef[ustride2      + tid];
      c[0][BS+tid] = coef[ustride2 + BS + tid];
      u_re   += a[ 2]*c[tid][0];
      u_im   += a[ 2]*c[tid][1];
      du_re  += a[ 6]*c[tid][0];
      du_im  += a[ 6]*c[tid][1];
      d2u_re += a[10]*c[tid][0];
      d2u_im += a[10]*c[tid][1];
      c[0][tid]    = coef[ustride3      + tid];
      c[0][BS+tid] = coef[ustride3 + BS + tid];
      u_re   += a[ 3]*c[tid][0];
      u_im   += a[ 3]*c[tid][1];
      du_re  += a[ 7]*c[tid][0];
      du_im  += a[ 7]*c[tid][1];
      d2u_re += a[11]*c[tid][0];
      d2u_im += a[11]*c[tid][1];
      // du_re  *= myOrbital.spline_dr_inv;
      // du_im  *= myOrbital.spline_dr_inv;
      // d2u_re *= myOrbital.spline_dr_inv * myOrbital.spline_dr_inv;
      // d2u_im *= myOrbital.spline_dr_inv * myOrbital.spline_dr_inv;
      // Now accumulate values
      v_re    +=   u_re * Ylm[lm][0] - u_im * Ylm[lm][1];
      v_im    +=   u_re * Ylm[lm][1] + u_im * Ylm[lm][0];
      g_rhat_re +=  du_re * Ylm[lm][0] - du_im * Ylm[lm][1];
      g_rhat_im +=  du_re * Ylm[lm][1] + du_im * Ylm[lm][0];
      g_thetahat_re += rInv*(u_re*dTheta[lm][0] - u_im*dTheta[lm][1]);
      g_thetahat_im += rInv*(u_re*dTheta[lm][1] + u_im*dTheta[lm][0]);
      g_phihat_re += rInv*sinthetaInv*
                     (u_re*dPhi[lm][0] - u_im*dPhi[lm][1]);
      g_phihat_im += rInv*sinthetaInv*
                     (u_re*dPhi[lm][1] + u_im*dPhi[lm][0]);
      lap_re += (Ylm[lm][0] * (lpref[lm]*u_re + d2u_re + 2.0f*rInv*du_re) -
                 Ylm[lm][1] * (lpref[lm]*u_im + d2u_im + 2.0f*rInv*du_im));
      lap_im += (Ylm[lm][1] * (lpref[lm]*u_re + d2u_re + 2.0f*rInv*du_re) +
                 Ylm[lm][0] * (lpref[lm]*u_im + d2u_im + 2.0f*rInv*du_im));
    }
    __shared__ T val[BS][2], grad[BS][3][2], lap[BS][2];
    val[tid][0] = phase_re *   v_re - phase_im *   v_im;
    val[tid][1] = phase_re *   v_im + phase_im *   v_re;
    lap[tid][0] = phase_re * lap_re - phase_im * lap_im;
    lap[tid][1] = phase_re * lap_im + phase_im * lap_re;
    for (int i=0; i<3; i++)
    {
      grad[tid][i][0] =
        (phase_re*g_rhat_re     - phase_im*g_rhat_im    )*rhat[i] +
        (phase_re*g_thetahat_re - phase_im*g_thetahat_im)*thetahat[i] +
        (phase_re*g_phihat_re   - phase_im*g_phihat_im  )*phihat[i];
      grad[tid][i][1] =
        (phase_im*g_rhat_re     + phase_re*g_rhat_im    )*rhat[i] +
        (phase_im*g_thetahat_re + phase_re*g_thetahat_im)*thetahat[i] +
        (phase_im*g_phihat_re   + phase_re*g_phihat_im  )*phihat[i];
    }
    // Now serialize to output buffers.
    int iend = min (BS, N-block*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == 0)
      {
        outval[outIndex]     = val[i][0];
        outgrad[outIndex][0] = grad[i][0][0];
        outgrad[outIndex][1] = grad[i][1][0];
        outgrad[outIndex][2] = grad[i][2][0];
        outlap[outIndex]     = lap[i][0];
      }
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        int off = outBlock++*BS+tid;
        myVal[off] = outval[tid];
        myGradLapl[0*row_stride+off] = outgrad[tid][0];
        myGradLapl[1*row_stride+off] = outgrad[tid][1];
        myGradLapl[2*row_stride+off] = outgrad[tid][2];
        myGradLapl[3*row_stride+off] = outlap[tid];
        outIndex = 0;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == 0)
        {
          outval[outIndex]     = val[i][1];
          outgrad[outIndex][0] = grad[i][0][1];
          outgrad[outIndex][1] = grad[i][1][1];
          outgrad[outIndex][2] = grad[i][2][1];
          outlap[outIndex]     = lap[i][1];
        }
        outIndex++;
        __syncthreads();
      }
      if (outIndex == BS)
      {
        int off = outBlock++*BS+tid;
        myVal[off] = outval[tid];
        myGradLapl[0*row_stride+off] = outgrad[tid][0];
        myGradLapl[1*row_stride+off] = outgrad[tid][1];
        myGradLapl[2*row_stride+off] = outgrad[tid][2];
        myGradLapl[3*row_stride+off] = outlap[tid];
        outIndex = 0;
      }
      __syncthreads();
    }
    // int off = block*BS + tid;
    // if (off < N) {
    //   myVal[off] = sign * val;
    //   myGradLapl[0*row_stride+off] = sign *
    // 	(g_rhat*rhat[0]+g_thetahat*thetahat[0]+g_phihat*phihat[0]);
    //   myGradLapl[1*row_stride+off] = sign *
    // 	(g_rhat*rhat[1]+g_thetahat*thetahat[1]+g_phihat*phihat[1]);
    //   myGradLapl[2*row_stride+off] = sign *
    // 	(g_rhat*rhat[2]+g_thetahat*thetahat[2]+g_phihat*phihat[2]);
    //   myGradLapl[3*row_stride+off] = sign * lap;
    // }
  }
  if (tid < outIndex)
  {
    int off = outBlock*BS+tid;
    myVal[off] = outval[tid];
    myGradLapl[0*row_stride+off] = outgrad[tid][0];
    myGradLapl[1*row_stride+off] = outgrad[tid][1];
    myGradLapl[2*row_stride+off] = outgrad[tid][2];
    myGradLapl[3*row_stride+off] = outlap[tid];
  }
}


template<typename T,int BS,int LMAX> __global__ void
evaluateHybridPolyComplexToReal_kernel
(HybridJobType *job_types, T* rhats,
 T **Ylm_complex, T **dYlm_dTheta, T **dYlm_dphi,
 AtomicOrbitalCuda<T> *orbitals, HybridData<T> *data,
 T *k_reduced, int *make2copies, T **vals, T **grad_lapl,
 int row_stride, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != ATOMIC_POLY_JOB)
    return;
  __shared__ T *myYlm, *mydTheta, *mydPhi, *myCoefs, *myVal, *myGradLapl;
  __shared__ HybridData<T> myData;
  __shared__ AtomicOrbitalCuda<T> myOrbital;
  __shared__ T k_red[BS][3];
  const int data_size = (sizeof(HybridData<T>)     +3)/sizeof(T);
  const int orb_size  = (sizeof(AtomicOrbitalCuda<T>)+3)/sizeof(T);
  if (tid < data_size)
    ((T*)&myData)[tid]    = ((T*)&(data[blockIdx.x]))[tid];
  if (tid < orb_size)
    ((T*)&myOrbital)[tid] = ((T*)(&orbitals[myData.ion]))[tid];
  if (tid == 0)
  {
    myYlm      = Ylm_complex[blockIdx.x];
    mydTheta   = dYlm_dTheta[blockIdx.x];
    mydPhi     = dYlm_dphi[blockIdx.x];
    myVal      = vals[blockIdx.x];
    myGradLapl = grad_lapl[blockIdx.x];
    myCoefs    = myOrbital.poly_coefs;
  }
  __shared__ T rhat[3], thetahat[3], phihat[3];
  __shared__ T sintheta, cosphi, sinphi, rInv, sinthetaInv;
  if (tid < 3)
    rhat[tid] = rhats[3*blockIdx.x+tid];
  if (tid ==0)
  {
    rInv = 1.0f/myData.dist;
    sintheta = sqrtf(1.0-rhat[2]*rhat[2]);
    sinthetaInv = 1.0/sintheta;
    cosphi = rhat[0]*sinthetaInv;
    sinphi = rhat[1]*sinthetaInv;
    thetahat[0] = rhat[2]*cosphi;
    thetahat[1] = rhat[2]*sinphi;
    thetahat[2] = -sintheta;
    phihat[0]   = -sinphi;
    phihat[1]   = cosphi;
    phihat[2]   = 0.0f;
  }
  __syncthreads();
  // Compute polynomial basis functions
  // Note maximum polynomial order of 16
  __shared__ T polyfuncs[16][3];
  if (tid < 16)
  {
    polyfuncs[tid][0] = pow(myData.dist,(T)tid);
    polyfuncs[tid][1] = (T)tid    *polyfuncs[tid][0] / myData.dist;
    polyfuncs[tid][2] = (T)(tid-1)*polyfuncs[tid][1] / myData.dist;
  }
  __syncthreads();
  __shared__ T Ylm[(LMAX+1)*(LMAX+1)][2],
             dTheta[(LMAX+1)*(LMAX+1)][2], dPhi[(LMAX+1)*(LMAX+1)][2],
             lpref[(LMAX+1)*(LMAX+1)];
  int numlm = (myOrbital.lMax+1)*(myOrbital.lMax+1);
  int Yblocks = (2*numlm+BS-1)/BS;
  for (int ib=0; ib<Yblocks; ib++)
    if (ib*BS + tid < 2*numlm)
    {
      Ylm[0][ib*BS+tid]    = myYlm[ib*BS + tid];
      dTheta[0][ib*BS+tid] = mydTheta[ib*BS + tid];
      dPhi[0][ib*BS+tid]   = mydPhi[ib*BS + tid];
    }
  for (int l=0; l<=myOrbital.lMax; l++)
    if (tid < 2*l+1)
    {
      int lm = l*l + tid;
      lpref[lm] = -rInv*rInv*(T)(l*(l+1));
    }
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int outIndex=0, outBlock=0;
  __shared__ T outval[BS], outgrad[BS][3], outlap[BS];
  for (int block=0; block<numBlocks; block++)
  {
    T phase_re, phase_im;
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS + tid;
      if (off < 3*N)
        k_red[0][i*BS+tid] = k_reduced[off];
    }
    __shared__ T m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    sincos(-(k_red[tid][0]*myData.img[0]+
             k_red[tid][1]*myData.img[1]+
             k_red[tid][2]*myData.img[2]),
           &phase_im, &phase_re);
    T *c0 =  myCoefs + 2*block*BS + tid;
    T v_re         =T(), v_im         =T();
    T g_rhat_re    =T(), g_rhat_im    =T(),
      g_thetahat_re=T(), g_thetahat_im=T(),
      g_phihat_re  =T(), g_phihat_im  =T(),
      lap_re       =T(), lap_im       =T();
    __shared__ T c[BS][2];
    for (int lm=0; lm<numlm; lm++)
    {
      T u_re=T(), u_im=T(), du_re=T(), du_im=T(),
        d2u_re=T(), d2u_im=T();
      T *coef = c0 + lm*myOrbital.lm_stride;
      for (int n=0; n<=myOrbital.poly_order; n++)
      {
        c[0][tid]    = coef[n*myOrbital.poly_stride];
        c[0][BS+tid] = coef[n*myOrbital.poly_stride+BS];
        __syncthreads();
        u_re   += polyfuncs[n][0]*c[tid][0];
        u_im   += polyfuncs[n][0]*c[tid][1];
        du_re  += polyfuncs[n][1]*c[tid][0];
        du_im  += polyfuncs[n][1]*c[tid][1];
        d2u_re += polyfuncs[n][2]*c[tid][0];
        d2u_im += polyfuncs[n][2]*c[tid][1];
        __syncthreads();
      }
      // Now accumulate values
      v_re    +=   u_re * Ylm[lm][0] - u_im * Ylm[lm][1];
      v_im    +=   u_re * Ylm[lm][1] + u_im * Ylm[lm][0];
      g_rhat_re +=  du_re * Ylm[lm][0] - du_im * Ylm[lm][1];
      g_rhat_im +=  du_re * Ylm[lm][1] + du_im * Ylm[lm][0];
      g_thetahat_re += rInv*(u_re*dTheta[lm][0] - u_im*dTheta[lm][1]);
      g_thetahat_im += rInv*(u_re*dTheta[lm][1] + u_im*dTheta[lm][0]);
      g_phihat_re += rInv*sinthetaInv*
                     (u_re*dPhi[lm][0] - u_im*dPhi[lm][1]);
      g_phihat_im += rInv*sinthetaInv*
                     (u_re*dPhi[lm][1] + u_im*dPhi[lm][0]);
      lap_re += (Ylm[lm][0] * (lpref[lm]*u_re + d2u_re + 2.0f*rInv*du_re) -
                 Ylm[lm][1] * (lpref[lm]*u_im + d2u_im + 2.0f*rInv*du_im));
      lap_im += (Ylm[lm][1] * (lpref[lm]*u_re + d2u_re + 2.0f*rInv*du_re) +
                 Ylm[lm][0] * (lpref[lm]*u_im + d2u_im + 2.0f*rInv*du_im));
    }
    __shared__ T val[BS][2], grad[BS][3][2], lap[BS][2];
    val[tid][0] = phase_re *   v_re - phase_im *   v_im;
    val[tid][1] = phase_re *   v_im + phase_im *   v_re;
    lap[tid][0] = phase_re * lap_re - phase_im * lap_im;
    lap[tid][1] = phase_re * lap_im + phase_im * lap_re;
    for (int i=0; i<3; i++)
    {
      grad[tid][i][0] =
        (phase_re*g_rhat_re     - phase_im*g_rhat_im    )*rhat[i] +
        (phase_re*g_thetahat_re - phase_im*g_thetahat_im)*thetahat[i] +
        (phase_re*g_phihat_re   - phase_im*g_phihat_im  )*phihat[i];
      grad[tid][i][1] =
        (phase_im*g_rhat_re     + phase_re*g_rhat_im    )*rhat[i] +
        (phase_im*g_thetahat_re + phase_re*g_thetahat_im)*thetahat[i] +
        (phase_im*g_phihat_re   + phase_re*g_phihat_im  )*phihat[i];
    }
    __syncthreads();
    // Now serialize to output buffers.
    int iend = min (BS, N-block*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == 0)
      {
        outval[outIndex]     = val[i][0];
        outgrad[outIndex][0] = grad[i][0][0];
        outgrad[outIndex][1] = grad[i][1][0];
        outgrad[outIndex][2] = grad[i][2][0];
        outlap[outIndex]     = lap[i][0];
      }
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        int off = outBlock++*BS+tid;
        myVal[off] = outval[tid];
        myGradLapl[0*row_stride+off] = outgrad[tid][0];
        myGradLapl[1*row_stride+off] = outgrad[tid][1];
        myGradLapl[2*row_stride+off] = outgrad[tid][2];
        myGradLapl[3*row_stride+off] = outlap[tid];
        outIndex = 0;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == 0)
        {
          outval[outIndex]     = val[i][1];
          outgrad[outIndex][0] = grad[i][0][1];
          outgrad[outIndex][1] = grad[i][1][1];
          outgrad[outIndex][2] = grad[i][2][1];
          outlap[outIndex]     = lap[i][1];
        }
        outIndex++;
        __syncthreads();
      }
      if (outIndex == BS)
      {
        int off = outBlock++*BS+tid;
        myVal[off] = outval[tid];
        myGradLapl[0*row_stride+off] = outgrad[tid][0];
        myGradLapl[1*row_stride+off] = outgrad[tid][1];
        myGradLapl[2*row_stride+off] = outgrad[tid][2];
        myGradLapl[3*row_stride+off] = outlap[tid];
        outIndex = 0;
      }
      __syncthreads();
    }
  }
  if (tid < outIndex)
  {
    int off = outBlock*BS+tid;
    myVal[off] = outval[tid];
    myGradLapl[0*row_stride+off] = outgrad[tid][0];
    myGradLapl[1*row_stride+off] = outgrad[tid][1];
    myGradLapl[2*row_stride+off] = outgrad[tid][2];
    myGradLapl[3*row_stride+off] = outlap[tid];
  }
}


template<typename T> void
evaluateHybridSplineComplexToReal
(HybridJobType *job_types, T *rhats,
 T **Ylm, T **dYlm_dTheta, T **dYlm_dphi,
 AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int *make2copies,
 T **vals, T **grad_lapl,
 int row_stride, int N, int numWalkers, int lMax)
{
  const int BS=32;
  dim3 dimGrid(numWalkers);
  dim3 dimBlock(BS);
  if (lMax == 0)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,0><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 1)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,1><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 2)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,2><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 3)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,3><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 4)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,4><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 5)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,5><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 6)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,6><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 7)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,7><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 8)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,8><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
  else if (lMax == 9)
  {
    evaluateHybridSplineComplexToReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
    evaluateHybridPolyComplexToReal_kernel<T,BS,9><<<dimGrid,dimBlock>>>
    (job_types, rhats, Ylm, dYlm_dTheta, dYlm_dphi, orbitals,
     data, k_reduced, make2copies, vals, grad_lapl, row_stride, N);
  }
}


template<typename T, int LMAX, int BS> __global__ void
CalcYlmComplex (T *rhats, HybridJobType  *job_types,
                T **Ylm_ptr, T **dYlm_dtheta_ptr, T **dYlm_dphi_ptr, int N)
{
  const T fourPiInv = 0.0795774715459477f;
  int tid = threadIdx.x;
  const int numlm = (LMAX+1)*(LMAX+1);
  __shared__ T* Ylm[BS], *dtheta[BS], *dphi[BS];
  if (blockIdx.x*BS+tid < N)
  {
    Ylm[tid]    = Ylm_ptr[blockIdx.x*BS+tid];
    dtheta[tid] = dYlm_dtheta_ptr[blockIdx.x*BS+tid];
    dphi[tid]   = dYlm_dphi_ptr[blockIdx.x*BS+tid];
  }
  __shared__ T rhat[BS][3];
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x + i)*BS + tid;
    if (off < 3*N)
      rhat[0][i*BS+tid] = rhats[off];
  }
  __syncthreads();
  T costheta = rhat[tid][2];
  T sintheta = sqrt(1.0f-costheta*costheta);
  T cottheta = costheta/sintheta;
  T cosphi, sinphi;
  cosphi=rhat[tid][0]/sintheta;
  sinphi=rhat[tid][1]/sintheta;
  __shared__ T phi[BS];
  phi[tid] = atan2f(sinphi, cosphi);
  __shared__ T XlmVec[BS][(LMAX+1)*(LMAX+2)/2],
             dXlmVec[BS][(LMAX+1)*(LMAX+2)/2];
  // Create a map from lm linear index to l and m
  __shared__ int l_lm[numlm];
  __shared__ int m_lm[numlm];
  __shared__ T  floatm_lm[numlm];
  int off=0;
  for (int l=0; l<=LMAX; l++)
  {
    if (tid < 2*l+1)
    {
      l_lm[off+tid] = l;
      m_lm[off+tid] = tid-l;
      floatm_lm[off+tid] = (T)(tid-l);
    }
    off += 2*l+1;
  }
  T lsign = 1.0f;
  T dl = 0.0f;
  for (int l=0; l<=LMAX; l++)
  {
    int index=l*(l+3)/2;
    XlmVec[tid][index]  = lsign;
    dXlmVec[tid][index] = dl * cottheta * XlmVec[tid][index];
    T dm = dl;
    for (int m=l; m>0; m--, index--)
    {
      T tmp = sqrt((dl+dm)*(dl-dm+1.0f));
      XlmVec[tid][index-1]  =
        -(dXlmVec[tid][index] + dm*cottheta*XlmVec[tid][index])/ tmp;
      dXlmVec[tid][index-1] =
        (dm-1.0f)*cottheta*XlmVec[tid][index-1] + XlmVec[tid][index]*tmp;
      dm -= 1.0f;
    }
    index = l*(l+1)/2;
    T sum = XlmVec[tid][index] * XlmVec[tid][index];
    for (int m=1; m<=l; m++)
      sum += 2.0f*XlmVec[tid][index+m]*XlmVec[tid][index+m];
    // Now, renormalize the Ylms for this l
    T norm = sqrt((2.0f*dl+1.0f)*fourPiInv / sum);
    for (int m=0; m<=l; m++)
    {
      XlmVec[tid][index+m]  *= norm;
      dXlmVec[tid][index+m] *= norm;
    }
    lsign *= -1.0f;
    dl += 1.0f;
  }
  __syncthreads();
  // Multiply by azimuthal phase and store in Ylm
  int end = min (N-blockIdx.x*BS, BS);
  int nb = ((LMAX+1)*(LMAX+1)+BS-1)/BS;
  __shared__ T outbuff[3][BS][2];
  for (int i=0; i < end; i++)
  {
    // __shared__ T sincosphi[2*LMAX+1][2];
    // if (tid < LMAX
    for (int block=0; block<nb; block++)
    {
      int lm = block*BS + tid;
      if (lm < numlm)
      {
        int l = l_lm[lm];
        int m = m_lm[lm];
        T fm = floatm_lm[lm];
        T re, im;
        sincos(fm*phi[i], &im, &re);
        // Switch sign if m<0 and it's odd
        if (m<0 && (m&1))
        {
          re *= -1.0f;
          im *= -1.0f;
        }
        int off = ((l*(l+1))>>1) + abs(m);
        // Ylm
        outbuff[0][tid][0] =     re *  XlmVec[i][off];
        outbuff[0][tid][1] =     im *  XlmVec[i][off];
        // dYlm_dtheta
        outbuff[1][tid][0] =     re * dXlmVec[i][off];
        outbuff[1][tid][1] =     im * dXlmVec[i][off];
        // dYlm_dphi
        outbuff[2][tid][0] = -fm*im *  XlmVec[i][off];
        outbuff[2][tid][1] =  fm*re *  XlmVec[i][off];
      }
      __syncthreads();
      // Now write back to global mem with coallesced writes
      int off = 2*block*BS + tid;
      if (off < 2*numlm)
      {
        Ylm[i][off]    = outbuff[0][0][tid];
        dtheta[i][off] = outbuff[1][0][tid];
        dphi[i][off]   = outbuff[2][0][tid];
      }
      off += BS;
      if (off < 2*numlm)
      {
        Ylm[i][off]    = outbuff[0][0][tid+BS];
        dtheta[i][off] = outbuff[1][0][tid+BS];
        dphi[i][off]   = outbuff[2][0][tid+BS];
      }
    }
  }
  // complex<T> e2imphi (1.0, 0.0);
  // complex<T> eye(0.0, 1.0);
  // for (int m=0; m<=l; m++) {
  //   Ylm[l*(l+1)+m]  =  XlmVec[tid][l+m]*e2imphi;
  //   Ylm[l*(l+1)-m]  =  XlmVec[tid][l-m]*conj(e2imphi);
  //   dYlm_dphi[l*(l+1)+m ]  =  (double)m * eye *XlmVec[tid][l+m]*e2imphi;
  //   dYlm_dphi[l*(l+1)-m ]  = -(double)m * eye *XlmVec[tid][l-m]*conj(e2imphi);
  //   dYlm_dtheta[l*(l+1)+m] = dXlmVec[tid][l+m]*e2imphi;
  //   dYlm_dtheta[l*(l+1)-m] = dXlmVec[tid][l-m]*conj(e2imphi);
  //   e2imphi *= e2iphi;
  // }
  // dl += 1.0;
  // lsign *= -1.0;
  // YlmTimer.stop();
}


template<typename T, int LMAX, int BS> __global__ void
CalcYlmReal (T *rhats, HybridJobType* job_type,
             T **Ylm_ptr, T **dYlm_dtheta_ptr, T **dYlm_dphi_ptr, int N)
{
  const T fourPiInv = 0.0795774715459477f;
  int tid = threadIdx.x;
  const int numlm = (LMAX+1)*(LMAX+2)/2;
  __shared__ T* Ylm[BS], *dtheta[BS], *dphi[BS];
  if (blockIdx.x*BS+tid < N)
  {
    Ylm[tid]    = Ylm_ptr[blockIdx.x*BS+tid];
    dtheta[tid] = dYlm_dtheta_ptr[blockIdx.x*BS+tid];
    dphi[tid]   = dYlm_dphi_ptr[blockIdx.x*BS+tid];
  }
  __shared__ T rhat[BS][3];
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x + i)*BS + tid;
    if (off < 3*N)
      rhat[0][i*BS+tid] = rhats[off];
  }
  __syncthreads();
  T costheta = rhat[tid][2];
  T sintheta = sqrt(1.0f-costheta*costheta);
  T cottheta = costheta/sintheta;
  T cosphi, sinphi;
  cosphi=rhat[tid][0]/sintheta;
  sinphi=rhat[tid][1]/sintheta;
  __shared__ T phi[BS];
  phi[tid] = atan2f(sinphi, cosphi);
  __shared__ T XlmVec[BS][numlm],
             dXlmVec[BS][numlm];
  // Create a map from lm linear index to l and m
  __shared__ int l_lm[numlm];
  __shared__ int m_lm[numlm];
  __shared__ T  floatm_lm[numlm];
  int off=0;
  for (int l=0; l<=LMAX; l++)
  {
    if (tid < l+1)
    {
      l_lm[off+tid] = l;
      m_lm[off+tid] = tid;
      floatm_lm[off+tid] = (T)tid;
    }
    off += l+1;
  }
  T lsign = 1.0f;
  T dl = 0.0f;
  for (int l=0; l<=LMAX; l++)
  {
    int index=l*(l+3)/2;
    XlmVec[tid][index]  = lsign;
    dXlmVec[tid][index] = dl * cottheta * XlmVec[tid][index];
    T dm = dl;
    for (int m=l; m>0; m--, index--)
    {
      T tmp = sqrt((dl+dm)*(dl-dm+1.0f));
      XlmVec[tid][index-1]  =
        -(dXlmVec[tid][index] + dm*cottheta*XlmVec[tid][index])/ tmp;
      dXlmVec[tid][index-1] =
        (dm-1.0f)*cottheta*XlmVec[tid][index-1] + XlmVec[tid][index]*tmp;
      dm -= 1.0f;
    }
    index = l*(l+1)/2;
    T sum = XlmVec[tid][index] * XlmVec[tid][index];
    for (int m=1; m<=l; m++)
      sum += 2.0f*XlmVec[tid][index+m]*XlmVec[tid][index+m];
    // Now, renormalize the Ylms for this l
    T norm = sqrt((2.0f*dl+1.0f)*fourPiInv / sum);
    for (int m=0; m<=l; m++)
    {
      XlmVec[tid][index+m]  *= norm;
      dXlmVec[tid][index+m] *= norm;
    }
    lsign *= -1.0f;
    dl += 1.0f;
  }
  __syncthreads();
  // Multiply by azimuthal phase and store in Ylm
  int end = min (N-blockIdx.x*BS, BS);
  int nb = (numlm+BS-1)/BS;
  __shared__ T outbuff[3][2*BS];
  for (int i=0; i < end; i++)
  {
    for (int block=0; block<nb; block++)
    {
      int lm = block*BS + tid;
      if (lm < numlm)
      {
        int l = l_lm[lm];
        int m = m_lm[lm];
        T fm = floatm_lm[lm];
        T re, im;
        sincos(fm*phi[i], &im, &re);
        int off = ((l*(l+1))>>1) + m;
        int iplus = l*(l+1)+m;
        int iminus = l*(l+1)-m;
        // Ylm
        outbuff[0][iplus] =     re *  XlmVec[i][off];
        // dYlm_dtheta
        outbuff[1][iplus] =     re * dXlmVec[i][off];
        // dYlm_dphi
        outbuff[2][iplus] = -fm*im *  XlmVec[i][off];
        if (m != 0)
        {
          outbuff[0][iminus] =     im *  XlmVec[i][off];
          outbuff[1][iminus] =     im * dXlmVec[i][off];
          outbuff[2][iminus] =  fm*re *  XlmVec[i][off];
        }
      }
      __syncthreads();
      // Now write back to global mem with coallesced writes
      int off = block*BS + tid;
      if (off < (LMAX+1)*(LMAX+1))
      {
        Ylm[i][off]    = outbuff[0][tid];
        dtheta[i][off] = outbuff[1][tid];
        dphi[i][off]   = outbuff[2][tid];
      }
      off += BS;
      if (off < (LMAX+1)*(LMAX+1))
      {
        Ylm[i][off]    = outbuff[0][tid+BS];
        dtheta[i][off] = outbuff[1][tid+BS];
        dphi[i][off]   = outbuff[2][tid+BS];
      }
    }
  }
}


template<typename T> void
CalcYlmRealCuda (T *rhats, HybridJobType *job_type,
                 T **Ylm_ptr, T **dYlm_dtheta_ptr, T **dYlm_dphi_ptr,
                 int lMax, int N)
{
  const int BS=32;
  int Nblocks = (N+BS-1)/BS;
  dim3 dimGrid(Nblocks);
  dim3 dimBlock(BS);
  if (lMax == 0)
    return;
  else if (lMax == 1)
    CalcYlmReal<T,1,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 2)
    CalcYlmReal<T,2,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 3)
    CalcYlmReal<T,3,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 4)
    CalcYlmReal<T,4,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 5)
    CalcYlmReal<T,5,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 6)
    CalcYlmReal<T,6,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 7)
    CalcYlmReal<T,7,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 8)
    CalcYlmReal<T,8,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
}


template<typename T> void
CalcYlmComplexCuda (T *rhats, HybridJobType *job_type,
                    T **Ylm_ptr, T **dYlm_dtheta_ptr, T **dYlm_dphi_ptr,
                    int lMax, int N)
{
  const int BS=32;
  int Nblocks = (N+BS-1)/BS;
  dim3 dimGrid(Nblocks);
  dim3 dimBlock(BS);
  if (lMax == 0)
    return;
  else if (lMax == 1)
    CalcYlmComplex<T,1,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 2)
    CalcYlmComplex<T,2,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 3)
    CalcYlmComplex<T,3,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 4)
    CalcYlmComplex<T,4,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 5)
    CalcYlmComplex<T,5,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 6)
    CalcYlmComplex<T,6,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 7)
    CalcYlmComplex<T,7,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
  else if (lMax == 8)
    CalcYlmComplex<T,8,BS><<<dimGrid,dimBlock>>>
    (rhats,job_type,Ylm_ptr,dYlm_dtheta_ptr,dYlm_dphi_ptr,N);
}


template<typename T, int LMAX, int BS> __global__ void
CalcYlmComplex (T *rhats, HybridJobType *job_types, T **Ylm_ptr, int N)
{
  const T fourPiInv = 0.0795774715459477f;
  int tid = threadIdx.x;
  const int numlm = (LMAX+1)*(LMAX+1);
  __shared__ T* Ylm[BS];
  if (blockIdx.x*BS+tid < N)
    Ylm[tid]    = Ylm_ptr[blockIdx.x*BS+tid];
  __shared__ T rhat[BS][3];
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x + i)*BS + tid;
    if (off < 3*N)
      rhat[0][i*BS+tid] = rhats[off];
  }
  __syncthreads();
  T costheta = rhat[tid][2];
  T sintheta = sqrt(1.0f-costheta*costheta);
  T cottheta = costheta/sintheta;
  T cosphi, sinphi;
  cosphi=rhat[tid][0]/sintheta;
  sinphi=rhat[tid][1]/sintheta;
  __shared__ T phi[BS];
  phi[tid] = atan2f(sinphi, cosphi);
  __shared__ T XlmVec[BS][(LMAX+1)*(LMAX+2)/2],
             dXlmVec[BS][(LMAX+1)*(LMAX+2)/2];
  // Create a map from lm linear index to l and m
  __shared__ int l_lm[numlm];
  __shared__ int m_lm[numlm];
  __shared__ T  floatm_lm[numlm];
  int off=0;
  for (int l=0; l<=LMAX; l++)
  {
    if (tid < 2*l+1)
    {
      l_lm[off+tid] = l;
      m_lm[off+tid] = tid-l;
      floatm_lm[off+tid] = (T)(tid-l);
    }
    off += 2*l+1;
  }
  T lsign = 1.0f;
  T dl = 0.0f;
  for (int l=0; l<=LMAX; l++)
  {
    int index=l*(l+3)/2;
    XlmVec[tid][index]  = lsign;
    dXlmVec[tid][index] = dl * cottheta * XlmVec[tid][index];
    T dm = dl;
    for (int m=l; m>0; m--, index--)
    {
      T tmp = sqrt((dl+dm)*(dl-dm+1.0f));
      XlmVec[tid][index-1]  =
        -(dXlmVec[tid][index] + dm*cottheta*XlmVec[tid][index])/ tmp;
      dXlmVec[tid][index-1] =
        (dm-1.0f)*cottheta*XlmVec[tid][index-1] + XlmVec[tid][index]*tmp;
      dm -= 1.0f;
    }
    index = l*(l+1)/2;
    T sum = XlmVec[tid][index] * XlmVec[tid][index];
    for (int m=1; m<=l; m++)
      sum += 2.0f*XlmVec[tid][index+m]*XlmVec[tid][index+m];
    // Now, renormalize the Ylms for this l
    T norm = sqrt((2.0f*dl+1.0f)*fourPiInv / sum);
    for (int m=0; m<=l; m++)
    {
      XlmVec[tid][index+m]  *= norm;
      dXlmVec[tid][index+m] *= norm;
    }
    lsign *= -1.0f;
    dl += 1.0f;
  }
  __syncthreads();
  // Multiply by azimuthal phase and store in Ylm
  int end = min (N-blockIdx.x*BS, BS);
  int nb = ((LMAX+1)*(LMAX+1)+BS-1)/BS;
  __shared__ T outbuff[BS][2];
  for (int i=0; i < end; i++)
  {
    // __shared__ T sincosphi[2*LMAX+1][2];
    // if (tid < LMAX
    for (int block=0; block<nb; block++)
    {
      int lm = block*BS + tid;
      if (lm < numlm)
      {
        int l = l_lm[lm];
        int m = m_lm[lm];
        T fm = floatm_lm[lm];
        T re, im;
        sincos(fm*phi[i], &im, &re);
        // Switch sign if m<0 and it's odd
        if (m<0 && (m&1))
        {
          re *= -1.0f;
          im *= -1.0f;
        }
        int off = ((l*(l+1))>>1) + abs(m);
        // Ylm
        outbuff[tid][0] =     re *  XlmVec[i][off];
        outbuff[tid][1] =     im *  XlmVec[i][off];
      }
      __syncthreads();
      // Now write back to global mem with coallesced writes
      int off = 2*block*BS + tid;
      if (off < 2*numlm)
        Ylm[i][off]    = outbuff[0][tid];
      off += BS;
      if (off < 2*numlm)
        Ylm[i][off]    = outbuff[0][tid+BS];
    }
  }
}


template<typename T, int LMAX, int BS> __global__ void
CalcYlmReal (T *rhats, HybridJobType *job_types, T **Ylm_ptr, int N)
{
  const T fourPiInv = 0.0795774715459477f;
  int tid = threadIdx.x;
  const int numlm = (LMAX+1)*(LMAX+2)/2;
  __shared__ T* Ylm[BS];
  if (blockIdx.x*BS+tid < N)
    Ylm[tid]    = Ylm_ptr[blockIdx.x*BS+tid];
  __shared__ T rhat[BS][3];
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x + i)*BS + tid;
    if (off < 3*N)
      rhat[0][i*BS+tid] = rhats[off];
  }
  __syncthreads();
  T costheta = rhat[tid][2];
  T sintheta = sqrt(1.0f-costheta*costheta);
  T cottheta = costheta/sintheta;
  T cosphi, sinphi;
  cosphi=rhat[tid][0]/sintheta;
  sinphi=rhat[tid][1]/sintheta;
  __shared__ T phi[BS];
  phi[tid] = atan2f(sinphi, cosphi);
  __shared__ T XlmVec[BS][numlm],
             dXlmVec[BS][numlm];
  // Create a map from lm linear index to l and m
  __shared__ int l_lm[numlm];
  __shared__ int m_lm[numlm];
  __shared__ T  floatm_lm[numlm];
  int off=0;
  for (int l=0; l<=LMAX; l++)
  {
    if (tid < l+1)
    {
      l_lm[off+tid] = l;
      m_lm[off+tid] = tid;
      floatm_lm[off+tid] = (T)tid;
    }
    off += l+1;
  }
  T lsign = 1.0f;
  T dl = 0.0f;
  for (int l=0; l<=LMAX; l++)
  {
    int index=l*(l+3)/2;
    XlmVec[tid][index]  = lsign;
    dXlmVec[tid][index] = dl * cottheta * XlmVec[tid][index];
    T dm = dl;
    for (int m=l; m>0; m--, index--)
    {
      T tmp = sqrt((dl+dm)*(dl-dm+1.0f));
      XlmVec[tid][index-1]  =
        -(dXlmVec[tid][index] + dm*cottheta*XlmVec[tid][index])/ tmp;
      dXlmVec[tid][index-1] =
        (dm-1.0f)*cottheta*XlmVec[tid][index-1] + XlmVec[tid][index]*tmp;
      dm -= 1.0f;
    }
    index = l*(l+1)/2;
    T sum = XlmVec[tid][index] * XlmVec[tid][index];
    for (int m=1; m<=l; m++)
      sum += 2.0f*XlmVec[tid][index+m]*XlmVec[tid][index+m];
    // Now, renormalize the Ylms for this l
    T norm = sqrt((2.0f*dl+1.0f)*fourPiInv / sum);
    for (int m=0; m<=l; m++)
    {
      XlmVec[tid][index+m]  *= norm;
      dXlmVec[tid][index+m] *= norm;
    }
    lsign *= -1.0f;
    dl += 1.0f;
  }
  __syncthreads();
  // Multiply by azimuthal phase and store in Ylm
  int end = min (N-blockIdx.x*BS, BS);
  int nb = (numlm+BS-1)/BS;
  __shared__ T outbuff[2*BS];
  for (int i=0; i < end; i++)
  {
    for (int block=0; block<nb; block++)
    {
      int lm = block*BS + tid;
      if (lm < numlm)
      {
        int l = l_lm[lm];
        int m = m_lm[lm];
        T fm = floatm_lm[lm];
        T re, im;
        sincos(fm*phi[i], &im, &re);
        int off = ((l*(l+1))>>1) + m;
        int iplus = l*(l+1)+m;
        int iminus = l*(l+1)-m;
        // Ylm
        outbuff[iplus] =     re *  XlmVec[i][off];
        if (m != 0)
          outbuff[iminus] =     im *  XlmVec[i][off];
      }
      __syncthreads();
      // Now write back to global mem with coallesced writes
      int off = block*BS + tid;
      if (off < (LMAX+1)*(LMAX+1))
        Ylm[i][off]    = outbuff[tid];
      off += BS;
      if (off < (LMAX+1)*(LMAX+1))
        Ylm[i][off]    = outbuff[tid+BS];
    }
  }
}

// YingWai: these seems to be unused  (Oct 1, 15)
/*
template<typename T> void
CalcYlmRealCuda (T *rhats, HybridJobType *job_types,
                 T **Ylm_ptr, int lMax, int N)
{
  const int BS=32;
  int Nblocks = (N+BS-1)/BS;
  dim3 dimGrid(Nblocks);
  dim3 dimBlock(BS);
  if (lMax == 0)
    return;
  else if (lMax == 1)
    CalcYlmReal<T,1,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 2)
    CalcYlmReal<T,2,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 3)
    CalcYlmReal<T,3,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 4)
    CalcYlmReal<T,4,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 5)
    CalcYlmReal<T,5,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 6)
    CalcYlmReal<T,6,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 7)
    CalcYlmReal<T,7,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
  else if (lMax == 8)
    CalcYlmReal<T,8,BS><<<dimGrid,dimBlock>>>(rhats,job_types,Ylm_ptr,N);
}

template<typename T> void
CalcYlmComplexCuda (T *rhats, HybridJobType *job_type,
                    T **Ylm_ptr, int lMax, int N)
{
  const int BS=32;
  int Nblocks = (N+BS-1)/BS;
  dim3 dimGrid(Nblocks);
  dim3 dimBlock(BS);
  if (lMax == 0)
    return;
  else if (lMax == 1)
    CalcYlmComplex<T,1,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 2)
    CalcYlmComplex<T,2,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 3)
    CalcYlmComplex<T,3,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 4)
    CalcYlmComplex<T,4,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 5)
    CalcYlmComplex<T,5,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 6)
    CalcYlmComplex<T,6,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 7)
    CalcYlmComplex<T,7,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
  else if (lMax == 8)
    CalcYlmComplex<T,8,BS><<<dimGrid,dimBlock>>>(rhats,job_type,Ylm_ptr,N);
}

void dummy_float()
{
  float *rhats(0), **Ylm_ptr(0), **dYlm_dtheta_ptr(0), **dYlm_dphi_ptr(0);
  HybridJobType* job_types(0);
  CalcYlmRealCuda(rhats,    job_types, Ylm_ptr, dYlm_dtheta_ptr, dYlm_dphi_ptr, 1, 1);
  CalcYlmComplexCuda(rhats, job_types, Ylm_ptr, dYlm_dtheta_ptr, dYlm_dphi_ptr, 1, 1);
  CalcYlmRealCuda(rhats,    job_types, Ylm_ptr, 1, 1);
  CalcYlmComplexCuda(rhats, job_types, Ylm_ptr, 1, 1);
}
*/

template<typename T, int BS> __global__ void
evaluate3DSplineReal_kernel (HybridJobType *job_types, T *pos, T *k_reduced,
                             typename cudaTypeTraits<T>::realType3 drInv, T *coefs, 
                             uint3 dim, uint3 strides,
                             T *Linv, T **vals, T **grad_lapl,
                             int row_stride, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.y];
  __syncthreads();
  if (myjob != BSPLINE_3D_JOB)
    return;
  int ir    = blockIdx.y;
  int off   = blockIdx.x*BS+threadIdx.x;
  __shared__ T *myval, *mygrad_lapl;
  __shared__ T r[3], u[3], img[3];
  __shared__ T k_red[BS][3];
  __shared__ T G[3][3], GGt[3][3];
  int i0 = tid/3;
  int i1 = tid - 3*i0;
  if (tid < 9)
  {
    G[0][tid] = Linv[tid];
    GGt[i0][i1] = (G[0][i0]*G[0][i1] +
                   G[1][i0]*G[1][i1] +
                   G[2][i0]*G[2][i1]);
  }
  if (tid == 0)
  {
    myval  = vals[ir];
    mygrad_lapl = grad_lapl[ir];
  }
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x+i)*BS+tid;
    if (off < 3*N)
      k_red[0][i*BS+tid] = k_reduced[off];
  }
  if (tid < 3)
  {
    r[tid] = pos[3*ir+tid];
    u[tid] = G[0][tid]*r[0] + G[1][tid]*r[1] + G[2][tid]*r[2];
    img[tid] = floor(u[tid]);
    u[tid] -= img[tid];
  }
  __syncthreads();
  T sign = cos(-(k_red[tid][0]*img[0]+
                 k_red[tid][1]*img[1]+
                 k_red[tid][2]*img[2]));
  // copysign(1.0f,sin(-(r[0]*kp[tid][0] +
  // 		   r[1]*kp[tid][1] +
  // 		   r[2]*kp[tid][2])));
  __syncthreads();
  int3 index;
  typename cudaTypeTraits<T>::realType3 t;
  T s, sf;
  typename cudaTypeTraits<T>::realType4 tp[3];
  s = u[0] * drInv.x;
  sf = floor(s);
  index.x = min(max(0,(int)sf), dim.x-1);
  t.x = s - sf;
  s = u[1] * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);
  t.y = s - sf;
  s = u[2] * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
  t.z = s - sf;
  tp[0] = cudaMakeType4<T>(t.x*t.x*t.x, t.x*t.x, t.x, 1.0);
  tp[1] = cudaMakeType4<T>(t.y*t.y*t.y, t.y*t.y, t.y, 1.0);
  tp[2] = cudaMakeType4<T>(t.z*t.z*t.z, t.z*t.z, t.z, 1.0);
  // First 4 of a are value, second 4 are derivative, last four are
  // second derivative.
  __shared__ T a[12], b[12], c[12];
  if (tid < 12)
  {
    a[tid] = Acuda[4*tid+0]*tp[0].x + Acuda[4*tid+1]*tp[0].y + Acuda[4*tid+2]*tp[0].z + Acuda[4*tid+3]*tp[0].w;
    b[tid] = Acuda[4*tid+0]*tp[1].x + Acuda[4*tid+1]*tp[1].y + Acuda[4*tid+2]*tp[1].z + Acuda[4*tid+3]*tp[1].w;
    c[tid] = Acuda[4*tid+0]*tp[2].x + Acuda[4*tid+1]*tp[2].y + Acuda[4*tid+2]*tp[2].z + Acuda[4*tid+3]*tp[2].w;
  }
  __syncthreads();
  __shared__ T abc[640];
  int i = (tid>>4)&3;
  int j = (tid>>2)&3;
  int k = (tid & 3);
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
  T v = 0.0, g0=0.0,  g1=0.0, g2=0.0,
        h00=0.0, h01=0.0, h02=0.0, h11=0.0, h12=0.0, h22=0.0;
  int n = 0;
  T *b0 = coefs + index.x*strides.x + index.y*strides.y + index.z*strides.z + off;
  if (off < N)
  {
    for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
      {
        T *base = b0 + i*strides.x + j*strides.y;
        for (int k=0; k<4; k++)
        {
          T c  = base[k*strides.z];
          v   += abc[n+  0] * c;
          g0  += abc[n+ 64] * c;
          g1  += abc[n+128] * c;
          g2  += abc[n+192] * c;
          h00 += abc[n+256] * c;
          h01 += abc[n+320] * c;
          h02 += abc[n+384] * c;
          h11 += abc[n+448] * c;
          h12 += abc[n+512] * c;
          h22 += abc[n+576] * c;
          n += 1;
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
    myval[off] = sign*v;
  }
  if (off < N)
  {
    // Store gradients back to global memory
    mygrad_lapl[off+0*row_stride] = sign*(G[0][0]*g0 + G[0][1]*g1 + G[0][2]*g2);
    mygrad_lapl[off+1*row_stride] = sign*(G[1][0]*g0 + G[1][1]*g1 + G[1][2]*g2);
    mygrad_lapl[off+2*row_stride] = sign*(G[2][0]*g0 + G[2][1]*g1 + G[2][2]*g2);
    // Store laplacians back to global memory
    // Hessian = H00 H01 H02 H11 H12 H22
    // Matrix = [0 1 2]
    //          [1 3 4]
    //          [2 4 5]
    // laplacian = Trace(GGt*Hessian)
    mygrad_lapl[off+3*row_stride] = sign *
                                    (GGt[0][0]*h00 + GGt[1][0]*h01 + GGt[2][0]*h02 +
                                     GGt[0][1]*h01 + GGt[1][1]*h11 + GGt[2][1]*h12 +
                                     GGt[0][2]*h02 + GGt[1][2]*h12 + GGt[2][2]*h22);
  }
}


template<typename T> void
evaluate3DSplineReal (HybridJobType *job_types, T *pos, T *kpoints_reduced,
                      typename SplineTraits<T,3>::CudaSplineType *multispline, T *Linv,
                      T **vals, T **grad_lapl,
                      int row_stride, int N, int numWalkers)
{
  const int BS=64;
  dim3 dimGrid((N+BS-1)/BS,numWalkers);
  dim3 dimBlock(BS);
  evaluate3DSplineReal_kernel<T,BS><<<dimGrid,dimBlock>>>
  (job_types, pos, kpoints_reduced, multispline->gridInv, multispline->coefs,
   multispline->dim, multispline->stride, Linv, vals, grad_lapl,
   row_stride, N);
}


template<typename T, int BS> __global__ void
evaluate3DSplineReal_kernel (HybridJobType *job_types, T *pos, T *kpoints_reduced,
                             typename cudaTypeTraits<T>::realType3 drInv, T *coefs,
                             uint3 dim, uint3 strides,
                             T *Linv, T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.y];
  __syncthreads();
  if (myjob != BSPLINE_3D_JOB)
    return;
  int ir    = blockIdx.y;
  int off   = blockIdx.x*BS+threadIdx.x;
  __shared__ T *myval;
  __shared__ T r[3], u[3], img[3];
  __shared__ T k_red[BS][3];
  __shared__ T G[3][3];
  T sign;
  if (tid < 9)
    G[0][tid] = Linv[tid];
  if (tid == 0)
    myval  = vals[ir];
  for (int i=0; i<3; i++)
  {
    int off = (3*blockIdx.x+i)*BS+tid;
    if (off < 3*N)
      k_red[0][i*BS+tid] = kpoints_reduced[off];
  }
  if (tid < 3)
  {
    r[tid] = pos[3*ir+tid];
    u[tid] = G[0][tid]*r[0] + G[1][tid]*r[1] + G[2][tid]*r[2];
    img[tid] = floor(u[tid]);
    u[tid] -= img[tid];
  }
  __syncthreads();
  sign = cos(-(k_red[tid][0]*img[0]+
               k_red[tid][1]*img[1]+
               k_red[tid][2]*img[2]));
  // sign = copysign(1.0f,sin(-(r[0]*kp[tid][0] +
  //                            r[1]*kp[tid][1] +
  //                            r[2]*kp[tid][2])));
  // for (int i=0; i<3; i++) {
  //   int off = (3*blockIdx.x+i)*BS+tid;
  //   if (off < 3*N)
  //     kp[0][i*BS+tid] = kpoints[off];
  // }
  __syncthreads();
  int3 index;
  typename cudaTypeTraits<T>::realType3 t;
  T s, sf;
  typename cudaTypeTraits<T>::realType4 tp[3];
  s = u[0] * drInv.x;
  sf = floor(s);
  index.x = min(max(0,(int)sf), dim.x-1);
  t.x = s - sf;
  s = u[1] * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);
  t.y = s - sf;
  s = u[2] * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
  t.z = s - sf;
  tp[0] = cudaMakeType4<T>(t.x*t.x*t.x, t.x*t.x, t.x, 1.0);
  tp[1] = cudaMakeType4<T>(t.y*t.y*t.y, t.y*t.y, t.y, 1.0);
  tp[2] = cudaMakeType4<T>(t.z*t.z*t.z, t.z*t.z, t.z, 1.0);
  // First 4 of a are value, second 4 are derivative, last four are
  // second derivative.
  __shared__ T a[4], b[4], c[4];
  if (tid < 4)
  {
    a[tid] = Acuda[4*tid+0]*tp[0].x + Acuda[4*tid+1]*tp[0].y + Acuda[4*tid+2]*tp[0].z + Acuda[4*tid+3]*tp[0].w;
    b[tid] = Acuda[4*tid+0]*tp[1].x + Acuda[4*tid+1]*tp[1].y + Acuda[4*tid+2]*tp[1].z + Acuda[4*tid+3]*tp[1].w;
    c[tid] = Acuda[4*tid+0]*tp[2].x + Acuda[4*tid+1]*tp[2].y + Acuda[4*tid+2]*tp[2].z + Acuda[4*tid+3]*tp[2].w;
  }
  __syncthreads();
  __shared__ T abc[64];
  int i = (tid>>4)&3;
  int j = (tid>>2)&3;
  int k = (tid & 3);
  if (tid < 64)
    abc[tid] = a[i]*b[j]*c[k];
  __syncthreads();
  if (off < N)
  {
    T val = 0.0;
    for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
      {
        T *base = coefs + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
        for (int k=0; k<4; k++)
          val += abc[16*i+4*j+k] * base[off+k*strides.z];
      }
    }
    myval[off] = sign * val;
  }
}


template<typename T> void
evaluate3DSplineReal (HybridJobType *job_types, T *pos, T *kpoints,
                      typename SplineTraits<T,3>::CudaSplineType *multispline, T *Linv,
                      T **vals, int N, int numWalkers)
{
  const int BS=64;
  dim3 dimGrid((N+BS-1)/BS,numWalkers);
  dim3 dimBlock(BS);
  evaluate3DSplineReal_kernel<T,BS><<<dimGrid,dimBlock>>>
  (job_types, pos, kpoints, multispline->gridInv, multispline->coefs,
   multispline->dim, multispline->stride, Linv, vals, N);
}


/***************************************************
// 3D B-spline complex-to-real evaluation functions //
***************************************************/

template<typename T, int BS> __global__ void
evaluate3DSplineComplexToReal_kernel
(HybridJobType *job_types, T *pos, T *kpoints, int *make2copies,
 typename cudaTypeTraits<T>::realType3 drInv, T *coefs, uint3 dim, uint3 strides,
 T *Linv, T **vals, T **grad_lapl,
 int row_stride, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != BSPLINE_3D_JOB)
    return;
  int ir    = blockIdx.x;
  __shared__ T *myVal, *myGradLapl;
  __shared__ T r[3], u[3];
  __shared__ T G[3][3], GGt[3][3];
  int i0 = tid/3;
  int i1 = tid - 3*i0;
  if (tid < 9)
  {
    G[0][tid] = Linv[tid];
    GGt[i0][i1] = (G[0][i0]*G[0][i1] +
                   G[1][i0]*G[1][i1] +
                   G[2][i0]*G[2][i1]);
  }
  if (tid == 0)
  {
    myVal      =      vals[ir];
    myGradLapl = grad_lapl[ir];
  }
  if (tid < 3)
  {
    r[tid] = pos[3*ir+tid];
    u[tid] = G[0][tid]*r[0] + G[1][tid]*r[1] + G[2][tid]*r[2];
    u[tid] -= floor(u[tid]);
  }
  __syncthreads();
  int3 index;
  typename cudaTypeTraits<T>::realType3 t;
  T s, sf;
  typename cudaTypeTraits<T>::realType4 tp[3];
  s = u[0] * drInv.x;
  sf = floor(s);
  index.x = min(max(0,(int)sf), dim.x-1);
  t.x = s - sf;
  s = u[1] * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);
  t.y = s - sf;
  s = u[2] * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
  t.z = s - sf;
  tp[0] = cudaMakeType4<T>(t.x*t.x*t.x, t.x*t.x, t.x, 1.0f);
  tp[1] = cudaMakeType4<T>(t.y*t.y*t.y, t.y*t.y, t.y, 1.0f);
  tp[2] = cudaMakeType4<T>(t.z*t.z*t.z, t.z*t.z, t.z, 1.0f);
  // First 4 of a are value, second 4 are derivative, last four are
  // second derivative.
  __shared__ T a[12], b[12], c[12];
  if (tid < 12)
  {
    a[tid] = Acuda[4*tid+0]*tp[0].x + Acuda[4*tid+1]*tp[0].y + Acuda[4*tid+2]*tp[0].z + Acuda[4*tid+3]*tp[0].w;
    b[tid] = Acuda[4*tid+0]*tp[1].x + Acuda[4*tid+1]*tp[1].y + Acuda[4*tid+2]*tp[1].z + Acuda[4*tid+3]*tp[1].w;
    c[tid] = Acuda[4*tid+0]*tp[2].x + Acuda[4*tid+1]*tp[2].y + Acuda[4*tid+2]*tp[2].z + Acuda[4*tid+3]*tp[2].w;
  }
  __syncthreads();
  __shared__ T abc[640];
  int i = (tid>>4)&3;
  int j = (tid>>2)&3;
  int k = (tid & 3);
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
  int numBlocks = (N+BS-1)/BS;
  int outIndex=0, outBlock=0;
  __shared__ T outval[BS], outgrad[BS][3], outlap[BS];
  //////////////////////
  // Outer block loop //
  //////////////////////
  for (int block=0; block<numBlocks; block++)
  {
    __shared__ T kp[BS][3];
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS+tid;
      if (off < 3*N)
        kp[0][i*BS+tid] = kpoints[off];
    }
    __shared__ int m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    T phase_re, phase_im;
    sincos(-(r[0]*kp[tid][0] +
             r[1]*kp[tid][1] +
             r[2]*kp[tid][2]), &phase_im, &phase_re);;
    T v=0.0f, g0=0.0f, g1=0.0f, g2=0.0f,
          h00=0.0f, h01=0.0f, h02=0.0f, h11=0.0f, h12=0.0f, h22=0.0f;
    __shared__ T val[BS][2], grad[BS][2][3], lapl[BS][2];
    int off   = 2*block*BS+threadIdx.x;
    int n = 0;
    T *b0 = coefs + index.x*strides.x + index.y*strides.y +
                index.z*strides.z + off;
    if (off < 2*N)
    {
      for (int i=0; i<4; i++)
      {
        for (int j=0; j<4; j++)
        {
          T *base = b0 + i*strides.x + j*strides.y;
          for (int k=0; k<4; k++)
          {
            T c  = base[k*strides.z];
            v   += abc[n+  0] * c;
            g0  += abc[n+ 64] * c;
            g1  += abc[n+128] * c;
            g2  += abc[n+192] * c;
            h00 += abc[n+256] * c;
            h01 += abc[n+320] * c;
            h02 += abc[n+384] * c;
            h11 += abc[n+448] * c;
            h12 += abc[n+512] * c;
            h22 += abc[n+576] * c;
            n += 1;
          }
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
    val[0][tid]    = v;
    grad[0][tid][0] = G[0][0]*g0 + G[0][1]*g1 + G[0][2]*g2;
    grad[0][tid][1] = G[1][0]*g0 + G[1][1]*g1 + G[1][2]*g2;
    grad[0][tid][2] = G[2][0]*g0 + G[2][1]*g1 + G[2][2]*g2;
    lapl[0][tid]    = (GGt[0][0]*h00 + GGt[1][0]*h01 + GGt[2][0]*h02 +
                       GGt[0][1]*h01 + GGt[1][1]*h11 + GGt[2][1]*h12 +
                       GGt[0][2]*h02 + GGt[1][2]*h12 + GGt[2][2]*h22);
    v=g0=g1=g2=h00=h01=h02=h11=h12=h22=0.0f;
    off   = (2*block+1)*BS+threadIdx.x;
    n = 0;
    b0 = coefs + index.x*strides.x + index.y*strides.y + index.z*strides.z + off;
    if (off < 2*N)
    {
      for (int i=0; i<4; i++)
      {
        for (int j=0; j<4; j++)
        {
          T *base = b0 + i*strides.x + j*strides.y;
          for (int k=0; k<4; k++)
          {
            T c  = base[k*strides.z];
            v   += abc[n+  0] * c;
            g0  += abc[n+ 64] * c;
            g1  += abc[n+128] * c;
            g2  += abc[n+192] * c;
            h00 += abc[n+256] * c;
            h01 += abc[n+320] * c;
            h02 += abc[n+384] * c;
            h11 += abc[n+448] * c;
            h12 += abc[n+512] * c;
            h22 += abc[n+576] * c;
            n += 1;
          }
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
    val[0][tid+BS]    = v;
    grad[0][tid+BS][0] = G[0][0]*g0 + G[0][1]*g1 + G[0][2]*g2;
    grad[0][tid+BS][1] = G[1][0]*g0 + G[1][1]*g1 + G[1][2]*g2;
    grad[0][tid+BS][2] = G[2][0]*g0 + G[2][1]*g1 + G[2][2]*g2;
    lapl[0][tid+BS]    = (GGt[0][0]*h00 + GGt[1][0]*h01 + GGt[2][0]*h02 +
                          GGt[0][1]*h01 + GGt[1][1]*h11 + GGt[2][1]*h12 +
                          GGt[0][2]*h02 + GGt[1][2]*h12 + GGt[2][2]*h22);
    __syncthreads();
    T re, im;
    // Add phase contribution to laplacian
    T k2 = (kp[tid][0]*kp[tid][0] +
                kp[tid][1]*kp[tid][1] +
                kp[tid][2]*kp[tid][2]);
    re = lapl[tid][0] - k2*val[tid][0] +
         2.0f*(kp[tid][0]*grad[tid][1][0]+
               kp[tid][1]*grad[tid][1][1]+
               kp[tid][2]*grad[tid][1][2]);
    im = lapl[tid][1] - k2*val[tid][1] -
         2.0f*(kp[tid][0]*grad[tid][0][0]+
               kp[tid][1]*grad[tid][0][1]+
               kp[tid][2]*grad[tid][0][2]);
    lapl[tid][0] = re*phase_re - im*phase_im;
    lapl[tid][1] = re*phase_im + im*phase_re;
    // Gradient = e^(-ikr)*(-i*u*k + gradu)
    for (int dim=0; dim<3; dim++)
    {
      re = grad[tid][0][dim] + kp[tid][dim]*val[tid][1];
      im = grad[tid][1][dim] - kp[tid][dim]*val[tid][0];
      grad[tid][0][dim] = phase_re*re - phase_im*im;
      grad[tid][1][dim] = phase_re*im + phase_im*re;
    }
    re = val[tid][0];
    im = val[tid][1];
    val[tid][0] = phase_re*re - phase_im*im;
    val[tid][1] = phase_re*im + phase_im*re;
    __syncthreads();
    // Now serialize to output buffers.
    int iend = min (BS, N-block*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == 0)
      {
        outval[outIndex]     = val[i][0];
        outgrad[outIndex][0] = grad[i][0][0];
        outgrad[outIndex][1] = grad[i][0][1];
        outgrad[outIndex][2] = grad[i][0][2];
        outlap[outIndex]     = lapl[i][0];
      }
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        int off = outBlock*BS+tid;
        myVal[off] = outval[tid];
        myGradLapl[0*row_stride+off] = outgrad[tid][0];
        myGradLapl[1*row_stride+off] = outgrad[tid][1];
        myGradLapl[2*row_stride+off] = outgrad[tid][2];
        myGradLapl[3*row_stride+off] = outlap[tid];
        outBlock++;
        outIndex = 0;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == 0)
        {
          outval[outIndex]     = val[i][1];
          outgrad[outIndex][0] = grad[i][1][0];
          outgrad[outIndex][1] = grad[i][1][1];
          outgrad[outIndex][2] = grad[i][1][2];
          outlap[outIndex]     = lapl[i][1];
        }
        outIndex++;
      }
      __syncthreads();
      if (outIndex == BS)
      {
        int off = outBlock*BS+tid;
        myVal[off]                   = outval[tid];
        myGradLapl[0*row_stride+off] = outgrad[tid][0];
        myGradLapl[1*row_stride+off] = outgrad[tid][1];
        myGradLapl[2*row_stride+off] = outgrad[tid][2];
        myGradLapl[3*row_stride+off] = outlap[tid];
        outIndex = 0;
        outBlock++;
      }
      __syncthreads();
    }
    __syncthreads();
  }
  __syncthreads();
  if (tid < outIndex)
  {
    int off = outBlock*BS+tid;
    myVal[off] = outval[tid];
    myGradLapl[0*row_stride+off] = outgrad[tid][0];
    myGradLapl[1*row_stride+off] = outgrad[tid][1];
    myGradLapl[2*row_stride+off] = outgrad[tid][2];
    myGradLapl[3*row_stride+off] = outlap[tid];
  }
}


template<typename T> void
evaluate3DSplineComplexToReal
(HybridJobType *job_types, T *pos, T *kpoints, int *make2copies,
 typename SplineTraits<std::complex<T>,3>::CudaSplineType *multispline, T *Linv,
 T **vals, T **grad_lapl,
 int row_stride, int N, int numWalkers)
{
  const int BS=64;
  dim3 dimGrid(numWalkers);
  dim3 dimBlock(BS);
  evaluate3DSplineComplexToReal_kernel<T,BS><<<dimGrid,dimBlock>>>
  (job_types, pos, kpoints, make2copies,
   multispline->gridInv, (T*) multispline->coefs,
   multispline->dim, multispline->stride, Linv, vals, grad_lapl,
   row_stride, N);
}


template<typename T, int BS> __global__ void
evaluate3DSplineComplexToReal_kernel
(HybridJobType *job_types, T *pos, T *kpoints, int *make2copies,
 typename cudaTypeTraits<T>::realType3 drInv, T *coefs, uint3 dim, uint3 strides,
 T *Linv, T **vals, int N)
{
  int tid = threadIdx.x;
  __shared__ HybridJobType myjob;
  if (tid == 0) myjob = job_types[blockIdx.x];
  __syncthreads();
  if (myjob != BSPLINE_3D_JOB)
    return;
  int ir    = blockIdx.x;
  __shared__ T *myVal;
  __shared__ T r[3], u[3];
  __shared__ T G[3][3];
  if (tid < 9)
    G[0][tid] = Linv[tid];
  if (tid == 0)
    myVal      =      vals[ir];
  if (tid < 3)
  {
    r[tid] = pos[3*ir+tid];
    u[tid] = G[0][tid]*r[0] + G[1][tid]*r[1] + G[2][tid]*r[2];
    u[tid] -= floor(u[tid]);
  }
  __syncthreads();
  int3 index;
  typename cudaTypeTraits<T>::realType3 t;
  T s, sf;
  typename cudaTypeTraits<T>::realType4 tp[3];
  s = u[0] * drInv.x;
  sf = floor(s);
  index.x = min(max(0,(int)sf), dim.x-1);
  t.x = s - sf;
  s = u[1] * drInv.y;
  sf = floor(s);
  index.y = min(max(0,(int)sf), dim.y-1);
  t.y = s - sf;
  s = u[2] * drInv.z;
  sf = floor(s);
  index.z = min(max(0,(int)sf), dim.z-1);
  t.z = s - sf;
  tp[0] = cudaMakeType4<T>(t.x*t.x*t.x, t.x*t.x, t.x, 1.0f);
  tp[1] = cudaMakeType4<T>(t.y*t.y*t.y, t.y*t.y, t.y, 1.0f);
  tp[2] = cudaMakeType4<T>(t.z*t.z*t.z, t.z*t.z, t.z, 1.0f);
  // First 4 of a are value, second 4 are derivative, last four are
  // second derivative.
  __shared__ T a[4], b[4], c[4];
  if (tid < 4)
  {
    a[tid] = Acuda[4*tid+0]*tp[0].x + Acuda[4*tid+1]*tp[0].y + Acuda[4*tid+2]*tp[0].z + Acuda[4*tid+3]*tp[0].w;
    b[tid] = Acuda[4*tid+0]*tp[1].x + Acuda[4*tid+1]*tp[1].y + Acuda[4*tid+2]*tp[1].z + Acuda[4*tid+3]*tp[1].w;
    c[tid] = Acuda[4*tid+0]*tp[2].x + Acuda[4*tid+1]*tp[2].y + Acuda[4*tid+2]*tp[2].z + Acuda[4*tid+3]*tp[2].w;
  }
  __syncthreads();
  __shared__ T abc[64];
  int i = (tid>>4)&3;
  int j = (tid>>2)&3;
  int k = (tid & 3);
  if (tid < 64)
    abc[(16*i+4*j+k)]   = a[i]*b[j]*c[k];
  __syncthreads();
  int numBlocks = (N+BS-1)/BS;
  int outIndex=0, outBlock=0;
  __shared__ T outval[BS];
  //////////////////////
  // Outer block loop //
  //////////////////////
  for (int block=0; block<numBlocks; block++)
  {
    __shared__ T kp[BS][3];
    for (int i=0; i<3; i++)
    {
      int off = (3*block+i)*BS+tid;
      if (off < 3*N)
        kp[0][i*BS+tid] = kpoints[off];
    }
    __shared__ int m2c[BS];
    if (block*BS+tid < N)
      m2c[tid] = make2copies[block*BS+tid];
    __syncthreads();
    T phase_re, phase_im;
    sincos(-(r[0]*kp[tid][0] +
             r[1]*kp[tid][1] +
             r[2]*kp[tid][2]), &phase_im, &phase_re);;
    T v=0.0f;
    __shared__ T val[BS][2];
    int off   = 2*block*BS+threadIdx.x;
    int n = 0;
    T *b0 = coefs + index.x*strides.x + index.y*strides.y +
                index.z*strides.z + off;
    if (off < 2*N)
    {
      for (int i=0; i<4; i++)
      {
        for (int j=0; j<4; j++)
        {
          T *base = b0 + i*strides.x + j*strides.y;
          for (int k=0; k<4; k++)
          {
            T c  = base[k*strides.z];
            v   += abc[n+  0] * c;
            n += 1;
          }
        }
      }
    }
    val[0][tid]    = v;
    v=0.0f;
    off   = (2*block+1)*BS+threadIdx.x;
    n = 0;
    b0 = coefs + index.x*strides.x + index.y*strides.y + index.z*strides.z + off;
    if (off < 2*N)
    {
      for (int i=0; i<4; i++)
      {
        for (int j=0; j<4; j++)
        {
          T *base = b0 + i*strides.x + j*strides.y;
          for (int k=0; k<4; k++)
          {
            T c  = base[k*strides.z];
            v   += abc[n+  0] * c;
            n += 1;
          }
        }
      }
    }
    val[0][tid+BS]    = v;
    __syncthreads();
    T re, im;
    // Add phase contribution to laplacian
    re = val[tid][0];
    im = val[tid][1];
    val[tid][0] = phase_re*re - phase_im*im;
    val[tid][1] = phase_re*im + phase_im*re;
    __syncthreads();
    // Now serialize to output buffers.
    int iend = min (BS, N-block*BS);
    for (int i=0; i<iend; i++)
    {
      if (tid == 0)
      {
        outval[outIndex]     = val[i][0];
      }
      outIndex++;
      __syncthreads();
      if (outIndex == BS)
      {
        int off = outBlock*BS+tid;
        myVal[off] = outval[tid];
        outBlock++;
        outIndex = 0;
      }
      __syncthreads();
      if (m2c[i])
      {
        if (tid == 0)
        {
          outval[outIndex]     = val[i][1];
        }
        outIndex++;
      }
      __syncthreads();
      if (outIndex == BS)
      {
        int off = outBlock*BS+tid;
        myVal[off]                   = outval[tid];
        outIndex = 0;
        outBlock++;
      }
      __syncthreads();
    }
    __syncthreads();
  }
  __syncthreads();
  if (tid < outIndex)
  {
    int off = outBlock*BS+tid;
    myVal[off] = outval[tid];
  }
}


template<typename T> void
evaluate3DSplineComplexToReal
(HybridJobType *job_types, T *pos, T *kpoints, int *make2copies,
 typename SplineTraits<std::complex<T>,3>::CudaSplineType *multispline, T *Linv,
 T **vals, int N, int numWalkers)
{
  const int BS=64;
  dim3 dimGrid(numWalkers);
  dim3 dimBlock(BS);
  evaluate3DSplineComplexToReal_kernel<T,BS><<<dimGrid,dimBlock>>>
  (job_types, pos, kpoints, make2copies,
   multispline->gridInv, (T*)multispline->coefs,
   multispline->dim, multispline->stride, Linv, vals, N);
}

/**********************************************************************************
// The followings are the explicit instantiations for the above template functions //
**********************************************************************************/

// MakeHyridJobList
template void MakeHybridJobList <float>
(float*, int, float*, float*, float*, int, float*, float*,
 HybridJobType*, float*, HybridData<float>*);

template void MakeHybridJobList <double>
(double*, int, double*, double*, double*, int, double*, double*,
 HybridJobType*, double*, HybridData<double>*);

// evaluateHybridSplineReal
template void evaluateHybridSplineReal <float>
(HybridJobType*, float**, AtomicOrbitalCuda<float>*,
 HybridData<float>*, float*, float**, int, int, int);

template void evaluateHybridSplineReal <double>
(HybridJobType*, double**, AtomicOrbitalCuda<double>*,
 HybridData<double>*, double*, double**, int, int, int);

// evaluateHybridSplineReal
template void evaluateHybridSplineReal <float>
(HybridJobType*, float*, float**, float**, float**,
 AtomicOrbitalCuda<float>*, HybridData<float>*, float*,
 float**, float**, int, int, int, int);

template void evaluateHybridSplineReal <double>
(HybridJobType*, double*, double**, double**, double**,
 AtomicOrbitalCuda<double>*, HybridData<double>*, double*,
 double**, double**, int, int, int, int);

// evaluateHybridSplineComplexToReal
template void evaluateHybridSplineComplexToReal <float>
(HybridJobType*, float**, AtomicOrbitalCuda<float>*,
 HybridData<float>*, float*, int*, float**, int, int, int);

template void evaluateHybridSplineComplexToReal <double>
(HybridJobType*, double**, AtomicOrbitalCuda<double>*,
 HybridData<double>*, double*, int*, double**, int, int, int);

// evaluateHybridSplineComplexToReal
template void evaluateHybridSplineComplexToReal <float>
(HybridJobType*, float*, float**, float**, float**,
 AtomicOrbitalCuda<float>*, HybridData<float>*,
 float*, int*, float**, float**, int, int, int, int);

template void evaluateHybridSplineComplexToReal <double>
(HybridJobType*, double*, double**, double**, double**,
 AtomicOrbitalCuda<double>*, HybridData<double>*,
 double*, int*, double**, double**, int, int, int, int);

// evaluate3DSplineReal
template void evaluate3DSplineReal <float>
(HybridJobType*, float*, float*,
 typename SplineTraits<float,3>::CudaSplineType*,
 float*, float**, int, int);

template void evaluate3DSplineReal <double>
(HybridJobType*, double*, double*,
 typename SplineTraits<double,3>::CudaSplineType*,
 double*, double**, int, int);

// evaluate3DSplineReal
template void evaluate3DSplineReal <float>
(HybridJobType*, float*, float*,
 typename SplineTraits<float,3>::CudaSplineType*,
 float*, float**, float**, int, int, int);

template void evaluate3DSplineReal <double>
(HybridJobType*, double*, double*,
 typename SplineTraits<double,3>::CudaSplineType*,
 double*, double**, double**, int, int, int);

// evaluate3DSplineComplexToReal
template void evaluate3DSplineComplexToReal <float>
(HybridJobType*, float*, float*, int*,
 typename SplineTraits<std::complex<float>,3>::CudaSplineType*,
 float*, float**, int, int);

template void evaluate3DSplineComplexToReal <double>
(HybridJobType*, double*, double*, int*,
 typename SplineTraits<std::complex<double>,3>::CudaSplineType*,
 double*, double**, int, int);

// evaluate3DSplineComplexToReal
template void evaluate3DSplineComplexToReal <float>
(HybridJobType*, float*, float*, int*,
 typename SplineTraits<std::complex<float>,3>::CudaSplineType*,
 float*, float**, float**, int, int, int);

template void evaluate3DSplineComplexToReal <double>
(HybridJobType*, double*, double*, int*,
 typename SplineTraits<std::complex<double>,3>::CudaSplineType*,
 double*, double**, double**, int, int, int);

// CalcYlmRealCuda
template void CalcYlmRealCuda <float>
(float*, HybridJobType*, float**, float**, float**, int, int);

template void CalcYlmRealCuda <double>
(double*, HybridJobType*, double**, double**, double**, int, int);

// CalcYlmComplexCuda
template void CalcYlmComplexCuda <float> 
(float*, HybridJobType*, float**, float**, float**, int, int); 

template void CalcYlmComplexCuda <double> 
(double*, HybridJobType*, double**, double**, double**, int, int); 


////////////////////
// TEST GPU Y_lm  //
////////////////////


#ifdef TEST_GPU_YLM

class Vec3
{
private:
  double r[3];
public:
  inline double  operator[](int i) const
  {
    return r[i];
  }
  inline double& operator[](int i)
  {
    return r[i];
  }
  Vec3(double x, double y, double z)
  {
    r[0]=x;
    r[1]=y;
    r[2]=z;
  }
  Vec3() { }
};


// Fast implementation
// See Geophys. J. Int. (1998) 135,pp.307-309
void
CalcYlm (Vec3 rhat,
         vector<complex<double> > &Ylm,
         vector<complex<double> > &dYlm_dtheta,
         vector<complex<double> > &dYlm_dphi,
         int lMax)
{
  const double fourPiInv = 0.0795774715459477;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cottheta = costheta/sintheta;
  double cosphi, sinphi;
  cosphi=rhat[0]/sintheta;
  sinphi=rhat[1]/sintheta;
  complex<double> e2iphi(cosphi, sinphi);
  double lsign = 1.0;
  double dl = 0.0;
  double XlmVec[2*lMax+1], dXlmVec[2*lMax+1];
  for (int l=0; l<=lMax; l++)
  {
    XlmVec[2*l]  = lsign;
    dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
    XlmVec[0]    = lsign*XlmVec[2*l];
    dXlmVec[0]   = lsign*dXlmVec[2*l];
    double dm = dl;
    double msign = lsign;
    for (int m=l; m>0; m--)
    {
      double tmp = std::sqrt((dl+dm)*(dl-dm+1.0));
      XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
      dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
      // Copy to negative m
      XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
      dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
      msign *= -1.0;
      dm -= 1.0;
    }
    double sum = 0.0;
    for (int m=-l; m<=l; m++)
      sum += XlmVec[l+m]*XlmVec[l+m];
    // Now, renormalize the Ylms for this l
    double norm = std::sqrt((2.0*dl+1.0)*fourPiInv / sum);
    for (int m=-l; m<=l; m++)
    {
      XlmVec[l+m]  *= norm;
      dXlmVec[l+m] *= norm;
    }
    // Multiply by azimuthal phase and store in Ylm
    complex<double> e2imphi (1.0, 0.0);
    complex<double> eye(0.0, 1.0);
    for (int m=0; m<=l; m++)
    {
      Ylm[l*(l+1)+m]  =  XlmVec[l+m]*e2imphi;
      Ylm[l*(l+1)-m]  =  XlmVec[l-m]*conj(e2imphi);
      dYlm_dphi[l*(l+1)+m ]  =  (double)m * eye *XlmVec[l+m]*e2imphi;
      dYlm_dphi[l*(l+1)-m ]  = -(double)m * eye *XlmVec[l-m]*conj(e2imphi);
      dYlm_dtheta[l*(l+1)+m] = dXlmVec[l+m]*e2imphi;
      dYlm_dtheta[l*(l+1)-m] = dXlmVec[l-m]*conj(e2imphi);
      e2imphi *= e2iphi;
    }
    dl += 1.0;
    lsign *= -1.0;
  }
}

// Fast implementation
// See Geophys. J. Int. (1998) 135,pp.307-309
void
CalcYlm (Vec3 rhat,
         vector<double> &Ylm,
         vector<double> &dYlm_dtheta,
         vector<double> &dYlm_dphi,
         int lMax)
{
  const double fourPiInv = 0.0795774715459477;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cottheta = costheta/sintheta;
  double cosphi, sinphi;
  cosphi=rhat[0]/sintheta;
  sinphi=rhat[1]/sintheta;
  complex<double> e2iphi(cosphi, sinphi);
  double lsign = 1.0;
  double dl = 0.0;
  double XlmVec[2*lMax+1], dXlmVec[2*lMax+1];
  for (int l=0; l<=lMax; l++)
  {
    XlmVec[2*l]  = lsign;
    dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
    XlmVec[0]    = lsign*XlmVec[2*l];
    dXlmVec[0]   = lsign*dXlmVec[2*l];
    double dm = dl;
    double msign = lsign;
    for (int m=l; m>0; m--)
    {
      double tmp = std::sqrt((dl+dm)*(dl-dm+1.0));
      XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
      dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
      // Copy to negative m
      XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
      dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
      msign *= -1.0;
      dm -= 1.0;
    }
    double sum = 0.0;
    for (int m=-l; m<=l; m++)
      sum += XlmVec[l+m]*XlmVec[l+m];
    // Now, renormalize the Ylms for this l
    double norm = std::sqrt((2.0*dl+1.0)*fourPiInv / sum);
    for (int m=-l; m<=l; m++)
    {
      XlmVec[l+m]  *= norm;
      dXlmVec[l+m] *= norm;
    }
    // Multiply by azimuthal phase and store in Ylm
    Ylm[l*(l+1)]         =  XlmVec[l];
    dYlm_dphi[l*(l+1) ]  = 0.0;
    dYlm_dtheta[l*(l+1)] = dXlmVec[l];
    complex<double> e2imphi = e2iphi;
    for (int m=1; m<=l; m++)
    {
      Ylm[l*(l+1)+m]         =  XlmVec[l+m]*e2imphi.real();
      Ylm[l*(l+1)-m]         =  XlmVec[l+m]*e2imphi.imag();
      dYlm_dphi[l*(l+1)+m ]  = -(double)m * XlmVec[l+m] *e2imphi.imag();
      dYlm_dphi[l*(l+1)-m ]  =  (double)m * XlmVec[l+m] *e2imphi.real();
      dYlm_dtheta[l*(l+1)+m] = dXlmVec[l+m]*e2imphi.real();
      dYlm_dtheta[l*(l+1)-m] = dXlmVec[l+m]*e2imphi.imag();
      e2imphi *= e2iphi;
    }
    dl += 1.0;
    lsign *= -1.0;
  }
}


#include <stdlib.h>

void TestYlmComplex()
{
  int numr = 1000;
  const int BS=32;
  float *rhat_device, *Ylm_device, *dtheta_device, *dphi_device;
  float **Ylm_ptr, **dtheta_ptr, **dphi_ptr;
  const int lmax = 5;
  const int numlm = (lmax+1)*(lmax+1);
  cudaMalloc ((void**)&rhat_device, 3*sizeof(float)*numr);
  cudaMalloc ((void**)&Ylm_device, 2*numlm*sizeof(float)*numr);
  cudaMalloc ((void**)&dtheta_device, 2*numlm*sizeof(float)*numr);
  cudaMalloc ((void**)&dphi_device, 2*numlm*sizeof(float)*numr);
  cudaMalloc ((void**)&Ylm_ptr,    numr*sizeof(float*));
  cudaMalloc ((void**)&dtheta_ptr, numr*sizeof(float*));
  cudaMalloc ((void**)&dphi_ptr,   numr*sizeof(float*));
  float *Ylm_host[numr], *dtheta_host[numr], *dphi_host[numr];
  float rhost[3*numr];
  vector<Vec3> rlist;
  for (int i=0; i<numr; i++)
  {
    Vec3 r;
    r[0] = 2.0*drand48()-1.0;
    r[1] = 2.0*drand48()-1.0;
    r[2] = 2.0*drand48()-1.0;
    double nrm = 1.0/std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    r[0] *= nrm;
    r[1] *= nrm;
    r[2] *= nrm;
    rlist.push_back(r);
    rhost[3*i+0]=r[0];
    rhost[3*i+1]=r[1];
    rhost[3*i+2]=r[2];
    Ylm_host[i] = Ylm_device+2*i*numlm;
    dtheta_host[i] = dtheta_device+2*i*numlm;
    dphi_host[i]   = dphi_device + 2*i*numlm;
  }
  cudaMemcpyAsync(rhat_device, rhost, 3*numr*sizeof(float),  cudaMemcpyHostToDevice);
  cudaMemcpyAsync(Ylm_ptr, Ylm_host, numr*sizeof(float*),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(dtheta_ptr, dtheta_host, numr*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(dphi_ptr,  dphi_host, numr*sizeof(float*), cudaMemcpyHostToDevice);
  dim3 dimBlock(BS);
  dim3 dimGrid((numr+BS-1)/BS);
  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++)
  {
    CalcYlmComplex<float,5,BS><<<dimGrid,dimBlock>>>
    (rhat_device, Ylm_ptr, dtheta_ptr, dphi_ptr, numr);
  }
  cudaThreadSynchronize();
  end = clock();
  fprintf (stderr, "Ylm rate = %1.8f\n",
           10000*numr/((double)(end-start)/(double)CLOCKS_PER_SEC));
  complex<float> Ylm[numr*numlm], dtheta[numr*numlm], dphi[numr*numlm];
  cudaMemcpy(Ylm, Ylm_device, 2*numr*numlm*sizeof(float),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(dtheta, dtheta_device, 2*numr*numlm*sizeof(float),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(dphi, dphi_device, 2*numr*numlm*sizeof(float),
             cudaMemcpyDeviceToHost);
  int n = 999;
  vector<complex<double> > Ylm_cpu(numlm), dtheta_cpu(numlm), dphi_cpu(numlm);
  CalcYlm (rlist[n], Ylm_cpu, dtheta_cpu, dphi_cpu, lmax);
  fprintf (stderr, "Ylm:\n");
  for (int lm=0; lm<numlm; lm++)
  {
    fprintf(stderr, "%12.7f %12.7f   %12.7f %12.7f  %3.0f %3.0f\n",
            Ylm_cpu[lm].real(), Ylm_cpu[lm].imag(),
            Ylm[lm+n*numlm].real(), Ylm[lm+n*numlm].imag(),
            Ylm_cpu[lm].real()/Ylm[lm+n*numlm].real(),
            Ylm_cpu[lm].imag()/Ylm[lm+n*numlm].imag());
  }
  fprintf (stderr, "dtheta:\n");
  for (int lm=0; lm<numlm; lm++)
  {
    fprintf(stderr, "%12.6f %12.6f   %12.6f %12.6f  %3.0f %3.0f\n",
            dtheta_cpu[lm].real(), dtheta_cpu[lm].imag(),
            dtheta[lm+n*numlm].real(), dtheta[lm+n*numlm].imag(),
            dtheta_cpu[lm].real()/dtheta[lm+n*numlm].real(),
            dtheta_cpu[lm].imag()/dtheta[lm+n*numlm].imag());
  }
  fprintf (stderr, "dphi:\n");
  for (int lm=0; lm<numlm; lm++)
  {
    fprintf(stderr, "%12.6f %12.6f   %12.6f %12.6f  %3.0f %3.0f\n",
            dphi_cpu[lm].real(), dphi_cpu[lm].imag(),
            dphi[lm+n*numlm].real(), dphi[lm+n*numlm].imag(),
            dphi_cpu[lm].real()/dphi[lm+n*numlm].real(),
            dphi_cpu[lm].imag()/dphi[lm+n*numlm].imag());
  }
}


void TestYlmReal()
{
  int numr = 1000;
  const int BS=32;
  float *rhat_device, *Ylm_device, *dtheta_device, *dphi_device;
  float **Ylm_ptr, **dtheta_ptr, **dphi_ptr;
  const int lmax = 5;
  const int numlm = (lmax+1)*(lmax+1);
  int block_size = ((numlm+15)/16)*16;
  cudaMalloc ((void**)&rhat_device,   3*sizeof(float)*numr);
  cudaMalloc ((void**)&Ylm_device,    block_size*sizeof(float)*numr);
  cudaMalloc ((void**)&dtheta_device, block_size*sizeof(float)*numr);
  cudaMalloc ((void**)&dphi_device,   block_size*sizeof(float)*numr);
  cudaMalloc ((void**)&Ylm_ptr,       numr*sizeof(float*));
  cudaMalloc ((void**)&dtheta_ptr,    numr*sizeof(float*));
  cudaMalloc ((void**)&dphi_ptr,      numr*sizeof(float*));
  float *Ylm_host[numr], *dtheta_host[numr], *dphi_host[numr];
  float rhost[3*numr];
  vector<Vec3> rlist;
  for (int i=0; i<numr; i++)
  {
    Vec3 r;
    r[0] = 2.0*drand48()-1.0;
    r[1] = 2.0*drand48()-1.0;
    r[2] = 2.0*drand48()-1.0;
    double nrm = 1.0/std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    r[0] *= nrm;
    r[1] *= nrm;
    r[2] *= nrm;
    rlist.push_back(r);
    rhost[3*i+0]=r[0];
    rhost[3*i+1]=r[1];
    rhost[3*i+2]=r[2];
    Ylm_host[i]    = Ylm_device    + i*block_size;
    dtheta_host[i] = dtheta_device + i*block_size;
    dphi_host[i]   = dphi_device   + i*block_size;
  }
  cudaMemcpyAsync(rhat_device, rhost, 3*numr*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(Ylm_ptr, Ylm_host, numr*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(dtheta_ptr, dtheta_host, numr*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(dphi_ptr,  dphi_host, numr*sizeof(float*), cudaMemcpyHostToDevice);
  dim3 dimBlock(BS);
  dim3 dimGrid((numr+BS-1)/BS);
  clock_t start, end;
  start = clock();
  for (int i=0; i<10000; i++)
  {
    CalcYlmReal<float,lmax,BS><<<dimGrid,dimBlock>>>
    (rhat_device, Ylm_ptr, dtheta_ptr, dphi_ptr, numr);
  }
  cudaThreadSynchronize();
  end = clock();
  fprintf (stderr, "Ylm rate = %1.8f\n",
           10000*numr/((double)(end-start)/(double)CLOCKS_PER_SEC));
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf (stderr, "CUDA error in CalcYlmReal:\n  %s\n",
             cudaGetErrorString(err));
    abort();
  }
  float Ylm[numr*block_size], dtheta[numr*block_size], dphi[numr*block_size];
  cudaMemcpy(Ylm, Ylm_device, numr*block_size*sizeof(float),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(dtheta, dtheta_device, numr*block_size*sizeof(float),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(dphi, dphi_device, numr*block_size*sizeof(float),
             cudaMemcpyDeviceToHost);
  int n = 999;
  vector<double> Ylm_cpu(numlm), dtheta_cpu(numlm), dphi_cpu(numlm);
  CalcYlm (rlist[n], Ylm_cpu, dtheta_cpu, dphi_cpu, lmax);
  fprintf (stderr, "Ylm:\n");
  for (int lm=0; lm<numlm; lm++)
  {
    fprintf(stderr, "%12.7f %12.7f %3.0f\n",
            Ylm_cpu[lm],
            Ylm[lm+n*block_size],
            Ylm_cpu[lm]/Ylm[lm+n*block_size]);
  }
  fprintf (stderr, "dtheta:\n");
  for (int lm=0; lm<numlm; lm++)
  {
    fprintf(stderr, "%12.6f %12.6f %3.0f \n",
            dtheta_cpu[lm],
            dtheta[lm+n*block_size],
            dtheta_cpu[lm]/dtheta[lm+n*block_size]);
  }
  fprintf (stderr, "dphi:\n");
  for (int lm=0; lm<numlm; lm++)
  {
    fprintf(stderr, "%12.6f %12.6f %3.0f\n",
            dphi_cpu[lm],
            dphi[lm+n*block_size],
            dphi_cpu[lm]/dphi[lm+n*block_size]);
  }
}



main()
{
  TestYlmComplex();
  TestYlmReal();
}
#endif
