//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef ATOMIC_ORBITAL_CUDA_H
#define ATOMIC_ORBITAL_CUDA_H

#include "einspline/multi_bspline.h"
#include "einspline/multi_bspline_create_cuda.h"
#include <cuda.h>

// type traits for cuda variable types
template<typename T> struct cudaTypeTraits;

template<> struct cudaTypeTraits<float>
{
  typedef float1 realType1; 
  typedef float2 realType2; 
  typedef float3 realType3; 
  typedef float4 realType4; 
};

template<> struct cudaTypeTraits<double>
{ 
  typedef double1 realType1;
  typedef double2 realType2;
  typedef double3 realType3;
  typedef double4 realType4;
};


// type traits for CudaSplineType similar to MultiOrbitalTraits
template<typename StorageType, int dim> struct SplineTraits {}; 

// 2D
template<> struct SplineTraits<float,2>
{
  typedef multi_UBspline_2d_s_cuda CudaSplineType;
};

template<> struct SplineTraits<double,2>
{
  typedef multi_UBspline_2d_d_cuda CudaSplineType;
};

template<> struct SplineTraits<std::complex<float>,2>
{
  typedef multi_UBspline_2d_c_cuda CudaSplineType;
};

template<> struct SplineTraits<std::complex<double>,2>
{
  typedef multi_UBspline_2d_z_cuda CudaSplineType;
};

// 3D
template<> struct SplineTraits<float,3>
{
  typedef multi_UBspline_3d_s_cuda CudaSplineType;
};

template<> struct SplineTraits<double,3>
{
  typedef multi_UBspline_3d_d_cuda CudaSplineType;
};

template<> struct SplineTraits<std::complex<float>,3>
{
  typedef multi_UBspline_3d_c_cuda CudaSplineType;
};

template<> struct SplineTraits<std::complex<double>,3>
{
  typedef multi_UBspline_3d_z_cuda CudaSplineType;
};


typedef enum { BSPLINE_3D_JOB, ATOMIC_POLY_JOB, ATOMIC_SPLINE_JOB } HybridJobType;


template<typename T>
class HybridData
{
public:
  // Integer specifying which image of the ion this electron is near
  T img[3];
  // Minimum image distance to the ion;
  T dist;
  // The index the ion this electron is near
  int ion;
  int lMax;
  int PAD[2];
};


template<typename T>
class AtomicOrbitalCuda
{
public:
  int lMax, spline_stride, lm_stride;
  int poly_order, poly_stride;
  T spline_dr_inv;
  T *spline_coefs, *poly_coefs;
  int PAD[6];
};


void init_atomic_cuda();


template<typename T> void
MakeHybridJobList (T* elec_list, int num_elecs, T* ion_list,
                   T* poly_radii, T* spline_radii,
                   int num_ions, T *L, T *Linv,
                   HybridJobType *job_list, T *rhat_list,
                   HybridData<T> *data_list);

template<typename T> void
evaluateHybridSplineReal (HybridJobType *job_types, T **Ylm_real,
                          AtomicOrbitalCuda<T> *orbitals,
                          HybridData<T> *data, T *k_reduced,
                          T **vals, int N, int numWalkers, int lMax);

template<typename T> void
evaluateHybridSplineReal (HybridJobType *job_types, T *rhats,
                          T **Ylm_real, T **dYlm_dTheta, T **dYlm_dphi,
                          AtomicOrbitalCuda<T> *orbitals,
                          HybridData<T> *data, T *k_reduced,
                          T **vals, T **grad_lapl,
                          int row_stride, int N, int numWalkers, int lMax);

template<typename T> void
evaluateHybridSplineComplexToReal
(HybridJobType *job_types, T **Ylm_real,
 AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int *make2copies,
 T **vals, int N, int numWalkers, int lMax);

template<typename T> void
evaluateHybridSplineComplexToReal
(HybridJobType *job_types, T *rhats,
 T **Ylm, T **dYlm_dTheta, T **dYlm_dphi,
 AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int *make2copies,
 T **vals, T **grad_lapl,
 int row_stride, int N, int numWalkers, int lMax);

template<typename T> void
evaluate3DSplineReal (HybridJobType *job_types, T *pos, T *kpoints,
                      typename SplineTraits<T,3>::CudaSplineType *multispline, T *Linv,
                      T **vals, int N, int numWalkers);

template<typename T> void
evaluate3DSplineReal (HybridJobType *job_types, T *pos, T *kpoints,
                      typename SplineTraits<T,3>::CudaSplineType *multispline, T *Linv,
                      T **vals, T **grad_lapl,
                      int row_stride, int N, int numWalkers);

template<typename T> void
evaluate3DSplineComplexToReal
(HybridJobType *job_types, T *pos, T *kpoints, int *make2copies,
 typename SplineTraits<std::complex<T>,3>::CudaSplineType *multispline, T *Linv,
 T **vals, int N, int numWalkers);

template<typename T> void
evaluate3DSplineComplexToReal
(HybridJobType *job_types, T *pos, T *kpoints, int *make2copies,
 typename SplineTraits<std::complex<T>,3>::CudaSplineType *multispline, T *Linv,
 T **vals, T **grad_lapl,
 int row_stride, int N, int numWalkers);

template<typename T> void
CalcYlmRealCuda (T *rhats,  HybridJobType *job_type,
                 T **Ylm_ptr, T **dYlm_dtheta_ptr, T **dYlm_dphi_ptr,
                 int lMax, int N);

template<typename T> void
CalcYlmComplexCuda (T *rhats,  HybridJobType *job_type,
                    T **Ylm_ptr, T **dYlm_dtheta_ptr, T **dYlm_dphi_ptr,
                    int lMax, int N);


// YingWai: these seems to be unused  (Oct 1, 15)
/*
template<typename T> void
evaluateHybridSplineComplexToRealNLPP
(HybridJobType *job_types,
 T **Ylm_real, AtomicOrbitalCuda<T> *orbitals,
 HybridData<T> *data, T *k_reduced, int* make2copies,
 T **vals, int N, int numWalkers, int lMax);

template<typename T> void
CalcYlmRealCuda (T *rhats,  HybridJobType *job_type,
                 T **Ylm_ptr, int lMax, int N);

template<typename T> void
CalcYlmComplexCuda (T *rhats,  HybridJobType *job_type,
                    T **Ylm_ptr, int lMax, int N);
*/

#endif
