#ifndef ATOMIC_ORBITAL_CUDA_H
#define ATOMIC_ORBITAL_CUDA_H

#include <einspline/multi_bspline.h>
#include <einspline/multi_bspline_create_cuda.h>


typedef enum { BSPLINE_3D_JOB, ATOMIC_POLY_JOB, ATOMIC_SPLINE_JOB } HybridJobType;

struct HybridDataFloat
{
  // Integer specifying which image of the ion this electron is near
  float img[3];
  // Minimum image distance to the ion;
  float dist;
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

void
MakeHybridJobList (float* elec_list, int num_elecs, float* ion_list,
                   float* poly_radii, float* spline_radii,
                   int num_ions, float *L, float *Linv,
                   HybridJobType *job_list, float *rhat_list,
                   HybridDataFloat *data_list);

void
evaluateHybridSplineReal (HybridJobType *job_types, float **Ylm_real,
                          AtomicOrbitalCuda<float> *orbitals,
                          HybridDataFloat *data, float *k_reduced,
                          float **vals, int N, int numWalkers, int lMax);
void
evaluateHybridSplineReal (HybridJobType *job_types, float *rhats,
                          float **Ylm_real, float **dYlm_dTheta, float **dYlm_dphi,
                          AtomicOrbitalCuda<float> *orbitals,
                          HybridDataFloat *data, float *k_reduced,
                          float **vals, float **grad_lapl,
                          int row_stride, int N, int numWalkers, int lMax);

void
evaluateHybridSplineComplexToReal
(HybridJobType *job_types, float **Ylm_real,
 AtomicOrbitalCuda<float> *orbitals,
 HybridDataFloat *data, float *k_reduced, int *make2copies,
 float **vals, int N, int numWalkers, int lMax);

void
evaluateHybridSplineComplexToRealNLPP
(HybridJobType *job_types,
 float **Ylm_real, AtomicOrbitalCuda<float> *orbitals,
 HybridDataFloat *data, float *k_reduced, int* make2copies,
 float **vals, int N, int numWalkers, int lMax);

void
evaluateHybridSplineComplexToReal
(HybridJobType *job_types, float *rhats,
 float **Ylm, float **dYlm_dTheta, float **dYlm_dphi,
 AtomicOrbitalCuda<float> *orbitals,
 HybridDataFloat *data, float *k_reduced, int *make2copies,
 float **vals, float **grad_lapl,
 int row_stride, int N, int numWalkers, int lMax);



void
evaluate3DSplineReal (HybridJobType *job_types, float *pos, float *kpoints,
                      multi_UBspline_3d_s_cuda *multispline, float *Linv,
                      float **vals, int N, int numWalkers);

void
evaluate3DSplineReal (HybridJobType *job_types, float *pos, float *kpoints,
                      multi_UBspline_3d_s_cuda *multispline, float *Linv,
                      float **vals, float **grad_lapl,
                      int row_stride, int N, int numWalkers);

void
evaluate3DSplineComplexToReal
(HybridJobType *job_types, float *pos, float *kpoints, int *make2copies,
 multi_UBspline_3d_c_cuda *multispline, float *Linv,
 float **vals, int N, int numWalkers);

void
evaluate3DSplineComplexToReal
(HybridJobType *job_types, float *pos, float *kpoints, int *make2copies,
 multi_UBspline_3d_c_cuda *multispline, float *Linv,
 float **vals, float **grad_lapl,
 int row_stride, int N, int numWalkers);



void CalcYlmRealCuda (float *rhats,  HybridJobType *job_type,
                      float **Ylm_ptr, float **dYlm_dtheta_ptr, float **dYlm_dphi_ptr,
                      int lMax, int N);

void CalcYlmComplexCuda (float *rhats,  HybridJobType *job_type,
                         float **Ylm_ptr, float **dYlm_dtheta_ptr, float **dYlm_dphi_ptr,
                         int lMax, int N);

void CalcYlmRealCuda (float *rhats,  HybridJobType *job_type,
                      float **Ylm_ptr, int lMax, int N);

void CalcYlmComplexCuda (float *rhats,  HybridJobType *job_type,
                         float **Ylm_ptr, int lMax, int N);

#endif
