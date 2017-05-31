//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef CUDA_DETERMINANT_UPDATE_H
#define CUDA_DETERMINANT_UPDATE_H

#include <complex>

struct updateJob
{
  void *A, *Ainv, *newRow, *AinvDelta, *AinvColk, *gradLapl, *newGradLapl, *dummy;
  int iat;
};



void
update_inverse_cuda(float *A_g[], float *Ainv_g[], float *u_g[],
                    float *Ainv_delta_g[], float *Ainv_colk_g[],
                    int N, int rowstride, int iat, int numWalkers);

void
update_inverse_cuda(updateJob jobList[], double dummy,
                    int N, int rowstride, int numWalkers);
void
update_inverse_cuda(updateJob jobList[], float dummy,
                    int N, int rowstride, int numWalkers);

void
update_inverse_cuda(double *A_g[], double *Ainv_g[], double *u_g[],
                    double *Ainv_delta_g[], double *Ainv_colk_g[],
                    int N, int rowstride, int iat, int numWalkers);


/////////////////////////////////////////////////
// New version with fewer PCI transfers needed //
/////////////////////////////////////////////////
// These are for nonlocal updates, in which the electron
// is different for each block
void
update_inverse_cuda(float **data, int iat[],
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);

void
update_inverse_cuda(double **data, int iat[],
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);

#ifdef QMC_COMPLEX
void
update_inverse_cuda(std::complex<float> **data, int iat[],
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);

void
update_inverse_cuda(std::complex<double> **data, int iat[],
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);
#endif


// These are for single-particle move updates.  Each
// walker is moving the same electron.
void
update_inverse_cuda(float **data, int iat,
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);

void
update_inverse_cuda(double **data, int iat,
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);

#ifdef QMC_COMPLEX
void
update_inverse_cuda(std::complex<float>** data, int iat,
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);

void
update_inverse_cuda(std::complex<double>** data, int iat,
                    int A_off, int Ainv_off, int newRow_off,
                    int AinvDelta_off, int AinvColk_off,
                    int N, int rowstride, int numWalkers);
#endif


void
determinant_ratios_cuda (float *Ainv_list[], float *new_row_list[],
                         float *ratios, int N, int row_stride, int iat,
                         int numWalkers);

void
determinant_ratios_cuda (double *Ainv_list[], double *new_row_list[],
                         double *ratios, int N, int row_stride, int iat,
                         int numWalkers);

#ifdef QMC_COMPLEX
void
determinant_ratios_cuda (std::complex<float> *Ainv_list[], std::complex<float> *new_row_list[],
                         std::complex<float> *ratios, int N, int row_stride, int iat,
                         int numWalkers);

void
determinant_ratios_cuda (std::complex<double> *Ainv_list[], std::complex<double> *new_row_list[],
                         std::complex<double> *ratios, int N, int row_stride, int iat,
                         int numWalkers);
#endif


void
calc_many_ratios (float *Ainv_list[], float *new_row_list[],
                  float* ratio_list[], int num_ratio_list[],
                  int N, int row_stride, int elec_list[],
                  int numWalkers);

void
calc_many_ratios (double *Ainv_list[], double *new_row_list[],
                  double* ratio_list[], int num_ratio_list[],
                  int N, int row_stride, int elec_list[],
                  int numWalkers);

#ifdef QMC_COMPLEX
void
calc_many_ratios (std::complex<float> *Ainv_list[], std::complex<float> *new_row_list[],
                  std::complex<float>* ratio_list[], int num_ratio_list[],
                  int N, int row_stride, int elec_list[],
                  int numWalkers);

void
calc_many_ratios (std::complex<double> *Ainv_list[], std::complex<double> *new_row_list[],
                  std::complex<double>* ratio_list[], int num_ratio_list[],
                  int N, int row_stride, int elec_list[],
                  int numWalkers);
#endif

void
scale_grad_lapl(float *grad_list[], float *hess_list[],
                float *grad_lapl_list[], float Linv[], int N,
                int num_walkers);

void
calc_grad_lapl (float *Ainv_list[], float *grad_lapl_list[],
                float *out_list[], int N, int row_stride, int num_mats);

void
calc_grad_lapl (double *Ainv_list[], double *grad_lapl_list[],
                double *out_list[], int N, int row_stride, int num_mats);


#ifdef QMC_COMPLEX
void
calc_grad_lapl (std::complex<float> *Ainv_list[], std::complex<float> *grad_lapl_list[],
                std::complex<float> *out_list[], int N, int row_stride, int num_mats);

void
calc_grad_lapl (std::complex<double> *Ainv_list[], std::complex<double> *grad_lapl_list[],
                std::complex<double> *out_list[], int N, int row_stride, int num_mats);
#endif


void
multi_copy (float *dest[], float *src[], int len, int num);
void
multi_copy (double *dest[], double *src[], int len, int num);

void
multi_copy (float *buff[], int dest_off, int src_off, int len, int num);
void
multi_copy (double *buff[], int dest_off, int src_off, int len, int num);


#ifdef QMC_COMPLEX
void
multi_copy (std::complex<float> *dest[], std::complex<float> *src[], int len, int num);
void
multi_copy (std::complex<double> *dest[], std::complex<double> *src[], int len, int num);

void
multi_copy (std::complex<float> *buff[], int dest_off, int src_off, int len, int num);
void
multi_copy (std::complex<double> *buff[], int dest_off, int src_off, int len, int num);
#endif

//YingWai: the following two are repeated   (Dec 28, 15)
/*
void
calc_many_ratios (float *Ainv_list[], float *new_row_list[],
                  float* ratio_list[], int num_ratio_list[],
                  int N, int row_stride, int elec_list[],
                  int numWalkers);

void
calc_many_ratios (double *Ainv_list[], double *new_row_list[],
                  double* ratio_list[], int num_ratio_list[],
                  int N, int row_stride, int elec_list[],
                  int numWalkers);
*/

void
determinant_ratios_grad_lapl_cuda
(float *Ainv_list[], float *new_row_list[],
 float *grad_lapl_list[], float ratios_grad_lapl[],
 int N, int row_stride, int iat, int numWalkers);

void
determinant_ratios_grad_lapl_cuda
(double *Ainv_list[], double *new_row_list[],
 double *grad_lapl_list[], double ratios_grad_lapl[],
 int N, int row_stride, int iat, int numWalkers);

#ifdef QMC_COMPLEX
void
determinant_ratios_grad_lapl_cuda
(std::complex<float> *Ainv_list[], std::complex<float> *new_row_list[],
 std::complex<float> *grad_lapl_list[], std::complex<float> ratios_grad_lapl[],
 int N, int row_stride, int iat, int numWalkers);

void
determinant_ratios_grad_lapl_cuda
(std::complex<double> *Ainv_list[], std::complex<double> *new_row_list[],
 std::complex<double> *grad_lapl_list[], std::complex<double> ratios_grad_lapl[], 
 int N, int row_stride, int iat, int numWalkers);
#endif

void
determinant_ratios_grad_lapl_cuda
(float *Ainv_list[], float *new_row_list[], float *grad_lapl_list[],
 float ratios_grad_lapl[], int N, int row_stride,
 int iat_list[], int numWalkers);

void
determinant_ratios_grad_lapl_cuda
(double *Ainv_list[], double *new_row_list[], double *grad_lapl_list[],
 double ratios_grad_lapl[], int N, int row_stride,
 int iat_list[], int numWalkers);

#ifdef QMC_COMPLEX
void
determinant_ratios_grad_lapl_cuda
(std::complex<float> *Ainv_list[], std::complex<float> *new_row_list[],
 std::complex<float> *grad_lapl_list[], std::complex<float> ratios_grad_lapl[], 
 int N, int row_stride, int iat_list[], int numWalkers);

void
determinant_ratios_grad_lapl_cuda
(std::complex<double> *Ainv_list[], std::complex<double> *new_row_list[], 
 std::complex<double> *grad_lapl_list[], std::complex<double> ratios_grad_lapl[], 
 int N, int row_stride, int iat_list[], int numWalkers);
#endif

void
calc_gradient (float *Ainv_list[], float *grad_lapl_list[],
               float grad[], int N, int row_stride, int elec,
               int numWalkers);

void
calc_gradient (double *Ainv_list[], double *grad_lapl_list[],
               double grad[], int N, int row_stride, int elec,
               int numWalkers);

#ifdef QMC_COMPLEX
void
calc_gradient (std::complex<float> *Ainv_list[], std::complex<float> *grad_lapl_list[],
               std::complex<float> grad[], int N, int row_stride, int elec,
               int numWalkers);

void
calc_gradient (std::complex<double> *Ainv_list[], std::complex<double> *grad_lapl_list[],
               std::complex<double> grad[], int N, int row_stride, int elec,
               int numWalkers);
#endif

#endif
