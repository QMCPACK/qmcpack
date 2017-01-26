//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef BSPLINE_JASTROW_CUDA_H
#define BSPLINE_JASTROW_CUDA_H

#include <complex>
#include "NLjobGPU.h"

///////////////////////
// Two-Body routines //
///////////////////////

void
two_body_sum (float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
              float spline_coefs[], int numCoefs, float rMax,
              float sum[], int numWalkers);

void
two_body_sum (double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
              double spline_coefs[], int numCoefs, double rMax,
              double sum[], int numWalkers);

void
two_body_ratio (float *R[], int first, int last,
                float Rnew[], int inew,
                float spline_coefs[], int numCoefs, float rMax,
                float sum[], int numWalkers);

void
two_body_ratio (double *R[], int first, int last,
                double Rnew[], int inew,
                double spline_coefs[], int numCoefs, double rMax,
                double sum[], int numWalkers);

void
two_body_ratio_grad(float *R[], int first, int last,
                    float  Rnew[], int inew,
                    float spline_coefs[], int numCoefs, float rMax,
                    bool zero, float ratio_grad[], int numWalkers);

void
two_body_ratio_grad(double *R[], int first, int last,
                    double  Rnew[], int inew,
                    double spline_coefs[], int numCoefs, double rMax,
                    bool zero,
                    double ratio_grad[], int numWalkers);

void
two_body_NLratios(NLjobGPU<float> jobs[], int first, int last,
                  float* spline_coefs[], int numCoefs[], float rMax[],
                  int numjobs);

void
two_body_NLratios(NLjobGPU<double> jobs[], int first, int last,
                  double* spline_coefs[], int numCoefs[], double rMax[],
                  int numjobs);


void
two_body_update(float *R[], int N, int iat, int numWalkers);

void
two_body_update(double *R[], int N, int iat, int numWalkers);


void
two_body_grad_lapl(float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                   float spline_coefs[], int numCoefs, float rMax,
                   float gradLapl[], int row_stride, int numWalkers);

void
two_body_grad_lapl(double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                   double spline_coefs[], int numCoefs, double rMax,
                   double gradLapl[], int row_stride, int numWalkers);


void
two_body_gradient (float *R[], int first, int last, int iat,
                   float spline_coefs[], int numCoefs, float rMax,
                   bool zeroOut, float grad[], int numWalkers);

void
two_body_gradient (double *R[], int first, int last, int iat,
                   double spline_coefs[], int numCoefs, double rMax,
                   bool zeroOut, double grad[], int numWalkers);

void
two_body_derivs(float *R[], float *gradLogPsi[], int e1_first, int e1_last,
                int e2_first, int e2_last,
                int numCoefs, float rMax,
                float *derivs[], int numWalkers);
void
two_body_derivs(double *R[], double *gradLogPsi[], int e1_first, int e1_last,
                int e2_first, int e2_last,
                int numCoefs, double rMax,
                double *derivs[], int numWalkers);


#ifdef QMC_COMPLEX
void
two_body_derivs(float *R[], std::complex<float> *gradLogPsi[], int e1_first, int e1_last,
                int e2_first, int e2_last,
                int numCoefs, float rMax,
                float *derivs[], int numWalkers);
void
two_body_derivs(double *R[], std::complex<double> *gradLogPsi[], int e1_first, int e1_last,
                int e2_first, int e2_last,
                int numCoefs, double rMax,
                double *derivs[], int numWalkers);
#endif

///////////////////////
// One-Body routines //
///////////////////////

void
one_body_sum (float C[], float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
              float spline_coefs[], int numCoefs, float rMax,
              float sum[], int numWalkers);

void
one_body_sum (double C[], double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
              double spline_coefs[], int numCoefs, double rMax,
              double sum[], int numWalkers);

void
one_body_ratio (float C[], float *R[], int first, int last,
                float Rnew[], int inew,
                float spline_coefs[], int numCoefs, float rMax,
                float sum[], int numWalkers);

void
one_body_ratio (double C[], double *R[], int first, int last,
                double Rnew[], int inew,
                double spline_coefs[], int numCoefs, double rMax,
                double sum[], int numWalkers);

void
one_body_ratio_grad (float C[], float *R[], int first, int last,
                     float Rnew[], int inew,
                     float spline_coefs[], int numCoefs, float rMax,
                     bool zero, float ratio_grad[], int numWalkers);
void
one_body_ratio_grad (double C[], double *R[], int first, int last,
                     double Rnew[], int inew,
                     double spline_coefs[], int numCoefs, double rMax,
                     bool zero, double ratio_grad[], int numWalkers);

void
one_body_NLratios(NLjobGPU<float> jobs[], float C[], int first, int last,
                  float spline_coefs[], int numCoefs, float rMax, int numjobs);

void
one_body_NLratios(NLjobGPU<float> jobs[], float C[], int first, int last,
                  float spline_coefs[], int numCoefs, float rMax,
                  int numjobs);


void
one_body_NLratios(NLjobGPU<double> jobs[], double C[], int first, int last,
                  double spline_coefs[], int numCoefs, double rMax, int numjobs);

void
one_body_update(float *R[], int N, int iat, int numWalkers);

void
one_body_update(double *R[], int N, int iat, int numWalkers);


void
one_body_grad_lapl(float C[], float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                   float spline_coefs[], int numCoefs, float rMax,
                   float gradLapl[], int row_stride, int numWalkers);

void
one_body_grad_lapl(double C[], double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                   double spline_coefs[], int numCoefs, double rMax,
                   double gradLapl[], int row_stride, int numWalkers);

void
one_body_gradient (float *Rlist[], int iat, float C[], int first, int last,
                   float spline_coefs[], int num_coefs, float rMax,
                   bool zeroSum, float grad[], int numWalkers);

void
one_body_gradient (double *Rlist[], int iat, double C[], int first, int last,
                   double spline_coefs[], int num_coefs, double rMax,
                   bool zeroSum, double grad[], int numWalkers);


void
one_body_derivs(float C[], float *R[], float *gradLogPsi[],
                int cfirst, int clast,
                int efirst, int elast,
                int numCoefs, float rMax,
                float *derivs[], int numWalkers);

void
one_body_derivs(double C[], double *R[], double *gradLogPsi[],
                int cfirst, int clast,
                int efirst, int elast,
                int numCoefs, double rMax,
                double *derivs[], int numWalkers);


#ifdef QMC_COMPLEX
void
one_body_derivs(float C[], float *R[], std::complex<float> *gradLogPsi[],
                int cfirst, int clast,
                int efirst, int elast,
                int numCoefs, float rMax,
                float *derivs[], int numWalkers);

void
one_body_derivs(double C[], double *R[], std::complex<double> *gradLogPsi[],
                int cfirst, int clast,
                int efirst, int elast,
                int numCoefs, double rMax,
                double *derivs[], int numWalkers);

#endif

#endif
