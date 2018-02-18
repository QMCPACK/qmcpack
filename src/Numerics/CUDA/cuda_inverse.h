//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//		              Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//
// File created by:  Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef CUDA_INVERSE_H
#define CUDA_INVERSE_H

#include <complex>
#include <cublas_v2.h>

//////////////////////////////
// Single / mixed precision //
//////////////////////////////
void 
cublas_inverse (cublasHandle_t handle,
                float *Alist_d[], float *Ainvlist_d[],
                float *AWorkList_d[], float *AinvWorkList_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats,
                bool useHigherPrecision = true);

//////////////////////
// Double precision //
//////////////////////
void
cublas_inverse (cublasHandle_t handle,
                double *Alist_d[], double *Ainvlist_d[],
                double *AWorklist_d[], double *AinvWorklist_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats, 
                bool useHigherPrecision = true);


//////////////////////////////////////
// Complex single / mixed precision //
//////////////////////////////////////
void 
cublas_inverse (cublasHandle_t handle,
                std::complex<float> *Alist_d[], std::complex<float> *Ainvlist_d[],
                std::complex<float> *AWorkList_d[], std::complex<float> *AinvWorkList_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats,
                bool useHigherPrecision = true);

//////////////////////////////
// Complex double precision //
//////////////////////////////
void
cublas_inverse (cublasHandle_t handle,
                std::complex<double> *Alist_d[], std::complex<double> *Ainvlist_d[],
                std::complex<double> *AWorklist_d[], std::complex<double> *AinvWorklist_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats, 
                bool useHigherPrecision = true);

#endif
