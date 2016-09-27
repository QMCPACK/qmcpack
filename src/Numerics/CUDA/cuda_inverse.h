//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//		      Jeremy McMinnis, jmcminis@gmail.com, Navar Inc. 
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//
// File created by:  Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef CUDA_INVERSE_H
#define CUDA_INVERSE_H

#include <cublas_v2.h>

//////////////////////////////
// Single / mixed precision //
//////////////////////////////
void 
cublas_inverse (cublasHandle_t handle,
                float *Alist_d[], float *Ainvlist_d[],
                float *AWorkList_d[], float *AinvWorkList_d[],
                int N, int rowStride, int numMats,
                bool useHigherPrecision = true);

//////////////////////
// Double precision //
//////////////////////
void
cublas_inverse (cublasHandle_t handle,
                double *Alist_d[], double *Ainvlist_d[],
                double *AWorklist_d[], double *AinvWorklist_d[],
                int N, int rowStride, int numMats, 
                bool useHigherPrecision = true);

#endif
