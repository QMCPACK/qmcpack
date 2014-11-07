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
