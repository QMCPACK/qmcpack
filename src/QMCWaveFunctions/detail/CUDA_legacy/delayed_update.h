#ifndef DELAYED_UPDATE_H
#define DELAYED_UPDATE_H

#include "config.h"
#ifndef QMC_CUDA2HIP
#include <complex>
#include <cublas_v2.h>
#else
#include <hipblas.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif

//////////////////////////////
// Single / mixed precision //
//////////////////////////////

void cublas_lemma_mats(cublasHandle_t handle,
                       float* AinvList_d[],
                       float* U_d[],
                       float* lemma_d[],
                       float* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride);

void cublas_ainv_row(cublasHandle_t handle,
                     float* AinvkList_d[],
                     float* AWorkList_d[],
                     float* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride);

void cublas_smw_update(cublasHandle_t handle,
                       float* AinvkList_d[],
                       float* AinvList_d[],
                       float* AinvUList_d[],
                       float* AWorkList_d[],
                       float* lemma_inv[],
                       float* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride);

void update_onemove(float* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num);

void multi_row_copy(float* dest[], float* src[], int len, int offset, int rows, int stride, int num);

void calc_lemma_column(float* ainv[],
                       float* newrow[],
                       float* lemma[],
                       float* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num);

void copy_update_block(float* lemma_lu[],
                       float* lemma[],
                       float* ainv_work[],
                       float* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num);

void copy_delayed(float* lemma_lu[],
                  float* lemma[],
                  float* ainv_row[],
                  float* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num);

void calc_gradlapl_and_collect(float* lemma_lu[],
                               float* Ainv_row[],
                               float* GL_col[],
                               float ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num);

void calc_gradient_delayed(float* Ainv_row[], float* GL_col[], float ratios[], int N, int rowstride, int num);

//////////////////////
// Double precision //
//////////////////////

void cublas_lemma_mats(cublasHandle_t handle,
                       double* AinvList_d[],
                       double* U_d[],
                       double* lemma_d[],
                       double* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride);

void cublas_ainv_row(cublasHandle_t handle,
                     double* AinvkList_d[],
                     double* AWorkList_d[],
                     double* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride);

void cublas_smw_update(cublasHandle_t handle,
                       double* AinvkList_d[],
                       double* AinvList_d[],
                       double* AinvUList_d[],
                       double* AWorkList_d[],
                       double* lemma_inv[],
                       double* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride);

void update_onemove(double* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num);

void multi_row_copy(double* dest[], double* src[], int len, int offset, int rows, int stride, int num);

void calc_lemma_column(double* ainv[],
                       double* newrow[],
                       double* lemma[],
                       double* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num);

void copy_update_block(double* lemma_lu[],
                       double* lemma[],
                       double* ainv_work[],
                       double* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num);

void copy_delayed(double* lemma_lu[],
                  double* lemma[],
                  double* ainv_row[],
                  double* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num);

void calc_gradlapl_and_collect(double* lemma_lu[],
                               double* Ainv_row[],
                               double* GL_col[],
                               double ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num);

void calc_gradient_delayed(double* Ainv_row[], double* GL_col[], double ratios[], int N, int rowstride, int num);


#ifdef QMC_COMPLEX

//////////////////////////////////////
// Complex single / mixed precision //
//////////////////////////////////////

void cublas_lemma_mats(cublasHandle_t handle,
                       std::complex<float>* AinvList_d[],
                       std::complex<float>* U_d[],
                       std::complex<float>* lemma_d[],
                       std::complex<float>* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride);

void cublas_ainv_row(cublasHandle_t handle,
                     std::complex<float>* AinvkList_d[],
                     std::complex<float>* AWorkList_d[],
                     std::complex<float>* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride);

void cublas_smw_update(cublasHandle_t handle,
                       std::complex<float>* AinvkList_d[],
                       std::complex<float>* AinvList_d[],
                       std::complex<float>* AinvUList_d[],
                       std::complex<float>* AWorkList_d[],
                       std::complex<float>* lemma_inv[],
                       std::complex<float>* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride);

void update_onemove(std::complex<float>* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num);

void multi_row_copy(std::complex<float>* dest[],
                    std::complex<float>* src[],
                    int len,
                    int offset,
                    int rows,
                    int stride,
                    int num);

void calc_lemma_column(std::complex<float>* ainv[],
                       std::complex<float>* newrow[],
                       std::complex<float>* lemma[],
                       std::complex<float>* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num);

void copy_update_block(std::complex<float>* lemma_lu[],
                       std::complex<float>* lemma[],
                       std::complex<float>* ainv_work[],
                       std::complex<float>* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num);

void copy_delayed(std::complex<float>* lemma_lu[],
                  std::complex<float>* lemma[],
                  std::complex<float>* ainv_row[],
                  std::complex<float>* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num);

void calc_gradlapl_and_collect(std::complex<float>* lemma_lu[],
                               std::complex<float>* Ainv_row[],
                               std::complex<float>* GL_col[],
                               std::complex<float> ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num);

void calc_gradient_delayed(std::complex<float>* Ainv_row[],
                           std::complex<float>* GL_col[],
                           std::complex<float> ratios[],
                           int N,
                           int rowstride,
                           int num);

//////////////////////////////
// Complex double precision //
//////////////////////////////

void cublas_lemma_mats(cublasHandle_t handle,
                       std::complex<double>* AinvList_d[],
                       std::complex<double>* U_d[],
                       std::complex<double>* lemma_d[],
                       std::complex<double>* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride);

void cublas_ainv_row(cublasHandle_t handle,
                     std::complex<double>* AinvkList_d[],
                     std::complex<double>* AWorkList_d[],
                     std::complex<double>* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride);

void cublas_smw_update(cublasHandle_t handle,
                       std::complex<double>* AinvkList_d[],
                       std::complex<double>* AinvList_d[],
                       std::complex<double>* AinvUList_d[],
                       std::complex<double>* AWorkList_d[],
                       std::complex<double>* lemma_inv[],
                       std::complex<double>* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride);

void update_onemove(std::complex<double>* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num);

void multi_row_copy(std::complex<double>* dest[],
                    std::complex<double>* src[],
                    int len,
                    int offset,
                    int rows,
                    int stride,
                    int num);

void calc_lemma_column(std::complex<double>* ainv[],
                       std::complex<double>* newrow[],
                       std::complex<double>* lemma[],
                       std::complex<double>* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num);

void copy_update_block(std::complex<double>* lemma_lu[],
                       std::complex<double>* lemma[],
                       std::complex<double>* ainv_work[],
                       std::complex<double>* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num);

void copy_delayed(std::complex<double>* lemma_lu[],
                  std::complex<double>* lemma[],
                  std::complex<double>* ainv_row[],
                  std::complex<double>* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num);

void calc_gradlapl_and_collect(std::complex<double>* lemma_lu[],
                               std::complex<double>* Ainv_row[],
                               std::complex<double>* GL_col[],
                               std::complex<double> ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num);

void calc_gradient_delayed(std::complex<double>* Ainv_row[],
                           std::complex<double>* GL_col[],
                           std::complex<double> ratios[],
                           int N,
                           int rowstride,
                           int num);

#endif // QMC_COMPLEX

#endif
