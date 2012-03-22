#ifndef MULTI_BSPLINE_STRUCTS_CUDA_H
#define MULTI_BSPLINE_STRUCTS_CUDA_H

#define SPLINE_BLOCK_SIZE 64

#include "bspline_base.h"
#include "bspline_base_cuda.h"



////////
// 1D //
////////
typedef struct
{
  float *coefs;
  uint stride;
  float gridInv;
  int num_splines;
  int dim;
} multi_UBspline_1d_s_cuda;

typedef struct
{
  complex_float *coefs;
  uint stride;
  float gridInv;
  int num_splines;
  int dim;
} multi_UBspline_1d_c_cuda;

typedef struct
{
  double *coefs;
  uint stride;
  double gridInv;
  int num_splines;
  int dim;
} multi_UBspline_1d_d_cuda;

typedef struct
{
  complex_double *coefs;
  uint stride;
  double gridInv;
  int num_splines;
  int dim;
} multi_UBspline_1d_z_cuda;

////////
// 2D //
////////

typedef struct
{
  float *coefs;
  uint2 stride;
  float2 gridInv;
  int num_splines;
} multi_UBspline_2d_s_cuda;

typedef struct
{
  complex_float *coefs;
  uint2 stride;
  float2 gridInv;
  int num_splines;
} multi_UBspline_2d_c_cuda;

typedef struct
{
  double *coefs;
  uint2 stride;
  double gridInv[2];
  int num_splines;
} multi_UBspline_2d_d_cuda;

typedef struct
{
  complex_double *coefs;
  uint2 stride;
  double gridInv[2];
  int num_splines;
} multi_UBspline_2d_z_cuda;

////////
// 3D //
////////

typedef struct
{
  float *coefs;
  uint3 stride;
  float3 gridInv;
  uint3 dim;
  int num_splines;
} multi_UBspline_3d_s_cuda;

typedef struct
{
  complex_float *coefs;
  uint3 stride;
  float3 gridInv;
  uint3 dim;
  int num_splines;
} multi_UBspline_3d_c_cuda;

typedef struct
{
  double *coefs;
  uint3 stride;
  double3 gridInv;
  uint3 dim;
  int num_splines;
} multi_UBspline_3d_d_cuda;

typedef struct
{
  complex_double *coefs;
  uint3 stride;
  double3 gridInv;
  uint3 dim;
  int num_splines;
} multi_UBspline_3d_z_cuda;



#endif
