//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef BSPLINE_STRUCTS_CUDA_H
#define BSPLINE_STRUCTS_CUDA_H

#include "bspline_base_cuda.h"

#define SPLINE_BLOCK_SIZE 64

////////
// 2D //
////////

#if CUDA_VERSION < 3000
typedef struct
{
  double x,y,z;
} double3;

typedef struct
{
  double x,y,z,w;
} double4;
#endif

typedef struct
{
  float *coefs;
  uint2 stride;
  float2 gridInv;
} UBspline_2d_s_cuda;

typedef struct
{
  float *coefs_real, *coefs_imag;
  uint2 stride;
  float2 gridInv;
} UBspline_2d_c_cuda;

typedef struct
{
  double *coefs;
  uint2 stride;
  double gridInv[2];
} UBspline_2d_d_cuda;

typedef struct
{
  complex_double *coefs;
  uint2 stride;
  double gridInv[2];
} UBspline_2d_z_cuda;

////////
// 3D //
////////

typedef struct
{
  float *coefs;
  uint3 stride;
  float3 gridInv;
  uint3 dim;
} UBspline_3d_s_cuda;

typedef struct
{
  complex_float *coefs;
  uint3 stride;
  float3 gridInv;
  uint3 dim;
} UBspline_3d_c_cuda;

typedef struct
{
  double *coefs;
  uint3 stride;
  double3 gridInv;
  uint3 dim;
} UBspline_3d_d_cuda;

typedef struct
{
  complex_double *coefs;
  uint3 stride;
  double3 gridInv;
  uint3 dim;
} UBspline_3d_z_cuda;



#endif
