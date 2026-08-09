/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <stdio.h>
#include "multi_bspline_create.h"
#include "multi_bspline_copy.h"
#include "multi_bspline_eval_d.h"
#include "bspline_create.h"

#ifdef __cplusplus
extern "C" {
#endif

  void copy_UBspline_3d_d(multi_UBspline_3d_d* multi, int i
      , const UBspline_3d_d* single, const int* offset, const int* N)
  {
    intptr_t z_stride=multi->z_stride;

    for(int ix=0; ix<N[0]; ++ix)
      for(int iy=0; iy<N[1]; ++iy)
      {
        intptr_t out=ix*multi->x_stride+iy*multi->y_stride+i;
        intptr_t in =(ix+offset[0])*single->x_stride+(iy+offset[1])*single->y_stride+offset[2];
        for(int iz=0; iz<N[2]; ++iz, in++, out +=z_stride)
        {
          multi->coefs[out]=single->coefs[in];
        }
      }
  }

  void copy_UBspline_3d_d_s(multi_UBspline_3d_s* multi, int i
      , const UBspline_3d_d* single, const int* offset, const int* N)
  {
    intptr_t x_stride_in=single->x_stride;
    intptr_t y_stride_in=single->y_stride;
    intptr_t x_stride_out=multi->x_stride;
    intptr_t y_stride_out=multi->y_stride;
    intptr_t z_stride_out=multi->z_stride;
    intptr_t offset0=(intptr_t)offset[0];
    intptr_t offset1=(intptr_t)offset[1];
    intptr_t offset2=(intptr_t)offset[2];
    const intptr_t istart=(intptr_t)i;
    const intptr_t n0=N[0],n1=N[1],n2=N[2];
    for(intptr_t ix=0; ix<n0; ++ix)
      for(intptr_t iy=0; iy<n1; ++iy)
      {
        float* restrict out=multi->coefs+ix*x_stride_out+iy*y_stride_out+istart;
        const double* restrict in =single->coefs+(ix+offset0)*x_stride_in+(iy+offset1)*y_stride_in+offset2;
        for(intptr_t iz=0; iz<n2; ++iz)
        {
          out[iz*z_stride_out]=(float)in[iz];
        }
      }
//
//    for(int ix=0; ix<N[0]; ++ix)
//      for(int iy=0; iy<N[1]; ++iy)
//      {
//        intptr_t out=ix*multi->x_stride+iy*multi->y_stride+i;
//        intptr_t in =(ix+offset[0])*single->x_stride+(iy+offset[1])*single->y_stride+offset[2];
//        for(int iz=0; iz<N[2]; ++iz, in++, out +=z_stride)
//        {
//          multi->coefs[out]=(float)single->coefs[in];
//        }
//      }
  }

  // Create 3D uniform single-precision, real Bspline
  multi_UBspline_3d_s *
    copy_multi_UBspline_3d_s (multi_UBspline_3d_s* spline)
    {
      multi_UBspline_3d_s *clone=create_multi_UBspline_3d_s(
          spline->x_grid, spline->y_grid, spline->z_grid
          ,spline->xBC, spline->yBC, spline->zBC
          ,spline->num_splines);
      memcpy(clone->coefs,spline->coefs,sizeof(float)*clone->coefs_size);
      return clone;
    }

  // just checking not used

  void copy_UBspline_1d_d(multi_UBspline_1d_d* multi, int i
      , const UBspline_1d_d* single, const int offset, const int N)
  {
    //fprintf(stdout,"debug xstride %ld i %d N %d \n", multi->x_stride, i, offset);
    for(int ix=0; ix<N; ++ix)
    {
      intptr_t out=ix*multi->x_stride+i;
      intptr_t in =ix+offset;
      multi->coefs[out]=single->coefs[in];
    }
  }

  void copy_UBspline_1d_d_s(multi_UBspline_1d_s* multi, int i
      , const UBspline_1d_d* single, const int offset, const int N)
  {
    for(int ix=0; ix<N; ++ix)
    {
      intptr_t out=ix*multi->x_stride+i;
      intptr_t in =ix+offset;
      multi->coefs[out]=(float)single->coefs[in];
    }
  }

#ifdef __cplusplus
}
#endif
