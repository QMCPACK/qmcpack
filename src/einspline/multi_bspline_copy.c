/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <stdio.h>
#include "multi_bspline_create.h"
#include "multi_bspline_copy.h"
#include "bspline_create.h"

#ifdef __cplusplus
extern "C" {
#endif



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
