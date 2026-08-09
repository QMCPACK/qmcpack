/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_COPY_H
#define MULTI_BSPLINE_COPY_H

#include "bspline_base.h"
#include "bspline_structs.h"
#include "multi_bspline_structs.h"

#ifdef __cplusplus
extern "C"
{
#endif

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////              Spline copy functions                 ////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  void copy_UBspline_1d_d(multi_UBspline_1d_d* multi,
                          int i,
                          const UBspline_1d_d* single,
                          const int offset,
                          const int N);

  void copy_UBspline_1d_d_s(multi_UBspline_1d_s* multi,
                            int i,
                            const UBspline_1d_d* single,
                            const int offset,
                            const int N);

#ifdef __cplusplus
}
#endif

#endif
