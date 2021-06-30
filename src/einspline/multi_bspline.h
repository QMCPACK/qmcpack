/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_H
#define MULTI_BSPLINE_H

#include "bspline_base.h"
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////           Bspline structure definitions            ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
#include "multi_bspline_structs.h"

// Currently, some of the single-precision routines use SSE2 instructions
#include "multi_bspline_eval_s.h"
// #include "multi_bspline_eval_c.h"
#include "multi_bspline_eval_d.h"
#include "multi_bspline_eval_z.h"

#include "bspline_create.h"
#include "multi_bspline_create.h"
#endif
