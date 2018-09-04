//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


///\file CUDABsplineConversions.h 
#ifndef QMCPLUSPLUS_CUDABSPLINECONVERSIONS_H
#define QMCPLUSPLUS_CUDABSPLINECONVERSIONS_H

// Cuda vector types
#include <vector_types.h>

#include "qmc_common.h"
#include "einspline/bspline_base_cuda.h"
#include "einspline/multi_bspline.h"
#include "einspline/multi_bspline_create_cuda.h"

namespace qmcplusplus
{
namespace cudasoatemp
{

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_d *in,
    multi_UBspline_3d_s_cuda* &out)
{
  out = create_multi_UBspline_3d_s_cuda_conv (in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_d *in,
    multi_UBspline_3d_d_cuda * &out)
{
  out = create_multi_UBspline_3d_d_cuda(in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_c_cuda* &out)
{
  out = create_multi_UBspline_3d_c_cuda_conv (in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_z_cuda * &out)
{
  out = create_multi_UBspline_3d_z_cuda(in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_d_cuda * &out)
{

  APP_ABORT("Attempted to convert complex CPU spline into a real GPU spline.\n");
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_s_cuda * &out)
{
  APP_ABORT("Attempted to convert complex CPU spline into a real GPU spline.\n");
}

}
}

#endif //QMCPLUSPLUS_CUDABSPLINECONVERSIONS_H
