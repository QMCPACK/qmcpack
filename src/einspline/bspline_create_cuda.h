#ifndef BSPLINE_CREATE_CUDA_H
#define BSPLINE_CREATE_CUDA_H

#include "bspline_structs_cuda.h"

extern "C" UBspline_3d_s_cuda*
create_UBspline_3d_s_cuda (UBspline_3d_s* spline);

extern "C" UBspline_3d_s_cuda*
create_UBspline_3d_s_cuda_conv (UBspline_3d_d* spline);


extern "C" UBspline_3d_c_cuda*
create_UBspline_3d_c_cuda (UBspline_3d_c* spline);

extern "C" UBspline_3d_c_cuda*
create_UBspline_3d_c_cuda_conv (UBspline_3d_z* spline);


extern "C" UBspline_3d_d_cuda*
create_UBspline_3d_d_cuda (UBspline_3d_d* spline);

extern "C" UBspline_3d_z_cuda*
create_UBspline_3d_z_cuda (UBspline_3d_z* spline);

#endif
