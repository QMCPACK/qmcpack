#ifndef MULTI_BSPLINE_CREATE_CUDA_H
#define MULTI_BSPLINE_CREATE_CUDA_H

#include "multi_bspline_structs_cuda.h"


////////
// 1D //
////////
extern "C" multi_UBspline_1d_s_cuda*
create_multi_UBspline_1d_s_cuda (multi_UBspline_1d_s* spline);

extern "C" multi_UBspline_1d_s_cuda*
create_multi_UBspline_1d_s_cuda_conv (multi_UBspline_1d_d* spline);


extern "C" multi_UBspline_1d_c_cuda*
create_multi_UBspline_1d_c_cuda (multi_UBspline_1d_c* spline);

extern "C" multi_UBspline_1d_c_cuda*
create_multi_UBspline_1d_c_cuda_conv (multi_UBspline_1d_z* spline);


extern "C" multi_UBspline_1d_d_cuda*
create_multi_UBspline_1d_d_cuda (multi_UBspline_1d_d* spline);

extern "C" multi_UBspline_1d_z_cuda*
create_multi_UBspline_1d_z_cuda (multi_UBspline_1d_z* spline);

////////
// 2D //
////////
extern "C" multi_UBspline_2d_s_cuda*
create_multi_UBspline_2d_s_cuda (multi_UBspline_2d_s* spline);

extern "C" multi_UBspline_2d_s_cuda*
create_multi_UBspline_2d_s_cuda_conv (multi_UBspline_2d_d* spline);


extern "C" multi_UBspline_2d_c_cuda*
create_multi_UBspline_2d_c_cuda (multi_UBspline_2d_c* spline);

extern "C" multi_UBspline_2d_c_cuda*
create_multi_UBspline_2d_c_cuda_conv (multi_UBspline_2d_z* spline);


extern "C" multi_UBspline_2d_d_cuda*
create_multi_UBspline_2d_d_cuda (multi_UBspline_2d_d* spline);

extern "C" multi_UBspline_2d_z_cuda*
create_multi_UBspline_2d_z_cuda (multi_UBspline_2d_z* spline);




////////
// 3D //
////////

extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda (multi_UBspline_3d_s* spline);

extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda_conv (multi_UBspline_3d_d* spline);


extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda (multi_UBspline_3d_c* spline);

extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda_conv (multi_UBspline_3d_z* spline);


extern "C" multi_UBspline_3d_d_cuda*
create_multi_UBspline_3d_d_cuda (multi_UBspline_3d_d* spline);

extern "C" multi_UBspline_3d_z_cuda*
create_multi_UBspline_3d_z_cuda (multi_UBspline_3d_z* spline);

#endif
