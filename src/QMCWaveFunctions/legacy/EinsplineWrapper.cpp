//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "Configuration.h"
#include "QMCWaveFunctions/EinsplineWrapper.h"

namespace qmcplusplus
{

#if OHMMS_DIM==2
////////////////////////////////////////////////////////////
//                         2D                             //
////////////////////////////////////////////////////////////

//////////////////////
// Uniform routines //
//////////////////////

//////////
// Real //
//////////
UBsplineClass_2d_d::UBsplineClass_2d_d(TinyVector<double,2> start,
                                       TinyVector<double,2> end,
                                       BCtype_d xBC, BCtype_d yBC,
                                       Array<double,2> &data) :
  Spline(NULL)
{
  if (Spline != NULL)
    destroy_Bspline (Spline);
  Ugrid xgrid, ygrid;
  xgrid.start = start[0];
  xgrid.end = end[0];
  xgrid.num =data.size(0);
  ygrid.start = start[1];
  ygrid.end = end[1];
  ygrid.num =data.size(1);
  Spline = create_UBspline_2d_d (xgrid, ygrid, xBC, yBC, data.data());
}

double
UBsplineClass_2d_d::operator()(TinyVector<double,2> r)
{
  double val;
  eval_UBspline_2d_d (Spline, r[0], r[1], &val);
  return val;
}


void
UBsplineClass_2d_d::evaluate (TinyVector<double,2> r,
                              double &val, TinyVector<double,2> &grad)
{
  eval_UBspline_2d_d_vg (Spline, r[0], r[1], &val, &grad[0]);
}


void
UBsplineClass_2d_d::evaluate (TinyVector<double,2> r, double &val,
                              TinyVector<double,2> &grad, double &lapl)
{
  eval_UBspline_2d_d_vgl (Spline, r[0], r[1], &val,
                          &grad[0], &lapl);
}

void
UBsplineClass_2d_d::evaluate (TinyVector<double,2> r, double &val,
                              TinyVector<double,2> &grad,
                              Tensor<double,2> &hess)
{
  eval_UBspline_2d_d_vgh (Spline, r[0], r[1], &val, &grad[0],
                          &hess(0,0));
}

/////////////
// Complex //
/////////////
UBsplineClass_2d_z::UBsplineClass_2d_z(TinyVector<double,2> start,
                                       TinyVector<double,2> end,
                                       BCtype_z xBC, BCtype_z yBC,
                                       Array<std::complex<double>,2> &data) :
  Spline(NULL)
{
  if (Spline != NULL)
    destroy_Bspline (Spline);
  Ugrid xgrid, ygrid;
  xgrid.start = start[0];
  xgrid.end = end[0];
  xgrid.num =data.size(0);
  ygrid.start = start[1];
  ygrid.end = end[1];
  ygrid.num =data.size(1);
  Spline = create_UBspline_2d_z (xgrid, ygrid, xBC, yBC, data.data());
}

complex<double>
UBsplineClass_2d_z::operator()(TinyVector<double,2> r)
{
  std::complex<double> val;
  eval_UBspline_2d_z (Spline, r[0], r[1], &val);
  return val;
}


void
UBsplineClass_2d_z::evaluate (TinyVector<double,2> r,
                              std::complex<double> &val,
                              TinyVector<std::complex<double>,2> &grad)
{
  eval_UBspline_2d_z_vg (Spline, r[0], r[1], &val, &grad[0]);
}


void
UBsplineClass_2d_z::evaluate (TinyVector<double,2> r, std::complex<double> &val,
                              TinyVector<std::complex<double>,2> &grad,
                              std::complex<double> &lapl)
{
  eval_UBspline_2d_z_vgl (Spline, r[0], r[1], &val,
                          &grad[0], &lapl);
}

void
UBsplineClass_2d_z::evaluate (TinyVector<double,2> r, std::complex<double> &val,
                              TinyVector<std::complex<double>,2> &grad,
                              Tensor<std::complex<double>,2> &hess)
{
  eval_UBspline_2d_z_vgh (Spline, r[0], r[1], &val, &grad[0],
                          &hess(0,0));
}


/////////////////////////
// Nonuniform routines //
/////////////////////////

//////////
// Real //
//////////
NUBsplineClass_2d_d::NUBsplineClass_2d_d(TinyVector<double,2> start,
    TinyVector<double,2> end,
    TinyVector<double,2> ratio,
    BCtype_d xBC, BCtype_d yBC,
    Array<double,2> &data) :
  Spline(NULL)
{
  xGrid = create_center_grid (start[0], end[0], ratio[0], data.size(0));
  yGrid = create_center_grid (start[1], end[1], ratio[1], data.size(1));
  Spline = create_NUBspline_2d_d (xGrid, yGrid, xBC, yBC, data.data());
}

double
NUBsplineClass_2d_d::operator()(TinyVector<double,2> r)
{
  double val;
  eval_NUBspline_2d_d (Spline, r[0], r[1], &val);
  return val;
}

void
NUBsplineClass_2d_d::evaluate (TinyVector<double,2> r,
                               double &val, TinyVector<double,2> &grad)
{
  eval_NUBspline_2d_d_vg (Spline, r[0], r[1], &val, &grad[0]);
}


void
NUBsplineClass_2d_d::evaluate (TinyVector<double,2> r, double &val,
                               TinyVector<double,2> &grad, double &lapl)
{
  eval_NUBspline_2d_d_vgl (Spline, r[0], r[1], &val,
                           &grad[0], &lapl);
}

void
NUBsplineClass_2d_d::evaluate (TinyVector<double,2> r, double &val,
                               TinyVector<double,2> &grad,
                               Tensor<double,2> &hess)
{
  eval_NUBspline_2d_d_vgh (Spline, r[0], r[1], &val, &grad[0],
                           &hess(0,0));
}

/////////////
// Complex //
/////////////
NUBsplineClass_2d_z::NUBsplineClass_2d_z(TinyVector<double,2> start,
    TinyVector<double,2> end,
    TinyVector<double,2> ratio,
    BCtype_z xBC, BCtype_z yBC,
    Array<std::complex<double>,2> &data)
  : Spline(NULL)
{
  xGrid = create_center_grid (start[0], end[0], ratio[0], data.size(0));
  yGrid = create_center_grid (start[1], end[1], ratio[1], data.size(1));
  Spline = create_NUBspline_2d_z (xGrid, yGrid, xBC, yBC, data.data());
}

complex<double>
NUBsplineClass_2d_z::operator()(TinyVector<double,2> r)
{
  std::complex<double> val;
  eval_NUBspline_2d_z (Spline, r[0], r[1], &val);
  return val;
}

void
NUBsplineClass_2d_z::evaluate (TinyVector<double,2> r, std::complex<double> &val,
                               TinyVector<std::complex<double>,2> &grad)
{
  eval_NUBspline_2d_z_vg (Spline, r[0], r[1], &val, &grad[0]);
}

void
NUBsplineClass_2d_z::evaluate (TinyVector<double,2> r, std::complex<double> &val,
                               TinyVector<std::complex<double>,2> &grad,
                               std::complex<double> &lapl)
{
  eval_NUBspline_2d_z_vgl (Spline, r[0], r[1], &val,
                           &grad[0], &lapl);
}

void
NUBsplineClass_2d_z::evaluate (TinyVector<double,2> r, std::complex<double> &val,
                               TinyVector<std::complex<double>,2> &grad,
                               Tensor<std::complex<double>,2> &hess)
{
  eval_NUBspline_2d_z_vgh (Spline, r[0], r[1], &val, &grad[0],
                           &hess(0,0));
}

#elif OHMMS_DIM==3
#if !defined(__xlC__)
////////////////////////////////////////////////////////////
//                         3D                             //
////////////////////////////////////////////////////////////

//////////////////////
// Uniform routines //
//////////////////////

//////////
// Real //
//////////
UBsplineClass_3d_d::UBsplineClass_3d_d(TinyVector<double,3> start,
                                       TinyVector<double,3> end,
                                       BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
                                       Array<double,3> &data) :
  Spline(NULL)
{
  if (Spline != NULL)
    destroy_Bspline (Spline);
  Ugrid xgrid, ygrid, zgrid;
  xgrid.start = start[0];
  xgrid.end = end[0];
  xgrid.num =data.size(0);
  ygrid.start = start[1];
  ygrid.end = end[1];
  ygrid.num =data.size(1);
  zgrid.start = start[2];
  zgrid.end = end[2];
  zgrid.num =data.size(2);
  Spline = create_UBspline_3d_d (xgrid, ygrid, zgrid, xBC, yBC, zBC, data.data());
}


double
UBsplineClass_3d_d::operator()(TinyVector<double,3> r)
{
  double val;
  eval_UBspline_3d_d (Spline, r[0], r[1], r[2], &val);
  return val;
}


void
UBsplineClass_3d_d::evaluate (TinyVector<double,3> r,
                              double &val, TinyVector<double,3> &grad)
{
  eval_UBspline_3d_d_vg (Spline, r[0], r[1], r[2], &val, &grad[0]);
}


void
UBsplineClass_3d_d::evaluate (TinyVector<double,3> r, double &val,
                              TinyVector<double,3> &grad, double &lapl)
{
  eval_UBspline_3d_d_vgl (Spline, r[0], r[1], r[2], &val,
                          &grad[0], &lapl);
}

void
UBsplineClass_3d_d::evaluate (TinyVector<double,3> r, double &val,
                              TinyVector<double,3> &grad,
                              Tensor<double,3> &hess)
{
  eval_UBspline_3d_d_vgh (Spline, r[0], r[1], r[2], &val, &grad[0],
                          &hess(0,0));
}

/////////////
// Complex //
/////////////
UBsplineClass_3d_z::UBsplineClass_3d_z(TinyVector<double,3> start,
                                       TinyVector<double,3> end,
                                       BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
                                       Array<std::complex<double>,3> &data) :
  Spline(NULL)
{
  if (Spline != NULL)
    destroy_Bspline (Spline);
  Ugrid xgrid, ygrid, zgrid;
  xgrid.start = start[0];
  xgrid.end = end[0];
  xgrid.num =data.size(0);
  ygrid.start = start[1];
  ygrid.end = end[1];
  ygrid.num =data.size(1);
  zgrid.start = start[2];
  zgrid.end = end[2];
  zgrid.num =data.size(2);
  Spline = create_UBspline_3d_z (xgrid, ygrid, zgrid, xBC, yBC, zBC, data.data());
}

complex<double>
UBsplineClass_3d_z::operator()(TinyVector<double,3> r)
{
  std::complex<double> val;
  eval_UBspline_3d_z (Spline, r[0], r[1], r[2], &val);
  return val;
}


void
UBsplineClass_3d_z::evaluate (TinyVector<double,3> r,
                              std::complex<double> &val,
                              TinyVector<std::complex<double>,3> &grad)
{
  eval_UBspline_3d_z_vg (Spline, r[0], r[1], r[2], &val, &grad[0]);
}


void
UBsplineClass_3d_z::evaluate (TinyVector<double,3> r, std::complex<double> &val,
                              TinyVector<std::complex<double>,3> &grad,
                              std::complex<double> &lapl)
{
  eval_UBspline_3d_z_vgl (Spline, r[0], r[1], r[2], &val,
                          &grad[0], &lapl);
}

void
UBsplineClass_3d_z::evaluate (TinyVector<double,3> r, std::complex<double> &val,
                              TinyVector<std::complex<double>,3> &grad,
                              Tensor<std::complex<double>,3> &hess)
{
  eval_UBspline_3d_z_vgh (Spline, r[0], r[1], r[2], &val, &grad[0],
                          &hess(0,0));
}


////////////////////////////////////
// Nonuniform evaluation routines //
////////////////////////////////////

//////////
// Real //
//////////
NUBsplineClass_3d_d::NUBsplineClass_3d_d(TinyVector<double,3> start,
    TinyVector<double,3> end,
    TinyVector<double,3> ratio,
    BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
    Array<double,3> &data)
{
  xGrid = create_center_grid (start[0], end[0], ratio[0], data.size(0));
  yGrid = create_center_grid (start[1], end[1], ratio[1], data.size(1));
  zGrid = create_center_grid (start[2], end[2], ratio[2], data.size(2));
  Spline = create_NUBspline_3d_d (xGrid, yGrid, zGrid,
                                  xBC, yBC, zBC, data.data());
}

double
NUBsplineClass_3d_d::operator()(TinyVector<double,3> r)
{
  double val;
  eval_NUBspline_3d_d (Spline, r[0], r[1], r[2], &val);
  return val;
}


void
NUBsplineClass_3d_d::evaluate (TinyVector<double,3> r,
                               double &val, TinyVector<double,3> &grad)
{
  eval_NUBspline_3d_d_vg (Spline, r[0], r[1], r[2], &val, &grad[0]);
}


void
NUBsplineClass_3d_d::evaluate (TinyVector<double,3> r, double &val,
                               TinyVector<double,3> &grad, double &lapl)
{
  eval_NUBspline_3d_d_vgl (Spline, r[0], r[1], r[2], &val,
                           &grad[0], &lapl);
}

void
NUBsplineClass_3d_d::evaluate (TinyVector<double,3> r, double &val,
                               TinyVector<double,3> &grad,
                               Tensor<double,3> &hess)
{
  eval_NUBspline_3d_d_vgh (Spline, r[0], r[1], r[2], &val, &grad[0],
                           &hess(0,0));
}

/////////////
// Complex //
/////////////
NUBsplineClass_3d_z::NUBsplineClass_3d_z(TinyVector<double,3> start,
    TinyVector<double,3> end,
    TinyVector<double,3> ratio,
    BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
    Array<std::complex<double>,3> &data)
{
  xGrid = create_center_grid (start[0], end[0], ratio[0], data.size(0));
  yGrid = create_center_grid (start[1], end[1], ratio[1], data.size(1));
  zGrid = create_center_grid (start[2], end[2], ratio[2], data.size(2));
  Spline = create_NUBspline_3d_z (xGrid, yGrid, zGrid,
                                  xBC, yBC, zBC, data.data());
}


complex<double>
NUBsplineClass_3d_z::operator()(TinyVector<double,3> r)
{
  std::complex<double> val;
  eval_NUBspline_3d_z (Spline, r[0], r[1], r[2], &val);
  return val;
}

void
NUBsplineClass_3d_z::evaluate (TinyVector<double,3> r, std::complex<double> &val,
                               TinyVector<std::complex<double>,3> &grad)
{
  eval_NUBspline_3d_z_vg (Spline, r[0], r[1], r[2], &val, &grad[0]);
}


void
NUBsplineClass_3d_z::evaluate (TinyVector<double,3> r, std::complex<double> &val,
                               TinyVector<std::complex<double>,3> &grad,
                               std::complex<double> &lapl)
{
  eval_NUBspline_3d_z_vgl (Spline, r[0], r[1], r[2], &val,
                           &grad[0], &lapl);
}

void
NUBsplineClass_3d_z::evaluate (TinyVector<double,3> r, std::complex<double> &val,
                               TinyVector<std::complex<double>,3> &grad,
                               Tensor<std::complex<double>,3> &hess)
{
  eval_NUBspline_3d_z_vgh (Spline, r[0], r[1], r[2], &val, &grad[0],
                           &hess(0,0));
}
#endif
#endif
}
