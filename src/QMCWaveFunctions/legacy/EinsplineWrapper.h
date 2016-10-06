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
    
    
#ifndef EINSPLINE_WRAPPER_H
#define EINSPLINE_WRAPPER_H

#include <config.h>
#include <einspline/bspline.h>
#include <einspline/nubspline.h>
#include "OhmmsPETE/OhmmsArray.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/TinyVectorOps.h"

namespace qmcplusplus
{
#if   OHMMS_DIM==2
////////////////////////////////////////////////////////////
//                         2D                             //
////////////////////////////////////////////////////////////
class BsplineClass_2d_d
{
public:
  virtual double operator()(TinyVector<double,2> r)                 = 0;
  virtual void evaluate (TinyVector<double,2> r,
                         double &val, TinyVector<double,2> &grad)   = 0;
  virtual void evaluate (TinyVector<double,2> r, double &val,
                         TinyVector<double,2> &grad, double &lapl)  = 0;
  virtual void evaluate (TinyVector<double,2> r, double &val,
                         TinyVector<double,2> &grad,
                         Tensor<double,2> &hess)                    = 0;
};

class UBsplineClass_2d_d : public BsplineClass_2d_d
{
private:
  UBspline_2d_d *Spline;
public:
  double operator()(TinyVector<double,2> r);
  void evaluate (TinyVector<double,2> r,
                 double &val, TinyVector<double,2> &grad);
  void evaluate (TinyVector<double,2> r, double &val,
                 TinyVector<double,2> &grad, double &lapl);
  void evaluate (TinyVector<double,2> r, double &val,
                 TinyVector<double,2> &grad,
                 Tensor<double,2> &hess);
  UBsplineClass_2d_d (TinyVector<double,2> start,
                      TinyVector<double,2> end,
                      BCtype_d xBC, BCtype_d yBC,
                      Array<double,2> &data);
};

class NUBsplineClass_2d_d : public BsplineClass_2d_d
{
private:
  NUBspline_2d_d *Spline;
  NUgrid *xGrid, *yGrid, *zGrid;
public:
  double operator()(TinyVector<double,2> r);
  void evaluate (TinyVector<double,2> r,
                 double &val, TinyVector<double,2> &grad);
  void evaluate (TinyVector<double,2> r, double &val,
                 TinyVector<double,2> &grad, double &lapl);
  void evaluate (TinyVector<double,2> r, double &val,
                 TinyVector<double,2> &grad,
                 Tensor<double,2> &hess);
  NUBsplineClass_2d_d (TinyVector<double,2> start,
                       TinyVector<double,2> end,
                       TinyVector<double,2> ratio,
                       BCtype_d xBC, BCtype_d yBC,
                       Array<double,2> &data);
};


/////////////////////////////
// Complex Bspline classes //
/////////////////////////////
class BsplineClass_2d_z
{
public:
  virtual std::complex<double> operator()(TinyVector<double,2> r)                 = 0;
  virtual void evaluate (TinyVector<double,2> r,
                         std::complex<double> &val,
                         TinyVector<std::complex<double>,2> &grad)   = 0;
  virtual void evaluate (TinyVector<double,2> r, std::complex<double> &val,
                         TinyVector<std::complex<double>,2> &grad,
                         std::complex<double> &lapl)  = 0;
  virtual void evaluate (TinyVector<double,2> r, std::complex<double> &val,
                         TinyVector<std::complex<double>,2> &grad,
                         Tensor<std::complex<double>,2> &hess)                    = 0;
};

class UBsplineClass_2d_z : public BsplineClass_2d_z
{
private:
  UBspline_2d_z *Spline;
public:
  std::complex<double> operator()(TinyVector<double,2> r);
  void evaluate (TinyVector<double,2> r,
                 std::complex<double> &val, TinyVector<std::complex<double>,2> &grad);
  void evaluate (TinyVector<double,2> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,2> &grad,
                 std::complex<double> &lapl);
  void evaluate (TinyVector<double,2> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,2> &grad,
                 Tensor<std::complex<double>,2> &hess);
  UBsplineClass_2d_z (TinyVector<double,2> start,
                      TinyVector<double,2> end,
                      BCtype_z xBC, BCtype_z yBC,
                      Array<std::complex<double>,2> &data);
};

class NUBsplineClass_2d_z : public BsplineClass_2d_z
{
private:
  NUBspline_2d_z *Spline;
  NUgrid *xGrid, *yGrid, *zGrid;
public:
  std::complex<double> operator()(TinyVector<double,2> r);
  void evaluate (TinyVector<double,2> r,
                 std::complex<double> &val, TinyVector<std::complex<double>,2> &grad);
  void evaluate (TinyVector<double,2> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,2> &grad, std::complex<double> &lapl);
  void evaluate (TinyVector<double,2> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,2> &grad,
                 Tensor<std::complex<double>,2> &hess);
  NUBsplineClass_2d_z (TinyVector<double,2> start,
                       TinyVector<double,2> end,
                       TinyVector<double,2> ratio,
                       BCtype_z xBC, BCtype_z yBC,
                       Array<std::complex<double>,2> &data);
};

#elif   OHMMS_DIM==3

////////////////////////////////////////////////////////////
//                         3D                             //
////////////////////////////////////////////////////////////
class BsplineClass_3d_d
{

public:
  virtual double operator()(TinyVector<double,3> r)                 = 0;
  virtual void evaluate (TinyVector<double,3> r,
                         double &val, TinyVector<double,3> &grad)   = 0;
  virtual void evaluate (TinyVector<double,3> r, double &val,
                         TinyVector<double,3> &grad, double &lapl)  = 0;
  virtual void evaluate (TinyVector<double,3> r, double &val,
                         TinyVector<double,3> &grad,
                         Tensor<double,3> &hess)                    = 0;
};
class UBsplineClass_3d_d : public BsplineClass_3d_d
{
private:
  UBspline_3d_d *Spline;
public:
  double operator()(TinyVector<double,3> r);
  void evaluate (TinyVector<double,3> r,
                 double &val, TinyVector<double,3> &grad);
  void evaluate (TinyVector<double,3> r, double &val,
                 TinyVector<double,3> &grad, double &lapl);
  void evaluate (TinyVector<double,3> r, double &val,
                 TinyVector<double,3> &grad,
                 Tensor<double,3> &hess);
  UBsplineClass_3d_d (TinyVector<double,3> start,
                      TinyVector<double,3> end,
                      BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
                      Array<double,3> &data);
};

class NUBsplineClass_3d_d : public BsplineClass_3d_d
{
private:
  NUBspline_3d_d *Spline;
  NUgrid *xGrid, *yGrid, *zGrid;
public:
  double operator()(TinyVector<double,3> r);
  void evaluate (TinyVector<double,3> r,
                 double &val, TinyVector<double,3> &grad);
  void evaluate (TinyVector<double,3> r, double &val,
                 TinyVector<double,3> &grad, double &lapl);
  void evaluate (TinyVector<double,3> r, double &val,
                 TinyVector<double,3> &grad,
                 Tensor<double,3> &hess);
  NUBsplineClass_3d_d (TinyVector<double,3> start,
                       TinyVector<double,3> end,
                       TinyVector<double,3> ratio,
                       BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
                       Array<double,3> &data);
};


/////////////////////////////
// Complex Bspline classes //
/////////////////////////////
class BsplineClass_3d_z
{
public:
  virtual std::complex<double> operator()(TinyVector<double,3> r)                 = 0;
  virtual void evaluate (TinyVector<double,3> r,
                         std::complex<double> &val,
                         TinyVector<std::complex<double>,3> &grad)   = 0;
  virtual void evaluate (TinyVector<double,3> r, std::complex<double> &val,
                         TinyVector<std::complex<double>,3> &grad,
                         std::complex<double> &lapl)  = 0;
  virtual void evaluate (TinyVector<double,3> r, std::complex<double> &val,
                         TinyVector<std::complex<double>,3> &grad,
                         Tensor<std::complex<double>,3> &hess)                    = 0;
};

class UBsplineClass_3d_z : public BsplineClass_3d_z
{
private:
  UBspline_3d_z *Spline;
public:
  std::complex<double> operator()(TinyVector<double,3> r);
  void evaluate (TinyVector<double,3> r,
                 std::complex<double> &val, TinyVector<std::complex<double>,3> &grad);
  void evaluate (TinyVector<double,3> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,3> &grad,
                 std::complex<double> &lapl);
  void evaluate (TinyVector<double,3> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,3> &grad,
                 Tensor<std::complex<double>,3> &hess);
  UBsplineClass_3d_z (TinyVector<double,3> start,
                      TinyVector<double,3> end,
                      BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
                      Array<std::complex<double>,3> &data);
};

class NUBsplineClass_3d_z : public BsplineClass_3d_z
{
private:
  NUBspline_3d_z *Spline;
  NUgrid *xGrid, *yGrid, *zGrid;
public:
  std::complex<double> operator()(TinyVector<double,3> r);
  void evaluate (TinyVector<double,3> r,
                 std::complex<double> &val, TinyVector<std::complex<double>,3> &grad);
  void evaluate (TinyVector<double,3> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,3> &grad, std::complex<double> &lapl);
  void evaluate (TinyVector<double,3> r, std::complex<double> &val,
                 TinyVector<std::complex<double>,3> &grad,
                 Tensor<std::complex<double>,3> &hess);
  NUBsplineClass_3d_z (TinyVector<double,3> start,
                       TinyVector<double,3> end,
                       TinyVector<double,3> ratio,
                       BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
                       Array<std::complex<double>,3> &data);
};
#endif

}
#endif
