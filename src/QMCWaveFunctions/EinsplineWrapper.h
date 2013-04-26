//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: esler@uiuc.edu                                     //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////
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
  virtual complex<double> operator()(TinyVector<double,2> r)                 = 0;
  virtual void evaluate (TinyVector<double,2> r,
                         complex<double> &val,
                         TinyVector<complex<double>,2> &grad)   = 0;
  virtual void evaluate (TinyVector<double,2> r, complex<double> &val,
                         TinyVector<complex<double>,2> &grad,
                         complex<double> &lapl)  = 0;
  virtual void evaluate (TinyVector<double,2> r, complex<double> &val,
                         TinyVector<complex<double>,2> &grad,
                         Tensor<complex<double>,2> &hess)                    = 0;
};

class UBsplineClass_2d_z : public BsplineClass_2d_z
{
private:
  UBspline_2d_z *Spline;
public:
  complex<double> operator()(TinyVector<double,2> r);
  void evaluate (TinyVector<double,2> r,
                 complex<double> &val, TinyVector<complex<double>,2> &grad);
  void evaluate (TinyVector<double,2> r, complex<double> &val,
                 TinyVector<complex<double>,2> &grad,
                 complex<double> &lapl);
  void evaluate (TinyVector<double,2> r, complex<double> &val,
                 TinyVector<complex<double>,2> &grad,
                 Tensor<complex<double>,2> &hess);
  UBsplineClass_2d_z (TinyVector<double,2> start,
                      TinyVector<double,2> end,
                      BCtype_z xBC, BCtype_z yBC,
                      Array<complex<double>,2> &data);
};

class NUBsplineClass_2d_z : public BsplineClass_2d_z
{
private:
  NUBspline_2d_z *Spline;
  NUgrid *xGrid, *yGrid, *zGrid;
public:
  complex<double> operator()(TinyVector<double,2> r);
  void evaluate (TinyVector<double,2> r,
                 complex<double> &val, TinyVector<complex<double>,2> &grad);
  void evaluate (TinyVector<double,2> r, complex<double> &val,
                 TinyVector<complex<double>,2> &grad, complex<double> &lapl);
  void evaluate (TinyVector<double,2> r, complex<double> &val,
                 TinyVector<complex<double>,2> &grad,
                 Tensor<complex<double>,2> &hess);
  NUBsplineClass_2d_z (TinyVector<double,2> start,
                       TinyVector<double,2> end,
                       TinyVector<double,2> ratio,
                       BCtype_z xBC, BCtype_z yBC,
                       Array<complex<double>,2> &data);
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
  virtual complex<double> operator()(TinyVector<double,3> r)                 = 0;
  virtual void evaluate (TinyVector<double,3> r,
                         complex<double> &val,
                         TinyVector<complex<double>,3> &grad)   = 0;
  virtual void evaluate (TinyVector<double,3> r, complex<double> &val,
                         TinyVector<complex<double>,3> &grad,
                         complex<double> &lapl)  = 0;
  virtual void evaluate (TinyVector<double,3> r, complex<double> &val,
                         TinyVector<complex<double>,3> &grad,
                         Tensor<complex<double>,3> &hess)                    = 0;
};

class UBsplineClass_3d_z : public BsplineClass_3d_z
{
private:
  UBspline_3d_z *Spline;
public:
  complex<double> operator()(TinyVector<double,3> r);
  void evaluate (TinyVector<double,3> r,
                 complex<double> &val, TinyVector<complex<double>,3> &grad);
  void evaluate (TinyVector<double,3> r, complex<double> &val,
                 TinyVector<complex<double>,3> &grad,
                 complex<double> &lapl);
  void evaluate (TinyVector<double,3> r, complex<double> &val,
                 TinyVector<complex<double>,3> &grad,
                 Tensor<complex<double>,3> &hess);
  UBsplineClass_3d_z (TinyVector<double,3> start,
                      TinyVector<double,3> end,
                      BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
                      Array<complex<double>,3> &data);
};

class NUBsplineClass_3d_z : public BsplineClass_3d_z
{
private:
  NUBspline_3d_z *Spline;
  NUgrid *xGrid, *yGrid, *zGrid;
public:
  complex<double> operator()(TinyVector<double,3> r);
  void evaluate (TinyVector<double,3> r,
                 complex<double> &val, TinyVector<complex<double>,3> &grad);
  void evaluate (TinyVector<double,3> r, complex<double> &val,
                 TinyVector<complex<double>,3> &grad, complex<double> &lapl);
  void evaluate (TinyVector<double,3> r, complex<double> &val,
                 TinyVector<complex<double>,3> &grad,
                 Tensor<complex<double>,3> &hess);
  NUBsplineClass_3d_z (TinyVector<double,3> start,
                       TinyVector<double,3> end,
                       TinyVector<double,3> ratio,
                       BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
                       Array<complex<double>,3> &data);
};
#endif

}
#endif
