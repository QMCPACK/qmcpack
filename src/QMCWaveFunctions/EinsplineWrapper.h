//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: esler@uiuc.edu                                     //
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)             //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////
#ifndef EINSPLINE_WRAPPER_H
#define EINSPLINE_WRAPPER_H

#include <einspline/bspline.h>
#include <einspline/nubspline.h>
#include "Numerics/HDFNumericAttrib.h"

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
		       TinyVector<double,3> end, double ratio,
		       BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
		       Array<double,3> &data);
};

#endif
