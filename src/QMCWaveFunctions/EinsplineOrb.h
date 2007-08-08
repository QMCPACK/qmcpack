//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)             //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EINSPLINE_ORB_H
#define QMCPLUSPLUS_EINSPLINE_ORB_H

#include <einspline/bspline.h>
#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus {

  template<typename T, int N>
  class EinsplineOrb : public QMCTraits
  {
  public:
    PosType Center;
    RealType Radius, Energy;
    bool Localized;
    
    inline void evaluate (PosType r, T &psi);
    inline void evaluate (PosType r, T &psi, TinyVector<T,N> &grad, T &lapl);  
    inline void evaluate (PosType r, T &psi, TinyVector<T,N> &grad,
			  Tensor<T,N> &hess);
    void Read (hid_t h5file, string baseName);
  };
  
  
  
  template<>
  class EinsplineOrb<double,2> : public QMCTraits
  {
  public:
    PosType Center;
    RealType Radius, Energy;
    bool Localized;
    UBspline_2d_d *Spline;
    
    inline void evaluate (PosType r, double &psi) 
    {
      eval_UBspline_2d_d (Spline, r[0], r[1], &psi);
    }
    inline void evaluate (PosType r, double &psi, TinyVector<double,2> &grad,
			  double &lapl)
    {
      eval_UBspline_2d_d_vgl (Spline, r[0], r[1], &psi, &(grad[0]), &lapl);
    }
    inline void evaluate (PosType r, double &psi, TinyVector<double,2> &grad,
			  Tensor<double,2> &hess)
    {
      eval_UBspline_2d_d_vgh (Spline, r[0], r[1], &psi, &(grad[0]), &(hess(0,0)));
    }
    void read (hid_t h5file, string baseName);
    EinsplineOrb() : Center(PosType()), Radius(0.0), Energy(0.0), Localized(false)
    {
    }
  };
  
  template<>
  class EinsplineOrb<complex<double>,2> : public QMCTraits
  {
  public:
    PosType Center;
    RealType Radius, Energy;
    bool Localized;
    UBspline_2d_z *Spline;
    
    inline void evaluate (PosType r, complex<double> &psi) 
    {
      eval_UBspline_2d_z (Spline, r[0], r[1], &psi);
    }
    inline void evaluate (PosType r, complex<double> &psi, 
			  TinyVector<complex<double>,3> &grad,
			  complex<double> &lapl)
    {
      eval_UBspline_2d_z_vgl (Spline, r[0], r[1],
			      &psi, &(grad[0]), &lapl);
    }
    inline void evaluate (PosType r, complex<double> &psi, 
			  TinyVector<complex<double>,2> &grad,
			  Tensor<complex<double>,2> &hess)
    {
      eval_UBspline_2d_z_vgh (Spline, r[0], r[1],
			      &psi, &(grad[0]), &(hess(0,0)));
    }
    void read (hid_t h5file, string baseName);

    EinsplineOrb() : Center(PosType()), Radius(0.0), Energy(0.0), Localized(false)
    {
    }

  };
  
  
  
  template<>
  class EinsplineOrb<double,3> : public QMCTraits
  {
  public:
    PosType Center;
    RealType Radius, Energy;
    bool Localized;
    UBspline_3d_d *Spline;
    
    inline void evaluate (PosType r, double &psi) 
    {
      eval_UBspline_3d_d (Spline, r[0], r[1], r[2], &psi);
    }
    inline void evaluate (PosType r, double &psi, TinyVector<double,3> &grad,
			  double &lapl)
    {
      eval_UBspline_3d_d_vgl (Spline, r[0], r[1], r[2],
			      &psi, &(grad[0]), &lapl);
    }
    inline void evaluate (PosType r, double &psi, TinyVector<double,3> &grad,
			  Tensor<double,3> &hess)
    {
      eval_UBspline_3d_d_vgh (Spline, r[0], r[1], r[2], 
			      &psi, &(grad[0]), &(hess(0,0)));
    }
    void read (hid_t h5file, string baseName);

    EinsplineOrb() : Center(PosType()), Radius(0.0), Energy(0.0), Localized(false)
    {
    }
  };
  
  template<>
  class EinsplineOrb<complex<double>,3> : public QMCTraits
  {
  public:
    PosType Center;
    RealType Radius, Energy;
    bool Localized;
    UBspline_3d_z *Spline;
    
    inline void evaluate (PosType r, complex<double> &psi) 
    {
      eval_UBspline_3d_z (Spline, r[0], r[1], r[2], &psi);
    }
    inline void evaluate (PosType r, complex<double> &psi, 
			  TinyVector<complex<double>,3> &grad,
			  complex<double> &lapl)
    {
      eval_UBspline_3d_z_vgl (Spline, r[0], r[1], r[2],
			      &psi, &(grad[0]), &lapl);
    }
    inline void evaluate (PosType r, complex<double> &psi, 
			  TinyVector<complex<double>,3> &grad,
			  Tensor<complex<double>,3> &hess)
    {
      eval_UBspline_3d_z_vgh (Spline, r[0], r[1], r[2], 
			      &psi, &(grad[0]), &(hess(0,0)));
    }
    void read (hid_t h5file, string baseName);

    EinsplineOrb() : Center(PosType()), Radius(0.0), Energy(0.0), Localized(false)
    {
    }
  };
}

#endif
