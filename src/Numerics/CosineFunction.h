//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_SINECOSINEFUNCTION_H
#define OHMMS_SINECOSINEFUNCTION_H
#include <math.h>

template<class T, class PT>
struct Sine3D
{
  typedef T value_type;
  typedef PT pos_type;

  Sine3D(value_type kx = 0, value_type ky =0, value_type kz = 0)
  {
    set(kx,ky,kz);
  }

  void set(value_type kx, value_type ky, value_type kz)
  {
    const double twopi = 2.0*M_PI;
    Kx = twopi*(kx+0.5);
    Ky = twopi*(ky+0.5);
    Kz = twopi*(kz+0.5);
    Knorm2 = Kx*Kx+Ky*Ky+Kz*Kz;
  }

  inline value_type evaluate(const PT& r)
  {
    return sin(Kx*r[0])*sin(Ky*r[1])*sin(Kz*r[2]);
  }

  inline value_type evaluate(const PT& r, PT& gr, value_type& lap)
  {
    value_type v;
    gr = gradient(r,v);
    lap = -Knorm2*v;
    return v;
  }


  inline PT gradient(const PT& r, value_type& v)
  {
    v = evaluate(r);
    return PT(Kx*cos(Kx*r[0])*sin(Ky*r[1])*sin(Kz*r[2]),
              Ky*sin(Kx*r[0])*cos(Ky*r[1])*sin(Kz*r[2]),
              Kz*sin(Kx*r[0])*sin(Ky*r[1])*cos(Kz*r[2]));
  }


  // update internal variables and return the funcional value at r
  inline value_type laplacian(const PT& r)
  {
    Lap = laplacian(r,Val,Grad);
    return Val;
  }

  inline value_type operator()(const PT& r)
  {
    return evaluate(r);
  }

  value_type Kx, Ky, Kz, Knorm2;
  value_type Val, Lap;
  pos_type   Grad;
};

template<class T, class PT>
struct Cosine3D
{
  typedef T value_type;
  typedef PT pos_type;

  Cosine3D(value_type kx = 0, value_type ky =0, value_type kz = 0)
  {
    set(kx,ky,kz);
  }

  void set(value_type kx, value_type ky, value_type kz)
  {
    const double twopi = 2.0*M_PI;
    Kx = twopi*(kx+0.5);
    Ky = twopi*(ky+0.5);
    Kz = twopi*(kz+0.5);
    Knorm2 = Kx*Kx+Ky*Ky+Kz*Kz;
  }

  inline value_type evaluate(const PT& r)
  {
    return cos(Kx*r[0])*cos(Ky*r[1])*cos(Kz*r[2]);
  }

  inline PT gradient(const PT& r, value_type& v)
  {
    v = evaluate(r);
    return PT(-Kx*sin(Kx*r[0])*cos(Ky*r[1])*cos(Kz*r[2]),
              -Ky*cos(Kx*r[0])*sin(Ky*r[1])*cos(Kz*r[2]),
              -Kz*cos(Kx*r[0])*cos(Ky*r[1])*sin(Kz*r[2]));
  }

  inline value_type laplacian(const PT& r, value_type& v, PT& gr)
  {
    gr = gradient(r,v);
    return -Knorm2*v;
  }

  // update internal variables
  inline value_type laplacian(const PT& r)
  {
    Lap = laplacian(r,Val,Grad);
    return Val;
  }

  inline value_type operator()(const PT& r)
  {
    return evaluate(r);
  }

  value_type Kx, Ky, Kz, Knorm2;
  value_type Val, Lap;
  pos_type   Grad;
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
