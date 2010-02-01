//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Kenneth P. Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/* @file stdfunc.h
 * @brief helper functions to support missing std functions
 */
#ifndef QMCPLUSPLUS_STDLIB_PORT_H
#define QMCPLUSPLUS_STDLIB_PORT_H
#include <config.h>
#include <cmath>
#include <cstdlib>
#ifndef TWOPI
#ifndef M_PI
#define TWOPI 6.2831853071795862
#else
#define TWOPI (2*M_PI)
#endif /* M_PI */
#endif /* TWOPI */

#if !defined(HAVE_STD_ROUND)
template<typename T> inline T round(T x)
{
  T dmy;
  x=modf(x,&dmy); 
  return x-static_cast<int>(x*2.0);
}
#endif

#if !defined(HAVE_SINCOS)
#include <cmath>
template<typename T> 
inline void sincos(T a, T* restrict s, T*  restrict c)
{
  *s=std::sin(a);
  *c=std::cos(a);
}
#endif

#if defined(HAVE_MKL_VML)
#include <mkl_vml_functions.h>

inline void vec_sqrt(int n, const double* restrict in, double* restrict out)
{
  vdSqrt(n,in,out);
}

inline void vec_inv(int n, const double* restrict in, double* restrict out)
{
  vdInv(n,in,out);
}

inline void vec_inv_sqrt(int n, const double* restrict in, double* restrict out)
{
  vdInvSqrt(n,in,out);
}

inline void vec_sqrt(int n, const float* restrict in, float* restrict out)
{
  vsSqrt(n,in,out);
}

#else
#if defined(HAVE_MASSV)
#include <mass.h>
#include <massv.h>
inline void vec_sqrt(int n, double* restrict in, double* restrict out)
{
    vsqrt(out,in,&n);
}
inline void vec_sqrt(int n, float* restrict in, float* restrict out)
{
    vssqrt(out,in,&n);
}

inline void vec_inv_sqrt(int n, double* restrict in, double* restrict out)
{
  vrsqrt(out,in,&n);
}

inline void vec_inv_sqrt(int n, float* restrict in, float* restrict out)
{
  vsrsqrt(out,in,&n);
}
#else
template<typename T>
inline void vec_sqrt(int n, const T* restrict in, T* restrict out)
{
  for(int i=0; i<n; ++i) out[i]=sqrt(in[i]);
}
#endif
template<typename T>
inline void vec_inv(int n, const T* restrict in, T* restrict out)
{
  for(int i=0; i<n; ++i) out[i]=1.0/in[i];
}

#endif

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3310 $   $Date: 2008-10-29 19:21:31 -0500 (Wed, 29 Oct 2008) $
 * $Id: stdfunc.h 3310 2008-10-30 00:21:31Z jnkim $ 
 ***************************************************************************/
