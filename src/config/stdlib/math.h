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

#if !defined(HAVE_STD_ROUND)
template<typename T> inline T round(T x)
{
  T dmy;
  x=modf(x,&dmy); 
  return x-static_cast<int>(x*2.0);
}
#endif

#if !defined(HAVE_SINCOS)
template<typename T> 
inline void sincos(T a, T* restrict s, T*  restrict c)
{
  *s=std::sin(a);
  *c=std::cos(a);
}
#endif

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3310 $   $Date: 2008-10-29 19:21:31 -0500 (Wed, 29 Oct 2008) $
 * $Id: stdfunc.h 3310 2008-10-30 00:21:31Z jnkim $ 
 ***************************************************************************/
