//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <complex>
/** dummy trait class **/
template <class T> struct QMCDataTrait {};

/** specialization for double **/
template <>
struct QMCDataTrait<double> {
  typedef double real_type;
  typedef double value_type;
};

/** specialization for complex<double> **/
template <>
struct QMCDataTrait<std::complex<double> > {
  typedef double real_type;
  typedef std::complex<double> value_type;
};
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
