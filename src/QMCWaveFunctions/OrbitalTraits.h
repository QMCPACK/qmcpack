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
#ifndef QMCPLUSPLUS_ORBITALTRAITS_H
#define QMCPLUSPLUS_ORBITALTRAITS_H
#include <complex>
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus {

  inline double real(double a) {
    return a;
  }
  
  inline TinyVector<double,3> real(const TinyVector<double,3>& a) { 
    return a;
  }

  inline TinyVector<double,3> real(const TinyVector<complex<double>,3>& a) { 
    return TinyVector<double,3>(a[0].real(),a[1].real(),a[2].real());
  }

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
