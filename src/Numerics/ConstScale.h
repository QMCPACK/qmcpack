//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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

#ifndef OHMMS_CONSTSCALINGSKMODEL_H
#define OHMMS_CONSTSCALINGSKMODEL_H

// SK<T,L1,L2,MID> where MID = 0 for no scaling
// dum scaling function which returns a constant function
struct ConstScale {
  double C;
  ConstScale(double c=1.0):C(c){ }
  ~ConstScale() { }
  inline double operator()(double r) { return C;}
  inline double operator()(double r, double& vr) { vr = C/r; return 0.0e0;}
  inline double operator()(double r, double& vr, double& dvr) {
    vr = C/r/r; dvr = 0.0e0; return 0.0e0;
  }
};

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * OHMMS_VERSION_ID: $Id$ 
 ***************************************************************************/
