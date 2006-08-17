//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
/**@file NonLocalTOperator.h
 *@brief Declaration of NonLocalTOperator
 */
#ifndef QMCPLUSPLUS_NONLOCALTRANSITIONOPERATOR_H
#define QMCPLUSPLUS_NONLOCALTRANSITIONOPERATOR_H

#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {

  struct NonLocalTOperator {

    typedef NonLocalData::RealType RealType;
    typedef NonLocalData::PosType PosType;

    RealType Tau;
    RealType Alpha;
    RealType Gamma;
    RealType plusFactor;
    RealType minusFactor;

    vector<NonLocalData> Txy;

    NonLocalTOperator();
    inline int size() const { return Txy.size();}

    inline int id(int ibar) const {
      return Txy[ibar].PID;
    }

    inline PosType delta(int ibar) const {
      return Txy[ibar].Delta;
    }

    bool put(xmlNodePtr cur);
    void reserve(int n);
    void reset();
    int selectMove(RealType prob);
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

