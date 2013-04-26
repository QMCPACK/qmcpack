//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_WFMC_UPDATEALL_H
#define QMCPLUSPLUS_WFMC_UPDATEALL_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{



class WFMCUpdateAllWithReweight: public QMCUpdateBase
{

public:
  int Elength,Eindex ;

  /// Constructor.
  WFMCUpdateAllWithReweight(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg, int length, int index);

  ///destructor
  ~WFMCUpdateAllWithReweight();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  /// Copy Constructor (disabled)
  WFMCUpdateAllWithReweight(const WFMCUpdateAllWithReweight& a): QMCUpdateBase(a) {}
  /// Copy operator (disabled).
  WFMCUpdateAllWithReweight& operator=(const WFMCUpdateAllWithReweight&)
  {
    return *this;
  }
};


}

#endif
/***************************************************************************
 * $RCSfile: VMCUpdateAll.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCUpdateAll.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
