//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
