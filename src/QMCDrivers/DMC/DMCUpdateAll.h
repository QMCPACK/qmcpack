//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_DMC_ALLPARTICLE_UPDATE_H
#define QMCPLUSPLUS_DMC_ALLPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus
{

class DMCUpdateAllWithRejection: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdateAllWithRejection(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdateAllWithRejection();

  //void initWalkers(WalkerIter_t it, WalkerIter_t it_end);
  //void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);
  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
private:
  /// Copy Constructor (disabled)
  DMCUpdateAllWithRejection(const DMCUpdateAllWithRejection& a): QMCUpdateBase(a) {}
  /// Copy operator (disabled).
  DMCUpdateAllWithRejection& operator=(const DMCUpdateAllWithRejection&)
  {
    return *this;
  }
};

class DMCUpdateAllWithKill: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdateAllWithKill(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdateAllWithKill();

  //void initWalkers(WalkerIter_t it, WalkerIter_t it_end);
  //void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);
  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
private:
  /// Copy Constructor (disabled)
  DMCUpdateAllWithKill(const DMCUpdateAllWithKill& a): QMCUpdateBase(a) {}
  /// Copy operator (disabled).
  DMCUpdateAllWithKill& operator=(const DMCUpdateAllWithKill&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCUpdateAll.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
