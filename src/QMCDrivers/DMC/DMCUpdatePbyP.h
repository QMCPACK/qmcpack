//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus
{

class DMCUpdatePbyPWithRejection: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdatePbyPWithRejection(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                             QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdatePbyPWithRejection();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  std::vector<NewTimer*> myTimers;
};

class DMCUpdatePbyPWithRejectionFast: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                                 QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdatePbyPWithRejectionFast();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  std::vector<NewTimer*> myTimers;
};


class DMCUpdatePbyPWithKill: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdatePbyPWithKill(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                        QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdatePbyPWithKill();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:

  /// Copy Constructor (disabled)
  DMCUpdatePbyPWithKill(const DMCUpdatePbyPWithKill& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  DMCUpdatePbyPWithKill& operator=(const DMCUpdatePbyPWithKill&)
  {
    return *this;
  }
  std::vector<NewTimer*> myTimers;

};

class RNDMCUpdatePbyPFast: public QMCUpdateBase
{

public:

  RNDMCUpdatePbyPFast(MCWalkerConfiguration& w, MCWalkerConfiguration& wg, TrialWaveFunction& psi, TrialWaveFunction& guide,
                      QMCHamiltonian& h, RandomGenerator_t& rg);

  ~RNDMCUpdatePbyPFast();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

//     void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//     , std::vector<MCWalkerConfiguration*>& wclone
//     , std::vector<QMCHamiltonian*>& hclone
//     , std::vector<RandomGenerator_t*>& rng
//     , std::vector<RealType>& ratio_i_0);

private:
  MCWalkerConfiguration W_G;
  std::vector<NewTimer*> myTimers;
  int maxS;
  RealType efn;
  int estimateCrossings, maxcopy;

};


class RNDMCUpdatePbyPCeperley: public QMCUpdateBase
{

public:

  /// Constructor.
  RNDMCUpdatePbyPCeperley(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                          QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~RNDMCUpdatePbyPCeperley();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
  void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

private:
  std::vector<NewTimer*> myTimers;
  int maxS;
  RealType efn;
  int estimateCrossings;
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCUpdatePbyP.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
