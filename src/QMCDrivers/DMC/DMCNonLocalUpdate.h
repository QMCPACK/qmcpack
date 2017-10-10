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
    
    


#ifndef QMCPLUSPLUS_DMC_NONLOCAL_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_DMC_NONLOCAL_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus
{

class DMCNonLocalUpdate: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCNonLocalUpdate(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                    QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCNonLocalUpdate();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  /// Copy Constructor (disabled)
  DMCNonLocalUpdate(const DMCNonLocalUpdate& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  DMCNonLocalUpdate& operator=(const DMCNonLocalUpdate&)
  {
    return *this;
  }

};

class DMCNonLocalUpdatePbyP: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCNonLocalUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                        QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCNonLocalUpdatePbyP();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  /// Copy Constructor (disabled)
  DMCNonLocalUpdatePbyP(const DMCNonLocalUpdatePbyP& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  DMCNonLocalUpdatePbyP& operator=(const DMCNonLocalUpdatePbyP&)
  {
    return *this;
  }
  std::vector<NewTimer*> myTimers;

};

class DMCNonLocalUpdatePbyPFast: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCNonLocalUpdatePbyPFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                            QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCNonLocalUpdatePbyPFast();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  /// Copy Constructor (disabled)
  DMCNonLocalUpdatePbyPFast(const DMCNonLocalUpdatePbyP& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  DMCNonLocalUpdatePbyPFast& operator=(const DMCNonLocalUpdatePbyP&)
  {
    return *this;
  }
  std::vector<NewTimer*> myTimers;

};

}

#endif
