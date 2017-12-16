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

  void advanceWalker(Walker_t& thisWalker, bool recompute);

#if (__cplusplus >= 201103L)
  DMCUpdateAllWithRejection(const DMCUpdateAllWithRejection& a)=delete;
  DMCUpdateAllWithRejection& operator=(const DMCUpdateAllWithRejection&)=delete;
#else
private:
  /// Copy Constructor (disabled)
  DMCUpdateAllWithRejection(const DMCUpdateAllWithRejection& a): QMCUpdateBase(a) {}
  /// Copy operator (disabled).
  DMCUpdateAllWithRejection& operator=(const DMCUpdateAllWithRejection&)
  {
    return *this;
  }
#endif
};

class DMCUpdateAllWithKill: public QMCUpdateBase
{

public:

  /// Constructor.
  DMCUpdateAllWithKill(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  ~DMCUpdateAllWithKill();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

#if (__cplusplus >= 201103L)
  DMCUpdateAllWithKill(const DMCUpdateAllWithKill& a)=delete;
  DMCUpdateAllWithKill& operator=(const DMCUpdateAllWithKill&)=delete;
#else
private:
  /// Copy Constructor (disabled)
  DMCUpdateAllWithKill(const DMCUpdateAllWithKill& a): QMCUpdateBase(a) {}
  /// Copy operator (disabled).
  DMCUpdateAllWithKill& operator=(const DMCUpdateAllWithKill&)
  {
    return *this;
  }
#endif
};
}

#endif
