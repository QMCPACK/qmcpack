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
#ifndef QMCPLUSPLUS_VMC_UPDATEALL_H
#define QMCPLUSPLUS_VMC_UPDATEALL_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCUpdateAll: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi,
               QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdateAll();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//       void advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone, vector<RandomGenerator_t*>& rng, vector<RealType>& c_i);
//       void estimateNormWalkers(vector<TrialWaveFunction*>& pclone
//     , vector<MCWalkerConfiguration*>& wclone
//     , vector<QMCHamiltonian*>& hclone
//     , vector<RandomGenerator_t*>& rng
//     , vector<RealType>& ratio_i_0);

private:
  /// Copy Constructor (disabled)
  VMCUpdateAll(const VMCUpdateAll& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  VMCUpdateAll& operator=(const VMCUpdateAll&)
  {
    return *this;
  }
};

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
 */
class VMCUpdateAllWithDrift: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                        QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdateAllWithDrift();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//       void advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone, vector<RandomGenerator_t*>& rng, vector<RealType>& c_i);

  RealType advanceWalkerForEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios);

private:
  /// Copy Constructor (disabled)
  VMCUpdateAllWithDrift(const VMCUpdateAllWithDrift& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  VMCUpdateAllWithDrift& operator=(const VMCUpdateAllWithDrift&)
  {
    return *this;
  }
};


/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
//   class VMCUpdateAllSampleRN: public QMCUpdateBase
//     {
//     public:
//       /// Constructor.
//       VMCUpdateAllSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide,
//                             QMCHamiltonian& h, RandomGenerator_t& rg);
//
//       ~VMCUpdateAllSampleRN();
//
//       void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//       void setLogEpsilon(RealType eps) { logEpsilon=eps; }
//
//     private:
//       /// Copy Constructor (disabled)
//       VMCUpdateAllSampleRN(const VMCUpdateAllSampleRN& a): QMCUpdateBase(a) { }
//       /// Copy operator (disabled).
//       VMCUpdateAllSampleRN& operator=(const VMCUpdateAllSampleRN&)
//       {
//         return *this;
//       }
//       RealType logEpsilon;
//     };

}

#endif
/***************************************************************************
 * $RCSfile: VMCUpdateAll.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCUpdateAll.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
