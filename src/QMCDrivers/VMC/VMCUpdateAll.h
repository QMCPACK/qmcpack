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

  void advanceWalker(Walker_t& thisWalker, bool recompute);
//       void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i);
//       void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//     , std::vector<MCWalkerConfiguration*>& wclone
//     , std::vector<QMCHamiltonian*>& hclone
//     , std::vector<RandomGenerator_t*>& rng
//     , std::vector<RealType>& ratio_i_0);

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

  void advanceWalker(Walker_t& thisWalker, bool recompute);

  RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);

private:
  /// Copy Constructor (disabled)
  VMCUpdateAllWithDrift(const VMCUpdateAllWithDrift& a): QMCUpdateBase(a) { }
  /// Copy operator (disabled).
  VMCUpdateAllWithDrift& operator=(const VMCUpdateAllWithDrift&)
  {
    return *this;
  }
};

}

#endif
