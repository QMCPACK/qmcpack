//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RMC_PARTICLEBYPARTICLE_UPDATE_H
#define QMCPLUSPLUS_RMC_PARTICLEBYPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the RMC algorithm using all electron moves
 */
class RMCUpdatePbyPWithDrift : public QMCUpdateBase
{
public:
  /// Constructor.
  RMCUpdatePbyPWithDrift(MCWalkerConfiguration& w,
                         TrialWaveFunction& psi,
                         QMCHamiltonian& h,
                         RandomBase<FullPrecRealType>& rg,
                         std::vector<int> act,
                         std::vector<int> tp);
  ~RMCUpdatePbyPWithDrift() override;

  enum
  {
    SYM_ACTION,
    DMC_ACTION
  };

  void advanceWalkersVMC();
  void advanceWalkersRMC();
  void advanceWalker(Walker_t& thisWalker, bool recompute) override;
  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure) override;
  void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end) override;
  void initWalkers(WalkerIter_t it, WalkerIter_t it_end) override;
  void accumulate(WalkerIter_t it, WalkerIter_t it_end);

  bool put(xmlNodePtr cur) override;

private:
  /// Copy Constructor (disabled)
  RMCUpdatePbyPWithDrift(const RMCUpdatePbyPWithDrift&) = delete;
  /// Copy operator (disabled).
  RMCUpdatePbyPWithDrift& operator=(const RMCUpdatePbyPWithDrift&) = delete;

  std::vector<int> Action, TransProb;
  bool scaleDrift;
  IndexType actionType;

  NewTimer& advance_timer_;
  NewTimer& movepbyp_timer_;
  NewTimer& update_mbo_timer_;
  NewTimer& energy_timer_;

  IndexType vmcSteps;
  IndexType equilSteps;
  IndexType vmcToDoSteps;
  IndexType equilToDoSteps;
};


} // namespace qmcplusplus

#endif
