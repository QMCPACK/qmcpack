//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RMC_UPDATEALL_H
#define QMCPLUSPLUS_RMC_UPDATEALL_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the RMC algorithm using all electron moves
 */
class RMCUpdateAllWithDrift : public QMCUpdateBase
{
public:
  /// Constructor.

  enum
  {
    SYM_ACTION,
    DMC_ACTION
  };

  RMCUpdateAllWithDrift(MCWalkerConfiguration& w,
                        TrialWaveFunction& psi,
                        QMCHamiltonian& h,
                        RandomBase<FullPrecRealType>& rg,
                        std::vector<int> act,
                        std::vector<int> tp);
  ~RMCUpdateAllWithDrift() override;
  void advanceWalker(Walker_t& thisWalker, bool recompute) override;
  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure) override;
  void advanceWalkersVMC();
  void advanceWalkersRMC();
  void checkReptile(WalkerIter_t it, WalkerIter_t it_end);
  void initWalkers(WalkerIter_t it, WalkerIter_t it_end) override;

  void accumulate(WalkerIter_t it, WalkerIter_t it_end);

  bool put(xmlNodePtr cur) override;

private:
  /// Copy Constructor (disabled)
  RMCUpdateAllWithDrift(const RMCUpdateAllWithDrift&) = delete;
  /// Copy operator (disabled).
  RMCUpdateAllWithDrift& operator=(const RMCUpdateAllWithDrift&) = delete;

  std::vector<int> Action, TransProb;

  bool scaleDrift;
  IndexType actionType;

  IndexType vmcSteps;
  IndexType equilSteps;
  IndexType vmcToDoSteps;
  IndexType equilToDoSteps;
};


} // namespace qmcplusplus

#endif
