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
class DMCUpdateAllWithRejection : public QMCUpdateBase
{
public:
  /// Constructor.
  DMCUpdateAllWithRejection(MCWalkerConfiguration& w,
                            TrialWaveFunction& psi,
                            QMCHamiltonian& h,
                            RandomBase<FullPrecRealType>& rg);
  ///destructor
  ~DMCUpdateAllWithRejection() override;

  void advanceWalker(Walker_t& thisWalker, bool recompute) override;

  /// Copy Constructor (disabled)
  DMCUpdateAllWithRejection(const DMCUpdateAllWithRejection&) = delete;
  /// Copy operator (disabled).
  DMCUpdateAllWithRejection& operator=(const DMCUpdateAllWithRejection&) = delete;
};

class DMCUpdateAllWithKill : public QMCUpdateBase
{
public:
  /// Constructor.
  DMCUpdateAllWithKill(MCWalkerConfiguration& w,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       RandomBase<FullPrecRealType>& rg);
  ///destructor
  ~DMCUpdateAllWithKill() override;

  void advanceWalker(Walker_t& thisWalker, bool recompute) override;

  /// Copy Constructor (disabled)
  DMCUpdateAllWithKill(const DMCUpdateAllWithKill&) = delete;
  /// Copy operator (disabled).
  DMCUpdateAllWithKill& operator=(const DMCUpdateAllWithKill&) = delete;
};
} // namespace qmcplusplus

#endif
