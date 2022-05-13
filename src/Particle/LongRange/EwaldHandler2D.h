//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EWALD2D_HANDLER_H
#define QMCPLUSPLUS_EWALD2D_HANDLER_H

#include "Particle/LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
/* LR breakup for the standard Ewald method in 2D
 * 1/r = \left[ \sum_r erfc(\alpha*r)/r + erf(\alpha*r)/r \right]
 *     = \left[ \sum_r erfc(\alpha*r)/r \right] +
 *       \left[ \sum_k 2\pi/k * erfc(k/(2*\alpha)) \right]
 */
class EwaldHandler2D : public LRHandlerBase
{
public:
  EwaldHandler2D(ParticleSet& ref, mRealType kc_in=-1.0);

  // copy constructor
  LRHandlerBase* makeClone(ParticleSet& ref) const override { return new EwaldHandler2D(*this); }

  // short-range part
  inline mRealType evaluate(mRealType r, mRealType rinv) const override { return erfc(alpha*r) * rinv; }
  inline mRealType evaluateLR_r0() const override { return 2.0 * alpha / std::sqrt(M_PI); }

  // long-range part
  inline mRealType evaluateLR(mRealType r) const override { return erf(alpha*r) / r; }
  inline mRealType evaluateSR_k0() const override { return 2.0 * std::sqrt(M_PI) / (alpha*area); }
  void fillFk(const KContainer& KList);

  // begin required overrides
  inline mRealType srDf(mRealType r, mRealType rinv) const override
  {
    throw std::runtime_error("2D Ewald sr derivative not implemented");
  }
  inline mRealType lrDf(mRealType r) const override
  {
    throw std::runtime_error("2D Ewald lr derivative not implemented");
  }
  virtual mRealType evaluate_vlr_k(mRealType k) const override
  {
    throw std::runtime_error("2D Ewald vlr_k not implemented");
  }
  void initBreakup(ParticleSet& ref) override {}
  void Breakup(ParticleSet& ref, mRealType rs_in) override { initBreakup(ref); }
  void resetTargetParticleSet(ParticleSet& ref) override {}
  // overrides end
private:
  mRealType alpha;
  mRealType area;
};
} // qmcplusplus
#endif
