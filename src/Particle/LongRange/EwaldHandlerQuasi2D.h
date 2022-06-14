//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EWALD_QUASI2D_HANDLER_H
#define QMCPLUSPLUS_EWALD_QUASI2D_HANDLER_H

#include "Particle/LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
/* LR breakup for the standard Ewald method in quasi-2D (slab geometry)
 */
class EwaldHandlerQuasi2D : public LRHandlerBase
{
public:
  EwaldHandlerQuasi2D(ParticleSet& ref, mRealType kc_in=-1.0);

  // copy constructor
  LRHandlerBase* makeClone(ParticleSet& ref) const override { return new EwaldHandlerQuasi2D(*this); }

  // short-range part
  inline mRealType evaluate(mRealType r, mRealType rinv) const override { return erfc(alpha*r) * rinv; }
  inline mRealType evaluateLR_r0() const override { return 2.0 * alpha / std::sqrt(M_PI); }

  // long-range part
  inline mRealType evaluateLR(mRealType r) const override { return erf(alpha*r) / r; }
  inline mRealType evaluateSR_k0() const override { return slab_vsr_k0(0)/area; }
  void fillFk(const KContainer& KList);

  // z-dependent long-range part
  mRealType evaluate_slab(pRealType z,
                          const std::vector<int>& kshell,
                          const pRealType* restrict rk1_r,
                          const pRealType* restrict rk1_i,
                          const pRealType* restrict rk2_r,
                          const pRealType* restrict rk2_i) const override;

  // begin required overrides
  inline mRealType srDf(mRealType r, mRealType rinv) const override
  {
    throw std::runtime_error("Quasi2D Ewald sr derivative not implemented");
  }
  inline mRealType lrDf(mRealType r) const override
  {
    throw std::runtime_error("Quasi2D Ewald lr derivative not implemented");
  }
  virtual mRealType evaluate_vlr_k(mRealType k) const override
  {
    throw std::runtime_error("Quasi2D Ewald vlr_k not implemented");
  }
  void initBreakup(ParticleSet& ref) override {}
  void Breakup(ParticleSet& ref, mRealType rs_in) override { initBreakup(ref); }
  void resetTargetParticleSet(ParticleSet& ref) override {}
  // overrides end
private:
  mRealType alpha;
  mRealType area;
  ///store |k|
  std::vector<mRealType> kmags;
  // f(z, k) = e^{z k} erfc(k/(2\alpha)+\alpha*z)
  mRealType slabFunc(mRealType z, mRealType k) const;
  // log( f(z, k) )
  mRealType slabLogf(mRealType z, mRealType k) const;
  mRealType slab_vsr_k0(mRealType z) const;
};
} // qmcplusplus
#endif
