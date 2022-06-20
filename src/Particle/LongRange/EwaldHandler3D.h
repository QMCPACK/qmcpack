//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_EWALD_HANLDER3D_H
#define QMCPLUSPLUS_EWALD_HANLDER3D_H

#include "LongRange/LRHandlerBase.h"
#include "coulomb_types.h"

namespace qmcplusplus
{
/* LR breakup for the standard Ewald method in 3D
 *
 */
class EwaldHandler3D : public LRHandlerBase
{
public:
  ///type of supercell
  int SuperCellEnum;
  /// Related to the Gaussian width: \f$ v_l = v(r)erf(\sigma r)\f$
  mRealType Sigma;
  ///Volume of the supercell
  mRealType Volume;
  ///Area of the supercell: always z is the slab direction
  mRealType Area;

  TinyVector<mRealType, 4> PreFactors;
  ///store |k|
  std::vector<mRealType> kMag;
  /// Constructor
  EwaldHandler3D(ParticleSet& ref, mRealType kc_in = -1.0) : LRHandlerBase(kc_in)
  {
    LRHandlerBase::ClassName = "EwaldHandler3D";
    Sigma = LR_kc = ref.getLattice().LR_kc;
  }

  /** "copy" constructor
   * @param aLR LRHandlerTemp
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  EwaldHandler3D(const EwaldHandler3D& aLR, ParticleSet& ref);

  LRHandlerBase* makeClone(ParticleSet& ref) const override { return new EwaldHandler3D(*this, ref); }

  void initBreakup(ParticleSet& ref) override;

  void Breakup(ParticleSet& ref, mRealType rs_in) override { initBreakup(ref); }

  void resetTargetParticleSet(ParticleSet& ref) override {}

  inline mRealType evaluate(mRealType r, mRealType rinv) const override { return erfc(r * Sigma) * rinv; }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline mRealType evaluateLR(mRealType r) const override { return erf(r * Sigma) / r; }

  inline mRealType evaluateSR_k0() const override
  {
    mRealType v0 = M_PI / Sigma / Sigma / Volume;
    return v0;
  }

  mRealType evaluate_vlr_k(mRealType k) const override;

  mRealType evaluateLR_r0() const override { return 2.0 * Sigma / std::sqrt(M_PI); }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline mRealType srDf(mRealType r, mRealType rinv) const override
  {
    return -2.0 * Sigma * std::exp(-Sigma * Sigma * r * r) / (std::sqrt(M_PI) * r) - erfc(Sigma * r) * rinv * rinv;
  }

  /**  evaluate the first derivative of the long range part (in real space) at r
   *
   * @param r  radius
   */
  inline mRealType lrDf(mRealType r) const override
  {
    mRealType rinv = 1.0 / r;
    return 2.0 * Sigma * std::exp(-Sigma * Sigma * r * r) / (std::sqrt(M_PI) * r) - erf(Sigma * r) * rinv * rinv;
  }

  void fillFk(const KContainer& KList);

  void fillYkgstrain(const KContainer& KList)
  {
    Fkgstrain.resize(KList.kpts_cart.size());
    const std::vector<int>& kshell(KList.kshell);
    MaxKshell = kshell.size() - 1;
    for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
    {
      mRealType uk = evalYkgstrain(std::sqrt(KList.ksq[ki]));
      while (ki < KList.kshell[ks + 1] && ki < Fkgstrain.size())
        Fkgstrain[ki++] = uk;
    }
  }

  void filldFk_dk(const KContainer& KList)
  {
    dFk_dstrain.resize(KList.kpts_cart.size());

    for (int ki = 0; ki < dFk_dstrain.size(); ki++)
    {
      dFk_dstrain[ki] = evaluateLR_dstrain(KList.kpts_cart[ki], std::sqrt(KList.ksq[ki]));
    }
  }

  //This returns the stress derivative of Fk, except for the explicit volume dependence.  The explicit volume dependence is factored away into V.
  inline SymTensor<mRealType, OHMMS_DIM> evaluateLR_dstrain(TinyVector<pRealType, OHMMS_DIM> k,
                                                            pRealType kmag) const override
  {
    SymTensor<mRealType, OHMMS_DIM> deriv_tensor = 0;

    for (int dim1 = 0; dim1 < OHMMS_DIM; dim1++)
      for (int dim2 = dim1; dim2 < OHMMS_DIM; dim2++)
      {
        deriv_tensor(dim1, dim2) =
            -evaldYkgstrain(kmag) * k[dim1] * k[dim2] / kmag; //- evaldFk_dk(kmag)*k[dim1]*k[dim2]/kmag ;

        if (dim1 == dim2)
          deriv_tensor(dim1, dim2) -= evalYkgstrain(kmag); //+ derivconst;
      }


    return deriv_tensor;
  }


  inline SymTensor<mRealType, OHMMS_DIM> evaluateSR_dstrain(TinyVector<pRealType, OHMMS_DIM> r,
                                                            pRealType rmag) const override
  {
    SymTensor<mRealType, OHMMS_DIM> deriv_tensor = 0;

    mRealType Sr_r = srDf(rmag, 1.0 / mRealType(rmag)) / mRealType(rmag);

    for (int dim1 = 0; dim1 < OHMMS_DIM; dim1++)
    {
      for (int dim2 = dim1; dim2 < OHMMS_DIM; dim2++)
      {
        deriv_tensor(dim1, dim2) = r[dim1] * r[dim2] * Sr_r;
      }
    }

    return deriv_tensor;
  }

  inline SymTensor<mRealType, OHMMS_DIM> evaluateSR_k0_dstrain() const override
  {
    mRealType v0 = -M_PI / Sigma / Sigma / Volume;
    SymTensor<mRealType, OHMMS_DIM> stress;
    for (int i = 0; i < OHMMS_DIM; i++)
      stress(i, i) = v0;

    return stress;
  }

  inline mRealType evaluateLR_r0_dstrain(int i, int j) const { return 0.0; }

  inline SymTensor<mRealType, OHMMS_DIM> evaluateLR_r0_dstrain() const override
  {
    SymTensor<mRealType, OHMMS_DIM> stress;
    return stress;
  }

private:
  inline mRealType evalYkgstrain(mRealType k) const
  {
    mRealType norm  = 4.0 * M_PI / Volume;
    mRealType denom = 4.0 * Sigma * Sigma;
    mRealType k2    = k * k;
    return norm * std::exp(-k2 / denom) / k2;
  }

  inline mRealType evaldYkgstrain(mRealType k) const
  {
    mRealType norm   = 4.0 * M_PI / Volume;
    mRealType denom  = 4.0 * Sigma * Sigma;
    mRealType sigma2 = Sigma * Sigma;
    mRealType k2     = k * k;
    return -norm * std::exp(-k2 / denom) * (denom + k2) / (k * k2 * 2.0 * sigma2);
  }
};
} // namespace qmcplusplus
#endif
