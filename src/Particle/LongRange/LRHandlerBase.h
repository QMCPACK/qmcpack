//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file LRHandlerBase.h
 * @brief Define LRHandlerBase and DummyLRHandler<typename Func>
 */
#ifndef QMCPLUSPLUS_LRHANLDERBASE_AND_DUMMY_H
#define QMCPLUSPLUS_LRHANLDERBASE_AND_DUMMY_H

#include "coulomb_types.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus
{
/** base class for LRHandlerTemp<FUNC,BASIS> and DummyLRHanlder<typename Func>
 */
struct LRHandlerBase
{
  DECLARE_COULOMB_TYPES

  /// Maxkimum Kshell for the given Kc
  int MaxKshell;
  /// Maximum k cutoff
  mRealType LR_kc;
  /// Maximum r cutoff
  mRealType LR_rc;
  ///Fourier component for all the k-point
  Vector<mRealType> Fk;
  ///Fourier component of the LR part, fit to optimize the gradients.
  Vector<mRealType> Fkg;
  ///Fourier component of the LR part of strain tensor, by optimized breakup.
  std::vector<SymTensor<mRealType, OHMMS_DIM>> dFk_dstrain;
  ///Vector of df_k/dk, fit as to optimize strains.
  Vector<mRealType> Fkgstrain;
  ///Fourier component for each k-shell
  Vector<mRealType> Fk_symm;

  ///Fourier component for each k-shell
  ///Coefficient
  std::vector<mRealType> coefs;
  ///Coefficient for gradient fit.
  std::vector<mRealType> gcoefs;
  ///Coefficient for strain fit.
  std::vector<mRealType> gstraincoefs;

  virtual mRealType evaluate_vlr_k(mRealType k) const = 0;


  //constructor
  explicit LRHandlerBase(mRealType kc) : MaxKshell(0), LR_kc(kc), LR_rc(0), ClassName("LRHandlerBase") {}

  // virtual destructor
  virtual ~LRHandlerBase() = default;

  //return r cutoff
  inline mRealType get_rc() const { return LR_rc; }
  //return k cutoff
  inline mRealType get_kc() const { return LR_kc; }

  inline mRealType evaluate_w_sk(const std::vector<int>& kshell, const pRealType* restrict sk) const
  {
    mRealType vk = 0.0;
    for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
    {
      mRealType u = 0;
      for (; ki < kshell[ks + 1]; ki++)
        u += (*sk++);
      vk += Fk_symm[ks] * u;
    }
    return vk;
  }

  /** evaluate \f$\sum_k F_{k} \rho^1_{-{\bf k}} \rho^2_{\bf k}\f$
   * @param kshell degeneracies of the vectors
   * @param rk1_r/i starting address of \f$\rho^1_{{\bf k}}\f$
   * @param rk2_r/i starting address of \f$\rho^2_{{\bf k}}\f$
   *
   * Valid for the strictly ordered k and \f$F_{k}\f$.
   */
  inline mRealType evaluate(const std::vector<int>& kshell,
                            const pRealType* restrict rk1_r,
                            const pRealType* restrict rk1_i,
                            const pRealType* restrict rk2_r,
                            const pRealType* restrict rk2_i) const
  {
    mRealType vk = 0.0;
    for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
    {
      mRealType u = 0;
      for (; ki < kshell[ks + 1]; ki++)
        u += ((*rk1_r++) * (*rk2_r++) + (*rk1_i++) * (*rk2_i++));
      vk += Fk_symm[ks] * u;
    }
    return vk;
  }

  /** Evaluate the long-range potential with the open BC for the D-1 direction */
  virtual mRealType evaluate_slab(pRealType z,
                                  const std::vector<int>& kshell,
                                  const pRealType* restrict rk1_r,
                                  const pRealType* restrict rk1_i,
                                  const pRealType* restrict rk2_r,
                                  const pRealType* restrict rk2_i) const
  {
    return 0.0;
  }

  /** evaluate \f$\sum_k F_{k} \rho^1_{-{\bf k}} \rho^2_{\bf k}\f$
   * and \f$\sum_k F_{k} \rho^1_{-{\bf k}} \rho^2_{\bf k}\f$
   * @param kshell degeneracies of the vectors
   * @param rk1 starting address of \f$\rho^1_{{\bf k}}\f$
   * @param rk2 starting address of \f$\rho^2_{{\bf k}}\f$
   *
   * Valid for the strictly ordered k and \f$F_{k}\f$.
   */
  inline void evaluateGrad(const ParticleSet& A,
                           const ParticleSet& B,
                           int specB,
                           std::vector<pRealType>& Zat,
                           std::vector<TinyVector<pRealType, OHMMS_DIM>>& grad1) const
  {
    const Matrix<pRealType>& e2ikrA_r = A.getSK().eikr_r;
    const Matrix<pRealType>& e2ikrA_i = A.getSK().eikr_i;
    const pRealType* rhokB_r          = B.getSK().rhok_r[specB];
    const pRealType* rhokB_i          = B.getSK().rhok_i[specB];
    const std::vector<PosType>& kpts  = A.getSimulationCell().getKLists().kpts_cart;
    for (int ki = 0; ki < Fk.size(); ki++)
    {
      PosType k = kpts[ki];
      for (int iat = 0; iat < Zat.size(); iat++)
      {
        grad1[iat] -= Zat[iat] * k * Fkg[ki] * (e2ikrA_r(iat, ki) * rhokB_i[ki] - e2ikrA_i(iat, ki) * rhokB_r[ki]);
      }
    }
  }

  ///FIX_PRECISION
  inline SymTensor<pRealType, OHMMS_DIM> evaluateStress(const std::vector<int>& kshell,
                                                        const pRealType* rhokA_r,
                                                        const pRealType* rhokA_i,
                                                        const pRealType* rhokB_r,
                                                        const pRealType* rhokB_i) const
  {
    SymTensor<pRealType, OHMMS_DIM> stress;
    for (int ki = 0; ki < dFk_dstrain.size(); ki++)
    {
      stress += (rhokA_r[ki] * rhokB_r[ki] + rhokA_i[ki] * rhokB_i[ki]) * dFk_dstrain[ki];
    }

    return stress;
  }

  /** evaluate \f$ v_{s}(k=0) = \frac{4\pi}{V}\int_0^{r_c} r^2 v_s(r) dr \f$
   */
  virtual mRealType evaluateSR_k0() const { return 0.0; }
  /** evaluate \f$ v_s(r=0) \f$ for the self-interaction term
   */
  virtual mRealType evaluateLR_r0() const { return 0.0; }

  ///These functions return the strain derivatives of all corresponding quantities
  /// in total energy.  See documentation (forthcoming).
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateLR_r0_dstrain() const { return 0; };
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateSR_k0_dstrain() const { return 0; };
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateLR_dstrain(TinyVector<pRealType, OHMMS_DIM> k, pRealType kmag) const
  {
    return 0;
  };
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateSR_dstrain(TinyVector<pRealType, OHMMS_DIM> r, pRealType rmag) const
  {
    return 0;
  };

  virtual void initBreakup(ParticleSet& ref)              = 0;
  virtual void Breakup(ParticleSet& ref, mRealType rs_in) = 0;
  virtual void resetTargetParticleSet(ParticleSet& ref)   = 0;

  virtual mRealType evaluate(mRealType r, mRealType rinv) const = 0;
  virtual mRealType evaluateLR(mRealType r) const               = 0;
  virtual mRealType srDf(mRealType r, mRealType rinv) const     = 0;

  virtual mRealType lrDf(mRealType r) const
  {
    APP_ABORT("Error: lrDf(r) is not implemented in " + ClassName + "\n");
    return 0.0;
  };

  /** make clone */
  virtual LRHandlerBase* makeClone(ParticleSet& ref) const = 0;

protected:
  std::string ClassName;
};

/** LRHandler without breakup.
 *
 * The template parameter Func should impelement operator()(kk) which
 * returns the Fourier component of a long-range function. Here kk
 * is \f$|{\bf k}|^2\f$.
 */
template<class Func>
struct DummyLRHandler : public LRHandlerBase
{
  Func myFunc;

  DummyLRHandler(mRealType kc = -1.0) : LRHandlerBase(kc) {}

  DummyLRHandler(const DummyLRHandler& aLR) : LRHandlerBase(aLR) {}

  void initBreakup(ParticleSet& ref) override
  {
    mRealType norm = 4.0 * M_PI / ref.getLattice().Volume;
    mRealType kcsq = LR_kc * LR_kc;
    auto& KList(ref.getSimulationCell().getKLists());
    int maxshell = KList.kshell.size() - 1;
    const auto& kk(KList.ksq);
    int ksh = 0, ik = 0;
    while (ksh < maxshell)
    {
      if (kk[ik] > kcsq)
        break; //exit
      ik = KList.kshell[++ksh];
    }
    MaxKshell = ksh;
    Fk_symm.resize(MaxKshell);
    Fk.resize(KList.kpts_cart.size());
    for (ksh = 0, ik = 0; ksh < MaxKshell; ksh++, ik++)
    {
      mRealType v  = norm * myFunc(kk[KList.kshell[ksh]]); //rpa=u0/kk[ik];
      Fk_symm[ksh] = v;
      for (; ik < KList.kshell[ksh + 1]; ik++)
        Fk[ik] = v;
    }
  }

  mRealType evaluate_vlr_k(mRealType k) const override { return 0.0; }
  mRealType evaluate(mRealType r, mRealType rinv) const override { return 0.0; }
  mRealType evaluateLR(mRealType r) const override { return 0.0; }
  mRealType srDf(mRealType r, mRealType rinv) const override { return 0.0; }
  void Breakup(ParticleSet& ref, mRealType rs_in) override {}
  void resetTargetParticleSet(ParticleSet& ref) override {}
  LRHandlerBase* makeClone(ParticleSet& ref) const override { return new DummyLRHandler<Func>(LR_kc); }
};

} // namespace qmcplusplus
#endif
