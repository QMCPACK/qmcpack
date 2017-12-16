//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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
  mRealType  LR_kc;
  /// Maximum r cutoff
  mRealType  LR_rc;
  ///Fourier component for all the k-point
  Vector<mRealType> Fk;
  ///Fourier component of the LR part, fit to optimize the gradients.
  Vector<mRealType> Fkg;
  ///Fourier component of the LR part of strain tensor, by optimized breakup.  
  std::vector< SymTensor<mRealType, OHMMS_DIM> > dFk_dstrain;
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

  
  //constructor
  explicit LRHandlerBase(mRealType kc):LR_kc(kc) {}

  // virtual destructor
  virtual ~LRHandlerBase() {}

  //return r cutoff
  inline mRealType get_rc() const
  {
    return LR_rc;
  }
  //return k cutoff
  inline mRealType get_kc() const
  {
    return LR_kc;
  }

  /** evaluate \f$\sum_k F_{k} \rho^1_{-{\bf k} \rho^2_{\bf k}\f$
   * @param kshell degeneracies of the vectors
   * @param rk1 starting address of \f$\rho^1_{{\bf k}\f$
   * @param rk2 starting address of \f$\rho^2_{{\bf k}\f$
   *
   * Valid for the strictly ordered k and \f$F_{k}\f$.
   */
  inline mRealType evaluate(const std::vector<int>& kshell
                           , const pComplexType* restrict rk1, const pComplexType* restrict rk2)
  {
    mRealType vk=0.0;
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType u=0;
      for(; ki<kshell[ks+1]; ki++,rk1++,rk2++)
        u += ((*rk1).real()*(*rk2).real()+(*rk1).imag()*(*rk2).imag());
      vk += Fk_symm[ks]*u;
    }
    return vk;
  }

  inline mRealType evaluate(const std::vector<int>& kshell
                           , const pRealType* restrict rk1_r, const pRealType* restrict rk1_i
                           , const pRealType* restrict rk2_r, const pRealType* restrict rk2_i)
  {
    mRealType vk=0.0;
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType u=0;
      for(; ki<kshell[ks+1]; ki++)
        u += ((*rk1_r++)*(*rk2_r++)+(*rk1_i++)*(*rk2_i++));
      vk += Fk_symm[ks]*u;
    }
    return vk;
  }

  /** Evaluate the long-range potential with the open BC for the D-1 direction */
  virtual  mRealType evaluate_slab(pRealType z, const std::vector<int>& kshell
                                  , const pComplexType* restrict eikr_i, const pComplexType* restrict eikr_j)
  {
    return 0.0;
  }

  inline mRealType evaluate(const std::vector<int>& kshell
                           , int iat, const pComplexType* restrict rk2, ParticleSet &P)
  {
    mRealType vk=0.0;
#if !defined(USE_REAL_STRUCT_FACTOR)
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType u=0;
      for(; ki<kshell[ks+1]; ki++,rk2++)
      {
        pComplexType eikr=P.SK->eikr(iat,ki);
        u += eikr.real()*(*rk2).real()+eikr.imag()*(*rk2).imag();
      }
      vk += Fk_symm[ks]*u;
    }
#endif
    return vk;
  }

  inline mRealType evaluate(const std::vector<int>& kshell
                           , int iat, const pRealType* restrict rk2_r, const pRealType* restrict rk2_i, ParticleSet &P)
  {
    mRealType vk=0.0;
#if defined(USE_REAL_STRUCT_FACTOR)
    const pRealType* restrict eikr_r=P.SK->eikr_r[iat];
    const pRealType* restrict eikr_i=P.SK->eikr_i[iat];
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType u=0;
      for(; ki<kshell[ks+1]; ki++)
        u += eikr_r[ki]*(*rk2_r++)+eikr_i[ki]*(*rk2_i++);
      vk += Fk_symm[ks]*u;
    }
#endif
    return vk;
  }

  /** evaluate \f$\sum_k F_{k} \rho^1_{-{\bf k} \rho^2_{\bf k}\f$
   * and \f$\sum_k F_{k} \rho^1_{-{\bf k} \rho^2_{\bf k}\f$
   * @param kshell degeneracies of the vectors
   * @param rk1 starting address of \f$\rho^1_{{\bf k}\f$
   * @param rk2 starting address of \f$\rho^2_{{\bf k}\f$
   *
   * Valid for the strictly ordered k and \f$F_{k}\f$.
   */
  inline void evaluateGrad(const ParticleSet &A, const ParticleSet &B,
                           int specB, std::vector<pRealType> &Zat,
                           std::vector<TinyVector<pRealType,OHMMS_DIM> > &grad1)
  {
#if !defined(USE_REAL_STRUCT_FACTOR)
    const Matrix<pComplexType> &e2ikrA = A.SK->eikr;
    const pComplexType *rhokB = B.SK->rhok[specB];
    const std::vector<PosType> &kpts = A.SK->KLists.kpts_cart;
    for (int ki=0; ki<Fk.size(); ki++)
    {
      PosType k = kpts[ki];
      for (int iat=0; iat<Zat.size(); iat++)
      {
        grad1[iat] -= Zat[iat]*k*Fkg[ki]*
                      (e2ikrA(iat,ki).real()*rhokB[ki].imag() - e2ikrA(iat,ki).imag()*rhokB[ki].real());
      }
    }
#else
    
    const Matrix<pRealType> &e2ikrA_r = A.SK->eikr_r;
    const Matrix<pRealType> &e2ikrA_i = A.SK->eikr_i;
    const pRealType *rhokB_r = B.SK->rhok_r[specB];
    const pRealType *rhokB_i = B.SK->rhok_i[specB];
    const std::vector<PosType> &kpts = A.SK->KLists.kpts_cart;
    for (int ki=0; ki<Fk.size(); ki++)
    {
      PosType k = kpts[ki];
      for (int iat=0; iat<Zat.size(); iat++)
      {
        grad1[iat] -= Zat[iat]*k*Fkg[ki]*
                      (e2ikrA_r(iat,ki)*rhokB_i[ki] - e2ikrA_i(iat,ki)*rhokB_r[ki]);
      }
    }
#endif
  }
  
  ///FIX_PRECISION
  inline SymTensor<pRealType, OHMMS_DIM> evaluateStress(const std::vector<int>& kshell, 
      const pRealType *rhokA_r, const pRealType *rhokA_i, 
      const pRealType *rhokB_r, const pRealType *rhokB_i)
  { 
    SymTensor<pRealType, OHMMS_DIM> stress;
    for (int ki=0; ki<dFk_dstrain.size(); ki++)
    { 
      stress += (rhokA_r[ki]*rhokB_r[ki] + rhokA_i[ki]*rhokB_i[ki]) * dFk_dstrain[ki];
    }

    return stress;
  }  
  
  ///FIX_PRECISION
  inline SymTensor<pRealType, OHMMS_DIM> evaluateStress(const std::vector<int>& kshell, 
      const pComplexType * rhokA, const pComplexType * rhokB)
  { 
    SymTensor<pRealType, OHMMS_DIM> stress;
    for (int ki=0; ki<dFk_dstrain.size(); ki++)
    {
      stress += (rhokA[ki].real()*rhokB[ki].real() + rhokA[ki].imag()*rhokB[ki].imag()) * dFk_dstrain[ki];
    }
    return stress;
  }  

  /** evaluate \f$ v_{s}(k=0) = \frac{4\pi}{V}\int_0^{r_c} r^2 v_s(r) dr \f$
   */
  virtual mRealType evaluateSR_k0()
  {
    return 0.0;
  }
  /** evaluate \f$ v_s(r=0) \f$ for the self-interaction term
   */
  virtual mRealType evaluateLR_r0()
  {
    return 0.0;
  }
  
  ///These functions return the strain derivatives of all corresponding quantities
  /// in total energy.  See documentation (forthcoming).  
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateLR_r0_dstrain(){return 0;};
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateSR_k0_dstrain(){return 0;};
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateLR_dstrain(TinyVector<pRealType, OHMMS_DIM> k, pRealType kmag){return 0;};
  virtual SymTensor<mRealType, OHMMS_DIM> evaluateSR_dstrain(TinyVector<pRealType, OHMMS_DIM> r, pRealType rmag){return 0;};

  virtual void initBreakup(ParticleSet& ref)=0;
  virtual void Breakup(ParticleSet& ref, mRealType rs_in)=0;
  virtual void resetTargetParticleSet(ParticleSet& ref)=0;

  virtual mRealType evaluate(mRealType r, mRealType rinv)=0;
  virtual mRealType evaluateLR(mRealType r)=0;
  virtual mRealType srDf(mRealType r, mRealType rinv)=0;

  /** make clone */
  virtual LRHandlerBase* makeClone(ParticleSet& ref)=0;

};

/** LRHandler without breakup.
 *
 * The template parameter Func should impelement operator()(kk) which
 * returns the Fourier component of a long-range function. Here kk
 * is \f$|{\bf k}|^2\f$.
 */
template<class Func>
struct DummyLRHandler: public LRHandlerBase
{

  Func myFunc;

  DummyLRHandler(mRealType kc=-1.0): LRHandlerBase(kc)
  {}

  DummyLRHandler(const DummyLRHandler& aLR): LRHandlerBase(aLR)
  {}

  void initBreakup(ParticleSet& ref)
  {
    mRealType kcsq=LR_kc*LR_kc;
    KContainer& KList(ref.SK->KLists);
    int maxshell=KList.kshell.size()-1;
    const KContainer::SContainer_t& kk(KList.ksq);
    int ksh=0,ik=0;
    while(ksh<maxshell)
    {
      if(kk[ik]>kcsq)
        break; //exit
      ik=KList.kshell[++ksh];
    }
    MaxKshell=ksh;
    Fk_symm.resize(MaxKshell);
    Fk.resize(KList.kpts_cart.size());
    mRealType u0 = 4.0*M_PI/ref.Lattice.Volume;
    for(ksh=0,ik=0; ksh<MaxKshell; ksh++, ik++)
    {
      mRealType v=u0*myFunc(kk[ik]);//rpa=u0/kk[ik];
      Fk_symm[ksh]=v;
      for(; ik<KList.kshell[ksh+1]; ik++)
        Fk[ik]=v;
    }
  }

  mRealType evaluate(mRealType r, mRealType rinv)
  {
    return 0.0;
  }
  mRealType evaluateLR(mRealType r)
  {
    return 0.0;
  }
  mRealType srDf(mRealType r, mRealType rinv)
  {
    return 0.0;
  }
};

}
#endif
