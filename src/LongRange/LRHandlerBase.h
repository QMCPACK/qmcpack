//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file LRHandlerBase.h
 * @brief Define LRHandlerBase and DummyLRHandler<typename Func>
 */
#ifndef QMCPLUSPLUS_LRHANLDERBASE_AND_DUMMY_H
#define QMCPLUSPLUS_LRHANLDERBASE_AND_DUMMY_H

#include "LongRange/StructFact.h"
namespace qmcplusplus {

  /** base class for LRHandlerTemp<FUNC,BASIS> and DummyLRHanlder<typename Func>
   */
  struct LRHandlerBase: public QMCTraits {

    /// Maxkimum Kshell for the given Kc
    IndexType MaxKshell;
    /// Maximum k cutoff 
    RealType  LR_kc;
    /// Maximum r cutoff 
    RealType  LR_rc;
    ///Fourier component for all the k-point
    Vector<RealType> Fk; 
    ///Fourier component for each k-shell
    Vector<RealType> Fk_symm; 
    ///Coefficient
    Vector<RealType> coefs; 

    //constructor
    LRHandlerBase(RealType kc):LR_kc(kc){}
    //return r cutoff
    inline RealType get_rc() const { return LR_rc;}
    //return k cutoff
    inline RealType get_kc() const { return LR_kc;}

    /** evaluate \f$\sum_k F_{k} \rho^1_{-{\bf k} \rho^2_{\bf k}\f$
     * @param kshell degeneracies of the vectors
     * @param rk1 starting address of \f$\rho^1_{{\bf k}\f$
     * @param rk2 starting address of \f$\rho^2_{{\bf k}\f$
     * 
     * Valid for the strictly ordered k and \f$F_{k}\f$. 
     */
    inline RealType evaluate(const vector<int>& kshell
        , const ComplexType* restrict rk1, const ComplexType* restrict rk2) 
    {
      RealType vk=0.0;
      for(int ks=0,ki=0; ks<MaxKshell; ks++) {
	RealType u=0;
	for(;ki<kshell[ks+1]; ki++,rk1++,rk2++)
	  u += ((*rk1).real()*(*rk2).real()+(*rk1).imag()*(*rk2).imag());
        vk += Fk_symm[ks]*u;

      }
      //for(int ki=0; ki<Fk.size(); ki++) {
      //  //vk += (rk1[ki]*rk2[minusk[ki]]).real()*Fk[ki];
      //  vk += (rk1[ki].real()*rk2[ki].real()+rk1[ki].imag()*rk2[ki].imag())*Fk[ki];
      //} //ki
      return vk;
    }

    /** Evaluate the long-range potential with the open BC for the D-1 direction */
    virtual  RealType evaluate_slab(RealType z, const vector<int>& kshell
        , const ComplexType* restrict eikr_i, const ComplexType* restrict eikr_j) 
    {
      return 0.0;
    }
    
    inline RealType evaluate(const vector<int>& kshell
        , int iat, const ComplexType* restrict rk2, ParticleSet &P) 
    {
      RealType vk=0.0;
      for(int ks=0,ki=0; ks<MaxKshell; ks++) {
        RealType u=0;
        for(;ki<kshell[ks+1]; ki++,rk2++){
          ComplexType eikr=P.SK->eikr(iat,ki);
          u += eikr.real()*(*rk2).real()+eikr.imag()*(*rk2).imag();
        }
        vk += Fk_symm[ks]*u;
      }
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
         int specB, vector<RealType> &Zat,
         vector<TinyVector<RealType,DIM> > &grad1) 
     {
       const Matrix<ComplexType> &e2ikrA = A.SK->eikr;
       const ComplexType *rhokB = B.SK->rhok[specB];
       const vector<PosType> &kpts = A.SK->KLists.kpts_cart;
       for (int ki=0; ki<Fk.size(); ki++) 
       {
         TinyVector<RealType,DIM> k = kpts[ki];
         for (int iat=0; iat<Zat.size(); iat++) 
         {
           grad1[iat] -= Zat[iat]*k*Fk[ki]*
             (e2ikrA(iat,ki).real()*rhokB[ki].imag() - e2ikrA(iat,ki).imag()*rhokB[ki].real());
         }
       }
     }

    /** evaluate \f$ v_{s}(k=0) = \frac{4\pi}{V}\int_0^{r_c} r^2 v_s(r) dr \f$
     */
    virtual RealType evaluateSR_k0() { return 0.0;}
    /** evaluate \f$ v_s(r=0) \f$ for the self-interaction term
     */
    virtual RealType evaluateLR_r0() { return 0.0;}

    virtual void initBreakup(ParticleSet& ref)=0;
    virtual void Breakup(ParticleSet& ref, RealType rs_in)=0;
    virtual void resetTargetParticleSet(ParticleSet& ref)=0;

    virtual RealType evaluate(RealType r, RealType rinv)=0;
    virtual RealType evaluateLR(RealType r)=0;
    virtual RealType srDf(RealType r, RealType rinv)=0;

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
    struct DummyLRHandler: public LRHandlerBase {

      Func myFunc;

      DummyLRHandler(RealType kc=-1.0): LRHandlerBase(kc)
        {}

      DummyLRHandler(const DummyLRHandler& aLR): LRHandlerBase(aLR)
        {}

      void initBreakup(ParticleSet& ref)
      {
        RealType kcsq=LR_kc*LR_kc;

        KContainer& KList(ref.SK->KLists);
        int maxshell=KList.kshell.size()-1;
        const KContainer::SContainer_t& kk(KList.ksq);

        int ksh=0,ik=0;
        while(ksh<maxshell)
        {
          if(kk[ik]>kcsq) break; //exit
          ik=KList.kshell[++ksh];
        }
        MaxKshell=ksh;

        Fk_symm.resize(MaxKshell);
        Fk.resize(KList.kpts_cart.size());

        RealType u0 = 4.0*M_PI/ref.Lattice.Volume;
        for(ksh=0,ik=0; ksh<MaxKshell; ksh++, ik++)
        {
          RealType v=u0*myFunc(kk[ik]);//rpa=u0/kk[ik];
          Fk_symm[ksh]=v;
          for(; ik<KList.kshell[ksh+1]; ik++) Fk[ik]=v;
        }
      }

      RealType evaluate(RealType r, RealType rinv) { return 0.0;}
      RealType evaluateLR(RealType r){return 0.0;}
      RealType srDf(RealType r, RealType rinv){ return 0.0;}
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2031 $   $Date: 2007-05-31 16:20:53 -0500 (Thu, 31 May 2007) $
 * $Id: LRHandlerTemp.h 2031 2007-05-31 21:20:53Z jnkim $
 ***************************************************************************/
