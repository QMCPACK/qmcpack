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
    Vector<RealType> coefs; 
    Vector<RealType> Fk; 
    Vector<RealType> Fk_symm; 

    //constructor
    LRHandlerBase(RealType kc):LR_kc(kc){}
    //copy constructor
    LRHandlerBase(const LRHandlerBase& rh):
      MaxKshell(rh.MaxKshell), LR_kc(rh.LR_kc), LR_rc(rh.LR_rc), 
    coefs(rh.coefs), Fk(rh.Fk), Fk_symm(rh.Fk_symm)
    {}

    //return r cutoff
    inline RealType get_rc() const { return LR_rc;}
    //return k cutoff
    inline RealType get_kc() const { return LR_kc;}
    virtual void initBreakup(ParticleSet& ref)=0;
    virtual void Breakup(ParticleSet& ref, RealType rs_in)=0;

    virtual RealType evaluate(RealType r, RealType rinv)=0;
    virtual RealType evaluateLR(RealType r)=0;
    virtual RealType srDf(RealType r, RealType rinv)=0;
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
