//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#define QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus {

  struct RPAConstraints: public OrbitalConstraintsBase {
    typedef RPAJastrow<RealType> FuncType;
    bool IgnoreSpin;
    RealType Rs;
    string ID;
    vector<FuncType*> FuncList;

    ~RPAConstraints();

    RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    inline void apply() {
      for(int i=0; i<FuncList.size(); i++) {
        FuncList[i]->reset(Rs);
      }
    }

    void addOptimizables(VarRegistry<RealType>& outVars) {
      outVars.add(ID,&Rs,1);
    }

    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);
    inline void addTwoBodyPart(ComboOrbital* jcombo) {
      OrbitalBase* j2 = createTwoBody();
      if (j2) jcombo->Psi.push_back(j2);
    }
    bool put(xmlNodePtr cur);
  };
     
  struct RPAPBCConstraints: public OrbitalConstraintsBase {
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    
    bool IgnoreSpin;
    RealType Rs;
    string ID;

    RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    ~RPAPBCConstraints();
    void apply();
    void addOptimizables(VarRegistry<RealType>& outVars);
    OrbitalBase* createSRTwoBody();
    OrbitalBase* createLRTwoBody();
    void addTwoBodyPart(ComboOrbital* jcombo);
    OrbitalBase* createOneBody(ParticleSet& source);
    bool put(xmlNodePtr cur);
     
  private:
  };

  template<class T>
  struct ShortRangePartAdapter : OptimizableFunctorBase<T> {
  private:
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    HandlerType* handler;
  public:
    typedef typename OptimizableFunctorBase<T>::real_type real_type;  
    
    explicit ShortRangePartAdapter(HandlerType* inhandler) {
      handler = inhandler;
    }
    inline real_type evaluate(real_type r) { return f(r); }
    inline real_type f(real_type r) { return handler->evaluate(r, 1.0/r); }
    inline real_type df(real_type r) {
 /*
      // for now just do this numerically (very crudely)
      real_type dr = 1e-8;
      real_type fr = handler->evaluate(r, 1.0/r);
      real_type frPlusDr = handler->evaluate(r+dr, 1.0/(r+dr));
      return (frPlusDr - fr) / dr;
  */    
       
      return handler->srDf(r, 1.0/r);
    }
    void reset() { }
    bool put(xmlNodePtr cur) {return true;}
    void addOptimizables( VarRegistry<real_type>& vlist){}
  };


  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
