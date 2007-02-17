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
    string ID_Rs;
    vector<FuncType*> FuncList;

    ~RPAConstraints();

    RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    void addOptimizables(OptimizableSetType& outVars); 
    /** update the optimizable variables
     */
    void resetParameters(OptimizableSetType& optVariables);


    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);
    void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
    bool put(xmlNodePtr cur);
  };
     
  struct RPAPBCConstraints: public OrbitalConstraintsBase {
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    
    bool IgnoreSpin;
    RealType Rs;
    string ID_Rs;

    RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    ~RPAPBCConstraints();
    void addOptimizables(OptimizableSetType& outVars);
    void resetParameters(OptimizableSetType& optVariables);
    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);
    void addExtra2ComboOrbital(ComboOrbital* jcombo);

    bool put(xmlNodePtr cur);
    OrbitalBase* createSRTwoBody();
    OrbitalBase* createLRTwoBody();
     
  private:
  };

  template<class T>
  struct ShortRangePartAdapter : OptimizableFunctorBase<T> {
  private:
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    HandlerType* handler;
  public:
    typedef typename OptimizableFunctorBase<T>::real_type real_type;  
    typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;  
    
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
    void resetParameters(OptimizableSetType& optVariables) { }
    bool put(xmlNodePtr cur) {return true;}
    void addOptimizables(OptimizableSetType& vlist){}
  };


  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
