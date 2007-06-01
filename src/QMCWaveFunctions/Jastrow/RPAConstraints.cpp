//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#include "QMCWaveFunctions/Jastrow/RPAConstraints.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRJastrowSingleton.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus {

  template<class T>
  struct ShortRangePartAdapter : OptimizableFunctorBase<T> {
  private:
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    T Uconst;
    HandlerType* myHandler;
  public:
    typedef typename OptimizableFunctorBase<T>::real_type real_type;  
    typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;  
    
    explicit ShortRangePartAdapter(HandlerType* inhandler): Uconst(0) {
      myHandler = inhandler;
    }
    inline void setRmax(real_type rm) { Uconst=myHandler->evaluate(rm,1.0/rm);}
    inline real_type evaluate(real_type r) { return f(r); }
    inline real_type f(real_type r) 
    { 
      return myHandler->evaluate(r, 1.0/r)-Uconst; 
    }
    inline real_type df(real_type r) 
    {
      return myHandler->srDf(r, 1.0/r);
    }
    void resetParameters(OptimizableSetType& optVariables) { }
    bool put(xmlNodePtr cur) {return true;}
    void addOptimizables(OptimizableSetType& vlist){}
  };


  ////////////////////////////////////////
  //RPAConstraints definitions
  ////////////////////////////////////////
  RPAConstraints::RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
    {
      JComponent.set(TWOBODY);
    }

  RPAConstraints::~RPAConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool RPAConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
    if(vit == inVars.end()) {
      vit=inVars.find("rs");
      if(vit == inVars.end()) {
        ID_Rs="rs"; Rs=-1.0;
      }
    }
    else 
    {
      ID_Rs=(*vit).second.first; 
      Rs=(*vit).second.second;
    }
    return true;
  }

  void RPAConstraints::addOptimizables(OptimizableSetType& outVars) 
  {
    outVars[ID_Rs]=Rs;
  }

  void RPAConstraints::resetParameters(OptimizableSetType& optVariables) 
  {
    OptimizableSetType::iterator it(optVariables.find(ID_Rs));
    if(it != optVariables.end()) 
    {
      Rs=(*it).second;
      for(int i=0; i<FuncList.size(); i++) {
        FuncList[i]->reset(Rs);
      }
    }
  }

  OrbitalBase* RPAConstraints::createTwoBody() {

    if(Rs<0) {
      if(targetPtcl.Lattice.SuperCellEnum) {
        Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/static_cast<RealType>(targetPtcl.getTotalNum()),1.0/3.0);
      } else { 
        Rs=1.0;
      }
    }

    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(targetPtcl);
//    if(IgnoreSpin) {
      FuncType *func=new FuncType(false);
      func->reset(Rs);
      for(int i=0; i<4; i++) {
        J2->addFunc(func);
      }
      FuncList.push_back(func);
      app_log() << "  RPAConstraints::Adding Spin-independent RPA Two-Body Jastrow " << endl;
      app_log() << "    Rs=" << Rs << "  F=" << 1.0/func->B << endl;
//    } else {
//      FuncType *funcUU=new FuncType(true);
//      FuncType *funcUD=new FuncType(false);
//      funcUU->reset(Rs);
//      funcUD->reset(Rs);
//
//      J2->addFunc(funcUU);
//      J2->addFunc(funcUD);
//      J2->addFunc(funcUD);
//      J2->addFunc(funcUU);
//
//      FuncList.push_back(funcUU);
//      FuncList.push_back(funcUD);
//      app_log() << "  RPAConstraints::Adding Spin-dependent RPA Two-Body Jastrow " << endl;
//      app_log() << "    Rs=" << Rs << "  F(uu)=" << 1.0/funcUU->B << "  F(ud)= " << 1.0/funcUD->B << endl;
//    }
    return J2;
  }

  OrbitalBase* RPAConstraints::createOneBody(ParticleSet& source) {
    return 0;
  }

  ////////////////////////////////////////
  //RPAPBCConstraints definitions
  ////////////////////////////////////////
  RPAPBCConstraints::RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin), Rs(-1.0), Kc(-1.0),
  myHandler(0), LongRangeRPA(0)
    {
      JComponent.set(MULTIPLE);
      JComponent.set(TWOBODY);
      JComponent.set(LONGRANGE);
    }

  RPAPBCConstraints::~RPAPBCConstraints() 
  {  
    //may delete LongRangeRPA
  }

  bool RPAPBCConstraints::put(xmlNodePtr cur) {
    string useL="yes";
    string useS="yes";
    ParameterSet params;
    params.add(useL,"longrange","string");
    params.add(useS,"shortrange","string");
    params.add(Kc,"kc","double");
    params.put(cur);

    DropLongRange = (useL=="no");
    DropShortRange = (useS=="no");
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
    if(vit == inVars.end()) {
      vit=inVars.find("rs");
      if(vit == inVars.end()) {
        ID_Rs="rs"; Rs=-1.0;
      }
    }
    else
    {
      ID_Rs=(*vit).second.first; 
      Rs=(*vit).second.second;
    }
    return true;
  }

  void RPAPBCConstraints::addOptimizables(OptimizableSetType& outVars) {
    //potentially add Rcut
    outVars[ID_Rs]=Rs;
    if(LongRangeRPA) LongRangeRPA->put(NULL,outVars);
  }

  void RPAPBCConstraints::resetParameters(OptimizableSetType& optVariables) 
  { 
    OptimizableSetType::iterator it(optVariables.find(ID_Rs));
    if(it != optVariables.end()) Rs=(*it).second;
  }

  /** create TwoBody Jastrow using LR breakup method
   */
  OrbitalBase* RPAPBCConstraints::createTwoBody() 
  {
    if(Rs<0) {
      if(targetPtcl.Lattice.SuperCellEnum) {
        Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/static_cast<RealType>(targetPtcl.getTotalNum()),1.0/3.0);
      } else {
        Rs=1.0;
      }
    }
    if(Kc<0) Kc=1.0/std::sqrt(Rs);

    app_log() << "    RPAPBCConstraints::addTwoBodyPart Rs = " << Rs <<  "  Kc= " << Kc << endl;

    if(myHandler==0)
        myHandler = LRJastrowSingleton::getHandler(targetPtcl,Kc);

    if(DropShortRange)
    {
      DropLongRange=true;//prevent adding twice
      app_log() << "    Disable Short-Range RPA. Return the Long-Range RPA." << endl;
      if(LongRangeRPA==0)
      {
        LongRangeRPA = new LRTwoBodyJastrow(targetPtcl, myHandler);
      }
      return LongRangeRPA;
    }
    else
    {
      typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
      typedef CubicSplineSingle<RealType,SplineEngineType> FuncType;
      typedef LinearGrid<RealType> GridType;

      RealType Rcut = myHandler->Basis.get_rc()-0.1;
      myGrid = new GridType;
      int npts=static_cast<int>(Rcut/0.01)+1;
      myGrid->set(0,Rcut,npts);

      //create the numerical functor
      FuncType* nfunc = new FuncType;
      ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(myHandler);
      sra->setRmax(Rcut);
      nfunc->initialize(sra, myGrid);

#if !defined(HAVE_MPI)
      ofstream fout("rpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
          << myHandler->evaluate(r,1.0/r) << " " 
          << myHandler->evaluateLR(r) << endl;
      }
#endif

      TwoBodyJastrowOrbital<FuncType> *J2 = new TwoBodyJastrowOrbital<FuncType>(targetPtcl);
      for (int i=0; i<4; i++) J2->addFunc(nfunc);
      return J2;
    }
  }
     
  //OrbitalBase* RPAPBCConstraints::createLRTwoBody() {
  //  if(LongRangeRPA==0)
  //    LongRangeRPA = new LRTwoBodyJastrow(targetPtcl, myHandler);
  //  return LongRangeRPA;
  //}


  void RPAPBCConstraints::addExtra2ComboOrbital(ComboOrbital* jcombo)
  {
    if(DropLongRange)
    {
      app_log() << "    Disable Long-Range RPA Jastrow " << endl;
    }
    else
    {
      app_log() << "    Add Long-Range RPA Jastrow " << endl;
      if(LongRangeRPA==0)
      {
        //HandlerType* handler = LRJastrowSingleton::getHandler(targetPtcl,1.0/std::sqrt(Rs));
        LongRangeRPA = new LRTwoBodyJastrow(targetPtcl, myHandler);
      }
      jcombo->Psi.push_back(LongRangeRPA);
    }
  }
     
  OrbitalBase* RPAPBCConstraints::createOneBody(ParticleSet& source) {
    return 0;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
