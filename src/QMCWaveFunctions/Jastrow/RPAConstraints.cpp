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
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRHandlerTemp.h"
//#include "QMCWaveFunctions/DiffOrbitalBase.h"
//#include "LongRange/DummyLRHandler.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  /** Functor which return \f$frac{Rs}{k^2 (k^2+(1/Rs)^2)}\f$
   */
  template<typename T>
    struct RPA0
    {
      T Rs;
      T OneOverRsSq;
      RPA0(T rs=1):Rs(rs){OneOverRsSq=1.0/(Rs*Rs);}
      inline T operator()(T kk) 
      {
        T k2=std::sqrt(kk);
        return Rs/(k2*(k2+OneOverRsSq));
        //return (-0.5+0.5*std::pow(1.0+12.0*OneOverRs3/kk/kk,0.5));
      }
    };


  ///initialize the static member
  map<string,RPAPBCConstraints::HandlerType*> RPAPBCConstraints::handlerSet;

  ////////////////////////////////////////
  //RPAPBCConstraints definitions
  ////////////////////////////////////////
  RPAPBCConstraints::RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),
  IgnoreSpin(nospin), DropLongRange(false), DropShortRange(false), OwnHandler(false),
  LongRangeForm(0),
  Rs(-1.0), Kc(-1.0), Rcut(0.0), ID_Rs("rs"),
  myHandler(0), realHandler(0), kspaceHandler(0), 
  LongRangeRPA(0), ShortRangeRPA(0)
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

    //capture attribute jastrow/@name
    MyName="Jee";
    OhmmsAttributeSet a;
    a.add(MyName,"name");
    a.put(cur);

    //xmlNodePtr tcur = cur->children;
    string useL="yes";
    string useS="yes";
    ParameterSet params;
    string rpafunc="breakup";
    params.add(useL,"longrange","string");
    params.add(useS,"shortrange","string");
    params.add(Kc,"kc","double");
    params.add(rpafunc,"function","string");
    params.put(cur);

    //THIS IS NOT GOOD for consistency.
    //set the form of long-range term
    if(rpafunc == "breakup") 
      LongRangeForm=USE_BREAKUP;
    else if(rpafunc == "rpa") 
      LongRangeForm=USE_RPA;
    
    app_log() <<endl<<"   LongRangeForm is "<<rpafunc<<endl;
    
    DropLongRange = (useL == "no");
    DropShortRange = (useS=="no");

    bool success=getVariables(cur);

    //<parameter name="rs" id="rs_id"></parameter>
    map<string,pair<string,RealType> >::iterator vit(inVars.find("rs"));
    if(vit == inVars.end()) 
    {
      ID_Rs="rs"; 
      Rs=-1.0;
    }
    else
    {
      ID_Rs=(*vit).second.first; 
      Rs=(*vit).second.second;
    }

    app_log() << "    Rs can be optimized using ID=" << ID_Rs << endl;
    RealType tlen = std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/ static_cast<RealType>(targetPtcl.getTotalNum()) ,1.0/3.0);
    
    if(Rs<0) {
      if(targetPtcl.Lattice.SuperCellEnum) {
        Rs=tlen;
      } else {
        Rs=1.0;
      }
    }
    if(Kc<0) Kc=1.0/std::sqrt(tlen);

    //do not this again
    if(myHandler) return true;

    map<string,HandlerType*>::iterator hit(handlerSet.find(MyName));
    if(hit != handlerSet.end())//reuse the handler
    {
      OwnHandler=false;
      myHandler=realHandler=(*hit).second;
      return true;
    }

    app_log() << "    RPAPBCConstraints::addTwoBodyPart Rs = " << Rs <<  "  Kc= " << Kc << endl;

    OwnHandler=true;
    realHandler= new LRHandlerTemp<RPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    realHandler->Breakup(targetPtcl,Rs);

    //For derivatives, commented out 2008-04-07
    //derivHandler= new LRHandlerTemp<DerivRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    //derivHandler->Breakup(targetPtcl,Rs);

    //assign realHandler to myHandler
    myHandler=realHandler;

    //add the handler to the set so that we can reuse it
    handlerSet[MyName]=myHandler;
    string dname="d"+MyName;
    //handlerSet[dname]=derivHandler;

    //JK 2008-01-08: Disable overwriting kspace handler
    ////show how to overwrite  the kspace
    //if(LongRangeForm!=USE_BREAKUP && kspaceHandler==0)
    //{
    //  if(LongRangeForm == USE_RPA)
    //  { //example of rpa
    //    std::cout<<endl<<"   using RPA for the k-space part"<<endl;
    //    DummyLRHandler<RPA0<RealType> > *dummy= new DummyLRHandler<RPA0<RealType> >(Kc);
    //    dummy->myFunc=RPA0<RealType>(Rs);
    //    dummy->initBreakup(targetPtcl);
    //    myHandler=kspaceHandler=dummy;
    //  }
    //}

    return true;
  }

  void RPAPBCConstraints::addOptimizables(OptimizableSetType& outVars) {
    //outVars.addVariable(ID_Rs,Rs);
    //potentially add Rcut
    //outVars[ID_Rs]=Rs;
    if(LongRangeRPA) LongRangeRPA->put(NULL,outVars);
  }

  void RPAPBCConstraints::resetParameters(OptimizableSetType& optVariables) 
  { 
    
    OptimizableSetType::iterator it(optVariables.find(ID_Rs));
    if(it != optVariables.end()){ 
      //targetPtcl.Lattice.LR_rc *= ((*it).second) / rs;
      Rs=(*it).second;
      //std::cout<<"resetting a parameter: "<<(*it).first<<" to: "<<(*it).second<<endl;
      //myHandler->Breakup(targetPtcl,(*it).second);
      if(OwnHandler) realHandler->Breakup(targetPtcl,Rs);
      //JK: Don't know why we need this
      //myHandler=realHandler;

      //2008-04-07 Longrange is not committed yet.
      //realHandler->resetTargetParticleSet(targetPtcl,Rs);
      if(LongRangeRPA) LongRangeRPA->resetParameters(optVariables);
    
      //reset the numerical functor
      //nfunc->initialize(SRA, myGrid);

#if !defined(HAVE_MPI)
      ofstream fout("rpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
          << realHandler->evaluate(r,1.0/r) << " " 
          << realHandler->evaluateLR(r) << endl;
      }
#endif
    }
  }

  /** create TwoBody Jastrow using LR breakup method
   */
  OrbitalBase* RPAPBCConstraints::createTwoBody() 
  {
    if(DropShortRange)
    {
      app_log() << "    Disable Short-Range RPA. Return the Long-Range RPA." << endl;
      if(LongRangeRPA==0) 
        LongRangeRPA = new LRTwoBodyJastrow(targetPtcl, myHandler);
      return LongRangeRPA;
    }
    else
    {
      //return existing ShortRangeRPA
      if(ShortRangeRPA) return ShortRangeRPA;

      //short-range uses realHandler
      Rcut = realHandler->get_rc()-0.1;
      myGrid = new GridType;
      int npts=static_cast<int>(Rcut/0.01)+1;
      myGrid->set(0,Rcut,npts);

      //create the numerical functor
      nfunc = new FuncType;
      SRA = new ShortRangePartAdapter<RealType>(myHandler);
      SRA->setRmax(Rcut);
      nfunc->initialize(SRA, myGrid);

#if !defined(HAVE_MPI)
      static  int counter=0;
      char fname[32];
      sprintf(fname,"%s.%d.dat",MyName.c_str(),counter++);
      ofstream fout(fname);
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
          << realHandler->evaluate(r,1.0/r) << " " 
          << realHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType> *j2 = new TwoBodyJastrowOrbital<FuncType>(targetPtcl);
      j2->addFunc("rpa",0,0,nfunc);

    //For derivatives, commented out 2008-04-07
//      //create the numerical functor for the derivatives
//      ndfunc = new FuncType;
//      dSRA = new ShortRangePartAdapter<RealType>(derivHandler);
//      dSRA->setRmax(Rcut);
//      ndfunc->initialize(dSRA, myGrid);

//      TwoBodyJastrowOrbital<FuncType> *dj2 = new TwoBodyJastrowOrbital<FuncType>(targetPtcl);
//      dj2->addFunc("rpa",0,0,ndfunc);
//
//#if defined(ENABLE_SMARTPOINTER)
//      if(!dPsi.get())
//        dPsi= DiffOrbitalBasePtr(new AnalyticDiffOrbital(targetPtcl.getTotalNum(),dj2));
//#else
//      if(dPsi==0) 
//        dPsi= new AnalyticDiffOrbital(targetPtcl.getTotalNum(),dj2);
//#endif
//      else
//        dPsi->addOrbital(dj2);

//      LRTwoBodyJastrow *dLong=new LRTwoBodyJastrow(targetPtcl,derivHandler);
//      dLong->put(NULL,targetPsi.VarList);//resize
//      dPsi->addOrbital(dLong);
//
//      //dPsi->resetTargetParticleSet(targetPtcl);
//      int irpa=targetPsi.VarList.addVariable(ID_Rs,Rs);
//      dPsi->setBounds(irpa);

#if !defined(HAVE_MPI)
      ofstream dfout("drpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        //dfout << r << "   " << ndfunc->evaluate(r) << "   "
          //<< derivHandler->evaluate(r,1.0/r) << " " 
          //<< derivHandler->evaluateLR(r) << " "
          //<< nfunc->evaluate(r)-ndfunc->evaluate(r)
          //<< endl;
      }
#endif

      return ShortRangeRPA=j2;
    }
  }

  void RPAPBCConstraints::addExtra2ComboOrbital(ComboOrbital* jcombo)
  {
    //if DropShortRange is true, LongRangeRPA is already used
    if(DropShortRange) return;
    if(DropLongRange) return;

    if(LongRangeRPA==0)
      LongRangeRPA = new LRTwoBodyJastrow(targetPtcl, myHandler);
    jcombo->Psi.push_back(LongRangeRPA);
  }


  OrbitalBase* RPAPBCConstraints::createOneBody(ParticleSet& source) 
  {
    return 0;
  }

 
//  ////////////////////////////////////////
//  //RPAConstraints definitions
//  ////////////////////////////////////////
//  RPAConstraints::RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
//    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
//    {
//      JComponent.set(TWOBODY);
//    }
//
//  RPAConstraints::~RPAConstraints() {
//    delete_iter(FuncList.begin(), FuncList.end());
//  }
//
//  bool RPAConstraints::put(xmlNodePtr cur) {
//    bool success=getVariables(cur);
//    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
//    if(vit == inVars.end()) {
//      vit=inVars.find("rs");
//      if(vit == inVars.end()) {
//        ID_Rs="rs"; Rs=-1.0;
//      }
//    }
//    else 
//    {
//      ID_Rs=(*vit).second.first; 
//      Rs=(*vit).second.second;
//    }
//    return true;
//  }
//
//  void RPAConstraints::addOptimizables(OptimizableSetType& outVars) 
//  {
//    outVars[ID_Rs]=Rs;
//  }
//
//  void RPAConstraints::resetParameters(OptimizableSetType& optVariables) 
//  {
//    OptimizableSetType::iterator it(optVariables.find(ID_Rs));
//    if(it != optVariables.end()) 
//    {
//      Rs=(*it).second;
//      for(int i=0; i<FuncList.size(); i++) {
//        FuncList[i]->reset(Rs);
//      }
//    }
//  }
//
//  OrbitalBase* RPAConstraints::createTwoBody() {
//
//    if(Rs<0) {
//      if(targetPtcl.Lattice.SuperCellEnum) {
//        Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/static_cast<RealType>(targetPtcl.getTotalNum()),1.0/3.0);
//      } else { 
//        Rs=1.0;
//      }
//    }
//
//    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
//    JeeType *J2 = new JeeType(targetPtcl);
////    if(IgnoreSpin) {
//      FuncType *func=new FuncType(false);
//      func->reset(Rs);
//      for(int i=0; i<4; i++) {
//        J2->addFunc(func);
//      }
//      FuncList.push_back(func);
//      app_log() << "  RPAConstraints::Adding Spin-independent RPA Two-Body Jastrow " << endl;
//      app_log() << "    Rs=" << Rs << "  F=" << 1.0/func->B << endl;
////    } else {
////      FuncType *funcUU=new FuncType(true);
////      FuncType *funcUD=new FuncType(false);
////      funcUU->reset(Rs);
////      funcUD->reset(Rs);
////
////      J2->addFunc(funcUU);
////      J2->addFunc(funcUD);
////      J2->addFunc(funcUD);
////      J2->addFunc(funcUU);
////
////      FuncList.push_back(funcUU);
////      FuncList.push_back(funcUD);
////      app_log() << "  RPAConstraints::Adding Spin-dependent RPA Two-Body Jastrow " << endl;
////      app_log() << "    Rs=" << Rs << "  F(uu)=" << 1.0/funcUU->B << "  F(ud)= " << 1.0/funcUD->B << endl;
////    }
//    return J2;
//  }
//
//  OrbitalBase* RPAConstraints::createOneBody(ParticleSet& source) {
//    return 0;
//  }
//
     
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
