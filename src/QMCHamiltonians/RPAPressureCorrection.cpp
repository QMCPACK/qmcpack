//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"


#include "Configuration.h"
#include "QMCHamiltonians/RPAPressureCorrection.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "LongRange/LRHandlerBase.h" 
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/RPAConstraints.h"
#include "Optimize/VarList.h"
#include "ParticleBase/ParticleAttribOps.h"




namespace qmcplusplus {

  typedef QMCHamiltonianBase::Return_t Return_t;
    
  ///initialize the static member
  map<string,RPAPressureCorrection::HandlerType*> RPAPressureCorrection::handlerSet;
  
  void RPAPressureCorrection::resetTargetParticleSet(ParticleSet& P) { };
    
  Return_t RPAPressureCorrection::evaluate(ParticleSet& P) {
    vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
    tValue=0.0;
    Value=0.0;
    dG = 0.0;
    dL = 0.0;
    while(dit != dit_end) {
      tValue += (*dit)-> evaluateLog(P,dG,dL);
      ++dit;
    }
    
    ZVCorrection =  -1.0 * drsdV * (0.5*Sum(dL)+Dot(dG,P.G));

//     tValue is d_r log dPsi
    Value = ZVCorrection;
    
    return 0.0;
  }

  Return_t RPAPressureCorrection::evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
    return evaluate(P);
  }
    
    
  bool RPAPressureCorrection::put(xmlNodePtr cur, ParticleSet& P) {


    MyName = "ZVZB";
    xmlNodePtr tcur = cur->children;
      
    double RPAKCut= -1.0;
    Rs=-1.0;
    //alpha is the derivative of r_s^star with respect to r_s. It could be optimized.
    string drsstardrs("alpha");
    string RPAPCorr("ZB");
    string RPAPfunc("RPA_LR");
    ParameterSet nattrib;
    OhmmsAttributeSet attrib;
    attrib.add(RPAPCorr,"etype" );
    attrib.add(RPAPfunc,"functor" );
    attrib.put(cur);
    nattrib.add(RPAKCut,"kc","double");
    nattrib.add(Rs,"Rs","double");
    nattrib.add(drsStar,"alpha","double");
    nattrib.put(cur);
//       while(tcur != NULL) {
//         nattrib.put(tcur);
//         tcur=tcur->next;
//  
//     if (realHandler) return true;
    map<string,HandlerType*>::iterator hit(handlerSet.find(MyName));
    if(hit != handlerSet.end())//reuse the handler
    {
      OwnHandler=false;
      realHandler=(*hit).second;
//       return true;
    }
    
    
    
    
    RealType tlen=std::pow(0.75/M_PI*P.Lattice.Volume/static_cast<RealType>(P.getTotalNum()),1.0/3.0);
    
    if(Rs<0.0) Rs=tlen;
    
    if(RPAKCut<0.0) {
      RPAKCut=P.Lattice.LR_kc;
    };
    
    drsdV= tlen/(3.0* P.Lattice.Volume)*drsStar;

    
    app_log() <<"    RPAPressureCorrection::createDeriv "<<endl<<
                "       Rs = " << Rs <<endl<<
                "    dRs^* = " << drsStar <<endl<<
                "       Kc = " << RPAKCut << endl;
    OwnHandler=true;

    if (RPAPfunc=="RPA"){
      app_log() <<"      and long range RPA wavefunction derivative "<<endl;
      if (OwnHandler){
      realHandler= new LRRPAHandlerTemp<DerivRPABreakup<Return_t>,LPQHIBasis >(P,RPAKCut);
      realHandler-> Breakup(P,Rs);
      //add the handler to the set so that we can reuse it
      handlerSet[MyName]=realHandler;
      };
      LRTwoBodyJastrow* LongRangeRPA = new LRTwoBodyJastrow(P,realHandler);
//       LongRangeRPA->resetByHandler();

      VarRegistry<RealType> tempO;
      LongRangeRPA->put(0,tempO);
      dPsi.push_back(LongRangeRPA);

      //short-range uses realHandler
      RealType Rcut = realHandler->get_rc()-0.1;
      GridType* myGrid = new GridType;
      int npts=static_cast<int>(Rcut/0.01)+1;
      myGrid->set(0,Rcut,npts);

      //create the numerical functor
      FuncType* nfunc = new FuncType;
      ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(realHandler);
      sra->setRmax(Rcut);
      nfunc->initialize(sra, myGrid);
#if !defined(HAVE_MPI)
      ofstream fout("drpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
            << realHandler->evaluate(r,1.0/r) << " " 
            << realHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j2= new TwoBodyJastrowOrbital<FuncType>(P);
      /*      j2->init(P);*/
//       for (int i=0; i<4; i++) j2->addFunc(nfunc);
//       for (int i=0; i<4; i++) j2->addFunc("ZVZB",0,0,nfunc);
      j2->addFunc("ZVZB",0,0,nfunc);
      dPsi.push_back(j2);
    }
    else if (RPAPfunc=="BREAKUP_LR" ){
      app_log() <<"      and long range breakup derivative"<<endl;
      if (OwnHandler){

      realHandler= new LRHandlerTemp<DerivYukawaBreakup<Return_t>,LPQHIBasis >(P,RPAKCut);
      realHandler-> Breakup(P,Rs);
    //add the handler to the set so that we can reuse it
      handlerSet[MyName]=realHandler;
      };
      LRTwoBodyJastrow* LongRangeRPA = new LRTwoBodyJastrow(P,realHandler);
//       LongRangeRPA->resetByHandler();

      VarRegistry<RealType> tempO;
      LongRangeRPA->put(0,tempO);
      dPsi.push_back(LongRangeRPA);
    }
    else if (RPAPfunc=="NUMERIC"){
      app_log() <<"and numeric Derivative.  NOT IMPLEMENTED YET!!!"<<endl;
    }
    else if (RPAPfunc=="Yukawa"){
      app_log() <<"and short/long range Yukawa derivative"<<endl;
      if (OwnHandler){
      realHandler= new LRHandlerTemp<DerivYukawaBreakup<Return_t>,LPQHIBasis >(P,RPAKCut);
      realHandler-> Breakup(P,Rs);
    //add the handler to the set so that we can reuse it
      handlerSet[MyName]=realHandler;
      };
      LRTwoBodyJastrow* LongRangeRPA = new LRTwoBodyJastrow(P,realHandler);
//       LongRangeRPA->resetByHandler();

      VarRegistry<RealType> tempO;
      LongRangeRPA->put(0,tempO);
      dPsi.push_back(LongRangeRPA);
      
      


      //short-range uses realHandler
      RealType Rcut = realHandler->get_rc()-0.1;
      GridType* myGrid = new GridType;
      int npts=static_cast<int>(Rcut/0.01)+1;
      myGrid->set(0,Rcut,npts);

      //create the numerical functor
      FuncType* nfunc = new FuncType;
      ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(realHandler);
      sra->setRmax(Rcut);
      nfunc->initialize(sra, myGrid);
#if !defined(HAVE_MPI)
      ofstream fout("drpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
            << realHandler->evaluate(r,1.0/r) << " " 
            << realHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j2= new TwoBodyJastrowOrbital<FuncType>(P);
      /*      j2->init(P);*/
//       for (int i=0; i<4; i++) j2->addFunc(nfunc);
//       for (int i=0; i<4; i++) j2->addFunc("ZVZB",0,0,nfunc);
      j2->addFunc("ZVZB",0,0,nfunc);
      dPsi.push_back(j2);
    }
    else if (RPAPfunc=="Yukawa_TEST"){
      app_log() <<"and TESTING!!!  "<<endl;
      HandlerType* realHandler;
      realHandler= new LRHandlerTemp<DerivYukawaBreakup<Return_t>,LPQHIBasis >(P,RPAKCut);
      realHandler-> initBreakup(P);

      LRTwoBodyJastrow* LongRangeRPA = new LRTwoBodyJastrow(P,realHandler);
//       LongRangeRPA->resetByHandler();

      VarRegistry<RealType> tempO;
      LongRangeRPA->put(0,tempO);
        dPsi.push_back(LongRangeRPA);
//       double dPsiKVal = (*dPsi.begin())->evaluateLog(P,P.G,P.L);
//       cout<<"point-1 "<<(*dPsi.begin())->evaluateLog(P,P.G,P.L)<<endl;



      //short-range uses realHandler
      RealType Rcut = realHandler->get_rc()-0.1;
      GridType* myGrid = new GridType;
      int npts=static_cast<int>(Rcut/0.01)+1;
      myGrid->set(0,Rcut,npts);

      //create the numerical functor
      FuncType* nfunc = new FuncType;
      ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(realHandler);
      sra->setRmax(Rcut);
      nfunc->initialize(sra, myGrid);
#if !defined(HAVE_MPI)
      ofstream fout("drpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
            << realHandler->evaluate(r,1.0/r) << " " 
            << realHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j2= new TwoBodyJastrowOrbital<FuncType>(P);
/*      j2->init(P);*/
//       for (int i=0; i<4; i++) j2->addFunc(nfunc);
      j2->addFunc("ZVZB",0,0,nfunc);
      dPsi.push_back(j2);
      
      double tValue=0.0;
      vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
      while(dit != dit_end) {
        tValue += (*dit)-> evaluateLog(P,P.G,P.L);
        std::cout<<"  "<<tValue<<endl;
        ++dit;
      }
//       double dPsiRVal = tValue-dPsiKVal;
      cout<<"point0 "<<tValue<<endl;
      
      ParticleSet* tempP = new ParticleSet(P);

      
      HandlerType* tempHandler= new LRHandlerTemp<YukawaBreakup<Return_t>,LPQHIBasis >(*tempP,RPAKCut);
      tempHandler-> initBreakup(*tempP);
      LRTwoBodyJastrow* tempLongRangeRPA = new LRTwoBodyJastrow(*tempP,tempHandler);
      tempLongRangeRPA->put(0,tempO);
      
      vector<OrbitalBase*> tempPsi;
        tempPsi.push_back(tempLongRangeRPA);
      
//       double tempPsiKVal = (*tempPsi.begin())->evaluateLog(*tempP,(*tempP).G,(*tempP).L);

      
      //short-range uses realHandler
      RealType Rcut2 = tempHandler->get_rc()-0.1;
      GridType* myGrid2 = new GridType;
      int npts2=static_cast<int>(Rcut2/0.01)+1;
      myGrid2->set(0,Rcut2,npts2);

      //create the numerical functor
      FuncType* nfunc2 = new FuncType;
      ShortRangePartAdapter<RealType>* sra2 = new ShortRangePartAdapter<RealType>(tempHandler);
      sra2->setRmax(Rcut2);
      nfunc2->initialize(sra2, myGrid2);
#if !defined(HAVE_MPI)
      ofstream fout2("drpa2.short.dat");
      for (int i = 0; i < myGrid2->size(); i++) {
        RealType r=(*myGrid2)(i);
        fout2 << r << "   " << nfunc2->evaluate(r) << "   "
            << tempHandler->evaluate(r,1.0/r) << " " 
            << tempHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j3= new TwoBodyJastrowOrbital<FuncType>(*tempP);
/*      j3->init(*tempP);*/
//       for (int i=0; i<4; i++) j3->addFunc(nfunc2);
      j3->addFunc("ZVZB",0,0,nfunc2);
      tempPsi.push_back(j3);
      
      cout<<"point0.5 "<<endl;
      double mValue=0.0;
      vector<OrbitalBase*>::iterator mit(tempPsi.begin()), mit_end(tempPsi.end());
      while(mit != mit_end) {
        mValue += (*mit)-> evaluateLog(*tempP,(*tempP).G,(*tempP).L);
        ++mit;
      }
      
//       double tempPsiRVal = mValue-tempPsiKVal;
//       double lowV = (*tempPsi.begin())->evaluateLog(*tempP,(*tempP).G,(*tempP).L);
      cout<<"point1 "<< mValue <<endl;
      
      ParticleSet* tempP2 = new ParticleSet(P);
      /*      (*tempP2).initParticleSet();*/
      P.convert2Unit(P.R,(*tempP).R);
      double scale = 1.0005;
      //RPAKCut=1.0/(std::sqrt(Rs)*scale); 
      (*tempP2).Lattice*=scale;
      (*tempP2).convert2Cart((*tempP).R,(*tempP2).R);
      (*tempP2).Lattice.SetLRCutoffs();
      (*tempP2).SK->UpdateNewCell((*tempP2).Lattice.LR_kc);
      (*tempP2).update();
      //(*tempP2).Lattice.LR_rc=P.Lattice.LR_rc;
      
      (*tempP2).Lattice.print(std::cout);
      std::cout<<(*tempP2).R<<endl;
      std::cout<<(*tempP2).Lattice.LR_kc<<endl;
      std::cout<<(*tempP).Lattice.Volume<<" <-p1,p2-> "<<(*tempP2).Lattice.Volume<<endl;
      std::cout<<" R_s for 1: "<<std::pow(3.0/4.0/M_PI*(*tempP).Lattice.Volume/static_cast<RealType>((*tempP).getTotalNum()),1.0/3.0)<<endl;
      std::cout<<" R_s for 2: "<<std::pow(3.0/4.0/M_PI*(*tempP2).Lattice.Volume/static_cast<RealType>((*tempP2).getTotalNum()),1.0/3.0)<<endl;
      std::cout<<Dot((*tempP2).R,(*tempP2).R)<<"  "<<Dot(P.R,P.R)<<endl;
          
      HandlerType* tempHandler2= new LRHandlerTemp<YukawaBreakup<Return_t>,LPQHIBasis >(*tempP2,RPAKCut);
      tempHandler2-> initBreakup(*tempP2);
      LRTwoBodyJastrow* tempLongRangeRPA2 = new LRTwoBodyJastrow(*tempP2,tempHandler2);
      tempLongRangeRPA2->put(0,tempO);
      
      vector<OrbitalBase*> tempPsi2;
      tempPsi2.push_back(tempLongRangeRPA2);
      
      //RealType Rcut3 = tempHandler2->get_rc()-0.1;
      RealType Rcut3 = Rcut2;
      GridType* myGrid3 = new GridType;
      int npts3=static_cast<int>(Rcut3/0.01)+1;
      myGrid3->set(0,Rcut3,npts3);

      //create the numerical functor
      FuncType* nfunc3 = new FuncType;
      ShortRangePartAdapter<RealType>* sra3 = new ShortRangePartAdapter<RealType>(tempHandler2);
      sra3->setRmax(Rcut3);
      nfunc3->initialize(sra3, myGrid3);
#if !defined(HAVE_MPI)
      ofstream fout3("drpa3.short.dat");
      for (int i = 0; i < myGrid3->size(); i++) {
        RealType r=(*myGrid3)(i);
        fout3 << r << "   " << nfunc3->evaluate(r) << "   "
            << tempHandler2->evaluate(r,1.0/r) << " " 
            << tempHandler2->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j4= new TwoBodyJastrowOrbital<FuncType>(*tempP2);
/*      j4->init(*tempP2);*/
//       for (int i=0; i<4; i++) j4->addFunc(nfunc3);
      j4->addFunc("ZVZB",0,0,nfunc3);
      tempPsi2.push_back(j4);
      
      cout<<"point1.5 "<<endl;
      double lValue=0.0;
      vector<OrbitalBase*>::iterator lit(tempPsi2.begin()), lit_end(tempPsi2.end());
      while(lit != lit_end) {
        lValue += (*lit)-> evaluateLog(*tempP2,(*tempP2).G,(*tempP2).L);
        ++lit;
      }
//       double lowV = (*tempPsi.begin())->evaluateLog(*tempP,(*tempP).G,(*tempP).L);
      cout<<"point1 "<< lValue <<endl;
      
      
      
/*      
      double highV = (*tempPsi2.begin())->evaluateLog(*tempP2,(*tempP2).G,(*tempP2).L);
      cout<<"point2 "<<highV<<endl;*/


      std::cout<<(lValue-mValue)/((scale-1.0)*Rs)<<" is numeric, "<<tValue<<" is analytic. "<<endl;
    }
    else if (RPAPfunc=="RPA_TEST"){
      app_log() <<"and TESTING!!!  "<<endl;
      HandlerType* realHandler;
      realHandler= new LRRPAHandlerTemp<DerivRPABreakup<Return_t>,LPQHIBasis >(P,RPAKCut);
      realHandler-> initBreakup(P);

      LRTwoBodyJastrow* LongRangeRPA = new LRTwoBodyJastrow(P,realHandler);
//       LongRangeRPA->resetByHandler();

      VarRegistry<RealType> tempO;
      LongRangeRPA->put(0,tempO);
      dPsi.push_back(LongRangeRPA);
//       double dPsiKVal = (*dPsi.begin())->evaluateLog(P,P.G,P.L);
//       cout<<"point-1 "<<(*dPsi.begin())->evaluateLog(P,P.G,P.L)<<endl;



      //short-range uses realHandler
      RealType Rcut = realHandler->get_rc()-0.1;
      GridType* myGrid = new GridType;
      int npts=static_cast<int>(Rcut/0.01)+1;
      myGrid->set(0,Rcut,npts);

      //create the numerical functor
      FuncType* nfunc = new FuncType;
      ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(realHandler);
      sra->setRmax(Rcut);
      nfunc->initialize(sra, myGrid);
#if !defined(HAVE_MPI)
      ofstream fout("drpa.short.dat");
      for (int i = 0; i < myGrid->size(); i++) {
        RealType r=(*myGrid)(i);
        fout << r << "   " << nfunc->evaluate(r) << "   "
            << realHandler->evaluate(r,1.0/r) << " " 
            << realHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j2= new TwoBodyJastrowOrbital<FuncType>(P);
      /*      j2->init(P);*/
//       for (int i=0; i<4; i++) j2->addFunc(nfunc);
      j2->addFunc("ZVZB",0,0,nfunc);
      dPsi.push_back(j2);
      
      double tValue=0.0;
      vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
      while(dit != dit_end) {
        tValue += (*dit)-> evaluateLog(P,P.G,P.L);
        std::cout<<"  "<<tValue<<endl;
        ++dit;
      }
//       double dPsiRVal = tValue-dPsiKVal;
      cout<<"point0 "<<tValue<<endl;
      
      ParticleSet* tempP = new ParticleSet(P);

      
      HandlerType* tempHandler= new LRRPAHandlerTemp<RPABreakup<Return_t>,LPQHIBasis >(*tempP,RPAKCut);
      tempHandler-> initBreakup(*tempP);
      LRTwoBodyJastrow* tempLongRangeRPA = new LRTwoBodyJastrow(*tempP,tempHandler);
      tempLongRangeRPA->put(0,tempO);
      
      vector<OrbitalBase*> tempPsi;
      tempPsi.push_back(tempLongRangeRPA);
      
//       double tempPsiKVal = (*tempPsi.begin())->evaluateLog(*tempP,(*tempP).G,(*tempP).L);

      
      //short-range uses realHandler
      RealType Rcut2 = tempHandler->get_rc()-0.1;
      GridType* myGrid2 = new GridType;
      int npts2=static_cast<int>(Rcut2/0.01)+1;
      myGrid2->set(0,Rcut2,npts2);

      //create the numerical functor
      FuncType* nfunc2 = new FuncType;
      ShortRangePartAdapter<RealType>* sra2 = new ShortRangePartAdapter<RealType>(tempHandler);
      sra2->setRmax(Rcut2);
      nfunc2->initialize(sra2, myGrid2);
#if !defined(HAVE_MPI)
      ofstream fout2("drpa2.short.dat");
      for (int i = 0; i < myGrid2->size(); i++) {
        RealType r=(*myGrid2)(i);
        fout2 << r << "   " << nfunc2->evaluate(r) << "   "
            << tempHandler->evaluate(r,1.0/r) << " " 
            << tempHandler->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j3= new TwoBodyJastrowOrbital<FuncType>(*tempP);
      /*      j3->init(*tempP);*/
//       for (int i=0; i<4; i++) j3->addFunc(nfunc2);
      j3->addFunc("ZVZB",0,0,nfunc2);
      tempPsi.push_back(j3);
      
      cout<<"point0.5 "<<endl;
      double mValue=0.0;
      vector<OrbitalBase*>::iterator mit(tempPsi.begin()), mit_end(tempPsi.end());
      while(mit != mit_end) {
        mValue += (*mit)-> evaluateLog(*tempP,(*tempP).G,(*tempP).L);
        ++mit;
      }
      
//       double tempPsiRVal = mValue-tempPsiKVal;
//       double lowV = (*tempPsi.begin())->evaluateLog(*tempP,(*tempP).G,(*tempP).L);
      cout<<"point1 "<< mValue <<endl;
      
      ParticleSet* tempP2 = new ParticleSet(P);
      /*      (*tempP2).initParticleSet();*/
      P.convert2Unit(P.R,(*tempP).R);
      double scale = 1.001;
      //RPAKCut=1.0/(std::sqrt(Rs)*scale); 
      (*tempP2).Lattice*=scale;
      (*tempP2).convert2Cart((*tempP).R,(*tempP2).R);
      (*tempP2).Lattice.SetLRCutoffs();
      (*tempP2).SK->UpdateNewCell((*tempP2).Lattice.LR_kc);
      (*tempP2).update();
      //(*tempP2).Lattice.LR_rc=P.Lattice.LR_rc;
      
      (*tempP2).Lattice.print(std::cout);
      std::cout<<(*tempP2).R<<endl;
      std::cout<<(*tempP2).Lattice.LR_kc<<endl;
      std::cout<<(*tempP).Lattice.Volume<<" <-p1,p2-> "<<(*tempP2).Lattice.Volume<<endl;
      std::cout<<" R_s for 1: "<<std::pow(3.0/4.0/M_PI*(*tempP).Lattice.Volume/static_cast<RealType>((*tempP).getTotalNum()),1.0/3.0)<<endl;
      std::cout<<" R_s for 2: "<<std::pow(3.0/4.0/M_PI*(*tempP2).Lattice.Volume/static_cast<RealType>((*tempP2).getTotalNum()),1.0/3.0)<<endl;
      std::cout<<Dot((*tempP2).R,(*tempP2).R)<<"  "<<Dot(P.R,P.R)<<endl;
          
      HandlerType* tempHandler2= new LRRPAHandlerTemp<RPABreakup<Return_t>,LPQHIBasis >(*tempP2,RPAKCut);
      tempHandler2-> initBreakup(*tempP2);
      LRTwoBodyJastrow* tempLongRangeRPA2 = new LRTwoBodyJastrow(*tempP2,tempHandler2);
      tempLongRangeRPA2->put(0,tempO);
      
      vector<OrbitalBase*> tempPsi2;
      tempPsi2.push_back(tempLongRangeRPA2);
      
      //RealType Rcut3 = tempHandler2->get_rc()-0.1;
      RealType Rcut3 = Rcut2;
      GridType* myGrid3 = new GridType;
      int npts3=static_cast<int>(Rcut3/0.01)+1;
      myGrid3->set(0,Rcut3,npts3);

      //create the numerical functor
      FuncType* nfunc3 = new FuncType;
      ShortRangePartAdapter<RealType>* sra3 = new ShortRangePartAdapter<RealType>(tempHandler2);
      sra3->setRmax(Rcut3);
      nfunc3->initialize(sra3, myGrid3);
#if !defined(HAVE_MPI)
      ofstream fout3("drpa3.short.dat");
      for (int i = 0; i < myGrid3->size(); i++) {
        RealType r=(*myGrid3)(i);
        fout3 << r << "   " << nfunc3->evaluate(r) << "   "
            << tempHandler2->evaluate(r,1.0/r) << " " 
            << tempHandler2->evaluateLR(r) << endl;
      }
#endif
      TwoBodyJastrowOrbital<FuncType>* j4= new TwoBodyJastrowOrbital<FuncType>(*tempP2);
      /*      j4->init(*tempP2);*/
//       for (int i=0; i<4; i++) j4->addFunc(nfunc3);
      j4->addFunc("ZVZB",0,0,nfunc3);
      tempPsi2.push_back(j4);
      
      cout<<"point1.5 "<<endl;
      double lValue=0.0;
      vector<OrbitalBase*>::iterator lit(tempPsi2.begin()), lit_end(tempPsi2.end());
      while(lit != lit_end) {
        lValue += (*lit)-> evaluateLog(*tempP2,(*tempP2).G,(*tempP2).L);
        ++lit;
      }
//       double lowV = (*tempPsi.begin())->evaluateLog(*tempP,(*tempP).G,(*tempP).L);
      cout<<"point1 "<< lValue <<endl;
      
      
      
/*      
      double highV = (*tempPsi2.begin())->evaluateLog(*tempP2,(*tempP2).G,(*tempP2).L);
      cout<<"point2 "<<highV<<endl;*/


      std::cout<<(lValue-mValue)/((scale-1.0)*Rs)<<" is numeric, "<<tValue<<" is analytic. "<<endl;
    }

    return true;
  }


  QMCHamiltonianBase* RPAPressureCorrection::makeClone(ParticleSet& P, TrialWaveFunction& psi)
  {
    return new RPAPressureCorrection(P);
  }
    
    
    
}


/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
* $Id: BareKineticEnergy.h 1581 2007-01-04 16:02:14Z jnkim $ 
***************************************************************************/

