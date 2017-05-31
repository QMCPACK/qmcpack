//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/RQMCMultiple.h"
#include "QMCDrivers/MultiChain.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Estimators/CSPolymerEstimator.h"
#include "Estimators/ComboPolymerEstimator.h"
#include "Estimators/MJPolymerEstimator.h"
#include "Estimators/HFPolymerEstimator.h"
#include "Estimators/HFDHE2PolymerEstimator.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
RQMCMultiple::RQMCMultiple(MCWalkerConfiguration& w,
                           TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),
  ReptileLength(21),
  NumTurns(0), Reptile(0), NewBead(0),
  multiEstimator(0), MSS(1.0)
{
  RootName = "rmc";
  QMCType ="RQMCMultiple";
  m_param.add(ReptileLength,"chains","int");
  QMCDriverMode.set(QMC_MULTIPLE,1);
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
  add_H_and_Psi(&h,&psi);
  nObs = h.sizeOfObservables();
  EvalInterval=1;
  //     SpeciesSet tspecies(W.getSpeciesSet());
  //     RealType mass = tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]));
  //     if (mass < 1e-12) {
//       mass=1.0;
//       tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]))=1.0;
//     }
//     oneovermass = 1.0/mass;
//     RealType oneoversqrtmass = std::sqrt(oneovermass);
//     Tau=brancher->getTau();
//     m_oneover2tau = 0.5*mass/Tau;
//     m_sqrttau = std::sqrt(Tau*oneovermass);
}

RQMCMultiple::~RQMCMultiple()
{
  if(Reptile)
    delete Reptile;
  delete NewBead;
}

/** main function to perform RQMC
*/
bool RQMCMultiple::run()
{
  if(MyCounter==0)
    initReptile();
  Reptile->open(RootName);
  multiEstimator->setTau(Tau);
  //multiEstimator->initialize(Reptile, Directionless, Tau, nSteps);
  Estimators->start(nBlocks,true);
  IndexType block = 0;
  do
  {
    //Loop over Blocks
    IndexType step = 0;
    NumTurns = 0;
    Estimators->startBlock(nSteps);
    do
    {
      //Loop over steps
      moveReptile();
      Reptile->Age +=1;
      step++;
      CurrentStep++;
      if (CurrentStep%EvalInterval ==0)
        Estimators->accumulate(W);
    }
    while(step<nSteps);
//     Estimators->accumulate(W);
    multiEstimator->evaluateDiff();
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAccept = 0;
    nReject = 0;
    block++;
    //       Reptile->record();
    recordBlock(block);
  }
  while(block<nBlocks);
  Estimators->stop();
  Reptile->record();
  Reptile->close();
  return finalize(nBlocks);
}

/** initialize Reptile
*
* The actions are
* - resize any initernal array
* - create Reptile if never initialized
*   -- if a previous configuration file is not "invalid", initialize the beads
* - set reference properties for the first run
*/
void RQMCMultiple::initReptile()
{
  //     m_oneover2tau=0.5/Tau;
  //     m_sqrttau=std::sqrt(Tau);
  //     Tauoverm = Tau/MSS;
  //     sqrtTauoverm = std::sqrt(Tauoverm);
  //Resize working arrays
  resizeArrays(Psi1.size());
  //Initial direction of growth. To be read if calculation is restarted.
  int InitialGrowthDirection(0);
  //Build NewBead. This takes care of a bunch of resizing operation and properties of the starting bead
  NewBead=new Bead(**W.begin() , scaleBeadDrift);
  //     NewBead->ScaleDrift = ((scaleBeadDrift=="true")||(scaleBeadDrift=="yes"))
  W.R = NewBead->R;
  W.update();
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    NewBead->Properties(ipsi,LOGPSI) = Psi1[ipsi]->evaluateLog(W);
    NewBead->Properties(ipsi,SIGN) = Psi1[ipsi]->getPhase();
    RealType eloc= H1[ipsi]->evaluate(W);
    NewBead->Properties(ipsi,LOCALENERGY)= eloc;
    H1[ipsi]->saveProperty(NewBead->getPropertyBase(ipsi));
    //*(NewBead->Gradients[ipsi])=W.G;
    Copy(W.G,*(NewBead->Gradients[ipsi]));
    NewBead->BeadSignWgt[ipsi]=1;
    NewBead->getScaledDriftSingle(branchEngine->LogNorm,Tauoverm,ipsi);
    NewBead->Action(ipsi,MinusDirection)= 0.5*m_oneover2tau*Dot(*NewBead->DriftVectors[ipsi],*NewBead->DriftVectors[ipsi]);
    NewBead->Action(ipsi,PlusDirection)=NewBead->Action(ipsi,MinusDirection);
    NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    NewBead->stepmade=0;
    NewBead->deltaRSquared[0]=0.0;
    NewBead->deltaRSquared[1]=0.0;
    NewBead->deltaRSquared[2]=0.0;
  }
  NewBead->getDrift(branchEngine->LogNorm);
  //     NewBead->TransProb[MinusDirection]=(0.5*Tau)*Dot(NewBead->Drift,NewBead->Drift) ;
  NewBead->TransProb[MinusDirection]=m_oneover2tau*Dot(NewBead->Drift,NewBead->Drift) ;
  NewBead->TransProb[PlusDirection]=NewBead->TransProb[MinusDirection];
  //Reptile is made up by replicating the first walker. To be read if restarted.
  //if(Reptile == 0) Reptile=new MultiChain(*W.begin(),ReptileLength,InitialGrowthDirection,nPsi);
  bool restartmode = false;
  bool reuseReptile= true;
  Reptile=W.getPolymer();
  if(Reptile==0)
  {
    reuseReptile=false;
    Reptile=new MultiChain(NewBead,ReptileLength,InitialGrowthDirection,nPsi);
    if(h5FileRoot.size() && (h5FileRoot != "invalid"))
    {
      app_log() << "Reading the previous multi-chain configurations" << std::endl;
      restartmode = Reptile->read(h5FileRoot);
    }
    W.setPolymer(Reptile);
  }
  multiEstimator->setPolymer(Reptile);
  if(!restartmode)
  {
    if(reuseReptile)
      checkReptileProperties();
    else
      setReptileProperties();
  }
} // END OF InitReptile

void RQMCMultiple::resizeArrays(int n)
{
  nPsi = n;
  nptcl=W.G.size();
  gRand.resize(nptcl);
  NewGlobalAction.resize(n);
  NewGlobalSignWgt.resize(n);
  DeltaG.resize(n);
  WeightSign.resize(n);
  ////Register properties for each walker
  //for(int ipsi=0; ipsi<nPsi; ipsi++) H1[ipsi]->add2WalkerProperty(W);
  //resize Walker::Properties to hold everything
  W.resetWalkerProperty(nPsi);
  H1[0]->setPrimary(true);
  for(int ipsi=1; ipsi<nPsi; ipsi++)
  {
    H1[ipsi]->setPrimary(false);
  }
}

//Set initial properties assuming all beads are occupying the same position
void RQMCMultiple::setReptileProperties()
{
  MultiChain::iterator bead(Reptile->begin());
  RealType spring_norm( -1.5e0 * std::log(4*std::acos(0.e0)) * (*bead)->Drift.size() * Reptile->Last );
  //Assign Reference Sign as the majority Sign
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Reptile->setRefSign(ipsi,(*bead)->Properties(ipsi,SIGN));
  //Compute the Global Action
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType LinkAction( 2*((*bead)->Action(ipsi,MinusDirection)+ (*bead)->Action(ipsi,Directionless) ));
    Reptile->GlobalAction[ipsi]= 2*(*bead)->Properties(ipsi,LOGPSI) + spring_norm
                                 - Reptile->Last*LinkAction - branchEngine->LogNorm[ipsi] ;
  }
  //Compute Global Sign weight (need to be initialized somewhere)
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Reptile->GlobalSignWgt[ipsi] = ReptileLength;
  //Compute reference action
  RealType RefAction(-1.0e20);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RefAction= std::max(RefAction,Reptile->GlobalAction[ipsi]);
  }
  //Compute Total Weight
  Reptile->GlobalWgt=0.0e0;
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType DeltaAction(Reptile->GlobalAction[ipsi]-RefAction);
    if(DeltaAction > -30)
      Reptile->GlobalWgt += std::exp(DeltaAction);
  }
  Reptile->GlobalWgt=std::log(Reptile->GlobalWgt)+RefAction;
  //Compute Umbrella Weight
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType DeltaAction(Reptile->GlobalAction[ipsi]-Reptile->GlobalWgt);
    if(DeltaAction > -30)
      Reptile->UmbrellaWeight[ipsi] = std::exp(DeltaAction);
    else
      Reptile->UmbrellaWeight[ipsi] = 0.0e0;
  }
}

void RQMCMultiple::checkReptileProperties()
{
  //Temporary vector
  std::vector<int> SumSign;
  SumSign.resize(nPsi);
  ///Assign a bunch of useful pointers
  MultiChain::iterator first_bead(Reptile->begin()), bead_end(Reptile->end());
  MultiChain::iterator bead(first_bead),last_bead(bead_end-1);
  ///Loop over beads to initialize action and WF
  // just for debugging
  int beadCount = 0;
  while(bead != bead_end)
  {
    ///Pointer to the current walker
    Bead& curW(**bead);
    //Do not re-evaluate the Properties
    ////Copy to W (ParticleSet) to compute distances, Psi and H
    W.R=curW.R;
    W.update();
    ///loop over WF to compute contribution to the action and WF
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      //Compute Energy and Psi and save in curW
      curW.Properties(ipsi,LOGPSI) = Psi1[ipsi]->evaluateLog(W);
      IndexType BeadSign = Reptile->getSign(curW.Properties(ipsi,SIGN) = Psi1[ipsi]->getPhase());
      RealType eloc= H1[ipsi]->evaluate(W);
      curW.Properties(ipsi,LOCALENERGY)= eloc;
      H1[ipsi]->saveProperty(curW.getPropertyBase(ipsi));
      //*curW.Gradients[ipsi]=W.G;
      Copy(W.G,*curW.Gradients[ipsi]);
      ///Initialize Kinetic Action
      RealType KinActMinus=0.0;
      RealType KinActPlus=0.0;
      curW.getScaledDrift(branchEngine->LogNorm,Tauoverm);
      // Compute contribution to the Action in the MinusDirection
      if(bead!=first_bead)
        //forward action
      {
        Bead& prevW(**(bead-1));
        deltaR=prevW.R-curW.R - (*curW.DriftVectors[ipsi]);
        KinActMinus=Dot(deltaR,deltaR);
      }
      // Compute contribution to the Action in the PlusDirection
      if(bead!=last_bead)
        //backward action
      {
        Bead& nextW(**(bead+1));
        deltaR=nextW.R-curW.R - (*curW.DriftVectors[ipsi]);
        KinActPlus=Dot(deltaR,deltaR);
      }
      //Compute the Total "Sign" of the Reptile
      SumSign[ipsi] += BeadSign;
      // Save them in curW
      curW.Action(ipsi,MinusDirection)=0.5*m_oneover2tau*KinActMinus;
      curW.Action(ipsi,PlusDirection)=0.5*m_oneover2tau*KinActPlus;
      curW.Action(ipsi,Directionless)=0.5*Tau*eloc;
    }
    ++bead;
    beadCount++;
  }// End Loop over beads
  //Assign Reference Sign as the majority Sign
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    if(SumSign[ipsi]>0)
      Reptile->RefSign[ipsi]=1;
    else
      Reptile->RefSign[ipsi]=-1;
  }
  //Compute Sign-weight for each bead
  bead=first_bead;
  while(bead != bead_end)
  {
    Bead& curW(**bead);
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      int BeadSign = Reptile->getSign(curW.Properties(ipsi,SIGN));
      curW.BeadSignWgt[ipsi]=std::abs((BeadSign+Reptile->RefSign[ipsi])/2);
    }
    ++bead;
  }
  // Compute initial drift for each bead
  bead=first_bead;
  while(bead != bead_end)
  {
    Bead& curW(**bead);
    //       curW.getScaledDrift(branchEngine->LogNorm,Tau);
    curW.getScaledDrift(branchEngine->LogNorm,Tauoverm);
    ++bead;
  }
  //Compute initial transition Probability within the chain
  bead=first_bead;
  while(bead != bead_end)
  {
    Bead& curW(**bead);
    RealType TrProbMinus=0.0;
    RealType TrProbPlus=0.0;
    // Compute contribution to the Transition Prob in the MinusDirection
    if(bead!=first_bead)
      //forward action
    {
      Bead& prevW(**(bead-1));
      deltaR=prevW.R-curW.R - curW.Drift;
      TrProbMinus=Dot(deltaR,deltaR);
    }
    // Compute contribution to the Transition Prob in the PlusDirection
    if(bead!=last_bead)
      //backward action
    {
      Bead& nextW(**(bead+1));
      deltaR=nextW.R-curW.R - curW.Drift;
      TrProbPlus=Dot(deltaR,deltaR);
    }
    curW.TransProb[MinusDirection]=TrProbMinus*m_oneover2tau ;
    curW.TransProb[PlusDirection]=TrProbPlus*m_oneover2tau;
    ++bead;
  }
  //Compute the Global Action
  bead=first_bead;
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    Reptile->GlobalAction[ipsi]=(*bead)->Properties(ipsi,LOGPSI);
    //app_log() << " WF : " << ipsi << " " << Reptile->GlobalAction[ipsi] << std::endl;
  }
  while(bead != last_bead)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      Reptile->GlobalAction[ipsi]-=
        ( (*bead)->Action(ipsi,PlusDirection) + (*(bead+1))->Action(ipsi,MinusDirection)+
          (*bead)->Action(ipsi,Directionless) + (*(bead+1))->Action(ipsi,Directionless)   );
    }
    bead++;
  }
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    Reptile->GlobalAction[ipsi]+=(*bead)->Properties(ipsi,LOGPSI);
    //app_log() << " WF : " << ipsi << " " << Reptile->GlobalAction[ipsi] << std::endl;
  }
  //Compute Global Sign weight (need to be initialized somewhere)
  bead=first_bead;
  while(bead != bead_end)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      Reptile->GlobalSignWgt[ipsi] += (*bead)->BeadSignWgt[ipsi];
    }
    ++bead;
  }
  //Compute reference action
  RealType RefAction(-10.0e20);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    WeightSign[ipsi]=std::max(0,Reptile->GlobalSignWgt[ipsi]-Reptile->Last);
    if(WeightSign[ipsi])
      RefAction= std::max(RefAction,Reptile->GlobalAction[ipsi]);
  }
  //Compute Total Weight
  Reptile->GlobalWgt=0.0e0;
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType DeltaAction(Reptile->GlobalAction[ipsi]-RefAction);
    if((WeightSign[ipsi]>0) && (DeltaAction > -30))
      Reptile->GlobalWgt += std::exp(DeltaAction);
  }
  Reptile->GlobalWgt=std::log(Reptile->GlobalWgt)+RefAction;
  //Compute Umbrella Weight
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType DeltaAction(Reptile->GlobalAction[ipsi]-Reptile->GlobalWgt);
    if((WeightSign[ipsi]>0) && (DeltaAction > -30))
      Reptile->UmbrellaWeight[ipsi] = std::exp(DeltaAction);
    else
      Reptile->UmbrellaWeight[ipsi] = 0.0e0;
    app_log() << "GA " << ipsi <<  " : " << Reptile->GlobalAction[ipsi] << std::endl;
    app_log() << "UW " << ipsi <<  " : " << Reptile->UmbrellaWeight[ipsi] << std::endl;
  }
  app_log() << "GW " <<  " : " << Reptile->GlobalWgt << std::endl;
}



void RQMCMultiple::recordBlock(int block)
{
  Reptile->record();
  //Write stuff
  //TEST CACHE
  //Estimators->report(CurrentStep);
  //TEST CACHE
  //app_error() << " BROKEN RQMCMultiple::recordBlock(int block) HDFWalkerOutput as 2007-04-16 " << std::endl;
  //HDFWalkerOutput WO(RootName,false,0);
  //WO.get(W);
  //WO.write(*branchEngine);
  //Reptile->write(WO.getFileID());
}

//   bool RQMCMultiple::put(xmlNodePtr q)
//   {
//     //Using enumeration
//     //MinusDirection=0;
//     //PlusDirection=1;
//     //Directionless=2;
//     nPsi=H1.size();
//     if(branchEngine->LogNorm.size()!=nPsi)
//     {
//       branchEngine->LogNorm.resize(nPsi);
//       for(int i=0; i<nPsi; i++) branchEngine->LogNorm[i]=0.e0;
//     }
//
//     Estimators = branchEngine->getEstimatorManager();
//     if(Estimators == 0)
//     {
//       Estimators = new EstimatorManager(myComm);
//       multiEstimator = new CSPolymerEstimator(H,nPsi);
//       //multiEstimator = new RQMCEstimator(H,nPsi);
//       Estimators->add(multiEstimator,Estimators->MainEstimatorName);
//       branchEngine->setEstimatorManager(Estimators);
//     }
//     return true;
//   }
bool RQMCMultiple::put(xmlNodePtr q)
{
  std::string observ("NONE");
  maxTouch=100;
  OhmmsAttributeSet attrib;
  ParameterSet nattrib;
  scaleBeadDrift = "no";
  attrib.add(EvalInterval,"EvaluateInterval");
  nattrib.add(MSS,"mass","double" );
  nattrib.add(maxTouch,"maxTouch","int" );
  attrib.add(scaleBeadDrift,"scaleDrift");
  attrib.put(q);
  nattrib.put(q);
  if ((scaleBeadDrift=="yes")||(scaleBeadDrift=="true"))
    app_log()<<"  Using Scaled drift."<< std::endl;
  else
    app_log()<<"  Drift set by gradient."<< std::endl;
  app_log()<<"  Evaluating observals every "<<EvalInterval<<" steps."<< std::endl;
  nPsi=H1.size();
  if(branchEngine->LogNorm.size()!=nPsi)
  {
    branchEngine->LogNorm.resize(nPsi);
    for(int i=0; i<nPsi; i++)
      branchEngine->LogNorm[i]=0.e0;
  }
  Estimators = branchEngine->getEstimatorManager();
  if(multiEstimator == 0)
  {
    if( Estimators)
      delete Estimators;
    Estimators = new EstimatorManager(myComm);
    xmlNodePtr q2(q);
    q2=q2->children;
    bool FoundEstimatorBlock=false;
//     while(q2!=NULL){
//       std::string cname((const char*)q2->name);
//       std::string observ="NONE";
//       int fq(1);
//       OhmmsAttributeSet Ettrib;
//       Ettrib.add(observ,"observables" );
//       Ettrib.add(fq,"freq" );
//       Ettrib.put(q2);
//       if (cname=="Estimators"){
// // 	if (observ=="NONE"){
// // 	  app_log()<<"Using normal Observables"<< std::endl;
// // 	  multiEstimator = new CSPolymerEstimator(H,nPsi);
// // 	  FoundEstimatorBlock=true;
// // 	} else
// //       if (observ=="ZVZB"){
// //           app_log()<<"Using ZVZB observables"<< std::endl;
// //           //         multiEstimator = new MJPolymerEstimator(H,nPsi);
// //           MJPolymerEstimator* MJp = new MJPolymerEstimator(H,nPsi);
// //           MJp->setpNorm(1.0/( W.Lattice.DIM *  W.Lattice.Volume));
// //           MJp->setrLen(ReptileLength);
// //
// //           multiEstimator = MJp;
// //           FoundEstimatorBlock=true;
// //       }
// //       if (observ=="HF"){
// //           app_log()<<"Using HF observables"<< std::endl;
// //           HFPolymerEstimator* MJp = new HFPolymerEstimator(H,nPsi);
// //           MJp->add_HF_Observables(ReptileLength,fq);
// //           MJp->setTau(Tau);
// //           multiEstimator = MJp;
// //           FoundEstimatorBlock=true;
// //       }
// //       else if (observ=="HFDHE2"){
// //
// // 	app_log()<<"TRUNC is BROKEN Using HFDHE2 observables"<< std::endl;
// // 	HFDHE2PolymerEstimator* HFp = new HFDHE2PolymerEstimator(H,nPsi);
// // 	HFp->setrLen(ReptileLength);
// // 	HFp->setpNorm(1.0/( W.Lattice.DIM *  W.Lattice.Volume));
// // 	multiEstimator = HFp;
// // 	FoundEstimatorBlock=true;
// //       } else if (observ=="Combo"){
// // 	app_log()<<"Using Combo observables"<< std::endl;
// // 	ComboPolymerEstimator* HFp = new ComboPolymerEstimator(H,nPsi);
// // 	HFp->put(q2,W,ReptileLength);
// // 	multiEstimator = HFp;
// // 	FoundEstimatorBlock=true;
// //       }
//     }
//     q2=q2->next;
//   }
    if (!FoundEstimatorBlock)
    {
      app_log()<<"Using normal Observables"<< std::endl;
      multiEstimator = new CSPolymerEstimator(H,nPsi);
    }
    Estimators->add(multiEstimator,Estimators->MainEstimatorName);
    branchEngine->setEstimatorManager(Estimators);
  }
  SpeciesSet tspecies(W.getSpeciesSet());
  RealType mass = tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]));
  if (mass < 1e-12)
  {
    mass=1.0;
    tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]))=1.0;
  }
  RealType oneovermass=1.0/mass;
  RealType oneoversqrtmass = std::sqrt(oneovermass);
  Tauoverm = (Tau*oneovermass);
  m_oneover2tau = 0.5/Tauoverm;
  m_sqrttau = std::sqrt(Tauoverm);
  return true;
}

void RQMCMultiple::moveReptile()
{
  //     m_oneover2tau=0.5/Tau;
  //     m_sqrttau=std::sqrt(Tau);
  //     Tauoverm = Tau/MSS;
  //     sqrtTauoverm = std::sqrt(Tauoverm);
  int ihead,inext,itail;
  //Depending on the growth direction initialize growth variables
  if(Reptile->GrowthDirection==MinusDirection)
  {
    forward=MinusDirection;
    backward=PlusDirection;
    ihead = 0;
    itail = Reptile->Last;
    inext = itail-1;
  }
  else
  {
    forward=PlusDirection;
    backward=MinusDirection;
    ihead = Reptile->Last;
    itail = 0;
    inext = 1;
  }
  //Point head to the growing end of the reptile
  Bead *head,*tail,*next;
  head = (*Reptile)[ihead];
  tail=(*Reptile)[itail];
  //Treat the case with one bead differently
  if(ReptileLength==1)
    next=NewBead;
  else
    next=(*Reptile)[inext];
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandom(gRand);
  tail->timesTouched +=1;
  bool FORCE=0;
  if (tail->timesTouched > maxTouch)
  {
    FORCE=1;
    app_log()<<"Forcing Accept!!"<< std::endl;
  };
  //     //new position,\f$R_{yp} = R_{xp} + \sqrt(2}*\Delta + tau*D_{xp}\f$
  //     W.R = head->R + m_sqrttau*gRand + Tau*head->Drift;
  //     //Save Transition Probability
  //     head->TransProb[forward]=0.5*Dot(gRand,gRand);
  //     //Save position in NewBead
  //     NewBead->R=W.R;
  //     //update the distance table associated with W
  //     //DistanceTable::update(W);
  //     W.update();
  //
  //     //Compute the "diffusion-drift"\f$=(R_{yp}-R{xp})/tau\f$
  //     deltaR= NewBead->R - head->R;
  //
  //     //Compute HEAD action in forward direction
  //     for(int ipsi=0; ipsi<nPsi; ipsi++) {
//       gRand = deltaR-Tau*(*head->Gradients[ipsi]);
//       head->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(gRand,gRand);
//     }
//Just making sure this is set
  NewBead->stepmade=Reptile->Age;
  NewBead->timesTouched=0;
//     //new position,\f$R_{yp} = R_{xp} + \sqrt(2}*\Delta + tau*D_{xp}\f$
//     W.R = head->R + m_sqrttau*gRand + Tau*head->Drift;
//     //Save Transition Probability
//     head->TransProb[forward]=0.5*Dot(gRand,gRand);
//     head->getScaledDrift(branchEngine->LogNorm,Tauoverm);
//     deltaR =  m_sqrttau*gRand + head->Drift;
  deltaR =  m_sqrttau*gRand + head->Drift;
  W.R = head->R + deltaR;
//Save Transition Probability
  head->TransProb[forward]=0.5*Dot(gRand,gRand);
  head->deltaRSquared[forward]= Dot(deltaR,deltaR);
  NewBead->deltaRSquared[backward]= Dot(deltaR,deltaR);
  NewBead->deltaRSquared[forward]=0.0;
//     ParticleSet::ParticlePos_t DR1 = head->Gradients[0];
  ParticleSet::ParticlePos_t DR1 = head->Drift;
  DR1 *= 1.0/Tauoverm;
  head->deltaRSquared[Directionless]= Dot(DR1,DR1);
//Save position in NewBead
  NewBead->R=W.R;
//update the distance table associated with W
//DistanceTable::update(W);
  W.update();
//Compute the "diffusion-drift"\f$=(R_{yp}-R{xp})/tau\f$
  deltaR= NewBead->R - head->R;
//Compute HEAD action in forward direction
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    //       gRand = deltaR-Tau*(*head->Gradients[ipsi]);
    gRand = deltaR - *head->DriftVectors[ipsi];
    //       head->Action(ipsi,forward)= 0.5*m_oneover2tau*Dot(gRand,gRand);
    head->Action(ipsi,forward)= 0.5*m_oneover2tau*Dot(gRand,gRand);
  }
  bool FAIL=0;
//evaluate all relevant quantities in the new position
  int totbeadwgt(0);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType* restrict NewBeadProp=NewBead->getPropertyBase(ipsi);
    RealType* restrict HeadProp=head->getPropertyBase(ipsi);
    //evaluate Psi and H
    NewBeadProp[LOGPSI]=Psi1[ipsi]->evaluateLog(W);
    NewBeadProp[SIGN]=Psi1[ipsi]->getPhase();
    RealType eloc=NewBeadProp[LOCALENERGY]= H1[ipsi]->evaluate(W);
    H1[ipsi]->auxHevaluate(W);
    //Save properties
    H1[ipsi]->saveProperty(NewBeadProp);
    //*(NewBead->Gradients[ipsi])=W.G;
    Copy(W.G,*(NewBead->Gradients[ipsi]));
    //       NewBead->deltaRSquared[Directionless]= Dot(W.G,W.G);
    //       NewBead->getScaledDrift(branchEngine->LogNorm,Tau);
    //       NewBead->getScaledDrift(branchEngine->LogNorm,Tauoverm);
    NewBead->getScaledDriftSingle(branchEngine->LogNorm,Tauoverm,ipsi);
    ParticleSet::ParticlePos_t DR2 = NewBead->Drift;
    DR2 *= 1.0/Tauoverm;
    NewBead->deltaRSquared[Directionless]= Dot(DR2,DR2);
    //Compute the backward part of the Kinetic action
    //gRand=deltaR+Tau*W.G;
    PAOps<RealType,DIM>::axpy(1.0,*NewBead->DriftVectors[ipsi],deltaR,gRand);
    //       NewBead->Action(ipsi,backward)= 0.5*m_oneover2tau*Dot(gRand,gRand);
    NewBead->Action(ipsi,backward)= 0.5*m_oneover2tau*Dot(gRand,gRand);
    NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
    NewBead->BeadSignWgt[ipsi]=beadwgt;
    totbeadwgt+=beadwgt;
    //       RealType* restrict NewBeadProp=NewBead->getPropertyBase(ipsi);
    //
    //       //evaluate Psi and H
    //       NewBeadProp[LOGPSI]=Psi1[ipsi]->evaluateLog(W);
    //       NewBeadProp[SIGN]=Psi1[ipsi]->getPhase();
    //       RealType eloc=NewBeadProp[LOCALENERGY]= H1[ipsi]->evaluate(W);
    //
    //       //Save properties
    //       H1[ipsi]->saveProperty(NewBeadProp);
    //       //*(NewBead->Gradients[ipsi])=W.G;
    //       Copy(W.G,*(NewBead->Gradients[ipsi]));
    //
    //       //Compute the backward part of the Kinetic action
    //       //gRand=deltaR+Tau*W.G;
    //       PAOps<RealType,DIM>::axpy(Tau,W.G,deltaR,gRand);
    //       NewBead->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(gRand,gRand);
    //
    //       NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    //       int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
    //       NewBead->BeadSignWgt[ipsi]=beadwgt;
    //       totbeadwgt+=beadwgt;
//    if ((NewBeadProp[LOCALENERGY] <= Ecut) && (NewBeadProp[LOCALENERGY] < HeadProp[LOCALENERGY]) ) {
//      FAIL=1;
//    }
  }
  if (!FAIL)
  {
    //Compute Drift and TransProb here. This could be done (and was originally done)
    //after acceptance if # of beads is greater than one but needs to be done here
    //to make the one-bead case working ok.
    //       NewBead->getScaledDrift(branchEngine->LogNorm,Tau);
    NewBead->getDrift(branchEngine->LogNorm);
    gRand=deltaR+ NewBead->Drift;
    //       NewBead->TransProb[backward]=m_oneover2tau*Dot(gRand,gRand);
    NewBead->TransProb[backward]=m_oneover2tau*Dot(gRand,gRand);
    RealType AcceptProb(-1.0),NewGlobalWgt(0.0);
    RealType RefAction(-1.0e20);
    if(totbeadwgt!=0)
    {
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        DeltaG[ipsi]=
          - head->Action(ipsi,forward)       - NewBead->Action(ipsi,backward)
          + tail->Action(ipsi,forward)       + next->Action(ipsi,backward)
          - head->Action(ipsi,Directionless) - NewBead->Action(ipsi,Directionless)
          + tail->Action(ipsi,Directionless) + next->Action(ipsi,Directionless)
          - head->Properties(ipsi,LOGPSI)    + NewBead->Properties(ipsi,LOGPSI)
          - tail->Properties(ipsi,LOGPSI)    + next->Properties(ipsi,LOGPSI);
        NewGlobalAction[ipsi]=Reptile->GlobalAction[ipsi]+DeltaG[ipsi];
        //Compute the new sign
        NewGlobalSignWgt[ipsi]=Reptile->GlobalSignWgt[ipsi]+
                               NewBead->BeadSignWgt[ipsi]-tail->BeadSignWgt[ipsi];
        //Weight: 1 if all beads have the same sign. 0 otherwise.
        WeightSign[ipsi]=std::max(0,NewGlobalSignWgt[ipsi]-Reptile->Last);
        // Assign Reference Action
        if(WeightSign[ipsi]>0)
          RefAction= std::max(RefAction,NewGlobalAction[ipsi]);
      }
      //Compute Log of global Wgt
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        RealType DeltaAction(NewGlobalAction[ipsi]-RefAction);
        if((WeightSign[ipsi]>0) && (DeltaAction > -30.0))
          NewGlobalWgt+=std::exp(DeltaAction);
      }
      NewGlobalWgt=std::log(NewGlobalWgt)+RefAction;
      AcceptProb=std::exp(NewGlobalWgt - Reptile->GlobalWgt + head->TransProb[forward] - next->TransProb[backward]);
    }
    //FORCE acceptance if tail bead is the stuck one.
    if((Random() < AcceptProb)||(FORCE))
    {
      //Update Reptile information
      Reptile->GlobalWgt=NewGlobalWgt;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        Reptile->GlobalAction[ipsi]=NewGlobalAction[ipsi];
        Reptile->GlobalSignWgt[ipsi]=NewGlobalSignWgt[ipsi];
        RealType DeltaAction(NewGlobalAction[ipsi]-NewGlobalWgt);
        if((WeightSign[ipsi]>0) && (DeltaAction > -30.0))
          Reptile->UmbrellaWeight[ipsi]=std::exp(DeltaAction);
        else
          Reptile->UmbrellaWeight[ipsi]=0.0e0;
      }
      //Add NewBead to the Polymer.
      if(Reptile->GrowthDirection==MinusDirection)
      {
        Reptile->push_front(NewBead);
        NewBead=tail;
        Reptile->pop_back();
      }
      else
      {
        Reptile->push_back(NewBead);
        NewBead=tail;
        Reptile->pop_front();
      }
      ++nAccept;
    }
    else
    {
      ++nReject;
      ++NumTurns;
      Reptile->flip();
    }
  }
  else
  {
    ++nReject;
    ++NumTurns;
    Reptile->flip();
  }
}
}
