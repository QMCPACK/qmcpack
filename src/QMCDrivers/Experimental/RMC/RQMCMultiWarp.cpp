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
    
    


#include "QMCDrivers/RQMCMultiWarp.h"
#include "QMCDrivers/MultiChain.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCApp/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus
{
RQMCMultiWarp::RQMCMultiWarp(MCWalkerConfiguration& w,
                             TrialWaveFunction& psi, QMCHamiltonian& h,
                             ParticleSetPool& ptclPool):
  QMCDriver(w,psi,h), PtclPool(ptclPool),
  ReptileLength(21),
  NumTurns(0), Reptile(0), NewBead(0)
{
  RootName = "rmc-warp";
  QMCType ="rmc-warp";
  m_param.add(ReptileLength,"chains","int");
  refSetName="invalid";
  m_param.add(refSetName,"reference","str");
  QMCDriverMode.set(QMC_MULTIPLE,1);
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
  LOGJACOB=w.addProperty("logjacob");
  add_H_and_Psi(&h,&psi);
}

RQMCMultiWarp::~RQMCMultiWarp()
{
  delete NewBead;
}

/** initialize Reptile
 *
 * The actions are
 * - resize any initernal array
 * - create Reptile if never initialized
 *   -- if a previous configuration file is not "invalid", initialize the beads
 * - set reference properties for the first run
 */
void RQMCMultiWarp::initReptile()
{
  //Initial direction of growth. To be read if calculation is restarted.
  int InitialGrowthDirection(0);
  //Build NewBead. This takes care of a bunch of resizing operation and properties of the starting bead
  //Resize working arrays
  resizeArrays(Psi1.size());
  NewBead=new Bead(**W.begin());
  W.R = NewBead->R;
  W.update();
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Jacobian[ipsi]=1.0;
  for(int iptcl=0; iptcl< nptcl; iptcl++)
  {
    PtclWarp.warp_one(iptcl,0);
    //Save particle position
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      WW[ipsi]->R[iptcl]=W.R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
      Jacobian[ipsi]*=PtclWarp.get_Jacobian(iptcl,ipsi);
    }
  }
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    WW[ipsi]->update();
    NewBead->Properties(ipsi,LOGPSI) = Psi1[ipsi]->evaluateLog(*WW[ipsi]);
    NewBead->Properties(ipsi,SIGN) = Psi1[ipsi]->getPhase();
    NewBead->Properties(ipsi,LOGJACOB)=std::log(std::abs(Jacobian[ipsi]));
    RealType eloc= H1[ipsi]->evaluate(*WW[ipsi]);
    NewBead->Properties(ipsi,LOCALENERGY)= eloc;
    H1[ipsi]->saveProperty(NewBead->getPropertyBase(ipsi));
    //*(NewBead->Gradients[ipsi])=WW[ipsi]->G;
    Copy(WW[ipsi]->G,*(NewBead->Gradients[ipsi]));
    NewBead->BeadSignWgt[ipsi]=1;
    NewBead->Action(ipsi,MinusDirection)=0.25*Tau*Dot(*NewBead->Gradients[ipsi],*NewBead->Gradients[ipsi]);
    NewBead->Action(ipsi,PlusDirection)=NewBead->Action(ipsi,MinusDirection);
    NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
  }
  NewBead->getDrift(Jacobian,PtclWarp,branchEngine->LogNorm);
  NewBead->TransProb[MinusDirection]=(0.5*Tau)*Dot(NewBead->Drift,NewBead->Drift) ;
  NewBead->TransProb[PlusDirection]=NewBead->TransProb[MinusDirection];
  //Reptile is made up by replicating the NewBead. To be read if restarted.
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
  if(restartmode)
  {
    int ihead;
    if(Reptile->GrowthDirection==MinusDirection)
      ihead = 0;
    else
      ihead = Reptile->Last;
    Bead *head;
    head = (*Reptile)[ihead];
    W.R=head->R;
    W.update();
    for(int iptcl=0; iptcl< nptcl; iptcl++)
    {
      PtclWarp.warp_one(iptcl,0);
      for(int ipsi=0; ipsi<nPsi; ipsi++)
        WW[ipsi]->R[iptcl]=head->R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
    }
    //Save initial warped position for move
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Warped_R[ipsi]=WW[ipsi]->R;
  }
  else
  {
    if(reuseReptile)
      checkReptileProperties();
    else
      setReptileProperties();
  }
} // END OF InitReptile


void RQMCMultiWarp::resizeArrays(int n)
{
  nPsi = n;
  nptcl=W.G.size();
  gRand.resize(nptcl);
  NewGlobalAction.resize(n);
  NewGlobalSignWgt.resize(n);
  DeltaG.resize(n);
  WeightSign.resize(n);
  Jacobian.resize(n);
  for (int ipsi=0; ipsi<n; ipsi++)
  {
    Warped_R.push_back(ParticlePos_t(nptcl));
    Warped_deltaR.push_back(ParticlePos_t(nptcl));
  }
  //Register properties for each walker
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
void RQMCMultiWarp::setReptileProperties()
{
  MultiChain::iterator bead(Reptile->begin());
  //Assign Reference Sign as the majority Sign
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Reptile->setRefSign(ipsi,(*bead)->Properties(ipsi,SIGN));
  //Reptile->RefSign[ipsi]=(*bead)->Properties(ipsi,SIGN);
  RealType spring_norm( -1.5e0 * std::log(4*std::acos(0.e0)) * (*bead)->Drift.size() * Reptile->Last );
  //Compute the Global Action
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType LinkAction( 2*((*bead)->Action(ipsi,MinusDirection)+ (*bead)->Action(ipsi,Directionless) ));
    Reptile->GlobalAction[ipsi]=2*(*bead)->Properties(ipsi,LOGPSI) + spring_norm
                                + ReptileLength*(*bead)->Properties(ipsi,LOGJACOB) - Reptile->Last*LinkAction
                                - branchEngine->LogNorm[ipsi] ;
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

void RQMCMultiWarp::checkReptileProperties()
{
  //Temporary vector
  typedef MCWalkerConfiguration::ParticlePos_t ParticlePos_t;
  std::vector<ParticlePos_t> Warped_R_prev,Warped_R_next,Warped_R_curr;
  std::vector<int> SumSign(nPsi);
  std::vector<RealType>Jacobian_next(nPsi);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    Warped_R_prev.push_back(ParticlePos_t(nptcl));
    Warped_R_next.push_back(ParticlePos_t(nptcl));
    Warped_R_curr.push_back(ParticlePos_t(nptcl));
  }
  ///Assign a bunch of useful pointers
  MultiChain::iterator first_bead(Reptile->begin()), bead_end(Reptile->end());
  MultiChain::iterator bead(first_bead),last_bead(bead_end-1);
  ///Loop over beads to initialize action and WF
  W.R=(*first_bead)->R;
  W.update();
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Jacobian[ipsi]=1.e0;
  for(int iptcl=0; iptcl< nptcl; iptcl++)
  {
    PtclWarp.warp_one(iptcl,0);
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      (Warped_R_curr[ipsi])[iptcl]=(*first_bead)->R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
      Jacobian[ipsi]*=PtclWarp.get_Jacobian(iptcl,ipsi);
    }
  }
  //cout << "<<<<<<<<< setRefProperties():Step1 " << std::endl;
  //Save initial warped position for next move
  //If growth is in Plus direction do it at the end
  if(Reptile->GrowthDirection==MinusDirection)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Warped_R[ipsi]=Warped_R_curr[ipsi];
  }
  //cout << "<<<<<<<<< setRefProperties():Step2 " << std::endl;
  while(bead != bead_end)
  {
    ///Pointer to the current walker
    Bead& curW(**bead);
    //Do not re-evaluate the Properties
    ////Copy to W (ParticleSet) to compute distances, Psi and H
    if(bead!=first_bead)
    {
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        Warped_R_prev[ipsi]=Warped_R_curr[ipsi];
        Warped_R_curr[ipsi]=Warped_R_next[ipsi];
        Jacobian[ipsi]=Jacobian_next[ipsi];
      }
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Jacobian_next[ipsi]=1.e0;
    if(bead!=last_bead)
    {
      Bead& nextW(**(bead+1));
      W.R=nextW.R;
      W.update();
      for(int iptcl=0; iptcl< nptcl; iptcl++)
      {
        PtclWarp.warp_one(iptcl,0);
        //Save particle position
        for(int ipsi=0; ipsi<nPsi; ipsi++)
        {
          (Warped_R_next[ipsi])[iptcl]=nextW.R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
          Jacobian_next[ipsi]*=PtclWarp.get_Jacobian(iptcl,ipsi);
        }
      }
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      WW[ipsi]->R = Warped_R_curr[ipsi];
      WW[ipsi]->update();
      //Compute Energy and Psi and save in curW
      curW.Properties(ipsi,LOGPSI) = Psi1[ipsi]->evaluateLog(*WW[ipsi]);
      RealType BeadSign = Reptile->getSign(curW.Properties(ipsi,SIGN) = Psi1[ipsi]->getPhase());
      //RealType BeadSign = curW.Properties(ipsi,SIGN) = Psi1[ipsi]->getPhase();
      RealType eloc= H1[ipsi]->evaluate(*WW[ipsi]);
      curW.Properties(ipsi,LOCALENERGY)= eloc;
      H1[ipsi]->saveProperty(curW.getPropertyBase(ipsi));
      //*curW.Gradients[ipsi]=WW[ipsi]->G;
      Copy(WW[ipsi]->G,*curW.Gradients[ipsi]);
      curW.Properties(ipsi,LOGJACOB)=std::log(std::abs(Jacobian[ipsi]));
      ///Initialize Kinetic Action
      RealType KinActMinus=0.0;
      RealType KinActPlus=0.0;
      // Compute contribution to the Action in the MinusDirection
      if(bead!=first_bead)
        //forward action
      {
        deltaR = Warped_R_prev[ipsi] - Warped_R_curr[ipsi] - Tau*(*curW.Gradients[ipsi]);
        KinActMinus=Dot(deltaR,deltaR);
      }
      // Compute contribution to the Action in the PlusDirection
      if(bead!=last_bead)
        //backward action
      {
        deltaR = Warped_R_next[ipsi] - Warped_R_curr[ipsi] - Tau*(*curW.Gradients[ipsi]);
        KinActPlus=Dot(deltaR,deltaR);
      }
      //Compute the Total "Sign" of the Reptile
      SumSign[ipsi] += int(BeadSign);
      // Save them in curW
      curW.Action(ipsi,MinusDirection)=0.25/Tau*KinActMinus;
      curW.Action(ipsi,PlusDirection)=0.25/Tau*KinActPlus;
      curW.Action(ipsi,Directionless)=0.5*Tau*eloc;
    }
    ++bead;
  }// End Loop over beads
  //cout << "<<<<<<<<< setRefProperties():Step3 " << std::endl;
  //Store warped position for next move if direction is +
  //otherwise this is done at the beginning of the bead loop
  if(Reptile->GrowthDirection==PlusDirection)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Warped_R[ipsi]=Warped_R_curr[ipsi];
  }
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
      //int BeadSign = int(curW.Properties(ipsi,SIGN));
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
    W.R=curW.R;
    W.update();
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Jacobian[ipsi]=1.0e0;
    //Warp, compute Jacobian and Jacobian drift
    for(int iptcl=0; iptcl< nptcl; iptcl++)
    {
      PtclWarp.warp_one(iptcl,0);
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        WW[ipsi]->R[iptcl]=(*first_bead)->R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
        Jacobian[ipsi]*=PtclWarp.get_Jacobian(iptcl,ipsi);
      }
    }
    curW.getDrift(Jacobian,PtclWarp,branchEngine->LogNorm);
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
      deltaR=prevW.R-curW.R - Tau*curW.Drift;
      TrProbMinus=Dot(deltaR,deltaR);
    }
    // Compute contribution to the Transition Prob in the PlusDirection
    if(bead!=last_bead)
      //backward action
    {
      Bead& nextW(**(bead+1));
      deltaR=nextW.R-curW.R - Tau*curW.Drift;
      TrProbPlus=Dot(deltaR,deltaR);
    }
    curW.TransProb[MinusDirection]=TrProbMinus*0.5/Tau ;
    curW.TransProb[PlusDirection]=TrProbPlus*0.5/Tau;
    ++bead;
  }
  //Compute the Global Action
  bead=first_bead;
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Reptile->GlobalAction[ipsi]=(*bead)->Properties(ipsi,LOGPSI);
  while(bead != last_bead)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      RealType LinkAction( (*bead)->Action(ipsi,PlusDirection) +
                           (*(bead+1))->Action(ipsi,MinusDirection)+
                           (*bead)->Action(ipsi,Directionless) +
                           (*(bead+1))->Action(ipsi,Directionless) );
      //Reptile->GlobalAction[ipsi]-=( (*bead)->Action(ipsi,PlusDirection) +
      //                               (*(bead+1))->Action(ipsi,MinusDirection)+
      //                               (*bead)->Action(ipsi,Directionless) +
      //                               (*(bead+1))->Action(ipsi,Directionless) );
      Reptile->GlobalAction[ipsi]+=((*bead)->Properties(ipsi,LOGJACOB)-LinkAction);
      //cout << "Psi: " << ipsi << " LA " << LinkAction << "   LJ " <<
      //        (*bead)->Properties(ipsi,LOGJACOB) << std::endl;
    }
    ++bead;
  }
  RealType spring_norm( -1.5e0 * std::log(4*std::acos(0.e0)) * (*bead)->Drift.size() * Reptile->Last );
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Reptile->GlobalAction[ipsi]+= ( (*bead)->Properties(ipsi,LOGPSI) +
                                    (*bead)->Properties(ipsi,LOGJACOB) -
                                    branchEngine->LogNorm[ipsi] + spring_norm );
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
  RealType RefAction(-1.0e20);
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
  }
}

bool RQMCMultiWarp::run()
{
  //TEST CACHE
  //Estimators->reportHeader(AppendRun);
  Estimators->start(nBlocks);
  //TEST CACHE
  //Check if we need to update the norm of the wave functions
  if(MyCounter==0)
    initReptile();
  Reptile->open(RootName);
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  std::vector<RealType> AveEloc,AveWeight;
  AveEloc.resize(nPsi);
  AveWeight.resize(nPsi);
  std::string filename(RootName);
  filename=RootName+".Eloc.dat";
  std::ofstream *OutEnergy;
  OutEnergy=new std::ofstream(filename.c_str());
  //PolymerEstimator reptileReport(*Reptile,nPsi);
  //reptileReport.resetReportSettings(RootName);
  //accumulate configuration: probably need to reorder
  //HDFWalkerOutput WO(RootName);
  std::vector<RealType>Norm(nPsi);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Norm[ipsi]=0.0;
  RealType oneoversteps=1.0/static_cast<RealType>(nSteps);
  do
    //Loop over Blocks
  {
    IndexType step = 0;
    NumTurns = 0;
    Estimators->startBlock(nSteps);
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      AveEloc[ipsi]=0.0;
      AveWeight[ipsi]=0.0;
    }
    do
      //Loop over steps
    {
      moveReptile();
      step++;
      CurrentStep++;
      //Copy the front and back to W to take average report
      //W.copyWalkerRefs(Reptile->front(),Reptile->back());
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        double WeightedEloc=Reptile->UmbrellaWeight[ipsi]*
                            ( Reptile->front()->Action(ipsi,Directionless)
                              +Reptile->back()->Action(ipsi,Directionless) );
        AveEloc[ipsi]+=WeightedEloc;
        AveWeight[ipsi]+=Reptile->UmbrellaWeight[ipsi];
      }
      Estimators->accumulate(W);
      //use the first energy for the branch
      //branchEngine->accumulate(AveEloc[0],1.0);
      //reptileReport.accumulate();
    }
    while(step<nSteps);
    //multiEstimator->evaluateDiff();
    RealType acceptedR = static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject);
    Estimators->stopBlock(acceptedR);
    *OutEnergy << block << " " ;
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      AveEloc[ipsi]/=(AveWeight[ipsi]*Tau+numeric_limits<RealType>::epsilon());
      *OutEnergy << AveEloc[ipsi] << " ";
      *OutEnergy << AveWeight[ipsi]/nSteps << " ";
    }
    *OutEnergy << acceptedR << std::endl;
    OutEnergy->flush();
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    nAccept = 0;
    nReject = 0;
    block++;
    Reptile->record();
    //recordBlock(block);
  }
  while(block<nBlocks);
  delete OutEnergy;
  //Need MPI-IO
  app_log() << "ratio = "
            << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
            << std::endl;
  Reptile->close();
  return finalize(nBlocks);
}

void RQMCMultiWarp::recordBlock(int block)
{
  //Write stuff
  //TEST CACHE
  //Estimators->report(CurrentStep);
  //TEST CACHE
  //app_error() << " BROKEN RQMCMultiWarp::recordBlock(int block) HDFWalkerOutput as 2007-04-16 " << std::endl;
  //HDFWalkerOutput WO(RootName,false,0);
  //WO.get(W);
  //WO.write(*branchEngine);
  //Reptile->write(WO.getFileID());
  //Reptile->write(RootName);
}

bool RQMCMultiWarp::put(xmlNodePtr q)
{
  MinusDirection=0;
  PlusDirection=1;
  Directionless=2;
  if(WW.empty())
  {
    W.clearDistanceTables();
  }
  //qmcsystem
  std::vector<ParticleSet*> ionSets;
  DistanceTableData* dtableReference;
  xmlNodePtr cur=q->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "qmcsystem")
    {
      OhmmsAttributeSet aAttrib;
      std::string source_name("i");
      aAttrib.add(source_name,"source");
      aAttrib.put(cur);
      ParticleSet* ions=PtclPool.getParticleSet(source_name);
      if(ions)
        ionSets.push_back(ions);
      //string source_name((const char*)xmlGetProp(cur,(const xmlChar*)"source"));
      //ionSets.push_back(PtclPool.getParticleSet(source_name));
    }
    cur=cur->next;
  }
  ParticleSet* p(0);
  if(refSetName!="invalid")
  {
    p=PtclPool.getParticleSet(refSetName);
    if(p==0)
    {
      std::cout << "The specified reference cannot be found. Stop." << std::endl;
      APP_ABORT("RQMCMultiWarp::put");
    }
  }
  else
  {
    refSetName=ionSets[0]->getName().c_str();
    p=PtclPool.getParticleSet(refSetName);
  }
  dtableReference=DistanceTable::add(*p,W);
  /*vector<DistanceTableData*> dtableList;
  std::string target_name(W.getName());
  xmlNodePtr cur=q->children;
  while(cur != NULL) {
    std::string cname((const char*)(cur->name));
    if(cname == "qmcsystem") {
      std::string source_name((const char*)xmlGetProp(cur,(const xmlChar*)"source"));
      dtableList.push_back(DistanceTable::getTable(source_name.c_str(),target_name.c_str()));
    }
    cur=cur->next;
  }*/
  PtclWarp.initialize(ionSets,dtableReference);
  //PtclWarp.initialize(dtableList);
  nPsi=Psi1.size();
  if(WW.empty())
  {
    //WW.push_back(&W);
    char newname[128];
    //for(int ipsi=1; ipsi<nPsi; ipsi++){
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      sprintf(newname,"%s%d", W.getName().c_str(),ipsi);
      ParticleSet* pclone=PtclPool.getParticleSet(newname);
      if(pclone == 0)
      {
        app_log() << "  Cloning particle set in VMCMultipleWarp " << newname << std::endl;
        pclone=new ParticleSet(W);
        pclone->setName(newname);
        PtclPool.addParticleSet(pclone);
      }
      else
      {
        app_log() << "  Cloned particle exists " << newname << std::endl;
      }
      WW.push_back(pclone);
      Psi1[ipsi]->resetTargetParticleSet(*WW[ipsi]);
      H1[ipsi]->resetTargetParticleSet(*WW[ipsi]);
    }
  }
  /*if(WW.empty()){
    WW.push_back(&W);
    char newname[128];
    for(int ipsi=1; ipsi<nPsi; ipsi++){
  sprintf(newname,"%s%d", W.getName().c_str(),ipsi);
      ParticleSet* pclone=PtclPool.getParticleSet(newname);
      if(pclone == 0) {
        app_log() << "  Cloning particle set in VMCMultipleWarp " << newname << std::endl;
        pclone=new ParticleSet(W);
        pclone->setName(newname);
        PtclPool.addParticleSet(pclone);
      } else {
        app_log() << "  Cloned particle exists " << newname << std::endl;
      }
  //Correct copy constructor????????
  WW.push_back(pclone);
  WW[ipsi]=pclone;
  Psi1[ipsi]->resetTargetParticleSet(*WW[ipsi]);
      H1[ipsi]->resetTargetParticleSet(*WW[ipsi]);
    }
  }*/
  if(branchEngine->LogNorm.size()==0)
  {
    branchEngine->LogNorm.resize(nPsi);
    for(int i=0; i<nPsi; i++)
      branchEngine->LogNorm[i]=0.e0;
  }
  return true;
}

void RQMCMultiWarp::moveReptile()
{
  //Used several times
  m_oneover2tau=0.5/Tau;
  m_sqrttau=std::sqrt(Tau);
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
  //new position,\f$R_{yp} = R_{xp} + \sqrt(2}*\Delta + tau*D_{xp}\f$
  W.R = head->R + m_sqrttau*gRand + Tau*head->Drift;
  //Save Transition Probability
  head->TransProb[forward]=0.5*Dot(gRand,gRand);
  //Save position in NewBead
  NewBead->R=W.R;
  //update the distance table associated with W
  //DistanceTable::update(W);
  W.update();
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Jacobian[ipsi]=1.e0;
  for(int iptcl=0; iptcl< nptcl; iptcl++)
  {
    PtclWarp.warp_one(iptcl,0);
    //Save particle position
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      //cout << PtclWarp.get_displacement(iptcl,ipsi) << std::endl;
      WW[ipsi]->R[iptcl]=W.R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
      Jacobian[ipsi]*=PtclWarp.get_Jacobian(iptcl,ipsi);
    }
  }
  //Save the Log of the jacobian
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    NewBead->Properties(ipsi,LOGJACOB)=std::log(std::abs(Jacobian[ipsi]));
  //Compute deltaR : necessary to compute transition probability
  deltaR= NewBead->R - head->R;
  //Compute Warped_deltaR : necessary to compute actions
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    Warped_deltaR[ipsi]=WW[ipsi]->R - Warped_R[ipsi];
  }
  //Compute HEAD action in forward direction
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    gRand = Warped_deltaR[ipsi]-Tau*(*head->Gradients[ipsi]);
    //head->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(gRand,gRand);
    head->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(gRand,gRand);
  }
  //evaluate all relevant quantities in the new position
  int totbeadwgt(0);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    WW[ipsi]->update();
    RealType* restrict NewBeadProp=NewBead->getPropertyBase(ipsi);
    //evaluate Psi and H
    NewBeadProp[LOGPSI]=Psi1[ipsi]->evaluateLog(*WW[ipsi]);
    NewBeadProp[SIGN]=Psi1[ipsi]->getPhase();
    RealType eloc=NewBeadProp[LOCALENERGY]= H1[ipsi]->evaluate(*WW[ipsi]);
    //Save properties
    H1[ipsi]->saveProperty(NewBeadProp);
    //*(NewBead->Gradients[ipsi])= WW[ipsi]->G;
    Copy(WW[ipsi]->G,*(NewBead->Gradients[ipsi]));
    //Compute the backward part of the Kinetic action
    //gRand=Warped_deltaR[ipsi]+Tau*WW[ipsi]->G;
    PAOps<RealType,DIM>::axpy(Tau,WW[ipsi]->G,Warped_deltaR[ipsi],gRand);
    NewBead->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(gRand,gRand);
    NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
    //int beadwgt=std::abs( ( int(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
    NewBead->BeadSignWgt[ipsi]=beadwgt;
    totbeadwgt+=beadwgt;
  }
  //Compute Drift and TransProb here. This could be done after acceptance
  //if # of beads is greater than 1 but needs to be done here to make the
  //one-bead case working ok.
  NewBead->getDrift(Jacobian,PtclWarp,branchEngine->LogNorm);
  gRand=deltaR+Tau*NewBead->Drift;
  NewBead->TransProb[backward]=m_oneover2tau*Dot(gRand,gRand);
  RealType AcceptProb(-1.0),NewGlobalWgt(0.0);
  RealType RefAction(-1.0e20);
  if(totbeadwgt!=0)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      DeltaG[ipsi]= - head->Action(ipsi,forward)       - NewBead->Action(ipsi,backward)
                    + tail->Action(ipsi,forward)       + next->Action(ipsi,backward)
                    - head->Action(ipsi,Directionless) - NewBead->Action(ipsi,Directionless)
                    + tail->Action(ipsi,Directionless) + next->Action(ipsi,Directionless)
                    - head->Properties(ipsi,LOGPSI)    + NewBead->Properties(ipsi,LOGPSI)
                    - tail->Properties(ipsi,LOGPSI)    + next->Properties(ipsi,LOGPSI)
                    - tail->Properties(ipsi,LOGJACOB)  + NewBead->Properties(ipsi,LOGJACOB);
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
  if(Random() < AcceptProb)
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
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Warped_R[ipsi]=WW[ipsi]->R;
    ++nAccept;
  }
  else
  {
    ++nReject;
    ++NumTurns;
    //Evaluate initial position with warp
    W.R=tail->R;
    W.update();
    for(int iptcl=0; iptcl< nptcl; iptcl++)
    {
      PtclWarp.warp_one(iptcl,0);
      //Save particle position
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        (Warped_R[ipsi])[iptcl]=tail->R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
      }
    }
    Reptile->flip();
  }
}
}
