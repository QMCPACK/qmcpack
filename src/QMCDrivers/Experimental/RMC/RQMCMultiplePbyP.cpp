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
    
    


#include "QMCDrivers/RQMCMultiplePbyP.h"
#include "QMCDrivers/MultiChain.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Estimators/CSPolymerEstimator.h"
#include "Estimators/MJPolymerEstimator.h"
#include "Estimators/HFDHE2PolymerEstimator.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/Bead_ParticleSet.h"
#include "LongRange/StructFact.h"
namespace qmcplusplus
{
RQMCMultiplePbyP::RQMCMultiplePbyP(MCWalkerConfiguration& w,
                                   TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h),
  ReptileLength(21),
  NumTurns(0), Reptile(0), NewBead(0),
  multiEstimator(0), MSS(1.0)
{
  RootName = "rmc";
  QMCType ="RQMCMultiplePbyP";
  m_param.add(ReptileLength,"chains","int");
  MaxLevel=-1;
  m_param.add(MaxLevel,"MaxLevel","int");
  QMCDriverMode.set(QMC_MULTIPLE,1);
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
  add_H_and_Psi(&h,&psi);
  nObs = h.sizeOfObservables();
}

RQMCMultiplePbyP::~RQMCMultiplePbyP()
{
  if(Reptile)
    delete Reptile;
  delete NewBead;
}

/** main function to perform RQMC
 */
bool RQMCMultiplePbyP::run()
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
      //	moveReptile();
      moveReptile_bisection();
      Estimators->accumulate(W);
      if (block>3)
      {
        moveReptile_displace();
        Estimators->accumulate(W);
      }
      moveReptile_bisection_end();
      for (int i=0; i<ReptileLength; i++)
      {
        (*Reptile)[i]->stepmade=(*Reptile)[i]->ptclAge[0];
        for (int ptcl=0; ptcl<W.R.size(); ptcl++)
        {
          if ((*Reptile)[i]->ptclAge[ptcl]<=(*Reptile)[i]->stepmade)
          {
            (*Reptile)[i]->stepmade=(*Reptile)[i]->ptclAge[ptcl];
            //              std::cerr <<(*Reptile)[i]->ptclAge[ptcl]<<" "<<ptcl<<" "<<i<< std::endl;
          }
        }
      }
      Estimators->accumulate(W);
      Reptile->Age +=1;
      int min=(*Reptile)[0]->stepmade;
      for (int i=0; i<ReptileLength; i++)
        if ((*Reptile)[i]->stepmade<min)
          min=(*Reptile)[i]->stepmade;
      std::cerr <<Reptile->Age-min<< std::endl;
      step++;
      CurrentStep++;
    }
    while(step<nSteps);
    multiEstimator->evaluateDiff();
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAccept = 0;
    nReject = 0;
    block++;
    Reptile->record();
    //recordBlock(block);
  }
  while(block<nBlocks);
  Estimators->stop();
  Reptile->close();
  return finalize(nBlocks);
}

void RQMCMultiplePbyP::initReptile()
{
  //Resize working arrays
  resizeArrays(Psi1.size());
  //Initial direction of growth. To be read if calculation is restarted.
  int InitialGrowthDirection(0);
  //Build NewBead. This takes care of a bunch of resizing operation and properties of the starting bead
  NewBead=new Bead(**W.begin());
  Buffer_t& w_buffer(NewBead->DataSet);
  {
    W.R = NewBead->R;
    w_buffer.clear();
    W.registerData(w_buffer);
    //W.update();
    //    W.registerData(w_buffer);
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      RealType logpsi=Psi1[ipsi]->registerData(W,w_buffer);
      RealType eloc=H1[ipsi]->registerData(W,w_buffer);
      NewBead->Properties(ipsi,LOGPSI) =logpsi;
      NewBead->Properties(ipsi,SIGN) = Psi1[ipsi]->getPhase();
      ///        RealType eloc= H1[ipsi]->evaluate(W);
      NewBead->Properties(ipsi,LOCALENERGY)= eloc;
      //      Psi1[ipsi]->registerbData(W,w_buffer);
      H1[ipsi]->saveProperty(NewBead->getPropertyBase(ipsi));
      //*(NewBead->Gradients[ipsi])=W.G;
      Copy(W.G,*(NewBead->Gradients[ipsi]));
      //      Copy(W.L,*(NewBead->Laplacians[ipsi]));
      *(NewBead->Laplacians[ipsi])=W.L;
      NewBead->BeadSignWgt[ipsi]=1;
      NewBead->Action(ipsi,MinusDirection)= 0.25*Tau*Dot(*NewBead->Gradients[ipsi],*NewBead->Gradients[ipsi]);// not hack any more
      //      NewBead->Action(ipsi,MinusDirection)= 0.0;
      NewBead->Action(ipsi,PlusDirection)=NewBead->Action(ipsi,MinusDirection);
      NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    }
    NewBead->getDrift(branchEngine->LogNorm);
    NewBead->TransProb[MinusDirection]=(0.5*Tau)*Dot(NewBead->Drift,NewBead->Drift) ;
    NewBead->TransProb[PlusDirection]=NewBead->TransProb[MinusDirection];
  }
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
  InitBisection(w_buffer);
} // END OF InitReptile




void RQMCMultiplePbyP::InitBisection(Buffer_t& w_buffer)
{
  grand_transProb.resize(W.G.size());//dothis better
  lambda=1.0/(2.0*MSS);
  assert(lambda==0.5); //6.059;
  std::cout <<"My max level is "<<MaxLevel<< std::endl;
  assert(MaxLevel!=-1);
  //    MaxLevel=2;
  int num_bisection_slices=std::pow(2.0,(double)MaxLevel);
  tempReptile.resize(num_bisection_slices+1);
  tempReptile2.resize(ReptileLength);
  psiReptile.resize(num_bisection_slices+1);
  hReptile.resize(num_bisection_slices+1);
  psiReptile2.resize(ReptileLength);
  hReptile2.resize(ReptileLength);
  //    for(int i=0; i<=num_bisection_slices;i++){
  for(int i=0; i<ReptileLength; i++)
  {
    char pname[16];
    sprintf(pname,"Bead_PS.%s.c%i",W.getName().c_str(),i);
    tempReptile2[i]=new Bead_ParticleSet(W,nPsi);
    tempReptile2[i]->setName(pname);
    Bead_ParticleSet &tempBPS(*tempReptile2[i]);
    //      sprintf(pname,"Bead_PS_H.%s.c%i",W.getName().c_str(),i);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile2[i].push_back(Psi1[ipsi]->makeClone(tempBPS));
      hReptile2[i].push_back(H1[ipsi]->makeClone(tempBPS,*psiReptile2[i][ipsi]));
      //	hReptile[i][ipsi]->setName(pname);//new
    }
  }
  for(int i=0; i<=num_bisection_slices; i++)
  {
    char pname[16];
    sprintf(pname,"Bead_PS.%s.c%i",W.getName().c_str(),i);
    tempReptile[i]=tempReptile2[i];
    psiReptile[i].resize(nPsi);
    hReptile[i].resize(nPsi);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile[i][ipsi]=psiReptile2[i][ipsi];
      hReptile[i][ipsi]=hReptile2[i][ipsi];
    }
  }
  RegisterMe();
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::BisectionMove"));   //0
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::GrowMove"));        //1
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::PsiCalc"));         //2
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::HCalc"));           //3
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::BisectionStartup"));//4
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::BisectionMoves"));  //5
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::other"));           //6
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::highLevel"));       //7
  myTimers.push_back(new NewTimer("RQMCMultiplePbyP::sumProperties"));  //8
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[2]);
  TimerManager.addTimer(myTimers[3]);
  TimerManager.addTimer(myTimers[4]);
  TimerManager.addTimer(myTimers[5]);
  TimerManager.addTimer(myTimers[6]);
  TimerManager.addTimer(myTimers[7]);
  TimerManager.addTimer(myTimers[8]);
}


void RQMCMultiplePbyP::RegisterMe()
{
  //resize tempReptile info
  for (int i=0; i<tempReptile2.size(); ++i)
  {
    Walker_t::Buffer_t tbuffer;
    tempReptile2[i]->registerData(tbuffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile2[i][ipsi]->registerData(*tempReptile2[i],tbuffer);
      RealType eloc=hReptile2[i][ipsi]->registerData(*tempReptile2[i],tbuffer);
    }
  }
}

//   void RQMCMultiplePbyP::InitBisection(Buffer_t& w_buffer)
//   {

//     grand_transProb.resize(W.G.size());//dothis better
//     lambda=1.0/(2.0*MSS);
//     assert(lambda==0.5); //6.059;
//     std::cout <<"My max level is "<<MaxLevel<< std::endl;
//     assert(MaxLevel!=-1);
//     //    MaxLevel=2;
//     int num_bisection_slices=std::pow(2.0,(double)MaxLevel);
//     tempReptile.resize(num_bisection_slices+1);
//     tempReptile2.resize(ReptileLength);
//     //    tempReptile_slow.resize(num_bisection_slices+1);
//     psiReptile.resize(num_bisection_slices+1);
//     hReptile.resize(num_bisection_slices+1);
//     for(int i=0; i<=num_bisection_slices;i++){

//       char pname[16];
//       sprintf(pname,"Bead_PS.%s.c%i",W.getName().c_str(),i);
//       tempReptile[i]=new Bead_ParticleSet(W,nPsi);
//       tempReptile[i]->setName(pname);
//       Bead_ParticleSet &tempBPS(*tempReptile[i]);
//       //      sprintf(pname,"Bead_PS_H.%s.c%i",W.getName().c_str(),i);
//       for (int ipsi=0;ipsi<nPsi;ipsi++){
// 	psiReptile[i].push_back(Psi1[ipsi]->makeClone(tempBPS));
// 	hReptile[i].push_back(H1[ipsi]->makeClone(tempBPS,*psiReptile[i][ipsi]));
// 	//	hReptile[i][ipsi]->setName(pname);//new
//       }
//     }
//     RegisterMe();

//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::BisectionMove"));   //0
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::GrowMove"));        //1
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::PsiCalc"));         //2
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::HCalc"));           //3
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::BisectionStartup"));//4
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::BisectionMoves"));  //5
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::other"));           //6
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::highLevel"));       //7
//     myTimers.push_back(new NewTimer("RQMCMultiplePbyP::sumProperties"));  //8
//     TimerManager.addTimer(myTimers[1]);
//     TimerManager.addTimer(myTimers[0]);
//     TimerManager.addTimer(myTimers[2]);
//     TimerManager.addTimer(myTimers[3]);
//     TimerManager.addTimer(myTimers[4]);
//     TimerManager.addTimer(myTimers[5]);
//     TimerManager.addTimer(myTimers[6]);
//     TimerManager.addTimer(myTimers[7]);
//     TimerManager.addTimer(myTimers[8]);

//   }

//   void RQMCMultiplePbyP::RegisterMe()
//   {

//     //resize tempReptile info
//     for (int i=0;i<tempReptile.size();++i){
//       Walker_t::Buffer_t tbuffer;
//       tempReptile[i]->registerData(tbuffer);
//       for (int ipsi=0;ipsi<nPsi;ipsi++){

// 	psiReptile[i][ipsi]->registerData(*tempReptile[i],tbuffer);
// 	RealType eloc=hReptile[i][ipsi]->registerData(*tempReptile[i],tbuffer);
//       }
//     }

//   }

void
RQMCMultiplePbyP::ActionChange_displace(bool oldData)
{
  RealType AddPsi=0.0;
  int toMult;
  if (oldData)
    toMult=1;
  else
    toMult=-1;
  for (int ipsi=0; ipsi<nPsi; ipsi++)
  {
    for (int i=0; i<tempReptile2.size()-1; i++)
    {
      DeltaG[ipsi]+= tempReptile2[i]->Action(ipsi,forward)*toMult;
      AddPsi+= tempReptile2[i]->Action(ipsi,forward)*toMult;
    }
    for (int i=1; i<tempReptile2.size(); i++)
    {
      DeltaG[ipsi]+= tempReptile2[i]->Action(ipsi,backward)*toMult;
      AddPsi+= tempReptile2[i]->Action(ipsi,backward)*toMult;
    }
    for (int i=0; i<tempReptile2.size(); i++)
    {
      DeltaG[ipsi]+= 2.0*tempReptile2[i]->Action(ipsi,Directionless)*toMult;
      //	AddPsi+= 2.0*tempReptile2[i]->Action(ipsi,Directionless)*toMult;
    }
    for (int i=0; i<tempReptile2.size(); i++)
    {
      DeltaSign[ipsi] -= tempReptile2[i]->BeadSignWgt[ipsi]*toMult;
    }
    DeltaG[ipsi]-=tempReptile2[0]->Properties(ipsi,LOGPSI)*toMult;
    //      AddPsi-=tempReptile2[0]->Properties(ipsi,LOGPSI)*toMult;
    DeltaG[ipsi]-=tempReptile2[0]->Action(ipsi,Directionless)*toMult; //correct the *2
    //      AddPsi-=tempReptile2[0]->Action(ipsi,Directionless)*toMult; //correct the *2
    DeltaG[ipsi]-=tempReptile2[tempReptile2.size()-1]->Properties(ipsi,LOGPSI)*toMult;
    //      AddPsi-=tempReptile2[tempReptile2.size()-1]->Properties(ipsi,LOGPSI)*toMult;
    DeltaG[ipsi]-=tempReptile2[tempReptile2.size()-1]->Action(ipsi,Directionless)*toMult; //correct the *2
    //      AddPsi-=tempReptile2[tempReptile2.size()-1]->Action(ipsi,Directionless)*toMult; //correct the *2
  }
  //    std::cerr <<"ADDPsi is "<<AddPsi<< std::endl;
  //    std::cerr <<"Delta G is "<<DeltaG<< std::endl;
}



void RQMCMultiplePbyP::moveReptile_displace()
{
  myTimers[0]->start();
  int oldForward=forward;
  int oldBackward=backward;
  forward=PlusDirection;
  backward=MinusDirection;
  myTimers[4]->start();
  startSlice=0;
  endSlice=ReptileLength-1;
  for (int i=startSlice; i<=endSlice; i++)
    tempReptile2[i-startSlice]->CopyFromBead(*(*Reptile)[i],psiReptile2[i-startSlice]);
  for (int i=startSlice,ii=0; i<=endSlice; i++,ii++)
  {
    //      tempReptile2[ii]->R=(*Reptile)[i]->R;
    Walker_t& thisWalker(*((*Reptile)[i]));
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    w_buffer.rewind();
    //a
    tempReptile2[ii]->updateBuffer(thisWalker,w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      RealType logPsi=psiReptile2[ii][ipsi]->updateBuffer(*tempReptile2[ii],w_buffer,true);
      hReptile2[ii][ipsi]->updateBuffer(*tempReptile2[ii],w_buffer);
      assert(logPsi=tempReptile2[ii]->getPropertyBase(ipsi)[LOGPSI]);
      Copy(tempReptile2[ii]->G,*(tempReptile2[ii]->Gradients[ipsi]));
      *(tempReptile2[ii]->Laplacians[ipsi])=tempReptile2[ii]->L;
    }
    w_buffer.rewind();
//a
    tempReptile2[ii]->copyFromBuffer(w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile2[ii][ipsi]->copyFromBuffer(*tempReptile2[ii],w_buffer);
      hReptile2[ii][ipsi]->copyFromBuffer(*tempReptile2[ii],w_buffer);
    }
  }
  myTimers[4]->stop();
  ///the invariant here is that you have the laplacians
  ///and gradients that stay correct because you will need
  ///them from previous configurations
  int totalNumParticles=W.getTotalNum();
  for (int k=0; k<1; k++)
    for (int movePtcl=0; movePtcl<totalNumParticles; movePtcl++)
    {
      bool shouldReject=false;
      myTimers[7]->start();
      myTimers[5]->start();
      DeltaG=0.0;
      DeltaSign=0.0;
      //old is added
      ///calculate the old action
      ActionChange_displace(true);
      ///save the old data
      for (int i=0; i<tempReptile2.size(); i++)
        tempReptile2[i]->SaveOldData();
      RealType sigma=3.0;
      ParticleSet::ParticlePos_t toMove_rand(1);
      for (int dim=0; dim<DIM; dim++)
        toMove_rand[0][dim]=sigma*(Random()-0.5);
      //	std::cerr <<toMove_rand[0][0]<< std::endl;
      for (int i=0; i<tempReptile2.size(); i++)
        tempReptile2[i]->makeMove(movePtcl,toMove_rand[0]);
      for (int i=0; i<tempReptile2.size(); i++)
      {
        std::vector<TrialWaveFunction*> &Psi_curr(psiReptile2[i]);
        std::vector<QMCHamiltonian*> &H_curr(hReptile2[i]);
        for (int ipsi=0; ipsi<nPsi; ipsi++)
        {
          RealType* restrict NewBeadProp=tempReptile2[i]->getPropertyBase(ipsi);
          tempReptile2[i]->SetGradientAndLaplacian(ipsi);
          RealType logPsi;
          RealType eloc;
          RealType oldLogPsi=tempReptile2[i]->getPropertyBase(ipsi)[LOGPSI];
          myTimers[2]->start();
          RealType ratio=Psi_curr[ipsi]->ratio(*tempReptile2[i],movePtcl,dG,dL);
          myTimers[2]->stop();
          logPsi=std::log(ratio)+oldLogPsi;
          NewBeadProp[LOGPSI]=logPsi;
          myTimers[3]->start();
          tempReptile2[i]->dG=dG;
          tempReptile2[i]->dL=dL;
          eloc=NewBeadProp[LOCALENERGY]=hReptile2[i][ipsi]->evaluatePbyP(*tempReptile2[i],movePtcl);
          NewBeadProp[LOCALPOTENTIAL]=hReptile2[i][ipsi]->getLocalPotential();
          hReptile2[i][ipsi]->saveProperty(NewBeadProp);
          myTimers[3]->stop();
          tempReptile2[i]->G=tempReptile2[i]->G+dG;
          tempReptile2[i]->L=tempReptile2[i]->L+dL;
          Copy(tempReptile2[i]->G,*(tempReptile2[i]->Gradients[ipsi]));
          *(tempReptile2[i]->Laplacians[ipsi])=tempReptile2[i]->L;
          if (ratio<0)
            NewBeadProp[SIGN]=std::abs(NewBeadProp[SIGN]-M_PI);
          myTimers[8]->start();
          tempReptile2[i]->Action(ipsi,Directionless)=0.5*Tau*eloc;
          tempReptile2[i]->getDrift(branchEngine->LogNorm);//a
          if (i!=0)
          {
            deltaR=(*tempReptile2[i-1]).R-(*tempReptile2[i]).R; //PBC???
            W.applyMinimumImage(deltaR);
            grand_transProb=deltaR-Tau*tempReptile2[i]->Drift;//a
            tempReptile2[i]->TransProb[MinusDirection]=Dot(grand_transProb,grand_transProb)*0.5/Tau;//a
            deltaR=deltaR-Tau*tempReptile2[i]->G; //not hack any more
            tempReptile2[i]->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(deltaR,deltaR);
            //	      std::cerr <<"backward: "<<tempReptile2[i]->Action(ipsi,backward)<<" ";
            //	      std::cerr <<tempReptile2[i]->Action_saved(ipsi,backward)<< std::endl;
            //	      std::cerr <<deltaR<< std::endl;
          }
          if (i!=tempReptile2.size()-1)
          {
            deltaR=(*tempReptile2[i+1]).R-(*tempReptile2[i]).R;
            W.applyMinimumImage(deltaR);
            grand_transProb=deltaR-Tau*tempReptile2[i]->Drift;//a
            tempReptile2[i]->TransProb[PlusDirection]=Dot(grand_transProb,grand_transProb)*0.5/Tau;//a
            deltaR=deltaR-Tau*tempReptile2[i]->G; //not hack any more
            tempReptile2[i]->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(deltaR,deltaR);
            //	      std::cerr <<"forward: "<<tempReptile2[i]->Action(ipsi,forward)<<" ";
            //	      std::cerr <<tempReptile2[i]->Action_saved(ipsi,forward)<< std::endl;
            //	      std::cerr <<deltaR<< std::endl;
          }
          //checks to see if the refsign and you are the same
          //if refSign==newBeadProb[sign] then beadwgt=1  else beadwgt=0
          int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
          tempReptile2[i]->BeadSignWgt[ipsi]=beadwgt;
          if (NewBeadProp[LOCALENERGY]<= KEcut) // && (Reptile->Age>100))
            shouldReject=false;
        }
        myTimers[8]->stop();
      }
      myTimers[6]->start();
      ActionChange_displace(false);
      /// Here we should check to see if the Hamiltonian for all the reptiles
      /// has managed to go over the nodes
      bool someHamiltonian_ok=false;
      for (int ipsi=0; ipsi<nPsi; ipsi++)
      {
        bool myHamiltonian_ok=true;
        for (int slice=0; slice<tempReptile2.size(); slice++)
        {
          myHamiltonian_ok = tempReptile2[slice]->BeadSignWgt[ipsi]==1 && myHamiltonian_ok;
        }
        someHamiltonian_ok = someHamiltonian_ok || myHamiltonian_ok;
      }
      ///If some hamiltonian has overlap then we are ok
      RealType AcceptProb(-1.0),NewGlobalWgt(0.0);
      RealType RefAction(-1.0e20);
      if (someHamiltonian_ok)
        //	  std::cerr <<"NOT OVER NODES"<< std::endl;
        if (shouldReject)
          std::cerr <<"REJECTING DUE TO KE CUTOFF"<< std::endl;
      if (someHamiltonian_ok && !shouldReject)
      {
        for (int ipsi=0; ipsi<nPsi; ipsi++)
        {
          NewGlobalAction[ipsi]=Reptile->GlobalAction[ipsi]+DeltaG[ipsi];
          //Compute the new sign
          NewGlobalSignWgt[ipsi]=Reptile->GlobalSignWgt[ipsi]+
                                 DeltaSign[ipsi];
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
        ///These forward and backward probs must be changed!
        AcceptProb=std::exp(NewGlobalWgt - Reptile->GlobalWgt );
        //	  std::cerr <<"My accept prob is "<<AcceptProb<< std::endl;
        //	  std::cerr <<NewGlobalWgt - Reptile->GlobalWgt - logSampleProb + oldLogSampleProb<<" "<<rdist2_new-rdist2_old<< std::endl;
      }
      if(Random() < AcceptProb )
      {
        //	  std::cerr <<"ACCEPT"<< std::endl;
        for (int i=0; i<tempReptile2.size(); i++)
        {
          std::vector<TrialWaveFunction*> &Psi_curr(psiReptile2[i]);
          tempReptile2[i]->acceptMove(movePtcl);
          (*Reptile)[i+startSlice]->ptclAge[movePtcl]=Reptile->Age;
          for (int ipsi=0; ipsi<nPsi; ipsi++)
          {
            Psi_curr[ipsi]->acceptMove(*(ParticleSet*)tempReptile2[i],movePtcl);
            hReptile2[i][ipsi]->acceptMove(movePtcl);
          }
        }
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
        ++nAccept;
      }
      else
      {
        //	  std::cerr <<"REJECT;"<< std::endl;
        for (int i=0; i<tempReptile2.size(); i++)
        {
          std::vector<TrialWaveFunction*> &Psi_curr(psiReptile2[i]);
          tempReptile2[i]->rejectMove(movePtcl);
          tempReptile2[i]->RestoreOldData();
          for (int ipsi=0; ipsi<nPsi; ipsi++)
          {
            Psi_curr[ipsi]->rejectMove(movePtcl);
            hReptile2[i][ipsi]->rejectMove(movePtcl);
          }
        }
        //	  ++nAccept;
        ++nReject;
      }
      myTimers[6]->stop();
      myTimers[7]->stop();
    }
  for (int i=startSlice; i<=endSlice; i++)
  {
    tempReptile2[i-startSlice]->CopyToBead(*(*Reptile)[i],psiReptile2[i-startSlice]);
  }
  for (int i=startSlice,ii=0; i<=endSlice; i++,ii++)
  {
    Walker_t& thisWalker(*((*Reptile)[i]));
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    w_buffer.rewind();
    tempReptile2[ii]->copyToBuffer(w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile2[ii][ipsi]->evaluateLog(*tempReptile2[ii],w_buffer);
      //psiReptile2[ii][ipsi]->evaluate(*tempReptile2[ii],w_buffer);
      hReptile2[ii][ipsi]->evaluate(*tempReptile2[ii],w_buffer);
    }
  }
  oldForward=forward;
  oldBackward=backward;
  myTimers[0]->stop();
}



void RQMCMultiplePbyP::ChooseSlices(int ReptileLength,int &startSlice,int &endSlice)
{
  int numSlice=1;
  for (int i=0; i<MaxLevel; i++)
    numSlice*=2;
  startSlice=std::floor(Random()*(ReptileLength-numSlice));
  endSlice=startSlice+numSlice;
  assert(endSlice<ReptileLength);
}

RQMCMultiplePbyP::RealType RQMCMultiplePbyP::LogSampleProb(std::vector<Bead_ParticleSet*> &tempReptile,
    int startSlice, int endSlice,
    std::vector<int> &particles, int level)
{
  int NumImage = 0;
  PosType rpp;
  int skip = 1<<(level);
  RealType logSampleProb=0.0;
  RealType levelTau = Tau*skip*0.5;
  for (int ptclIndex=0; ptclIndex<particles.size(); ptclIndex++)
  {
    int ptcl = particles[ptclIndex];
    RealType sigma2=(1.0*lambda*levelTau);
    RealType sigma=std::sqrt(sigma2);
    RealType prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=0; slice<tempReptile.size()-1; slice+=skip)
    {
      PosType r = tempReptile[slice]->R[ptcl];  // Path(slice,ptcl);
      PosType rdiff = tempReptile[slice+skip]->R[ptcl]-tempReptile[slice]->R[ptcl];
      W.Lattice.applyMinimumImage(rdiff);
      PosType rpp=tempReptile[slice+(skip>>1)]->R[ptcl];
      W.Lattice.applyMinimumImage(r);
      W.Lattice.applyMinimumImage(rpp);
      ///We've ignored boundary conditions here (well we think this is fixed but we're not sure)
      PosType rbar=r + 0.5*rdiff;
      W.Lattice.applyMinimumImage(rbar);
      PosType Delta= rpp - rbar;
      W.Lattice.applyMinimumImage(Delta);
      RealType GaussProd=1.0;
      for (int dim=0; dim<DIM; dim++)
      {
        RealType GaussSum = 0.0;
        for (int image=-NumImage; image <= NumImage; image++)
        {
          RealType dist = Delta[dim]+(RealType)image*W.Lattice.R(dim,dim);
          GaussSum += std::exp(-0.5*dist*dist/sigma2);
        }
        GaussProd *= GaussSum;
      }
      logSampleProb += prefactorOfSampleProb + std::log(GaussProd);
    }
  }
  return logSampleProb;
}

RQMCMultiplePbyP::RealType
RQMCMultiplePbyP::BisectionSample(std::vector<Bead_ParticleSet*> &tempReptile,
                                  int particleToMove,
                                  int level)
{
  //    std::cerr <<"Bisection sample called"<< std::endl;
  int NumImage = 0;
  //WARNING: This is masking higher level gRand on purpose
  ///SHOULD JSUT RENAME
  ParticleSet::ParticlePos_t gRand;
  gRand.resize(1);
  PosType rpp;
  int skip = 1<<(level);
  RealType logNewSampleProb=0.0;
  RealType levelTau = Tau*skip*0.5;
  int  ptcl=particleToMove;
  RealType sigma2=(1.0*lambda*levelTau);
  RealType sigma=std::sqrt(sigma2);
  RealType prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
  for (int slice=0; slice<tempReptile.size()-1; slice+=skip)
  {
    PosType r = tempReptile[slice]->R[ptcl];
    W.Lattice.applyMinimumImage(r);
    PosType rpp;
    PosType rdiff=tempReptile[slice+skip]->R[ptcl]-tempReptile[slice]->R[ptcl];
    W.Lattice.applyMinimumImage(rdiff);
    PosType rbar = r + 0.5*rdiff;
    PosType newDelta;
    makeGaussRandom(gRand);
    newDelta=sigma*gRand[0];
    W.Lattice.applyMinimumImage(newDelta);
    RealType GaussProd=1.0;
    for (int dim=0; dim<DIM; dim++)
    {
      RealType GaussSum = 0.0;
      for (int image=-NumImage; image <= NumImage; image++)
      {
        RealType dist = newDelta[dim]+(RealType)image*W.Lattice.R(dim,dim);
        GaussSum += std::exp(-0.5*dist*dist/sigma2);
      }
      GaussProd *= GaussSum;
    }
    logNewSampleProb += prefactorOfSampleProb + std::log(GaussProd);
    rpp=rbar+newDelta;
    PosType toMove= rpp-tempReptile[slice+(skip>>1)]->R[ptcl];
    //      std::cerr <<slice+(skip>>1)<<" "<<toMove<< std::endl;
    ///Here we've stored the new position in the path
    //      tempReptile[slice+(skip>>1)]->R[ptcl]=rpp;
    W.Lattice.applyMinimumImage(toMove);
    //      std::cerr <<toMove<< std::endl;
    PosType newPlace=tempReptile[slice+(skip>>1)]->makeMove(ptcl,toMove);
    //tempReptile[slice+(skip>>1)]->R[ptcl]=rpp;
    //      assert(newPlace==rpp);
    //      tempReptile[0]->makeMove(ptcl,toMove);
  }
  return logNewSampleProb;
}


void
RQMCMultiplePbyP::ActionChange_wf(bool oldData,int sliceToThrow)
{
  int toMult;
  RealType AddPsi=0.0;
  if (oldData)
    toMult=1;
  else
    toMult=-1;
  for (int ipsi=0; ipsi<nPsi; ipsi++)
  {
    for (int i=0; i<tempReptile.size()-1; i++)
    {
      DeltaG[ipsi]+= tempReptile[i]->Action(ipsi,forward)*toMult;
      AddPsi+= tempReptile[i]->Action(ipsi,forward)*toMult;
    }
    for (int i=1; i<tempReptile.size(); i++)
    {
      DeltaG[ipsi]+= tempReptile[i]->Action(ipsi,backward)*toMult;
      AddPsi+= tempReptile[i]->Action(ipsi,backward)*toMult;
    }
    for (int i=0; i<tempReptile.size(); i++)
    {
      DeltaG[ipsi]+= 2.0*tempReptile[i]->Action(ipsi,Directionless)*toMult;
    }
    for (int i=0; i<tempReptile.size(); i++)
    {
      DeltaSign[ipsi] -= tempReptile[i]->BeadSignWgt[ipsi]*toMult;
    }
    DeltaG[ipsi]-=tempReptile[sliceToThrow]->Properties(ipsi,LOGPSI)*toMult;
    DeltaG[ipsi]-=tempReptile[sliceToThrow]->Action(ipsi,Directionless)*toMult; //correct the *2
  }
  //    std::cerr <<"Add psi for wf is "<<AddPsi<< std::endl;
  //    std::cerr <<"Delta G for wf is "<<DeltaG<< std::endl;
}



void
RQMCMultiplePbyP::PutInBox (PosType &v)
{
  std::cerr <<"My box is of size "<< W.Lattice.R(0,0)<<" "<<W.Lattice.R(1,1)<<" "<<W.Lattice.R(2,2)<<" "<<W.Lattice.R(1,2)<< std::endl;;
  for (int i=0; i<DIM; i++)
  {
    double boxLength=W.Lattice.R(i,i);
    double n = -std::floor(v(i)*(1.0/boxLength)+0.5);
    v(i) += n*boxLength;
  }
}




void RQMCMultiplePbyP::moveReptile_bisection_end()
{
  myTimers[0]->start();
  int oldForward=forward;
  int oldBackward=backward;
  forward=PlusDirection;
  backward=MinusDirection;
  //    int startSlice;
  //    int endSlice;
  //    int sliceToThrow;
  int numSlice=1;
  for (int i=0; i<MaxLevel; i++)
    numSlice*=2;
  if (Random()>0.5)
  {
    startSlice=0;
    endSlice=startSlice+numSlice;
    sliceToThrow=0;
  }
  else
  {
    endSlice=ReptileLength-1;
    startSlice=endSlice-numSlice;
    sliceToThrow=tempReptile.size()-1;
  }
  myTimers[4]->start();
  for (int i=startSlice; i<=endSlice; i++)
    tempReptile[i-startSlice]->CopyFromBead(*(*Reptile)[i],psiReptile[i-startSlice]);
  for (int i=startSlice,ii=0; i<=endSlice; i++,ii++)
  {
    //      tempReptile[ii]->R=(*Reptile)[i]->R;
    Walker_t& thisWalker(*((*Reptile)[i]));
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    w_buffer.rewind();
    //a
    tempReptile[ii]->updateBuffer(thisWalker,w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      RealType logPsi=psiReptile[ii][ipsi]->updateBuffer(*tempReptile[ii],w_buffer,true);
      hReptile[ii][ipsi]->updateBuffer(*tempReptile[ii],w_buffer);
      assert(logPsi=tempReptile[ii]->getPropertyBase(ipsi)[LOGPSI]);
      Copy(tempReptile[ii]->G,*(tempReptile[ii]->Gradients[ipsi]));
      *(tempReptile[ii]->Laplacians[ipsi])=tempReptile[ii]->L;
    }
    w_buffer.rewind();
//a
    tempReptile[ii]->copyFromBuffer(w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile[ii][ipsi]->copyFromBuffer(*tempReptile[ii],w_buffer);
      hReptile[ii][ipsi]->copyFromBuffer(*tempReptile[ii],w_buffer);
    }
  }
  myTimers[4]->stop();
  ///the invariant here is that you have the laplacians
  ///and gradients that stay correct because you will need
  ///them from previous configurations
  int totalNumParticles=W.getTotalNum();
  for (int k=0; k<1; k++)
    for (int movePtcl=0; movePtcl<totalNumParticles; movePtcl++)
    {
      bool shouldReject=false;
      myTimers[7]->start();
      myTimers[5]->start();
      DeltaG=0.0;
      DeltaSign=0.0;
      //old is added
      ///calculate the old action
      ActionChange_wf(true,sliceToThrow);
      ///save the old data
      for (int i=0; i<tempReptile.size(); i++)
        tempReptile[i]->SaveOldData();
      std::vector<int> particleToMove_vector(1,movePtcl);
      RealType oldLogSampleProb=0.0;
      for (int level=MaxLevel; level>0; level--)
        oldLogSampleProb+=LogSampleProb(tempReptile,
                                        startSlice,endSlice,
                                        particleToMove_vector,level);
      PosType old_rdiff=tempReptile[tempReptile.size()-1]->R[movePtcl]-tempReptile[0]->R[movePtcl];
      W.Lattice.applyMinimumImage(old_rdiff);
      int skip = 1<<(MaxLevel);
      RealType levelTau = Tau*skip*0.5;
      RealType sigma=std::sqrt(4.0*lambda*levelTau);
      PosType rand_toss;
      ParticleSet::ParticlePos_t toThrow_rand(1);
      makeGaussRandom(toThrow_rand);
      toThrow_rand=toThrow_rand*sigma;
      //	for (int dim=0;dim<DIM;dim++)
      //	  rand_toss[dim]=sigma*(Random()-0.5);
      int not_sliceToThrow;
      if (sliceToThrow==0)
        not_sliceToThrow=tempReptile.size()-1;
      else
        not_sliceToThrow=0;
      PosType desiredMove=toThrow_rand[0]+tempReptile[not_sliceToThrow]->R[movePtcl];
      PosType toMove=desiredMove-tempReptile[sliceToThrow]->R[movePtcl];
      //	tempReptile[sliceToThrow]->makeMove(movePtcl,rand_toss);
      PosType actualMove=tempReptile[sliceToThrow]->makeMove(movePtcl,toMove);
      RealType logSampleProb=0.0;
      for (int level=MaxLevel; level>0; level--)
        logSampleProb+=BisectionSample(tempReptile,movePtcl,level);
      myTimers[5]->stop();
      PosType new_rdiff=tempReptile[tempReptile.size()-1]->R[movePtcl]-tempReptile[0]->R[movePtcl];
      for (int i=0; i<tempReptile.size(); i++)
      {
        std::vector<TrialWaveFunction*> &Psi_curr(psiReptile[i]);
        std::vector<QMCHamiltonian*> &H_curr(hReptile[i]);
        for (int ipsi=0; ipsi<nPsi; ipsi++)
        {
          RealType* restrict HeadProp=(*Reptile)[0]->getPropertyBase(ipsi);
          //	    myTimers[2]->start();
          RealType* restrict NewBeadProp=tempReptile[i]->getPropertyBase(ipsi);
          tempReptile[i]->SetGradientAndLaplacian(ipsi);
          RealType logPsi;
          RealType eloc;
          if (i==sliceToThrow || (i!=0 && i!=tempReptile.size()-1))
          {
            RealType oldLogPsi=tempReptile[i]->getPropertyBase(ipsi)[LOGPSI];
            myTimers[2]->start();
            RealType ratio=Psi_curr[ipsi]->ratio(*tempReptile[i],movePtcl,dG,dL);
            myTimers[2]->stop();
            logPsi=std::log(ratio)+oldLogPsi;
            NewBeadProp[LOGPSI]=logPsi;
            myTimers[3]->start();
            tempReptile[i]->dG=dG;
            tempReptile[i]->dL=dL;
            //	      W.R=(tempReptile[i]->R);
            //	      W.G=tempReptile[i]->G+dG;
            //	      W.L=tempReptile[i]->L+dL;
            //	      W.update();
            //	      eloc=NewBeadProp[LOCALENERGY]=H1[ipsi]->evaluate(W);
            eloc=NewBeadProp[LOCALENERGY]=hReptile[i][ipsi]->evaluatePbyP(*tempReptile[i],movePtcl);
            hReptile[i][ipsi]->saveProperty(NewBeadProp);
            //	      H1[ipsi]->saveProperty(NewBeadProp);
            myTimers[3]->stop();
            //	      myTimers[2]->start();
            tempReptile[i]->G=tempReptile[i]->G+dG;
            tempReptile[i]->L=tempReptile[i]->L+dL;
            Copy(tempReptile[i]->G,*(tempReptile[i]->Gradients[ipsi]));
            *(tempReptile[i]->Laplacians[ipsi])=tempReptile[i]->L;
            if (ratio<0)
              NewBeadProp[SIGN]=std::abs(NewBeadProp[SIGN]-M_PI);
          }
          else
          {
            logPsi=tempReptile[i]->getPropertyBase(ipsi)[LOGPSI];
            eloc=NewBeadProp[LOCALENERGY];
          }
          //	    myTimers[2]->stop();
          myTimers[8]->start();
          tempReptile[i]->Action(ipsi,Directionless)=0.5*Tau*eloc;
          tempReptile[i]->getDrift(branchEngine->LogNorm);//a
          if (i!=0)
          {
            deltaR=(*tempReptile[i-1]).R-(*tempReptile[i]).R; //PBC???
            W.applyMinimumImage(deltaR);
            grand_transProb=deltaR-Tau*tempReptile[i]->Drift;//a
            tempReptile[i]->TransProb[MinusDirection]=Dot(grand_transProb,grand_transProb)*0.5/Tau;//a
            deltaR=deltaR-Tau*tempReptile[i]->G; //not hack any more
            tempReptile[i]->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(deltaR,deltaR);
          }
          if (i!=tempReptile.size()-1)
          {
            deltaR=(*tempReptile[i+1]).R-(*tempReptile[i]).R;
            W.applyMinimumImage(deltaR);
            grand_transProb=deltaR-Tau*tempReptile[i]->Drift;//a
            tempReptile[i]->TransProb[PlusDirection]=Dot(grand_transProb,grand_transProb)*0.5/Tau;//a
            deltaR=deltaR-Tau*tempReptile[i]->G; //not hack any more
            tempReptile[i]->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(deltaR,deltaR);
          }
          //checks to see if the refsign and you are the same
          //if refSign==newBeadProb[sign] then beadwgt=1  else beadwgt=0
          int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
          tempReptile[i]->BeadSignWgt[ipsi]=beadwgt;
          if (NewBeadProp[LOCALENERGY]<= KEcut) // && (Reptile->Age>100))
            shouldReject=false;
        }
        myTimers[8]->stop();
        //	    checkBeadInfo(i);
      }
      myTimers[6]->start();
      ActionChange_wf(false,sliceToThrow);
      /// Here we should check to see if the Hamiltonian for all the reptiles
      /// has managed to go over the nodes
      bool someHamiltonian_ok=false;
      for (int ipsi=0; ipsi<nPsi; ipsi++)
      {
        bool myHamiltonian_ok=true;
        for (int slice=0; slice<tempReptile.size(); slice++)
        {
          myHamiltonian_ok = tempReptile[slice]->BeadSignWgt[ipsi]==1 && myHamiltonian_ok;
        }
        someHamiltonian_ok = someHamiltonian_ok || myHamiltonian_ok;
      }
      ///If some hamiltonian has overlap then we are ok
      RealType AcceptProb(-1.0),NewGlobalWgt(0.0);
      RealType RefAction(-1.0e20);
      if (someHamiltonian_ok && !shouldReject)
      {
        for (int ipsi=0; ipsi<nPsi; ipsi++)
        {
          NewGlobalAction[ipsi]=Reptile->GlobalAction[ipsi]+DeltaG[ipsi];
          //Compute the new sign
          NewGlobalSignWgt[ipsi]=Reptile->GlobalSignWgt[ipsi]+
                                 DeltaSign[ipsi];
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
        RealType rdist2_old=dot(old_rdiff,old_rdiff)/(2.0*sigma*sigma);
        RealType rdist2_new=dot(toThrow_rand[0],toThrow_rand[0])/(2.0*sigma*sigma);
        ///These forward and backward probs must be changed!
        AcceptProb=std::exp(NewGlobalWgt - Reptile->GlobalWgt - logSampleProb + oldLogSampleProb +
                            +rdist2_new-rdist2_old);
        //	  std::cerr <<"My accept prob is "<<AcceptProb<< std::endl;
        //	  std::cerr <<NewGlobalWgt - Reptile->GlobalWgt - logSampleProb + oldLogSampleProb<<" "<<rdist2_new-rdist2_old<< std::endl;
      }
      //	std::cerr <<"ACCEpt prob is "<<AcceptProb<< std::endl;
      if(Random() < AcceptProb )
      {
        //	  std::cerr <<"ACCEPT end"<< std::endl;
        for (int i=0; i<tempReptile.size(); i++)
        {
          std::vector<TrialWaveFunction*> &Psi_curr(psiReptile[i]);
          if (i==sliceToThrow || (i!=0 && i!=tempReptile.size()-1))
          {
            tempReptile[i]->acceptMove(movePtcl);
            (*Reptile)[i+startSlice]->ptclAge[movePtcl]=Reptile->Age;
            for (int ipsi=0; ipsi<nPsi; ipsi++)
            {
              Psi_curr[ipsi]->acceptMove(*(ParticleSet*)tempReptile[i],movePtcl);
              hReptile[i][ipsi]->acceptMove(movePtcl);
            }
          }
        }
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
        //	  ++nAccept;
      }
      else
      {
        //	  std::cerr <<"REJECT end"<< std::endl;
        for (int i=0; i<tempReptile.size(); i++)
        {
          std::vector<TrialWaveFunction*> &Psi_curr(psiReptile[i]);
          if (i==sliceToThrow || (i!=0 && i!=tempReptile.size()-1))
            tempReptile[i]->rejectMove(movePtcl);
          tempReptile[i]->RestoreOldData();
          if (i==sliceToThrow || (i!=0 && i!=tempReptile.size()-1))
            for (int ipsi=0; ipsi<nPsi; ipsi++)
            {
              Psi_curr[ipsi]->rejectMove(movePtcl);
              hReptile[i][ipsi]->rejectMove(movePtcl);
            }
        }
        //	  ++nAccept;
        //	  ++nReject;
      }
      myTimers[6]->stop();
      myTimers[7]->stop();
    }
  for (int i=startSlice; i<=endSlice; i++)
  {
    tempReptile[i-startSlice]->CopyToBead(*(*Reptile)[i],psiReptile[i-startSlice]);
  }
  for (int i=startSlice,ii=0; i<=endSlice; i++,ii++)
  {
    Walker_t& thisWalker(*((*Reptile)[i]));
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    w_buffer.rewind();
    tempReptile[ii]->copyToBuffer(w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile[ii][ipsi]->evaluateLog(*tempReptile[ii],w_buffer);
      //psiReptile[ii][ipsi]->evaluate(*tempReptile[ii],w_buffer);
      hReptile[ii][ipsi]->evaluate(*tempReptile[ii],w_buffer);
    }
  }
  oldForward=forward;
  oldBackward=backward;
  myTimers[0]->stop();
}


void RQMCMultiplePbyP::moveReptile_bisection()
{
  myTimers[0]->start();
  int oldForward=forward;
  int oldBackward=backward;
  forward=PlusDirection;
  backward=MinusDirection;
  //setup the walkers intitially
  //    int startSlice;
  //    int endSlice;
  ChooseSlices(ReptileLength,startSlice,endSlice);
  bool isFirst=startSlice==0;
  bool isLast=endSlice==ReptileLength-1;
  myTimers[4]->start();
  for (int i=startSlice; i<=endSlice; i++)
    tempReptile[i-startSlice]->CopyFromBead(*(*Reptile)[i],psiReptile[i-startSlice]);
  for (int i=startSlice,ii=0; i<=endSlice; i++,ii++)
  {
    //      tempReptile[ii]->R=(*Reptile)[i]->R;
    Walker_t& thisWalker(*((*Reptile)[i]));
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    //       w_buffer.rewind();
    //       tempReptile[ii]->updateBuffer(thisWalker,w_buffer);
    //       for (int ipsi=0;ipsi<nPsi;ipsi++){
    //	 RealType logPsi=psiReptile[ii][ipsi]->updateBuffer(*tempReptile[ii],w_buffer,true);
    //	 assert(logPsi=tempReptile[ii]->getPropertyBase(ipsi)[LOGPSI]);
    //         Copy(tempReptile[ii]->G,*(tempReptile[ii]->Gradients[ipsi]));
    //	 *(tempReptile[ii]->Laplacians[ipsi])=tempReptile[ii]->L;
    //       }
    w_buffer.rewind();
    tempReptile[ii]->copyFromBuffer(w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile[ii][ipsi]->copyFromBuffer(*tempReptile[ii],w_buffer);
      hReptile[ii][ipsi]->copyFromBuffer(*tempReptile[ii],w_buffer);
    }
  }
  myTimers[4]->stop();
  ///the invariant here is that you have the laplacians
  ///and gradients that stay correct because you will need
  ///them from previous configurations
  int totalNumParticles=W.getTotalNum();
  for (int k=0; k<1; k++)
    for (int movePtcl=0; movePtcl<totalNumParticles; movePtcl++)
    {
      bool shouldReject=false;
      myTimers[7]->start();
      myTimers[5]->start();
      DeltaG=0.0;
      DeltaSign=0.0;
      //old is added
      ///calculate the old action
      if (isFirst)
        ActionChange_wf(true,0);
      else
        if (isLast)
          ActionChange_wf(true,tempReptile.size()-1);
        else
          ActionChange(true);
      ///save the old data
      for (int i=0; i<tempReptile.size(); i++)
        tempReptile[i]->SaveOldData();
      std::vector<int> particleToMove_vector(1,movePtcl);
      RealType oldLogSampleProb=0.0;
      for (int level=MaxLevel; level>0; level--)
        oldLogSampleProb+=LogSampleProb(tempReptile,
                                        startSlice,endSlice,
                                        particleToMove_vector,level);
      RealType logSampleProb=0.0;
      for (int level=MaxLevel; level>0; level--)
        logSampleProb+=BisectionSample(tempReptile,movePtcl,level);
      myTimers[5]->stop();
      for (int i=0; i<tempReptile.size(); i++)
      {
        //	  std::cerr <<i<< std::endl;
        std::vector<TrialWaveFunction*> &Psi_curr(psiReptile[i]);
        std::vector<QMCHamiltonian*> &H_curr(hReptile[i]);
        for (int ipsi=0; ipsi<nPsi; ipsi++)
        {
          RealType* restrict HeadProp=(*Reptile)[0]->getPropertyBase(ipsi);
          //	    myTimers[2]->start();
          RealType* restrict NewBeadProp=tempReptile[i]->getPropertyBase(ipsi);
          tempReptile[i]->SetGradientAndLaplacian(ipsi);
          RealType logPsi;
          RealType eloc;
          if (i!=0 && i!=tempReptile.size()-1)
          {
            RealType oldLogPsi=tempReptile[i]->getPropertyBase(ipsi)[LOGPSI];
            myTimers[2]->start();
            RealType ratio=Psi_curr[ipsi]->ratio(*tempReptile[i],movePtcl,tempReptile[i]->dG,tempReptile[i]->dL);
            myTimers[2]->stop();
            logPsi=std::log(ratio)+oldLogPsi;
            NewBeadProp[LOGPSI]=logPsi;
            myTimers[3]->start();
            eloc=NewBeadProp[LOCALENERGY]= hReptile[i][ipsi]->evaluatePbyP(*tempReptile[i],movePtcl);
            hReptile[i][ipsi]->saveProperty(NewBeadProp);
            myTimers[3]->stop();
            NewBeadProp[LOGPSI]=logPsi;
            tempReptile[i]->G=tempReptile[i]->G+tempReptile[i]->dG;
            tempReptile[i]->L=tempReptile[i]->L+tempReptile[i]->dL;
            Copy(tempReptile[i]->G,*(tempReptile[i]->Gradients[ipsi]));
            *(tempReptile[i]->Laplacians[ipsi])=tempReptile[i]->L;
            if (ratio<0)
              NewBeadProp[SIGN]=std::abs(NewBeadProp[SIGN]-M_PI);
          }
          else
          {
            logPsi=tempReptile[i]->getPropertyBase(ipsi)[LOGPSI];
            eloc=NewBeadProp[LOCALENERGY];
          }
          myTimers[8]->start();
          tempReptile[i]->Action(ipsi,Directionless)=0.5*Tau*eloc;
          tempReptile[i]->getDrift(branchEngine->LogNorm);//a
          if (i!=0)
          {
            deltaR=(*tempReptile[i-1]).R-(*tempReptile[i]).R; //PBC???
            tempReptile[i]->applyMinimumImage(deltaR);
            grand_transProb=deltaR-Tau*tempReptile[i]->Drift;//a
            tempReptile[i]->TransProb[MinusDirection]=Dot(grand_transProb,grand_transProb)*0.5/Tau;//a
            deltaR=deltaR-Tau*tempReptile[i]->G; //not hack any more
            tempReptile[i]->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(deltaR,deltaR);
          }
          if (i!=tempReptile.size()-1)
          {
            deltaR=(*tempReptile[i+1]).R-(*tempReptile[i]).R;
            tempReptile[i]->applyMinimumImage(deltaR);
            grand_transProb=deltaR-Tau*tempReptile[i]->Drift;//a
            tempReptile[i]->TransProb[PlusDirection]=Dot(grand_transProb,grand_transProb)*0.5/Tau;//a
            deltaR=deltaR-Tau*tempReptile[i]->G; //not hack any more
            tempReptile[i]->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(deltaR,deltaR);
          }
          //checks to see if the refsign and you are the same
          //if refSign==newBeadProb[sign] then beadwgt=1  else beadwgt=0
          int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
          tempReptile[i]->BeadSignWgt[ipsi]=beadwgt;
          if (NewBeadProp[LOCALENERGY] <= KEcut) // && (Reptile->Age>100))
            shouldReject=false;
        }
        myTimers[8]->stop();
        //	    checkBeadInfo(i);
      }
      myTimers[6]->start();
      if (isFirst)
        ActionChange_wf(false,0);
      else
        if (isLast)
          ActionChange_wf(false,tempReptile.size()-1);
        else
          ActionChange(false);
      /// Here we should check to see if the Hamiltonian for all the reptiles
      /// has managed to go over the nodes
      bool someHamiltonian_ok=false;
      for (int ipsi=0; ipsi<nPsi; ipsi++)
      {
        bool myHamiltonian_ok=true;
        for (int slice=0; slice<tempReptile.size(); slice++)
        {
          myHamiltonian_ok = tempReptile[slice]->BeadSignWgt[ipsi]==1 && myHamiltonian_ok;
        }
        someHamiltonian_ok = someHamiltonian_ok || myHamiltonian_ok;
      }
      ///If some hamiltonian has overlap then we are ok
      RealType AcceptProb(-1.0),NewGlobalWgt(0.0);
      RealType RefAction(-1.0e20);
      if (someHamiltonian_ok && !shouldReject)
      {
        for (int ipsi=0; ipsi<nPsi; ipsi++)
        {
          NewGlobalAction[ipsi]=Reptile->GlobalAction[ipsi]+DeltaG[ipsi];
          //Compute the new sign
          NewGlobalSignWgt[ipsi]=Reptile->GlobalSignWgt[ipsi]+
                                 DeltaSign[ipsi];
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
        ///These forward and backward probs must be changed!
        AcceptProb=std::exp(NewGlobalWgt - Reptile->GlobalWgt - logSampleProb + oldLogSampleProb);
      }
      if(Random() < AcceptProb )
      {
        //	  std::cerr <<"ACCEPT"<< std::endl;
        for (int i=0; i<tempReptile.size(); i++)
        {
          std::vector<TrialWaveFunction*> &Psi_curr(psiReptile[i]);
          if (i!=0 && i!=tempReptile.size()-1)
          {
            tempReptile[i]->acceptMove(movePtcl);
            for (int ipsi=0; ipsi<nPsi; ipsi++)
            {
              Psi_curr[ipsi]->acceptMove(*(ParticleSet*)tempReptile[i],movePtcl);
              hReptile[i][ipsi]->acceptMove(movePtcl);
            }
          }
        }
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
        //aa	  ++nAccept;
      }
      else
      {
        //	  std::cerr <<"REJECT"<< std::endl;
        for (int i=0; i<tempReptile.size(); i++)
        {
          std::vector<TrialWaveFunction*> &Psi_curr(psiReptile[i]);
          if (i!=0 && i!=tempReptile.size()-1)
            tempReptile[i]->rejectMove(movePtcl);
          tempReptile[i]->RestoreOldData();
          if (i!=0 && i!=tempReptile.size()-1)
            for (int ipsi=0; ipsi<nPsi; ipsi++)
            {
              Psi_curr[ipsi]->rejectMove(movePtcl);
              hReptile[i][ipsi]->rejectMove(movePtcl);
            }
        }
        //bb	  ++nReject;
      }
// 	for (int i=0;i<tempReptile.size();i++)
// 	  checkBeadInfo(i);
      myTimers[6]->stop();
      myTimers[7]->stop();
    }
  //    myTimers[6]->start();
  for (int i=startSlice; i<=endSlice; i++)
  {
    tempReptile[i-startSlice]->CopyToBead(*(*Reptile)[i],psiReptile[i-startSlice]);
  }
  for (int i=startSlice,ii=0; i<=endSlice; i++,ii++)
  {
    Walker_t& thisWalker(*((*Reptile)[i]));
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    w_buffer.rewind();
    tempReptile[ii]->copyToBuffer(w_buffer);
    for (int ipsi=0; ipsi<nPsi; ipsi++)
    {
      psiReptile[ii][ipsi]->evaluateLog(*tempReptile[ii],w_buffer);
      //psiReptile[ii][ipsi]->evaluate(*tempReptile[ii],w_buffer);
      hReptile[ii][ipsi]->evaluate(*tempReptile[ii],w_buffer);
    }
  }
  //myTimers[6]->stop();
  oldForward=forward;
  oldBackward=backward;
  myTimers[0]->stop();
}



/** initialize Reptile
 *
 * The actions are
 * - resize any initernal array
 * - create Reptile if never initialized
 *   -- if a previous configuration file is not "invalid", initialize the beads
 * - set reference properties for the first run
 */
void RQMCMultiplePbyP::initReptile_new()
{
  std::cerr <<"INIT REPTILE BEING CALLED"<< std::endl;
//     m_oneover2tau=0.5/Tau;
//     m_sqrttau=std::sqrt(Tau);
//     Tauoverm = Tau/MSS;
//     sqrtTauoverm = std::sqrt(Tauoverm);
  //Resize working arrays
  resizeArrays(Psi1.size());
  //Initial direction of growth. To be read if calculation is restarted.
  int InitialGrowthDirection(0);
  //Build NewBead. This takes care of a bunch of resizing operation and properties of the starting bead
  NewBead=new Bead(**W.begin());
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
    *(NewBead->Laplacians[ipsi])=W.L;
    NewBead->BeadSignWgt[ipsi]=1;
    NewBead->getScaledDrift(branchEngine->LogNorm,Tauoverm);
    NewBead->Action(ipsi,MinusDirection)= 0.5*m_oneover2tau*Dot(*NewBead->DriftVectors[ipsi],*NewBead->DriftVectors[ipsi]);
    NewBead->Action(ipsi,PlusDirection)=NewBead->Action(ipsi,MinusDirection);
    NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    NewBead->stepmade=0;
    NewBead->deltaRSquared[0]=0.0;
    NewBead->deltaRSquared[1]=0.0;
    NewBead->deltaRSquared[2]=0.0;
  }
  //     NewBead->TransProb[MinusDirection]=(0.5*Tau)*Dot(NewBead->Drift,NewBead->Drift) ;
  NewBead->TransProb[MinusDirection]=m_oneover2tau*Dot(NewBead->Drift,NewBead->Drift) ;
  NewBead->TransProb[PlusDirection]=NewBead->TransProb[MinusDirection];
  std::cerr <<"STARTINFO: "<<NewBead->TransProb[forward]<<" "<<NewBead->Drift<< std::endl;
  std::cerr <<"STARTINFO: "<<NewBead->TransProb[backward]<< std::endl;
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
  ///THis does not belong here!!
//     multiEstimator->pnorm = 1.0/( W.Lattice.DIM *  W.Lattice.Volume);
  if(!restartmode)
  {
    if(reuseReptile)
      checkReptileProperties();
    else
      setReptileProperties();
  }
} // END OF InitReptile

void RQMCMultiplePbyP::resizeArrays(int n)
{
  nPsi = n;
  nptcl=W.G.size();
  gRand.resize(nptcl);
  NewGlobalAction.resize(n);
  NewGlobalSignWgt.resize(n);
  DeltaG.resize(n);
  DeltaSign.resize(n);
  WeightSign.resize(n);
  dG.resize(nptcl);
  dL.resize(nptcl);
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
void RQMCMultiplePbyP::setReptileProperties()
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

void RQMCMultiplePbyP::checkBeadInfo(int i,bool dontDie)
{
  //    std::cerr <<"Checking bead info "<<i<< std::endl;
  W.R=tempReptile[i]->R;
  W.update();
  //    std::cerr <<"R is "<<W.R<< std::endl;
  Psi.evaluateLog(W);
  //    std::cerr <<"G is "<<W.G<< std::endl;
  tempReptile[i]->SetGradientAndLaplacian(0);
  if (!dontDie)
  {
    bool noGood=false;
    for (int k=0; k<W.R.size(); k++)
      for (int j=0; j<W.R[0].size(); j++)
      {
        if (!(std::abs(tempReptile[i]->G[k][j]-W.G[k][j])<1e-6))
        {
          std::cerr <<tempReptile[i]->G[k][j]<<" "<<W.G[k][j]<< std::endl;
          std::cerr <<tempReptile[i]->R<< std::endl;
          std::cerr <<W.Lattice.R(0,0)<< std::endl;
          noGood=true;
          exit(1);
        }
      }
    for (int k=0; k<W.R.size(); k++)
      for (int j=0; j<W.R[0].size(); j++)
      {
        //	  std::cerr <<tempReptile[i]->G[k][j]<<" "<<W.G[k][j]<< std::endl;
        assert(std::abs(tempReptile[i]->G[k][j]-W.G[k][j])<1e-6);
      }
  }
  //    std::cerr <<"tempReptile G is "<<tempReptile[i]->G<< std::endl;
  //    W.getDrift(branchEngine->LogNorm);
  //    Bead_ParticleSet &prev(tempReptile[i]);
  ParticleSet::ParticlePos_t deltaR;
  deltaR.resize(W.R.size());
  ParticleSet::ParticlePos_t deltaR_forward;
  ParticleSet::ParticlePos_t deltaR_backward;
  deltaR_forward.resize(W.R.size());
  deltaR_backward.resize(W.R.size());
  if (i!=tempReptile.size()-1)
    deltaR_forward=tempReptile[i+1]->R-tempReptile[i]->R;
  W.applyMinimumImage(deltaR_forward);
  if (i!=0)
    deltaR_backward=tempReptile[i-1]->R-tempReptile[i]->R;
  W.applyMinimumImage(deltaR_backward);
  for (int ipsi=0; ipsi<nPsi; ipsi++)
  {
    if (i!=0)
    {
      //	ParticleSet::ParticlePos_t deltaR=deltaR_backward+W.G;
      PAOps<RealType,DIM>::axpy(-Tau,W.G,deltaR_backward,deltaR);
      std::cerr <<tempReptile[i]->Action(ipsi,MinusDirection)<<" "<<0.5*m_oneover2tau*Dot(deltaR,deltaR)<< std::endl;
      assert(std::abs(tempReptile[i]->Action(ipsi,MinusDirection)-0.5*m_oneover2tau*Dot(deltaR,deltaR))<1e-6);
      std::cerr <<"Broken: "<<tempReptile[i]->TransProb[MinusDirection]<<" "<<Dot(deltaR,deltaR)*0.5/Tau<< std::endl;
      assert(std::abs(tempReptile[i]->TransProb[MinusDirection]-Dot(deltaR,deltaR)*0.5/Tau)<1e-6);
    }
    if (i!=tempReptile.size()-1)
    {
      //	ParticleSet::ParticlePos_t deltaR=deltaR_forward-Tau*W.G;
      PAOps<RealType,DIM>::axpy(-Tau,W.G,deltaR_forward,deltaR);
      //	std::cerr <<i<<" "<<tempReptile[i]->Action(ipsi,PlusDirection)<<" "<<0.5*m_oneover2tau*Dot(deltaR,deltaR)<< std::endl;
      //	std::cerr <<i<<" "<<tempReptile[i]->Action(ipsi,MinusDirection)<< std::endl;
      //	std::cerr <<" "<<tempReptile[i]->Action(ipsi,PlusDirection)<<" "<<0.5*m_oneover2tau*Dot(deltaR,deltaR)<<" "<<i<<" "<<Tau<< std::endl;
      assert(std::abs(tempReptile[i]->Action(ipsi,PlusDirection)-0.5*m_oneover2tau*Dot(deltaR,deltaR))<1e-6);
      //	deltaR=deltaR_forward-Tau*W.G; //W.Drift;
      //	std::cerr <<" "<<tempReptile[i]->TransProb[PlusDirection]<<" "<<Dot(deltaR,deltaR)*0.5/Tau<<" "<<i<<" "<<Tau<< std::endl;
      assert(std::abs(tempReptile[i]->TransProb[PlusDirection]-Dot(deltaR,deltaR)*0.5/Tau)<1e-6);
    }
    RealType eloc= H1[ipsi]->evaluate(W);
    std::cerr <<tempReptile[i]->Action(ipsi,Directionless)<<" "<<0.5*Tau*eloc<< std::endl;
    assert(std::abs(tempReptile[i]->Action(ipsi,Directionless)-0.5*Tau*eloc)<1e-6);
  }
  //    std::cerr <<"Done checking bead info"<< std::endl;
}

void RQMCMultiplePbyP::checkReptileProperties()
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
      *(NewBead->Laplacians[ipsi])=W.L;
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
    //cout << " WF : " << ipsi << " " << Reptile->GlobalAction[ipsi] << std::endl;
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
    //cout << " WF : " << ipsi << " " << Reptile->GlobalAction[ipsi] << std::endl;
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



void RQMCMultiplePbyP::recordBlock(int block)
{
  Reptile->record();
  //Write stuff
  //TEST CACHE
  //Estimators->report(CurrentStep);
  //TEST CACHE
  //app_error() << " BROKEN RQMCMultiplePbyP::recordBlock(int block) HDFWalkerOutput as 2007-04-16 " << std::endl;
  //HDFWalkerOutput WO(RootName,false,0);
  //WO.get(W);
  //WO.write(*branchEngine);
  //Reptile->write(WO.getFileID());
}

//   bool RQMCMultiplePbyP::put(xmlNodePtr q)
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
bool RQMCMultiplePbyP::put(xmlNodePtr q)
{
  std::cerr <<"PUT BEING CALLED"<< std::endl;
  std::string observ("NONE");
  OhmmsAttributeSet attrib;
  ParameterSet nattrib;
  attrib.add(observ,"observables" );
  nattrib.add(MSS,"mass","double" );
  nattrib.add(KEcut,"KEcut","double" );
  attrib.put(q);
  nattrib.put(q);
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
    if (observ=="NONE")
    {
      std::cout <<"Using normal Observables"<< std::endl;
      multiEstimator = new CSPolymerEstimator(H,nPsi);
    }
    else
      if (observ=="ZVZB")
      {
        std::cout <<"Using ZVZB observables"<< std::endl;
//         multiEstimator = new MJPolymerEstimator(H,nPsi);
        MJPolymerEstimator* MJp = new MJPolymerEstimator(H,nPsi);
        MJp->setpNorm(1.0/( W.Lattice.DIM *  W.Lattice.Volume));
        multiEstimator = MJp;
      }
      else
        if (observ=="HFDHE2")
        {
          std::cout <<"Using HFDHE2 observables"<< std::endl;
          HFDHE2PolymerEstimator* HFp = new HFDHE2PolymerEstimator(H,nPsi);
          HFp->setpNorm(1.0/( W.Lattice.DIM *  W.Lattice.Volume));
          multiEstimator = HFp;
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
//     std::cout <<"  Mass for Propagator is: "<<1.0/oneovermass<< std::endl;
//     std::cout <<"  m_over2t: "<<m_oneover2tau<< std::endl;
  return true;
}


void RQMCMultiplePbyP::moveReptile()
{
  myTimers[1]->start();
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
  //Compute the "diffusion-drift"\f$=(R_{yp}-R{xp})/tau\f$
  deltaR= NewBead->R - head->R;
  ParticleSet::ParticlePos_t newGrowth;
  ParticleSet::ParticlePos_t oldGrowth;
  newGrowth.resize(deltaR.size());
  oldGrowth.resize(deltaR.size());
  newGrowth=  gRand;
  //    oldGrowth= tail->R-next->R;
  //Compute HEAD action in forward direction
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    gRand = deltaR-Tau*(*head->Gradients[ipsi]);
    head->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(gRand,gRand);
  }
  //evaluate all relevant quantities in the new position
  int totbeadwgt(0);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    RealType* restrict NewBeadProp=NewBead->getPropertyBase(ipsi);
    //evaluate Psi and H
    NewBeadProp[LOGPSI]=Psi1[ipsi]->evaluateLog(W);
    NewBeadProp[SIGN]=Psi1[ipsi]->getPhase();
    RealType eloc=NewBeadProp[LOCALENERGY]= H1[ipsi]->evaluate(W);
    //Save properties
    H1[ipsi]->saveProperty(NewBeadProp);
    //*(NewBead->Gradients[ipsi])=W.G;
    Copy(W.G,*(NewBead->Gradients[ipsi]));
    ///hack      *(NewBead->Laplacians[ipsi])=W.L;
    //      Copy(W.L,*(NewBead->Laplacians[ipsi]));
    //Compute the backward part of the Kinetic action
    //gRand=deltaR+Tau*W.G;
    PAOps<RealType,DIM>::axpy(Tau,W.G,deltaR,gRand);
    NewBead->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(gRand,gRand);
    NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    int beadwgt=std::abs( ( Reptile->getSign(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
    NewBead->BeadSignWgt[ipsi]=beadwgt;
    totbeadwgt+=beadwgt;
  }
  //Compute Drift and TransProb here. This could be done (and was originally done)
  //after acceptance if # of beads is greater than one but needs to be done here
  //to make the one-bead case working ok.
  NewBead->getDrift(branchEngine->LogNorm);
  gRand=deltaR+Tau*NewBead->Drift;
  NewBead->TransProb[backward]=m_oneover2tau*Dot(gRand,gRand);
  oldGrowth= std::sqrt(m_oneover2tau*2.0)*gRand;
  //    std::cerr <<"old growth: "<<oldGrowth<<" "<< std::endl;
  //    std::cerr <<"new growth: "<<newGrowth<< std::endl;
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
    RealType dist2_new=0; //=Dot(newGrowth,newGrowth);
    RealType dist2_old=0; //=Dot(oldGrowth,oldGrowth);
    for (int i=0; i<newGrowth.size(); i++)
    {
      dist2_new+=std::log(dot(newGrowth[i],newGrowth[i]));
      dist2_old+=std::log(dot(oldGrowth[i],oldGrowth[i]));
    }
    AcceptProb=std::exp(NewGlobalWgt - Reptile->GlobalWgt + head->TransProb[forward] - next->TransProb[backward]); // - (dist2_new) + dist2_old);
  }
  //    std::cerr <<"Accept Prob is " <<AcceptProb<< std::endl;
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
    //HACK!
//       Walker_t::Buffer_t& w_buffer(NewBead->DataSet);
//       w_buffer.rewind();
//        W.copyToBuffer(w_buffer);
//        for (int ipsi=0;ipsi<nPsi;ipsi++)
// 	 Psi1[ipsi]->evaluate(W,w_buffer);
//       //PUT COMMENT BACK IN TO DEHACK
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
  myTimers[1]->stop();
}


void RQMCMultiplePbyP::ActionChange(bool oldData)
{
  RealType AddPsi=0.0;
  int toMult;
  if (oldData)
    toMult=1;
  else
    toMult=-1;
  for (int ipsi=0; ipsi<nPsi; ipsi++)
  {
    for (int i=0; i<tempReptile.size()-1; i++)
    {
      DeltaG[ipsi]+= tempReptile[i]->Action(ipsi,forward)*toMult;
      AddPsi+= tempReptile[i]->Action(ipsi,forward)*toMult;
      //      std::cerr <<"CHANGE: "<<i<<" "<<tempReptile[i]->Action(0,1)*toMult<< std::endl;
    }
    for (int i=1; i<tempReptile.size(); i++)
    {
      DeltaG[ipsi]+= tempReptile[i]->Action(ipsi,backward)*toMult;
      AddPsi+= tempReptile[i]->Action(ipsi,backward)*toMult;
      //      std::cerr <<"CHANGE2: "<<i<<" "<<tempReptile[i]->Action(ipsi,forward)*toMult<< std::endl;
    }
    for (int i=1; i<tempReptile.size()-1; i++)
    {
      DeltaG[ipsi]+= 2.0*tempReptile[i]->Action(ipsi,Directionless)*toMult;
    }
    for (int i=0; i<tempReptile.size(); i++)
    {
      DeltaSign[ipsi] -= tempReptile[i]->BeadSignWgt[ipsi]*toMult;
    }
  }
  //  std::cerr <<"AddPsi for ActionChange is "<<AddPsi<< std::endl;
  //  std::cerr <<"Delta G for ACtionchange is "<<DeltaG<< std::endl;
}



}
