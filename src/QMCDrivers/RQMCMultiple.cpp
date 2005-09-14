//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/RQMCMultiple.h"
#include "QMCDrivers/MultiChain.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"
#include "OhmmsPETE/OhmmsMatrix.h"
namespace ohmmsqmc { 
  RQMCMultiple::RQMCMultiple(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), 
    ReptileLength(21),
    NumTurns(0), Reptile(0)
  { 
    RootName = "rmc-multi";
    QMCType ="rmc-multi";
    m_param.add(ReptileLength,"chains","int");
    QMCDriverMode.set(QMC_MULTIPLE,1);
    //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
    add_H_and_Psi(&h,&psi);
  }

  RQMCMultiple::~RQMCMultiple() {
    if(Reptile) delete Reptile;
    delete NewBead;
  }


  void RQMCMultiple::resizeArrays(int n) {

    nPsi = n;
    nptcl=W.G.size();
    gRand.resize(nptcl);
    NewGlobalAction.resize(n);
    NewGlobalSignWgt.resize(n);
    DeltaG.resize(n);
    WeightSign.resize(n);

    //Register properties for each walker
    for(int ipsi=0; ipsi<nPsi; ipsi++) H1[ipsi]->add2WalkerProperty(W);

    //resize Walker::Properties to hold everything
    W.resetWalkerProperty(nPsi);

    H1[0]->setPrimary(true);
    for(int ipsi=1; ipsi<nPsi; ipsi++) {
      H1[ipsi]->setPrimary(false);
    }
  }


  ///initialize Reptile
  void RQMCMultiple::initReptile() {

    //Resize working arrays
    resizeArrays(Psi1.size());

    //Temporary vector
    vector<int> SumSign;
    SumSign.resize(nPsi);

    //Initial direction of growth. To be read if calculation is restarted.
    int InitialGrowthDirection(0);

    //Reptile is made up by replicating the first walker. To be read if restarted.
    if(Reptile == 0) Reptile=new MultiChain(*W.begin(),ReptileLength,InitialGrowthDirection,nPsi);

    //Build NewBead. This takes care of a bunch of resizing operation.
    NewBead=new Bead(**W.begin());

    ///Assign a bunch of useful pointers
    MultiChain::iterator first_bead(Reptile->begin()), bead_end(Reptile->end());
    MultiChain::iterator bead(first_bead),last_bead(bead_end-1);

    ///Loop over beads to initialize action and WF
    while(bead != bead_end){

      ///Pointer to the current walker
      Bead& curW(**bead);

      ///Copy to W (ParticleSet) to compute distances, Psi and H
      W.R=curW.R;

      ///Compute Distances
      DistanceTable::update(W);

      ///loop over WF to compute contribution to the action and WF
      for(int ipsi=0; ipsi<nPsi; ipsi++) {

	///Compute Energy and Psi and save in curW
	curW.Properties(ipsi,LOGPSI) = Psi1[ipsi]->evaluateLog(W);
	RealType BeadSign = curW.Properties(ipsi,SIGN) = Psi1[ipsi]->getSign();
	RealType eloc= H1[ipsi]->evaluate(W);
	curW.Properties(ipsi,LOCALENERGY)= eloc;
	H1[ipsi]->saveProperty(curW.getPropertyBase(ipsi));
	*curW.Gradients[ipsi]=W.G;

	///Initialize Kinetic Action
	RealType KinActMinus=0.0;
	RealType KinActPlus=0.0;

	// Compute contribution to the Action in the MinusDirection
	if(bead!=first_bead){//forward action
	  Bead& prevW(**(bead-1));
	  deltaR=prevW.R-curW.R - Tau*W.G;
	  KinActMinus=Dot(deltaR,deltaR);
	}

	// Compute contribution to the Action in the PlusDirection
	if(bead!=last_bead){//backward action
	  Bead& nextW(**(bead+1));
	  deltaR=nextW.R-curW.R - Tau*W.G;
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

    //Assign Reference Sign as the majority Sign
    for(int ipsi=0; ipsi<nPsi; ipsi++) {
      if(SumSign[ipsi]>0)
       Reptile->RefSign[ipsi]=1;
      else Reptile->RefSign[ipsi]=-1;
    }

    //Compute Sign-weight for each bead
    bead=first_bead;
    while(bead != bead_end){
      Bead& curW(**bead);
      for(int ipsi=0; ipsi<nPsi; ipsi++) {
        int BeadSign = int(curW.Properties(ipsi,SIGN));
        curW.BeadSignWgt[ipsi]=abs((BeadSign+Reptile->RefSign[ipsi])/2);
      }
      ++bead;
    }

    // Compute initial drift for each bead
    bead=first_bead;
    while(bead != bead_end) {
      Bead& curW(**bead);
      curW.Drift=curW.BeadSignWgt[0]*(*curW.Gradients[0]);
      RealType denom=curW.BeadSignWgt[0];
      for(int ipsi=1; ipsi<nPsi; ipsi++) {
        RealType wgtpsi=curW.BeadSignWgt[ipsi]*
          exp(2.0*(curW.Properties(ipsi,LOGPSI)-curW.Properties(0,LOGPSI)));
        curW.Drift += (wgtpsi*(*curW.Gradients[ipsi]));
        denom += wgtpsi;
      }
      denom=1.0/denom;
      curW.Drift = curW.Drift*denom;
      ++bead;
    }

    //Compute initial transition Probability within the chain
    bead=first_bead;
    while(bead != bead_end) {
      Bead& curW(**bead);
      RealType TrProbMinus=0.0;
      RealType TrProbPlus=0.0;
      // Compute contribution to the Transition Prob in the MinusDirection
      if(bead!=first_bead){//forward action
        Bead& prevW(**(bead-1));
        deltaR=prevW.R-curW.R - Tau*curW.Drift;
        TrProbMinus=Dot(deltaR,deltaR);
      }
      
      // Compute contribution to the Transition Prob in the PlusDirection
      if(bead!=last_bead){//backward action
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
    for(int ipsi=0; ipsi<nPsi; ipsi++) Reptile->GlobalAction[ipsi]=(*bead)->Properties(ipsi,LOGPSI);
    while(bead != last_bead){
      for(int ipsi=0; ipsi<nPsi; ipsi++){
	Reptile->GlobalAction[ipsi]-=( (*bead)->Action(PlusDirection) + (*(bead+1))->Action(MinusDirection)+
                                       (*bead)->Action(Directionless) + (*(bead+1))->Action(Directionless)   );
      } 
      bead++;
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++) Reptile->GlobalAction[ipsi]+=(*bead)->Properties(ipsi,LOGPSI);

    //Compute Global Sign weight (need to be initialized somewhere)
    bead=first_bead;
    while(bead != bead_end){
      for(int ipsi=0; ipsi<nPsi; ipsi++){ 
        Reptile->GlobalSignWgt[ipsi] += (*bead)->BeadSignWgt[ipsi];
      }
      ++bead;
    }

    //Compute reference action
    RealType RefAction(-10.0e20);
    for(int ipsi=0; ipsi<nPsi; ipsi++){
      WeightSign[ipsi]=std::max(0,Reptile->GlobalSignWgt[ipsi]-Reptile->Last);
      if(WeightSign[ipsi])RefAction=max(RefAction,Reptile->GlobalAction[ipsi]);
    }

    //Compute Total Weight
    Reptile->GlobalWgt=0.0e0;
    for(int ipsi=0; ipsi<nPsi; ipsi++){
      RealType DeltaAction(Reptile->GlobalAction[ipsi]-RefAction);
      if((WeightSign[ipsi]>0) && (DeltaAction > -30)) Reptile->GlobalWgt += exp(DeltaAction);
    }
    Reptile->GlobalWgt=log(Reptile->GlobalWgt)+RefAction;

    //Compute Umbrella Weight 
    for(int ipsi=0; ipsi<nPsi; ipsi++){
      RealType DeltaAction(Reptile->GlobalAction[ipsi]-Reptile->GlobalWgt);
      if((WeightSign[ipsi]>0) && (DeltaAction > -30)) Reptile->UmbrellaWeight[ipsi] = exp(DeltaAction);
      else Reptile->UmbrellaWeight[ipsi] = 0.0e0;
    }
 
  } // END OF InitReptile





  bool RQMCMultiple::run() { 

    Estimators->reportHeader();

    cout << "Initializing Reptile ... " ;
    initReptile();
    cout << "Reptile has been initialized" << endl;

    IndexType block = 0;
    Pooma::Clock timer;
    IndexType accstep=0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;

    int ipsi(0);

    std::vector<double>AveEloc,AveWeight;
    AveEloc.resize(nPsi);
    AveWeight.resize(nPsi);
    std::string filename(RootName);
    filename=RootName+".Eloc.dat";
    ofstream *OutEnergy;
    OutEnergy=new ofstream(filename.c_str());
    cout << "Energies written on " << filename << endl;

    //PolymerEstimator reptileReport(*Reptile,nPsi);
    //reptileReport.resetReportSettings(RootName);

    //accumulate configuration: probably need to reorder
    HDFWalkerOutput WO(RootName);

    RealType oneoversteps=1.0/static_cast<RealType>(nSteps);

    cout << "Start Looping over Blocks " << endl;
    do { //Loop over Blocks

      IndexType step = 0;
      timer.start();
      NumTurns = 0;

      for(int ipsi=0; ipsi<nPsi; ipsi++){
	AveEloc[ipsi]=0.0;
	AveWeight[ipsi]=0.0;
      }

      do { //Loop over steps

	moveReptile();
	step++; accstep++;

	//Copy the front and back to W to take average report
	//W.copyWalkerRefs(Reptile->front(),Reptile->back());
	for(int ipsi=0; ipsi<nPsi; ipsi++){
	  double WeightedEloc=Reptile->UmbrellaWeight[ipsi]*
	    ( Reptile->front()->Action(ipsi,Directionless)
	      +Reptile->back()->Action(ipsi,Directionless) );
	  AveEloc[ipsi]+=WeightedEloc;
	  AveWeight[ipsi]+=Reptile->UmbrellaWeight[ipsi];
	}
	Estimators->accumulate(W);

	//reptileReport.accumulate();

      } while(step<nSteps);

      *OutEnergy << block << " " ;
      for(int ipsi=0; ipsi<nPsi; ipsi++){
	AveEloc[ipsi]/=(AveWeight[ipsi]*Tau);
	*OutEnergy << AveEloc[ipsi] << " ";
	*OutEnergy << AveWeight[ipsi]/nSteps << " ";
      }
      *OutEnergy << endl;
      OutEnergy->flush();


      nAcceptTot += nAccept;
      nRejectTot += nReject;

      RealType acceptedR = static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject); 
      Estimators->flush();
      Estimators->setColumn(AcceptIndex,acceptedR);
      Estimators->report(accstep);

      //reptileReport.report(accstep);

      //change NumCuts to make accstep ~ 50%
      LogOut->getStream() 
	<< "Block " << block << " " 
	<< timer.cpu_time() << " " << NumTurns << endl;

      nAccept = 0; nReject = 0;
      block++;

      if(pStride) WO.get(W);

    } while(block<nBlocks);

    delete OutEnergy;

    LogOut->getStream() 
      << "ratio = " 
      << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << endl;

    Estimators->finalize();
    return true;
  }





  bool RQMCMultiple::put(xmlNodePtr q){
    //nothing to do yet
    MinusDirection=0;
    PlusDirection=1;
    Directionless=2;
    return true;
  }






  void RQMCMultiple::moveReptile(){

    //Used several times
    m_oneover2tau=0.5/Tau; 
    m_sqrttau=sqrt(Tau);
    int ihead,inext,itail;

    //Depending on the growth direction initialize growth variables
    if(Reptile->GrowthDirection==MinusDirection){
      forward=MinusDirection;
      backward=PlusDirection;
      ihead = 0; 
      itail = Reptile->Last;
      inext = itail-1;
    }else{
      forward=PlusDirection;
      backward=MinusDirection;
      ihead = Reptile->Last;
      itail = 0;
      inext = 1;
    }

    //Point head to the growing end of the reptile
    Bead *head,*tail,*next;
    head = (*Reptile)[ihead];
    tail=(*Reptile)[itail]; next=(*Reptile)[inext];
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(gRand);

    //new position,\f$R_{yp} = R_{xp} + \sqrt(2}*\Delta + tau*D_{xp}\f$
    W.R = head->R + m_sqrttau*gRand + Tau*head->Drift;
    //Save Transition Probability
    head->TransProb[forward]=0.5*Dot(gRand,gRand);
    //Save position in NewBead
    NewBead->R=W.R; 
    //update the distance table associated with W
    DistanceTable::update(W);

    //Compute the "diffusion-drift"\f$=(R_{yp}-R{xp})/tau\f$
    deltaR= NewBead->R - head->R;

    //Compute HEAD action in forward direction
    for(int ipsi=0; ipsi<nPsi; ipsi++) {
      gRand = deltaR-Tau*(*head->Gradients[ipsi]);
      head->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(gRand,gRand);
    }

    //evaluate all relevant quantities in the new position
    int totbeadwgt(0);
    for(int ipsi=0; ipsi<nPsi; ipsi++) {

      RealType* restrict NewBeadProp=NewBead->getPropertyBase(ipsi);

      //evaluate Psi and H
      NewBeadProp[LOGPSI]=Psi1[ipsi]->evaluateLog(W);
      NewBeadProp[SIGN]=Psi1[ipsi]->getSign();
      RealType eloc=NewBeadProp[LOCALENERGY]= H1[ipsi]->evaluate(W);

      //Save properties
      H1[ipsi]->saveProperty(NewBeadProp);
      *(NewBead->Gradients[ipsi])=W.G; 

      //Compute the backward part of the Kinetic action
      gRand=deltaR+Tau*W.G;
      NewBead->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(gRand,gRand);

      NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
      int beadwgt=abs( ( int(NewBeadProp[SIGN])+Reptile->RefSign[ipsi] )/2 );
      NewBead->BeadSignWgt[ipsi]=beadwgt;
      totbeadwgt+=beadwgt;
    }

    RealType AcceptProb(-1.0),NewGlobalWgt(0.0);
    RealType RefAction(-10.0e20);
    if(totbeadwgt!=0){
      for(int ipsi=0; ipsi<nPsi; ipsi++) {

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
        if(WeightSign[ipsi]>0)RefAction=max(RefAction,NewGlobalAction[ipsi]);
      }
  
      //Compute Log of global Wgt
      for(int ipsi=0; ipsi<nPsi; ipsi++) {
        RealType DeltaAction(NewGlobalAction[ipsi]-RefAction);
        if((WeightSign[ipsi]>0) && (DeltaAction > -30.0)) NewGlobalWgt+=exp(DeltaAction);
      }
      NewGlobalWgt=log(NewGlobalWgt)+RefAction;

      AcceptProb=exp(NewGlobalWgt - Reptile->GlobalWgt + head->TransProb[forward] - next->TransProb[backward]);
    }

    if(Random() < AcceptProb){

      //Update Reptile information
      Reptile->GlobalWgt=NewGlobalWgt;
      for(int ipsi=0; ipsi<nPsi; ipsi++) {
	Reptile->GlobalAction[ipsi]=NewGlobalAction[ipsi];
        Reptile->GlobalSignWgt[ipsi]=NewGlobalSignWgt[ipsi];
        RealType DeltaAction(NewGlobalAction[ipsi]-NewGlobalWgt);
        if((WeightSign[ipsi]>0) && (DeltaAction > -30.0))
	  Reptile->UmbrellaWeight[ipsi]=exp(DeltaAction);
        else Reptile->UmbrellaWeight[ipsi]=0.0e0;
      }

      //Compute Drift and TransProb in NewBead position
      NewBead->Drift=NewBead->BeadSignWgt[0]*(*NewBead->Gradients[0]);
      RealType denom=NewBead->BeadSignWgt[0];
      for(int ipsi=1; ipsi<nPsi; ipsi++) {
        RealType wgtpsi=NewBead->BeadSignWgt[ipsi]*
	  exp(2.0*(NewBead->Properties(ipsi,LOGPSI)-NewBead->Properties(0,LOGPSI)));
        NewBead->Drift += (wgtpsi*(*NewBead->Gradients[ipsi]));
        denom += wgtpsi;
      }
      denom=1.0/denom;
      NewBead->Drift = denom*NewBead->Drift;
      gRand=deltaR+Tau*NewBead->Drift;
      NewBead->TransProb[backward]=m_oneover2tau*Dot(gRand,gRand);

      //Add NewBead to the Polymer.
      if(Reptile->GrowthDirection==MinusDirection){
	Reptile->push_front(NewBead);
	NewBead=tail;
	Reptile->pop_back();
      }else{
	Reptile->push_back(NewBead);
	NewBead=tail;
	Reptile->pop_front();
      }
      ++nAccept;
    } else {
      ++nReject; 
      ++NumTurns; 
      Reptile->flip();
    }
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
