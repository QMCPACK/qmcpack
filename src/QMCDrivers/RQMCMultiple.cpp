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
    WReptile.resize(n);
    logpsi.resize(n);
    gRand.resize(nptcl);
    NewTotalSign.resize(n);
    WeightSign.resize(n);
    RefSign.resize(n); RefSign=1;

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
	Reptile->TotalSign[ipsi] += int(BeadSign);

	// Save them in curW
	curW.Action(ipsi,MinusDirection)=0.25/Tau*KinActMinus;
	curW.Action(ipsi,PlusDirection)=0.25/Tau*KinActPlus;
	curW.Action(ipsi,Directionless)=0.5*Tau*eloc;
      }

      ++bead;
    }// End Loop over beads

    //Loop over Links
    bead=first_bead;
    //Initialize Reptile weights 
    for(int ipsi=0; ipsi<nPsi; ipsi++) WReptile[ipsi]=(*bead)->Properties(ipsi,LOGPSI);
    while(bead != last_bead){
      //Some the link action for the current link
      for(int ipsi=0; ipsi<nPsi; ipsi++){
	WReptile[ipsi]-=( (*bead)->Action(PlusDirection)+(*bead)->Action(Directionless)+
	    (*(bead+1))->Action(MinusDirection)+(*(bead+1))->Action(Directionless) );
      } 
      bead++;
    }
    //Finalize Reptile weights 
    for(int ipsi=0; ipsi<nPsi; ipsi++) WReptile[ipsi]+=(*bead)->Properties(ipsi,LOGPSI);

    //Weight: 1 if all beads have the same sign. 0 otherwise.
    for(int ipsi=0; ipsi<nPsi; ipsi++){ 
      //Compute the reference sign. This should be read if restarted.
      if(Reptile->TotalSign[ipsi]<0)RefSign[ipsi]=-1;
      WeightSign[ipsi]=std::max(0,RefSign[ipsi]*Reptile->TotalSign[ipsi]-Reptile->Last);
    }

    // Compute all kinds of Relative Weights 
    if(nPsi>1){
      RealType wdiff=WReptile[1]-WReptile[0];
      Reptile->LogRatioActionIJ=wdiff;
      Reptile->SumRatio = WeightSign[0]+WeightSign[1]*exp(wdiff);
      Reptile->UmbrellaWeight[0]=WeightSign[0]/Reptile->SumRatio;
      Reptile->UmbrellaWeight[1]=WeightSign[1]*(1.0-Reptile->UmbrellaWeight[0]);
    }else{
      Reptile->SumRatio = 1.0;
      Reptile->UmbrellaWeight[0]=1.0;
    }

    // Compute initial drift for each bead
    bead=first_bead;
    while(bead != bead_end) {
      Bead& curW(**bead);
      if(nPsi>1){
	RealType w0(std::max( 0.0,curW.Properties(0,SIGN) ) ),
		 w1(std::max( 0.0,curW.Properties(1,SIGN) ) );
	///Compute (Psi1)^2/(Psi0)^2 at the walker position
	RealType wgtsign=w0+w1*exp(2.0*(curW.Properties(1,LOGPSI)-curW.Properties(0,LOGPSI)));
	/// Convert to (Psi0)^2/(Psi0^2+Psi1^2)
	wgtsign=w0/wgtsign;
	/// Compute Drift as weighted average
	curW.Drift = wgtsign*(*curW.Gradients[0]) + (1.0-wgtsign)*(*curW.Gradients[1]);
      }else{
	/// Drift is just the gradient if only 1 WF
	curW.Drift = *curW.Gradients[0];
      }
      ++bead;
    }
  } // END OF InitReptile





  bool RQMCMultiple::run() { 

    Estimators->reportHeader();

    initReptile();

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

    //PolymerEstimator reptileReport(*Reptile,nPsi);
    //reptileReport.resetReportSettings(RootName);

    //accumulate configuration: probably need to reorder
    HDFWalkerOutput WO(RootName);

    RealType oneoversteps=1.0/static_cast<RealType>(nSteps);

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
      }
      if(nPsi>1) *OutEnergy << AveEloc[1]-AveEloc[0];
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
    Bead *head,*tail,*next;
    int ihead,inext,itail;
    //MultiChain::iterator anchor(Reptile->begin()),tail(Reptile->end()),next(tail-1);

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
    head = (*Reptile)[ihead];
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(gRand);

    //new position,\f$R_{yp} = R_{xp} + \sqrt(2}*\Delta + tau*D_{xp}\f$
    //W.R = head->R + sqrt(Tau)*gRand + Tau*head->Drift;
    W.R = head->R + m_sqrttau*gRand + Tau*head->Drift;
    //Save position in NewBead
    NewBead->R=W.R; 
    //update the distance table associated with W
    DistanceTable::update(W);

    //Compute the "diffusion-drift"\f$=(R_{yp}-R{xp})/tau\f$
    deltaR= NewBead->R - head->R;

    //Compute the probability for this move to happen, see DMC for the factors
    RealType LogForwardProb=-0.5*Dot(gRand,gRand);

    for(int ipsi=0; ipsi<nPsi; ipsi++) {
      gRand = deltaR-Tau*(*head->Gradients[ipsi]);
      head->Action(ipsi,forward)=0.5*m_oneover2tau*Dot(gRand,gRand);
    }

    //evaluate all relevant quantities in the new position
    for(int ipsi=0; ipsi<nPsi; ipsi++) {

      RealType* restrict NewBeadProp=NewBead->getPropertyBase(ipsi);

      //evaluate Psi and H
      logpsi[ipsi]=NewBeadProp[LOGPSI]=Psi1[ipsi]->evaluateLog(W);
      NewBeadProp[SIGN]=Psi1[ipsi]->getSign();
      RealType eloc=NewBeadProp[LOCALENERGY]= H1[ipsi]->evaluate(W);

      //Save properties
      H1[ipsi]->saveProperty(NewBeadProp);
      *(NewBead->Gradients[ipsi])=W.G; 

      //Compute the backward part of the Kinetic action
      gRand=deltaR+Tau*W.G;
      NewBead->Action(ipsi,backward)=0.5*m_oneover2tau*Dot(gRand,gRand);

      NewBead->Action(ipsi,Directionless)=0.5*Tau*eloc;
    }

    //Compute the Drift in new position (could be done after accepting)
    if(nPsi>1){
      RealType w0( std::max(0.0,NewBead->Properties(0,SIGN))),
	       w1( std::max(0.0,NewBead->Properties(1,SIGN)));
      ///Compute (Psi1)^2/(Psi0)^2 at the new position
      double wgtsign=w0+w1*exp(2.0*(logpsi[1]-logpsi[0]));
      /// Convert to (Psi0)^2/(Psi0^2+Psi1^2)
      wgtsign=w0/wgtsign;
      /// Compute Drift as weighted average
      NewBead->Drift = wgtsign*(*NewBead->Gradients[0])+(1.0-wgtsign)*(*NewBead->Gradients[1]);
    }else{
      /// Drift is Gradient if only 1 Psi
      NewBead->Drift=*NewBead->Gradients[0];
    }

    tail=(*Reptile)[itail]; next=(*Reptile)[inext];
    gRand = tail->R - next->R - Tau*next->Drift;
    RealType LogBackwardProb=-m_oneover2tau*Dot(gRand,gRand);

    // Evaluate the probability ratio between direct and inverse move
    RealType LogRatioTransProbability = LogBackwardProb-LogForwardProb;
    for(int ipsi=0; ipsi<nPsi; ipsi++) {

      RealType DeltaAction, DeltaEdge;
      // Initialize with the contribution of the new link
      DeltaAction = ( head->Action(ipsi,forward)+head->Action(ipsi,Directionless)+
	  NewBead->Action(ipsi,backward)+NewBead->Action(ipsi,Directionless) );
      // Subtract the contribution of the old link
      DeltaAction-= ( tail->Action(ipsi,forward)+tail->Action(ipsi,Directionless)+
	  next->Action(ipsi,backward)+next->Action(ipsi,Directionless) );
      /// Initialize with the new Wave funtion contribution at the edges
      DeltaEdge=NewBead->Properties(ipsi,LOGPSI)+next->Properties(ipsi,LOGPSI);
      /// Subtract the old wave funtion contribution at the edges
      DeltaEdge-= ( head->Properties(ipsi,LOGPSI)+tail->Properties(ipsi,LOGPSI) );
      //Global Weight
      WReptile[ipsi]= DeltaEdge-DeltaAction;
      //Compute the new sign
      NewTotalSign[ipsi]=Reptile->TotalSign[ipsi]+
	int(NewBead->Properties(ipsi,SIGN)-tail->Properties(ipsi,SIGN));
      //Weight: 1 if all beads have the same sign. 0 otherwise.
      WeightSign[ipsi]=std::max(0,RefSign[ipsi]*NewTotalSign[ipsi]-Reptile->Last);
    }


    //Move is completed. Compute acceptance probability. 
    if(nPsi>1){
      NewLogRatioAction=Reptile->LogRatioActionIJ+WReptile[1]-WReptile[0];
      NewSumRatio=WeightSign[0]+WeightSign[1]*exp(NewLogRatioAction);
      accept=exp(WReptile[0]+LogRatioTransProbability)*NewSumRatio/Reptile->SumRatio;
    }else{
      accept=WeightSign[0]*exp(WReptile[0]+LogRatioTransProbability);
    }

    //Update Reptile information
    if(Random() < accept){

      if(nPsi>1){
	Reptile->LogRatioActionIJ=NewLogRatioAction;
	Reptile->SumRatio=NewSumRatio;
	Reptile->UmbrellaWeight[0]=WeightSign[0]/Reptile->SumRatio;
	Reptile->UmbrellaWeight[1]=WeightSign[1]*(1.0-Reptile->UmbrellaWeight[0]);
	Reptile->TotalSign=NewTotalSign;
      }else{
	Reptile->SumRatio = 1.0;
	Reptile->UmbrellaWeight[0]=1.0;
      }

      //Add the NewBead to the Polymer.
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
