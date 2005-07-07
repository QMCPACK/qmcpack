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
#include "QMCDrivers/PolymerChain.h"
#include "QMCDrivers/PolymerEstimator.h"
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
    UseBounce(false),
    ReptileLength(21),
    NumCuts(1),
    NumTurns(0), Reptile(0)
  { 
    RootName = "rmc-multi";
    QMCType ="rmc-multi";
    m_param.add(ReptileLength,"chains","int");
    m_param.add(NumCuts,"cuts","int");
    m_param.add(UseBounce,"bounce","int");
    //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
    add_H_and_Psi(&h,&psi);
  }

  RQMCMultiple::~RQMCMultiple() {
    if(Reptile) delete Reptile;
  }


  void RQMCMultiple::resizeArrays(int n) {

    nPsi = n;
    int npair=n*(n-1)/2;
    TotalSign.resize(n);
    SignWeight.resize(n);
    beadSignWeight.resize(n);
    sumratio.resize(n);
    SumRatioAction.resize(n);
    WReptile.resize(n);
    logpsi.resize(n);
    LogRatioActionIJ.resize(npair);

    DiffusionDrift.resize(W.getTotalNum());
    LocalDrift.resize(W.getTotalNum());

    //Register properties for each walker
    for(int ipsi=0; ipsi<nPsi; ipsi++) H1[ipsi]->add2WalkerProperty(W);

    //Add forward/backward KineticAction ParticleSet::PropertyList 
    MinusKineticAction=W.addProperty("MinusKineticAction");
    PlusKineticAction=W.addProperty("PlusKineticAction");

    //resize Walker::Properties to hold everything
    W.resetWalkerProperty(nPsi);
  }

  ///initialize Reptile
  void RQMCMultiple::initReptile() {

    resizeArrays(Psi1.size());

    //overwrite the number of cuts for Bounce algorithm
    if(UseBounce) NumCuts = 1;
    LOGMSG("Moving " << NumCuts << " for each reptation step");
    if(Reptile == 0) {
      //Reptile is made up by replicating the first walker
      Reptile=new PolymerChain(*W.begin(),ReptileLength,NumCuts);
      //resize internal arrays to handle multiple H/Psi pairs 
      Reptile->resizeArrays(nPsi);
    }

    PolymerChain::iterator first_bead(Reptile->begin()), bead_end(Reptile->end());
    PolymerChain::iterator bead(first_bead),last_bead(bead_end-1);
    
    RealType oneovertau=1.0/Tau;
    WReptile=0.0;
    Reptile->TotalSign=0;

    ///Initialize each bead of the Reptile
    //This is general: it works also if beads are not on top of each other
    while(bead != bead_end){

      Walker_t& curW(**bead);

      W.R=curW.R;
      DistanceTable::update(W);

      for(int ipsi=0; ipsi<nPsi; ipsi++) {

        curW.Properties(ipsi,LOGPSI) = Psi1[ipsi]->evaluateLog(W);
        curW.Properties(ipsi,SIGN) = Psi1[ipsi]->getSign();
        RealType eloc= H1[ipsi]->evaluate(W);
        curW.Properties(ipsi,LOCALENERGY)= eloc;

        //update the Walker's properties
        H1[ipsi]->saveProperty(curW.getPropertyBase(ipsi));

          // Sign initialization
        Reptile->TotalSign[ipsi]+=static_cast<int>(Psi1[ipsi]->getSign());

        //ID is used to assign Gradients
        *Reptile->Gradients(curW.ID,ipsi)=W.G;

        RealType KinActMinus=0.0;
        RealType KinActPlus=0.0;

        // initialize Action info for each bead
        if(bead!=first_bead){//forward action
          Walker_t& prevW(**(bead-1));
          deltaR=oneovertau*(prevW.R-curW.R) - W.G;
          KinActMinus=Dot(deltaR,deltaR);
        }

        if(bead!=last_bead){//backward action
          Walker_t& nextW(**(bead+1));
          deltaR=oneovertau*(nextW.R-curW.R) - W.G;
          KinActPlus=Dot(deltaR,deltaR);
        } 

        curW.Properties(ipsi,MinusKineticAction)=KinActMinus;
        curW.Properties(ipsi,PlusKineticAction)=KinActPlus;

        // Construct Weight of the reptile for WF ipsi
        WReptile[ipsi] -= (Tau*(eloc+0.25*(KinActMinus+KinActPlus)));
      }

      ++bead;
    }// End Loop over beads

    //Check the sign-related properties with the first and last beads
    //Add the WF contribution to the action and correct the end beads
    for(int ipsi=0; ipsi<nPsi; ipsi++) {
      RealType* firstProp=(*first_bead)->getPropertyBase(ipsi);
      RealType* lastProp=(*last_bead)->getPropertyBase(ipsi);
      int totsign=Reptile->TotalSign[ipsi];
      //WReptile can be problematic: LOGPSI = log(abs(psi)) if first and last have different signs
      WReptile[ipsi]+=(firstProp[LOGPSI]+lastProp[LOGPSI]); 
      WReptile[ipsi]+=(0.5*Tau*(firstProp[LOCALENERGY]+lastProp[LOCALENERGY]));

      SignWeight[ipsi]=abs(totsign/ReptileLength);
      SumRatioAction[ipsi] = SignWeight[ipsi];

      //Compute reference sign for the polymer the same thing
      //Reptile->RefSigb[ipsi]=totsign/abs(totsign)
      Reptile->RefSign[ipsi] = (totsign<0)?-1:1;
    }

    int ij(0);
    for(int ipsi=0; ipsi<nPsi-1; ipsi++) {
      for(int jpsi=ipsi+1; jpsi<nPsi; jpsi++,ij++) {
        RealType wdiff=WReptile[jpsi]-WReptile[ipsi];
        Reptile->LogRatioActionIJ[ij]=wdiff;
        RealType r=exp(wdiff);
        SumRatioAction[ipsi]+=SignWeight[jpsi]*r;
        SumRatioAction[jpsi]+=SignWeight[ipsi]/r;
      }
    }

    Reptile->SumRatioPrimary=SumRatioAction[0];

    // Compute initial Umbrella Weight
    for(int ipsi=0; ipsi<nPsi-1; ipsi++) {
      Reptile->UmbrellaWeight[ipsi]=SignWeight[ipsi]/SumRatioAction[ipsi];
    }

    // Compute initial drift for each bead
    bead=first_bead;
    while(bead != bead_end) {
      Walker_t& curW(**bead);
      int totBeadWgt=0;
      for(int ipsi=0; ipsi<nPsi; ipsi++){
        totBeadWgt+=beadSignWeight[ipsi]=abs(curW.Properties(ipsi,SIGN)+Reptile->RefSign[ipsi])/2;
        sumratio[ipsi]=beadSignWeight[ipsi];
      }

      if(totBeadWgt==0) {
        ERRORMSG("The initial total bead weight is zero. Abort")
        OHMMS::Controller->abort();
      }

      for(int ipsi=0; ipsi<nPsi-1; ipsi++) {
        for(int jpsi=ipsi+1; jpsi<nPsi; jpsi++,ij++) {
          RealType r=exp(2.0*(curW.Properties(ipsi,LOGPSI)-curW.Properties(jpsi,LOGPSI)));
          sumratio[ipsi]+=beadSignWeight[jpsi]*r;
          sumratio[jpsi]+=beadSignWeight[ipsi]/r;
        }
      }
       
      int curID=curW.ID;
      RealType ratioM = beadSignWeight[0]/sumratio[0];
      LocalDrift = ratioM*(*Reptile->Gradients(curID,0));
      for(int ipsi=1; ipsi<nPsi; ipsi++) {
        ratioM = beadSignWeight[ipsi]/sumratio[ipsi];
        LocalDrift += ratioM*(*Reptile->Gradients(curID,ipsi));
      }
      curW.Drift= LocalDrift;
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

    PolymerEstimator reptileReport(*Reptile,nPsi);
    reptileReport.resetReportSettings(RootName);

    //accumulate configuration: probably need to reorder
    HDFWalkerOutput WO(RootName);

    RealType oneoversteps=1.0/static_cast<RealType>(nSteps);

    do { //Loop over Blocks

      IndexType step = 0;
      timer.start();
      NumTurns = 0;

      do { //Loop over steps

	moveReptile();
	step++; accstep++;

        //Copy the front and back to W to take average report
        W.copyWalkerRefs(Reptile->front(),Reptile->back());
	Estimators->accumulate(W);

	reptileReport.accumulate();

      } while(step<nSteps);

      nAcceptTot += nAccept;
      nRejectTot += nReject;

      RealType acceptedR = static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject); 
      Estimators->flush();
      Estimators->setColumn(AcceptIndex,acceptedR);
      Estimators->report(accstep);

      reptileReport.report(accstep);

      //change NumCuts to make accstep ~ 50%
      LogOut->getStream() 
	<< "Block " << block << " " 
	<< timer.cpu_time() << " " << NumTurns << " " << Reptile->getID() << endl;

      nAccept = 0; nReject = 0;
      block++;

      if(pStride) WO.get(W);

    } while(block<nBlocks);

    LogOut->getStream() 
      << "ratio = " 
      << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << endl;

    Estimators->finalize();
    return true;
  }


  bool RQMCMultiple::put(xmlNodePtr q){
    //nothing to do yet
    return true;
  }

  void RQMCMultiple::moveReptile(){

    RealType g = sqrt(Tau);
    RealType oneovertau=1.0/Tau;
    RealType oneover2tau=0.5*oneovertau;
    RealType tauover2 = 0.5*Tau;

    typedef MCWalkerConfiguration::Walker_t Walker_t;

    if(!UseBounce && Random()<0.5) {
      Reptile->flip(); 	  
      NumTurns++;
    }

    //temporary array for the move are filled.
    Walker_t* anchor = Reptile->makeEnds();

    //Transition is between x-y--X--xp to y--X--xp--yp:
    //tail-next---X--anchor to next--X--anchor-head
    //X is itself a Reptile made up of many beads.
    //The link x-y is deleted and xp-yp is created.
    //save the local energies of the anchor and tails
    //tails store the piece which is cut
    //heads the piece to be added
    //int forward=Reptile->MoveHead(int(!Reptile->MoveHead));
    //int backward(int(Reptile->MoveHead));
    int forward=MinusKineticAction;
    int backward=PlusKineticAction;
    if(!Reptile->MoveHead) {//moving backward, swap forward and backward
      forward=PlusKineticAction;
      backward=MinusKineticAction;
    } 

    NumCuts = Reptile->NumCuts;
    WReptile=0.0;
    TotalSign=Reptile->TotalSign;

    RealType LogRatioTransProbability = 0.0;

    //reference to Gradients of the moving parts
    Matrix<ParticleSet::ParticleGradient_t*>& headG(Reptile->HeadGradients);

    for(int i=0, iplus=1; i<NumCuts; i++,iplus++) {

      Walker_t* head=Reptile->heads[i]; //head walker (moving)

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);

      //new position,\f$R_{yp} = R_{xp} + \sqrt(2}*\Delta + tau*D_{xp}\f$
      W.R = anchor->R + g*deltaR + Tau*anchor->Drift;

      //Compute the "diffusion-drift"\f$=(R_{yp}-R{xp})/tau\f$
      DiffusionDrift= oneovertau*(W.R-anchor->R);

      //Compute the probability for this move to happen, see DMC for the factors
      RealType LogForwardProb=-0.5*Dot(deltaR,deltaR);

      for(int ipsi=0; ipsi<nPsi; ipsi++) {
        //\f$d\{\bf R\} = V_{d} - G_{0}\f$
        deltaR = DiffusionDrift-(*(headG(i,ipsi)));
        anchor->Properties(ipsi,forward)=Dot(deltaR,deltaR);
      }

      //update the distance table associated with W
      DistanceTable::update(W);
      int TotalBeadWgt=0;

      //evaluate Psi and H and update the weights 
      for(int ipsi=0; ipsi<nPsi; ipsi++) {

        RealType* restrict headProp=head->getPropertyBase(ipsi);//yp

        //evaluate Psi and H and save the properties
        logpsi[ipsi]=headProp[LOGPSI]=Psi1[ipsi]->evaluateLog(W);
        headProp[SIGN]=Psi1[ipsi]->getSign();
        headProp[LOCALENERGY]= H1[ipsi]->evaluate(W);
        H1[ipsi]->saveProperty(headProp);

        //Save gradient for the next 
        *headG(iplus,ipsi)=W.G; 

        //Compute the backward part of the Kinetic action
        deltaR=DiffusionDrift+W.G;

        //KinActBackward_yp[ipsi]=Dot(deltaR,deltaR);
        headProp[backward]=Dot(deltaR,deltaR);

        //beadSignWeight = 0 or 1, nothing else
        beadSignWeight[ipsi]=abs((headProp[SIGN]+Reptile->RefSign[ipsi])/2);
        sumratio[ipsi] = beadSignWeight[ipsi];
        TotalBeadWgt += beadSignWeight[ipsi];
      }

      //if(TotalBeadWgt == 0) all the Psi's change their signs
      if(TotalBeadWgt == 0) { //sign changed
	TotalSign=0; break;
      }

      for(int ipsi=0;ipsi<nPsi-1; ipsi++) {
        for(int jpsi=ipsi+1; jpsi<nPsi; jpsi++) {
          RealType ratioij=exp(2.0*logpsi(jpsi)-logpsi(ipsi));
          sumratio[ipsi]+=beadSignWeight[jpsi]*ratioij;
          sumratio[jpsi]+=beadSignWeight[ipsi]/ratioij;
        }
      }

      // Evaluate the drift in the new position
      //LocalDrift = (beadSignWeight[0]/sumratio[0])*head->Gradient[0];
      RealType r(beadSignWeight[0]/sumratio[0]);
      LocalDrift = r*(*headG(iplus,0));
      for(int ipsi=1;ipsi<nPsi; ipsi++) {
        r=beadSignWeight[ipsi]/sumratio[ipsi];
        LocalDrift += r*(*headG(iplus,ipsi));
      }
      head->R = W.R;
      head->Drift = LocalDrift;

      Walker_t* tail=Reptile->tails[i];  //tail walker 
      Walker_t* next=Reptile->tails[i+1];//inner walker to tail

      // \f$ R_{x}-R_{y}-tau*D_{y}\f$
      deltaR = (tail->R-next->R)-Tau*next->Drift;
      RealType LogBackwardProb=-oneover2tau*Dot(deltaR,deltaR);

      // Evaluate the probability ratio between direct and inverse move
      LogRatioTransProbability += (LogBackwardProb-LogForwardProb);
      for(int ipsi=0; ipsi<nPsi; ipsi++) {
        const RealType* restrict xp=anchor->getPropertyBase(ipsi);
        const RealType* restrict yp=head->getPropertyBase(ipsi);
        const RealType* restrict x =tail->getPropertyBase(ipsi);
        const RealType* restrict y =next->getPropertyBase(ipsi);
        //WReptile[I] is the ratio of the weight of I after
        //and before the move	
        //WReptile[ipsi] += 
        //  (-tauover2*(eloc_yp[ipsi]+eloc_xp[ipsi]-eloc_x[ipsi]-eloc_y[ipsi]
        //  	    +0.5*(KinActBackward_yp[ipsi]+KinActForward_xp[ipsi]
        //  	      -KinActForward_x[ipsi]-KinActBackward_y[ipsi]) ));
        //WReptile[ipsi] += (logpsi_yp[ipsi]-logpsi_xp[ipsi]-
        //    logpsi_x[ipsi]+logpsi_y[ipsi]);
        //TotalSign[ipsi] += (Sign_yp[ipsi]-Sign_x[ipsi]);
        WReptile[ipsi] 
          +=(-tauover2*(yp[LOCALENERGY]+xp[LOCALENERGY]-x[LOCALENERGY]-y[LOCALENERGY] 
                +0.5*(yp[backward]+xp[forward] -x[forward]-y[backward])));
        WReptile[ipsi] += (yp[LOGPSI]-xp[LOGPSI]-x[LOGPSI]+y[LOGPSI]);
        TotalSign[ipsi] += (yp[SIGN]-x[SIGN]);
      }

      //move the anchor and swap the local energies for WReptile
      anchor=head;
    }

    //Move is completed. Compute acceptance probability int/int is not a good idea
    for(int ipsi=0; ipsi<nPsi; ipsi++) {
      SignWeight[ipsi]=abs(TotalSign[ipsi]/ReptileLength);
      SumRatioAction[ipsi] = SignWeight[ipsi];
    }

    int ij(0);
    for(int ipsi=0; ipsi<nPsi-1; ipsi++) {
      for(int jpsi=ipsi+1; jpsi<nPsi; jpsi++,ij++) {
        RealType newaction=Reptile->LogRatioActionIJ[ij]+WReptile[jpsi]-WReptile[ipsi];
        LogRatioActionIJ[ij]=newaction;
        RealType r=exp(newaction);
        SumRatioAction[ipsi]+=SignWeight[jpsi]*r;
        SumRatioAction[jpsi]+=SignWeight[ipsi]/r;
      }
    }

    double accept = std::min( 1.0,SumRatioAction[0]/Reptile->SumRatioPrimary 
    	*exp(WReptile[0]+LogRatioTransProbability) );

    //Update Reptile information
    if(Random() < accept){
      Reptile->updateEnds();
      Reptile->LogRatioActionIJ=LogRatioActionIJ;
      Reptile->SumRatioPrimary=SumRatioAction[0];
      Reptile->TotalSign=TotalSign;
      for(int ipsi=0; ipsi<nPsi; ipsi++) 
        Reptile->UmbrellaWeight[ipsi]=SignWeight[ipsi]/SumRatioAction[ipsi];
      ++nAccept;
    } else {
      ++nReject; 
      if(UseBounce) {
	++NumTurns; Reptile->flip();
      }
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
