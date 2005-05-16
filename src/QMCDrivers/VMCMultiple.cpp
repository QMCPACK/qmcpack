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
#include "QMCDrivers/VMCMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

  /// Constructor.
  VMCMultiple::VMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h) { 
    RootName = "vmc";
    QMCType ="vmc";

    //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
    add_H_and_Psi(&h,&psi);
  }

  /** allocate internal data here before run() is called
   * @author SIMONE
   *
   * See QMCDriver::process
   */
  bool VMCMultiple::put(xmlNodePtr q){

    nPsi=Psi1.size();	

    if(IndexPsi.size() != nPsi) {
      IndexPsi.resize(nPsi,-1);//set it to -1 so that we add index only once
    }

    logpsi.resize(nPsi);
    dgrad.resize(nPsi);		
    lap.resize(nPsi);		
    sumratio.resize(nPsi);	
    invsumratio.resize(nPsi);	
    totinvsumratio.resize(nPsi);
    for(int ipsi=0; ipsi< Psi1.size(); ipsi++){
      char cname[6];
      int jpsi=sprintf(cname,"WPsi%d",ipsi);
      if(IndexPsi[ipsi]<0) IndexPsi[ipsi]=Estimators->addColumn(cname);
    }

    return true;
  }
  
  /** Run the VMCMultiple algorithm.
   *
   * Similar to VMC::run 
   */
  bool VMCMultiple::run() { 
    
    Estimators->reportHeader();

    //probably unnecessary
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    while(it != it_end) {
      (*it)->Properties(WEIGHT) = 1.0;
      ++it;
    }
    
    Estimators->reset();

    //create an output engine
    HDFWalkerOutput WO(RootName);
    
    IndexType block = 0;
    
    Pooma::Clock timer;
    
    double wh=0.0;
    IndexType accstep=0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    do {
      IndexType step = 0;
      timer.start();
      nAccept = 0; nReject=0;
      for(int ipsi=0; ipsi< nPsi; ipsi++)			//SIMONE
        totinvsumratio[ipsi] = 0.0;					//SIMONE
      do {
        advanceWalkerByWalker();
        step++;accstep++;
        Estimators->accumulate(W);
      } while(step<nSteps);
      
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      RealType BlockTotal=static_cast<RealType>(nAccept+nReject);				//SIMONE
      Estimators->flush();
      Estimators->setColumn(AcceptIndex,static_cast<RealType>(nAccept)/BlockTotal);	//SIMONE
      for(int ipsi=0; ipsi< nPsi; ipsi++)						//SIMONE
        Estimators->setColumn(IndexPsi[ipsi],totinvsumratio[ipsi]/BlockTotal);		//SIMONE
      Estimators->report(accstep);
      
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << endl;
      if(pStride) WO.get(W);
      nAccept = 0; nReject = 0; //SIMONE - Necessary?
      block++;
    } while(block<nBlocks);
    
    LogOut->getStream() 
      << "Ratio = " 
      << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << endl;
    
    if(!pStride) WO.get(W);
    
    Estimators->finalize();

    return true;
  }

  /**  Advance all the walkers one timstep. 
   */
  void 
  VMCMultiple::advanceWalkerByWalker() {
    
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    MCWalkerConfiguration::PropertyContainer_t Properties;
    
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 

    while(it != it_end) {
      
      //copy the properties of the working walker
      Properties = (*it)->Properties;           //SIMONE - Necessary?
      //save old local energy
      ValueType eold = Properties(LOCALENERGY); //SIMONE - Necessary?
      
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      
      W.R = g*deltaR + (*it)->R + (*it)->Drift;
      
      //update the distance table associated with W
      DistanceTable::update(W);
      
      //evaluate wave function
      //ValueType psi = Psi.evaluate(W);
      //Properties(LOGPSI) =log(fabs(psi));
      //Properties(PSI) = psi;
      //update the properties: note that we are getting 
      //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
      for(int ipsi=0; ipsi< nPsi;ipsi++) {			   //SIMONE
	logpsi[ipsi]=(Psi1[ipsi]->evaluateLog(W)); 		   //SIMONE
	lap[ipsi]=W.L;
	dgrad[ipsi]=W.G;					   //SIMONE
	sumratio[ipsi]=1.0;					   //SIMONE
      } 							   //SIMONE
      // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i	   //SIMONE 
      for(int ipsi=0; ipsi< nPsi; ipsi++) {			   //SIMONE
	for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++){     		   //SIMONE
	  RealType ratioij=exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));   //SIMONE
	  sumratio[ipsi] += ratioij;                               //SIMONE
	  sumratio[jpsi] += 1.0/ratioij;			   //SIMONE
	}                                                          //SIMONE
      }                                                            //SIMONE
      for(int ipsi=0; ipsi< nPsi; ipsi++)			   //SIMONE
	invsumratio[ipsi]=1.0/sumratio[ipsi];			   //SIMONE
      // Only these properties need to be updated                  //SIMONE
      // Using the sum of the ratio Psi^2[j]/Psi^2[iwref]	   //SIMONE
      // because these are number of order 1. Potentially	   //SIMONE
      // the sum of Psi^2[j] can get very big			   //SIMONE
      Properties(LOGPSI) =logpsi[0];		   		   //SIMONE
      Properties(SUMRATIO) = sumratio[0];			   //SIMONE
 
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      ValueType scale = Tau; // ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);	//SIMONE
      for(int ipsi=0; ipsi< nPsi ;ipsi++) {               		//SIMONE
        drift += dgrad[ipsi]*invsumratio[ipsi];                       	//SIMONE
      } 							    	//SIMONE
      drift *= scale; 						   	//SIMONE
      deltaR = (*it)->R - W.R - drift;
      RealType logGb = -oneover2tau*Dot(deltaR,deltaR);
      
      RealType g= Properties(SUMRATIO)/(*it)->Properties(SUMRATIO)*   		//SIMONE
	exp(logGb-logGf+2.0*(Properties(LOGPSI)-(*it)->Properties(LOGPSI)));	//SIMONE
      if(Random() > g) {
	(*it)->Properties(AGE)++;     
	++nReject; 
      } else {
	Properties(AGE) = 0;
	(*it)->R = W.R;
	(*it)->Drift = drift;
	(*it)->Properties = Properties;
	for(int ipsi=0; ipsi< nPsi; ipsi++){ 		             	//SIMONE
	  W.L = lap[ipsi];
	  W.G = dgrad[ipsi];					   //SIMONE
	  RealType eloc=H1[ipsi]->evaluate(W); 	   		   //SIMONE
	  H1[ipsi]->copy((*it)->getEnergyBase(ipsi),invsumratio[ipsi]);			//SIMONE
          //for(int jp=0; jp< np; jp++)					//SIMONE - define np
	  //  (*it)->PsiProperty[jp][ipsi] *= invsumratio[ipsi];		//SIMONE
	}								//SIMONE
	++nAccept;
	for(int ipsi=0; ipsi<nPsi; ipsi++)				//SIMONE
	  totinvsumratio[ipsi] += invsumratio[ipsi];			//SIMONE
      }
      ++it; 
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
