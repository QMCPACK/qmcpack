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
#include "QMC/VMCMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

//This macro is a temporary one to test the efficiency of all-walker moves. So far, not as much
//as to convince JK to implement it globally.
#define MOVE_ONE

  /// Constructor.
  VMCMultiple::VMCMultiple(MCWalkerConfiguration& w, 
	   TrialWaveFunction& psi, 
	   QMCHamiltonian& h, 
	   xmlNodePtr q): 
    QMCDriver(w,psi,h,q) { 
    RootName = "vmc";
    QMCType ="vmc";
    add_H_and_Psi(&h,&psi);
  }
  
  /** Run the VMCMultiple algorithm.
   * \param nblocks number of blocks
   * \param nsteps number of steps
   * \param tau the timestep
   *
   * Advance the walkers nblocks*nsteps timesteps. 
   * For each timestep:
   * <ul>
   * <li> Move all the particles of a walker.
   * <li> Calculate the properties for the new walker configuration.
   * <li> Accept/reject the new configuration.
   * <li> Accumulate the estimators.
   * </ul>
   * For each block:
   * <ul>
   * <li> Flush the estimators and print to file.
   * <li> (Optional) Print the ensemble of walker configurations.
   * </ul>
   *
   * Default mode: Print the ensemble of walker configurations 
   * at the end of the run.
   */
  bool VMCMultiple::run() { 
    
    DistanceTable::create(1);

    if(put(qmc_node)){

      H.setTau(0.0);
      
      //set the data members to start a new run
      getReady();
      
      //probably unnecessary
      MCWalkerConfiguration::iterator it = W.begin(); 
      MCWalkerConfiguration::iterator it_end = W.end(); 
      while(it != it_end) {
	(*it)->Properties(WEIGHT) = 1.0;
        ++it;
      }

      deltaR.resize(W.getTotalNum());
      drift.resize(W.getTotalNum());

      nPsi=Psi1.size();			//SIMONE
      IndexPsi.resize(nPsi);		//SIMONE
      logpsi.resize(nPsi);		//SIMONE
      dgrad.resize(nPsi);		//SIMONE
      lap.resize(nPsi);			//SIMONE
      sumratio.resize(nPsi);		//SIMONE
      invsumratio.resize(nPsi);		//SIMONE
      totinvsumratio.resize(nPsi);	//SIMONE
      
      Estimators.reset();

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
	  Estimators.accumulate(W);
	} while(step<nSteps);
	
	timer.stop();
	nAcceptTot += nAccept;
	nRejectTot += nReject;
	double BlockTotal=static_cast<double>(nAccept+nReject);				//SIMONE
	Estimators.flush();
	Estimators.setColumn(AcceptIndex,static_cast<double>(nAccept)/BlockTotal);	//SIMONE
	for(int ipsi=0; ipsi< nPsi; ipsi++)						//SIMONE
	  Estimators.setColumn(IndexPsi[ipsi],totinvsumratio[ipsi]/BlockTotal);		//SIMONE
	Estimators.report(accstep);
	
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
      
      Estimators.finalize();
      
      return true;
    } else {
      ERRORMSG("Error with Input")
      return false;
    }
  }


  bool 
  VMCMultiple::put(xmlNodePtr q){
    xmlNodePtr qsave=q;
    bool success = putQMCInfo(q);

    //check <estimator>'s
    success = Estimators.put(qsave);

    for(int ipsi=0; ipsi< Psi1.size(); ipsi++){
      char cname[6];
      int jpsi=sprintf(cname,"WPsi%d",ipsi);
      IndexPsi[ipsi]=Estimators.addColumn(cname);
    }
    return success;
  }

  /**  Advance all the walkers one timstep. 
   *
   Propose a move for each walker from its old 
   position \f${\bf R'}\f$ to a new position \f${\bf R}\f$
   \f[ 
   {\bf R'} + {\bf \chi} + 
   \tau {\bf v_{drift}}({\bf R'}) =  {\bf R},
   \f]
   where \f$ {\bf \chi} \f$ is a 3N-diminsional Gaussian 
   of mean zero and variance \f$ \tau \f$ and 
   \f$ {\bf v_{drift}} \f$ is the drift velocity
   \f[ 
   {\bf v_{drift}}({\bf R'}) = {\bf \nabla} 
   \ln |\Psi_T({\bf R'})| = \Psi_T({\bf R'})^{-1} 
   {\bf \nabla} \Psi_T({\bf R'}). 
   \f]
   Metropolis accept/reject with probability
   \f[
   P_{accept}(\mathbf{R'}\rightarrow\mathbf{R}) = 
   \min\left[1,\frac{G(\mathbf{R}\rightarrow\mathbf{R'})
   \Psi_T(\mathbf{R})^2}{G(\mathbf{R'}\rightarrow\mathbf{R})
   \Psi_T(\mathbf{R'})^2}\right],
   \f] 
   where \f$ G \f$ is the drift-diffusion Green's function 
   \f[
   G(\mathbf{R'} \rightarrow 
   \mathbf{R}) = (2\pi\tau)^{-3/2}\exp \left[ -
   (\mathbf{R}-\mathbf{R'}-\tau \mathbf{v_{drift}}
   (\mathbf{R'}))^2/2\tau \right].
   \f]
   If the move is accepted, update the walker configuration and
   properties.  For rejected moves, do not update except for the
   Age which needs to be incremented by one.
   *
   */

  /**@ingroup advanceWalkersVMCMultiple
   * \brief Loop through the walkers and advance one at a time.
   */
  
  void 
  VMCMultiple::advanceWalkerByWalker() {
    
    //Pooma::Clock timer;
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    MCWalkerConfiguration::PropertyContainer_t Properties;
    int nh = H.size()+1;  //SIMONE - What is this?
    
    MCWalkerConfiguration::iterator it = W.begin(); 
    MCWalkerConfiguration::iterator it_end = W.end(); 
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
