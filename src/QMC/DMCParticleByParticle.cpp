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
#include "QMC/DMCParticleByParticle.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMC/MolecuFixedNodeBranch.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

    /// Constructor.
    DMCParticleByParticle::DMCParticleByParticle(MCWalkerConfiguration& w, 
						 TrialWaveFunction& psi, 
						 QMCHamiltonian& h, 
						 xmlNodePtr q): 
	QMCDriver(w,psi,h,q) { 
	RootName = "dmc";
	QMCType ="dmc";
    }
  
    bool DMCParticleByParticle::run() { 

	W.setUpdateMode(MCWalkerConfiguration::Update_Particle);
	//DistanceTable::create(1);

	if(put(qmc_node)){
      
	    //set the data members to start a new run
	    //    getReady();
	    int PopIndex, E_TIndex;
	    Estimators.resetReportSettings(RootName);
	    AcceptIndex = Estimators.addColumn("AcceptRatio");
	    PopIndex = Estimators.addColumn("Population");
	    E_TIndex = Estimators.addColumn("E_T");
	    Estimators.reportHeader();

	    for(MCWalkerConfiguration::iterator it = W.begin(); 
		it != W.end(); ++it) {
		(*it)->Properties(WEIGHT) = 1.0;
		(*it)->Properties(MULTIPLICITY) = 1.0;
	    }

	    MolecuFixedNodeBranch<RealType> brancher(Tau,W.getActiveWalkers());
	    //initialize parameters for fixed-node branching
	    brancher.put(qmc_node,LogOut);
	
	    /*if VMC/DMC directly preceded DMC (Counter > 0) then
	      use the average value of the energy estimator for
	      the reference energy of the brancher*/
	    if(Counter) {
		RealType e_ref = W.getLocalEnergy();
		LOGMSG("Overwriting the reference energy by the local energy " << e_ref)  
		    brancher.setEguess(e_ref);
	    }


	    //going to add routines to calculate how much we need
	    bool require_register =  W.createAuxDataSet();
	    int iwalker=0;
	    MCWalkerConfiguration::iterator it = W.begin();
	    if(require_register) {
		while(it != W.end()) {
		    W.DataSet[iwalker]->rewind();
		    W.registerData(**it,*(W.DataSet[iwalker]));
		    Psi.registerData(W,*(W.DataSet[iwalker]));
		    it++;iwalker++;
		} 
	    }      

	    Estimators.reset();

	    //create an output engine
	    HDFWalkerOutput WO(RootName);
      
	    IndexType block = 0;
      
	    Pooma::Clock timer;
	    int Population = W.getActiveWalkers();
	    int tPopulation = W.getActiveWalkers();
	    RealType E_T = brancher.E_T;
	    RealType Eest = E_T;
	    RealType oneovertau = 1.0/Tau;
	    RealType oneover2tau = 0.5*oneovertau;
	    RealType g = sqrt(Tau);
      
	    MCWalkerConfiguration::PropertyContainer_t Properties;
	    ParticleSet::ParticlePos_t deltaR(W.getTotalNum());
	    ParticleSet::ParticlePos_t drift(W.getTotalNum());
	    ParticleSet::ParticleGradient_t G(W.getTotalNum()), dG(W.getTotalNum());
	    ParticleSet::ParticleLaplacian_t L(W.getTotalNum()), dL(W.getTotalNum());

	    IndexType accstep=0;
	    IndexType nAcceptTot = 0;
	    IndexType nRejectTot = 0;

	    drift = 0.0;

	    do {
		IndexType step = 0;
		timer.start();
		nAccept = 0; nReject=0;
		IndexType nAllRejected = 0;
		do {
		    Population = W.getActiveWalkers();
		    it = W.begin();	 
		    iwalker=0; 
		    while(it != W.end()) {


			(*it)->Properties(WEIGHT) = 1.0;
			(*it)->Properties(MULTIPLICITY) = 1.0;
	      
			//save old local energy
			ValueType eold = (*it)->Properties(LOCALENERGY);
			ValueType emixed = eold;  

			MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
			W.R = (*it)->R;
			w_buffer.rewind();
			W.copyFromBuffer(w_buffer);
			Psi.copyFromBuffer(W,w_buffer);

			ValueType psi_old = (*it)->Properties(PSI);
			ValueType psi = psi_old;
			//create a 3N-Dimensional Gaussian with variance=1
			makeGaussRandom(deltaR);
			bool moved = false;
			bool crossed = false;

			RealType reject;
			int nAcceptTemp = 0;
			int nRejectTemp = 0;

			for(int iat=0; iat<W.getTotalNum(); iat++) {

			    PosType dr = g*deltaR[iat]+(*it)->Drift[iat];
			    PosType Rp = (*it)->R[iat] + dr;
			    //PosType dr = g*deltaR[iat];
			    W.makeMove(iat,dr);

			    RealType ratio = Psi.ratio(W,iat,dG,dL);
			    RealType prob = 0.0;
			    if(ratio < 0.0) crossed = true;
			
				G = W.G+dG;
				RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));

				//ValueType vsq = Dot(G,G);
				//ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
				//	drift = scale*G;
				drift = Tau*G;

				deltaR[iat] = (*it)->R[iat] - Rp - drift[iat];
				RealType backwardGF = exp(-oneover2tau*dot(deltaR[iat],deltaR[iat]));
				/* green function
				   deltaR = -scale*G; 
				   dletaR[iat] -= dr; //subtract dr = (*it)->R(iat) - W.R(iat)
				   RealType backwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
				   //a better way to do this
				   RealType forwardGF = exp(-0.5*dot(dr,dr));
				   ValueType vsq = Dot(G,G);
				   ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
				   ValueType scale2 = scale*scale;
				   RealType dist = scale2*vsq+2.0*scale*dot(dr,G(iat))+dot(dr,dr);
				   RealType backwardGF = exp(-oneover2tau*dist);
				*/

				RealType ratio2 = pow(ratio,2);
				//set acceptance probability
				prob= min(backwardGF/forwardGF*ratio2,1.0);
				reject = 1.0 - prob;
				if(Random() < prob) {
				    moved = true;
				    ++nAcceptTemp;
				    W.acceptMove(iat);
				    //Psi.update(W,iat);
				    Psi.update2(W,iat);
				    W.G = G;
				    W.L += dL;
				    //Need to change the drift
				    (*it)->Drift = drift;
				} else {
				    ++nRejectTemp; 
				    Psi.restore(iat);
				}
			}//particles
			
			if(moved) {
			    //(*it)->Data.rewind();
			    //W.copyToBuffer((*it)->Data);
			    //psi = Psi.evaluate(W,(*it)->Data);
			    w_buffer.rewind();
			    W.copyToBuffer(w_buffer);
			    psi = Psi.evaluate(W,w_buffer);
			    
			    (*it)->R = W.R;
			    // (*it)->Drift = drift;
			    (*it)->Properties(AGE) = 0;
			    (*it)->Properties(PSISQ) = psi*psi;
			    (*it)->Properties(PSI) = psi;
			    (*it)->Properties(LOCALENERGY) = H.evaluate(W);
			    H.copy((*it)->getEnergyBase());
			    (*it)->Properties(LOCALPOTENTIAL) = H.getLocalPotential();
			    emixed += (*it)->Properties(LOCALENERGY);
			} else {
			    WARNMSG("All the particle moves are rejected.")
				nAllRejected++;
			    (*it)->Properties(AGE)++;
			    emixed += eold;
			}
			//calculate the weight and multiplicity
			reject = 0.0; //temporary
			ValueType M = brancher.branchGF(Tau,emixed*0.5,reject);

			if((*it)->Properties(AGE) > 3.0) M = min(0.5,M);
			if((*it)->Properties(AGE) > 0.9) M = min(1.0,M);
			(*it)->Properties(WEIGHT) = M; 
			(*it)->Properties(MULTIPLICITY) = M + Random();

			//node-crossing: kill it for the time being
			if(brancher(psi_old,(*it)->Properties(PSI))) {
			    cout << "crossed node" << endl;
			nAllRejected++;
			nAcceptTemp = 0;
			nRejectTemp += W.getTotalNum();
			// accepted=false;     
			(*it)->Properties(WEIGHT) = 0.0; 
			(*it)->Properties(MULTIPLICITY) = 0.0;
		    }

//			cout << "M = " << (*it)->Properties(MULTIPLICITY) << endl;
			nAccept += nAcceptTemp;
			nReject += nRejectTemp;
			it++; iwalker++;
		    }//walkers
		    step++;accstep++;
		    Estimators.accumulate(W);
//		    E_T = brancher.update(Population,Eest);
		    brancher.branch(accstep,W);
		} while(step<nSteps);
	
		timer.stop();
		nAcceptTot += nAccept;
		nRejectTot += nReject;
	
		Estimators.flush();
		Eest = Estimators.average(0);
	
		Estimators.setColumn(PopIndex,static_cast<double>(Population));
		Estimators.setColumn(E_TIndex,E_T);
		Estimators.setColumn(AcceptIndex,
				     static_cast<double>(nAccept)/static_cast<double>(nAccept+nReject));
		Estimators.report(accstep);
	
		LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
				    << static_cast<double>(nAllRejected)/static_cast<double>(step*W.getActiveWalkers()) << endl;
		if(pStride) WO.get(W);
		nAccept = 0; nReject = 0;
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
    DMCParticleByParticle::put(xmlNodePtr q){
	xmlNodePtr qsave=q;
	bool success = putQMCInfo(q);
	success = Estimators.put(qsave);
	return success;
    }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
