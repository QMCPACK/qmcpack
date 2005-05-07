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
#include "QMC/MolecuFixedNodeBranch.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

  /// Constructor.
  DMCParticleByParticle::DMCParticleByParticle(MCWalkerConfiguration& w, 
					       TrialWaveFunction& psi, 
					       QMCHamiltonian& h, 
					       xmlNodePtr q): 
    QMCDriver(w,psi,h,q),BranchInfo("default") { 
    RootName = "dmc";
    QMCType ="dmc";
  }
  
  void DMCParticleByParticle::setBranchInfo(const string& aname) {
    BranchInfo = aname;
  }
  bool DMCParticleByParticle::run() { 

    W.setUpdateMode(MCWalkerConfiguration::Update_Particle);
    //DistanceTable::create(1);

    if(put(qmc_node)){
      
      //set the data members to start a new run
      // getReady();
      MolecuFixedNodeBranch<RealType> brancher(Tau,W.getActiveWalkers());

      //initialize parameters for fixed-node branching
      brancher.put(qmc_node,LogOut);

      if(BranchInfo != "default")  brancher.read(BranchInfo);
      else {
        /*if VMC/DMC directly preceded DMC (Counter > 0) then
          use the average value of the energy estimator for
          the reference energy of the brancher*/
        if(Counter) {
          RealType e_ref = W.getLocalEnergy();
          LOGMSG("Overwriting the reference energy by the local energy " << e_ref)  
          brancher.setEguess(e_ref);
        }
      }
      
      int PopIndex, E_TIndex;
      Estimators.resetReportSettings(RootName);
      AcceptIndex = Estimators.addColumn("AcceptRatio");
      PopIndex = Estimators.addColumn("Population");
      E_TIndex = Estimators.addColumn("E_T");
      Estimators.reportHeader();
      
      MCWalkerConfiguration::iterator it(W.begin()); 
      MCWalkerConfiguration::iterator it_end(W.end()); 
      while(it!=it_end) {
        (*it)->Properties(WEIGHT) = 1.0;
        (*it)->Properties(MULTIPLICITY) = 1.0;
        ++it;
      }

      //going to add routines to calculate how much we need
      bool require_register =  W.createAuxDataSet();
      int iwalker=0;
      it = W.begin();
      it_end = W.end();
      if(require_register) {
        while(it != it_end) {
	  W.DataSet[iwalker]->rewind();
  	  W.registerData(**it,*(W.DataSet[iwalker]));
  	  Psi.registerData(W,*(W.DataSet[iwalker]));
	  ++it;++iwalker;
        } 
      }      

      Estimators.reset();
      
      IndexType block = 0;
      
      Pooma::Clock timer;
      int Population = W.getActiveWalkers();
      int tPopulation = W.getActiveWalkers();
      RealType Eest = brancher.E_T;
      RealType E_T = Eest;
      RealType oneovertau = 1.0/Tau;
      RealType oneover2tau = 0.5*oneovertau;
      RealType g = sqrt(Tau);
      
      MCWalkerConfiguration::PropertyContainer_t Properties;
      ParticleSet::ParticlePos_t deltaR(W.getTotalNum());
      ParticleSet::ParticleGradient_t G(W.getTotalNum()), dG(W.getTotalNum());
      ParticleSet::ParticleLaplacian_t L(W.getTotalNum()), dL(W.getTotalNum());

      IndexType accstep=0;
      IndexType nAcceptTot = 0;
      IndexType nRejectTot = 0;
      IndexType nat = W.getTotalNum();

      int ncross = 0;
      do {
	IndexType step = 0;
	timer.start();
	nAccept = 0; nReject=0;
	IndexType nAllRejected = 0;
	do {
	  Population = W.getActiveWalkers();

	  it = W.begin();	 
          it_end = W.end();

	  iwalker=0; 
	  while(it != it_end) {

	    MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);

	    (*it)->Properties(WEIGHT) = 1.0;
	    (*it)->Properties(MULTIPLICITY) = 1.0;
	    //save old local energy
	    ValueType eold((*it)->Properties(LOCALENERGY));
	    ValueType emixed(eold);

	    W.R = (*it)->R;
	    w_buffer.rewind();
	    W.copyFromBuffer(w_buffer);
	    Psi.copyFromBuffer(W,w_buffer);

            ValueType psi_old((*it)->Properties(SIGN));
	    ValueType psi(psi_old);

	    //create a 3N-Dimensional Gaussian with variance=1
	    makeGaussRandom(deltaR);
	    bool notcrossed(true);
            int nAcceptTemp(0);
            int nRejectTemp(0);

            int iat=0;
            while(notcrossed && iat<nat){

	      PosType dr(g*deltaR[iat]+(*it)->Drift[iat]);
	      PosType newpos(W.makeMove(iat,dr));
	      RealType ratio(Psi.ratio(W,iat,dG,dL));

	      if(ratio < 0.0) {//node is crossed, stop here
		notcrossed = false;
	      } else {
		G = W.G+dG;
		RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
		
		ValueType vsq = Dot(G,G);
		ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
		dr = (*it)->R[iat]-newpos-scale*G[iat]; 
		//dr = (*it)->R[iat]-newpos-Tau*G[iat]; 
		RealType logGb = -oneover2tau*dot(dr,dr);
		
		//RealType ratio2 = pow(ratio,2)
		RealType prob = std::min(1.0,pow(ratio,2)*exp(logGb-logGf));
		if(Random() < prob) { 
		  ++nAcceptTemp;
		  W.acceptMove(iat);
		  Psi.update2(W,iat);
		  W.G = G;
		  W.L += dL;
		  //  (*it)->Drift = Tau*G;
		  (*it)->Drift = scale*G;
		} else {
		  ++nRejectTemp; 
		  Psi.restore(iat);
		}
	      } 
              ++iat;
	    }

	    if(notcrossed) {
	      if(nAcceptTemp) {//need to overwrite the walker properties
		w_buffer.rewind();
		W.copyToBuffer(w_buffer);
		psi = Psi.evaluate(W,w_buffer);
		(*it)->R = W.R;
		(*it)->Properties(AGE) = 0;
                //This is not so useful: allow overflow/underflow
		(*it)->Properties(LOGPSI) = log(fabs(psi));
		(*it)->Properties(SIGN) = psi;
		(*it)->Properties(LOCALENERGY) = H.evaluate(W);
		H.copy((*it)->getEnergyBase());
		(*it)->Properties(LOCALPOTENTIAL) = H.getLocalPotential();
		emixed += (*it)->Properties(LOCALENERGY);
	      } else {
		//WARNMSG("All the particle moves are rejected.")
	        (*it)->Properties(AGE)++;
		++nAllRejected;
		emixed += eold;
	      }
	      
	      ValueType M = brancher.branchGF(Tau,emixed*0.5,0.0);
	      // if((*it)->Properties(AGE) > 3.0) M = min(0.5,M);
	      //persistent configurations
	      if((*it)->Properties(AGE) > 1.9) M = std::min(0.5,M);
	      if((*it)->Properties(AGE) > 0.9) M = std::min(1.0,M);
	      (*it)->Properties(WEIGHT) = M; 
	      (*it)->Properties(MULTIPLICITY) = M + Random();
	      nAccept += nAcceptTemp;
	      nReject += nRejectTemp;
	    } else {//set the weight and multiplicity to zero
	      (*it)->Properties(WEIGHT) = 0.0; 
	      (*it)->Properties(MULTIPLICITY) = 0.0;
	      nReject += W.getTotalNum();//not make sense
	    }

	    ++it; ++iwalker;
	  }

	  ++step;++accstep;
	  Estimators.accumulate(W);
	  Eest = brancher.update(Population,Eest); 
          //E_T = brancher.update(Population,Eest);
	  brancher.branch(accstep,W);
	} while(step<nSteps);
	
        //WARNMSG("The number of a complete rejectoin " << nAllRejected)
	timer.stop();
	nAcceptTot += nAccept;
	nRejectTot += nReject;
	
	Estimators.flush();
	
	Estimators.setColumn(PopIndex,static_cast<double>(Population));
	Estimators.setColumn(E_TIndex,Eest); //E_T);
	Estimators.setColumn(AcceptIndex,
			     static_cast<double>(nAccept)/static_cast<double>(nAccept+nReject));
	Estimators.report(accstep);
	
	Eest = Estimators.average(0);
	LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
			    << static_cast<double>(nAllRejected)/static_cast<double>(step*W.getActiveWalkers()) << endl;
	if(pStride) {
          //create an output engine
	  HDFWalkerOutput WO(RootName);
          WO.get(W); 
          brancher.write(WO.getGroupID());
        }
	nAccept = 0; nReject = 0;
	block++;
      } while(block<nBlocks);

      LogOut->getStream() 
	<< "Ratio = " 
	<< static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
	<< endl;
      
      if(!pStride) {
        //create an output engine
        HDFWalkerOutput WO(RootName);
        WO.get(W); 
        brancher.write(WO.getGroupID());
      }
      
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
