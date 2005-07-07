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
#include "QMCDrivers/DMCParticleByParticle.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/MolecuFixedNodeBranch.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

  /// Constructor.
  DMCParticleByParticle::DMCParticleByParticle(MCWalkerConfiguration& w, 
					       TrialWaveFunction& psi, 
					       QMCHamiltonian& h):
    QMCDriver(w,psi,h),
    PopIndex(-1), EtrialIndex(-1),
    BranchInfo("default") { 
    RootName = "dmc";
    QMCType ="dmc";
  }
  
  void DMCParticleByParticle::setBranchInfo(const string& aname) {
    BranchInfo = aname;
  }

  bool DMCParticleByParticle::run() { 

    //add columns
    IndexType PopIndex = Estimators->addColumn("Population");
    IndexType EtrialIndex = Estimators->addColumn("E_T");
    //write the header
    Estimators->reportHeader();

    MolecuFixedNodeBranch<RealType> brancher(Tau,W.getActiveWalkers());
    //initialize parameters for fixed-node branching
    brancher.put(qmcNode,LogOut);

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
    
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    while(it!=it_end) {
      (*it)->Weight= 1.0;
      (*it)->Multiplicity=1.0;
      ++it;
    }

    //going to add routines to calculate how much we need
    bool require_register =  W.createAuxDataSet();
    int iwalker=0;
    it = W.begin();
    it_end = W.end();
    typedef Walker_t::Buffer_t Buffer_t;

    if(require_register) {
      while(it != it_end) {
        (*it)->DataSet.rewind();
	W.registerData(**it,(*it)->DataSet);
	Psi.registerData(W,(*it)->DataSet);
        //W.DataSet[iwalker]->rewind();
        //W.registerData(**it,*(W.DataSet[iwalker]));
        //Psi.registerData(W,*(W.DataSet[iwalker]));
        ++it;++iwalker;
      } 
    }      

    Estimators->reset();
    
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

          //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
          Walker_t& thisWalker(**it);
          Buffer_t& w_buffer(thisWalker.DataSet);

          thisWalker.Weight = 1.0e0;
          thisWalker.Multiplicity=1.0e0;
          //save old local energy
          ValueType eold(thisWalker.Properties(LOCALENERGY));
          ValueType emixed(eold);

          W.R = thisWalker.R;
          w_buffer.rewind();
          W.copyFromBuffer(w_buffer);
          Psi.copyFromBuffer(W,w_buffer);

          ValueType psi_old(thisWalker.Properties(SIGN));
          ValueType psi(psi_old);

          //create a 3N-Dimensional Gaussian with variance=1
          makeGaussRandom(deltaR);
          bool notcrossed(true);
          int nAcceptTemp(0);
          int nRejectTemp(0);

          int iat=0;
          while(notcrossed && iat<nat){

            PosType dr(g*deltaR[iat]+thisWalker.Drift[iat]);
            PosType newpos(W.makeMove(iat,dr));
            RealType ratio(Psi.ratio(W,iat,dG,dL));

            if(ratio < 0.0) {//node is crossed, stop here
              notcrossed = false;
            } else {
      	      G = W.G+dG;
      	      RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
      	      
      	      ValueType vsq = Dot(G,G);
      	      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      	      dr = thisWalker.R[iat]-newpos-scale*G[iat]; 
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
      	        thisWalker.Drift = scale*G;
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
      	      thisWalker.R = W.R;

              RealType enew= H.evaluate(W);
              thisWalker.resetProperty(log(abs(psi)),psi,enew);
              H.saveProperty(thisWalker.getPropertyBase());
      	      emixed += enew;
      	      //thisWalker.Age = 0;
              ////This is not so useful: allow overflow/underflow
      	      //thisWalker.Properties(LOGPSI) = log(fabs(psi));
      	      //thisWalker.Properties(SIGN) = psi;
      	      //thisWalker.Properties(LOCALENERGY) = H.evaluate(W);
      	      //H.copy(thisWalker.getEnergyBase());
      	      //thisWalker.Properties(LOCALPOTENTIAL) = H.getLocalPotential();
      	      //emixed += thisWalker.Properties(LOCALENERGY);
           } else {
              thisWalker.Age++;
      	      ++nAllRejected;
      	      emixed += eold;
            }
            
            ValueType M = brancher.branchGF(Tau,emixed*0.5,0.0);
            // if((*it)->Properties(AGE) > 3.0) M = min(0.5,M);
            //persistent configurations
           
            //AGE and MULTIPLICITY is removed
            //if(thisWalker.Properties(AGE) > 1.9) M = std::min(0.5,M);
            //if(thisWalker.Properties(AGE) > 0.9) M = std::min(1.0,M);
            //thisWalker.Properties(WEIGHT) = M; 
            //thisWalker.Properties(MULTIPLICITY) = M + Random();
            //if(thisWalker.Age > 1) M = std::min(0.5,M);
            if(thisWalker.Age > 2) M = std::min(0.5,M);
            if(thisWalker.Age > 0) M = std::min(1.0,M);
            thisWalker.Weight = M; 
            thisWalker.Multiplicity=M + Random();
            nAccept += nAcceptTemp;
            nReject += nRejectTemp;
          } else {//set the weight and multiplicity to zero
            thisWalker.willDie();
            nReject += W.getTotalNum();//not make sense
          }

          ++it; ++iwalker;
        }

        ++step;++accstep;
        Estimators->accumulate(W);
        Eest = brancher.update(Population,Eest); 
        //E_T = brancher.update(Population,Eest);
        brancher.branch(accstep,W);
      } while(step<nSteps);
      
      //WARNMSG("The number of a complete rejectoin " << nAllRejected)
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      Estimators->flush();
      
      Estimators->setColumn(PopIndex,static_cast<RealType>(Population));
      Estimators->setColumn(EtrialIndex,Eest); //E_T);
      Estimators->setColumn(AcceptIndex,
      		     static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
      Estimators->report(accstep);
      
      Eest = Estimators->average(0);
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
      		    << static_cast<RealType>(nAllRejected)/static_cast<RealType>(step*W.getActiveWalkers()) << endl;
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
    
    Estimators->finalize();
    
    return true;
  }

  bool 
  DMCParticleByParticle::put(xmlNodePtr q){
    return true;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
