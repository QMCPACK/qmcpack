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
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

  /// Constructor.
  DMCParticleByParticle::DMCParticleByParticle(MCWalkerConfiguration& w, 
					       TrialWaveFunction& psi, 
					       QMCHamiltonian& h):
    QMCDriver(w,psi,h),
    PopIndex(-1), EtrialIndex(-1),
    BranchInfo("default"), branchEngine(0){ 
    RootName = "dmc";
    QMCType ="dmc";
  }
  
  DMCParticleByParticle::~DMCParticleByParticle() {
    if(branchEngine) delete branchEngine;
  }

  bool DMCParticleByParticle::run() { 

    //add columns
    IndexType PopIndex = Estimators->addColumn("Population");
    IndexType EtrialIndex = Estimators->addColumn("Etrial");
    //write the header
    Estimators->reportHeader();

    if(branchEngine == 0) {
      branchEngine=new BranchEngineType(Tau,W.getActiveWalkers());
      RealType e_ref = W.getLocalEnergy();
      branchEngine->setEguess(e_ref);
      branchEngine->put(qmcNode,LogOut);
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

    Estimators->reset();
    
    IndexType block = 0;
    
    Pooma::Clock timer;
    int Population = W.getActiveWalkers();
    int tPopulation = W.getActiveWalkers();
    RealType Eest = branchEngine->E_T;
    oneover2tau = 0.5/Tau;
    sqrttau = sqrt(Tau);
    
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

          thisWalker.Weight = 1.0e0;
          thisWalker.Multiplicity=1.0e0;

          W.R = thisWalker.R;
          DistanceTable::update(W);

          //save old local energy
          ValueType eold(thisWalker.Properties(LOCALENERGY));
          ValueType emixed(eold), enew(eold);

          //Evaluate again to updat G and L
          ValueType psi_old(Psi.evaluateLog(W));
          ValueType sign_old(Psi.getSign());
          G=W.G;
          L=W.L;

          //create a 3N-Dimensional Gaussian with variance=1
          makeGaussRandom(deltaR);
          bool notcrossed(true);
          int nAcceptTemp(0);
          int nRejectTemp(0);
          int iat=0;
          while(iat<nat) {//particle-by-particle move
            PosType dr(sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
            W.R[iat]+=dr;
            //cout << " New position " << thisWalker.R[iat] << g*deltaR[iat] << thisWalker.Drift[iat] << endl;
            DistanceTable::update(W);
            ValueType psi_new(Psi.evaluateLog(W));
            ValueType sign_new(Psi.getSign());

            if(sign_new*sign_old < 0.0) {//node is crossed reject the move
              //cout << accstep << " " << iwalker << " " << iat << " rejection due to node crosssing " << endl;
              ++nRejectTemp;
              //restore properties
              W.R[iat]=thisWalker.R[iat];
              W.G=G;
              W.L=L;
            } else {
      	      RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
      	      //ValueType vsq = Dot(G,G);
      	      //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
              ValueType scale=Tau;
      	      dr = thisWalker.R[iat]-W.R[iat]-scale*W.G[iat]; 
      	      RealType logGb = -oneover2tau*dot(dr,dr);
      	      RealType prob = std::min(1.0,
                  exp(2.0*(psi_new-psi_old)+logGb-logGf));
              //cout << accstep << " " <<iwalker << " " <<  iat << " " 
              //  << exp(2.0*(psi_new-psi_old)) << " " << exp(logGb-logGf)
              //  << " " << prob ;

      	      if(Random() < prob) { 
      	        ++nAcceptTemp;
                thisWalker.R[iat]=W.R[iat];
      	        thisWalker.Drift = scale*W.G;
                psi_old=psi_new;
                sign_old=sign_new;
                //cout << "  accepted " << endl;
      	      } else {
      	        ++nRejectTemp; 
                W.R[iat]=thisWalker.R[iat];
                W.G=G;
                W.L=L;
                //cout << "  rejected " << endl;
      	      }
            } 
            ++iat;
          }

          if(nAcceptTemp>0) {//need to overwrite the walker properties
      	    thisWalker.R = W.R;
            enew= H.evaluate(W);
            thisWalker.resetProperty(psi_old,sign_old,enew);
            H.saveProperty(thisWalker.getPropertyBase());
      	    emixed = (eold+enew)*0.5e0;

            ValueType M = branchEngine->branchGF(Tau,emixed,0.0);
            //if(thisWalker.Age > 3) M = std::min(0.5,M);
            //else if(thisWalker.Age > 0) M = std::min(1.0,M);
            thisWalker.Weight = M; 
            thisWalker.Multiplicity=M + Random();
          } else {
            thisWalker.Age++;
      	    ++nAllRejected;
            thisWalker.Multiplicity=0.0;
          }

          //accumulate the energy
          branchEngine->accumulate(emixed,thisWalker.Weight);
          //cout << "### mixed energy " << emixed << " " << thisWalker.Weight <<  " " << thisWalker.Multiplicity << endl;

          nAccept += nAcceptTemp;
          nReject += nRejectTemp;
          /*
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
      	      ValueType psi = Psi.evaluate(W,w_buffer);
      	      thisWalker.R = W.R;
              RealType enew= H.evaluate(W);
              thisWalker.resetProperty(log(abs(psi)),psi,enew);
              H.saveProperty(thisWalker.getPropertyBase());
      	      emixed += enew;
            } else {
              thisWalker.Age++;
      	      ++nAllRejected;
      	      emixed += eold;
            }
            
            ValueType M = branchEngine->branchGF(Tau,emixed*0.5,0.0);
            // if((*it)->Properties(AGE) > 3.0) M = min(0.5,M);
            //persistent configurations
            //AGE and MULTIPLICITY is removed
            //if(thisWalker.Properties(AGE) > 1.9) M = std::min(0.5,M);
            //if(thisWalker.Properties(AGE) > 0.9) M = std::min(1.0,M);
            //thisWalker.Properties(WEIGHT) = M; 
            //thisWalker.Properties(MULTIPLICITY) = M + Random();
            //if(thisWalker.Age > 1) M = std::min(0.5,M);
            if(thisWalker.Age > 2) M = std::min(0.5,M);
            else if(thisWalker.Age > 0) M = std::min(1.0,M);
            thisWalker.Weight = M; 
            thisWalker.Multiplicity=M + Random();

            //accumulate the energy
            branchEngine->accumulate(emixed*0.5,M);

            nAccept += nAcceptTemp;
            nReject += nRejectTemp;
          } else {//set the weight and multiplicity to zero
            thisWalker.willDie();
            nReject += W.getTotalNum();//not make sense
          }
        */

          ++it; ++iwalker;
        }

        ++step;++accstep;
        Estimators->accumulate(W);
        Eest = branchEngine->update(Population); 
        branchEngine->branch(accstep,W);
      } while(step<nSteps);
      
      //WARNMSG("The number of a complete rejectoin " << nAllRejected)
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      Estimators->flush();
      
      Estimators->setColumn(PopIndex,static_cast<RealType>(Population));

      Estimators->setColumn(EtrialIndex,Eest); 
      Estimators->setColumn(AcceptIndex,
      		     static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
      Estimators->report(accstep);
      
      Eest = Estimators->average(0);
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
      		    << static_cast<RealType>(nAllRejected)/static_cast<RealType>(step*W.getActiveWalkers()) << endl;

      if(pStride) { //create an output engine
        HDFWalkerOutput WO(RootName);
        WO.get(W); 
      }
      nAccept = 0; nReject = 0;
      block++;
    } while(block<nBlocks);

    LogOut->getStream() 
      << "Ratio = " 
      << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << endl;
    
    //branchEngine->updateBranchData();
    if(!pStride) { //create an output engine
      HDFWalkerOutput WO(RootName);
      WO.get(W); 
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
