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

    //this should be controlled by continue="yes/no" of qmc element
    branchEngine->flush(0);

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
        ++it;++iwalker;
      } 
    } else {
      while(it != it_end) {
        Buffer_t& w_buffer((*it)->DataSet);
      	w_buffer.rewind();
        W.updateBuffer(**it,w_buffer);
        ValueType logpsi=Psi.updateBuffer(W,w_buffer);
        RealType enew= H.evaluate(W);
        (*it)->resetProperty(logpsi,Psi.getSign(),enew);
        H.saveProperty((*it)->getPropertyBase());
        ++it;
      }
    }

    Estimators->reset();
    
    IndexType block = 0;
    
    Pooma::Clock timer;
    int Population = W.getActiveWalkers();
    int tPopulation = W.getActiveWalkers();
    RealType Eest = branchEngine->E_T;
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
          ValueType emixed(eold), enew(eold);

          W.R = thisWalker.R;
          w_buffer.rewind();
          W.copyFromBuffer(w_buffer);
          Psi.copyFromBuffer(W,w_buffer);

          //create a 3N-Dimensional Gaussian with variance=1
          makeGaussRandom(deltaR);
          bool notcrossed(true);
          int nAcceptTemp(0);
          int nRejectTemp(0);
          int iat=0;
         while(iat<nat) {//particle-by-particle move
           PosType dr(g*deltaR[iat]+thisWalker.Drift[iat]);
           PosType newpos(W.makeMove(iat,dr));

           RealType ratio=Psi.ratio(W,iat,dG,dL);

           if(ratio < 0.0) {//node is crossed reject the move
             ++nRejectTemp;
    	     Psi.restore(iat);
           } else {
      	     G = W.G+dG;
      	     RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
      	     //ValueType vsq = Dot(G,G);
      	     //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
             ValueType scale=Tau;
      	     dr = thisWalker.R[iat]-newpos-scale*G[iat]; 
      	     RealType logGb = -oneover2tau*dot(dr,dr);
      	     RealType prob = std::min(1.0,ratio*ratio*exp(logGb-logGf));
      	     if(Random() < prob) { 
      	       ++nAcceptTemp;
      	       W.acceptMove(iat);
      	       Psi.update2(W,iat);
      	       W.G = G;
      	       W.L += dL;
      	       thisWalker.Drift = scale*G;
      	     } else {
      	       ++nRejectTemp; 
      	       Psi.restore(iat);
      	     }
           } 

           ++iat;
         }

         if(nAcceptTemp>0) {//need to overwrite the walker properties
      	   w_buffer.rewind();
      	   W.copyToBuffer(w_buffer);
      	   ValueType psi = Psi.evaluate(W,w_buffer);
      	   thisWalker.R = W.R;
           enew= H.evaluate(W);
           thisWalker.resetProperty(log(abs(psi)),psi,enew);
           H.saveProperty(thisWalker.getPropertyBase());
      	   emixed = (eold+enew)*0.5e0;
         } else {
           thisWalker.Age++;
           ++nAllRejected;
         }

         ValueType M = branchEngine->branchGF(Tau,emixed,0.0);
         if(thisWalker.Age > 3) M = std::min(0.5,M);
         else if(thisWalker.Age > 0) M = std::min(1.0,M);
         thisWalker.Weight = M; 
         thisWalker.Multiplicity=M + Random();

         branchEngine->accumulate(emixed,thisWalker.Weight);//accumulate the energy
         // cout << "### mixed energy " << emixed << " " << thisWalker.Weight <<  " " << thisWalker.Multiplicity << endl;
         nAccept += nAcceptTemp;
         nReject += nRejectTemp;
        
         /** killing any walker with a node-crossing mode
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

      //re-evaluate the buffer every block
      it=W.begin(); 
      it_end=W.end(); 
      while(it != it_end) {
        Buffer_t& w_buffer((*it)->DataSet);
      	w_buffer.rewind();
        W.updateBuffer(**it,w_buffer);
        ValueType logpsi=Psi.updateBuffer(W,w_buffer);
        RealType enew= H.evaluate(W);
        (*it)->resetProperty(logpsi,Psi.getSign(),enew);
        H.saveProperty((*it)->getPropertyBase());
        ++it;
      }

      nAccept = 0; nReject = 0;
      block++;
    } while(block<nBlocks);

    LogOut->getStream() 
      << "Ratio = " 
      << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << endl;
    
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
