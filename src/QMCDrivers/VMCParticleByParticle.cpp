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
#include "QMCDrivers/VMCParticleByParticle.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

  /// Constructor.
  VMCParticleByParticle::VMCParticleByParticle(MCWalkerConfiguration& w, 
					       TrialWaveFunction& psi, 
					       QMCHamiltonian& h):
    QMCDriver(w,psi,h) { 
    RootName = "vmc";
    QMCType ="vmc";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
  }
  
  bool VMCParticleByParticle::run() { 


    Estimators->reportHeader(AppendRun);
    Estimators->reset();

    IndexType block = 0;
    
    Pooma::Clock timer;
    m_oneover2tau = 0.5/Tau;
    m_sqrttau = sqrt(Tau);

    G.resize(W.getTotalNum());
    dG.resize(W.getTotalNum());
    L.resize(W.getTotalNum());
    dL.resize(W.getTotalNum());
    
    nAcceptTot = 0;
    nRejectTot = 0;

    bool appendWalker=false;
    int now=0;

    do {
      IndexType step = 0;
      timer.start();
      nAccept = 0; nReject=0;
      nAllRejected = 0;
      do {
        advanceWalkerByWalker();
        ++step;++CurrentStep;
        Estimators->accumulate(W);
        if(CurrentStep%100 == 0) updateWalkers();
      } while(step<nSteps);
      
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      Estimators->flush();
      Estimators->setColumn(AcceptIndex,
      		   static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
      Estimators->report(CurrentStep);
      branchEngine->accumulate(Estimators->average(0),1.0);
      
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
      		    << static_cast<RealType>(nAllRejected)/static_cast<RealType>(step*W.getActiveWalkers()) << endl;

      if(pStride) {
        appendWalker=AppendRun || block>0;
        now=block;
      } 
      HDFWalkerOutput WO(RootName,appendWalker, now);
      WO.get(W);

      nAccept = 0; nReject = 0;
      ++block;

    } while(block<nBlocks);

    branchEngine->update(W.getActiveWalkers(), Estimators->average(0));

    LogOut->getStream() 
      << "Ratio = " 
      << static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot)
      << endl;
    
    int nconf= pStride ? block:1;
    HDFWalkerOutput WOextra(RootName,true,nconf);
    WOextra.write(*branchEngine);

    Estimators->finalize();
    return true;
  }

  void VMCParticleByParticle::advanceWalkerByWalker() {

    MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
    while(it != it_end) {

      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Buffer_t& w_buffer(thisWalker.DataSet);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      ValueType psi_old = thisWalker.Properties(SIGN);
      ValueType psi = psi_old;
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      bool moved = false;

      for(int iat=0; iat<W.getTotalNum(); iat++) {

        PosType dr = m_sqrttau*deltaR[iat]+thisWalker.Drift[iat];
        PosType newpos = W.makeMove(iat,dr);

        //RealType ratio = Psi.ratio(W,iat);
        RealType ratio = Psi.ratio(W,iat,dG,dL);

        G = W.G+dG;

        //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
        //dr = (*it)->R[iat]-newpos-Tau*G[iat]; 
        //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
        RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);

        ValueType vsq = Dot(G,G);
        ValueType scale = ((-1.0e0+sqrt(1.0e0+2.0e0*Tau*vsq))/vsq);
        dr = thisWalker.R[iat]-newpos-scale*G[iat];
        //dr = (*it)->R[iat]-newpos-Tau*G[iat]; 

        RealType logGb = -m_oneover2tau*dot(dr,dr);

        RealType prob = std::min(1.0e0,ratio*ratio*exp(logGb-logGf));
        //alternatively
        if(Random() < prob) { 
          moved = true;
          ++nAccept;
          W.acceptMove(iat);
          //Psi.update(W,iat);
          Psi.update2(W,iat);
          W.G = G;
          W.L += dL;
          //(*it)->Drift = Tau*G;
          thisWalker.Drift = scale*G;
        } else {
          ++nReject; 
          Psi.restore(iat);
        }
      }

      if(moved) {
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        psi = Psi.evaluate(W,w_buffer);

        thisWalker.R = W.R;
        RealType eloc=H.evaluate(W);
        thisWalker.resetProperty(log(abs(psi)),psi,eloc);
        H.saveProperty(thisWalker.getPropertyBase());
      }
      //Keep until everything is tested: debug routine
      // 	    if(moved) {
      // 	      //(*it)->Data.rewind();
      // 	      //W.copyToBuffer((*it)->Data);
      //   	      //psi = Psi.evaluate(W,(*it)->Data);
      // 	      DataSet[iwalker]->rewind();
      // 	      W.copyToBuffer(*DataSet[iwalker]);
      // 	      psi = Psi.evaluate(W,*DataSet[iwalker]);

      // 	      cout << "What is the ratio " << psi/psi_old << endl;
      // 	      // 		 std::cout << "Are they the same ratio=" << psi << endl;

      // 	      RealType please =  H.evaluate(W);

      // 	      /**@note Debug the particle-by-particle move */
      // 	      G = W.G;
      // 	      L = W.L;
      // 	      DistanceTable::update(W);
      // 	      ValueType psinew = Psi.evaluate(W);
      // 	      ValueType dsum=pow(psi-psinew,2);
      // 	      std::cout<< "Diff of updated and calculated gradients/laplacians" << endl;
      // 	      for(int iat=0; iat<W.getTotalNum(); iat++) {
      // 		dsum += pow(G[iat][0]-W.G[iat][0],2)+pow(G[iat][1]-W.G[iat][1],2)
      // 		  +pow(G[iat][2]-W.G[iat][2],2)+pow(L[iat]-W.L[iat],2);
      // 		std::cout << G[iat]-W.G[iat] << " " << L[iat]-W.L[iat] <<endl;
      // 	      }
      // 	      std::cout << "Are they the same ratio=" << psi << " eval=" << psinew << " " << sqrt(dsum) << endl;
      // 	      std::cout << "============================ " << endl;
      // 	      for(int i=0; i<W.getLocalNum(); i++) {
      // 		cout << W.G[i] << " " << W.L[i] << endl;
      // 	      }

      // 	      (*it)->Properties(PSISQ) = psi*psi;
      // 	      (*it)->Properties(PSI) = psi;
      // 	      (*it)->Properties(LOCALENERGY) = H.evaluate(W);
      // 	      std::cout << "Energy " << please << " " << (*it)->Properties(LOCALENERGY) 
      // 			<< endl;
      // 	      H.copy((*it)->getEnergyBase());
      // 	      ValueType vsq = Dot(W.G,W.G);
      // 	      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      // 	      (*it)->Drift = scale*W.G;
      // 	      (*it)->R =W.R;
      //	    } 
      else {
        //WARNMSG("All the particle moves are rejected.")
        ++nAllRejected;
      }
      ++it;
    }

  }
  bool 
  VMCParticleByParticle::put(xmlNodePtr q){
    //nothing to add
    return true;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
