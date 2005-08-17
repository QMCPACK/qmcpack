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
#include "Particle/DistanceTable.h"
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
  }
  
  bool VMCParticleByParticle::run() { 


    Estimators->reportHeader();
      
    //going to add routines to calculate how much we need
    bool require_register =  W.createAuxDataSet();
    int iwalker=0;
    MCWalkerConfiguration::iterator it(W.begin());
    MCWalkerConfiguration::iterator it_end(W.end());

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

    //create an output engine
    HDFWalkerOutput WO(RootName);
    
    IndexType block = 0;
    
    Pooma::Clock timer;

    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    MCWalkerConfiguration::PropertyContainer_t Properties;
    ParticleSet::ParticleGradient_t G(W.getTotalNum()), dG(W.getTotalNum());
    ParticleSet::ParticleLaplacian_t L(W.getTotalNum()), dL(W.getTotalNum());

    IndexType accstep=0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;

    //  ofstream fout("test.txt");
    do {
      IndexType step = 0;
      timer.start();
      nAccept = 0; nReject=0;
      IndexType nAllRejected = 0;
      do {
        it = W.begin();	 
        it_end = W.end();	 
        iwalker=0; 
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

            PosType dr = g*deltaR[iat]+thisWalker.Drift[iat];
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

            RealType logGb = -oneover2tau*dot(dr,dr);

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
            //thisWalker.Properties(SIGN) = psi;
            //thisWalker.Properties(LOCALENERGY) = H.evaluate(W);
            //H.copy(thisWalker.getEnergyBase());
            //thisWalker.Properties(LOCALPOTENTIAL) = H.getLocalPotential();
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
          ++it; ++iwalker;
        }
        ++step;++accstep;
        Estimators->accumulate(W);
      } while(step<nSteps);
      
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      Estimators->flush();
      Estimators->setColumn(AcceptIndex,
      		   static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
      Estimators->report(accstep);
      
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
      		    << static_cast<RealType>(nAllRejected)/static_cast<RealType>(step*W.getActiveWalkers()) << endl;
      if(pStride) WO.get(W);
      nAccept = 0; nReject = 0;
      ++block;
    } while(block<nBlocks);

    LogOut->getStream() 
      << "Ratio = " 
      << static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot)
      << endl;
    
    if(!pStride) WO.get(W);
    
    Estimators->finalize();
    
    return true;
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
