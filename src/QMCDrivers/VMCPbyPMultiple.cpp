//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
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
#include "QMCDrivers/VMCPbyPMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"
#include "Estimators/MultipleEnergyEstimator.h"

namespace ohmmsqmc { 

  /// Constructor.
  VMCPbyPMultiple::VMCPbyPMultiple(MCWalkerConfiguration& w, 
					       TrialWaveFunction& psi, 
					       QMCHamiltonian& h):
    QMCDriver(w,psi,h) { 
    RootName = "vmc";
    QMCType ="vmc";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    QMCDriverMode.set(QMC_MULTIPLE,1);
    add_H_and_Psi(&h,&psi);
  }

  VMCPbyPMultiple::~VMCPbyPMultiple() {
    for(int i=0; i<G.size(); i++) delete G[i];
    for(int i=0; i<dL.size(); i++) delete dL[i];
  }
  
  bool VMCPbyPMultiple::run() { 

    Estimators->reportHeader();
      
    //going to add routines to calculate how much we need
    bool require_register =  W.createAuxDataSet();

    multiEstimator->initialize(W,H1,Psi1,Tau,require_register);

    Estimators->reset();

    IndexType block = 0;
    
    Pooma::Clock timer;

    m_oneover2tau = 0.5/Tau;
    m_sqrttau = sqrt(Tau);
    RealType nPsi_minus_one = nPsi-1;

    ParticleSet::ParticleGradient_t dG(W.getTotalNum());

    IndexType accstep=0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;

    MCWalkerConfiguration::iterator it;
    MCWalkerConfiguration::iterator it_end(W.end());
    bool appendwalker=pStride>0;
    do {  //Blocks loop
      IndexType step = 0;
      timer.start();
      nAccept = 0; nReject=0;
      IndexType nAllRejected = 0;
      do {  //Steps loop
        it = W.begin();	 
        int iwalker=0; 
        while(it != it_end) {  //Walkers loop

          Walker_t& thisWalker(**it);
          Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

          W.R = thisWalker.R;
          w_buffer.rewind();

	  // Copy walker info in W
          W.copyFromBuffer(w_buffer);
          for(int ipsi=0; ipsi<nPsi; ipsi++){
	    // Copy wave function info in W and Psi1
	    Psi1[ipsi]->copyFromBuffer(W,w_buffer);  
	    Psi1[ipsi]->G=W.G;
	    Psi1[ipsi]->L=W.L;
          }

	  // Point to the correct walker in the ratioij buffer
	  RealType *ratioijPtr=multiEstimator->RatioIJ[iwalker];

          //This is not used
          //ValueType psi_old = thisWalker.Properties(SIGN);
          //ValueType psi = psi_old;
          //create a 3N-Dimensional Gaussian with variance=1
          makeGaussRandom(deltaR);
          bool moved = false;

          for(int iat=0; iat<W.getTotalNum(); iat++) {  //Particles loop

            PosType dr = m_sqrttau*deltaR[iat]+thisWalker.Drift[iat];
            PosType newpos = W.makeMove(iat,dr);

	    for(int ipsi=0; ipsi<nPsi; ipsi++){
	      // Compute ratios before and after the move
	      ratio[ipsi] = Psi1[ipsi]->ratio(W,iat,dG,*dL[ipsi]); 
              logpsi2[ipsi]=log(ratio[ipsi]*ratio[ipsi]);
	      // Compute Gradient in new position
              *G[ipsi]=Psi1[ipsi]->G + dG;
	      // Initialize: sumratio[i]=(Psi[i]/Psi[i])^2=1.0
	      sumratio[ipsi]=1.0;
	    }

	    // Compute new (Psi[i]/Psi[j])^2 and their sum
	    int indexij(0);
	    for(int ipsi=0; ipsi< nPsi_minus_one; ipsi++){
	      for(int jpsi=ipsi+1; jpsi < nPsi; jpsi++, indexij++){
		//RealType rji=ratio[jpsi]/ratio[ipsi];
		//rji = rji*rji*ratioijPtr[indexij]; 
                RealType rji=exp(logpsi2[jpsi]-logpsi2[ipsi])*ratioijPtr[indexij];
		ratioij[indexij]=rji;
		sumratio[ipsi] += rji;
		sumratio[jpsi] += 1.0/rji;
	      }
	    }

            RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
	    drift=0.0;
	    // Evaluate new Umbrella Weight and new drift
	    for(int ipsi=0; ipsi< nPsi; ipsi++){
	      invsumratio[ipsi]=1.0/sumratio[ipsi];
	      drift += invsumratio[ipsi]*(*G[ipsi]);
	    }
            ValueType vsq = Dot(drift,drift);
            ValueType scale = ((-1.0e0+sqrt(1.0e0+2.0e0*Tau*vsq))/vsq);
            //ValueType scale=Tau;
	    drift *= scale;
            dr = thisWalker.R[iat]-newpos-drift[iat];
            RealType logGb = -m_oneover2tau*dot(dr,dr);

	    // td = Target Density ratio
	    //RealType td=pow(ratio[0],2)*sumratio[0]/(*it)->Properties(SUMRATIO);
	    RealType td=ratio[0]*ratio[0]*sumratio[0]/(*it)->Multiplicity;
	    RealType prob = std::min(1.0,td*exp(logGb-logGf));

	    if(Random() < prob) { 
	      /* Electron move is accepted. Update:
		 -ratio (Psi[i]/Psi[j])^2 for this walker
		 -Gradient and laplacian for each Psi1[i]
		 -Drift
		 -buffered info for each Psi1[i]*/
	      moved = true;
	      ++nAccept;
	      W.acceptMove(iat);
	      // Update Buffer for (Psi[i]/Psi[j])^2 
	      std::copy(ratioij.begin(),ratioij.end(),ratioijPtr);
	      // Update Umbrella weight
	      UmbrellaWeight=invsumratio;

	      // Store sumratio for next Accept/Reject step to Multiplicity
	      //thisWalker.Properties(SUMRATIO)=sumratio[0];
	      thisWalker.Multiplicity=sumratio[0];
	      for(int ipsi=0; ipsi< nPsi; ipsi++){
		////Update local Psi1[i] buffer for the next move
		Psi1[ipsi]->update2(W,iat);  
		// Update G and L in Psi1[i]
		Psi1[ipsi]->G = *G[ipsi];
		Psi1[ipsi]->L += *dL[ipsi];
                thisWalker.Properties(ipsi,LOGPSI)+=log(abs(ratio[ipsi]));
	      }
	      // Update Drift
	      (*it)->Drift = drift;
	    } else {
	      ++nReject;
	      for(int ipsi=0; ipsi< nPsi; ipsi++)
		Psi1[ipsi]->restore(iat);
	    }
	  }

	  if(moved) {
	    /* The walker moved: Info are copied back to buffers:
	       -copy (Psi[i]/Psi[j])^2 to ratioijBuffer
	       -Gradient and laplacian for each Psi1[i]
	       -Drift
	       -buffered info for each Psi1[i]
	       Physical properties are updated */
            (*it)->Age=0;
	    (*it)->R = W.R;
	    w_buffer.rewind();
	    W.copyToBuffer(w_buffer);
	    for(int ipsi=0; ipsi< nPsi; ipsi++){
	      W.G=Psi1[ipsi]->G;
	      W.L=Psi1[ipsi]->L;
	      ValueType psi = Psi1[ipsi]->evaluate(W,w_buffer);
	      RealType et = H1[ipsi]->evaluate(W);

	      //multiEstimator->updateSample(iwalker,ipsi,et,UmbrellaWeight[ipsi]);
              //Properties is used for UmbrellaWeight and UmbrellaEnergy
              thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=UmbrellaWeight[ipsi];
              thisWalker.Properties(ipsi,LOCALENERGY)=et;
              H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
	    }
	  }
	  else {
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
	<< static_cast<RealType>(nAllRejected)/static_cast<RealType>(step*W.getActiveWalkers()) << 
	" nPsi " << nPsi << endl;

      HDFWalkerOutput WO(RootName,block&&appendwalker, block);
      WO.get(W);

      nAccept = 0; nReject = 0;
      ++block;

      //re-evaluate the ratio
      multiEstimator->initialize(W,H1,Psi1,Tau,false);
    } while(block<nBlocks);

    LogOut->getStream() 
      << "Ratio = " 
      << static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot)
      << endl;

    int nconf= appendwalker ? block:1;
    HDFWalkerOutput WOextra(RootName,true,nconf);
    WOextra.write(*branchEngine);

    Estimators->finalize();

    return true;
  }


  bool 
  VMCPbyPMultiple::put(xmlNodePtr q){
    nPsi=Psi1.size();
    resize(nPsi,W.getTotalNum());

    for(int ipsi=0; ipsi<nPsi; ipsi++) 
      H1[ipsi]->add2WalkerProperty(W);

    if(Estimators == 0) {
      Estimators = new ScalarEstimatorManager(H);
      multiEstimator = new MultipleEnergyEstimator(H,nPsi);
      Estimators->add(multiEstimator,"elocal");
    }

    H1[0]->setPrimary(true);
    for(int ipsi=1; ipsi<nPsi; ipsi++) {
      H1[ipsi]->setPrimary(false);
    }
    return true;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
