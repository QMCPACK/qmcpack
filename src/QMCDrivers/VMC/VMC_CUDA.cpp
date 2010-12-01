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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus { 

  /// Constructor.
  VMCcuda::VMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
		   QMCHamiltonian& h,WaveFunctionPool& ppool):
    QMCDriver(w,psi,h,ppool), myWarmupSteps(0), UseDrift("yes"),    nSubSteps(1),
    myPeriod4WalkerDump(0)
  { 
    RootName = "vmc";
    QMCType ="VMCcuda";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    QMCDriverMode.set(QMC_WARMUP,0);
    m_param.add(UseDrift,"useDrift","string"); 
    m_param.add(UseDrift,"usedrift","string");
    m_param.add(myWarmupSteps,"warmupSteps","int");
    m_param.add(myWarmupSteps,"warmupsteps","int");
    m_param.add(nTargetSamples,"targetWalkers","int");
    m_param.add(nSubSteps, "substeps", "int");
    m_param.add(nSubSteps, "subSteps", "int");
  }
  
  bool VMCcuda::checkBounds (vector<PosType> &newpos,
			     vector<bool> &valid)
  {
    for (int iw=0; iw<newpos.size(); iw++) {
      PosType red = W.Lattice.toUnit(newpos[iw]);
      valid[iw] = W.Lattice.isValid(red);
    }
  }

  bool VMCcuda::run() { 
    if (UseDrift == "yes")
      return runWithDrift();

    resetRun();
    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    IndexType updatePeriod= (QMCDriverMode[QMC_UPDATE_MODE]) 
      ? Period4CheckProperties 
      : (nBlocks+1)*nSteps;
    
    int nat = W.getTotalNum();
    int nw  = W.getActiveWalkers();
    
    vector<RealType>  LocalEnergy(nw);
    vector<PosType>   delpos(nw);
    vector<PosType>   newpos(nw);
    vector<ValueType> ratios(nw);
    vector<GradType>  oldG(nw), newG(nw);
    vector<ValueType> oldL(nw), newL(nw);
    vector<Walker_t*> accepted(nw);
    Matrix<ValueType> lapl(nw, nat);
    Matrix<GradType>  grad(nw, nat);
    double Esum;

    // First do warmup steps
    for (int step=0; step<myWarmupSteps; step++) {
      for(int iat=0; iat<nat; ++iat)  {
	//create a 3N-Dimensional Gaussian with variance=1
	makeGaussRandomWithEngine(delpos,Random);
	for(int iw=0; iw<nw; ++iw) {
	  PosType G = W[iw]->G[iat];
	  newpos[iw]=W[iw]->R[iat] + m_sqrttau*delpos[iw];
	  ratios[iw] = 1.0;
	}
	W.proposeMove_GPU(newpos, iat);
	    
	Psi.ratio(W,iat,ratios,newG, newL);
	    
	accepted.clear();

	vector<bool> acc(nw, true);
	if (W.UseBoundBox)
	  checkBounds (newpos, acc);

	for(int iw=0; iw<nw; ++iw) {
	  if(acc[iw] && ratios[iw]*ratios[iw] > Random()) {
	    accepted.push_back(W[iw]);
	    nAccept++;
	    W[iw]->R[iat] = newpos[iw];
	    acc[iw] = true;
	  }
	  else {
	    acc[iw] = false;
	    nReject++;
	  }
	}
	W.acceptMove_GPU(acc);
	if (accepted.size())
	  Psi.update(accepted,iat);
      }
    }

    do {
      IndexType step = 0;
      nAccept = nReject = 0;
      Esum = 0.0;
      Estimators->startBlock(nSteps);
      do
      {
        ++step;++CurrentStep;
	for (int isub=0; isub<nSubSteps; isub++) {
	  for(int iat=0; iat<nat; ++iat)  {
	    //create a 3N-Dimensional Gaussian with variance=1
	    makeGaussRandomWithEngine(delpos,Random);
	    for(int iw=0; iw<nw; ++iw) {
	      PosType G = W[iw]->G[iat];
	      newpos[iw]=W[iw]->R[iat] + m_sqrttau*delpos[iw];
	      ratios[iw] = 1.0;
	    }
	    W.proposeMove_GPU(newpos, iat);
	    
	    Psi.ratio(W,iat,ratios,newG, newL);
	    
	    accepted.clear();

	    vector<bool> acc(nw, true);
	    if (W.UseBoundBox)
	      checkBounds (newpos, acc);

	    for(int iw=0; iw<nw; ++iw) {
	      if(acc[iw] && ratios[iw]*ratios[iw] > Random()) {
		accepted.push_back(W[iw]);
		nAccept++;
		W[iw]->R[iat] = newpos[iw];
		acc[iw] = true;
	      }
	      else {
		acc[iw]=false;
		nReject++;
	      }
	    }
	    W.acceptMove_GPU(acc);
	    if (accepted.size())
	      Psi.update(accepted,iat);
	  }
	}
	Psi.gradLapl(W, grad, lapl);
	H.evaluate (W, LocalEnergy);
	if (myPeriod4WalkerDump && (CurrentStep % myPeriod4WalkerDump)==0) 
	  W.saveEnsemble();
	Estimators->accumulate(W);
      } while(step<nSteps);
      Psi.recompute(W);

      // vector<RealType> logPsi(W.WalkerList.size(), 0.0);
      // Psi.evaluateLog(W, logPsi);
      
      double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
      Estimators->stopBlock(accept_ratio);

      nAcceptTot += nAccept;
      nRejectTot += nReject;
      ++block;
      recordBlock(block);
    } while(block<nBlocks);

    //Mover->stopRun();

    //finalize a qmc section
    return finalize(block);
  }



  bool VMCcuda::runWithDrift() 
  { 
    resetRun();
    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    int nat = W.getTotalNum();
    int nw  = W.getActiveWalkers();
    
    vector<RealType>  LocalEnergy(nw), oldScale(nw), newScale(nw);
    vector<PosType>   delpos(nw);
    vector<PosType>   dr(nw);
    vector<PosType>   newpos(nw);
    vector<ValueType> ratios(nw), rplus(nw), rminus(nw);
    vector<PosType>  oldG(nw), newG(nw);
    vector<ValueType> oldL(nw), newL(nw);
    vector<Walker_t*> accepted(nw);
    Matrix<ValueType> lapl(nw, nat);
    Matrix<GradType>  grad(nw, nat);

    // First, do warmup steps
    for (int step=0; step<myWarmupSteps; step++) {
      for(int iat=0; iat<nat; iat++) {
	Psi.getGradient (W, iat, oldG);
	    
	//create a 3N-Dimensional Gaussian with variance=1
	makeGaussRandomWithEngine(delpos,Random);
	for(int iw=0; iw<nw; iw++) {
	  oldScale[iw] = getDriftScale(m_tauovermass,oldG[iw]);
	  dr[iw] = (m_sqrttau*delpos[iw]) + (oldScale[iw]*oldG[iw]);
	  newpos[iw]=W[iw]->R[iat] + dr[iw];
	  ratios[iw] = 1.0;
	}
	W.proposeMove_GPU(newpos, iat);
	    
	Psi.ratio(W,iat,ratios,newG, newL);
	    
	accepted.clear();
	vector<bool> acc(nw, true);
	if (W.UseBoundBox)
	  checkBounds (newpos, acc);

	for(int iw=0; iw<nw; ++iw) {
	  PosType drOld = 
	    newpos[iw] - (W[iw]->R[iat] + oldScale[iw]*oldG[iw]);
	  RealType logGf = -m_oneover2tau * dot(drOld, drOld);
	  newScale[iw]   = getDriftScale(m_tauovermass,newG[iw]);
	  PosType drNew  = 
	    (newpos[iw] + newScale[iw]*newG[iw]) - W[iw]->R[iat];

	  RealType logGb =  -m_oneover2tau * dot(drNew, drNew);
	  RealType x = logGb - logGf;
	  RealType prob = ratios[iw]*ratios[iw]*std::exp(x);

	  if(acc[iw] && Random() < prob) {
	    accepted.push_back(W[iw]);
	    nAccept++;
	    W[iw]->R[iat] = newpos[iw];
	    acc[iw] = true;
	  }
	  else {
	    acc[iw] = false;
	    nReject++;
	  }
	}
	W.acceptMove_GPU(acc);
	if (accepted.size())
	  Psi.update(accepted,iat);
      }
    }

    // Now do data collection steps
    do {
      IndexType step = 0;
      nAccept = nReject = 0;
      Estimators->startBlock(nSteps);
      do {
        step++;
	CurrentStep++;
	for (int isub=0; isub<nSubSteps; isub++) {
	  for(int iat=0; iat<nat; iat++) {
	    Psi.getGradient (W, iat, oldG);
	    
	    //create a 3N-Dimensional Gaussian with variance=1
	    makeGaussRandomWithEngine(delpos,Random);
	    for(int iw=0; iw<nw; iw++) {
	      oldScale[iw] = getDriftScale(m_tauovermass,oldG[iw]);
	      dr[iw] = (m_sqrttau*delpos[iw]) + (oldScale[iw]*oldG[iw]);
	      newpos[iw]=W[iw]->R[iat] + dr[iw];
	      ratios[iw] = 1.0;
	    }
	    W.proposeMove_GPU(newpos, iat);
	    
	    Psi.ratio(W,iat,ratios,newG, newL);
	    
	    accepted.clear();

	    vector<bool> acc(nw, true);
	    if (W.UseBoundBox)
	      checkBounds (newpos, acc);

	    for(int iw=0; iw<nw; ++iw) {
	      PosType drOld = 
		newpos[iw] - (W[iw]->R[iat] + oldScale[iw]*oldG[iw]);
	      // if (dot(drOld, drOld) > 25.0)
	      //   cerr << "Large drift encountered!  Old drift = " << drOld << endl;
	      RealType logGf = -m_oneover2tau * dot(drOld, drOld);
	      newScale[iw]   = getDriftScale(m_tauovermass,newG[iw]);
	      PosType drNew  = 
		(newpos[iw] + newScale[iw]*newG[iw]) - W[iw]->R[iat];
	      // if (dot(drNew, drNew) > 25.0)
	      //   cerr << "Large drift encountered!  Drift = " << drNew << endl;
	      RealType logGb =  -m_oneover2tau * dot(drNew, drNew);
	      RealType x = logGb - logGf;
	      RealType prob = ratios[iw]*ratios[iw]*std::exp(x);
	      
	      if(acc[iw] && Random() < prob) {
		accepted.push_back(W[iw]);
		nAccept++;
		W[iw]->R[iat] = newpos[iw];
		acc[iw] = true;
	      }
	      else {
		acc[iw] = false;
		nReject++;
	      }
	    }
	    W.acceptMove_GPU(acc);
	    if (accepted.size())
	      Psi.update(accepted,iat);
	  }
	  // cerr << "Rank = " << myComm->rank() <<
	  //   "  CurrentStep = " << CurrentStep << "  isub = " << isub << endl;
	}
	Psi.gradLapl(W, grad, lapl);
	H.evaluate (W, LocalEnergy);
	if (myPeriod4WalkerDump && (CurrentStep % myPeriod4WalkerDump)==0) 
	  W.saveEnsemble();
	Estimators->accumulate(W);
      } while(step<nSteps);
      Psi.recompute(W);
      
      double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
      Estimators->stopBlock(accept_ratio);

      nAcceptTot += nAccept;
      nRejectTot += nReject;
      ++block;
      recordBlock(block);
    } while(block<nBlocks);
    //finalize a qmc section
    if (!myComm->rank())
      gpu::cuda_memory_manager.report();
    return finalize(block);
  }





  void VMCcuda::resetRun()
  {
    SpeciesSet tspecies(W.getSpeciesSet());
    int massind=tspecies.addAttribute("mass");
    RealType mass = tspecies(massind,0);
    RealType oneovermass = 1.0/mass;
    RealType oneoversqrtmass = std::sqrt(oneovermass);
    m_oneover2tau = 0.5*mass/Tau;
    m_sqrttau = std::sqrt(Tau/mass);
    m_tauovermass = Tau/mass;

    // Compute the size of data needed for each walker on the GPU card
    PointerPool<Walker_t::cuda_Buffer_t > pool;
    
    Psi.reserve (pool);
    app_log() << "Each walker requires " << pool.getTotalSize() * sizeof(CudaRealType)
	      << " bytes in GPU memory.\n";

    // Now allocate memory on the GPU card for each walker
    // for (int iw=0; iw<W.WalkerList.size(); iw++) {
    //   Walker_t &walker = *(W.WalkerList[iw]);
    //   walker.resizeCuda(pool.getTotalSize());
    //   // pool.allocate(walker.cuda_DataSet);
    // }
    W.allocateGPU(pool.getTotalSize());
    app_log() << "Successfully allocated walkers.\n";
    W.copyWalkersToGPU();
    W.updateLists_GPU();
    vector<RealType> logPsi(W.WalkerList.size(), 0.0);
    //Psi.evaluateLog(W, logPsi);
    Psi.recompute(W, true);
    Estimators->start(nBlocks, true);

    // Compute sample dumping frequency
    if (nTargetSamples) {
      int nnodes = myComm->size();
      int nw = W.WalkerList.size();
      int samples_per_node = (nTargetSamples+nnodes-1)/nnodes; 
      int dumps_per_node   = (samples_per_node+nw-1) / nw;
      myPeriod4WalkerDump = Period4WalkerDump;
      app_log() << "  Dumping walker ensemble every " << myPeriod4WalkerDump
		<< " steps.\n";
    }

    W.clearEnsemble();
    int samples_this_node = nTargetSamples/myComm->size();
    if (nTargetSamples%myComm->size() > myComm->rank()) samples_this_node+=1;
    app_log() << "  Node zero will generate " << samples_this_node << " samples.\n";
    W.setNumSamples(samples_this_node);
  }

  bool 
  VMCcuda::put(xmlNodePtr q){
    //nothing to add
    return true;
  }
}

/***************************************************************************
 * $RCSfile: VMCParticleByParticle.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCParticleByParticle.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
