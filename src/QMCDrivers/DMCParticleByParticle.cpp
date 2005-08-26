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
    KillNodeCrossing(0),
    PopIndex(-1), EtrialIndex(-1),
    BranchInfo("default"), KillWalker("no"){ 
    RootName = "dmc";
    QMCType ="dmc";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    m_param.add(KillWalker,"killnode","string");
  }
  
  /// destructor
  DMCParticleByParticle::~DMCParticleByParticle() {
  }

  bool DMCParticleByParticle::run() { 

    KillNodeCrossing = (KillWalker == "yes");
    if(KillNodeCrossing) {
      LOGMSG("Walkers will be killed if a node crossing is detected.")
    } else {
      LOGMSG("Walkers will be kept even if a node crossing is detected.")
    }

    //add columns
    IndexType PopIndex = Estimators->addColumn("Population");
    IndexType EtrialIndex = Estimators->addColumn("Etrial");
    //write the header
    Estimators->reportHeader();

    //initialization of branchEngine has moved to QMCDriver
    //if(branchEngine == 0) {
    //  branchEngine=new BranchEngineType(Tau,W.getActiveWalkers());
    //  RealType e_ref = W.getLocalEnergy();
    //  branchEngine->setEguess(e_ref);
    //}
    //branchEngine->put(qmcNode,LogOut);
    ////this should be controlled by continue="yes/no" of qmc element
    //branchEngine->flush(0);
    //branchEngine->put(qmcNode,LogOut);

    /* This has moved to QMCDriver::initialize()
    //going to add routines to calculate how much we need
    bool require_register =  W.createAuxDataSet();
    //check if buffer is properly initialized
    if(require_register) {
      MCWalkerConfiguration::iterator it(W.begin()); 
      MCWalkerConfiguration::iterator it_end(W.end()); 
      while(it != it_end) {
        (*it)->Weight= 1.0;
        (*it)->Multiplicity=1.0;
        (*it)->DataSet.rewind();
        W.registerData(**it,(*it)->DataSet);
        Psi.registerData(W,(*it)->DataSet);
        ++it;
      } 
    } else {
      updateWalkers();
    }
    */

    Estimators->reset();
    
    IndexType block = 0;
    
    Pooma::Clock timer;
    RealType Eest = branchEngine->E_T;
    IndexType nat = W.getTotalNum();
    G.resize(nat);
    dG.resize(nat);
    L.resize(nat);
    dL.resize(nat);

    IndexType accstep=0;
    nAcceptTot = 0;
    nRejectTot = 0;
    m_oneover2tau = 0.5/Tau;
    m_sqrttau = sqrt(Tau);
    do {
      IndexType step = 0;
      timer.start();
      nAccept = 0; 
      nReject=0;
      nAllRejected = 0;
      nNodeCrossing=0;
      IndexType pop_acc=0; 
      do {
        pop_acc += W.getActiveWalkers();
        //default is killing
        if(KillNodeCrossing) 
          advanceKillNodeCrossing(nat);
        else
          advanceRejectNodeCrossing(nat);
        ++step; ++accstep;
        Estimators->accumulate(W);
        branchEngine->branch(accstep,W);
        Eest = branchEngine->update(W.getActiveWalkers(), Eest); 
        //every 100 step, we update the walkers
        if(accstep%100 == 0) updateWalkers();
      } while(step<nSteps);
      
      //WARNMSG("The number of a complete rejectoin " << nAllRejected)
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      //update estimator
      Estimators->flush();
      Estimators->setColumn(PopIndex,static_cast<RealType>(pop_acc)/static_cast<RealType>(nSteps));
      Estimators->setColumn(EtrialIndex,Eest); 
      Estimators->setColumn(AcceptIndex,
      		     static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
      Estimators->report(accstep);
      
      Eest = Estimators->average(0);
      RealType totmoves=1.0/static_cast<RealType>(step*W.getActiveWalkers());
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
      		    << static_cast<RealType>(nAllRejected)*totmoves << " Node crossing "
                    << static_cast<RealType>(nNodeCrossing)*totmoves << endl;

      if(pStride) { //create an output engine
        HDFWalkerOutput WO(RootName);
        WO.get(W); 
        WO.write(*branchEngine);
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
      WO.write(*branchEngine);
    }
    
    Estimators->finalize();
    
    return true;
  }

  /** advance all the walkers with killnode==yes
   * @param nat number of particles to move
   * 
   * When killnode==yes, any move resulting in node-crossing will lead to
   * the death of the walker.
   */
  void DMCParticleByParticle::advanceKillNodeCrossing(int nat) {

    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
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
      RealType rr_proposed=0.0;
      RealType rr_accepted=0.0;
      while(notcrossed && iat<nat){

        PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        PosType newpos(W.makeMove(iat,dr));
        RealType ratio(Psi.ratio(W,iat,dG,dL));

        if(ratio < 0.0) {//node is crossed, stop here
          notcrossed = false;
          Psi.restore(iat);
          ++nNodeCrossing;
        } else {
          RealType rr=dot(dr,dr);
          rr_proposed+=rr;
          G = W.G+dG;
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
          ValueType scale=Tau;
          //ValueType vsq = dot(G[iat],G[iat]);
          //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
          dr = thisWalker.R[iat]-newpos-scale*G[iat]; 
          RealType logGb = -m_oneover2tau*dot(dr,dr);
          RealType prob = std::min(1.0,ratio*ratio*exp(logGb-logGf));
          if(Random() < prob) { 
            ++nAcceptTemp;
            W.acceptMove(iat);
            Psi.update2(W,iat);
            W.G = G;
            W.L += dL;
            thisWalker.Drift = scale*G;
            rr_accepted+=rr;
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
          emixed = (eold+enew)*0.5;
          eold=enew;
        } else {
          thisWalker.Age++;
          ++nAllRejected;
        }

        RealType tau_eff=Tau;
        ValueType M = branchEngine->branchGF(tau_eff,emixed,0.0);
        if(thisWalker.Age > 3) M = std::min(0.5,M);
        else if(thisWalker.Age > 0) M = std::min(1.0,M);
        thisWalker.Weight = M; 
        thisWalker.Multiplicity=M + Random();

        //accumulate the energy
        nAccept += nAcceptTemp;
        nReject += nRejectTemp;
      } else {//set the weight and multiplicity to zero
        thisWalker.willDie();
        nReject += nat;
      }
      branchEngine->accumulate(eold,thisWalker.Weight);
      ++it; 
    }
  }

  /** advance all the walkers with killnode==no
   * @param nat number of particles to move
   * 
   * When killnode==no, any move resulting in node-crossing is treated
   * as a normal rejection.
   */
  void DMCParticleByParticle::advanceRejectNodeCrossing(int nat) {
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    int  iwalker=0; 
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

      RealType rr_proposed=0.0;
      RealType rr_accepted=0.0;
      while(iat<nat) {//particle-by-particle move
        PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        PosType newpos(W.makeMove(iat,dr));

        RealType ratio=Psi.ratio(W,iat,dG,dL);

        RealType rr=dot(dr,dr);
        rr_proposed+=rr;

        if(ratio < 0.0) {//node is crossed reject the move
          ++nRejectTemp;
          ++nNodeCrossing;
          Psi.restore(iat);
        } else {
          G = W.G+dG;
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
          RealType scale=Tau;
          //ValueType vsq = Dot(G,G);
          //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
          dr = thisWalker.R[iat]-newpos-scale*G[iat]; 
          RealType logGb = -m_oneover2tau*dot(dr,dr);
          RealType prob = std::min(1.0,ratio*ratio*exp(logGb-logGf));
          if(Random() < prob) { 
            ++nAcceptTemp;
            W.acceptMove(iat);
            Psi.update2(W,iat);
            W.G = G;
            W.L += dL;
            thisWalker.Drift = scale*G;
            rr_accepted+=rr;
          } else {
            ++nRejectTemp; 
            Psi.restore(iat);
          }
        } 

        ++iat;
      }

      if(nAcceptTemp>0) {//need to overwrite the walker properties
        thisWalker.R = W.R;
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        ValueType psi = Psi.evaluate(W,w_buffer);
        enew= H.evaluate(W);
        thisWalker.resetProperty(log(abs(psi)),psi,enew);
        H.saveProperty(thisWalker.getPropertyBase());
        emixed = (eold+enew)*0.5e0;
      } else {
        thisWalker.Age++;
        ++nAllRejected;
        rr_accepted=0.0;
      }

      ValueType M = branchEngine->branchGF(Tau*rr_accepted/rr_proposed,emixed,0.0);
      if(thisWalker.Age > 1) M = std::min(0.5,M);
      else if(thisWalker.Age > 0) M = std::min(1.0,M);
      thisWalker.Weight = M; 
      thisWalker.Multiplicity=M + Random();
      branchEngine->accumulate(emixed,thisWalker.Weight);//accumulate the energy
      nAccept += nAcceptTemp;
      nReject += nRejectTemp;
      ++it;
    }
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
