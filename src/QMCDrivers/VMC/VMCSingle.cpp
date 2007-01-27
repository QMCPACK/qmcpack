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
#include "QMCDrivers/VMC/VMCSingle.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"

namespace qmcplusplus { 

  /// Constructor.
  VMCSingle::VMCSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), Mover(0), UseDrift("yes") { 
    RootName = "vmc";
    QMCType ="VMPSingle";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    m_param.add(UseDrift,"useDrift","string"); m_param.add(UseDrift,"usedrift","string");
  }
  
  bool VMCSingle::run() { 
    if(Mover ==0)
    {
      if(QMCDriverMode[QMC_UPDATE_MODE])
      {
        if(UseDrift == "yes")
          Mover=new VMCUpdatePbyPWithDrift(W,Psi,H,Random);
        else
          Mover=new VMCUpdatePbyP(W,Psi,H,Random);
        Mover->resetRun(branchEngine);
        Mover->initWalkersForPbyP(W.begin(),W.end());
      }
      else
      {
        if(UseDrift == "yes")
          Mover=new VMCUpdateAllWithDrift(W,Psi,H,Random);
        else
          Mover=new VMCUpdateAll(W,Psi,H,Random);
        Mover->resetRun(branchEngine);
        Mover->initWalkers(W.begin(),W.end());
      }
    }

    Mover->put(qmcNode);
    Estimators->reportHeader(AppendRun);
    Estimators->reset();

    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;

    do {
      Estimators->startBlock();
      Mover->startBlock();
      IndexType step = 0;
      do
      {
        Mover->advanceWalkers(W.begin(),W.end());
        ++step;++CurrentStep;
        Estimators->accumulate(W);
      } while(step<nSteps);
      
      Estimators->stopBlock(Mover->acceptRatio());

      nAcceptTot += Mover->nAccept;
      nRejectTot += Mover->nReject;

      ++block;

      recordBlock(block);

      if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0) 
        Mover->updateWalkers(W.begin(),W.end());

    } while(block<nBlocks);

    //Need MPI-IO
    app_log() << "Ratio = " 
      << static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot)
      << endl;
    
    //finalize a qmc section
    return finalize(block);
  }

  bool 
  VMCSingle::put(xmlNodePtr q){
    //nothing to add
    return true;
  }
}

/***************************************************************************
 * $RCSfile: VMCParticleByParticle.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCParticleByParticle.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
