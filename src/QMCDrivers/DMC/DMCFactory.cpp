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
#include "QMCDrivers/DMC/DMCFactory.h" 
#include "QMCDrivers/DMC/DMCMoveAll.h" 
#include "QMCDrivers/DMC/DMCPbyP.h" 
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/DMC/DMCPbyPOMP.h" 
#endif
#include "Message/OpenMP.h"

namespace qmcplusplus {

  QMCDriver* DMCFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, HamiltonianPool& hpool) {
    
    //determine the type of move
    bool pbyp=(DMCType == "dmc-ptcl");
   // string reconfigure("no");
   // ParameterSet param;
   // param.add(reconfigure,"reconfiguration");
   // param.put(cur);
   // bool fixW = (reconfigure == "yes");
   QMCDriver* qmc=0;
    if(pbyp) {
#if defined(ENABLE_OPENMP)
      int np=omp_get_max_threads();
      if(np>1) {
        app_log() << "Creating DMCPbyPOpenMP for the qmc driver" << endl;
        DMCPbyPOMP* a = new DMCPbyPOMP(w,psi,h); 
        a->makeClones(hpool,np);
        qmc=a;
      } else {
        app_log() << "Creating DMCPbyP for the qmc driver" << endl;
        qmc = new DMCPbyP(w,psi,h); 
      }
#else
      app_log() << "Creating DMCPbyP for the qmc driver" << endl;
      qmc = new DMCPbyP(w,psi,h); 
#endif

    } else {
       qmc = new DMCMoveAll(w,psi,h); 
    }

    return qmc;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
