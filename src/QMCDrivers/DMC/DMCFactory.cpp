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
#include "QMCDrivers/DMC/DMCFactory.h" 
#include "QMCDrivers/DMC/DMC.h" 
#include "QMCDrivers/DMC/DMCPeta.h" 
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/DMC/DMCOMP.h" 
#endif
#include "Message/OpenMP.h"

//#define PETA_DMC_TEST

namespace qmcplusplus {

  QMCDriver* DMCFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, HamiltonianPool& hpool) {
    QMCDriver* qmc=0;
#if defined(PETA_DMC_TEST)
    qmc = new DMCPeta(w,psi,h);
#else

#if defined(ENABLE_OPENMP)
    int np=omp_get_max_threads();
    if(np>1) {
      app_log() << "Creating DMCMP for the qmc driver" << endl;
      qmc = new DMCOMP(w,psi,h,hpool);
    } else {
      app_log() << "Creating DMC for the qmc driver" << endl;
      qmc = new DMC(w,psi,h); 
    }
#else
    app_log() << "Creating DMC for the qmc driver" << endl;
    qmc = new DMC(w,psi,h); 
#endif

#endif
    qmc->setUpdateMode(PbyPUpdate);

    return qmc;
  }
}
/***************************************************************************
 * $RCSfile: DMCFactory.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
