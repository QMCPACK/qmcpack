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
#include "QMCDrivers/DMC/DMC.h" 
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/DMC/DMCOMP.h" 
#endif
#include "Message/OpenMP.h"

namespace qmcplusplus {

  QMCDriver* DMCFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomNumberControl& rc, HamiltonianPool& hpool) 
  {
    QMCDriver* qmc=0;
#if defined(ENABLE_OPENMP)
    int np=omp_get_max_threads();
    if(np>1) {
      app_log() << "Creating DMCMP for the qmc driver" << endl;
      qmc = new DMCOMP(w,psi,h,rc,hpool);
    } else {
      app_log() << "Creating DMC for the qmc driver" << endl;
      qmc = new DMC(w,psi,h,rc); 
    }
#else
    app_log() << "Creating DMC for the qmc driver" << endl;
    qmc = new DMC(w,psi,h,rc); 
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
