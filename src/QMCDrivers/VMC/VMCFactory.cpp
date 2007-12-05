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
#include "QMCDrivers/VMC/VMCFactory.h" 
#include "QMCDrivers/VMC/VMCSingle.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#endif
//#include "QMCDrivers/VMC/VMCMultiple.h"
//#include "QMCDrivers/VMC/VMCPbyPMultiple.h"
#if !defined(QMC_COMPLEX)
#include "QMCDrivers/VMC/VMCMultipleWarp.h"
#include "QMCDrivers/VMC/VMCPbyPMultiWarp.h"
#endif
#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#include "Message/OpenMP.h"

namespace qmcplusplus {

  QMCDriver* VMCFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h,RandomNumberControl& rc, 
      ParticleSetPool& ptclpool, HamiltonianPool& hpool) 
  {
    int np=omp_get_max_threads();

    //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
    QMCDriver* qmc=0;
    if(VMCMode == 0 || VMCMode == 1) //(0,0,0) (0,0,1)
    {
#if defined(ENABLE_OPENMP)
      if(np>1)
        qmc = new VMCSingleOMP(w,psi,h,rc,hpool);
      else
#endif
        qmc = new VMCSingle(w,psi,h,rc);
    } 
    //else if(VMCMode == 2) //(0,1,0)
    //{
    //  qmc = new VMCMultiple(w,psi,h);
    //} 
    //else if(VMCMode == 3) //(0,1,1)
    //{
    //  qmc = new VMCPbyPMultiple(w,psi,h);
    //} 
    else if(VMCMode ==2 || VMCMode ==3)
    {
      qmc = new CSVMC(w,psi,h,rc);
    }
#if !defined(QMC_COMPLEX)
    else if(VMCMode == 6) //(1,1,0)
    {
      qmc = new VMCMultipleWarp(w,psi,h, rc,ptclpool);
    } 
    else if(VMCMode == 7) //(1,1,1)
    {
      qmc = new VMCPbyPMultiWarp(w,psi,h,rc, ptclpool);
    }
#endif
    qmc->setUpdateMode(VMCMode&1);
    return qmc;
  }
}
/***************************************************************************
 * $RCSfile: DMCFactory.cpp,v $   $Author: jnkim $
 * $Revision: 1.3 $   $Date: 2006/04/05 00:49:59 $
 * $Id: DMCFactory.cpp,v 1.3 2006/04/05 00:49:59 jnkim Exp $ 
 ***************************************************************************/
