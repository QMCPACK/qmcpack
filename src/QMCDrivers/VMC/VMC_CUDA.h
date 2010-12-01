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
#ifndef QMCPLUSPLUS_VMC_CUDA_H
#define QMCPLUSPLUS_VMC_CUDA_H
#include "QMCDrivers/QMCDriver.h" 
namespace qmcplusplus {

  class QMCUpdateBase;

  /** @ingroup QMCDrivers  PbyP
   *@brief Implements the VMC algorithm using particle-by-particle move. 
   */
  class VMCcuda: public QMCDriver {
  public:
    /// Constructor.
    VMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool);
    bool run();
    bool runWithDrift();
    bool put(xmlNodePtr cur);
 
  private:
    /// tau/mass
    RealType m_tauovermass;
    /// Whether to use drift or not
    string UseDrift;
    ///number of warmup steps
    int myWarmupSteps, nSubSteps;
    ///period for walker dump
    int myPeriod4WalkerDump;
    /// Copy Constructor (disabled)
    VMCcuda(const VMCcuda& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    VMCcuda& operator=(const VMCcuda&) { return *this;}
    ///hide initialization from the main function
    bool checkBounds (vector<PosType> &newpos, vector<bool> &valid);

    void resetRun();
  };
}

#endif
/***************************************************************************
 * $RCSfile: VMCcuda.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCcuda.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $ 
 ***************************************************************************/
