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
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#include "QMCDrivers/QMCDriver.h" 
#include "QMCDrivers/DMCPbyPUpdate.h" 
namespace qmcplusplus {


  class HamiltonianPool;

  /** @ingroup QMCDrivers 
   *@brief A dummy QMCDriver for testing
   */
  class DMCPbyPOpenMP: public QMCDriver {
  public:

    /// Constructor.
    DMCPbyPOpenMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

    void makeClones(HamiltonianPool& hpool, int np=-1);

    bool run();
    bool put(xmlNodePtr cur);
 
  private:
    /// Copy Constructor (disabled)
    DMCPbyPOpenMP(const DMCPbyPOpenMP& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    DMCPbyPOpenMP& operator=(const DMCPbyPOpenMP&) { return *this;}

    int NumThreads;

    vector<DMCPbyPUpdate*> Movers;
    vector<ParticleSet*> wClones;
    vector<TrialWaveFunction*> psiClones;
    vector<QMCHamiltonian*> hClones;
    vector<RandomGenerator_t*> Rng;
    vector<int> wPerNode;

    void resetRun();
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
