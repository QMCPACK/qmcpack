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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_VMCRENYI_OMP_H
#define QMCPLUSPLUS_VMCRENYI_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCRenyiOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  VMCRenyiOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
              HamiltonianPool& hpool, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
  //inline vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  ///number of warmup steps
  int myWarmupSteps;
  ///number of RN warmup steps
  int myRNWarmupSteps;
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///option to enable/disable drift equation or RN for VMC
  string UseDrift;
  ///order or renyi to compute
  int EEN;
  ///patch geometry (sphere,...
  string computeEE;
  ///patch size (sphere:r^2,...
  RealType vsize;
  ///index of particles in region a and b for threads
  vector<vector<int> > region;
  vector<RealType> stats_vec;

  bool plotSwapAmplitude;
  int grid_spacing;

  fstream file_out;
  int Nmin,Nmax;
  stringstream ee_dat,pe_dat;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  VMCRenyiOMP(const VMCRenyiOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  VMCRenyiOMP& operator=(const VMCRenyiOMP&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile: VMCRenyiOMP.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCRenyiOMP.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
