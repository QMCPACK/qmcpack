//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  //inline std::vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  ///number of warmup steps
  int myWarmupSteps;
  ///number of RN warmup steps
  int myRNWarmupSteps;
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///option to enable/disable drift equation or RN for VMC
  std::string UseDrift;
  ///order or renyi to compute
  int EEN;
  ///patch geometry (sphere,...
  std::string computeEE;
  ///patch size (sphere:r^2,...
  RealType vsize;
  ///index of particles in region a and b for threads
  std::vector<std::vector<int> > region;
  std::vector<RealType> stats_vec;

  bool plotSwapAmplitude;
  int grid_spacing;

  fstream file_out;
  int Nmin,Nmax;
  std::stringstream ee_dat,pe_dat;
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
