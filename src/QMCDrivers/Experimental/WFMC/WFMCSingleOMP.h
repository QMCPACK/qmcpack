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
    
    


#ifndef QMCPLUSPLUS_WFMCSINGLE_OMP_H
#define QMCPLUSPLUS_WFMCSINGLE_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class WFMCSingleOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  WFMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                HamiltonianPool& hpool,WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
private:
  ///number of warmup steps
  int myWarmupSteps;
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///Index for Energy Index
  int Eindex;
  ///option to enable/disable drift equation for VMC
  std::string UseDrift,reweight;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  WFMCSingleOMP(const WFMCSingleOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  WFMCSingleOMP& operator=(const WFMCSingleOMP&)
  {
    return *this;
  }
};
}

#endif
