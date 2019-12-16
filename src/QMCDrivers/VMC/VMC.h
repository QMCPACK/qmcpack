//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_VMC_H
#define QMCPLUSPLUS_VMC_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMC : public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  VMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool, Communicate* comm);
  bool run();
  bool put(xmlNodePtr cur);
  QMCRunType getRunType() { return QMCRunType::VMC; }
  //inline std::vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  int prevSteps;
  int prevStepsBetweenSamples;

  ///Ways to set rn constant
  RealType logoffset, logepsilon;
  ///option to enable/disable drift equation or RN for VMC
  std::string UseDrift;
  ///turn on spin moves
  std::string SpinMoves;
  RealType SpinMass;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  VMC(const VMC&) = delete;
  /// Copy operator (disabled).
  VMC& operator=(const VMC&) = delete;
};
} // namespace qmcplusplus

#endif
