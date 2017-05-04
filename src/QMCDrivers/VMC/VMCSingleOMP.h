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
    
    


#ifndef QMCPLUSPLUS_VMCSINGLE_OMP_H
#define QMCPLUSPLUS_VMCSINGLE_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

class ParticleSet;

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCSingleOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  VMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
               HamiltonianPool& hpool, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
  //inline std::vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  int prevSteps;
  int prevStepsBetweenSamples;

  ///Ways to set rn constant
  RealType logoffset,logepsilon;
  ///coefficient for the laplacian term in the nodeless guiding function formula
  RealType NodelessEpsilon;
//  ///coefficient for the exponential damping in the nodeless guiding function formula
//  RealType NodelessAlpha;
//  RealType NodelessBeta;
  ///option to enable/disable drift equation or RN for VMC
  std::string UseDrift;
  ///pointer to the ion particle set needed by nodeless guiding
  ParticleSet * ionSetPtr;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  VMCSingleOMP(const VMCSingleOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  VMCSingleOMP& operator=(const VMCSingleOMP&)
  {
    return *this;
  }
};
}

#endif
