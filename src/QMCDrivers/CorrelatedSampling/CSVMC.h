//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/**@file CSVMC.h
 * @brief Definition of CSVMC
 */
#ifndef QMCPLUSPLUS_CS_VMCMULTIPLE_H
#define QMCPLUSPLUS_CS_VMCMULTIPLE_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

class CSEnergyEstimator;
class CSUpdateBase;

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class CSVMC: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  CSVMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, 
            HamiltonianPool& hpool, WaveFunctionPool& ppool);

  bool run();
  bool put(xmlNodePtr cur);

private:
  std::string UseDrift;
  int prevSteps;
  int prevStepsBetweenSamples;
  ///blocks over which normalization factors are accumulated
  int equilBlocks;
  /// Copy Constructor (disabled)
  CSVMC(const CSVMC& a): QMCDriver(a), CloneManager(a) { }
  /// Copy operator (disabled).
  CSVMC& operator=(const CSVMC&)
  {
    return *this;
  }
  
  void resetRun();
  

  CSEnergyEstimator *multiEstimator;
  CSUpdateBase* Mover;
};
}

#endif
