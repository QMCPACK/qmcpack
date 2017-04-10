//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

/**@file CSVMC.h
 * @brief Definition of CSVMC
 */
#ifndef QMCPLUSPLUS_CS_VMCMULTIPLE_H
#define QMCPLUSPLUS_CS_VMCMULTIPLE_H
#include "QMCDrivers/QMCDriver.h"
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
class CSVMC: public QMCDriver
{
public:
  /// Constructor.
  CSVMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool);

  bool run();
  bool put(xmlNodePtr cur);

private:
  ///blocks over which normalization factors are accumulated
  int equilBlocks;
  /// Copy Constructor (disabled)
  CSVMC(const CSVMC& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  CSVMC& operator=(const CSVMC&)
  {
    return *this;
  }

  CSEnergyEstimator *multiEstimator;
  CSUpdateBase* Mover;
};
}

#endif
