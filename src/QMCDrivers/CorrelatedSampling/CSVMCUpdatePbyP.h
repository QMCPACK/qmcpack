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
    
    


#ifndef QMCPLUSPLUS_CS_VMC_UPDATEPBYP_H
#define QMCPLUSPLUS_CS_VMC_UPDATEPBYP_H
#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers MultiplePsi ParticleByParticle
 * @brief Implements the VMC algorithm
 */
class CSVMCUpdatePbyP: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdatePbyP(MCWalkerConfiguration& w,
                   std::vector<TrialWaveFunction*>& psi,
                   std::vector<QMCHamiltonian*>& h,
                  RandomGenerator_t& rg);

  ~CSVMCUpdatePbyP();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
};
/** @ingroup QMCDrivers MultiplePsi ParticleByParticle
 *
 * @brief Implements the VMC algorithm with drift.  
 */
class CSVMCUpdatePbyPWithDriftFast: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w,
                   std::vector<TrialWaveFunction*>& psi,
                   std::vector<QMCHamiltonian*>& h,
                  RandomGenerator_t& rg);

  ~CSVMCUpdatePbyPWithDriftFast();

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
};
}

#endif
