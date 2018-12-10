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
    
    


/**@file CSVMCUpdateAll.h
 * @brief Definition of CSVMCUpdateAll
 */
#ifndef QMCPLUSPLUS_CS_VMC_UPDATEALL_H
#define QMCPLUSPLUS_CS_VMC_UPDATEALL_H

#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class CSVMCUpdateAll: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdateAll(MCWalkerConfiguration& w,  std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& h,
                 RandomGenerator_t& rg);

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
};


class CSVMCUpdateAllWithDrift: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdateAllWithDrift(MCWalkerConfiguration& w,  std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& h,
                 RandomGenerator_t& rg);

  void advanceWalker(Walker_t& thisWalker, bool recompute);

private:
};
}

#endif
