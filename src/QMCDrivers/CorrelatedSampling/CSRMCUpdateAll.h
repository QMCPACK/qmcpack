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
  CSVMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                 RandomGenerator_t& rg);

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
};
}

#endif
