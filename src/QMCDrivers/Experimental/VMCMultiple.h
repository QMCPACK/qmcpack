//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file VMCMultiple.h
 * @brief Definition of VMCMultiple
 */
#ifndef QMCPLUSPLUS_VMCMULTIPLE_H_PSI_H
#define QMCPLUSPLUS_VMCMULTIPLE_H_PSI_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class MultipleEnergyEstimator;

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class VMCMultiple: public QMCDriver
{
public:
  /// Constructor.
  VMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  void advanceWalkerByWalker();
  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMCMultiple(const VMCMultiple& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCMultiple& operator=(const VMCMultiple&)
  {
    return *this;
  }

  MultipleEnergyEstimator *multiEstimator;
  ///temporary storage
  int nPsi;
  ///number of blocks to compute the normalization factor
  int equilBlocks;
  std::vector<RealType> logpsi;
  std::vector<RealType> sumratio;
  std::vector<RealType> invsumratio;
  std::vector<RealType> Norm;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
