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
    
    
#ifndef QMCPLUSPLUS_VMC_WITHUPDATEENINGE_H
#define QMCPLUSPLUS_VMC_WITHUPDATEENINGE_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCSingle: public QMCDriver
{
public:
  /// Constructor.
  VMCSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);

private:
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///update engine
  QMCUpdateBase* Mover;
  ///option to enable/disable drift term
  std::string UseDrift;
  /// Copy Constructor (disabled)
  VMCSingle(const VMCSingle& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCSingle& operator=(const VMCSingle&)
  {
    return *this;
  }
  ///hide initialization from the main function
  void resetRun();
};
}

#endif
