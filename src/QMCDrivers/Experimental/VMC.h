//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_VMC_H
#define QMCPLUSPLUS_VMC_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers WalkerByWalker
 *@brief implements the VMC algorithm.
 */
class VMC: public QMCDriver
{
public:
  /// Constructor.
  VMC(MCWalkerConfiguration& w,
      TrialWaveFunction& psi,
      QMCHamiltonian& h);

  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMC(const VMC& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMC& operator=(const VMC&)
  {
    return *this;
  }

  void advanceWalkerByWalker();

};
}

#endif
