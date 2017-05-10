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
    
    
#ifndef QMCPLUSPLUS_VMCMoveAll_H
#define QMCPLUSPLUS_VMCMoveAll_H
#include "QMC/QMCDriver.h"
namespace qmcplusplus
{

/** Implements the VMCMoveAll algorithm.
 *
 * This class is just for a record. Not being used by qmcPlusPlus applications.
 * Possible that we can use this class for Vector machines!
 */
class VMCMoveAll: public QMCDriver
{
public:
  /// Constructor.
  VMCMoveAll(MCWalkerConfiguration& w,
             TrialWaveFunction& psi,
             QMCHamiltonian& h,
             xmlNodePtr q);

  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMCMoveAll(const VMCMoveAll& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCMoveAll& operator=(const VMCMoveAll&)
  {
    return *this;
  }

  ///temporary storage for drift
  ParticleSet::ParticlePos_t drift;

  ///temporary storage for random displacement
  ParticleSet::ParticlePos_t deltaR;

  void advanceAllWalkers();
};
}

#endif
