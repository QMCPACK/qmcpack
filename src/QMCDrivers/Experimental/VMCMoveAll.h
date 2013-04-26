//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMoveAll.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
