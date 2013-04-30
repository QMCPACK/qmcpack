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
  string UseDrift;
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
/***************************************************************************
 * $RCSfile: VMCSingle.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCSingle.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
