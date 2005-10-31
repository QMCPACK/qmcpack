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
#ifndef QMCPLUSPLUS_WAVEFUNCTIONTEST_H
#define QMCPLUSPLUS_WAVEFUNCTIONTEST_H

#include "QMC/QMCDriver.h" 
namespace qmcplusplus {

  /** Test the correctness of TrialWaveFunction for the values,
      gradients and laplacians
  */
  class WaveFunctionTester: public QMCDriver 
  {
  public:
    /// Constructor.
    WaveFunctionTester(MCWalkerConfiguration& w, 
		       TrialWaveFunction& psi, 
		       QMCHamiltonian& h, 
		       xmlNodePtr q);

    bool run();
    bool put(xmlNodePtr q);
  private:
    /// Copy Constructor (disabled)
    WaveFunctionTester(const WaveFunctionTester& a): QMCDriver(a) { }
    /// Copy Operator (disabled)
    WaveFunctionTester& operator=(const WaveFunctionTester&) { return *this;}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
