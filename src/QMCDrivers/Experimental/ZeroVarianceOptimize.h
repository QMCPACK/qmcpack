//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ZeroVarianceOptimize.h
 * @brief Definition of QMCDriver which performs zero-variance optimization.
 */
#ifndef QMCPLUSPLUS_ZEROVARIANCEOPTIMIZATION_H
#define QMCPLUSPLUS_ZEROVARIANCEOPTIMIZATION_H

#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"

namespace qmcplusplus
{

class DiffOrbitalBase;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization based on the zero-variance method
 *
 * To be derived from a clone manager
 */
class ZeroVarianceOptimize: public QMCDriver
{
public:

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::iterator WalkerIter_t;
  ///Constructor.
  ZeroVarianceOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                       QMCHamiltonian& h );//, HamiltonianPool& hpool);

  ///Destructor
  ~ZeroVarianceOptimize();
  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  ///set the wavefunction node
  void setWaveFunctionNode(xmlNodePtr cur)
  {
    wfNode=cur;
  }

private:
  ///Size of Hessian and Overalp Nopt+1
  int NoptPlusOne;
  ///total number of Warmup Blocks
  int WarmupBlocks;
  ///update engine
  QMCUpdateBase* Mover;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;
  ///option to use drift
  string UseDrift;

  ///these going to be matrix
  Vector <Vector<RealType> > dLogPsi;
  Vector <Vector<RealType> > dHPsi;
  Vector <Matrix<RealType> > Hessian;
  Vector <Matrix<RealType> > Overlap;


  ///hide initialization from the main function
  void resetRun();
  ///accumulate matrix elements
  void accumulate(WalkerIter_t it, WalkerIter_t it_end);

  ///Copy Constructor (disabled).
  ZeroVarianceOptimize(const ZeroVarianceOptimize& a): QMCDriver(a) { }
  ///Copy operator (disabled).
  ZeroVarianceOptimize& operator=(const ZeroVarianceOptimize&)
  {
    return *this;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: ZeroVarianceOptimize.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
