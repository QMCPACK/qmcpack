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
#ifndef QMCPLUSPLUS_REPTATION_H
#define QMCPLUSPLUS_REPTATION_H

#include "QMCDrivers/QMCDriver.h"
#include <deque>
namespace qmcplusplus
{


class PolymerChain;

/** @ingroup QMCDrivers
 * @brief Implements the RMC algorithm
 */
class ReptationMC: public QMCDriver
{

public:

  /// Constructor.
  ReptationMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  /// Destructor
  ~ReptationMC();

  bool run();
  bool put(xmlNodePtr q);

protected:

  typedef MCWalkerConfiguration::Walker_t Walker_t;

  ///boolean for using bounce algorithm. true, if bounce algorithm of D. Ceperley
  bool UseBounce;

  /** boolean for initialization
   *
   *\if true,
   *use clones for a chain.
   *\else
   *use drift-diffusion to form a chain
   *\endif
   */
  bool ClonePolymer;

  ///The length of polymers
  int PolymerLength;

  ///the number of the beads that will be cut
  int  NumCuts;

  ///the number of turns per block
  int NumTurns;

  ///PolymerChain
  PolymerChain* Reptile;

  ///move polymers
  void moveReptile();

  ///initialize polymers
  void initReptile();
private:

  /// Copy Constructor (disabled)
  ReptationMC(const ReptationMC& a): QMCDriver(a) { }

  /// Copy operator (disabled).
  ReptationMC& operator=(const ReptationMC&)
  {
    return *this;
  }


};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 759 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: ReptationMC.h 759 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
