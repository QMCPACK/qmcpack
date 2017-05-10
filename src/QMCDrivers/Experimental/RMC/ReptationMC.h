//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
