//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_REPMULTIPLE_H
#define QMCPLUSPLUS_REPMULTIPLE_H

#include "QMCDrivers/QMCDriver.h"
#include "OhmmsPETE/OhmmsVector.h"
//#include "Estimators/RQMCMultipleEstimator.h"

namespace qmcplusplus
{

class Bead;
class MultiChain;
class PolymerEstimator;

/** @ingroup QMCDrivers MultiplePsi
 * @brief Implements the RMC algorithm for energy differences
 */
class RQMCMultiple: public QMCDriver
{

public:

  enum { MinusDirection=0, PlusDirection=1, Directionless=2};

  /// Constructor.
  RQMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool);

  /// Destructor
  ~RQMCMultiple();

  bool run();
  bool put(xmlNodePtr q);

protected:

  /** boolean for initialization
   *
   *\if true,
   *use clones for a chain.
   *\else
   *use drift-diffusion to form a chain
   *\endif
   */
  ///The length of polymers
  int ReptileLength;

  PolymerEstimator *multiEstimator;

  ///
  int forward,backward,itail,inext;

  ///the number of turns per block
  int NumTurns;

  int nptcl,nObs;

  ///the number of H/Psi pairs
  int nPsi;
  ///use scaled drift "true" = yes (default) else no
  std::string scaleBeadDrift;
  ///The Reptile: a chain of beads
  MultiChain* Reptile;

  ///The new bead
  Bead *NewBead;

  ///move polymers
  void moveReptile();

  ///initialize polymers
  void initReptile();

  ///for the first run starting with a point, set reference properties (sign)
  void setReptileProperties();

  void checkReptileProperties();

  ///Working arrays
  Vector<RealType> NewGlobalAction,DeltaG;
  Vector<int> NewGlobalSignWgt,WeightSign;

  void resizeArrays(int n);

  ///overwrite recordBlock
  void recordBlock(int block);

private:

  /// Copy Constructor (disabled)
  RQMCMultiple(const RQMCMultiple& a): QMCDriver(a) { }

  /// Copy operator (disabled).
  RQMCMultiple& operator=(const RQMCMultiple&)
  {
    return *this;
  }

  ParticleSet::ParticlePos_t gRand;
  RealType MSS,Tauoverm,sqrtTauoverm;
  int maxTouch,EvalInterval;
};
}
#endif
