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
    
    



#ifndef QMCPLUSPLUS_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_POLYMERESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"

namespace qmcplusplus
{

class QMCHamiltonian;
class TrialWaveFunction;
class MultiChain;

struct PolymerEstimator: public ScalarEstimatorBase
{
  enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, LE_INDEX};

  ///The Reptile: a chain of beads
  MultiChain* Reptile;
  ///number of correlated systems
  int NumCopies;
  ///number of observables
  int NumObservables;
  ///Time step
  RealType Tau;
  ///1/Tau
  RealType OneOverTau;

  inline PolymerEstimator(): NumCopies(0), NumObservables(0), Reptile(0), Tau(0), OneOverTau(0) {}
//

  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  PolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0) {};

  PolymerEstimator(const PolymerEstimator& mest):
//      FirstIndex(mest.FirstIndex),LastIndex(mest.LastIndex),scalars(mest.scalars), scalars_saved(mest.scalars_saved),
    Reptile(mest.Reptile),NumCopies(mest.NumCopies),NumObservables(mest.NumObservables),Tau(mest.Tau),
    OneOverTau(mest.OneOverTau)
  {}

//     ScalarEstimatorBase* clone();

  inline RealType getUmbrellaWeight(int ipsi)
  {
    return scalars_saved[ipsi*LE_INDEX+WEIGHT_INDEX].result();
    //return d_data[ipsi*LE_INDEX+WEIGHT_INDEX];
  }

  inline void setTau(RealType dt)
  {
    Tau=dt;
    OneOverTau=1.0/dt;
  }

  inline void setPolymer(MultiChain* polymer)
  {
    Reptile=polymer;
  }

  void evaluateDiff() {};
  void setConstants(double x) {};
  //virtual ScalarEstimatorBase* clone()=0;
  //virtual void add2Record(RecordNamedProperty<RealType>& record) = 0;
  //virtual void accumulate(const Walker_t& awalker, RealType wgt) = 0;
  //virtual void accumulate(const MCWalkerConfiguration& W, WalkerIterator first, WalkerIterator last, RealType wgt) = 0;
};

}
#endif
