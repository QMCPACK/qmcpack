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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_CORRELATED_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_CORRELATED_POLYMERESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"

namespace qmcplusplus {

  class QMCHamiltonian;
  class TrialWaveFunction;
  class MultiChain;

  struct CSPolymerEstimator: public ScalarEstimatorBase {

    //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, PE_INDEX, KE_INDEX, LE_INDEX};
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

    /** constructor
     * @param h QMCHamiltonian to define the components
     * @param hcopy number of copies of QMCHamiltonians
     */
    CSPolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0); 
    CSPolymerEstimator(const CSPolymerEstimator& mest);
    ScalarEstimatorBase* clone();

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

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     *@param record storage of scalar records (name,value)
     */
    void add2Record(RecordNamedProperty<RealType>& record);

    inline  void accumulate(const Walker_t& awalker, RealType wgt) {}

    /** @warning Incomplete. Only to avoid compiler problems
     */
    //inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) { }

    void accumulate(WalkerIterator first, WalkerIterator last, RealType wgt);

    void evaluateDiff();
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: CSPolymerEstimator.h 1926 2007-04-20 17:30:26Z jnkim $ 
 ***************************************************************************/
