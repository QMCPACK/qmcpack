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
#ifndef QMCPLUSPLUS_LOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_LOCALENERGYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus {

  /** Class to accumulate the local energy and components
   *
   * Use Walker::Properties to accumulate Hamiltonian-related quantities.
   */
  class LocalEnergyEstimator: public ScalarEstimatorBase 
  {

    enum {ENERGY_INDEX, ENERGY2_INDEX, POTENTIAL_INDEX, LE_MAX};

    int FirstHamiltonian;
    int SizeOfHamiltonians;
    const QMCHamiltonian& refH;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    LocalEnergyEstimator(QMCHamiltonian& h);

    /** accumulation per walker
     * @param awalker current walker
     * @param wgt weight
     * 
     * Weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
     */
    inline void accumulate(const Walker_t& awalker, RealType wgt) 
    {
      const RealType* restrict ePtr = awalker.getPropertyBase();
      RealType wwght= wgt* awalker.Weight;
      scalars[0](ePtr[LOCALENERGY],wwght);
      scalars[1](ePtr[LOCALENERGY]*ePtr[LOCALENERGY],wwght);
      scalars[2](ePtr[LOCALPOTENTIAL],wwght);
      for(int target=3, source=FirstHamiltonian; target<scalars.size(); 
          ++target, ++source)
        scalars[target](ePtr[source],wwght);
    }

    /*@{*/
    inline void accumulate(const MCWalkerConfiguration& W
        , WalkerIterator first, WalkerIterator last, RealType wgt) 
    {
      for(; first != last; ++first) accumulate(**first,wgt);
    }
    void add2Record(RecordListType& record);
    void registerObservables(vector<observable_helper*>& h5dec, hid_t gid) {}
    ScalarEstimatorBase* clone();
    /*@}*/
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
