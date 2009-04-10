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
#ifndef QMCPLUSPLUS_FORWARDWALKING_ESTIMATOR_H
#define QMCPLUSPLUS_FORWARDWALKING_ESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus {

  /** Class to accumulate the local energy and components
   *
   * Use Walker::Properties to accumulate Hamiltonian-related quantities.
   */
  class ForwardWalkingEstimator: public ScalarEstimatorBase 
  {

    enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

    int FirstHamiltonian;
    int SizeOfHamiltonians;
    const QMCHamiltonian& refH;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    ForwardWalkingEstimator(QMCHamiltonian& h);

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
      RealType wwght2= wwght*wwght;
      scalars[0](ePtr[LOCALENERGY],wwght);
      scalars[1](ePtr[LOCALENERGY]*ePtr[LOCALENERGY],wwght2);
      scalars[2](ePtr[LOCALPOTENTIAL],wwght);
      scalars[3](ePtr[LOCALPOTENTIAL]*ePtr[LOCALPOTENTIAL],wwght2);
      for(int target=4, source=FirstHamiltonian; target<scalars.size(); target+=2, ++source)
      {
        scalars[target](ePtr[source],wwght);
        scalars[target+1](ePtr[source]*ePtr[source],wwght2);
      }
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
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3503 $   $Date: 2009-02-02 11:24:37 -0600 (Mon, 02 Feb 2009) $
 * $Id: ForwardWalkingEstimator.h 3503 2009-02-02 17:24:37Z jnkim $ 
 ***************************************************************************/
