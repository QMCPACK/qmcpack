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
#ifndef QMCPLUSPLUS_LOCALENERGY_ONLY_ESTIMATOR_H
#define QMCPLUSPLUS_LOCALENERGY_ONLY_ESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
namespace qmcplusplus {

  /** Estimator for local energy only
   */
  struct LocalEnergyOnlyEstimator: public ScalarEstimatorBase {

    inline LocalEnergyOnlyEstimator() 
    {
      scalars.resize(2);
      scalars_saved.resize(2);
    }

    ScalarEstimatorBase* clone()
    {
      return new LocalEnergyOnlyEstimator();
    }

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     * @param record storage of scalar records (name,value)
     */
    inline void add2Record(RecordListType& record) {
      FirstIndex = record.add("LocalEnergy");
      int s1=record.add("LocalPotential");
      LastIndex = FirstIndex+2;
      // int s2=record.add("KineticEnergy");
      //LastIndex = FirstIndex+3;
      clear();
    }

    inline void accumulate(const Walker_t& awalker, RealType wgt) 
    {
      scalars[0](awalker.Properties(LOCALENERGY),wgt);
      scalars[1](awalker.Properties(LOCALPOTENTIAL),wgt);
      //scalars[2](awalker.Properties(NUMPROPERTIES),wgt);
    }

    inline void accumulate(WalkerIterator first, WalkerIterator last, RealType wgt) 
    {
      while(first != last) 
      {
        accumulate(**first,wgt);++first;
      }
    }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
