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
namespace qmcplusplus {

  /** Estimator for local energy only
   */
  struct ForwardWalkingEstimator: public ScalarEstimatorBase {

    inline ForwardWalkingEstimator() 
    {
      scalars.resize(2);
      scalars_saved.resize(2);
    }

    inline void accumulate(const MCWalkerConfiguration& W
        , WalkerIterator first, WalkerIterator last, RealType wgt) 
    {
      for(; first != last; ++first)
      {
        scalars[0]((*first)->Properties(LOCALENERGY),wgt);
        scalars[1]((*first)->Properties(LOCALPOTENTIAL),wgt);
      }
    }

    void registerObservables(vector<observable_helper*>& h5dec, hid_t gid)
    {}
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
    ScalarEstimatorBase* clone()
    {
      return new ForwardWalkingEstimator();
    }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3503 $   $Date: 2009-02-02 11:24:37 -0600 (Mon, 02 Feb 2009) $
 * $Id: LocalEnergyOnlyEstimator.h 3503 2009-02-02 17:24:37Z jnkim $ 
 ***************************************************************************/
