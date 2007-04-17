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
#ifndef QMCPLUSPLUS_LOCALENERGY_ONLY_ESTIMATOR_H
#define QMCPLUSPLUS_LOCALENERGY_ONLY_ESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
namespace qmcplusplus {

  class LocalEnergyOnlyEstimator: public ScalarEstimatorBase {

    enum {ENERGY_INDEX, ENERGY_SQ_INDEX, LE_MAX};
    ///index of local energy
    int LocalEnergyIndex;
    ///sum of local energy 
    RealType eSum;
    ///sum of local energy squared
    RealType e2Sum;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    inline LocalEnergyOnlyEstimator(): eSum(0),e2Sum(0) {}

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     * @param record storage of scalar records (name,value)
     * @param msg buffer for message passing
     */
    inline void add2Record(RecordListType& record, BufferType& msg) {
      LocalEnergyIndex = record.add("LocalEnergy");
      int VarianceIndex = record.add("Variance");
      msg.add(eSum);
      msg.add(e2Sum);
    }

    inline void accumulate(const Walker_t& awalker, RealType wgt) {
      RealType e(awalker.Properties(LOCALENERGY));
      eSum += wgt*e;
      e2Sum += wgt*e*e;
    }

    inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) 
    {
      RealType e(awalker.Properties(LOCALENERGY));
      RealType wgt(awalker.Weight);
      eSum += wgt*e;
      e2Sum += wgt*e*e;
    }

    inline RealType accumulate(WalkerIterator first, WalkerIterator last) {
      //RealType deltaE = Href.getEnsembleAverage();
      int wsum=0;
      while(first != last) {
        RealType e((*first)->Properties(LOCALENERGY));
        eSum += e;
        e2Sum += e*e;
        ++first; wsum++;
      }
      return wsum;
      //elocal[ENERGY_INDEX] += static_cast<RealType>(wsum)*deltaE;
    }

    ///reset all the cumulative sums to zero
    inline void reset() { 
      eSum=0.0;
      e2Sum=0.0;
    }

    ///copy the value to a message buffer
    inline void copy2Buffer(BufferType& msg) {
      msg.put(eSum);
      msg.put(e2Sum);
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    inline void report(RecordListType& record, RealType wgtinv) 
    {
      record[LocalEnergyIndex] =d_average = eSum*wgtinv;
      record[LocalEnergyIndex+1] =d_variance = e2Sum*wgtinv-d_average*d_average;
      d_sum=eSum; eSum=0.0;
      d_sumsq=e2Sum; e2Sum=0.0;
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    inline void report(RecordListType& record, RealType wgtinv, BufferType& msg) {
      msg.put(eSum);
      msg.put(e2Sum);
      report(record,wgtinv);
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
