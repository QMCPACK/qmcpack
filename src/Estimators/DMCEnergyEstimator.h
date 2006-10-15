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
#ifndef QMCPLUSPLUS_DMC_LOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_DMC_LOCALENERGYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "Message/CommOperators.h"

namespace qmcplusplus {

  template<class T>
  struct DMCEnergyEstimator: public ScalarEstimatorBase<T> {

    enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX,LE_MAX};

    typedef typename ScalarEstimatorBase<T>::Walker_t Walker_t;
    typedef typename ScalarEstimatorBase<T>::WalkerIterator WalkerIterator;
    typedef typename ScalarEstimatorBase<T>::BufferType BufferType;

    using ScalarEstimatorBase<T>::CollectSum;
    using ScalarEstimatorBase<T>::b_average;
    using ScalarEstimatorBase<T>::b_variance;

    int averageIndex;
    int varianceIndex;
    TinyVector<std::string,LE_MAX> obName;
    TinyVector<T,LE_MAX> observables;

    /** constructor
     */
    DMCEnergyEstimator() { 
      obName[ENERGY_INDEX]="LocalEnergy";
      obName[ENERGY_SQ_INDEX]="Variance";
    }

    inline T getWeight() const { return observables[WEIGHT_INDEX];}

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     * @param record storage of scalar records (name,value)
     * @param msg buffer for message passing
     */
    void add2Record(RecordNamedProperty<T>& record, BufferType& msg) {
      averageIndex = record.add(obName[ENERGY_INDEX].c_str());
      varianceIndex= record.add(obName[ENERGY_SQ_INDEX].c_str());
      //add elocal to the message buffer
      msg.add(observables.begin(),observables.begin()+1);
    }

    inline void accumulate(const Walker_t& awalker, T wgt) {
      T e(awalker.Properties(LOCALENERGY));
      observables[ENERGY_INDEX] += wgt*e;
      observables[ENERGY_SQ_INDEX] += wgt*e*e;
      observables[WEIGHT_INDEX] += wgt;
    }

    inline void collect(WalkerIterator first, WalkerIterator last) {
      TinyVector<T,LE_MAX> temp;
      while(first != last) {
        T wgt((*first)->Weight);
        T e((*first)->Properties(LOCALENERGY));
        temp[ENERGY_INDEX] += wgt*e;
        temp[ENERGY_SQ_INDEX] += wgt*e*e;
        temp[WEIGHT_INDEX] += wgt;
        ++first;
      }
      gsum(temp,0);
      T wgtInv=1.0/temp[WEIGHT_INDEX];
      observables[ENERGY_INDEX] += temp[ENERGY_INDEX]*wgtInv;
      observables[ENERGY_SQ_INDEX] += temp[ENERGY_SQ_INDEX]*wgtInv;
      observables[WEIGHT_INDEX] += temp[WEIGHT_INDEX];
    }

    inline void accumulate(WalkerIterator first, WalkerIterator last, T wgtnorm) {
      while(first != last) {
        T wgt((*first)->Weight*wgtnorm);
        T e((*first)->Properties(LOCALENERGY));
        observables[ENERGY_INDEX] += wgt*e;
        observables[ENERGY_SQ_INDEX] += wgt*e*e;
        observables[WEIGHT_INDEX] += wgt;
        ++first;
      }
    }

    ///reset all the cumulative sums to zero
    inline void reset() { 
      observables=0.0;
    }

    ///copy the value to a message buffer
    inline void copy2Buffer(BufferType& msg) {
      msg.put(observables.begin(),observables.begin()+1);
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    inline void report(RecordNamedProperty<T>& record, T wgtinv, BufferType& msg) {
      b_average =  observables[ENERGY_INDEX]*wgtinv;
      b_variance = observables[ENERGY_SQ_INDEX]*wgtinv-b_average*b_average;
      record[averageIndex]=b_average;
      record[varianceIndex]=b_variance;
      observables=0.0;
    }
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
