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

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    inline LocalEnergyOnlyEstimator() 
    {
      d_data.resize(LE_MAX,0);
    }

    ScalarEstimatorBase* clone()
    {
      return new LocalEnergyOnlyEstimator();
    }

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     * @param record storage of scalar records (name,value)
     * @param msg buffer for message passing
     */
    inline void add2Record(RecordListType& record, BufferType& msg) {
      LocalEnergyIndex = record.add("LocalEnergy");
      int VarianceIndex = record.add("LocalEnergy2");
      msg.add(d_data[0]);
      msg.add(d_data[1]);
    }

    inline void accumulate(const Walker_t& awalker, RealType wgt) {
      RealType e(awalker.Properties(LOCALENERGY));
      d_data[0] += wgt*e;
      d_data[1] += wgt*e*e;
      d_wgt+=wgt;
    }

    inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) 
    {
      RealType e(awalker.Properties(LOCALENERGY));
      RealType wgt(awalker.Weight);
      d_data[0] += wgt*e;
      d_data[1] += wgt*e*e;
      d_wgt+=wgt;
    }

    inline void accumulate(WalkerIterator first, WalkerIterator last) {
      while(first != last) {
        RealType e((*first)->Properties(LOCALENERGY));
        d_data[0] += e;
        d_data[1] += e*e;
        ++first; 
        d_wgt+=1.0;
      }
    }

    void reset()
    {
      d_wgt=0.0;
      std::fill(d_data.begin(), d_data.end(),0.0);
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    inline void report(RecordListType& record, RealType wgtinv) 
    {
      record[LocalEnergyIndex] =d_average = d_data[0]*wgtinv;
      record[LocalEnergyIndex+1] =d_variance = d_data[1]*wgtinv-d_average*d_average;
      d_sum=d_data[0]; d_data[0]=0.0;
      d_sumsq=d_data[1]; d_data[1]=0.0;
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    inline void report(RecordListType& record, RealType wgtinv, BufferType& msg) {
      msg.put(d_data[0]);
      msg.put(d_data[1]);
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
