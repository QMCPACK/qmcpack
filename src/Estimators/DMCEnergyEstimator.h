//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_DMC_LOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_DMC_LOCALENERGYESTIMATOR_H
#include "Message/CommOperators.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** EnergyEstimator specialized for DMC
 *
 * This is a fake ScalarEstimatorBase<T> in that nothing is accumulated by
 * this object. However, it uses WalkerControlBase* which performs
 * the functions of an estimator and more.
 */
struct DMCEnergyEstimator: public ScalarEstimatorBase
{

  enum {ENERGY_INDEX=0, ENERGY_SQ_INDEX, WALKERSIZE_INDEX, WEIGHT_INDEX, EREF_INDEX, LE_MAX};

  ///WalkerControlBase* which actually performs accumulation
  WalkerControlBase* walkerControl;

  /** constructor
   */
  DMCEnergyEstimator():walkerControl(0)
  {
    //this handles three scalars: LocalEnergy, LocalEnergy2, Weight, Weight2
    d_data.resize(4);
  }

  ScalarEstimatorBase* clone()
  {
    DMCEnergyEstimator *d = new DMCEnergyEstimator();
    d->walkerControl=walkerControl;
    return d;
  }

  void setWalkerControl(WalkerControlBase* wc)
  {
    walkerControl=wc;
  }

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   * @param record storage of scalar records (name,value)
   * @param msg buffer for message passing
   *
   * Nothing to be added to msg.
   */
  void add2Record(RecordListType& record, BufferType& msg)
  {
    FirstIndex = record.add("LocalEnergy");
    int dummy=record.add("LocalEnergy2");
    dummy=record.add("Weight");
    dummy=record.add("Weight2");
    LastIndex = FirstIndex+4;
    //popIndex     = record.add("Population");
  }

  inline void accumulate(const Walker_t& awalker, RealType wgt)
  {
    app_error() << "Disabled DMCEnergyEstimator::accumulate(const Walker_t& awalker, T wgt) " << std::endl;
  }

  inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker)
  {
  }

  /** accumulate a Current counter for normalizations.
   *
   * walkerControl accumulate the local energy and weights
   */
  inline void accumulate(WalkerIterator first, WalkerIterator last)
  {
    RealType wgt(walkerControl->getCurrentValue(WEIGHT_INDEX));
    RealType norm=1.0/wgt;
    d_data[0] += norm*walkerControl->getCurrentValue(ENERGY_INDEX);
    d_data[1] += norm*walkerControl->getCurrentValue(ENERGY_SQ_INDEX);
    d_data[2] += wgt;
    d_data[3] += wgt*wgt;
    d_wgt += 1.0;
    //d_data[0]+= walkerControl->getCurrentValue(ENERGY_INDEX);
    //d_data[1]+= walkerControl->getCurrentValue(ENERGY_SQ_INDEX);
    ////d_data[2]+= walkerControl->getCurrentValue(WALKERSIZE_INDEX);
    //d_wgt += walkerControl->getCurrentValue(WEIGHT_INDEX);
  }

  ///reset all the cumulative sums to zero
  inline void reset()
  {
    d_wgt=0.0;
    std::fill(d_data.begin(), d_data.end(),0.0);
    walkerControl->reset();
  }

  /** calculate the averages and reset to zero
   * @param record a container class for storing scalar records (name,value)
   * @param wgtinv the inverse weight
   *
   * The weight is overwritten by 1.0/Current, since the normalized weight over
   * an ensemble, i.e., over the walkers, has been already
   * included by WalkerControlBase::branch operation.
   */
  inline void report(RecordListType& record, RealType wgtinv)
  {
    wgtinv=1.0/d_wgt;
    d_sum=walkerControl->getValue(ENERGY_INDEX);
    d_sumsq=walkerControl->getValue(ENERGY_SQ_INDEX);
    record[FirstIndex]= d_average =  wgtinv*d_sum;
    record[FirstIndex+1]= d_variance = wgtinv*d_sumsq-d_average*d_average;
    record[FirstIndex+2]= wgtinv*walkerControl->getValue(WALKERSIZE_INDEX);
    walkerControl->reset();
    walkerControl->setEnergyAndVariance(d_average,d_variance);
  }

  inline void report(RecordListType& record, RealType wgtinv, BufferType& msg)
  {
    report(record,wgtinv);
  }
};

}
#endif
