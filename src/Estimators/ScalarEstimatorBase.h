//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim 
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
#ifndef QMCPLUSPLUS_SCALAR_ESTIMATORBASE_H
#define QMCPLUSPLUS_SCALAR_ESTIMATORBASE_H
#include "Configuration.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {

  /** Abstract class for an estimator of a scalar operator.
   *
   * ScalarEstimators derived from ScalarEstimatorBase  implement three main functions
   * - reset : reset the internal values so that observables can be accumulated
   * - accumulate : measure and accumulate its value and the square of the value
   * - report : evaluate the block average and variance
   * ScalarEstimatorBase and its derived classes do not perform any I/O function.
   */
  struct ScalarEstimatorBase: public QMCTraits {

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef MCWalkerConfiguration::iterator WalkerIterator;
    typedef RecordNamedProperty<RealType>   RecordListType;
    typedef PooledData<RealType>            BufferType;

    ///sum of a scalar observable
    RealType d_sum;
    ///sum of a scalar observable squared
    RealType d_sumsq;
    ///average 
    RealType d_average;
    ///variance
    RealType d_variance;

    inline ScalarEstimatorBase(): 
      d_sum(RealType()), d_sumsq(RealType()), d_average(RealType()), d_variance(RealType()){}
    virtual ~ScalarEstimatorBase(){}

    ///retrun average
    inline RealType average() const { return d_average;}
    ///retrun variance
    inline RealType variance() const { return d_variance;}

    /** add the content of the scalar estimator to the record
     *\param record scalar data list 
     *
     *Each ScalarEstimatorBase object adds a number of scalar data to record.
     */
    virtual void add2Record(RecordNamedProperty<RealType>& record, BufferType& msg) = 0;

    /** a virtual function to accumulate expectation values
     *\param awalker a single walker
     *\param wgt the weight
     */
    virtual void accumulate(const Walker_t& awalker, RealType wgt) = 0;

    virtual RealType accumulate(WalkerIterator first, WalkerIterator last) = 0;
    virtual void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker)=0;

    /** a virtual function to report the scalar estimator
     *@param record 
     *@param wgtinv inverse of the weight
     *
     *Evalaute the block-average and flush the internal data for new averages
     */
    virtual void report(RecordListType&, RealType wgtinv) = 0;

    /** a virtual function to report the scalar estimator
     * @param record 
     * @param wgtinv inverse of the weight
     * @param msg temporary buffer for MPI
     *
     * Evalaute the block-average and flush the internal data for new averages
     */
    virtual void report(RecordListType&, RealType wgtinv, BufferType& msg) = 0;

    /** update the message buffer */
    virtual void copy2Buffer(BufferType& msg)=0;

    /// a virtual function to flush the internal data to start a new block average
    virtual void reset() = 0;
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
