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
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/RecordProperty.h"
#include "Estimators/accumulators.h"

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

    typedef accumulator_set<RealType> accumulator_type;
    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef MCWalkerConfiguration::iterator WalkerIterator;
    typedef RecordNamedProperty<RealType>   RecordListType;

    ///first index within an record of the first element handled by an object
    int FirstIndex;
    ///last index within an record of the first element handled by an object
    int LastIndex;
    ///scalars to be measured
    vector<accumulator_type> scalars; 
    ///scalars saved
    vector<accumulator_type> scalars_saved;

    inline ScalarEstimatorBase(): FirstIndex(0), LastIndex(0) {}

    /** copy constructor */
    ScalarEstimatorBase(const ScalarEstimatorBase& est):
      FirstIndex(est.FirstIndex),LastIndex(est.LastIndex),
    scalars(est.scalars), scalars_saved(est.scalars_saved)
    {}

    virtual ~ScalarEstimatorBase(){}

    ///retrun average of the
    inline RealType average(int i=0) const { return scalars_saved[i].mean();}
    ///return a variance
    inline RealType variance(int i=0) const { return scalars_saved[i].variance();}
    ///retrun mean and variance
    inline pair<RealType,RealType> operator[](int i) const
    { return scalars[i].mean_and_variance();}

    ///clear the scalars to collect
    inline void clear() 
    {
      for(int i=0; i<scalars.size(); i++) scalars[i].clear();
    }

    /** take block average and write to a common container */
    template<typename IT>
    inline void takeBlockAverage(IT first)
    {
      first += FirstIndex;
      for(int i=0; i<scalars.size(); i++)
      {  
        *first++ = scalars[i].mean(); 
        scalars_saved[i]=scalars[i]; //save current block
        scalars[i].clear();
      }
    }

    ///clone the object
    virtual ScalarEstimatorBase* clone()=0;

    /** add the content of the scalar estimator to the record
     * @param record scalar data list 
     *
     * Each ScalarEstimatorBase object adds 1 to many accumulator_type 
     */
    virtual void add2Record(RecordNamedProperty<RealType>& record) = 0;

    /** a virtual function to accumulate expectation values
     *\param awalker a single walker
     *\param wgt the weight
     */
    virtual void accumulate(const Walker_t& awalker, RealType wgt) = 0;
    /** a virtual function to accumulate expectatio values 
     * @param first iterator for the first walker
     * @param last iterafor for the last walker
     * @param wgt weight
     */
    virtual void accumulate(WalkerIterator first, WalkerIterator last, RealType wgt) = 0;

  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
