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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_VECTOR_ESTIMATO_IMPL_H
#define QMCPLUSPLUS_VECTOR_ESTIMATO_IMPL_H

namespace qmcplusplus {

  /** Class to manage observables of vector types
   *
   * Container and utility class for observables like gofr, sk.
   * Do not inherit from VectorEstimatorImpl but use it.
   */
  template<typename T>
  struct VectorEstimatorImpl {

    ///name of the observables
    std::string Name;
    /** current accumulative data 
     *
     * d_data.size() = 2 * size of vector elements
     * For the i-th element
     * - d_data[2*i] += \f$= \sum_{s} O_i^s\f$
     * - d_data[2*i+1] += \f$= \sum_{s} O_i^s O_i^s\f$
     */
    vector<T> d_data;
    ///d_sum[i] block sum
    vector<T> d_sum;
    ///d_sum2[i] block sum of squared
    vector<T> d_sum2;

    ///default constructor
    inline VectorEstimatorImpl(){}
    
    ///default constructor
    explicit inline VectorEstimatorImpl(int n)
    { resize(n);}

    ///copy constructor
    VectorEstimatorImpl(const VectorEstimatorImpl& est): d_data(est.d_data), 
    d_sum(est.d_sum), d_sum2(est.d_sum2)
    {}

    ///destructo
    ~VectorEstimatorImpl(){}

    /** resize the data containers
     * @param n number of elements of the vector observable
     */
    inline void resize(int n)
    {
      d_data.resize(2*n,T());
      d_sum.resize(n,T());
      d_sum2.resize(n,T());
    }

    inline void init()
    {
      std::fill(d_data.begin(),d_data.end(),T());
      std::fill(d_sum.begin(),d_sum.end(),T());
      std::fill(d_sum2.begin(),d_sum2.end(),T());
    }

    /// zero the active data
    inline void reset()
    {
      std::fill(d_data.begin(),d_data.end(),T());
    }

    /** accumulate expectation values
     * @param first1 vector data
     * @param first2 weight data
     */
    template<typename IT1, typename IT2>
    inline void accumulate(IT1 first1, IT2 first2)
    {
      typename vector<T>::iterator it(d_data.begin());
      typename vector<T>::iterator it_end(d_data.end());
      while(it != it_end)
      {
        T v=(*first1)*(*first2++);//w[i]*v[i]
        (*it++)+=v;
        (*it++)+=v*(*first1++);//w[i]*v[i]*v[i]
      }
    }

    /** accumulate expectation values
     *\param awalker a single walker
     *\param wgt the weight
     */
    template<typename IT>
    inline void accumulate(IT first, T wgt)
    {
      typename vector<T>::iterator it(d_data.begin());
      typename vector<T>::iterator it_end(d_data.end());
      while(it != it_end)
      {
        (*it++)+=wgt*(*first);
        (*it++)+=wgt*(*first)*(*first);
        ++first;
      }
    }

    /** accumulate expectation values
     *\param awalker a single walker
     *\param wgt the weight
     */
    template<typename IT>
    inline void accumulate(IT first)
    {
      typename vector<T>::iterator it(d_data.begin());
      typename vector<T>::iterator it_end(d_data.end());
      while(it != it_end)
      {
        (*it++)+=(*first);
        (*it++)+=(*first)*(*first);
        ++first;
      }
    }


    inline void takeBlockAverage(double wgtnorm)
    {
      typename vector<T>::const_iterator it(d_data.begin());
      typename vector<T>::const_iterator it_end(d_data.end());
      typename vector<T>::iterator sit(d_sum.begin());
      typename vector<T>::iterator s2it(d_sum2.begin());
      while(it != it_end)
      {
        (*sit++)=wgtnorm*(*it++);
        (*s2it++)=wgtnorm*(*it++);
      }
    }
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1931 $   $Date: 2007-04-21 11:56:27 -0500 (Sat, 21 Apr 2007) $
 * $Id: VectorEstimatorImpl.h 1931 2007-04-21 16:56:27Z jnkim $ 
 ***************************************************************************/
