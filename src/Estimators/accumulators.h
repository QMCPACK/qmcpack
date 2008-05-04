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
/** @file accumulators.h
 * @brief Define and declare accumulator_set
 *
 * A temporary implementation to handle scalar samples and will be replaced by
 * boost.Accumulator
 */
#ifndef QMCPLUSPLUS_ACCUMULATORS_H
#define QMCPLUSPLUS_ACCUMULATORS_H

/** generic accumulator of a scalar type **/
template<typename T>
class accumulator_set
{
  public:
    typedef T value_type;
    typedef T return_type;

    inline accumulator_set(): val(value_type()), valsq(value_type()),weight(value_type()){}

    /** add a sample */
    inline void operator()(value_type x) { 
      val+=x; valsq+=x*x;weight+=1.0; 
    }

    /** add a sample and  */
    inline void operator()(value_type x, value_type w) { 
      val+=x*w; valsq+=w*x*x;weight+=w; 
    }

    /** return the sum */
    inline return_type result() const { return val;}
    /** return the count
     *
     * Will return the sum of weights of each sample
     */ 
    inline return_type count() const { return weight;}

    inline pair<return_type,return_type> mean_and_variance() const { 
      value_type norm=1.0/weight;
      value_type avg=val*norm;
      return pair<return_type,return_type>(avg,norm*valsq-avg*avg);
    }

    inline return_type mean() const { return val/weight; }

    inline return_type variance() const { 
      value_type norm=1.0/weight;
      return norm*(valsq-val*val*norm);
    }

    inline void clear()
    {
      val=value_type();valsq=value_type();weight=value_type();
    }

  private:
    value_type weight;
    value_type val;
    value_type valsq;
    value_type dummy;//just for alignment
};

  template<typename ACC>
inline typename ACC::value_type mean(const ACC& ac)
{
  return ac.mean();
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1931 $   $Date: 2007-04-21 11:56:27 -0500 (Sat, 21 Apr 2007) $
 * $Id: ScalarEstimatorBase.h 1931 2007-04-21 16:56:27Z jnkim $ 
 ***************************************************************************/
