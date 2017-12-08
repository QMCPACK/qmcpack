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





/** @file accumulators.h
 * @brief Define and declare accumulator_set
 *
 * A temporary implementation to handle scalar samples and will be replaced by
 * boost.Accumulator
 */
#ifndef QMCPLUSPLUS_ACCUMULATORS_H
#define QMCPLUSPLUS_ACCUMULATORS_H

#include <config/stdlib/limits.h>
#include <iostream>
#include <type_traits>

/** generic accumulator of a scalar type
 *
 * To simplify i/o, the values are storged in contens
 */
template<typename T,
  typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
struct accumulator_set
{
  typedef T value_type;
  typedef T return_type;

  enum {VALUE=0,VALUESQ=1,WEIGHT=2,CAPACITY=4};
  T properties[CAPACITY];

  inline accumulator_set()
  {
    for(int i=0; i<CAPACITY; ++i)
      properties[i]=value_type();
  }

  /** add a sample */
  inline void operator()(value_type x)
  {
    properties[VALUE]  +=x;
    properties[VALUESQ]+=x*x;
    properties[WEIGHT] +=1.0;
  }

  /** add a sample and  */
  inline void operator()(value_type x, value_type w)
  {
    properties[VALUE]  +=w*x;
    properties[VALUESQ]+=w*x*x;
    properties[WEIGHT] +=w;
  }

  /** reset properties
   * @param v cummulative value
   * @param vv cummulative valuesq
   * @param w cummulative weight
   */
  inline void reset(value_type v, value_type vv, value_type w)
  {
    properties[VALUE]  =v;
    properties[VALUESQ]=vv;
    properties[WEIGHT] =w;
  }

  /** reset properties
   * @param v cummulative value
   * @param w cummulative weight
   */
  inline void reset(value_type v, value_type w)
  {
    properties[VALUE]  =v;
    properties[VALUESQ]=v*v;
    properties[WEIGHT] =w;
  }

  /** add a value but set the weight 1
   *
   * @todo Jeremy provides the reasonin of having this function. Suggest rename it to make the meaning clear.
   */
  inline void add(value_type x)
  {
    properties[VALUE]+=x;
    properties[VALUESQ] += x*x;
    properties[WEIGHT]=1;
  }

  /** return true if Weight!= 0 */
  inline bool good() const
  {
    return properties[WEIGHT]>0;
  }
  /** return true if Weight== 0 */
  inline bool bad() const
  {
    return iszero(properties[WEIGHT]);
  }

  /** return the sum */
  inline return_type result() const
  {
    return properties[VALUE];
  }

  /** return the sum of value squared */
  inline return_type result2() const
  {
    return properties[VALUESQ];
  }
  /** return the count
   *
   * Will return the sum of weights of each sample
   */
  inline return_type count() const
  {
    return properties[WEIGHT];
  }

  inline std::pair<return_type,return_type> mean_and_variance() const
  {
    value_type norm=1.0/properties[WEIGHT];
    value_type avg=properties[VALUE]*norm;
    return std::pair<return_type,return_type>(avg,norm*properties[VALUESQ]-avg*avg);
  }

  ///return the mean
  inline return_type mean() const
  {
    return good()?properties[VALUE]/properties[WEIGHT]:0.0;
  }

  ///return the mean of squared values
  inline return_type mean2() const
  {
    return good()?properties[VALUESQ]/properties[WEIGHT]:0.0;
  }

  inline return_type variance() const
  {
    if(iszero(properties[WEIGHT]))
      return std::numeric_limits<T>::max();
    value_type norm=1.0/properties[WEIGHT];
    return norm*(properties[VALUESQ]-properties[VALUE]*properties[VALUE]*norm);
  }

  inline void clear()
  {
    for(int i=0; i<CAPACITY; ++i)
      properties[i]=value_type();
  }
};

template<typename ACC>
inline typename ACC::value_type mean(const ACC& ac)
{
  return ac.mean();
}

template<typename T>
std::ostream& operator<<(std::ostream& os, accumulator_set<T>& rhs)
{
  os << "accumulator_set: "
     << " value = " << rhs.properties[0]
     << " value_sq = " << rhs.properties[1]
     << " weight = " << rhs.properties[2];
  return os;
}


#endif
