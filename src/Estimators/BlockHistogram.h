//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file BlockHistogram.h
 * @brief Define a accumulator whose average is evaluated for a moving block of  a fixed steps
 *
 * Use deque for the buffer but will use boost::circular_buffer with the update
 */
#ifndef QMCPLUSPLUS_BLOCKHISTOGRAM_H
#define QMCPLUSPLUS_BLOCKHISTOGRAM_H
#include <deque>

template<class T>
struct BlockHistogram
{
  int           maxSize;
  T             mySum;
  T             myWeightInv;
  std::deque<T> myData;

  //typedef for iterator
  typedef T         value_type;
  typedef T&        reference;
  typedef typename  std::deque<T>::size_type size_type;
  typedef typename  std::deque<T>::iterator iterator;
  typedef typename  std::deque<T>::const_iterator const_iterator;

  /** default constructor
   * @param n size of the samples defult=100
   */
  explicit inline BlockHistogram(size_type n=100): maxSize(n), mySum(T()), myWeightInv(1)
  {
  }

  /** resize the histogram
   */
  inline void resize(size_type n, const T& x)
  {
    myData.resize(n,x);
    mySum = n*x;
    myWeightInv=1/static_cast<T>(n);
  }

  inline void reserve(size_type n)
  {
    maxSize=n;
    this->clear();
  }

  inline void clear()
  {
    myData.clear();
    mySum = T();
    myWeightInv=1;
  }

  inline value_type mean() const
  {
    return mySum*myWeightInv;
  }

  inline value_type result() const
  {
    return mySum;
  }

  inline size_type size() const
  {
    return myData.size();
  }

  /** add a value x
   */
  inline void operator()(const T& x)
  {
    //add to the end
    myData.push_back(x);
    if(myData.size()<maxSize)
      mySum +=x;
    else
    {
      mySum += x - myData.front();
      myData.pop_front();
    }
    myWeightInv=1/static_cast<T>(myData.size());
  }

};
#endif
