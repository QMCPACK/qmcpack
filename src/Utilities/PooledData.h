//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_QMC_POOLEDDATA_H
#define OHMMS_QMC_POOLEDDATA_H
#include <vector>

/**class PooledData<T>
 */
template<class T>
class PooledData: public std::vector<T> {

  size_t Current;

public:
  
  //default constructor
  inline PooledData(): Current(0) { }

  //copy constructor
  inline PooledData(const PooledData& a): std::vector<T>(a) { }

  //assignement operator
  PooledData<T>& operator=(const PooledData<T>& a) {
    clear();
    if(a.size()) {
      insert(begin(), a.begin(), a.end());
    }
    return *this;
  }

  ///set the Current to zero
  inline void rewind() { Current = 0;}

  inline void add(T x) { push_back(x);}

  template<class _InputIterator>
  inline void add(_InputIterator first, _InputIterator last) {
    while(first != last) {
      push_back(*first++);
    }
  }

  inline void get(T& x) { x = operator[](Current++);}
  template<class _OutputIterator>
  inline void get(_OutputIterator first, _OutputIterator last) {
    while(first != last) {
      *first++ = operator[](Current++);
    }
  }
  
  inline void put(T x) { operator[](Current++) = x;}
  template<class _InputIterator>
  inline void put(_InputIterator first, _InputIterator last){
    while(first != last) {
      operator[](Current++) = *first++;
    }
  }


  inline void print(std::ostream& os) {
    std::copy(begin(), end(), ostream_iterator<T>(os," "));
  }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
