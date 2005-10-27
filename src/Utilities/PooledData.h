//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
    this->clear();
    if(a.size()) {
      this->insert(this->begin(), a.begin(), a.end());
    }
    return *this;
  }

  inline int current() const { return Current;}

  ///set the Current to zero
  inline void rewind() { Current = 0;}

  template<class T1>
  inline void add(T1 x) { Current++; this->push_back(static_cast<T>(x));}

  inline void add(T x) { Current++; this->push_back(x);}

  template<class _InputIterator>
  inline void add(_InputIterator first, _InputIterator last) {
    while(first != last) {
      Current++; this->push_back(*first++);
    }
  }

  template<class T1>
  inline void get(T1& x) {
    x = static_cast<T1>((*this)[Current++]);
  }

  inline void get(T& x) { x = (*this)[Current++];}

  template<class _OutputIterator>
  inline void get(_OutputIterator first, _OutputIterator last) {
    while(first != last) {
      *first++ = (*this)[Current++];
    }
  }
  
  inline void put(T x) { (*this)[Current++] = x;}

  template<class _InputIterator>
  inline void put(_InputIterator first, _InputIterator last){
    while(first != last) {
      (*this)[Current++] = *first++;
    }
  }


  /** return the address of the first element **/
  inline T* data() { return &((*this)[0]);}

  inline void print(std::ostream& os) {
    std::copy(this->begin(), this->end(), ostream_iterator<T>(os," "));
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m) {
    m.Pack(&((*this)[0]),this->size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m) {
    m.Unpack(&((*this)[0]),this->size());
    return m;
  }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
