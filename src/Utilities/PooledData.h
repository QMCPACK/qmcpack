//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/** @file PooledData.h
 * @brief Define a serialized buffer to store anonymous data
 *
 * Two implementations differ in how the iterators are handled.
 * USE_POOLED_DATA_ITERATOR selects an implementation using std::vector<T>::iterator 
 * Otherwise, the original implementation using an index Current is used.
 */
#ifndef QMCPLUSPLUS_POOLEDDATA_H
#define QMCPLUSPLUS_POOLEDDATA_H
#include <vector>
#include <complex>

#if defined(USE_POOLED_DATA_ITERATOR)
template<class T>
class PooledData: public std::vector<T> {

  typename std::vector<T>::iterator Anchor;

public:
  
  //default constructor
  inline PooledData() { }

  //copy constructor
  inline PooledData(const PooledData& a): std::vector<T>(a) { 
    Anchor=this->begin();
  }

  //constructor with a size
  explicit inline PooledData(int n): std::vector<T>(n)
  {
    Anchor=this->begin();
  }

  //assignement operator
  PooledData<T>& operator=(const PooledData<T>& a) 
  {
    this->clear();
    if(a.size()) {
      this->insert(this->begin(), a.begin(), a.end());
    }
    return *this;
  }

  ///return the location of the Anchor
  inline int current() const { return Anchor-this->begin();}

  ///set the Anchor at the first iterator
  inline void rewind() { Anchor=this->begin();}

  ///@{Add data to the pool
  template<class T1>
  inline void add(T1 x) { this->push_back(static_cast<T>(x));}

  inline void add(T x) { this->push_back(x);}

  inline void add(complex<T>& x) { 
    this->push_back(x.real()); 
    this->push_back(x.imag());
  }

  template<class _InputIterator>
    inline void add(_InputIterator first, _InputIterator last) {
      this->insert(this->end(),first,last);
    }

  inline void add(std::complex<T>* first,std::complex<T>* last) {
    while(first != last) {
      this->push_back((*first).real()); 
      this->push_back((*first).imag());
      ++first;
    }
  }
  ///@}

  ///@{Assign value to the arguments from the pool and advance the Anchor
  template<class T1>
  inline void get(T1& x) {
    x = static_cast<T1>(*Anchor++);
  }

  inline void get(T& x) { x = *Anchor++;} 

  inline void get(std::complex<T>& x) 
  { 
    x.real()=*Anchor++;
    x.imag()=*Anchor++;
  }

  template<class _ForwardIterator>
    inline void get(_ForwardIterator first, _ForwardIterator last) {
      typename std::vector<T>::iterator here=Anchor;
      Anchor += last-first;
      std::copy(here,Anchor,first);
    }

  inline void get(std::complex<T>* first, std::complex<T>* last)
  {
    while(first != last) {
      (*first).real()=*Anchor++;
      (*first).imag()=*Anchor++;
      ++first;
    }
  }
  ///@}
  
  ///@{Assign value from the arguments to the pool and advnace the Anchor
  /** Add  a value
   * @param x value to add
   */
  inline void put(T x) { *Anchor++ = x;}
  /** Add a complex<T>
   * @param x a complex value to add
   */
  inline void put(std::complex<T>& x) { 
    *Anchor++ = x.real();  
    *Anchor++ = x.imag();
  }

  /** Add values from [first,last)
   * @param first starting iterator
   * @param last ending iterator
   */
  template<class _FowardIterator>
  inline void put(_FowardIterator first, _FowardIterator last)
  {
    std::copy(first,last,Anchor);
    Anchor += last-first;
  }

  /** Add complex values from [first,last)
   * @param first starting iterator
   * @param last ending iterator
   */
  inline void put(std::complex<T>* first, std::complex<T>* last) 
  {
    while(first != last) {
      *Anchor++ = (*first).real();  
      *Anchor++ = (*first).imag();
      ++first;
    }
  }
  ///@}


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
#else
template<class T>
class PooledData: public std::vector<T> {

  int Current;

public:
  
  //default constructor
  inline PooledData(): Current(0) { }

  //copy constructor
  inline PooledData(const PooledData& a): std::vector<T>(a),Current(0) { }

  //constructor with a size
  explicit inline PooledData(int n): std::vector<T>(n) , Current(0){}

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

  inline void add(complex<T>& x) { Current+=2; 
    this->push_back(x.real()); 
    this->push_back(x.imag());
  }

  template<class _InputIterator>
  inline void add(_InputIterator first, _InputIterator last) {
    while(first != last) {
      Current++; this->push_back(*first++);
    }
  }

  inline void add(std::complex<T>* first,std::complex<T>* last) {
    while(first != last) {
      this->push_back((*first).real()); 
      this->push_back((*first).imag());
      Current+=2; ++first;
    }
  }

  template<class T1>
  inline void get(T1& x) {
    x = static_cast<T1>((*this)[Current++]);
  }

  inline void get(T& x) { x = (*this)[Current++];}

  inline void get(std::complex<T>& x) { 
    x=std::complex<T>((*this)[Current],(*this)[Current+1]); Current+=2;
  }

  template<class _OutputIterator>
  inline void get(_OutputIterator first, _OutputIterator last) {
    while(first != last) {
      *first++ = (*this)[Current++];
    }
  }

  inline void get(std::complex<T>* first, std::complex<T>* last){
    while(first != last) {
      (*first)=std::complex<T>((*this)[Current],(*this)[Current+1]);
      ++first; Current+=2;
    }
  }
  
  inline void put(T x) { (*this)[Current++] = x;}
  inline void put(std::complex<T>& x) { 
    (*this)[Current++] = x.real();  
    (*this)[Current++] = x.imag();
  }

  template<class _InputIterator>
  inline void put(_InputIterator first, _InputIterator last){
    while(first != last) {
      (*this)[Current++] = *first++;
    }
  }

  inline void put(std::complex<T>* first, std::complex<T>* last) {
    while(first != last) {
      (*this)[Current++] = (*first).real();
      (*this)[Current++] = (*first).imag();
      ++first;
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
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
