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
#include <limits>

//#define USE_POOLED_DATA_ITERATOR 1
//typedef double RealType;

#if defined(USE_POOLED_DATA_ITERATOR)

template<class T>
struct PooledData {
  
  //typedef for iterator
  typedef typename  std::vector<T>::size_type size_type;
  typedef typename  std::vector<T>::iterator iterator;
  typedef typename  std::vector<T>::const_iterator const_iterator;

  //default constructor
  inline PooledData() { Anchor=myData.begin();}

  ////copy constructor
  //inline PooledData(const PooledData& a): myData(a.myData) { 
  //  Anchor=myData.begin();
  //}

  //constructor with a size
  explicit inline PooledData(int n): myData(n)
  {
    Anchor=myData.begin();
  }

  ////assignement operator
  //PooledData<T>& operator=(const PooledData<T>& a) 
  //{
  //  myData.clear();
  //  if(a.size()) {
  //    myData.reserve(a.size());
  //    myData.insert(myData.begin(), a.begin(), a.end());
  //  }
  //  return *this;
  //}

  inline size_type size() const { return myData.size();}
  inline size_type capacity() const { return myData.capacity();}

  ///return the location of the Anchor
  inline int current() const { return Anchor-myData.begin();}

  ///set the Anchor at the first iterator
  inline void rewind() { Anchor=myData.begin();}

  inline iterator begin() { return myData.begin();}
  inline iterator end() { return myData.end();}
  inline const_iterator begin() const { return myData.begin();}
  inline const_iterator end() const { return myData.end();}

  inline void reserve(size_type n) { myData.reserve(n);}
  inline void clear() { myData.clear();Anchor=myData.begin();}

  ///@{Add data to the pool
  template<class T1>
  inline void add(T1 x) { myData.push_back(static_cast<T>(x));}

  inline void add(T x) { myData.push_back(x);}

  inline void add(complex<T>& x) { 
    myData.push_back(x.real()); 
    myData.push_back(x.imag());
  }

  template<class _InputIterator>
    inline void add(_InputIterator first, _InputIterator last) {
      myData.insert(myData.end(),first,last);
    }

  inline void add(std::complex<T>* first,std::complex<T>* last) {
    while(first != last) {
      myData.push_back((*first).real()); 
      myData.push_back((*first).imag());
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
      //typename std::vector<T>::iterator here=Anchor;
      size_type dn=last-first;
      std::copy(Anchor,Anchor+dn,first);
      Anchor += dn;
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
  inline T* data() { return &(myData[0]);}

  inline void print(std::ostream& os) {
    std::copy(myData.begin(), myData.end(), ostream_iterator<T>(os," "));
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m) {
    m.Pack(&(myData[0]),myData.size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m) {
    m.Unpack(&(myData[0]),myData.size());
    return m;
  }

  iterator Anchor;
  std::vector<T> myData;
};
#else
template<class T>
struct PooledData 
{
  typedef T value_type;
  typedef typename std::vector<T>::size_type size_type;

  size_type Current;
  std::vector<T> myData;
  
  ///default constructor
  inline PooledData(): Current(0) { }

  ///constructor with a size
  inline PooledData(size_type n):Current(0){ myData.resize(n,T());}

  ///return the size of the data
  inline size_type size() const { return myData.size();}

  //inline void clear() { std::vector<T>::clear(); Current=0;}
  inline size_type current() const { return Current;}

  /** set the Current to a cursor
   * @param cur locator to which Current is assigned
   */
  inline void rewind(size_type cur=0) { Current = cur;}

  ///return the starting iterator
  inline typename std::vector<T>::iterator begin() { return myData.begin();}
  ///return the ending iterator
  inline typename std::vector<T>::iterator end() { return myData.end();}

  /*@{ matching functions to std::vector functions */
  ///clear the data and set Current=0
  inline void clear() { myData.clear(); Current=0;}
  ///reserve the memory using vector<T>::reserve
  inline void reserve(size_type n) { myData.reserve(n); Current=0;}
  ///resize
  inline void resize(size_type n, T val=T()) { myData.resize(n,val); Current=0;}
  ///return i-th value
  inline T operator[](size_type i) const { return myData[i];}
  ///return i-th value to assign
  inline T& operator[](size_type i) { return myData[i];}
  /*@}*/

  template<class T1>
  inline void add(T1 x) { Current++; myData.push_back(static_cast<T>(x));}

  inline void add(T x) { Current++; myData.push_back(x);}

  inline void add(complex<T>& x) { 
    Current+=2; 
    myData.push_back(x.real()); 
    myData.push_back(x.imag());
  }

  template<class _InputIterator>
  inline void add(_InputIterator first, _InputIterator last) {
    Current += last-first;
    myData.insert(myData.end(),first,last);
    //while(first != last) {
    //  Current++; this->push_back(*first++);
    //}
  }

  inline void add(std::complex<T>* first,std::complex<T>* last) {
    size_type dn=2*(last-first);
    //TEMPORARY FIX FOR SGI-IA64 
#if defined(SGI_IA64)  || defined(PROFILING_ON)
    for(;first != last; ++first)
    {
      myData.push_back((*first).real()); 
      myData.push_back((*first).imag());
    }
#else
    myData.insert(myData.end(),&(first->real()),&(first->real())+dn);
#endif
    Current += dn;
  }

  template<class T1>
  inline void get(T1& x) {
    x = static_cast<T1>(myData[Current++]);
  }

  inline void get(T& x) { x = myData[Current++];}

  inline void get(std::complex<T>& x) 
  { 
    x=std::complex<T>(myData[Current],myData[Current+1]); Current+=2;
  }

  template<class _OutputIterator>
  inline void get(_OutputIterator first, _OutputIterator last) {
    size_type now=Current;
    Current+=last-first;
    std::copy(myData.begin()+now,myData.begin()+Current,first);
    //while(first != last) {
    //  *first++ = (*this)[Current++];
    //}
  }

  inline void get(std::complex<T>* first, std::complex<T>* last){
    while(first != last) {
      (*first)=std::complex<T>(myData[Current],myData[Current+1]);
      ++first; Current+=2;
    }
  }
  
  inline void put(T x) { myData[Current++] = x;}
  inline void put(std::complex<T>& x) { 
    myData[Current++] = x.real();  
    myData[Current++] = x.imag();
  }

  template<class _InputIterator>
  inline void put(_InputIterator first, _InputIterator last){
    std::copy(first,last,myData.begin()+Current);
    Current+=last-first;
    //while(first != last) {
    //  (*this)[Current++] = *first++;
    //}
  }

  inline void put(std::complex<T>* first, std::complex<T>* last) {
    while(first != last) {
      myData[Current++] = (*first).real();
      myData[Current++] = (*first).imag();
      ++first;
    }
  }


  /** return the address of the first element **/
  inline T* data() { return &(myData[0]);}

  inline void print(std::ostream& os) {
    std::copy(myData.begin(), myData.end(), ostream_iterator<T>(os," "));
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m) {
    m.Pack(&(myData[0]),myData.size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m) {
    m.Unpack(&(myData[0]),myData.size());
    return m;
  }
};
#endif

/** operator to check if two buffers are identical */
template<class T>
bool operator==(const PooledData<T>& a, const PooledData<T>& b)
{
  if(a.size() != b.size()) return false;
  //if(a.Current != b.Current) return false;
  for(typename PooledData<T>::size_type i=0; i<a.size(); ++i)
  {
    if(abs(a[i]-b[i])>numeric_limits<T>::epsilon()) return false;
  }
  return true;
}

/** operator to check if two buffers are different */
template<class T>
bool operator!=(const PooledData<T>& a, const PooledData<T>& b)
{
  if(a.size() != b.size()) return true;
  //if(a.Current != b.Current) return true;
  for(typename PooledData<T>::size_type i=0; i<a.size(); ++i)
  {
    if(abs(a[i]-b[i])>numeric_limits<T>::epsilon()) return true;
  }
  return false;
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
