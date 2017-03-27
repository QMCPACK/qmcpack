//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file PooledData.h
 * @brief Define a serialized buffer to store anonymous data
 *
 * JK: Remove iterator version on 2016-01-04
 */
#ifndef QMCPLUSPLUS_POOLEDDATA_H
#define QMCPLUSPLUS_POOLEDDATA_H
#include <vector>
#include <complex>
#include <limits>
#include <iterator>

using std::complex;

template<class T>
struct PooledData
{
  typedef T value_type;
  typedef OHMMS_PRECISION_FULL fp_value_type;
  typedef typename std::vector<T>::size_type size_type;

  size_type Current, Current_DP;
  std::vector<T> myData;
  /// only active when T!=fp_value_type
  std::vector<fp_value_type> myData_DP;

  ///default constructor
  inline PooledData(): Current(0), Current_DP(0) { }

  ///constructor with a size
  inline PooledData(size_type n, size_type n1=0): Current(0), Current_DP(0)
  {
    myData.resize(n,T());
    myData_DP.resize(n1,fp_value_type());
  }

  ///return the size of the data
  inline size_type byteSize() const
  {
    return sizeof(T)*myData.size() + sizeof(fp_value_type)*myData_DP.size();
  }

  ///return the size of the data
  inline size_type size() const
  {
    return myData.size();
  }

  ///return the size of the DP data
  inline size_type size_DP() const
  {
    return myData_DP.size();
  }

  //inline void clear() { std::vector<T>::clear(); Current=0;}
  inline size_type current() const
  {
    return Current;
  }

   inline size_type current_DP() const
  {
    return Current_DP;
  }

  /** set the Current to a cursor
   * @param cur locator to which Current is assigned
   */
  inline void rewind(size_type cur=0, size_type cur_DP=0)
  {
    Current = cur;
    Current_DP = cur_DP;
  }

  ///return the starting iterator
  inline typename std::vector<T>::iterator begin()
  {
    return myData.begin();
  }
  ///return the ending iterator
  inline typename std::vector<T>::iterator end()
  {
    return myData.end();
  }

  /*@{ matching functions to std::vector functions */
  ///clear the data and set Current=0
  inline void clear()
  {
    myData.clear();
    myData_DP.clear();
    Current=0;
    Current_DP=0;
  }
  ///reserve the memory using std::vector<T>::reserve
  inline void reserve(size_type n)
  {
    myData.reserve(n);
    Current=0;
  }
  ///resize
  inline void resize(size_type n, T val=T())
  {
    myData.resize(n,val);
    Current=0;
  }
  ///return i-th value
  inline T operator[](size_type i) const
  {
    return myData[i];
  }
  ///return i-th value to assign
  inline T& operator[](size_type i)
  {
    return myData[i];
  }
  /*@}*/

  template<class T1>
  inline void add(T1 x)
  {
    Current++;
    myData.push_back(static_cast<T>(x));
  }

  inline void add(T x)
  {
    Current++;
    myData.push_back(x);
  }

  inline void add(std::complex<T>& x)
  {
    Current+=2;
    myData.push_back(x.real());
    myData.push_back(x.imag());
  }

  template<class _InputIterator>
  inline void add(_InputIterator first, _InputIterator last)
  {
    Current += last-first;
    myData.insert(myData.end(),first,last);
  }

  inline void add(T* first, T* last)
  {
    Current += last-first;
    myData.insert(myData.end(),first,last);
  }

  template<class T1>
  inline void add(T1* first, T1* last)
  {
    Current_DP += last-first;
    myData_DP.insert(myData_DP.end(),first,last);
  }

  inline void add(std::complex<T>* first,std::complex<T>* last)
  {
    size_type dn=2*(last-first);
    T* t=reinterpret_cast<T*>(first);
    myData.insert(myData.end(),t,t+dn);
    Current += dn;
  }

  template<typename T1>
  inline void add(std::complex<T1>* first,std::complex<T1>* last)
  {
    size_type dn=2*(last-first);
    T1* t=reinterpret_cast<T1*>(first);
    myData_DP.insert(myData_DP.end(),t,t+dn);
    Current_DP += dn;
  }

  template<class T1>
  inline void get(T1& x)
  {
    x = static_cast<T1>(myData[Current++]);
  }

  inline void get(T& x)
  {
    x = myData[Current++];
  }

  inline void get(std::complex<T>& x)
  {
    x=std::complex<T>(myData[Current],myData[Current+1]);
    Current+=2;
  }

  template<class _OutputIterator>
  inline void get(_OutputIterator first, _OutputIterator last)
  {
    size_type now=Current;
    Current+=last-first;
    copy(myData.begin()+now,myData.begin()+Current,first);
  }

  inline void get(T* first, T* last)
  {
    size_type now=Current;
    Current+=last-first;
    std::copy(myData.begin()+now,myData.begin()+Current,first);
  }

  template<class T1>
  inline void get(T1* first, T1* last)
  {
    size_type now=Current_DP;
    Current_DP+=last-first;
    std::copy(myData_DP.begin()+now,myData_DP.begin()+Current_DP,first);
  }

  inline void get(std::complex<T>* first, std::complex<T>* last)
  {
    while(first != last)
    {
      (*first)=std::complex<T>(myData[Current],myData[Current+1]);
      ++first;
      Current+=2;
    }
  }

  template<typename T1>
  inline void get(std::complex<T1>* first, std::complex<T1>* last)
  {
    while(first != last)
    {
      (*first)=std::complex<T1>(myData_DP[Current_DP],myData_DP[Current_DP+1]);
      ++first;
      Current_DP+=2;
    }
  }

  inline void put(T x)
  {
    myData[Current++] = x;
  }
  inline void put(std::complex<T>& x)
  {
    myData[Current++] = x.real();
    myData[Current++] = x.imag();
  }

  template<class _InputIterator>
  inline void put(_InputIterator first, _InputIterator last)
  {
    copy(first,last,myData.begin()+Current);
    Current+=last-first;
  }

  inline void put(T* first, T* last)
  {
    std::copy(first,last,myData.begin()+Current);
    Current+=last-first;
  }

  template<class T1>
  inline void put(T1* first, T1* last)
  {
    std::copy(first,last,myData_DP.begin()+Current_DP);
    Current_DP+=last-first;
  }

  inline void put(std::complex<T>* first, std::complex<T>* last)
  {
    while(first != last)
    {
      myData[Current++] = (*first).real();
      myData[Current++] = (*first).imag();
      ++first;
    }
  }

  template<typename T1>
  inline void put(std::complex<T1>* first, std::complex<T1>* last)
  {
    while(first != last)
    {
      myData_DP[Current_DP++] = (*first).real();
      myData_DP[Current_DP++] = (*first).imag();
      ++first;
    }
  }


  /** return the address of the first element **/
  inline T* data()
  {
    return &(myData[0]);
  }

  /** return the address of the first DP element **/
  inline fp_value_type* data_DP()
  {
    return &(myData_DP[0]);
  }

  inline void print(std::ostream& os)
  {
    copy(myData.begin(), myData.end(), std::ostream_iterator<T>(os," "));
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m)
  {
    m.Pack(&(myData[0]),myData.size());
    m.Pack(&(myData_DP[0]),myData_DP.size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m)
  {
    m.Unpack(&(myData[0]),myData.size());
    m.Unpack(&(myData_DP[0]),myData_DP.size());
    return m;
  }

  inline PooledData<T>& operator += (const PooledData<T>& s)
  {
    for(int i=0; i<myData.size(); ++i)
      myData[i] += s[i];
    return *this;
  }
  inline PooledData<T>& operator *= (T scale)
  {
    for(int i=0; i<myData.size(); ++i)
      myData[i] *= scale;
    return *this;
  }
};

/** operator to check if two buffers are identical */
template<class T>
bool operator==(const PooledData<T>& a, const PooledData<T>& b)
{
  if(a.size() != b.size())
    return false;
  //if(a.Current != b.Current) return false;
  for(typename PooledData<T>::size_type i=0; i<a.size(); ++i)
  {
    if( std::abs(a[i]-b[i]) > std::numeric_limits<T>::epsilon() )
      return false;
  }
  return true;
}

/** operator to check if two buffers are different */
template<class T>
bool operator!=(const PooledData<T>& a, const PooledData<T>& b)
{
  if(a.size() != b.size())
    return true;
  //if(a.Current != b.Current) return true;
  for(typename PooledData<T>::size_type i=0; i<a.size(); ++i)
  {
    if( std::abs(a[i]-b[i]) > std::numeric_limits<T>::epsilon() )
      return true;
  }
  return false;
}
#endif
