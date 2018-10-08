//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_POOLEDMEMORY_H
#define QMCPLUSPLUS_POOLEDMEMORY_H

#include <complex>
#include <cstring>
#include "simd/allocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

/** Memory pool to manage arrays+scalars with alignment
 * @tparam T_scalar : type of the scalar data type, typically OHMMS_PRECISION_FULL
 * @tparam PageSize : page size in bytes, default=4096 (4K)
 * @tparam Alloc : allocator, Mallocator<T, PageSize>
 *
 * The bulk part is accessed directly by address
 * The scalar part works as PooledData, all the values are static_cast to T_scalar.
 */
template<typename T_scalar=OHMMS_PRECISION_FULL,
         size_t PageSize=4096,
         typename Alloc=Mallocator<char, PageSize> >
struct PooledMemory
{
  typedef char T;
  typedef T value_type;
  typedef typename Vector<T, Alloc>::size_type size_type;

  const int scalar_multiplier;
  size_type Current, Current_scalar;
  T_scalar *Scalar_ptr;
  Vector<T, Alloc> myData;

  ///default constructor
  inline PooledMemory():
    scalar_multiplier(sizeof(T_scalar)),
    Current(0), Scalar_ptr(nullptr), Current_scalar(0)
  {
  }

  ///copy constructor
  PooledMemory(const PooledMemory &in):
    scalar_multiplier(in.scalar_multiplier),
    myData(in.myData), Current(in.Current),
    Current_scalar(in.Current_scalar)
  {
    Scalar_ptr = reinterpret_cast<T_scalar*>(myData.data()+in.scalar_offset());
  }

  ///copy assign operator
  PooledMemory& operator=(const PooledMemory &in)
  {
    myData = in.myData;
    Current = in.Current;
    Current_scalar = in.Current_scalar;
    Scalar_ptr = reinterpret_cast<T_scalar*>(myData.data()+in.scalar_offset());
    return *this;
  }

  ///return the size of the data
  inline size_type byteSize() const
  {
    return sizeof(T)*myData.size();
  }

  ///return the size of the data
  inline size_type size() const
  {
    return myData.size();
  }

  //return the cursor
  inline size_type current() const
  {
    return Current;
  }

  //return the cursor of scalar
  inline size_type current_scalar() const
  {
    return Current_scalar;
  }

  /** set the cursors
   * @param cur locator to which Current is assigned
   * @param cur_scalar locator to which Current_scalar is assigned
   */
  inline void rewind(size_type cur=0, size_type cur_scalar=0)
  {
    Current = cur;
    Current_scalar = cur_scalar;
  }

  ///clear the data and set Current=0
  inline void clear()
  {
    myData.clear();
    Current=0;
    Scalar_ptr=nullptr;
    Current_scalar=0;
  }

  ///allocate the data
  inline void allocate()
  {
    myData.resize(Current+Current_scalar*scalar_multiplier);
    Scalar_ptr = reinterpret_cast<T_scalar*>(myData.data()+Current);
  }

  template<typename T1>
  inline void add(std::complex<T1>& x)
  {
    Current_scalar+=2;
  }

  template<typename T1>
  inline void add(T1& x)
  {
    Current_scalar++;
  }

  template<typename T1>
  inline void add(T1* first, T1* last)
  {
    constexpr int multiplier=sizeof(T1);
    Current += getAlignedSize<T>((last-first)*multiplier);
  }

  template<typename T1>
  inline void get(std::complex<T1>& x)
  {
    x = std::complex<T1>(static_cast<T1>(Scalar_ptr[Current_scalar]),
                         static_cast<T1>(Scalar_ptr[Current_scalar+1]));
    Current_scalar+=2;
  }

  template<typename T1>
  inline void get(T1& x)
  {
    x = static_cast<T1>(Scalar_ptr[Current_scalar++]);
  }

  template<typename T1>
  inline void get(T1* first, T1* last)
  {
    // for backward compatibility
    const size_t nbytes=(last-first)*sizeof(T1);
    std::memcpy(first, myData.data()+Current, nbytes);
    Current += getAlignedSize<T>(nbytes);
  }

  template<typename T1>
  inline T1* lendReference(size_type n)
  {
    constexpr int multiplier=sizeof(T1);
    T1 *ptr = reinterpret_cast<T1*>(myData.data()+Current);
    Current += getAlignedSize<T>(n*multiplier);
    return ptr;
  }

  inline void forward(size_type n)
  {
    Current += getAlignedSize<T>(n);
  }

  template<typename T1>
  inline void put(std::complex<T1>& x)
  {
    Scalar_ptr[Current_scalar++] = static_cast<T_scalar>(x.real());
    Scalar_ptr[Current_scalar++] = static_cast<T_scalar>(x.imag());
  }

  template<typename T1>
  inline void put(T1& x)
  {
    Scalar_ptr[Current_scalar++] = static_cast<T_scalar>(x);
  }

  template<typename T1>
  inline void put(T1* first, T1* last)
  {
    // for backward compatibility
    const size_t nbytes=(last-first)*sizeof(T1);
    std::memcpy(myData.data()+Current, first, nbytes);
    Current += getAlignedSize<T>(nbytes);
  }

  /** return the address of the first element **/
  inline T* data()
  {
    return myData.data();
  }

  /** return the address offset of the first scalar element **/
  inline size_type scalar_offset() const
  {
    return reinterpret_cast<T*>(Scalar_ptr) - myData.data();
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m)
  {
    m.Pack(myData.data(),myData.size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m)
  {
    m.Unpack(myData.data(),myData.size());
    return m;
  }

  template<typename T1>
  inline PooledMemory& operator<<(T1 &x)
  {
    put(x);
    return *this;
  }

  template<typename T1>
  inline PooledMemory& operator>>(T1 &x)
  {
    get(x);
    return *this;
  }

};

}
#endif
