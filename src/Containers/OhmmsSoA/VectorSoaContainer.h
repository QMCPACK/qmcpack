//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file VectorSoaContainer.h
 * Soa Container for D-dim vectors
 *
 * Alternative to Container<TinyVector<T,D> > to support SoA algorithms
 */
#ifndef QMCPLUSPLUS_VECTOR_SOA_H
#define QMCPLUSPLUS_VECTOR_SOA_H

#include <type_traits>
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/algorithm.hpp"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsSoA/PosTransformer.h"

namespace qmcplusplus
{
/** SoA adaptor class for Vector<TinyVector<T,D> >
 * @tparm T data type, float, double, complex<float>, complex<double>
 * @tparm Alloc memory allocator
 */
template<typename T, unsigned D, typename Alloc = aligned_allocator<T>>
struct VectorSoaContainer
{
  using AoSElement_t = TinyVector<T, D>;
  using Element_t    = T;

  ///default constructor
  VectorSoaContainer() : nLocal(0), nGhosts(0), nAllocated(0), myData(nullptr) {}

  /** constructor for a view, attach to pre-allocated data
   * @param ptr new myData
   * @param n new nLocal
   * @param n_padded new nGhosts
   */
  VectorSoaContainer(T* ptr, size_t n, size_t n_padded) : nLocal(n), nGhosts(n_padded), nAllocated(0), myData(ptr) {}

  ///destructor
  ~VectorSoaContainer() { free(); }

  ///default copy constructor
  VectorSoaContainer(const VectorSoaContainer& in) : nLocal(0), nGhosts(0), nAllocated(0), myData(nullptr)
  {
    resize(in.nLocal);
    std::copy_n(in.myData, nGhosts * D, myData);
  }

  ///default copy operator
  VectorSoaContainer& operator=(const VectorSoaContainer& in)
  {
    if (myData != in.myData)
    {
      resize(in.nLocal);
      std::copy_n(in.myData, nGhosts * D, myData);
    }
    return *this;
  }

  ///move constructor
  VectorSoaContainer(VectorSoaContainer&& in) noexcept
      : nLocal(in.nLocal), nGhosts(in.nGhosts), nAllocated(in.nAllocated), myData(std::move(in.myData))
  {
    in.myData     = nullptr;
    in.nAllocated = 0;
    in.nLocal     = 0;
    in.nGhosts    = 0;
  }

  /// constructor with size n without initialization
  explicit VectorSoaContainer(size_t n) : nLocal(0), nGhosts(0), nAllocated(0), myData(nullptr) { resize(n); }

  /** constructor with Vector<T1,D> */
  template<typename T1>
  VectorSoaContainer(const Vector<TinyVector<T1, D>>& in)
      : nLocal(0), nGhosts(0), nAllocated(0), myData(nullptr)
  {
    resize(in.size());
    copyIn(in);
  }

  template<typename T1>
  VectorSoaContainer& operator=(const Vector<TinyVector<T1, D>>& in)
  {
    resize(in.size());
    copyIn(in);
    return *this;
  }

  /** need A=0.0;
       */
  template<typename T1>
  VectorSoaContainer& operator=(T1 in)
  {
    std::fill(myData, myData + nGhosts * D, static_cast<T>(in));
    return *this;
  }

  /** resize myData
   * @param n nLocal
   *
   * nAllocated is used to ensure no memory leak
   */
  __forceinline void resize(size_t n)
  {
    static_assert(std::is_same<Element_t, typename Alloc::value_type>::value,
                  "VectorSoaContainer and Alloc data types must agree!");
    if (isRefAttached())
      throw std::runtime_error("Resize not allowed on VectorSoaContainer constructed by initialized memory.");

    size_t n_padded = getAlignedSize<T, Alloc::alignment>(n);

    if (n_padded * D > nAllocated)
    {
      if (nAllocated)
        mAllocator.deallocate(myData, nAllocated);
      nLocal     = n;
      nGhosts    = n_padded;
      nAllocated = nGhosts * D;
      myData     = mAllocator.allocate(nAllocated);
    }
    else
    {
      nLocal  = n;
      nGhosts = n_padded;
    }
  }

  ///clear
  inline void clear()
  {
    nLocal  = 0;
    nGhosts = 0;
  }

  /// free allocated memory and clear status variables
  __forceinline void free()
  {
    if (nAllocated)
      mAllocator.deallocate(myData, nAllocated);
    nLocal     = 0;
    nGhosts    = 0;
    nAllocated = 0;
    myData     = nullptr;
  }

  /** attach to pre-allocated data
   * @param n new nLocal
   * @param n_padded new nGhosts
   * @param ptr new myData
   *
   * To attach to existing memory, currently owned memory must be freed before calling attachReference
   */
  __forceinline void attachReference(size_t n, size_t n_padded, T* ptr)
  {
    if (nAllocated)
    {
      free();
      // This is too noisy right now.
      // std::cerr << "OhmmsVectorSoa attachReference called on previously allocated vector.\n" << std::endl;
      /// \todo return this when buffer system is simplified.
    }
    nLocal  = n;
    nGhosts = n_padded;
    myData  = ptr;
  }

  ///return the physical size
  __forceinline size_t size() const { return nLocal; }
  ///return the physical size
  __forceinline size_t capacity() const { return nGhosts; }

  /** AoS to SoA : copy from Vector<TinyVector<>>
       *
       * The same sizes are assumed.
       */
  template<typename T1>
  void copyIn(const Vector<TinyVector<T1, D>>& in)
  {
    //if(nLocal!=in.size()) resize(in.size());
    PosAoS2SoA(nLocal, D, reinterpret_cast<const T1*>(in.first_address()), D, myData, nGhosts);
  }

  /** SoA to AoS : copy to Vector<TinyVector<>>
       *
       * The same sizes are assumed.
       */
  template<typename T1>
  void copyOut(Vector<TinyVector<T1, D>>& out) const
  {
    PosSoA2AoS(nLocal, D, myData, nGhosts, reinterpret_cast<T1*>(out.first_address()), D);
  }

  /** return TinyVector<T,D>
       */
  __forceinline const AoSElement_t operator[](size_t i) const { return AoSElement_t(myData + i, nGhosts); }

  ///helper class for operator ()(size_t i) to assign a value
  struct Accessor
  {
    T* _base;
    size_t M;
    __forceinline Accessor(T* a, size_t ng) : _base(a), M(ng) {}
    template<typename T1>
    __forceinline Accessor& operator=(const TinyVector<T1, D>& rhs)
    {
#pragma unroll
      for (size_t i = 0; i < D; ++i)
        *(_base + M * i) = rhs[i];
      return *this;
    }

    /** assign value */
    template<typename T1>
    __forceinline Accessor& operator=(T1 rhs)
    {
#pragma unroll
      for (size_t i = 0; i < D; ++i)
        *(_base + M * i) = rhs;
      return *this;
    }
  };

  /** access operator for assignment of the i-th value
       *
       * Use for (*this)[i]=TinyVector<T,D>;
       */
  __forceinline Accessor operator()(size_t i) { return Accessor(myData + i, nGhosts); }
  ///return the base
  __forceinline T* data() { return myData; }
  ///return the base
  __forceinline const T* data() const { return myData; }
  /// return non_const data
  T* getNonConstData() const { return myData; }
  ///return the pointer of the i-th components
  __forceinline T* restrict data(size_t i) { return myData + i * nGhosts; }
  ///return the const pointer of the i-th components
  __forceinline const T* restrict data(size_t i) const { return myData + i * nGhosts; }
  ///return the end
  __forceinline T* end() { return myData + D * nGhosts; }
  ///return the end
  __forceinline const T* end() const { return myData + D * nGhosts; }


  ///return the base, device
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  __forceinline T* device_data() { return mAllocator.getDevicePtr(); }
  ///return the base, device
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  __forceinline const T* device_data() const { return mAllocator.getDevicePtr(); }
  ///return the pointer of the i-th components, device
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  __forceinline T* restrict device_data(size_t i) { return mAllocator.getDevicePtr() + i * nGhosts; }
  ///return the const pointer of the i-th components, device
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  __forceinline const T* restrict device_data(size_t i) const { return mAllocator.getDevicePtr() + i * nGhosts; }

private:
  /// number of elements
  size_t nLocal;
  /// number of elements + padded
  size_t nGhosts;
  /// number of elements allocated by myAlloc
  size_t nAllocated;
  /// pointer: what type????
  T* myData;
  /// allocator
  Alloc mAllocator;

  /// return true if memory is not owned by the container but from outside.
  inline bool isRefAttached() const { return nGhosts * D > nAllocated; }

  // We need this because of the hack propagation of the VectorSoaContainer allocator
  // to allow proper OhmmsVector based views of VectorSoAContainer elements with OMPallocator
  friend class qmcplusplus::Vector<T, Alloc>;
};

} // namespace qmcplusplus

#endif
