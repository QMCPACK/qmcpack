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
    
    



#ifndef QMCPLUSPLUS_WALKER_ALIGNED_ALLOCATOR_H
#define QMCPLUSPLUS_WALKER_ALIGNED_ALLOCATOR_H

#include <vector>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <OhmmsPETE/OhmmsVectorRef.h>

namespace qmcplusplus
{

#if defined(__INTEL_COMPILER)
  template<unsigned ALIGN>
  struct MMMalloc
  {
    static void* allocate(size_t n)
    {
      return  _mm_malloc(n,ALIGN);
    }

    static void deallocate(void*  ptr)
    {
      _mm_free(ptr);
    }
  };

  template<unsigned ALIGN> using StdAlloc=MMMalloc<ALIGN>;
#else
  template<unsigned ALIGN>
    struct StdAlloc
    {
      static void* allocate(size_t n)
      {
        return std::aligned_alloc(ALIGN,n);
      }

      static void deallocate(void*  ptr)
      {
        std::free(ptr);
      }
    };
#endif

  /** Allocator to manage scalars+arrays
   * @tparam T : type of the primary data type, typically OHMMS_PRECISION
   * @tparam PageSize : page size in bytes, default=4096
   * @tparam CacheLieSize : cacheline size in bytes, default=64
   * @tparam Alloc : allocator
   *
   */
  template<typename T, 
    size_t PageSize=4096, 
    size_t CacheLineSize=64,
    typename Alloc=StdAlloc<PageSize> >
  struct WalkerAllocator
  {
    ///counter for bulk
    size_t BulkCount=0; 
    ///counter for scalar
    size_t ScalarCount=0; 
    ///total bytes allocated 
    size_t TotCapacity=0;
    ///allocated pool: unique pointer
    void* m_pool=nullptr; 
    ///bulk data
    VectorRef<T> m_bulk;
    ///scalars all in double
    VectorRef<double> m_scalars;

    /** copy constructor */
    WalkerAllocator(const WalkerAllocator& rhs) : TotCapacity(rhs.TotCapacity)
    {
      if(TotCapacity>0)
        m_pool=Alloc::allocate(TotCapacity);
    }

    ~WalkerAllocator()
    {
      if(TotCapacity>0) Alloc::deallocate(m_pool);
    }

    inline size_t size() const { return TotCapacity;}
    //inline const T* data() const { return static_cast<const T*>(m_pool);}
    inline T* data() { return m_bulk.data();}
    inline T* address(size_t offset) {return m_bulk.data()+offset;}

    inline void reset() 
    {
      BulkCount=0;
      ScalarCount=0;
    }

    //return the starting address of n items
    inline size_t append(size_t n)
    {
      size_t cur=BulkCount;
      BulkCount += simd::getAligedSize<T>(n);
      return cur;
    }

    //add a scalar, no assinment, all converted into double
    template<typanem T1>
    inline size_t add(T1 a)
    {
      return ScalarCount++;
    }

    template<typanem T1>
    inline void put(T1 a)
    {
      m_scalars[ScalarCount++]=a;
    }

    template<typanem T1>
    inline void get(T1& a) const
    {
      a=m_scalars[ScalarCount++];
    }

    inline size_t add(const std::complex<T>& a)
    {
      ScalarCount+=2;
      return ScalarCount;
    }

    inline void put(const std::complex<T>& a)
    {
      m_scalars[ScalarCount++]=a.real();
      m_scalars[ScalarCount++]=a.imag();
    }

    inline void get(std::complex<T>& a) 
    {
      a=std::complex<T>(m_scalars[ScalarCount],m_scalars[ScalarCount+1]);
      ScalarCount+=2;
    }

    ///return the bytesize with alignment
    inline size_t getAlignedBytes(size_t n, size_t AL)
    {
      return ((n+AL-1)/AL)*AL;
    }

    inline void allocate()
    {
      //this cannot happen but just in case
      if(m_pool!=nullptr) Alloc::deallocate(m_pool);
      TotCapacity=getAlignedBytes(BulkCount*sizeof(T)+m_scalars.size()*sizeof(double),PageSize);
      m_pool=Alloc::allocate(TotCapacity);
      m_bulk.set(static_cast<T*>(m_pool));
      m_scalars.set(static_cast<double*>(m_pool+BulkCount*sizeof(T)));
    }

  };
}

#endif
