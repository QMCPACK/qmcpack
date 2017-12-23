// Copyright (c) 2017, NVIDIA CORPORATION. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of NVIDIA CORPORATION nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <memory>

namespace detail_uninit
{

template<typename T> struct alignment_of_impl;

template<typename T, std::size_t size_diff>
  struct helper
{
  static const std::size_t value = size_diff;
};

template<typename T>
  struct helper<T,0>
{
  static const std::size_t value = alignment_of_impl<T>::value;
};

template<typename T>
  struct alignment_of_impl
{
  struct big { T x; char c; };

  static const std::size_t value = helper<big, sizeof(big) - sizeof(T)>::value;
};
  
} // end detail_uninit

template<typename T>
  struct alignment_of
    : detail_uninit::alignment_of_impl<T>
{};

template<std::size_t Len, std::size_t Align>
  struct aligned_storage
{
  union type
  {
    unsigned char data[Len];
    struct __align__(Align) { } align;
  };
};

template<typename T>
  class uninitialized
{
  private:
    typename aligned_storage<sizeof(T), alignment_of<T>::value>::type storage;

    __host__ __device__ inline const T* ptr() const
    {
      return reinterpret_cast<const T*>(storage.data);
    }

    __host__ __device__ inline T* ptr()
    {
      return reinterpret_cast<T*>(storage.data);
    }

  public:
    // copy assignment
    __host__ __device__ inline uninitialized<T> &operator=(const T &other)
    {
      T& self = *this;
      self = other;
      return *this;
    }

    __host__ __device__ inline T& get()
    {
      return *ptr();
    }

    __host__ __device__ inline const T& get() const
    {
      return *ptr();
    }

    __host__ __device__ inline operator T& ()
    {
      return get();
    }

    __host__ __device__ inline operator const T&() const
    {
      return get();
    }

    inline __host__ __device__ void construct()
    {
      ::new(ptr()) T();
    }

    template<typename Arg>
    inline __host__ __device__ void construct(const Arg &a)
    {
      ::new(ptr()) T(a);
    }

    template<typename Arg1, typename Arg2>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2)
    {
      ::new(ptr()) T(a1,a2);
    }

    template<typename Arg1, typename Arg2, typename Arg3>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3)
    {
      ::new(ptr()) T(a1,a2,a3);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4)
    {
      ::new(ptr()) T(a1,a2,a3,a4);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4, const Arg5 &a5)
    {
      ::new(ptr()) T(a1,a2,a3,a4,a5);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4, const Arg5 &a5, const Arg6 &a6)
    {
      ::new(ptr()) T(a1,a2,a3,a4,a5,a6);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4, const Arg5 &a5, const Arg6 &a6, const Arg7 &a7)
    {
      ::new(ptr()) T(a1,a2,a3,a4,a5,a6,a7);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7, typename Arg8>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4, const Arg5 &a5, const Arg6 &a6, const Arg7 &a7, const Arg8 &a8)
    {
      ::new(ptr()) T(a1,a2,a3,a4,a5,a6,a7,a8);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7, typename Arg8, typename Arg9>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4, const Arg5 &a5, const Arg6 &a6, const Arg7 &a7, const Arg8 &a8, const Arg9 &a9)
    {
      ::new(ptr()) T(a1,a2,a3,a4,a5,a6,a7,a8,a9);
    }

    template<typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7, typename Arg8, typename Arg9, typename Arg10>
    inline __host__ __device__ void construct(const Arg1 &a1, const Arg2 &a2, const Arg3 &a3, const Arg4 &a4, const Arg5 &a5, const Arg6 &a6, const Arg7 &a7, const Arg8 &a8, const Arg9 &a9, const Arg10 &a10)
    {
      ::new(ptr()) T(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10);
    }

    inline __host__ __device__ void destroy()
    {
      T& self = *this;
      self.~T();
    }
};

