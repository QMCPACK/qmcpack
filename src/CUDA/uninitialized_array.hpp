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

#include "uninitialized.hpp"

template<typename T, std::size_t N>
  class uninitialized_array
{
  public:
    typedef T             value_type; 
    typedef T&            reference;
    typedef const T&      const_reference;
    typedef T*            pointer;
    typedef const T*      const_pointer;
    typedef pointer       iterator;
    typedef const_pointer const_iterator;
    typedef std::size_t   size_type;

    __forceinline__ __host__ __device__ iterator begin()
    {
      return data();
    }

    __forceinline__ __host__ __device__ const_iterator begin() const
    {
      return data();
    }

    __forceinline__ __host__ __device__ iterator end()
    {
      return begin() + size();
    }

    __forceinline__ __host__ __device__ const_iterator end() const
    {
      return begin() + size();
    }

    __forceinline__ __host__ __device__ const_iterator cbegin() const
    {
      return begin();
    }

    __forceinline__ __host__ __device__ const_iterator cend() const
    {
      return end();
    }

    __forceinline__ __host__ __device__ size_type size() const
    {
      return N;
    }

    __forceinline__ __host__ __device__ bool empty() const
    {
      return false;
    }

    __forceinline__ __host__ __device__ T* data()
    {
      return impl.get();
    }

    __forceinline__ __host__ __device__ const T* data() const
    {
      return impl.get();
    }

    // element access
    __forceinline__ __host__ __device__ reference operator[](size_type n)
    {
      return data()[n];
    }

    __forceinline__ __host__ __device__ const_reference operator[](size_type n) const
    {
      return data()[n];
    }

    __forceinline__ __host__ __device__ reference front()
    {
      return *data();
    }

    __forceinline__ __host__ __device__ const_reference front() const
    {
      return *data();
    }

    __forceinline__ __host__ __device__ reference back()
    {
      return data()[size() - size_type(1)];
    }

    __forceinline__ __host__ __device__ const_reference back() const
    {
      return data()[size() - size_type(1)];
    }

  private:
    uninitialized<T[N]> impl;
};

