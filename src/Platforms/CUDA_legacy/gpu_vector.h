//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
// Modifications Copyright (C) 2021 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef GPU_VECTOR_H
#define GPU_VECTOR_H

#include <malloc.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#ifdef QMC_CUDA
#ifndef QMC_CUDA2HIP
#include <cuda_runtime_api.h>
#else
#include <hip/hip_runtime.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif
#include "gpu_misc.h"
#endif

namespace gpu
{
struct gpu_mem_object
{
  size_t num_objects;
  size_t total_bytes;
  gpu_mem_object(size_t size) : num_objects(1), total_bytes(size) {}
  gpu_mem_object() : num_objects(0), total_bytes(0) {}
};


class cuda_memory_manager_type
{
private:
  std::map<std::string, gpu_mem_object> gpu_mem_map;
  std::map<void*, std::pair<std::string, size_t>> gpu_pointer_map;

public:
  void* allocate(size_t bytes, std::string name = "");
  void* allocate_managed(size_t bytes, std::string name = "", unsigned int flags = cudaMemAttachGlobal);

  void deallocate(void* p);

  void report();
};

extern cuda_memory_manager_type cuda_memory_manager;

template<typename T>
class host_vector;

template<typename T>
class device_vector
{
private:
  T* data_pointer;
  size_t current_size, alloc_size;
  std::string name;
  // True if the data was allocated by this vector.  False if we're
  // referencing memory
  bool own_data;
  // True if managed memory was requested using resize function, starts out false
  bool managedmem;
  // Flags for managed memory creation (defaults to cudaMemAttachGlobal) that can be set with set_managed_flags function
  unsigned int managed_flags;

public:
  using pointer = T*;

  void set_name(std::string n) { name = n; }

  void set_managed_flags(unsigned int flags) { managed_flags = flags; }

  inline device_vector()
      : data_pointer(NULL),
        current_size(0),
        alloc_size(0),
        own_data(true),
        managedmem(false),
        managed_flags(cudaMemAttachGlobal)
  {}

  inline device_vector(std::string myName)
      : data_pointer(NULL),
        current_size(0),
        alloc_size(0),
        name(myName),
        own_data(true),
        managedmem(false),
        managed_flags(cudaMemAttachGlobal)
  {}

  inline device_vector(size_t size)
      : data_pointer(NULL),
        current_size(0),
        alloc_size(0),
        own_data(true),
        managedmem(false),
        managed_flags(cudaMemAttachGlobal)
  {
    resize(size);
  }

  inline device_vector(std::string myName, size_t size)
      : name(myName),
        data_pointer(NULL),
        current_size(0),
        alloc_size(0),
        own_data(true),
        managedmem(false),
        managed_flags(cudaMemAttachGlobal)
  {
    resize(size);
  }

  inline device_vector(const host_vector<T>& vec);

  ~device_vector()
  {
    if (alloc_size > 0 && data_pointer && own_data)
      cuda_memory_manager.deallocate(data_pointer);
  }

  inline void reference(T* p, size_t size)
  {
    //       fprintf (stderr, "reference called for name=%s size=%ld ptr=%p\n",
    // 	       name.c_str(), size, p);
    if (own_data && alloc_size > 0)
      cuda_memory_manager.deallocate(data_pointer);
    data_pointer = p;
    current_size = size;
    alloc_size   = 0;
    own_data     = false;
    managedmem   = false;
  }


  inline T& operator[](size_t i) const { return data_pointer[i]; }


  inline void resize(size_t size, double reserve_factor = 1.0, bool managed = false)
  {
    if (!size)
    {
      current_size = 0;
      return;
    }
    size_t reserve_size = (size_t)std::ceil(reserve_factor * size);
    size_t byte_size    = sizeof(T) * reserve_size;
    bool error          = false;
    if (managed != managedmem)
    {
      if (managedmem)
      {
        if (alloc_size > 0) // Only trigger error message if memory is allocated
        {
          fprintf(stderr, "device_vector.resize from managed (%p) ", data_pointer);
          error = true;
        }
      }
      else
      {
        if (alloc_size != 0)
          fprintf(stderr, "device_vector.resize from non-managed to managed.\n");
      }
    }
    if ((size > alloc_size) || (alloc_size == 0))
    {
      if (own_data && (alloc_size > 0))
        cuda_memory_manager.deallocate(data_pointer);
      if (managed)
        data_pointer = (T*)cuda_memory_manager.allocate_managed(byte_size, name, managed_flags);
      else
        data_pointer = (T*)cuda_memory_manager.allocate(byte_size, name);
      alloc_size = reserve_size;
      own_data   = true;
      managedmem = managed;
    }
    current_size = size;
    if (error)
      fprintf(stderr, "to non-managed (%p).\n", data_pointer);
  }

  inline void clear()
  {
    if (alloc_size)
    {
      cuda_memory_manager.deallocate(data_pointer);
      data_pointer = NULL;
      current_size = alloc_size = 0;
    }
  }

  inline size_t size() const { return current_size; }


  inline device_vector& operator=(const device_vector<T>& vec)
  {
    if (this->size() != vec.size())
    {
      if (!own_data)
      {
        fprintf(stderr,
                "Assigning referenced GPU vector, but it has the "
                "wrong size.\n");
        fprintf(stderr, "Name = %s.  This size = %ld, vec size = %ld\n", name.c_str(), size(), vec.size());
        abort();
      }
      resize(vec.size(), 1.0, managedmem);
    }
#ifdef QMC_CUDA
    cudaMemcpyAsync(data_pointer, &(vec[0]), this->size() * sizeof(T), cudaMemcpyDeviceToDevice);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error in device_vector::operator=(device_vector):\n  %s\n", cudaGetErrorString(err));
      fprintf(stderr, "vec.size() = %ld\n", vec.size());
      abort();
    }
#endif
    return *this;
  }

  device_vector(const device_vector<T>& vec)
      : data_pointer(NULL),
        current_size(0),
        alloc_size(0),
        name(vec.name),
        own_data(true),
        managedmem(vec.managedmem),
        managed_flags(vec.managed_flags)
  {
    resize(vec.size(), 1.0, managedmem);
    // fprintf(stderr, "device_vector copy constructor called, name=%s.\n",
    // 	      name.c_str());
#ifdef QMC_CUDA
    if (this->size() != 0)
    {
      cudaMemcpy(data_pointer, &(vec[0]), vec.size() * sizeof(T), cudaMemcpyDeviceToDevice);
      cudaError_t err = cudaGetLastError();
      if (err != cudaSuccess)
      {
        fprintf(stderr, "CUDA error in device_vector::copy constructor:\n  %s\n", cudaGetErrorString(err));
        abort();
      }
    }
#endif
  }

  device_vector& operator=(const std::vector<T, std::allocator<T>>& vec)
  {
    if (this->size() != vec.size())
    {
      if (!own_data)
      {
        fprintf(stderr,
                "Assigning referenced GPU vector, but it has the "
                "wrong size.\n");
        // fprintf (stderr, "Name = %s.  This size = %ld, vec size = %ld\n",
        // 	   name.c_str(), size(), vec.size());
        abort();
      }
      resize(vec.size(), 1.0, managedmem);
    }
#ifdef QMC_CUDA
    cudaMemcpyAsync(data_pointer, &(vec[0]), this->size() * sizeof(T), cudaMemcpyHostToDevice);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error in device_vector::operator (%p)=(std::vector %p):\n  %s\n", data_pointer, &(vec[0]),
              cudaGetErrorString(err));
      abort();
    }
#endif
    return *this;
  }

  device_vector& operator=(const host_vector<T>& vec)
  {
    if (this->size() != vec.size())
    {
      if (!own_data)
      {
        fprintf(stderr,
                "Assigning referenced GPU vector, but it has the "
                "wrong size.\n");
        fprintf(stderr, "Name = %s.  This size = %ld, vec size = %ld\n", name.c_str(), size(), vec.size());
        abort();
      }
      resize(vec.size(), 1.0, managedmem);
    }
#ifdef QMC_CUDA
    cudaMemcpy(&((*this)[0]), &(vec[0]), vec.size() * sizeof(T), cudaMemcpyHostToDevice);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "In operator=, name=%s, size=%ld  vec.size()=%ld\n", name.c_str(), size(), vec.size());
      fprintf(stderr, "this pointer = %p  vec pointer=%p\n", data_pointer, &(vec[0]));
      fprintf(stderr, "CUDA error in device_vector::operator=(const host_vector<T> &vec) for %s:\n  %s\n", name.c_str(),
              cudaGetErrorString(err));
      abort();
    }
#endif
    return *this;
  }

  void asyncCopy(const host_vector<T>& vec)
  {
    if (this->size() != vec.size())
    {
      if (!own_data)
      {
        fprintf(stderr,
                "Assigning referenced GPU vector, but it has the "
                "wrong size.\n");
        fprintf(stderr, "Name = %s.  This size = %ld, vec size = %ld\n", name.c_str(), size(), vec.size());
        abort();
      }
      resize(vec.size(), 1.0, managedmem);
    }
#ifdef QMC_CUDA
    cudaMemcpyAsync(&((*this)[0]), &(vec[0]), vec.size() * sizeof(T), cudaMemcpyHostToDevice, kernelStream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "In operator=, name=%s, size=%ld  vec.size()=%ld\n", name.c_str(), size(), vec.size());
      fprintf(stderr, "this pointer = %p  vec pointer=%p\n", data_pointer, &(vec[0]));
      fprintf(stderr, "CUDA error in device_vector::asyncCopy(const host_vector<T> &vec) for %s:\n  %s\n", name.c_str(),
              cudaGetErrorString(err));
      abort();
    }
#endif
  }

  void asyncCopy(const T* vec_ptr, size_t len, size_t offset, size_t datalen)
  {
    if ((this->size() != len) || (this->size() < offset + datalen))
    {
      if (!own_data)
      {
        fprintf(stderr,
                "Assigning referenced GPU vector, but it has the "
                "wrong size.\n");
        fprintf(stderr, "Name = %s.  This size = %ld, vec size = %ld\n", name.c_str(), size(), len);
        abort();
      }
      if (len < offset + datalen)
      {
        fprintf(stderr, "Trying to write more than the length of the vector.\n");
        fprintf(stderr, "Name = %s.  This size = %ld, vec size = %ld, needed length = %ld\n", name.c_str(), size(), len,
                offset + datalen);
        abort();
      }
      resize(len);
    }
#ifdef QMC_CUDA
    cudaMemcpyAsync(&((*this)[offset]), vec_ptr, datalen * sizeof(T), cudaMemcpyHostToDevice, kernelStream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "In operator=, name=%s, size=%ld  vec.size()=%ld\n", name.c_str(), size(), len);
      fprintf(stderr, "this pointer = %p  vec pointer=%p\n", data_pointer, vec_ptr);
      fprintf(stderr, "CUDA error in device_vector::asyncCopy(const T* vec_ptr, len, offset, datalen) for %s:\n  %s\n",
              name.c_str(), cudaGetErrorString(err));
      abort();
    }
#endif
  }

  void asyncCopy(const std::vector<T, std::allocator<T>>& vec)
  {
    if (this->size() != vec.size())
    {
      if (!own_data)
      {
        fprintf(stderr,
                "Assigning referenced GPU vector, but it has the "
                "wrong size.\n");
        fprintf(stderr, "Name = %s.  This size = %ld, vec size = %ld\n", name.c_str(), size(), vec.size());
        abort();
      }
      resize(vec.size(), 1.0, managedmem);
    }
#ifdef QMC_CUDA
    cudaMemcpyAsync(&((*this)[0]), &(vec[0]), vec.size() * sizeof(T), cudaMemcpyHostToDevice, kernelStream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "In operator=, name=%s, size=%ld  vec.size()=%ld\n", name.c_str(), size(), vec.size());
      fprintf(stderr, "this pointer = %p  vec pointer=%p\n", data_pointer, &(vec[0]));
      fprintf(stderr,
              "CUDA error in device_vector::asyncCopy(const std::vector<T, std::allocator<T> > &vec) for %s:\n  %s\n",
              name.c_str(), cudaGetErrorString(err));
      abort();
    }
#endif
  }
  void copyFromGPU(std::vector<T, std::allocator<T>>& vec)
  {
    if (this->size() != vec.size())
    {
      vec.resize(size());
    }
#ifdef QMC_CUDA
    cudaMemcpy(&(vec[0]), &((*this)[0]), vec.size() * sizeof(T), cudaMemcpyDeviceToHost);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "In operator=, name=%s, size=%ld  vec.size()=%ld\n", name.c_str(), size(), vec.size());
      fprintf(stderr, "this pointer = %p  vec pointer=%p\n", data_pointer, &(vec[0]));
      fprintf(stderr,
              "CUDA error in device_vector::copyFromGPU(std::vector<T, std::allocator<T> > &vec) for %s:\n  %s\n",
              name.c_str(), cudaGetErrorString(err));
      abort();
    }
#endif
  }


  inline T* data() const { return data_pointer; }
};


template<typename T>
class host_vector
{
private:
  T* data;
  size_t current_size;
  size_t capacity;

public:
  host_vector()
  {
    data         = NULL;
    current_size = 0;
    capacity     = 0;
  }

  host_vector(const host_vector<T>& vec)
  {
    if (vec.size() != 0)
    {
      cudaHostAlloc((void**)&data, vec.size() * sizeof(T), 0);
      cudaMemcpy(data, vec.data, vec.size() * sizeof(T), cudaMemcpyHostToHost);
    }
    else
    {
      data = NULL;
    }
    current_size = vec.size();
    capacity     = current_size;
  }

  host_vector(int size)
  {
    data         = NULL;
    current_size = 0;
    capacity     = 0;
    resize(size);
  }

  ~host_vector()
  {
    if (data)
    {
      cudaFreeHost(data);
      data         = NULL;
      current_size = 0;
      capacity     = 0;
    }
  }

  host_vector(const device_vector<T>& vec)
  {
    data         = NULL;
    current_size = 0;
    capacity     = 0;
    resize(vec.size());
#ifdef QMC_CUDA
    cudaMemcpy(&(data[0]), &(vec[0]), current_size * sizeof(T), cudaMemcpyDeviceToHost);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error in host_vector::copy constructor():\n  %s\n", cudaGetErrorString(err));
      abort();
    }
#endif
  }


  host_vector& operator=(const host_vector<T>& vec)
  {
    if (this->size() != vec.size())
      this->resize(vec.size());
#ifdef QMC_CUDA
    cudaMemcpyAsync(&((*this)[0]), &(vec[0]), this->size() * sizeof(T), cudaMemcpyHostToDevice);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error in host_vector::operator=(host_vector):\n  %s\n", cudaGetErrorString(err));
      abort();
    }
#endif
    return *this;
  }

  host_vector& operator=(const device_vector<T>& vec)
  {
    if (this->size() != vec.size())
      this->resize(vec.size());
#ifdef QMC_CUDA
    cudaMemcpy(&((*this)[0]), &(vec[0]), this->size() * sizeof(T), cudaMemcpyDeviceToHost);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error in host_vector::operator=(device_vector %p):\n  %s\n", &(vec[0]),
              cudaGetErrorString(err));
      abort();
    }
#endif
    return *this;
  }

  void asyncCopy(const device_vector<T>& vec)
  {
    if (this->size() != vec.size())
      resize(vec.size());
#ifdef QMC_CUDA
    cudaMemcpyAsync(&((*this)[0]), &(vec[0]), this->size() * sizeof(T), cudaMemcpyDeviceToHost, memoryStream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA error in host_vector::asyncCopy:\n  %s\n", cudaGetErrorString(err));
      abort();
    }
#endif
  }

  inline size_t size() const { return current_size; }

  void reserve(size_t new_size)
  {
    if (new_size <= capacity)
      return;
    T* new_data;
    // QMCPACK often does repeated resizes like 256->257 then 257->258.
    // this anticipates the further resizes by pre-allocating an additional
    // 5% above what was requested.
    new_size = 1.05 * new_size;
    cudaHostAlloc((void**)&new_data, new_size * sizeof(T), 0);
    if (data != NULL)
    {
      cudaMemcpy(new_data, data, current_size * sizeof(T), cudaMemcpyHostToHost);
      cudaFreeHost(data);
      data = NULL;
    }
    data     = new_data;
    capacity = new_size;
  }

  T& operator[](const int n) { return data[n]; }
  const T& operator[](const int n) const
  {
    const T& a = data[n];
    return a;
  }

  inline void resize(size_t new_size)
  {
    if (new_size <= current_size)
    {
      current_size = new_size;
      if (new_size == 0)
      {
        clear();
      }
      return;
    }
    reserve(new_size);
    for (int i = current_size; i < new_size; ++i)
      data[i] = T();
    current_size = new_size;
  }

  inline void clear()
  {
    if (data != NULL)
    {
      //cudaFreeHost(data);
      //data = NULL;
      current_size = 0;
      //capacity = 0;
    }
  }

  inline void push_back(const T& x)
  {
    if (current_size < capacity)
    {
      data[current_size] = x;
      ++current_size;
      return;
    }
    reserve(2 * capacity + 1);
    push_back(x);
  }
};

template<typename T>
device_vector<T>::device_vector(const host_vector<T>& vec)
    : data_pointer(NULL), current_size(0), alloc_size(0), own_data(true), managedmem(false)
{
  this->resize(vec.size());
#ifdef QMC_CUDA
  cudaMemcpy(&((*this)[0]), &(vec[0]), this->size() * sizeof(T), cudaMemcpyDeviceToHost);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in host_vector::operator=() for %s:\n  %s\n", name.c_str(), cudaGetErrorString(err));
    abort();
  }
#endif
}

template<typename T>
class device_host_vector
{
public:
  host_vector<T> CPU;
  device_vector<T> GPU;

  device_host_vector() {}

  device_host_vector(size_t size) : CPU(size), GPU(size) {}

  inline void resize(size_t size)
  {
    CPU.resize(size);
    GPU.resize(size);
  }

  inline void host_to_device() { GPU = CPU; }

  inline void device_to_host() { CPU = GPU; }

  inline T operator[](size_t i) const { return CPU[i]; }

  inline T& operator[](size_t i) { return CPU[i]; }

  inline void fill(T val)
  {
    std::fill(CPU.begin(), CPU.end(), val);
    host_to_device();
  }

  inline T* gpu_data() { return GPU.data(); }
};
} // namespace gpu

#endif
