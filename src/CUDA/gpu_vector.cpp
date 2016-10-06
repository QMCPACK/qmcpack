//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "config.h"
#include "gpu_vector.h"
#include <cassert>


namespace gpu
{
cuda_memory_manager_type cuda_memory_manager;

void*
cuda_memory_manager_type::allocate(size_t bytes, std::string name)
{
  //    fprintf (stderr, "Allocating %ld bytes for %s\n", bytes, name.c_str());
  // Make sure size is a multiple of 16
  bytes = (((bytes+15)/16)*16);
  void *p;
  cudaMalloc ((void**)&p, bytes);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf (stderr, "Failed to allocate %ld bytes for object %s:\n  %s\n",
             bytes, name.c_str(), cudaGetErrorString(err));
    abort();
  }
  else
  {
    gpu_pointer_map[p] = std::pair<std::string,size_t> (name,bytes);
    std::map<std::string,gpu_mem_object>::iterator iter = gpu_mem_map.find(name);
    if (iter == gpu_mem_map.end())
    {
      gpu_mem_object obj(bytes);
      gpu_mem_map[name] = obj;
    }
    else
    {
      gpu_mem_object &obj = iter->second;
      obj.num_objects++;
      obj.total_bytes += bytes;
    }
  }
  return p;
}

void
cuda_memory_manager_type::deallocate (void *p)
{
  std::map<void*,std::pair<std::string,size_t> >::iterator piter;
  piter = gpu_pointer_map.find(p);
  if (piter == gpu_pointer_map.end())
  {
    fprintf (stderr, "Attempt to free a GPU pointer not in the map.\n");
    abort();
  }
  std::string name = piter->second.first;
  size_t bytes = piter->second.second;
  // fprintf (stderr, "Deallocating %s from GPU memory of size %ld at pointer %p.\n",
  //  	     name.c_str(), bytes, p);
  gpu_pointer_map.erase(piter);
  cudaFree (p);
  std::map<std::string,gpu_mem_object>::iterator iter = gpu_mem_map.find(name);
  if (iter == gpu_mem_map.end())
  {
    fprintf (stderr, "Error:  CUDA memory object %s not found "
             "in the memory map\n", name.c_str());
    abort();
  }
  gpu_mem_object &obj = iter->second;
  if (obj.num_objects == 1)
  {
    assert (bytes == obj.total_bytes);
    gpu_mem_map.erase(iter);
  }
  else
  {
    obj.num_objects--;
    obj.total_bytes -= bytes;
  }
}

void
cuda_memory_manager_type::report()
{
  fprintf (stderr, "\n\nGPU memory map contents:\n");
  fprintf (stderr, "Object name                                                    Num objects            Total bytes\n");
  std::map<std::string,gpu_mem_object>::iterator iter;
  size_t total_bytes = 0, total_num=0;
  for (iter=gpu_mem_map.begin(); iter != gpu_mem_map.end(); iter++)
  {
    fprintf (stderr, "%60s %8ld            %13ld\n",
             iter->first.c_str(),
             iter->second.num_objects,
             iter->second.total_bytes);
    total_bytes += iter->second.total_bytes;
    total_num += iter->second.num_objects;
  }
  fprintf (stderr, "%60s %8ld            %13ld\n", "Total",
           total_num, total_bytes);
}
}

