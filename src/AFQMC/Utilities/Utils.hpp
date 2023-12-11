////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_UTILITIES_UTILS_HPP
#define AFQMC_UTILITIES_UTILS_HPP

#include <numeric>
#include <stack>
#include <iostream>
#include <fstream>
#include <complex>
#include <list>

#include "Host/sysutil.h"
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/config.0.h"

#ifdef ENABLE_CUDA
#include "AFQMC/Memory/device_pointers.hpp"
#include "AFQMC/Memory/CUDA/cuda_arch.h"
#include "AFQMC/Numerics/device_kernels.hpp"
#include "cuda_runtime.h"
#elif ENABLE_HIP
#include "AFQMC/Memory/device_pointers.hpp"
#include "AFQMC/Memory/HIP/hip_arch.h"
#include "AFQMC/Numerics/device_kernels.hpp"
#include "hip/hip_runtime.h"
#endif

namespace qmcplusplus
{
namespace afqmc
{
inline bool file_exists(const std::string& name)
{
  std::ifstream f(name.c_str());
  return f.good();
}

template<class T, class Factory>
T get_parameter(Factory& F, std::string const& obj_name, std::string const& pname, T const& def)
{
  xmlNodePtr cur = F.getXML(obj_name);
  if (cur == nullptr)
    return def;
  else
  {
    T res(def);
    ParameterSet m_param;
    m_param.add(res, pname.c_str());
    m_param.put(cur);
    return res;
  }
}

template<class Iter, class Compare>
void parallel_inplace_merge(int np, int rk, Iter* beg, Iter* mid, Iter* end, MPI_Comm comm, Compare comp)
{
  if (np == 1)
  {
    std::inplace_merge(beg, mid, end, comp);
    return;
  }

  MPI_Barrier(comm);

  Iter *p1, *p2;
  if (std::distance(beg, mid) >= std::distance(mid, end))
  {
    p1      = beg + std::distance(beg, mid) / 2;
    auto it = std::lower_bound(mid, end, *p1, comp);
    p2      = &(*it);
  }
  else
  {
    p2      = mid + std::distance(mid, end) / 2;
    auto it = std::lower_bound(beg, mid, *p2, comp);
    p1      = &(*it);
  }

  MPI_Barrier(comm);
  if (rk == 0)
    std::rotate(p1, mid, p2);
  MPI_Barrier(comm);

  mid = p1 + std::distance(mid, p2);

  if (rk < np / 2)
    parallel_inplace_merge<Iter, Compare>(np / 2, rk, beg, p1, mid, comm, comp);
  else
    parallel_inplace_merge<Iter, Compare>(np / 2, rk - np / 2, mid, p2, end, comm, comp);
}

// Simple circular buffer object.
// Assumed to be of fixed size and should not be resized / reallocated.
// Access elements of the buffer through member function `current'.
template<class T>
class CircularBuffer
{
public:
  // default constructor.
  CircularBuffer() : nelem(0), head(0), buffer(0){};
  // construct 2D buffer. Useful for indexing buffer of vectors or objects.
  CircularBuffer(int nelements, int ncols)
  {
    buffer.resize(nelements, T(ncols));
    nelem = nelements;
    head  = 0;
  }
  // construct 1D buffer.
  CircularBuffer(int nelements)
  {
    buffer.resize(nelements);
    nelem = nelements;
    head  = 0;
  }
  // copy
  CircularBuffer& operator=(const CircularBuffer& cb)
  {
    if (this == &cb)
    {
      return *this;
    }
    else
    {
      nelem  = cb.nelem;
      head   = cb.head;
      buffer = cb.buffer;
    }
  }
  // destructor
  ~CircularBuffer(){};
  // get pointer to current entry (defined by head) in buffer.
  T* current() { return &buffer[head]; }
  // forward traversal.
  void increment() { head = (head + 1) % nelem; }
  // backward traversal.
  void decrement() { head = (head - 1 + nelem) % nelem; }
  // reset head to original location.
  void reset() { head = 0; }

private:
  // number of elements in buffer.
  int nelem;
  // current entry in buffer.
  int head;
  std::vector<T> buffer;
};

template<typename IType, typename integer>
void balance_partition_ordered_set(integer N, IType const* indx, std::vector<IType>& subsets)
{
  int64_t avg = 0;

  if (*(indx + N) == 0)
    exit(1);

  IType nsets = subsets.size() - 1;
  IType i0    = 0;
  IType iN    = N;
  while (*(indx + i0) == *(indx + i0 + 1))
    i0++;
  while (*(indx + iN - 1) == *(indx + iN))
    iN--;
  avg = static_cast<int64_t>(*(indx + iN)) - static_cast<int64_t>(*(indx + i0));
  avg /= nsets;

  // no c++14 :-(
  //    template<class Iter>
  //    auto partition = [=] (IType i0, IType iN, int n, Iter vals) {
  auto partition = [=](IType i0, IType iN, int n, typename std::list<IType>::iterator vals) {
    // finds optimal position for subsets[i]
    auto step = [=](IType i0, IType iN, IType& ik) {
      IType imin  = ik;
      ik          = i0 + 1;
      double v1   = double(std::abs(static_cast<int64_t>(*(indx + ik)) - static_cast<int64_t>(*(indx + i0)) - avg));
      double v2   = double(std::abs(static_cast<int64_t>(*(indx + iN)) - static_cast<int64_t>(*(indx + ik)) - avg));
      double vmin = v1 * v1 + v2 * v2;
      for (int k = i0 + 2, kend = iN; k < kend; k++)
      {
        v1       = double(std::abs(static_cast<int64_t>(*(indx + k)) - static_cast<int64_t>(*(indx + i0)) - avg));
        v2       = double(std::abs(static_cast<int64_t>(*(indx + iN)) - static_cast<int64_t>(*(indx + k)) - avg));
        double v = v1 * v1 + v2 * v2;
        if (v < vmin)
        {
          vmin = v;
          ik   = k;
        }
      }
      return ik != imin;
    };

    if (n == 2)
    {
      *vals = i0 + 1;
      step(i0, iN, *vals);
      return;
    }

    std::vector<IType> set(n + 1);
    set[0] = i0;
    set[n] = iN;
    for (int i = n - 1; i >= 1; i--)
      set[i] = iN + i - n;
    bool changed;
    do
    {
      changed = false;
      for (IType i = 1; i < n; i++)
        changed |= step(set[i - 1], set[i + 1], set[i]);
    } while (changed);

    std::copy_n(set.begin() + 1, n - 1, vals);

    return;
  };

  // dummy factorization
  std::stack<IType> factors;
  IType n0 = nsets;
  for (IType i = 2; i <= nsets; i++)
  {
    while (n0 % i == 0)
    {
      factors.push(i);
      n0 /= i;
    }
    if (n0 == 1)
      break;
  }
  assert(n0 == 1);

  std::list<IType> sets;
  sets.push_back(i0);
  sets.push_back(iN);

  while (factors.size() > 0)
  {
    auto ns = factors.top();
    factors.pop();

    // divide all current partitions into ns sub-partitions
    typename std::list<IType>::iterator it = sets.begin();
    it++;
    for (; it != sets.end(); it++)
    {
      typename std::list<IType>::iterator its = it;
      its--;
      auto i0 = *its;
      its     = sets.insert(it, std::size_t(ns - 1), i0 + 1);
      partition(i0, *it, ns, its);
    }
  }

  typename std::list<IType>::iterator it    = sets.begin();
  typename std::vector<IType>::iterator itv = subsets.begin();
  for (; itv < subsets.end(); itv++, it++)
    *itv = *it;

  return;
}

template<typename IType>
void balance_partition_ordered_set(std::vector<IType> const& indx, std::vector<IType>& subsets)
{
  if (indx.size() < 2 || subsets.size() < 2)
    return;
  balance_partition_ordered_set(indx.size() - 1, indx.data(), subsets);
}

template<typename IType>
void balance_partition_ordered_set_wcounts(std::vector<IType> const& counts, std::vector<IType>& subsets)
{
  if (counts.size() == 0 || subsets.size() < 2)
    return;
  subsets.resize(counts.size() + 1);
  std::vector<IType> indx(counts.size() + 1);
  IType cnt = 0;
  auto it   = indx.begin();
  *it++     = 0;
  for (auto& v : counts)
    *it++ = (cnt += v);
  balance_partition_ordered_set(counts.size(), indx.data(), subsets);
}

template<class Vec,
         class RandomNumberGenerator_,
         typename = typename std::enable_if_t<std::decay<Vec>::type::dimensionality == 1>>
void sampleGaussianFields(Vec&& V, RandomNumberGenerator_& rng)
{
  size_t n = V.size();
  for (int i = 0; i + 1 < n; i += 2)
  {
    RealType temp1 = 1 - 0.9999999999 * rng(), temp2 = rng();
    RealType mag = std::sqrt(-2.0 * std::log(temp1));
    V[i]         = mag * std::cos(6.283185306 * temp2);
    V[i + 1]     = mag * std::sin(6.283185306 * temp2);
  }
  if (n % 2 == 1)
  {
    RealType temp1 = 1 - 0.9999999999 * rng(), temp2 = rng();
    V[n - 1] = std::sqrt(-2.0 * std::log(temp1)) * std::cos(6.283185306 * temp2);
  }
}

template<class Mat,
         class RandomNumberGenerator_,
         typename = typename std::enable_if_t<(std::decay<Mat>::type::dimensionality > 1)>,
         typename = void>
void sampleGaussianFields(Mat&& M, RandomNumberGenerator_& rng)
{
  for (int i = 0, iend = M.size(0); i < iend; ++i)
    sampleGaussianFields(M[i], rng);
}

template<class T, class RandomNumberGenerator_>
void sampleGaussianFields_n(T* V, int n, RandomNumberGenerator_& rng)
{
  for (int i = 0; i + 1 < n; i += 2)
  {
    RealType temp1 = 1 - 0.9999999999 * rng(), temp2 = rng();
    RealType mag = std::sqrt(-2.0 * std::log(temp1));
    V[i]         = T(mag * std::cos(6.283185306 * temp2));
    V[i + 1]     = T(mag * std::sin(6.283185306 * temp2));
  }
  if (n % 2 == 1)
  {
    RealType temp1 = 1 - 0.9999999999 * rng(), temp2 = rng();
    V[n - 1] = T(std::sqrt(-2.0 * std::log(temp1)) * std::cos(6.283185306 * temp2));
  }
}

inline void memory_report()
{
  qmcplusplus::app_log() << "\n --> CPU Memory Available: " << (freemem() >> 20) << std::endl;
#ifdef ENABLE_CUDA
  size_t free_, tot_;
  cudaMemGetInfo(&free_, &tot_);
  qmcplusplus::app_log() << " --> GPU Memory Available,  Total in MB: " << (free_ >> 20) << " "
                         << (tot_ >> 20) << "\n"
                         << std::endl;
#elif ENABLE_HIP
  size_t free_, tot_;
  hipMemGetInfo(&free_, &tot_);
  qmcplusplus::app_log() << " --> GPU Memory Available,  Total in MB: " << (free_ >> 20) << " "
                         << (tot_ >> 20) << "\n"
                         << std::endl;
#endif
}

// TODO: FDM : why not use standard naming convention like arch::afqmc_rand_generator?
#if defined(ENABLE_CUDA)
template<class T, class Dummy>
void sampleGaussianFields_n(device::device_pointer<T> V, int n, Dummy& r)
{
  kernels::sampleGaussianRNG(to_address(V), n, arch::afqmc_curand_generator);
}
#elif defined(ENABLE_HIP)
template<class T, class Dummy>
void sampleGaussianFields_n(device::device_pointer<T> V, int n, Dummy& r)
{
  kernels::sampleGaussianRNG(to_address(V), n, arch::afqmc_rocrand_generator);
}
#endif

} // namespace afqmc

} // namespace qmcplusplus

#endif
