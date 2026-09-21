//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MPI_CONTAINER_PROXY_H
#define QMCPLUSPLUS_MPI_CONTAINER_PROXY_H

#include <stdexcept>

#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Pools/PooledData.h"

namespace qmcplusplus
{
template<class T>
struct scalar_traits
{
  enum
  {
    DIM = 1
  };
  using real_type = T;
  static inline T* get_address(T* a) { return a; }
};

template<typename T>
struct scalar_traits<std::complex<T>>
{
  enum
  {
    DIM = 2
  };
  using real_type = T;
  static inline T* get_address(std::complex<T>* a) { return reinterpret_cast<T*>(a); }
};

template<typename T, unsigned D>
struct scalar_traits<TinyVector<T, D>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM * D
  };
  using real_type = typename scalar_traits<T>::real_type;
  static inline real_type* get_address(TinyVector<T, D>* a) { return scalar_traits<T>::get_address(a->data()); }
};

template<typename T, unsigned D>
struct scalar_traits<Tensor<T, D>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM * D * D
  };
  using real_type = typename scalar_traits<T>::real_type;
  static inline real_type* get_address(Tensor<T, D>* a) { return scalar_traits<T>::get_address(a->data()); }
};


template<typename T>
struct container_proxy
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename scalar_traits<T>::real_type*;
  T& ref;
  inline container_proxy(T& a) : ref(a) {}
  inline size_t size() const { return DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(&ref); }
};

template<typename T, typename Alloc>
struct container_proxy<std::vector<T, Alloc>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  std::vector<T, Alloc>& ref;
  inline container_proxy(std::vector<T, Alloc>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};

template<typename Alloc>
struct container_proxy<std::vector<bool, Alloc>>
{
  enum
  {
    DIM = 1
  };
  using pointer = int*;
  std::vector<bool, Alloc>& ref;
  std::vector<int> my_copy;
  inline container_proxy(std::vector<bool, Alloc>& a) : ref(a)
  {
    my_copy.resize(a.size());
    copy(a.begin(), a.end(), my_copy.begin());
  }
  ~container_proxy() { copy(my_copy.begin(), my_copy.end(), ref.begin()); }
  inline size_t size() const { return my_copy.size(); }
  inline pointer data() { return &my_copy[0]; }
};

template<typename T>
struct container_proxy<PooledData<T>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  PooledData<T>& ref;
  inline container_proxy(PooledData<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};

template<typename T, typename Alloc>
struct container_proxy<Vector<T, Alloc>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  Vector<T, Alloc>& ref;
  inline container_proxy(Vector<T, Alloc>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};

template<typename T, typename Alloc>
struct container_proxy<Matrix<T, Alloc>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  Matrix<T, Alloc>& ref;
  inline container_proxy(Matrix<T, Alloc>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};

template<typename T, unsigned D, typename ALLOC>
struct container_proxy<Array<T, D, ALLOC>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  Array<T, D, ALLOC>& ref;
  inline container_proxy(Array<T, D, ALLOC>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};
} // namespace qmcplusplus
#endif
