//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MPI_CONTAINER_PROXY_H
#define QMCPLUSPLUS_MPI_CONTAINER_PROXY_H

#include <stdexcept>

#include "OhmmsPETE/Tensor.h"
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
  using real_type  = T;
  using value_type = T;
  static inline T* get_address(T* a) { return a; }
};

template<typename T>
struct scalar_traits<std::complex<T>>
{
  enum
  {
    DIM = 2
  };
  using real_type  = T;
  using value_type = std::complex<T>;
  static inline T* get_address(std::complex<T>* a) { return reinterpret_cast<T*>(a); }
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

template<typename T, unsigned D>
struct container_proxy<TinyVector<T, D>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM * D
  };
  using pointer = typename scalar_traits<T>::real_type*;
  TinyVector<T, D>& ref;
  inline container_proxy(TinyVector<T, D>& a) : ref(a) {}
  inline size_t size() const { return DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};

template<typename T, unsigned D>
struct container_proxy<Tensor<T, D>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM * D * D
  };
  using pointer = typename scalar_traits<T>::real_type*;
  Tensor<T, D>& ref;
  inline container_proxy(Tensor<T, D>& a) : ref(a) {}
  inline size_t size() const { return DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};

template<typename T>
struct container_proxy<std::vector<T>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  std::vector<T>& ref;
  inline container_proxy(std::vector<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_proxy<T>::DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(&ref[0]); }

  inline void resize(size_t n) { ref.resize(n); }

  template<typename I>
  inline void resize(I* n, int d)
  {
    size_t nt = n[0];
    for (int i = 1; i < d; ++i)
      nt *= n[i];
    ref.resize(nt);
  }
};

template<>
struct container_proxy<std::vector<bool>>
{
  enum
  {
    DIM = 1
  };
  using pointer = int*;
  std::vector<bool>& ref;
  std::vector<int> my_copy;
  inline container_proxy(std::vector<bool>& a) : ref(a)
  {
    my_copy.resize(a.size());
    copy(a.begin(), a.end(), my_copy.begin());
  }
  ~container_proxy() { copy(my_copy.begin(), my_copy.end(), ref.begin()); }
  inline size_t size() const { return my_copy.size(); }
  inline pointer data() { return &my_copy[0]; }
};

template<typename T, unsigned D>
struct container_proxy<std::vector<TinyVector<T, D>>>
{
  enum
  {
    DIM = D * scalar_traits<T>::DIM
  };
  using pointer   = typename container_proxy<T>::pointer;
  using data_type = std::vector<TinyVector<T, D>>;
  data_type& ref;
  inline container_proxy(data_type& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref[0].data()); }
};


template<typename T>
struct container_proxy<PooledData<T>>
{
  enum
  {
    DIM = 1
  };
  using pointer = typename container_proxy<T>::pointer;
  PooledData<T>& ref;
  inline container_proxy(PooledData<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_proxy<T>::DIM; }
  inline pointer data() { return ref.data(); }
  template<typename I>
  inline void resize(I* n)
  {
    ref.resize(static_cast<size_t>(n[0]));
  }
};

template<typename T>
struct container_proxy<Vector<T>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  Vector<T>& ref;
  inline container_proxy(Vector<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_proxy<T>::DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
  template<typename I>
  inline void resize(I* n)
  {
    ref.resize(static_cast<size_t>(n[0]));
  }
};

template<typename T>
struct container_proxy<Matrix<T>>
{
  enum
  {
    DIM = scalar_traits<T>::DIM
  };
  using pointer = typename container_proxy<T>::pointer;
  Matrix<T>& ref;
  inline container_proxy(Matrix<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size(); }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
  template<typename I>
  inline void resize(I* n, int d)
  {
    if (d != 2)
      throw std::runtime_error("OhmmsMatrix can only be resized with int[2].");
    ref.resize(n[0], n[1]);
  }
};

template<typename T, unsigned D>
struct container_proxy<Vector<TinyVector<T, D>>>
{
  enum
  {
    DIM = D * scalar_traits<T>::DIM
  };
  using pointer   = typename container_proxy<T>::pointer;
  using data_type = Vector<TinyVector<T, D>>;
  data_type& ref;
  inline container_proxy(data_type& a) : ref(a) {}
  inline size_t size() const { return ref.size() * DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref[0].data()); }
};

template<typename T, unsigned D>
struct container_proxy<Array<T, D>>
{
  using pointer = typename container_proxy<T>::pointer;
  Array<T, D>& ref;
  inline container_proxy(Array<T, D>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_proxy<T>::DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
};
/*
template<typename T, class Alloc>
struct container_proxy<boost::multi::array<T,2,Alloc> >
{
  enum {DIM=scalar_traits<T>::DIM};
  using pointer = typename container_proxy<T>::pointer;
  boost::multi::array<T,2,Alloc>& ref;
  inline container_proxy(boost::multi::array<T,2,Alloc>& a):ref(a) {}
  inline size_t size() const
  {
    return ref.num_elements()*DIM;
  }
  inline pointer data()
  {
    //using detail::to_address;
    //return scalar_traits<T>::get_address(to_address(ref.origin()));
    return scalar_traits<T>::get_address(std::addressof(*ref.origin()));
  }
  inline void resize(size_t n)
  {
    APP_ABORT(" Error: Can not resize container_proxy<boost::multi::array<T,D,Alloc> >. \n");
  }
  template<typename I>
  inline void resize(I* n, int d)
  {
    if(d < 2)
      APP_ABORT(" Error: Inconsistent dimension in container_proxy<boost::multi::array<T,D,Alloc> >::resize(I*,int). \n");
    ref.reextent({n[0],n[1]});
  }
};

template<typename T>
struct container_proxy<boost::multi::array_ref<T,2> >
{
  enum {DIM=scalar_traits<T>::DIM};
  using pointer = typename container_proxy<T>::pointer;
  boost::multi::array_ref<T,2>& ref;
  inline container_proxy(boost::multi::array_ref<T,2>& a):ref(a) {}
  inline size_t size() const
  {
    return ref.num_elements()*DIM;
  }
  inline pointer data()
  {
    //using detail::to_address;
    //return scalar_traits<T>::get_address(to_address(ref.origin()));
    return scalar_traits<T>::get_address(std::addressof(*ref.origin()));
  }
  inline void resize(size_t n)
  {
    APP_ABORT(" Error: Can not resize container_proxy<boost::multi::array_ref<T,D> >. \n");
  }

  template<typename I>
  inline void resize(I* n, int d)
  {
    APP_ABORT(" Error: Can not resize container_proxy<boost::multi::array_ref<T,D> >. \n");
  }

};
*/
} // namespace qmcplusplus
#endif
