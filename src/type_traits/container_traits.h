//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Lab
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Lab
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MPI_CONTAINER_TRAITS_H
#define QMCPLUSPLUS_MPI_CONTAINER_TRAITS_H

namespace qmcplusplus
{

template<typename CT>
struct container_traits;

template<typename T>
struct container_traits<std::vector<T>>
{
  enum
  {
    DIM = 1
  };
/*
  typedef typename container_traits<T>::pointer pointer;
  std::vector<T>& ref;
  inline container_traits(std::vector<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_traits<T>::DIM; }
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
*/
};

template<>
struct container_traits<std::vector<bool>>
{
  enum
  {
    DIM = 1
  };
/*
  typedef int* pointer;
  std::vector<bool>& ref;
  std::vector<int> my_copy;
  inline container_traits(std::vector<bool>& a) : ref(a)
  {
    my_copy.resize(a.size());
    copy(a.begin(), a.end(), my_copy.begin());
  }
  ~container_traits() { copy(my_copy.begin(), my_copy.end(), ref.begin()); }
  inline size_t size() const { return my_copy.size(); }
  inline pointer data() { return &my_copy[0]; }
*/
};

template<typename T>
struct container_traits<PooledData<T>>
{
  enum
  {
    DIM = 1
  };
/*
  typedef typename container_traits<T>::pointer pointer;
  PooledData<T>& ref;
  inline container_traits(PooledData<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_traits<T>::DIM; }
  inline pointer data() { return ref.data(); }
  template<typename I>
  inline void resize(I* n)
  {
    ref.resize(static_cast<size_t>(n[0]));
  }
*/
};

template<typename T>
struct container_traits<Vector<T>>
{
  enum
  {
    DIM = 1
  };
/*
  typedef typename container_traits<T>::pointer pointer;
  Vector<T>& ref;
  inline container_traits(Vector<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_traits<T>::DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
  template<typename I>
  inline void resize(I* n)
  {
    ref.resize(static_cast<size_t>(n[0]));
  }
*/
};

template<typename T>
struct container_traits<Matrix<T>>
{
  enum
  {
    DIM = 2
  };
/*
  typedef typename container_traits<T>::pointer pointer;
  Matrix<T>& ref;
  inline container_traits(Matrix<T>& a) : ref(a) {}
  inline size_t size() const { return ref.size(); }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
  template<typename I>
  inline void resize(I* n, int d)
  {
    if ( d != 2 )
      throw std::runtime_error("OhmmsMatrix can only be resized with int[2].");
    ref.resize(n[0],n[1]);
  }
*/
};

template<typename T, unsigned D>
struct container_traits<Array<T, D>>
{
  enum
  {
    DIM = D
  };
/*
  typedef typename container_traits<T>::pointer pointer;
  Array<T, D>& ref;
  inline container_traits(Array<T, D>& a) : ref(a) {}
  inline size_t size() const { return ref.size() * container_traits<T>::DIM; }
  inline pointer data() { return scalar_traits<T>::get_address(ref.data()); }
*/
};
} // namespace qmcplusplus
#endif
