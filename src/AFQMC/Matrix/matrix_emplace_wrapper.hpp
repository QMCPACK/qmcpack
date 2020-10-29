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

#ifndef MATRIX_EMPLACE_WRAPPER_HPP
#define MATRIX_EMPLACE_WRAPPER_HPP

#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
#include <mutex>

#include "mpi3/shm/mutex.hpp"

namespace csr
{
// keeps a local (to the core) buffer during matrix assembly to minimize impact
// of locks during emplace of shm matrix
// can writ a specialization for non-shm matrices later
template<class Matrix>
struct matrix_emplace_wrapper
{
  using value_type = typename Matrix::value_type;
  using index_type = typename Matrix::index_type;
  using shm_mutex  = boost::mpi3::shm::mutex;

public:
  matrix_emplace_wrapper() : M(nullptr), m(nullptr) {}

  matrix_emplace_wrapper(Matrix& mat_,
                         boost::mpi3::shared_communicator& node,
                         std::size_t sz = MAXIMUM_EMPLACE_BUFFER_SIZE)
      : M(std::addressof(mat_)), m(nullptr)
  {
    m = std::make_unique<shm_mutex>(node);
    buff.reserve(std::max(sz, std::size_t(0)));
  }

  // not sure this makes sense, but needed for TTI
  matrix_emplace_wrapper(matrix_emplace_wrapper const& other) : M(other.M), m(nullptr)
  {
    m = std::move(std::make_unique<shm_mutex>(other.m->scomm_));
    buff.reserve(other.buff.size());
  }
  matrix_emplace_wrapper operator=(matrix_emplace_wrapper const& other) = delete;
  matrix_emplace_wrapper(matrix_emplace_wrapper&& other)
  {
    APP_ABORT(" Error: matrix_emplace_wrapper move constructor has been disabled. \n");
  }
  matrix_emplace_wrapper operator=(matrix_emplace_wrapper&& other) = delete;

  std::array<std::size_t, 2> shape()
  {
    if (not M)
      return {0, 0};
    return M->shape();
  }

  template<typename Size>
  std::size_t size(Size d)
  {
    if (not M)
      return 0;
    return M->size(d);
  }

  template<typename integer_type>
  void reserve(integer_type sz)
  {
    if (M)
      M->reserve(sz);
  }

  template<typename integer_type>
  void reserve(std::vector<integer_type> const& sz)
  {
    if (M)
      M->reserve(sz);
  }

  void emplace(index_type i, index_type j, value_type v)
  {
    if (not M)
      return;
    if (buff.capacity() > 0)
    {
      buff.emplace_back(std::make_tuple(i, j, v));
      if (buff.size() == buff.capacity())
        push_buffer();
    }
    else
    {
      M->emplace({i, j}, v);
    }
  }

  template<class Pair = std::array<index_type, 2>>
  void emplace(Pair&& indices, value_type v)
  {
    if (not M)
      return;
    using std::get;
    if (buff.capacity() > 0)
    {
      buff.emplace_back(std::make_tuple(get<0>(indices), get<1>(indices), v));
      if (buff.size() == buff.capacity())
        push_buffer();
    }
    else
    {
      M->emplace({get<0>(indices), get<1>(indices)}, v);
    }
  }

  void emplace(std::tuple<index_type, index_type, value_type> const& v)
  {
    if (not M)
      return;
    using std::get;
    if (buff.capacity() > 0)
    {
      buff.emplace_back(v);
      if (buff.size() == buff.capacity())
        push_buffer();
    }
    else
    {
      M->emplace({get<0>(v), get<1>(v)}, get<2>(v));
    }
  }

  void emplace_back(index_type i, index_type j, value_type v)
  {
    if (not M)
      return;
    if (buff.capacity() > 0)
    {
      buff.emplace_back(std::make_tuple(i, j, v));
      if (buff.size() == buff.capacity())
        push_buffer();
    }
    else
    {
      M->emplace_back({i, j}, v);
    }
  }

  template<class Pair = std::array<index_type, 2>>
  void emplace_back(Pair&& indices, value_type v)
  {
    if (not M)
      return;
    using std::get;
    if (buff.capacity() > 0)
    {
      buff.emplace_back(std::make_tuple(get<0>(indices), get<1>(indices), v));
      if (buff.size() == buff.capacity())
        push_buffer();
    }
    else
    {
      M->emplace_back({get<0>(indices), get<1>(indices)}, v);
    }
  }

  void emplace_back(std::tuple<index_type, index_type, value_type> const& v)
  {
    if (not M)
      return;
    using std::get;
    if (buff.capacity() > 0)
    {
      buff.emplace_back(v);
      if (buff.size() == buff.capacity())
        push_buffer();
    }
    else
    {
      M->emplace_back({get<0>(v), get<1>(v)}, get<2>(v));
    }
  }

  void push_buffer()
  {
    using std::get;
    if (not M)
      return;
    if (buff.size() == 0)
      return;
    {
      std::lock_guard<shm_mutex> guard(*m);
      for (auto& t : buff)
        M->emplace({get<0>(t), get<1>(t)}, get<2>(t));
      buff.clear();
    }
  }

private:
  Matrix* M;
  std::vector<std::tuple<index_type, index_type, value_type>> buff;
  std::unique_ptr<shm_mutex> m;
};

} // namespace csr

#endif
