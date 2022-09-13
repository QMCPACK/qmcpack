////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by:
// Alfredo Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef SPARSE_CFO_MATRIX_HPP
#define SPARSE_CFO_MATRIX_HPP

#include <array>
#include <cassert>
#include <iostream>
//#include<cstddef>  // ptrdiff_t
#include <vector>
#include <tuple>
#include <numeric>
#include <memory>
#include <type_traits> // enable_if
#include <algorithm>
#include <utility>

#include "Configuration.h"
#include "AFQMC/Utilities/tuple_iterator.hpp"

#include "AFQMC/Memory/custom_pointers.hpp"

#include "mpi3/shared_communicator.hpp"
#include "mpi.h"

namespace ma
{
namespace sparse
{
using tp_ul_ul  = std::tuple<std::size_t, std::size_t>;
using size_type = std::size_t;
using qmcplusplus::afqmc::to_address;
//using difference_type     = std::ptrdiff_t;
//using index               = std::ptrdiff_t;

template<class Allocator>
struct null_is_root
{
  null_is_root(Allocator) {}
  bool root() { return true; }
  int size() { return 1; }
  int rank() { return 0; }
  void barrier(){};
};

struct is_root
{
  mpi3::shared_communicator& comm_;
  template<class Allocator>
  is_root(Allocator& a) : comm_(*a.commP_)
  {}
  bool root() { return comm_.root(); }
  int size() { return comm_.size(); }
  int rank() { return comm_.rank(); }
  void barrier() { comm_.barrier(); }
};

template<class ValTypePtr, class IndxTypePtr = int*, class IntTypePtr = size_type*>
class csr_matrix_ref
{
public:
  using value_type  = decltype(*std::declval<ValTypePtr>());
  using element     = decltype(*std::declval<ValTypePtr>());
  using element_ptr = ValTypePtr;
  using index_type  = decltype(*std::declval<IndxTypePtr>());
  using int_type    = decltype(*std::declval<IntTypePtr>());

protected:
  using this_t = csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>;
  size_type size1_;
  size_type size2_;
  size_type local_origin1_;
  size_type local_origin2_;
  size_type global_origin1_;
  size_type global_origin2_;
  size_type capacity_;
  ValTypePtr data_;
  IndxTypePtr jdata_;
  IntTypePtr pointers_begin_;
  IntTypePtr pointers_end_;
  // set object to null state
  void reset()
  {
    size1_ = size2_ = 0;
    local_origin1_ = local_origin2_ = 0;
    global_origin1_ = global_origin2_ = 0;
    capacity_                         = 0;
    data_                             = ValTypePtr(nullptr);
    jdata_                            = IndxTypePtr(nullptr);
    pointers_begin_                   = IntTypePtr(nullptr);
    pointers_end_                     = IntTypePtr(nullptr);
  }

public:
  using row_type                  = std::tuple<size_type, ValTypePtr, IndxTypePtr>;
  static const bool sparse        = true;
  static const int dimensionality = -2;
  csr_matrix_ref(std::tuple<size_type, size_type> const& arr,
                 std::tuple<size_type, size_type> const& local,
                 std::tuple<size_type, size_type> const& global,
                 size_type cap_,
                 ValTypePtr __data_,
                 IndxTypePtr __jdata_,
                 IntTypePtr __pointers_begin_,
                 IntTypePtr __pointers_end_)
      : size1_(std::get<0>(arr)),
        size2_(std::get<1>(arr)),
        local_origin1_(std::get<0>(local)),
        local_origin2_(std::get<1>(local)),
        global_origin1_(std::get<0>(global)),
        global_origin2_(std::get<1>(global)),
        capacity_(cap_),
        data_(__data_),
        jdata_(__jdata_),
        pointers_begin_(__pointers_begin_),
        pointers_end_(__pointers_end_)
  {}
  csr_matrix_ref(this_t const& other) = delete;
  csr_matrix_ref& operator=(this_t const& other) = delete;
  // pointer movement is handled by derived classes
  csr_matrix_ref(this_t&& other) = default; //:csr_matrix_ref() { }
  csr_matrix_ref& operator=(this_t&& other) = default;
  ~csr_matrix_ref() {}
  auto pointers_begin(size_type i = 0) const { return pointers_begin_ + i; }
  auto pointers_end(size_type i = 0) const { return pointers_end_ + i; }
  auto local_origin() const { return std::array<size_type, 2>{{local_origin1_, local_origin2_}}; }
  auto global_origin() const { return std::array<size_type, 2>{{global_origin1_, global_origin2_}}; }
  auto size() const { return size1_; }
  template<typename integer_type = size_type>
  // not callable from the GPU!!!
  auto capacity(integer_type i) const
  {
    if (not pointers_begin_)
      return size_type(0);
    return static_cast<size_type>(pointers_begin_[i + 1] - pointers_begin_[i]);
  }
  auto capacity() const
  {
    if (not pointers_begin_)
      return size_type(0);
    return capacity_;
  }
  auto num_non_zero_elements() const
  {
    size_type ret = 0;
    for (size_type i = 0; i != size(); ++i)
      ret += static_cast<size_type>(pointers_end_[i] - pointers_begin_[i]);
    return ret;
  }
  template<typename integer_type = size_type>
  auto num_non_zero_elements(integer_type i) const
  {
    assert(i >= 0 && i < size1_);
    return static_cast<size_type>(pointers_end_[i] - pointers_begin_[i]);
  }
  auto shape() const { return std::array<size_type, 2>{{size(), size2_}}; }
  auto sizes() const {return shape();}
  template<typename integer_type = size_type>
  auto size(integer_type d) const
  {
    assert(d == integer_type{0} || d == integer_type{1});
    return (d == integer_type{0}) ? size() : size2_;
  }
  auto non_zero_values_data(size_type i = 0) const
  {
    if (i == 0)
      return data_;
    else
      return data_ + (pointers_begin_[i] - pointers_begin_[0]);
  }
  auto non_zero_indices2_data(size_type i = 0) const
  {
    if (i == 0)
      return jdata_;
    else
      return jdata_ + (pointers_begin_[i] - pointers_begin_[0]);
  }
  auto sparse_row(int i) const
  {
    assert(i >= 0 && i < size1_);
    return std::make_tuple(size_type(pointers_end_[i] - pointers_begin_[i]), data_ + pointers_begin_[i],
                           jdata_ + pointers_begin_[i]);
  }
  auto release_capacity()
  {
    auto t    = capacity_;
    capacity_ = size_type(0);
    return t;
  }
  auto release_non_zero_values_data()
  {
    auto t = data_;
    data_  = ValTypePtr(nullptr);
    return t;
  }
  auto release_non_zero_indices2_data()
  {
    auto t = jdata_;
    jdata_ = IndxTypePtr(nullptr);
    return t;
  }
  auto release_pointer_begin()
  {
    auto t          = pointers_begin_;
    pointers_begin_ = IntTypePtr(nullptr);
    return t;
  }
  auto release_pointer_end()
  {
    auto t        = pointers_end_;
    pointers_end_ = IntTypePtr(nullptr);
    return t;
  }
  // ??? not sure how to do this better
  auto release_shape()
  {
    auto t = shape();
    size1_ = size2_ = size_type(0);
    return t;
  }
  auto release_local_origin()
  {
    auto t         = local_origin();
    local_origin1_ = local_origin2_ = size_type(0);
    return t;
  }
  auto release_global_origin()
  {
    auto t          = global_origin();
    global_origin1_ = global_origin2_ = size_type(0);
    return t;
  }
  friend decltype(auto) size(this_t const& s) { return s.size(); }
  friend decltype(auto) shape(this_t const& s) { return s.shape(); }
  //        friend auto index_bases(this_t const& s){return s.index_bases();}
  friend auto capacity(this_t const& s) { return s.capacity(); }
  friend auto num_non_zero_elements(this_t const& s) { return s.num_non_zero_elements(); }
  friend auto non_zero_values_data(this_t const& s) { return s.non_zero_values_data(); }
  friend auto non_zero_indices2_data(this_t const& s) { return s.non_zero_indices2_data(); }
  friend auto pointers_begin(this_t const& s) { return s.pointers_begin(); }
  friend auto pointers_end(this_t const& s) { return s.pointers_end(); }
};

// View of a sub matrix. Owns pointer_begin and pointer_end
// Don't know how to make this class derive from the current csr_matrix_ref,
// so I'm just making a new class
template<class ValTypePtr, class IndxTypePtr = int*, class IntTypePtr = int*>
class csr_matrix_view : public csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>
{
public:
  using value_type  = typename std::decay<decltype(*std::declval<ValTypePtr>())>::type;
  using element     = typename std::decay<decltype(*std::declval<ValTypePtr>())>::type;
  using element_ptr = ValTypePtr;
  using index_type  = typename std::decay<decltype(*std::declval<IndxTypePtr>())>::type;
  using int_type    = typename std::decay<decltype(*std::declval<IntTypePtr>())>::type;

protected:
  std::vector<int_type> ptrb;
  std::vector<int_type> ptre;
  using base   = csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>;
  using this_t = csr_matrix_view<ValTypePtr, IndxTypePtr, IntTypePtr>;
  csr_matrix_view() : base(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, nullptr, nullptr, nullptr, nullptr) {}

public:
  using row_type                  = std::tuple<size_type, ValTypePtr, IndxTypePtr>;
  static const bool sparse        = true;
  static const bool view          = true;
  static const int dimensionality = -2;
  csr_matrix_view(std::tuple<size_type, size_type> const& arr,
                  std::tuple<size_type, size_type> const& local,
                  std::tuple<size_type, size_type> const& global,
                  ValTypePtr __data_,
                  IndxTypePtr __jdata_,
                  std::vector<int_type>&& va,
                  std::vector<int_type>&& vb)
      : base(arr, local, global, 0, __data_, __jdata_, IntTypePtr(nullptr), IntTypePtr(nullptr)),
        ptrb(std::move(va)),
        ptre(std::move(vb))
  {
    assert(ptrb.size() == base::size1_ + 1);
    assert(ptre.size() == base::size1_);
    base::pointers_begin_ = ptrb.data();
    base::pointers_end_   = ptre.data();
    base::capacity_       = base::num_non_zero_elements();
  }
  csr_matrix_view(this_t const& other) = delete;
  csr_matrix_view& operator=(this_t const& other) = delete;
  csr_matrix_view(this_t&& other) : csr_matrix_view() { *this = std::move(other); }
  csr_matrix_view& operator=(this_t&& other)
  {
    if (this != std::addressof(other))
    {
      base::size1_          = std::exchange(other.size1_, size_type(0));
      base::size2_          = std::exchange(other.size2_, size_type(0));
      base::local_origin1_  = std::exchange(other.local_origin1_, size_type(0));
      base::local_origin2_  = std::exchange(other.local_origin2_, size_type(0));
      base::global_origin1_ = std::exchange(other.global_origin1_, size_type(0));
      base::global_origin2_ = std::exchange(other.global_origin2_, size_type(0));
      base::capacity_       = std::exchange(other.capacity_, size_type(0));
      base::data_           = std::exchange(other.data_, ValTypePtr(nullptr));
      base::jdata_          = std::exchange(other.jdata_, IndxTypePtr(nullptr));
      ptrb                  = std::move(other.ptrb);
      ptre                  = std::move(other.ptre);
      base::pointers_begin_ = ptrb.data();
      base::pointers_end_   = ptre.data();
      other.pointers_begin_ = IntTypePtr(nullptr);
      other.pointers_end_   = IntTypePtr(nullptr);
    }
    return *this;
  }
  ~csr_matrix_view() {}
};

template<class ValType,
         class IndxType       = int,
         class IntType        = size_type,
         class ValType_alloc  = std::allocator<ValType>,
         class IsRoot         = null_is_root<ValType_alloc>,
         class IndxType_alloc = typename ValType_alloc::template rebind<IndxType>::other,
         class IntType_alloc  = typename ValType_alloc::template rebind<IntType>::other>
class ucsr_matrix : public csr_matrix_ref<typename ValType_alloc::pointer,
                                          typename IndxType_alloc::pointer,
                                          typename IntType_alloc::pointer>
{
public:
  using value_type  = ValType;
  using element     = ValType;
  using element_ptr = typename ValType_alloc::pointer;
  using index_type  = IndxType;
  using int_type    = IntType;
  using alloc_type  = ValType_alloc;

protected:
  using this_t = ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot, IndxType_alloc, IntType_alloc>;
  //	using IndxType_alloc = typename ValType_alloc::template rebind<IndxType>::other;
  //	using IntType_alloc = typename ValType_alloc::template rebind<IntType>::other;
  using ValTypePtr  = typename ValType_alloc::pointer;
  using IndxTypePtr = typename IndxType_alloc::pointer;
  using IntTypePtr  = typename IntType_alloc::pointer;
  using Valloc_ts   = std::allocator_traits<ValType_alloc>;
  using Ialloc_ts   = std::allocator_traits<IndxType_alloc>;
  using Palloc_ts   = std::allocator_traits<IntType_alloc>;
  using base        = csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>;
  ValType_alloc Valloc_;
  IndxType_alloc Ialloc_;
  IntType_alloc Palloc_;
  void reset()
  {
    IsRoot r(Valloc_);
    // calling barrier for safety right now
    r.barrier();
    size_type tot_sz = base::capacity();
    if (r.root())
    {
      // turning off destroy until I figure out a solution for GPU
      /*
                        if(base::data_ && base::pointers_begin_ && base::pointers_end_) 
                                for(size_type i = 0; i != base::size1_; ++i)
                                        for(auto p = base::data_ + base::pointers_begin_[i]; 
                                                 p != base::data_ + base::pointers_end_[i]; ++p)
                                                Valloc_.destroy(to_address(p));
                        if(base::jdata_ && base::pointers_begin_ && base::pointers_end_) 
                                for(size_type i = 0; i != base::size1_; ++i)
                                        for(auto p = base::jdata_ + base::pointers_begin_[i]; 
                                                 p != base::jdata_ + base::pointers_end_[i]; ++p)
                                                Ialloc_.destroy(to_address(p));
                        if(base::pointers_begin_ && base::pointers_end_) {
                                for(size_type i = 0; i != base::size1_; ++i){
                                        Palloc_.destroy(std::addressof(base::pointers_begin_[i]));
                                        Palloc_.destroy(std::addressof(base::pointers_end_[i]));
                                }
                                Palloc_.destroy(std::addressof(base::pointers_begin_[base::size1_]));
                        }
*/
    }
    r.barrier();
    if (base::data_)
      Valloc_.deallocate(base::data_, tot_sz);
    if (base::jdata_)
      Ialloc_.deallocate(base::jdata_, tot_sz);
    if (base::pointers_begin_)
      Palloc_.deallocate(base::pointers_begin_, base::size1_ + 1);
    if (base::pointers_end_)
      Palloc_.deallocate(base::pointers_end_, base::size1_);
    base::reset();
    r.barrier();
  }

public:
  static const bool sparse        = true;
  static const int dimensionality = -2;
  static const bool sorted        = false;
  template<typename integer_type = size_type>
  ucsr_matrix(std::tuple<size_type, size_type> const& arr,
              std::tuple<size_type, size_type> const& global,
              integer_type nnzpr_unique = 0,
              ValType_alloc alloc       = ValType_alloc{})
      : csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>(arr,
                                                            tp_ul_ul{0, 0},
                                                            global,
                                                            std::get<0>(arr) * nnzpr_unique,
                                                            ValTypePtr(nullptr),
                                                            IndxTypePtr(nullptr),
                                                            IntTypePtr(nullptr),
                                                            IntTypePtr(nullptr)),
        Valloc_(alloc),
        Ialloc_(alloc),
        Palloc_(alloc)
  {
    base::data_           = Valloc_.allocate(std::get<0>(arr) * nnzpr_unique);
    base::jdata_          = Ialloc_.allocate(std::get<0>(arr) * nnzpr_unique);
    base::pointers_begin_ = Palloc_.allocate(std::get<0>(arr) + 1);
    base::pointers_end_   = Palloc_.allocate(std::get<0>(arr));

    IsRoot r(Valloc_);
    if (nnzpr_unique > 0)
    {
      if (r.root())
      {
        auto pb(to_address(base::pointers_begin_));
        auto pe(to_address(base::pointers_end_));
        for (size_type i = 0; i != base::size1_; ++i)
        {
          *(pb + i) = i * nnzpr_unique;
          *(pe + i) = i * nnzpr_unique;
        }
        *(pb + base::size1_) = base::size1_ * nnzpr_unique;
      }
    }
    else
    {
      using std::fill_n;
      fill_n(base::pointers_begin_, std::get<0>(arr) + 1, int_type(0));
      fill_n(base::pointers_end_, std::get<0>(arr), int_type(0));
    }
    r.barrier();
  }
  template<typename integer_type = size_type>
  ucsr_matrix(std::tuple<size_type, size_type> const& arr,
              std::tuple<size_type, size_type> const& global,
              std::vector<integer_type> const& nnzpr = std::vector<integer_type>(0),
              ValType_alloc const& alloc             = {})
      : csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>(arr,
                                                            tp_ul_ul{0, 0},
                                                            global,
                                                            0,
                                                            ValTypePtr(nullptr),
                                                            IndxTypePtr(nullptr),
                                                            IntTypePtr(nullptr),
                                                            IntTypePtr(nullptr)),
        Valloc_(alloc),
        Ialloc_(alloc),
        Palloc_(alloc)
  {
    size_type sz = size_type(std::accumulate(nnzpr.begin(), nnzpr.end(), integer_type(0)));
    if (std::get<0>(arr) == 0)
      sz = 0; // no rows, no capacity
    base::capacity_       = sz;
    base::data_           = Valloc_.allocate(sz);
    base::jdata_          = Ialloc_.allocate(sz);
    base::pointers_begin_ = Palloc_.allocate(std::get<0>(arr) + 1);
    base::pointers_end_   = Palloc_.allocate(std::get<0>(arr));

    assert(nnzpr.size() >= base::size1_);
    IsRoot r(Valloc_);
    if (sz > 0)
    {
      if (r.root())
      {
        IntType cnter(0);
        auto pb(to_address(base::pointers_begin_));
        auto pe(to_address(base::pointers_end_));
        for (size_type i = 0; i != base::size1_; ++i)
        {
          *(pb + i) = cnter;
          *(pe + i) = cnter;
          cnter += static_cast<IntType>(nnzpr[i]);
        }
        *(pb + base::size1_) = cnter;
      }
    }
    else
    {
      using std::fill_n;
      fill_n(base::pointers_begin_, std::get<0>(arr) + 1, int_type(0));
      fill_n(base::pointers_end_, std::get<0>(arr), int_type(0));
    }
    r.barrier();
  }
  ~ucsr_matrix() { reset(); }
  ucsr_matrix(this_t const& other)
      : csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>(std::tuple<size_type, size_type>{other.size1_,
                                                                                             other.size2_},
                                                            std::tuple<size_type, size_type>{other.local_origin1_,
                                                                                             other.local_origin2_},
                                                            std::tuple<size_type, size_type>{other.global_origin1_,
                                                                                             other.global_origin2_},
                                                            0,
                                                            ValTypePtr(nullptr),
                                                            IndxTypePtr(nullptr),
                                                            IntTypePtr(nullptr),
                                                            IntTypePtr(nullptr)),
        Valloc_(other.getAlloc()),
        Ialloc_(Valloc_),
        Palloc_(Valloc_)
  {
    base::capacity_       = other.capacity_;
    base::data_           = Valloc_.allocate(base::capacity_);
    base::jdata_          = Ialloc_.allocate(base::capacity_);
    base::pointers_begin_ = Palloc_.allocate(base::size1_ + 1);
    base::pointers_end_   = Palloc_.allocate(base::size1_);
    //                IsRoot r(Valloc_);
    //                if(r.root()){
    using std::copy_n;
    copy_n(other.data_, base::capacity_, base::data_);
    copy_n(other.jdata_, base::capacity_, base::jdata_);
    copy_n(other.pointers_begin_, base::size1_ + 1, base::pointers_begin_);
    copy_n(other.pointers_end_, base::size1_, base::pointers_end_);
    //                }
    //                r.barrier();
  }
  ucsr_matrix& operator=(this_t const& other)
  {
    if (this != std::addressof(other))
    {
      base::reset();
      base::size1_          = other.size1_;
      base::size2_          = other.size2_;
      base::local_origin1_  = other.local_origin1_;
      base::local_origin2_  = other.local_origin2_;
      base::global_origin1_ = other.global_origin1_;
      base::global_origin2_ = other.global_origin2_;
      base::capacity_       = other.capacity_;
      base::data_           = Valloc_.allocate(base::capacity_);
      base::jdata_          = Ialloc_.allocate(base::capacity_);
      base::pointers_begin_ = Palloc_.allocate(base::size1_ + 1);
      base::pointers_end_   = Palloc_.allocate(base::size1_);
      //                    IsRoot r(Valloc_);
      //                    if(r.root()){
      using std::copy_n;
      copy_n(other.data_, base::capacity_, base::data_);
      copy_n(other.jdata_, base::capacity_, base::jdata_);
      copy_n(other.pointers_begin_, base::size1_ + 1, base::pointers_begin_);
      copy_n(other.pointers_end_, base::size1_, base::pointers_end_);
      //                    }
      //                    r.barrier();
    }
    return *this;
  }
  ucsr_matrix(this_t&& other) : ucsr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, other.Valloc_)
  {
    *this = std::move(other);
  }
  // Instead of moving allocators, require they are the same right now
  ucsr_matrix& operator=(this_t&& other)
  {
    if (this != std::addressof(other))
    {
      if (Valloc_ != other.Valloc_ || Ialloc_ != other.Ialloc_ || Palloc_ != other.Palloc_)
        APP_ABORT(" Error: Can only move assign between ucsr_matrices with equivalent allocators. \n");
      reset();
      base::size1_          = std::exchange(other.size1_, size_type(0));
      base::size2_          = std::exchange(other.size2_, size_type(0));
      base::local_origin1_  = std::exchange(other.local_origin1_, size_type(0));
      base::local_origin2_  = std::exchange(other.local_origin2_, size_type(0));
      base::global_origin1_ = std::exchange(other.global_origin1_, size_type(0));
      base::global_origin2_ = std::exchange(other.global_origin2_, size_type(0));
      base::capacity_       = std::exchange(other.capacity_, size_type(0));
      base::data_           = std::exchange(other.data_, ValTypePtr(nullptr));
      base::jdata_          = std::exchange(other.jdata_, IndxTypePtr(nullptr));
      base::pointers_begin_ = std::exchange(other.pointers_begin_, IntTypePtr(nullptr));
      base::pointers_end_   = std::exchange(other.pointers_end_, IntTypePtr(nullptr));
    }
    return *this;
  }
  ValType_alloc getAlloc() const { return Valloc_; }
  //auto getAlloc() const { return Valloc_; }
  template<typename integer_type = size_type>
  void reserve(integer_type nnzpr_unique)
  {
    if (base::size1_ == 0)
      return;
    IntType minN = IntType(base::pointers_begin_[1] - base::pointers_begin_[0]);
    for (size_type i = 0; i != base::size1_; ++i)
      minN = std::min(minN, base::pointers_begin_[i + 1] - base::pointers_begin_[i]);
    if (static_cast<IntType>(nnzpr_unique) <= minN)
      return;
    this_t other(tp_ul_ul{base::size1_, base::size2_}, tp_ul_ul{base::global_origin1_, base::global_origin2_},
                 nnzpr_unique, Valloc_);
    if (base::capacity_ > 0)
    {
      IsRoot r(Valloc_);
      //                    if(r.root()){
      for (size_type i = 0; i < base::size1_; i++)
      {
        size_type disp = static_cast<size_type>(base::pointers_end_[i] - base::pointers_begin_[i]);
        using std::copy_n;
        copy_n(base::data_ + base::pointers_begin_[i], disp, other.data_ + other.pointers_begin_[i]);
        copy_n(base::jdata_ + base::pointers_begin_[i], disp, other.jdata_ + other.pointers_begin_[i]);
        if (r.root())
          other.pointers_end_[i] = other.pointers_begin_[i] + disp;
      }
      //                    }
      r.barrier();
    }
    *this = std::move(other);
  }
  template<typename integer_type>
  void reserve(std::vector<integer_type> const& nnzpr)
  {
    IsRoot r_(Valloc_);
    if (base::size1_ == 0)
      return;
    bool resz = false;
    assert(nnzpr.size() >= base::size1_);
    for (size_type i = 0; i < base::size1_; i++)
      if (static_cast<IntType>(nnzpr[i]) > base::pointers_begin_[i + 1] - base::pointers_begin_[i])
      {
        resz = true;
        break;
      }
    if (not resz)
      return;
    this_t other(tp_ul_ul{base::size1_, base::size2_}, tp_ul_ul{base::global_origin1_, base::global_origin2_}, nnzpr,
                 Valloc_);
    if (base::capacity_ > 0)
    {
      IsRoot r(Valloc_);
      //                    if(r.root()){
      for (size_type i = 0; i < base::size1_; i++)
      {
        size_type disp = static_cast<size_type>(base::pointers_end_[i] - base::pointers_begin_[i]);
        using std::copy_n;
        std::copy_n(base::data_ + base::pointers_begin_[i], disp, other.data_ + other.pointers_begin_[i]);
        std::copy_n(base::jdata_ + base::pointers_begin_[i], disp, other.jdata_ + other.pointers_begin_[i]);
        if (r.root())
          other.pointers_end_[i] = other.pointers_begin_[i] + disp;
      }
      //                    }
      r.barrier();
    }
    *this = std::move(other);
  }
  template<class Pair = std::array<IndxType, 2>, class... Args>
  void emplace(Pair&& indices, Args&&... args)
  {
    using std::get;
    assert(get<0>(indices) >= 0);
    assert(get<0>(indices) < base::size1_);
    if (base::pointers_end_[get<0>(indices)] < base::pointers_begin_[get<0>(indices) + 1])
    {
      ::new (static_cast<void*>(to_address(base::data_ + base::pointers_end_[get<0>(indices)])))
          value_type(std::forward<Args>(args)...);
      //Valloc_ts::construct(Valloc_,to_address(base::data_ + base::pointers_end_[get<0>(indices)]), std::forward<Args>(args)...);
      ::new (static_cast<void*>(to_address(base::jdata_ + base::pointers_end_[get<0>(indices)])))
          index_type(get<1>(indices));
      //Ialloc_ts::construct(Ialloc_,to_address(base::jdata_ + base::pointers_end_[get<0>(indices)]), get<1>(indices));
      ++base::pointers_end_[get<0>(indices)];
    }
    else
    {
      APP_ABORT(" Error: row size exceeded the maximum \n\n\n");
      throw std::out_of_range("row size exceeded the maximum");
    }
  }
  template<typename integer_type = IndxType, typename value_type = ValType>
  void emplace(std::tuple<integer_type, integer_type, value_type> const& val)
  {
    using std::get;
    emplace({get<0>(val), get<1>(val)}, static_cast<ValType>(get<2>(val)));
  }
  template<class Pair = std::array<IndxType, 2>, class... Args>
  void emplace_back(Pair&& indices, Args&&... args)
  {
    emplace(std::forward<Pair>(indices), std::forward<Args>(args)...);
  }
  template<typename integer_type = IndxType, typename value_type = ValType>
  void emplace_back(std::tuple<integer_type, integer_type, value_type> const& val)
  {
    using std::get;
    emplace({get<0>(val), get<1>(val)}, static_cast<ValType>(get<2>(val)));
  }

protected:
  struct row_reference
  {
    ucsr_matrix& self_;
    IndxType i_;
    struct element_reference
    {
      row_reference& self_;
      IndxType j_;
      template<class TT>
      element_reference&& operator=(TT&& tt) &&
      {
        self_.self_.emplace({{self_.i_, j_}}, std::forward<TT>(tt));
        return std::move(*this);
      }
    };
    using reference = element_reference;
    reference operator[](IndxType i) && { return reference{*this, i}; }
  };

public:
  using reference = row_reference;
  template<typename integer_type = size_type>
  reference operator[](integer_type i)
  {
    return reference{*this, static_cast<IndxType>(i)};
  }
  ValType_alloc& getValloc() { return Valloc_; }
  IndxType_alloc& getIalloc() { return Ialloc_; }
  IntType_alloc& getPalloc() { return Palloc_; }
};

template<class ValType,
         class IndxType       = int,
         class IntType        = size_type,
         class ValType_alloc  = std::allocator<ValType>,
         class IsRoot         = null_is_root<ValType_alloc>,
         class IndxType_alloc = typename ValType_alloc::template rebind<IndxType>::other,
         class IntType_alloc  = typename ValType_alloc::template rebind<IntType>::other>
class csr_matrix : public ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot, IndxType_alloc, IntType_alloc>
{
  using this_t = csr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot, IndxType_alloc, IntType_alloc>;

public:
  using base        = ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot, IndxType_alloc, IntType_alloc>;
  using value_type  = ValType;
  using element     = ValType;
  using element_ptr = typename ValType_alloc::pointer;
  using index_type  = IndxType;
  using int_type    = IntType;
  using alloc_type  = ValType_alloc;
  using ValTypePtr  = typename ValType_alloc::pointer;
  using IndxTypePtr = typename IndxType_alloc::pointer;
  using IntTypePtr  = typename IntType_alloc::pointer;
  static const bool sparse        = true;
  static const int dimensionality = -2;
  static const bool sorted        = true;
  template<typename integer_type = size_type>
  csr_matrix(std::tuple<size_type, size_type> const& arr,
             std::tuple<size_type, size_type> const& global,
             integer_type nnzpr_unique = 0,
             ValType_alloc alloc       = ValType_alloc{})
      : base(arr, global, nnzpr_unique, alloc)
  {}
  template<typename integer_type = size_type>
  csr_matrix(std::tuple<size_type, size_type> const& arr,
             std::tuple<size_type, size_type> const& global,
             std::vector<integer_type>& nnzpr = std::vector<integer_type>(0),
             ValType_alloc alloc              = ValType_alloc{})
      : base(arr, global, nnzpr, alloc)
  {}
  // for now, copy is inefficient, but keeping it simple
  csr_matrix(this_t const& csr)
      : base(std::tuple<size_type, size_type>{csr.size1_, csr.size2_},
             std::tuple<size_type, size_type>{csr.global_origin1_, csr.global_origin2_},
             0,
             csr.getAlloc())
  {
    *this = csr;
  }
  // right now, this routine is limited to transfers from host-host, or host-device.
  // Will fail if transferring device-to-host, since I need to use to_address on source
  template<class ValType_,
           class IndxType_,
           class IntType_,
           class ValType_alloc_,
           class IsRoot_,
           class IndxType_alloc_,
           class IntType_alloc_,
           typename = std::enable_if_t<not std::is_same<ValType_alloc_, ValType_alloc>::value>>
  csr_matrix(
      csr_matrix<ValType_, IndxType_, IntType_, ValType_alloc_, IsRoot_, IndxType_alloc_, IntType_alloc_> const& csr,
      ValType_alloc const& alloc = {})
      : //IntType_alloc_> const& csr, ValType_alloc const& alloc = {}):
        base(std::tuple<size_type, size_type>{0, 0}, std::tuple<size_type, size_type>{0, 0}, 0, alloc)
  {
    base::reset();
    auto shape_           = csr.shape();
    base::size1_          = shape_[0];
    base::size2_          = shape_[1];
    auto local_           = csr.local_origin();
    base::local_origin1_  = local_[0];
    base::local_origin2_  = local_[1];
    auto global_          = csr.global_origin();
    base::global_origin1_ = global_[0];
    base::global_origin2_ = global_[1];
    base::capacity_       = csr.capacity();
    base::data_           = base::Valloc_.allocate(base::capacity_);
    base::jdata_          = base::Ialloc_.allocate(base::capacity_);
    base::pointers_begin_ = base::Palloc_.allocate(base::size1_ + 1);
    base::pointers_end_   = base::Palloc_.allocate(base::size1_);
    //                IsRoot r(base::Valloc_);
    //                if(r.root()){
    using std::copy_n;
    copy_n(to_address(csr.non_zero_values_data()), base::capacity_, base::data_);
    copy_n(to_address(csr.non_zero_indices2_data()), base::capacity_, base::jdata_);
    copy_n(to_address(csr.pointers_begin()), base::size1_ + 1, base::pointers_begin_);
    copy_n(to_address(csr.pointers_end()), base::size1_, base::pointers_end_);
    //                }
    //                r.barrier();
  }

  csr_matrix& operator=(this_t const& csr)
  {
    if (this != std::addressof(csr))
    {
      base::reset();
      base::size1_          = csr.size1_;
      base::size2_          = csr.size2_;
      base::local_origin1_  = csr.local_origin1_;
      base::local_origin2_  = csr.local_origin2_;
      base::global_origin1_ = csr.global_origin1_;
      base::global_origin2_ = csr.global_origin2_;
      base::capacity_       = csr.capacity_;
      base::data_           = base::Valloc_.allocate(base::capacity_);
      base::jdata_          = base::Ialloc_.allocate(base::capacity_);
      base::pointers_begin_ = base::Palloc_.allocate(base::size1_ + 1);
      base::pointers_end_   = base::Palloc_.allocate(base::size1_);
      //                    IsRoot r(base::Valloc_);
      //                    if(r.root()){
      using std::copy_n;
      copy_n(csr.data_, base::capacity_, base::data_);
      copy_n(csr.jdata_, base::capacity_, base::jdata_);
      copy_n(csr.pointers_begin_, base::size1_ + 1, base::pointers_begin_);
      copy_n(csr.pointers_end_, base::size1_, base::pointers_end_);
      //                    }
      //                    r.barrier();
    }
    return *this;
  }
  template<class ValType_alloc_,
           class IsRoot_,
           class IndxType_alloc_,
           class IntType_alloc_,
           typename = std::enable_if_t<not std::is_same<ValType_alloc_, ValType_alloc>::value>>
  csr_matrix& operator=(
      csr_matrix<ValType, IndxType, IntType, ValType_alloc_, IsRoot_, IndxType_alloc_, IntType_alloc_> const& csr)
  {
    base::reset();
    auto shape_           = csr.shape();
    base::size1_          = shape_[0];
    base::size2_          = shape_[1];
    auto local_           = csr.local_origin();
    base::local_origin1_  = local_[0];
    base::local_origin2_  = local_[1];
    auto global_          = csr.global_origin();
    base::global_origin1_ = global_[0];
    base::global_origin2_ = global_[1];
    base::capacity_       = csr.capacity();
    base::data_           = base::Valloc_.allocate(base::capacity_);
    base::jdata_          = base::Ialloc_.allocate(base::capacity_);
    base::pointers_begin_ = base::Palloc_.allocate(base::size1_ + 1);
    base::pointers_end_   = base::Palloc_.allocate(base::size1_);
    //                IsRoot r(base::Valloc_);
    //                if(r.root()){
    using std::copy_n;
    copy_n(to_address(csr.non_zero_values_data()), base::capacity_, base::data_);
    copy_n(to_address(csr.non_zero_indices2_data()), base::capacity_, base::jdata_);
    copy_n(to_address(csr.pointers_begin()), base::size1_ + 1, base::pointers_begin_);
    copy_n(to_address(csr.pointers_end()), base::size1_, base::pointers_end_);
    //                }
    //                r.barrier();
    return *this;
  }

  csr_matrix& operator=(ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot> const& other)
  {
    base::reset();
    auto shape_           = other.shape();
    base::size1_          = shape_[0];
    base::size2_          = shape_[1];
    auto local_           = other.local_origin();
    base::local_origin1_  = local_[0];
    base::local_origin2_  = local_[1];
    auto global_          = other.global_origin();
    base::global_origin1_ = global_[0];
    base::global_origin2_ = global_[1];
    base::capacity_       = other.capacity();
    base::data_           = base::Valloc_.allocate(base::capacity_);
    base::jdata_          = base::Ialloc_.allocate(base::capacity_);
    base::pointers_begin_ = base::Palloc_.allocate(base::size1_ + 1);
    base::pointers_end_   = base::Palloc_.allocate(base::size1_);
    if (base::size1_ == 0 || base::capacity_ == 0)
      return *this;
    using qmcplusplus::make_paired_iterator;
    IsRoot r(base::Valloc_);
    //                if(r.root()){
    using std::copy_n;
    copy_n(other.non_zero_values_data(), base::capacity_, base::data_);
    copy_n(other.non_zero_indices2_data(), base::capacity_, base::jdata_);
    copy_n(other.pointers_begin(), base::size1_ + 1, base::pointers_begin_);
    copy_n(other.pointers_end(), base::size1_, base::pointers_end_);
    //                }
    //                r.barrier();
    for (size_type p = 0; p < base::size1_; p++)
    {
      if (p % static_cast<size_type>(r.size()) == static_cast<size_type>(r.rank()))
      {
        auto i1 = base::pointers_begin_[p];
        auto i2 = base::pointers_end_[p];
        std::sort(make_paired_iterator(to_address(base::jdata_ + i1), to_address(base::data_ + i1)),
                  make_paired_iterator(to_address(base::jdata_ + i2), to_address(base::data_ + i2)),
                  [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });
      }
    }
    r.barrier();
    return *this;
  }
  csr_matrix(this_t&& other) : csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, other.Valloc_)
  {
    *this = std::move(other);
  }
  csr_matrix(ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot> const& ucsr)
      : csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, ucsr.getAlloc())
  {
    *this = ucsr;
  }
  csr_matrix(ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot>&& ucsr)
      : csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, ucsr.getAlloc())
  {
    *this = std::move(ucsr);
  }
  csr_matrix& operator=(csr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot>&& other)
  {
    if (this != std::addressof(other))
    {
      if (base::Valloc_ != other.Valloc_ || base::Ialloc_ != other.Ialloc_ || base::Palloc_ != other.Palloc_)
        APP_ABORT(" Error: Can only move assign between csr_matrices with equivalent allocators. \n");
      base::reset();
      base::size1_          = std::exchange(other.size1_, size_type(0));
      base::size2_          = std::exchange(other.size2_, size_type(0));
      base::local_origin1_  = std::exchange(other.local_origin1_, size_type(0));
      base::local_origin2_  = std::exchange(other.local_origin2_, size_type(0));
      base::global_origin1_ = std::exchange(other.global_origin1_, size_type(0));
      base::global_origin2_ = std::exchange(other.global_origin2_, size_type(0));
      base::capacity_       = std::exchange(other.capacity_, size_type(0));
      base::data_           = std::exchange(other.data_, ValTypePtr(nullptr));
      base::jdata_          = std::exchange(other.jdata_, IndxTypePtr(nullptr));
      base::pointers_begin_ = std::exchange(other.pointers_begin_, IntTypePtr(nullptr));
      base::pointers_end_   = std::exchange(other.pointers_end_, IntTypePtr(nullptr));
    }
    return *this;
  }
  csr_matrix& operator=(ucsr_matrix<ValType, IndxType, IntType, ValType_alloc, IsRoot>&& other)
  {
    if (base::Valloc_ != other.getValloc() || base::Ialloc_ != other.getIalloc() || base::Palloc_ != other.getPalloc())
      APP_ABORT(" Error: Can only move assign between ucsr/csr_matrices with equivalent allocators. \n");

    base::reset();
    auto shape_           = other.release_shape();
    base::size1_          = shape_[0];
    base::size2_          = shape_[1];
    auto local_           = other.release_local_origin();
    base::local_origin1_  = local_[0];
    base::local_origin2_  = local_[1];
    auto global_          = other.release_global_origin();
    base::global_origin1_ = global_[0];
    base::global_origin2_ = global_[1];
    base::capacity_       = other.release_capacity();
    base::data_           = other.release_non_zero_values_data();
    base::jdata_          = other.release_non_zero_indices2_data();
    base::pointers_begin_ = other.release_pointer_begin();
    base::pointers_end_   = other.release_pointer_end();
    if (base::size1_ == 0 || base::capacity_ == 0)
      return *this;
    using qmcplusplus::make_paired_iterator;
    IsRoot r(base::Valloc_);
    for (size_type p = 0; p < base::size1_; p++)
    {
      if (p % static_cast<size_type>(r.size()) == static_cast<size_type>(r.rank()))
      {
        auto i1 = base::pointers_begin_[p];
        auto i2 = base::pointers_end_[p];
        std::sort(make_paired_iterator(to_address(base::jdata_ + i1), to_address(base::data_ + i1)),
                  make_paired_iterator(to_address(base::jdata_ + i2), to_address(base::data_ + i2)),
                  [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });
      }
    }
    r.barrier();
    return *this;
  }
  template<class Pair = std::array<IndxType, 2>, class... Args>
  void emplace(Pair&& indices, Args&&... args)
  {
    using std::get;
    assert(get<0>(indices) >= 0);
    assert(get<0>(indices) < base::size1_);
    if (base::pointers_end_[get<0>(indices)] < base::pointers_begin_[get<0>(indices) + 1])
    {
      auto loc = std::lower_bound(to_address(base::jdata_ + base::pointers_begin_[get<0>(indices)]),
                                  to_address(base::jdata_ + base::pointers_end_[get<0>(indices)]), get<1>(indices));
      size_type disp =
          std::distance(to_address(base::jdata_ + base::pointers_begin_[get<0>(indices)]), to_address(loc));
      size_type disp_ = std::distance(to_address(loc), to_address(base::jdata_ + base::pointers_end_[get<0>(indices)]));
      if (disp_ > 0 && *loc == get<1>(indices))
      {
        // value exists, construct in place
        //base::Valloc_ts::construct(base::Valloc_,to_address(base::data_ + base::pointers_begin_[get<0>(indices)] + disp), std::forward<Args>(args)...);
        ::new (static_cast<void*>(to_address(base::data_ + base::pointers_begin_[get<0>(indices)] + disp)))
            value_type(std::forward<Args>(args)...);
      }
      else
      {
        // new value, shift back and add in correct place
        if (disp_ > 0)
        {
          std::move_backward(to_address(base::data_ + base::pointers_begin_[get<0>(indices)] + disp),
                             to_address(base::data_ + base::pointers_end_[get<0>(indices)]),
                             to_address(base::data_ + base::pointers_end_[get<0>(indices)] + 1));
          std::move_backward(to_address(base::jdata_ + base::pointers_begin_[get<0>(indices)] + disp),
                             to_address(base::jdata_ + base::pointers_end_[get<0>(indices)]),
                             to_address(base::jdata_ + base::pointers_end_[get<0>(indices)] + 1));
        }
        ++base::pointers_end_[get<0>(indices)];
        //base::Valloc_ts::construct(base::Valloc_,to_address(base::data_+base::pointers_begin_[get<0>(indices)] + disp), std::forward<Args>(args)...);
        ::new (static_cast<void*>(to_address(base::data_ + base::pointers_begin_[get<0>(indices)] + disp)))
            value_type(std::forward<Args>(args)...);
        //base::Ialloc_ts::construct(base::Ialloc_, to_address(base::jdata_+base::pointers_begin_[get<0>(indices)] + disp), get<1>(indices));
        ::new (static_cast<void*>(to_address(base::jdata_ + base::pointers_begin_[get<0>(indices)] + disp)))
            index_type(get<1>(indices));
      }
    }
    else
      throw std::out_of_range("row size exceeded the maximum");
  }
  // new column index must be larger than all previous column indexes in the row
  // otherwise throws
  template<class Pair = std::array<IndxType, 2>, class... Args>
  void emplace_back(Pair&& indices, Args&&... args)
  {
    using std::get;
    assert(get<0>(indices) >= 0);
    assert(get<0>(indices) < base::size1_);
    if (base::pointers_end_[get<0>(indices)] < base::pointers_begin_[get<0>(indices) + 1])
    {
      // if row is empty or new column index is larger than last column in row
      if (base::pointers_begin_[get<0>(indices)] == base::pointers_end_[get<0>(indices)] or
          get<1>(indices) > base::jdata_[base::pointers_end_[get<0>(indices)] - 1])
      {
        //base::Valloc_ts::construct(base::Valloc_,to_address(base::data_+base::pointers_end_[get<0>(indices)]), std::forward<Args>(args)...);
        ::new (static_cast<void*>(to_address(base::data_ + base::pointers_end_[get<0>(indices)])))
            value_type(std::forward<Args>(args)...);
        //base::Ialloc_ts::construct(base::Ialloc_, to_address(base::jdata_+base::pointers_end_[get<0>(indices)]), get<1>(indices));
        ::new (static_cast<void*>(to_address(base::jdata_ + base::pointers_end_[get<0>(indices)])))
            index_type(get<1>(indices));
        ++base::pointers_end_[get<0>(indices)];
      }
      else // otherwise throw
        throw std::runtime_error("inconsistent column index in emplace_back");
    }
    else
      throw std::out_of_range("row size exceeded the maximum");
  }
  template<typename integer_type = IndxType, typename value_type = ValType>
  void emplace(std::tuple<integer_type, integer_type, value_type> const& val)
  {
    using std::get;
    emplace({get<0>(val), get<1>(val)}, static_cast<ValType>(get<2>(val)));
  }
  template<typename integer_type = IndxType, typename value_type = ValType>
  void emplace_back(std::tuple<integer_type, integer_type, value_type> const& val)
  {
    using std::get;
    emplace_back({get<0>(val), get<1>(val)}, static_cast<ValType>(get<2>(val)));
  }
  void remove_empty_spaces()
  {
    IsRoot r(base::Valloc_);
    if (r.root())
    {
      for (size_type i = 0; i < base::size1_ - 1; i++)
      {
        if (base::pointers_end_[i] == base::pointers_begin_[i + 1])
          continue;
        auto ni = static_cast<size_type>(base::pointers_end_[i + 1] - base::pointers_begin_[i + 1]);
        std::move(to_address(base::data_ + base::pointers_begin_[i + 1]),
                  to_address(base::data_ + base::pointers_end_[i + 1]),
                  to_address(base::data_ + base::pointers_end_[i]));
        std::move(to_address(base::jdata_ + base::pointers_begin_[i + 1]),
                  to_address(base::jdata_ + base::pointers_end_[i + 1]),
                  to_address(base::jdata_ + base::pointers_end_[i]));
        base::pointers_begin_[i + 1] = base::pointers_end_[i];
        base::pointers_end_[i + 1]   = base::pointers_begin_[i + 1] + ni;
      }
      base::pointers_begin_[base::size1_] = base::pointers_end_[base::size1_ - 1];
    }
    r.barrier();
  }

protected:
  struct row_reference
  {
    static const int dimensionality = 1;
    static const bool sparse        = true;
    this_t& self_;
    IndxType i_;
    struct element_reference
    {
      row_reference& self_;
      IndxType j_;
      template<class TT>
      element_reference&& operator=(TT&& tt) &&
      {
        self_.self_.emplace({{self_.i_, j_}}, std::forward<TT>(tt));
        return std::move(*this);
      }
    };
    using reference = element_reference;
    template<typename integer_type = size_type>
    reference operator[](integer_type i) &&
    {
      return reference{*this, static_cast<IndxType>(i)};
    }

    auto non_zero_values_data() const { return self_.non_zero_values_data(i_); }
    auto non_zero_indices2_data() const { return self_.non_zero_indices2_data(i_); }
    auto num_non_zero_elements() const { return size_type(self_.pointers_end_[i_] - self_.pointers_begin_[i_]); }
    auto capacity() const { return self_.capacity(i_); }
    auto shape() const { return std::array<size_type, 1>{{self_.size2_}}; }
    template<typename integer_type = size_type>
    auto size(integer_type d) const
    {
      return size_type{self_.size2_};
    }
  };

  struct const_row_reference
  {
    static const int dimensionality = 1;
    static const bool sparse        = true;
    const this_t& self_;
    const IndxType i_;

    auto non_zero_values_data() const { return self_.non_zero_values_data(i_); }
    auto non_zero_indices2_data() const { return self_.non_zero_indices2_data(i_); }
    auto num_non_zero_elements() const { return size_type(self_.pointers_end_[i_] - self_.pointers_begin_[i_]); }
    auto capacity() const { return self_.capacity(i_); }
    auto shape() const { return std::array<size_type, 1>{{self_.size2_}}; }
	auto sizes() const {return shape();}
    template<typename integer_type = size_type>
    auto size(integer_type d) const
    {
      return size_type{self_.size2_};
    }
  };

public:
  using reference       = row_reference;
  using const_reference = const_row_reference;
  template<class integer_type>
  reference operator[](integer_type i)
  {
    return reference{*this, static_cast<IndxType>(i)};
  }

  template<class integer_type = size_type>
  const_reference operator[](integer_type i) const
  {
    return const_reference{*this, static_cast<IndxType>(i)};
  }

  using sub_matrix = csr_matrix_ref<ValTypePtr, IndxTypePtr, IntTypePtr>;
  template<class integer_type>
  sub_matrix operator[](std::array<integer_type, 2> const& arr)
  {
    assert(arr[0] >= 0 && size_type(arr[1]) <= base::size1_);
    assert(arr[0] < arr[1]);
    size_type disp = static_cast<size_type>(base::pointers_begin_[arr[0]] - base::pointers_begin_[0]);
    size_type cap  = 0;
    for (size_type i = arr[0]; i < arr[1]; i++)
      cap += base::capacity(i);
    // careful!!! elements or row 'r' are always indexed with pointer_begin[r]-pointer_begin[0]
    return sub_matrix(tp_ul_ul{(arr[1] - arr[0]), base::size2_}, tp_ul_ul{arr[0], 0},
                      tp_ul_ul{base::global_origin1_ + arr[0], base::global_origin2_}, cap, base::data_ + disp,
                      base::jdata_ + disp, base::pointers_begin_ + static_cast<size_type>(arr[0]),
                      base::pointers_end_ + static_cast<size_type>(arr[0]));
  }

  template<typename IntT = int>
  using matrix_view = csr_matrix_view<ValTypePtr, IndxTypePtr, IntT*>;
  template<class IntT = int, class integer_type = size_type>
  matrix_view<IntT> operator[](std::array<integer_type, 4> const& arr)
  {
    // limited right now
    assert(base::capacity_ > 0);
    assert(base::size1_ > 0);
    assert(base::size2_ > 0);
    assert(arr[0] >= 0 && size_type(arr[1]) <= base::size1_);
    assert(arr[2] >= 0 && size_type(arr[3]) <= base::size2_);
    assert(arr[0] < arr[1]);
    assert(arr[2] < arr[3]);

    std::vector<IntT> ptrb, ptre;
    ptrb.reserve((arr[1] - arr[0]) + 1);
    ptre.reserve((arr[1] - arr[0]));

    if (arr[2] == 0 && arr[3] == base::size2_)
    {
      // no need to search
      std::size_t p0 = base::pointers_begin_[arr[0]];
      for (std::size_t r = arr[0]; r < arr[1]; r++)
      {
        std::size_t dbr = base::pointers_begin_[r] - p0;
        std::size_t der = base::pointers_end_[r] - p0;
        if (dbr >= static_cast<std::size_t>(std::numeric_limits<IntT>::max()))
          throw std::out_of_range("row size exceeded the maximum");
        if (der >= static_cast<std::size_t>(std::numeric_limits<IntT>::max()))
          throw std::out_of_range("row size exceeded the maximum");
        ptrb.emplace_back(static_cast<IntT>(dbr));
        ptre.emplace_back(static_cast<IntT>(der));
      }
      ptrb.emplace_back(ptre.back());

      // columns always begin in 0, since column values can't be shifted
      // but only columns in range [arr[2],arr[3]) are accessible/visible
      return matrix_view<IntT>(tp_ul_ul{(arr[1] - arr[0]), arr[3]},
                               tp_ul_ul{arr[0], arr[2]}, // local_origin2_ is the only way to know the "true" origin
                               tp_ul_ul{base::global_origin1_ + arr[0], base::global_origin2_},
                               base::non_zero_values_data(arr[0]), base::non_zero_indices2_data(arr[0]),
                               std::move(ptrb), std::move(ptre));
    }
    else
    {
      // reference position (wrt base::pointer_begin_[0])
      auto p0 = base::pointers_begin_[0];
      auto ref_col =
          std::lower_bound(to_address(base::jdata_ + base::pointers_begin_[arr[0]] - p0),
                           to_address(base::jdata_ + base::pointers_end_[arr[0]] - p0), static_cast<IndxType>(arr[2]));
      for (std::size_t r = arr[0]; r < arr[1]; r++)
      {
        auto br =
            std::lower_bound(to_address(base::jdata_ + base::pointers_begin_[r] - p0),
                             to_address(base::jdata_ + base::pointers_end_[r] - p0), static_cast<IndxType>(arr[2]));
        auto er =
            std::lower_bound(to_address(base::jdata_ + base::pointers_begin_[r] - p0),
                             to_address(base::jdata_ + base::pointers_end_[r] - p0), static_cast<IndxType>(arr[3]));
        std::size_t dbr = std::distance(ref_col, br);
        std::size_t der = std::distance(ref_col, er);
        if (dbr >= static_cast<std::size_t>(std::numeric_limits<IntT>::max()))
          throw std::out_of_range("row size exceeded the maximum");
        if (der >= static_cast<std::size_t>(std::numeric_limits<IntT>::max()))
          throw std::out_of_range("row size exceeded the maximum");
        ptrb.emplace_back(static_cast<IntT>(dbr));
        ptre.emplace_back(static_cast<IntT>(der));
      }
      ptrb.emplace_back(ptre.back());
      auto d0 = std::distance(to_address(base::jdata_), ref_col);

      // columns always begin in 0, since column values can't be shifted
      // but only columns in range [arr[2],arr[3]) are accessible/visible
      return matrix_view<IntT>(tp_ul_ul{(arr[1] - arr[0]), arr[3]},
                               tp_ul_ul{arr[0], arr[2]}, // local_origin2_ is the only way to know the "true" origin
                               tp_ul_ul{base::global_origin1_ + arr[0], base::global_origin2_}, base::data_ + d0,
                               base::jdata_ + d0, std::move(ptrb), std::move(ptre));
    }
  }
};

} // namespace sparse
} // namespace ma

#endif
