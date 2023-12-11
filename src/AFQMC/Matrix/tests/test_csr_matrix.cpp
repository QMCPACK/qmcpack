//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
//////////////////////////////////////////////////////////////////////////////////////

#undef NDEBUG

#include "catch.hpp"

#include <algorithm> // std::sort
#include <cassert>
#include <iostream>
#include <random>

#include "mpi3/shared_window.hpp"
#include "mpi3/shared_communicator.hpp"
#include "AFQMC/Memory/SharedMemory/shm_ptr_with_raw_ptr_dispatch.hpp"

#include "AFQMC/Matrix/csr_matrix.hpp"
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "AFQMC/Memory/custom_pointers.hpp"
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::get;
using tp_ul_ul = std::tuple<std::size_t, std::size_t>;

namespace mpi3 = boost::mpi3;

namespace qmcplusplus
{
template<typename Type, typename IndxType, typename IntType, class Alloc, class is_root>
void test_csr_matrix_shm_allocator(Alloc A, bool serial)
{
  auto world = boost::mpi3::environment::get_world_instance();
  mpi3::shared_communicator node(world.split_shared());

  using ucsr_matrix = ma::sparse::ucsr_matrix<Type, IndxType, IntType, Alloc, is_root>;
  using csr_matrix  = ma::sparse::csr_matrix<Type, IndxType, IntType, Alloc, is_root>;

  std::vector<Type> v_                  = {9, 10, 3, 1};
  auto itv                              = v_.begin();
  std::vector<IndxType> c_              = {2, 1, 1, 3};
  auto itc                              = c_.begin();
  std::vector<int> non_zero_per_row     = {2, 0, 1, 1};
  std::vector<int> max_non_zero_per_row = {5, 3, 4, 2};

  {
    ucsr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, 2, A);
    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();
  }

  {
    ucsr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    itv      = v_.begin();
    itc      = c_.begin();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();
  }

  {
    ucsr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, max_non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    itv      = v_.begin();
    itc      = c_.begin();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();
  }

  {
    ucsr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, 2, A);

    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    itv      = v_.begin();
    itc      = c_.begin();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();

    ucsr_matrix small0(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, 2, A);

    if (serial || node.rank() == 0)
      small0[3][3] = 1;
    if (serial || node.rank() == 0)
      small0[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small0[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small0[0][1] = 10;
    node.barrier();


    REQUIRE(small0.num_non_zero_elements() == 4);
    val = small0.non_zero_values_data();
    col = small0.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small0.size(0); i++)
    {
      for (auto it = small0.pointers_begin()[i]; it != small0.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    ucsr_matrix small2(std::move(small));
    REQUIRE(small2.num_non_zero_elements() == 4);
    val = small2.non_zero_values_data();
    col = small2.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small2.size(0); i++)
    {
      for (auto it = small2.pointers_begin()[i]; it != small2.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();

    ucsr_matrix small3(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, A);
    small3 = std::move(small2);
    REQUIRE(small3.num_non_zero_elements() == 4);
    val = small3.non_zero_values_data();
    col = small3.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small3.size(0); i++)
    {
      for (auto it = small3.pointers_begin()[i]; it != small3.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }


    // copy assignment
    small2 = small3;
    REQUIRE(small2.num_non_zero_elements() == 4);
    val = small2.non_zero_values_data();
    col = small2.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small2.size(0); i++)
    {
      for (auto it = small2.pointers_begin()[i]; it != small2.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    // copy constructor
    auto small9(small2);
    REQUIRE(small9.num_non_zero_elements() == 4);
    val = small9.non_zero_values_data();
    col = small9.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small9.size(0); i++)
    {
      for (auto it = small9.pointers_begin()[i]; it != small9.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();
    // ordered
    v_ = {10, 9, 3, 1};
    c_ = {1, 2, 1, 3};

    csr_matrix small10(small2);
    REQUIRE(small10.num_non_zero_elements() == 4);
    val = small10.non_zero_values_data();
    col = small10.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small10.size(0); i++)
    {
      for (auto it = small10.pointers_begin()[i]; it != small10.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    csr_matrix small4(std::move(small3));
    REQUIRE(small4.num_non_zero_elements() == 4);
    val = small4.non_zero_values_data();
    col = small4.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small4.size(0); i++)
    {
      for (auto it = small4.pointers_begin()[i]; it != small4.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    {
      using dev_csr_matrix = ma::sparse::csr_matrix<Type, IndxType, IntType, device::device_allocator<Type>>;
      dev_csr_matrix small11(small4);
      REQUIRE(small11.num_non_zero_elements() == 4);
      auto val_ = small11.non_zero_values_data();
      auto col_ = small11.non_zero_indices2_data();
      itv       = v_.begin();
      itc       = c_.begin();
      for (std::size_t i = 0; i < small11.size(0); i++)
      {
        for (IntType it = small11.pointers_begin()[i]; it != small11.pointers_end()[i]; it++)
        {
          REQUIRE(Type(val_[it]) == *(itv++));
          REQUIRE(col_[it] == *(itc++));
        }
      }

      // this is not a move! just making sure it makes a copy, otherwise
      // it will fail below
      dev_csr_matrix small12(std::move(small4));
      REQUIRE(small12.num_non_zero_elements() == 4);
      val_ = small12.non_zero_values_data();
      col_ = small12.non_zero_indices2_data();
      itv  = v_.begin();
      itc  = c_.begin();
      for (std::size_t i = 0; i < small12.size(0); i++)
      {
        for (IntType it = small12.pointers_begin()[i]; it != small12.pointers_end()[i]; it++)
        {
          REQUIRE(val_[it] == *(itv++));
          REQUIRE(col_[it] == *(itc++));
        }
      }
    }
#endif

    ucsr_matrix small5(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small5[3][3] = 1;
    if (serial || node.rank() == 0)
      small5[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small5[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small5[0][1] = 10;
    node.barrier();
    small4.reserve(max_non_zero_per_row);
    small4 = std::move(small5);
    REQUIRE(small4.num_non_zero_elements() == 4);
    val = small4.non_zero_values_data();
    col = small4.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small4.size(0); i++)
    {
      for (auto it = small4.pointers_begin()[i]; it != small4.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    ucsr_matrix small6(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, max_non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small6[3][3] = 1;
    if (serial || node.rank() == 0)
      small6[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small6[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small6[0][1] = 10;
    node.barrier();
    small4.reserve(100);
    small4 = std::move(small6);
    REQUIRE(small4.num_non_zero_elements() == 4);
    val = small4.non_zero_values_data();
    col = small4.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small4.size(0); i++)
    {
      for (auto it = small4.pointers_begin()[i]; it != small4.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    csr_matrix small7(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, max_non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small7[3][3] = 1;
    if (serial || node.rank() == 0)
      small7[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small7[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small7[0][1] = 10;
    node.barrier();
    small7.remove_empty_spaces();
    REQUIRE(small7.num_non_zero_elements() == 4);
    val = small7.non_zero_values_data();
    col = small7.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small7.size(0); i++)
    {
      if (i < small7.size(0) - 1)
        REQUIRE(small7.pointers_end()[i] - small7.pointers_begin()[i] == small7.capacity(i));
      for (auto it = small7.pointers_begin()[i]; it != small7.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    // copy assignment
    small7 = small4;
    REQUIRE(small7.num_non_zero_elements() == 4);
    val = small7.non_zero_values_data();
    col = small7.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small7.size(0); i++)
    {
      for (auto it = small7.pointers_begin()[i]; it != small7.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    // copy constructor
    auto small8(small4);
    REQUIRE(small8.num_non_zero_elements() == 4);
    val = small8.non_zero_values_data();
    col = small8.non_zero_indices2_data();
    itv = v_.begin();
    itc = c_.begin();
    for (std::size_t i = 0; i < small8.size(0); i++)
    {
      for (auto it = small8.pointers_begin()[i]; it != small8.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }
  }

  // ordered
  v_ = {10, 9, 3, 1};
  c_ = {1, 2, 1, 3};

  {
    csr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, 2, A);
    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    itv      = v_.begin();
    itc      = c_.begin();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();

    std::array<IndxType, 2> range = {0, 2};
    auto small2                   = small[range];
    REQUIRE(small2.num_non_zero_elements() == 2);
    val     = small2.non_zero_values_data();
    col     = small2.non_zero_indices2_data();
    itv     = v_.begin();
    itc     = c_.begin();
    auto i0 = small2.pointers_begin()[0];
    for (std::size_t i = 0; i < small2.size(0); i++)
    {
      for (auto it = small2.pointers_begin()[i]; it != small2.pointers_end()[i]; it++)
      {
        REQUIRE(val[it - i0] == *(itv++));
        REQUIRE(col[it - i0] == *(itc++));
      }
    }

    range       = {2, 4};
    auto small3 = small[range];
    REQUIRE(small3.num_non_zero_elements() == 2);
    val = small3.non_zero_values_data();
    col = small3.non_zero_indices2_data();
    itv = v_.begin() + 2;
    itc = c_.begin() + 2;
    i0  = small3.pointers_begin()[0];
    for (std::size_t i = 0; i < small3.size(0); i++)
    {
      for (auto it = small3.pointers_begin()[i]; it != small3.pointers_end()[i]; it++)
      {
        REQUIRE(val[it - i0] == *(itv++));
        REQUIRE(col[it - i0] == *(itc++));
      }
    }

    node.barrier();

    std::vector<Type> v__     = {10, 3};
    std::vector<IndxType> c__ = {1, 1};

    std::array<IndxType, 4> range2 = {0, 3, 0, 2};
    auto small4                    = small[range2];
    REQUIRE(small4.num_non_zero_elements() == 2);
    REQUIRE(small4.size(0) == 3);
    REQUIRE(small4.size(1) == 2);
    val     = small4.non_zero_values_data();
    col     = small4.non_zero_indices2_data();
    itv     = v__.begin();
    itc     = c__.begin();
    auto i1 = small4.pointers_begin()[0];
    for (std::size_t i = 0; i < small4.size(0); i++)
    {
      for (auto it = small4.pointers_begin()[i]; it != small4.pointers_end()[i]; it++)
      {
        REQUIRE(val[it - i1] == *(itv++));
        REQUIRE(col[it - i1] == *(itc++));
      }
    }

    node.barrier();

    v__ = {9, 1};
    c__ = {2, 3};

    range2      = {0, 4, 2, 4};
    auto small5 = small[range2];
    REQUIRE(small5.num_non_zero_elements() == 2);
    REQUIRE(small5.size(0) == 4);
    REQUIRE(small5.size(1) == 4);
    val = small5.non_zero_values_data();
    col = small5.non_zero_indices2_data();
    itv = v__.begin();
    itc = c__.begin();
    i1  = small5.pointers_begin()[0];
    for (std::size_t i = 0; i < small5.size(0); i++)
    {
      for (auto it = small5.pointers_begin()[i]; it != small5.pointers_end()[i]; it++)
      {
        REQUIRE(val[it - i1] == *(itv++));
        REQUIRE(col[it - i1] == *(itc++));
      }
    }
  }
  {
    csr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    itv      = v_.begin();
    itc      = c_.begin();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();
  }

  {
    csr_matrix small(tp_ul_ul{4, 4}, tp_ul_ul{0, 0}, max_non_zero_per_row, A);
    if (serial || node.rank() == 0)
      small[3][3] = 1;
    if (serial || node.rank() == 0)
      small[0][2] = 9;
    node.barrier();
    if (serial || node.rank() == node.size() - 1)
      small[2][1] = 3;
    if (serial || node.rank() == node.size() - 1)
      small[0][1] = 10;
    node.barrier();

    REQUIRE(small.num_non_zero_elements() == 4);
    auto val = small.non_zero_values_data();
    auto col = small.non_zero_indices2_data();
    itv      = v_.begin();
    itc      = c_.begin();
    for (std::size_t i = 0; i < small.size(0); i++)
    {
      for (auto it = small.pointers_begin()[i]; it != small.pointers_end()[i]; it++)
      {
        REQUIRE(val[it] == *(itv++));
        REQUIRE(col[it] == *(itc++));
      }
    }

    node.barrier();
  }
};

TEST_CASE("csr_matrix_serial", "[csr]")
{
  // serial
  {
    using Type    = double;
    using Alloc   = std::allocator<Type>;
    using is_root = ma::sparse::null_is_root<Alloc>;
    test_csr_matrix_shm_allocator<Type, int, std::size_t, Alloc, is_root>(Alloc(), true);
    test_csr_matrix_shm_allocator<Type, int, int, Alloc, is_root>(Alloc(), true);
  }
  {
    using Type    = std::complex<double>;
    using Alloc   = std::allocator<Type>;
    using is_root = ma::sparse::null_is_root<Alloc>;
    test_csr_matrix_shm_allocator<Type, int, std::size_t, Alloc, is_root>(Alloc(), true);
    test_csr_matrix_shm_allocator<Type, int, int, Alloc, is_root>(Alloc(), true);
  }
}

TEST_CASE("csr_matrix_shm", "[csr]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  mpi3::shared_communicator node(world.split_shared());

  {
    using Type    = double;
    using Alloc   = shm::allocator_shm_ptr_with_raw_ptr_dispatch<Type>;
    using is_root = ma::sparse::is_root;
    test_csr_matrix_shm_allocator<Type, int, std::size_t, Alloc, is_root>(Alloc(node), false);
    test_csr_matrix_shm_allocator<Type, int, int, Alloc, is_root>(Alloc(node), false);
  }

  {
    using Type    = std::complex<double>;
    using Alloc   = shm::allocator_shm_ptr_with_raw_ptr_dispatch<Type>;
    using is_root = ma::sparse::is_root;
    test_csr_matrix_shm_allocator<Type, int, std::size_t, Alloc, is_root>(Alloc(node), false);
    test_csr_matrix_shm_allocator<Type, int, int, Alloc, is_root>(Alloc(node), false);
  }
}


//#define TEST_CSR_LARGE_MEMORY
#ifdef TEST_CSR_LARGE_MEMORY
TEST_CASE("csr_matrix_shm_large_memory", "[csr]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  mpi3::shared_communicator node(world.split_shared());

  using Type    = std::complex<double>;
  using Alloc   = shm::allocator_shm_ptr_with_raw_ptr_dispatch<Type>;
  using is_root = ma::sparse::is_root;

  using ucsr_matrix = ma::sparse::ucsr_matrix<Type, int, size_t, Alloc, is_root>;
  using csr_matrix  = ma::sparse::csr_matrix<Type, int, size_t, Alloc, is_root>;

  Alloc A(node);
  ucsr_matrix umat({400000, 400000}, tp_ul_ul{0, 0}, 7000, A);
  world.barrier();
  if (node.root())
  {
    std::cout << " capacity: " << umat.capacity() << std::endl;
    umat.emplace({399999, 0}, Type(1));
    std::cout << " pbegin[399999]: " << umat.pointers_begin()[399999] << std::endl;
  }
  world.barrier();
  Alloc B(node);
  ucsr_matrix umat2({400000, 400000}, tp_ul_ul{0, 0}, 7000, B);
  world.barrier();
  if (node.root())
  {
    std::cout << " capacity: " << umat2.capacity() << std::endl;
    umat2.emplace({399999, 0}, Type(1));
    std::cout << " pbegin[399999]: " << umat2.pointers_begin()[399999] << std::endl;
  }
  world.barrier();
}
#endif


} // namespace qmcplusplus
