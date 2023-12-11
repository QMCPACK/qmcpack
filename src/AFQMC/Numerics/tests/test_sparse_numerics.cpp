//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Configuration.h"

// Always test the fallback code, regardless of MKL definition
#define MKL_INT int
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>

#undef APP_ABORT
#define APP_ABORT(x) \
  {                  \
    std::cout << x;  \
    throw;           \
  }

#include <iostream>
#include <vector>

#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Memory/buffer_managers.h"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

using boost::multi::array;
using boost::multi::array_ref;
using std::vector;
template<std::ptrdiff_t D>
using iextensions = typename boost::multi::iextensions<D>;

namespace qmcplusplus
{
// tests dispatching through ma_operations
template<class Allocator = std::allocator<double>>
void test_sparse_matrix_mult(Allocator const& alloc = {})
{
  using csr_matrix = ma::sparse::csr_matrix<double, int, int>;
  csr_matrix A_(std::tuple<std::size_t, std::size_t>{4, 4}, std::tuple<std::size_t, std::size_t>{0, 0}, 4);
  A_[3][3] = 1.;
  A_[2][1] = 3.;
  A_[0][1] = 9.;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  A_.remove_empty_spaces();
  ma::sparse::csr_matrix<double, int, int, Allocator> A(A_);
#else
  csr_matrix& A = A_;
#endif

  // matrix-matrix
  {
    vector<double> b = {1., 2., 1., 5., 2., 5., 8., 7., 1., 8., 9., 9., 4., 1., 2., 3.};
    array<double, 2, Allocator> B({4, 4}, alloc);
    using std::copy_n;
    copy_n(b.data(), b.size(), B.origin());
    REQUIRE(B.num_elements() == b.size());

    array<double, 2, Allocator> C({4, 4}, 0.0, alloc);
    REQUIRE(C.num_elements() == 16);

    ma::product(A, B, C); // C = A*B

    vector<double> c2 = {18., 45., 72., 63., 0., 0., 0., 0., 6., 15., 24., 21., 4., 1., 2., 3.};
    array_ref<double, 2> C2(c2.data(), {4, 4});
    REQUIRE(C2.num_elements() == c2.size());
    verify_approx(C, C2);

    using ma::T;
    ma::product(T(A), B, C); // D = T(A)*B

    vector<double> d2 = {0, 0, 0, 0, 12, 42, 36, 72, 0, 0, 0, 0, 4, 1, 2, 3};
    array_ref<double, 2> D2(d2.data(), {4, 4});
    REQUIRE(D2.num_elements() == d2.size());
    verify_approx(C, D2);
  }

  //// matrix-vector
  {
    vector<double> b = {1., 2., 1., 4.};
    array<double, 1, Allocator> B(iextensions<1u>{4}, alloc);
    using std::copy_n;
    copy_n(b.data(), b.size(), B.origin());
    REQUIRE(B.num_elements() == b.size());

    array<double, 1, Allocator> C(iextensions<1u>{4}, alloc);
    REQUIRE(C.num_elements() == 4);

    ma::product(A, B, C); // C = A*B

    vector<double> c2 = {18., 0., 6., 4.};
    array_ref<double, 1> C2(c2.data(), iextensions<1u>{4});
    REQUIRE(C2.num_elements() == c2.size());
    verify_approx(C, C2);

    using ma::T;
    ma::product(T(A), B, C); // D = T(A)*B

    vector<double> d2 = {0., 12., 0., 4.};
    array_ref<double, 1> D2(d2.data(), iextensions<1u>{4});
    REQUIRE(D2.num_elements() == d2.size());
    verify_approx(C, D2);
  }

#if !defined(ENABLE_CUDA) && !defined(ENABLE_HIP)
  // test that everything is fine after this
  A.remove_empty_spaces();
  // matrix-matrix
#endif

  {
    vector<double> b = {1., 2., 1., 5., 2., 5., 8., 7., 1., 8., 9., 9., 4., 1., 2., 3.};
    array<double, 2, Allocator> B({4, 4}, alloc);
    using std::copy_n;
    copy_n(b.data(), b.size(), B.origin());
    REQUIRE(B.num_elements() == b.size());

    array<double, 2, Allocator> C({4, 4}, alloc);
    REQUIRE(C.num_elements() == 16);

    ma::product(A, B, C); // C = A*B

    vector<double> c2 = {18., 45., 72., 63., 0., 0., 0., 0., 6., 15., 24., 21., 4., 1., 2., 3.};
    array_ref<double, 2> C2(c2.data(), {4, 4});
    REQUIRE(C2.num_elements() == c2.size());
    verify_approx(C, C2);

    using ma::T;
    ma::product(T(A), B, C); // D = T(A)*B
    vector<double> d2 = {0, 0, 0, 0, 12, 42, 36, 72, 0, 0, 0, 0, 4, 1, 2, 3};
    array_ref<double, 2> D2(d2.data(), {4, 4});
    REQUIRE(D2.num_elements() == d2.size());
    verify_approx(C, D2);
  }

  // matrix-vector
  {
    vector<double> b = {1., 2., 1., 4.};
    array<double, 1, Allocator> B(iextensions<1u>{4}, alloc);
    using std::copy_n;
    copy_n(b.data(), b.size(), B.origin());
    REQUIRE(B.num_elements() == b.size());

    array<double, 1, Allocator> C(iextensions<1u>{4}, alloc);
    REQUIRE(C.num_elements() == 4);

    ma::product(A, B, C); // C = A*B

    vector<double> c2 = {18., 0., 6., 4.};
    array_ref<double, 1> C2(c2.data(), iextensions<1u>{4});
    verify_approx(C, C2);

    using ma::T;
    ma::product(T(A), B, C); // D = T(A)*B

    vector<double> d2 = {0., 12., 0., 4.};
    array_ref<double, 1> D2(d2.data(), iextensions<1u>{4});
    REQUIRE(D2.num_elements() == d2.size());
    verify_approx(C, D2);
  }
}

TEST_CASE("sparse_ma_operations", "[matrix_operations]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
  using Alloc = device::device_allocator<double>;
  afqmc::setup_memory_managers(node, 10uL * 1024uL * 1024uL);
  test_sparse_matrix_mult<Alloc>();
#else
  afqmc::setup_memory_managers(node, 10uL * 1024uL * 1024uL);
#endif
  test_sparse_matrix_mult();
  afqmc::release_memory_managers();
}

} // namespace qmcplusplus
