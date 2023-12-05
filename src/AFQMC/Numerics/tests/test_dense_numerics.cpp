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


//#include "catch.hpp"

#include "catch.hpp"
#include "Configuration.h"

// Always test the fallback code, regardless of MKL definition
//#undef HAVE_MKL
#define MKL_INT int
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x) \
  {                  \
    std::cout << x;  \
    throw;           \
  }

#include <stdio.h>
#include <string>
#include <complex>
#include <vector>
#include "multi/array.hpp"
#include "multi/array_ref.hpp"

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "AFQMC/Memory/custom_pointers.hpp"
#endif

#include "AFQMC/Matrix/tests/matrix_helpers.h"

// Include the templates directly so all the needed types get instantiated
//  and so the undef of HAVE_MKL has an effect
#include "AFQMC/Numerics/ma_operations.hpp"

using boost::multi::array;
using boost::multi::array_ref;
using std::complex;
using std::cout;
using std::endl;
using std::string;
using std::vector;
template<std::ptrdiff_t D>
using iextensions = typename boost::multi::iextensions<D>;

namespace qmcplusplus
{
// tests dispatching through ma_operations
void test_dense_matrix_mult()
{
  {
    vector<double> m = {9., 24., 30., 4., 10., 12., 14., 16., 36.};
    array_ref<double, 2> M(m.data(), {3, 3});
    REQUIRE(M.num_elements() == m.size());
    vector<double> x = {1., 2., 3.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y(3);
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));

    using ma::T;
    ma::product(M, X, Y); // Y := M X

    vector<double> mx = {147., 60., 154.};
    array_ref<double, 1> MX(mx.data(), iextensions<1u>(mx.size()));
    verify_approx(MX, Y);
  }
  {
    vector<double> m = {9., 24., 30., 2., 4., 10., 12., 1., 14., 16., 36., 20.};
    array_ref<double, 2> M(m.data(), {3, 4});
    REQUIRE(M.num_elements() == m.size());
    vector<double> x = {1., 2., 3., 4.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y(3);
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));

    using ma::T;
    ma::product(M, X, Y); // Y := M X

    vector<double> mx = {155., 64., 234.};
    array_ref<double, 1> MX(mx.data(), iextensions<1u>(mx.size()));
    verify_approx(MX, Y);
  }
  {
    vector<double> m = {9., 24., 30., 2., 4., 10., 12., 1., 14., 16., 36., 20.};
    array_ref<double, 2> M(m.data(), {3, 4});
    REQUIRE(M.num_elements() == m.size());
    vector<double> x = {1., 2., 3.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y(4);
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));

    using ma::T;
    ma::product(T(M), X, Y); // Y := T(M) X

    vector<double> mx = {59., 92., 162., 64.};
    array_ref<double, 1> MX(mx.data(), iextensions<1u>(mx.size()));
    verify_approx(MX, Y);
  }
  {
    vector<double> m = {9., 24., 30., 9., 4., 10., 12., 7., 14., 16., 36., 1.};
    array_ref<double, 2> M(m.data(), {3, 4});
    vector<double> x = {1., 2., 3., 4.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y = {4., 5., 6.};
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));
    ma::product(M, X, Y); // y := M x

    vector<double> y2 = {183., 88., 158.};
    array_ref<double, 1> Y2(y2.data(), iextensions<1u>(y2.size()));
    verify_approx(Y, Y2);
  }

  {
    vector<double> m = {1., 2., 1., 2., 5., 8., 1., 8., 9.};
    array_ref<double, 2> M(m.data(), {3, 3});
    REQUIRE(ma::is_hermitian(M));
  }
  {
    vector<double> m = {
        1., 0., 2., 0., 1., 0., 2., 0., 5., 0., 8., -1., 1., 0., 8., 1., 9., 0.,
    };
    array_ref<complex<double>, 2> M(reinterpret_cast<complex<double>*>(m.data()), {3, 3});
    REQUIRE(ma::is_hermitian(M));
  }
  {
    vector<double> m = {1., 2., 1., 2., 5., 8., 1., 8., 9.};
    array_ref<double, 2> M(m.data(), {3, 3});
    REQUIRE(ma::is_hermitian(M));
  }
  {
    vector<double> a = {1., 0., 1., 3., 5., 8., 4., 8., 9.};
    array_ref<double, 2> A(a.data(), {3, 3});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {6., 2., 8., 9., 5., 5., 1., 7., 9.};
    array_ref<double, 2> B(b.data(), {3, 3});
    REQUIRE(B.num_elements() == b.size());

    vector<double> c(9);
    array_ref<double, 2> C(c.data(), {3, 3});
    REQUIRE(C.num_elements() == c.size());

    ma::product(A, B, C);

    vector<double> ab = {7., 9., 17., 71., 87., 121., 105., 111., 153.};
    array_ref<double, 2> AB(ab.data(), {3, 3});
    REQUIRE(AB.num_elements() == ab.size());

    verify_approx(C, AB);

    using ma::N;
    ma::product(N(A), N(B), C); // same as ma::product(A, B, C);
    verify_approx(C, AB);

    using ma::T;

    ma::product(T(A), B, C);
    vector<double> atb = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    array_ref<double, 2> AtB(atb.data(), {3, 3});
    verify_approx(C, AtB);

    ma::product(A, T(B), C);
    vector<double> abt = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
    array_ref<double, 2> ABt(abt.data(), {3, 3});
    verify_approx(C, ABt);

    ma::product(T(A), T(B), C);
    vector<double> atbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
    array_ref<double, 2> AtBt(atbt.data(), {3, 3});
    verify_approx(C, AtBt);

    using ma::H;
    ma::product(H(A), T(B), C);
    vector<double> ahbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
    array_ref<double, 2> AhBt(ahbt.data(), {3, 3});
    verify_approx(C, AhBt);

    ma::product(A, H(B), C);
    vector<double> abh = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
    array_ref<double, 2> ABh(abh.data(), {3, 3});
    verify_approx(C, ABh);
  }

  {
    vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    array_ref<double, 2> A(a.data(), {3, 3});
    REQUIRE(A.num_elements() == a.size());
    array<double, 2> B = A;
    ma::invert(B, 0.0);

    array<double, 2> Id({3, 3});
    ma::set_identity(Id);

    array<double, 2> Id2({3, 3});
    ma::product(A, B, Id2);

    verify_approx(Id, Id2);
  }
  {
    std::vector<double> WORK;
    array<double, 1> TAU(iextensions<1u>{3});

    vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    array_ref<double, 2> A(a.data(), {3, 3});
    REQUIRE(A.num_elements() == a.size());
    WORK.resize(std::max(ma::gelqf_optimal_workspace_size(A), ma::glq_optimal_workspace_size(A)));
    ma::gelqf(A, TAU, WORK);
    ma::glq(A, TAU, WORK);

    array<double, 2> Id({3, 3});
    ma::set_identity(Id);

    using ma::H;
    array<double, 2> Id2({3, 3});
    ma::product(H(A), A, Id2);

    verify_approx(Id, Id2);
  }
  {
    std::vector<double> WORK;
    array<double, 1> TAU(iextensions<1u>{4});

    vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129., 10., 23., 35.};
    array_ref<double, 2> A(a.data(), {4, 3});
    REQUIRE(A.num_elements() == a.size());
    WORK.resize(std::max(ma::gelqf_optimal_workspace_size(A), ma::glq_optimal_workspace_size(A)));
    ma::gelqf(A, TAU, WORK);
    ma::glq(A, TAU, WORK);

    array<double, 2> Id({3, 3});
    ma::set_identity(Id);

    using ma::H;
    array<double, 2> Id2({3, 3});
    ma::product(H(A), A, Id2);

    verify_approx(Id, Id2);
  }
  {
    std::vector<double> WORK;
    array<double, 1> TAU(iextensions<1u>{3});

    vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    array_ref<double, 2> A(a.data(), {3, 3});
    REQUIRE(A.num_elements() == a.size());
    WORK.resize(std::max(ma::geqrf_optimal_workspace_size(A), ma::gqr_optimal_workspace_size(A)));
    ma::geqrf(A, TAU, WORK);
    ma::gqr(A, TAU, WORK);

    array<double, 2> Id({3, 3});
    ma::set_identity(Id);

    using ma::H;
    array<double, 2> Id2({3, 3});
    ma::product(H(A), A, Id2);

    verify_approx(Id, Id2);
  }
  {
    std::vector<double> WORK;
    array<double, 1> TAU(iextensions<1u>{4});

    vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129., 10., 23., 35.};
    array_ref<double, 2> A(a.data(), {3, 4});
    REQUIRE(A.num_elements() == a.size());
    WORK.resize(std::max(ma::geqrf_optimal_workspace_size(A), ma::gqr_optimal_workspace_size(A)));
    ma::geqrf(A, TAU, WORK);
    ma::gqr(A, TAU, WORK);

    array<double, 2> Id({3, 3});
    ma::set_identity(Id);

    using ma::H;
    array<double, 2> Id2({3, 3});
    ma::product(A, H(A), Id2);

    verify_approx(Id, Id2);
  }
  {
    vector<double> a = {9., 24., 30., 45., 4., 10., 12., 12.};
    array_ref<double, 2> A(a.data(), {2, 4});
    vector<double> at = {9., 4., 24., 10., 30., 12., 45., 12.};
    array_ref<double, 2> AT(at.data(), {4, 2});
    array<double, 2> B({4, 2});
    ma::transpose(A, B);
    verify_approx(AT, B);
  }
  {
    using namespace std::complex_literals;
    vector<std::complex<double>> m_a = {1.90000, 1.40000 + 0.90000i, 0.40000 + 0.80000i, 1.40000 - 0.90000i,
                                        0.20000, 2.20000 + 0.60000i, 0.40000 - 0.80000i, 2.20000 - 0.60000i,
                                        0.60000};
    vector<std::complex<double>> m_b = {25.9622476651464 + 0.0000000000000i,  17.7794485121929 + 13.1574958765530i,
                                        11.2649352514491 + 16.4823940873968i, 17.7794485121928 - 13.1574958765530i,
                                        20.5657808536051 - 0.0000000000000i,  17.9925255171787 + 6.0065935802308i,
                                        11.2649352514491 - 16.4823940873968i, 17.9925255171787 - 6.0065935802308i,
                                        17.9429273455619 - 0.0000000000000i};

    array<std::complex<double>, 2> A({3, 3});
    array<std::complex<double>, 2> B({3, 3});

    for (int i = 0, k = 0; i < std::get<0>(A.sizes()); i++)
      for (int j = 0; j < std::get<1>(A.sizes()); j++, k++)
        A[i][j] = m_a[k];
    for (int i = 0, k = 0; i < std::get<0>(A.sizes()); i++)
      for (int j = 0; j < std::get<1>(A.sizes()); j++, k++)
        B[i][j] = m_b[k];

    array<std::complex<double>, 2> C = ma::exp(A);
    verify_approx(C, B);
  }
}

TEST_CASE("dense_ma_operations", "[matrix_operations]") { test_dense_matrix_mult(); }

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
template<class Allocator>
void test_dense_mat_vec_device(Allocator& alloc)
{
  using T = typename Allocator::value_type;

  //SECTION("mat_vec")
  {
    vector<T> m = {9., 24., 30., 4., 10., 12., 14., 16., 36.};
    vector<T> x = {1., 2., 3.};
    vector<T> y(3);

    array<T, 2, Allocator> M({3, 3}, alloc);
    array<T, 1, Allocator> X(iextensions<1u>(x.size()), alloc);
    array<T, 1, Allocator> Y(iextensions<1u>(y.size()), alloc);

    copy_n(m.data(), m.size(), M.origin());
    REQUIRE(M.num_elements() == m.size());
    copy_n(x.data(), x.size(), X.origin());
    REQUIRE(X.num_elements() == x.size());
    REQUIRE(Y.num_elements() == y.size());

    using ma::product;
    ma::product(M, X, Y); // Y := M X

    vector<T> mx = {147., 60., 154.};
    array_ref<T, 1> MX(mx.data(), iextensions<1u>(mx.size()));
    verify_approx(MX, Y);
  }

  //SECTION("mat_vec_rec")
  {
    vector<T> m = {9., 24., 30., 2., 4., 10., 12., 1., 14., 16., 36., 20.};
    vector<T> x = {1., 2., 3., 4.};
    vector<T> y(3);

    array<T, 2, Allocator> M({3, 4}, alloc);
    array<T, 1, Allocator> X(iextensions<1u>(x.size()), alloc);
    array<T, 1, Allocator> Y(iextensions<1u>(y.size()), alloc);

    copy_n(m.data(), m.size(), M.origin());
    REQUIRE(M.num_elements() == m.size());
    copy_n(x.data(), x.size(), X.origin());
    REQUIRE(X.num_elements() == x.size());
    REQUIRE(Y.num_elements() == y.size());

    ma::product(M, X, Y); // Y := M X

    vector<T> mx = {155., 64., 234.};
    array_ref<T, 1> MX(mx.data(), iextensions<1u>(mx.size()));
    verify_approx(MX, Y);
  }
  //SECTION("mat_vec_trans")
  {
    vector<T> m = {9., 24., 30., 2., 4., 10., 12., 1., 14., 16., 36., 20.};
    vector<T> x = {1., 2., 3.};
    vector<T> y(4);

    array<T, 2, Allocator> M({3, 4}, alloc);
    array<T, 1, Allocator> X(iextensions<1u>(x.size()), alloc);
    array<T, 1, Allocator> Y(iextensions<1u>(y.size()), alloc);

    copy_n(m.data(), m.size(), M.origin());
    REQUIRE(M.num_elements() == m.size());
    copy_n(x.data(), x.size(), X.origin());
    REQUIRE(X.num_elements() == x.size());
    REQUIRE(Y.num_elements() == y.size());

    ma::product(ma::T(M), X, Y); // Y := M X

    vector<T> mx = {59., 92., 162., 64.};
    array_ref<T, 1> MX(mx.data(), iextensions<1u>(mx.size()));
    verify_approx(MX, Y);
  }
  //SECTION("mat_vec_add")
  {
    vector<T> m = {9., 24., 30., 9., 4., 10., 12., 7., 14., 16., 36., 1.};
    vector<T> x = {1., 2., 3., 4.};
    vector<T> y = {4., 5., 6.};

    array<T, 2, Allocator> M({3, 4}, alloc);
    copy_n(m.data(), m.size(), M.origin());
    REQUIRE(M.num_elements() == m.size());

    array<T, 1, Allocator> X(iextensions<1u>(x.size()), alloc);
    copy_n(x.data(), x.size(), X.origin());
    REQUIRE(X.num_elements() == x.size());

    array<T, 1, Allocator> Y(iextensions<1u>(y.size()), alloc);
    REQUIRE(Y.num_elements() == y.size());

    ma::product(M, X, Y); // Y := M X

    vector<T> y2 = {183., 88., 158.};
    array_ref<T, 1> Y2(y2.data(), iextensions<1u>(y2.size()));
    verify_approx(Y, Y2);
  }
}

template<class Allocator>
void test_dense_mat_mul_device(Allocator& alloc)
{
  using T = typename Allocator::value_type;
  //SECTION("mat_herm")
  {
    vector<T> m = {1., 2., 1., 2., 5., 8., 1., 8., 9.};
    array<T, 2, Allocator> M({3, 3}, alloc);
    copy_n(m.data(), m.size(), M.origin());
    REQUIRE(M.num_elements() == m.size());
    //  REQUIRE( ma::is_hermitian(M) );
  }

  //SECTION("mat_herm_ref")
  {
    vector<T> m = {1., 2., 1., 2., 5., 8., 1., 8., 9.};
    array<T, 2, Allocator> M({3, 3}, alloc);
    copy_n(m.data(), m.size(), M.origin());
    REQUIRE(M.num_elements() == m.size());

    array_ref<T, 2, typename Allocator::pointer> Mref(M.origin(), M.extensions());
    // not yet implemented in GPU
    //    REQUIRE( ma::is_hermitian(Mref) );
  }

  //SECTION("mat_mat_op")
  {
    vector<T> a = {1., 0., 1., 3., 5., 8., 4., 8., 9.};
    vector<T> b = {6., 2., 8., 9., 5., 5., 1., 7., 9.};

    array<T, 2, Allocator> A({3, 3}, alloc);
    copy_n(a.data(), a.size(), A.origin());
    REQUIRE(A.num_elements() == a.size());

    array<T, 2, Allocator> B({3, 3}, alloc);
    copy_n(b.data(), b.size(), B.origin());
    REQUIRE(B.num_elements() == b.size());

    array<T, 2, Allocator> D({3, 3}, alloc);

    ma::product(A, B, D);

    vector<T> ab = {7., 9., 17., 71., 87., 121., 105., 111., 153.};
    array_ref<T, 2> AB(ab.data(), {3, 3});
    verify_approx(D, AB);

    ma::product(ma::N(A), ma::N(B), D);
    verify_approx(D, AB);

    ma::product(ma::T(A), B, D);
    vector<T> atb = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    array_ref<T, 2> AtB(atb.data(), {3, 3});
    verify_approx(D, AtB);

    ma::product(A, ma::T(B), D);
    vector<T> abt = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
    array_ref<T, 2> ABt(abt.data(), {3, 3});
    verify_approx(D, ABt);

    ma::product(ma::T(A), ma::T(B), D);
    vector<T> atbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
    array_ref<T, 2> AtBt(atbt.data(), {3, 3});
    verify_approx(D, AtBt);

    ma::product(ma::H(A), ma::T(B), D);
    vector<T> ahbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
    array_ref<T, 2> AhBt(ahbt.data(), {3, 3});
    verify_approx(D, AhBt);

    ma::product(A, ma::H(B), D);
    vector<T> abh = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
    array_ref<T, 2> ABh(abh.data(), {3, 3});
    verify_approx(D, ABh);
  }

  //SECTION("mat_mat_op_inv")
  {
    vector<T> a  = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    vector<T> id = {1., 0., 0., 0., 1., 0., 0., 0., 1.};

    array<T, 2, Allocator> A({3, 3}, alloc);
    copy_n(a.data(), a.size(), A.origin());
    REQUIRE(A.num_elements() == a.size());

    array<T, 2, Allocator> B({3, 3}, alloc);
    copy_n(a.data(), a.size(), B.origin());
    REQUIRE(B.num_elements() == a.size());

    array<T, 2, Allocator> I({3, 3}, alloc);

    ma::invert(B, 0.0);

    ma::product(A, B, I);

    array_ref<T, 2> Id2(id.data(), {3, 3});
    verify_approx(I, Id2);
  }
  //SECTION("mat_trans")
  {
    vector<T> a  = {9., 24., 30., 45., 4., 10., 12., 12.};
    vector<T> at = {9., 4., 24., 10., 30., 12., 45., 12.};
    array<T, 2, Allocator> A({2, 4}, alloc);
    copy_n(a.data(), a.size(), A.origin());
    REQUIRE(A.num_elements() == a.size());

    array<T, 2, Allocator> B({4, 2}, alloc);
    ma::transpose(A, B);

    array_ref<T, 2> AT(at.data(), {4, 2});
    verify_approx(AT, B);
  }
}

template<class Allocator>
void test_dense_gerf_gqr_device(Allocator& alloc)
{
  ////SECTION("mat_gerf_gqr")
  using T = typename Allocator::value_type;
  {
    vector<T> a  = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    vector<T> id = {1., 0., 0., 0., 1., 0., 0., 0., 1.};

    array<T, 2, Allocator> A({3, 3}, alloc);
    copy_n(a.data(), a.size(), A.origin());
    REQUIRE(A.num_elements() == a.size());

    array<T, 2, Allocator> Id({3, 3}, alloc);

    auto sz = std::max(ma::geqrf_optimal_workspace_size(A), ma::gqr_optimal_workspace_size(A));
    array<T, 1, Allocator> WORK(iextensions<1u>{sz}, alloc);
    array<T, 1, Allocator> TAU(iextensions<1u>{3}, alloc);

    ma::geqrf(A, TAU, WORK);
    ma::gqr(A, TAU, WORK);

    ma::product(ma::H(A), A, Id);

    array_ref<T, 2> Id2(id.data(), {3, 3});
    verify_approx(Id, Id2);
  }
  //SECTION("mat_gerf_gqr_product")
  {
    vector<T> a  = {37., 45., 59., 53., 81., 97., 87., 105., 129., 10., 23., 35.};
    vector<T> id = {1., 0., 0., 0., 1., 0., 0., 0., 1.};

    array<T, 2, Allocator> A({3, 4}, alloc);
    copy_n(a.data(), a.size(), A.origin());
    REQUIRE(A.num_elements() == a.size());

    array<T, 2, Allocator> Id({3, 3}, alloc);

    auto sz = std::max(ma::geqrf_optimal_workspace_size(A), ma::gqr_optimal_workspace_size(A));
    array<T, 1, Allocator> WORK(iextensions<1u>{sz}, alloc);
    array<T, 1, Allocator> TAU(iextensions<1u>{3}, alloc);

    ma::geqrf(A, TAU, WORK);
    ma::gqr(A, TAU, WORK);

    ma::product(A, ma::H(A), Id);

    array_ref<T, 2> Id2(id.data(), {3, 3});
    verify_approx(Id, Id2);
  }
}

template<class Allocator>
void test_dense_gerf_gqr_strided_device(Allocator& alloc)
{
  using T = typename Allocator::value_type;
  {
    vector<T> a  = {37., 45., 59., 53., 81., 97., 87., 105., 129., 10., 23., 35.};
    vector<T> id = {1., 0., 0., 0., 1., 0., 0., 0., 1.};

    array<T, 3, Allocator> A({2, 3, 4}, alloc);
    copy_n(a.data(), a.size(), A[0].origin());
    copy_n(a.data(), a.size(), A[1].origin());
    REQUIRE(A.num_elements() == 2 * a.size());

    auto sz = std::max(ma::geqrf_optimal_workspace_size(A[0]), ma::gqr_optimal_workspace_size(A[0]));
    array<T, 1, Allocator> WORK(iextensions<1u>{sz}, alloc);
    array<T, 2, Allocator> Id({3, 3}, alloc);
    using IAllocator = typename Allocator::template rebind<int>::other;
    array<int, 1, IAllocator> info(iextensions<1u>{2}, IAllocator{alloc});
    array<T, 2, Allocator> TAU({2, 4}, alloc);

    geqrfStrided(4, 3, A.origin(), 4, 12, TAU.origin(), sz, info.origin(), 2);
    gqrStrided(4, 3, 3, A.origin(), 4, 12, TAU.origin(), 4, WORK.origin(), sz, info.origin(), 2);
    for (int i = 0; i < 2; i++)
    {
      ma::product(A[i], ma::H(A[i]), Id);
      array_ref<T, 2> Id2(id.data(), {3, 3});
      verify_approx(Id, Id2);
    }
  }
}

template<class Allocator>
void test_dense_batched_gemm(Allocator& alloc)
{
  using T       = typename Allocator::value_type;
  using pointer = typename Allocator::pointer;
  {
    int nbatch                 = 3;
    array<T, 2, Allocator> a   = {{0.0, 1.0, 2.0}, {3.0, 4.0, 5.0}, {6.0, 7.0, 8.0}};
    array<T, 2, Allocator> b   = {{0.0, 1.0, 2.0}, {3.0, 4.0, 5.0}, {6.0, 7.0, 8.0}};
    array<T, 2, Allocator> res = {{15.0, 18.0, 21.0}, {42.0, 54.0, 66.0}, {69.0, 90.0, 111.0}};
    array<T, 3, Allocator> c({3, 3, 3}, 0.0, alloc);
    T alpha = 1.0;
    T beta  = 0.0;

    std::vector<pointer> A_array;
    std::vector<pointer> B_array;
    std::vector<pointer> C_array;
    for (int i = 0; i < nbatch; i++)
    {
      A_array.emplace_back(a.origin());
      B_array.emplace_back(b.origin());
      C_array.emplace_back(c[i].origin());
    }
    using ma::gemmBatched;
    gemmBatched('N', 'N', 3, 3, 3, alpha, A_array.data(), 3, B_array.data(), 3, beta, C_array.data(), 3, nbatch);
    for (int i = 0; i < nbatch; i++)
    {
      verify_approx(c[i], res);
    }
  }
}

template<class Allocator>
void test_dense_geqrf_getri_batched_device(Allocator& alloc)
{
  using T       = typename Allocator::value_type;
  using pointer = typename Allocator::pointer;
  // Test non packed data lda > M
  vector<T> a  = {37., 45., 59., 0, 53., 81., 97., 0.0, 87., 105., 129., 0.0};
  vector<T> a2 = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
  vector<T> id = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  array<T, 3, Allocator> A({2, 3, 4}, 0.0, alloc), Ai({2, 3, 3}, 0.0, alloc);

  std::vector<pointer> A_array, Ai_array;
  for (int i = 0; i < 2; i++)
  {
    copy_n(a.data(), a.size(), A[i].origin());
    A_array.emplace_back(A[i].origin());
    Ai_array.emplace_back(Ai[i].origin());
  }

  array<T, 1, Allocator> WORK(iextensions<1u>{9}, alloc);
  array<T, 2, Allocator> Id({3, 3}, alloc);
  using IAllocator = typename Allocator::template rebind<int>::other;
  array<int, 1, IAllocator> info(iextensions<1u>{2}, IAllocator{alloc});
  array<int, 1, IAllocator> piv(iextensions<1u>{2 * (3 + 1)}, IAllocator{alloc});
  array<T, 2, Allocator> B({3, 3}, 0.0, alloc), Bi({3, 3}, 0.0, alloc);
  array<int, 1, IAllocator> spiv(iextensions<1u>{4}, IAllocator{alloc});
  int status;
  //SECTION("getrf_batched")
  {
    using ma::getrfBatched;
    getrfBatched(3, A_array.data(), 4, ma::pointer_dispatch(piv.origin()), ma::pointer_dispatch(info.origin()), 2);
    copy_n(a2.data(), a2.size(), B.origin());
    using ma::getrf;
    getrf(3, 3, ma::pointer_dispatch(B.origin()), 3, ma::pointer_dispatch(spiv.data()), status,
          ma::pointer_dispatch(WORK.data()));
  }
  //SECTION("getri_batched")
  {
    using ma::getriBatched;
    getriBatched(3, A_array.data(), 4, ma::pointer_dispatch(piv.origin()), Ai_array.data(), 3,
                 ma::pointer_dispatch(info.origin()), 2);
    //SECTION("getri")
    {
      getri(3, ma::pointer_dispatch(B.origin()), 3, ma::pointer_dispatch(piv.origin()),
            ma::pointer_dispatch(WORK.origin()), 9, status);
      for (int i = 0; i < 2; i++)
      {
        verify_approx(Ai[i], B);
      }
    }
  }
  //SECTION("mat_inv")
  {
    using std::copy_n;
    copy_n(B.origin(), B.num_elements(), Bi.origin());
    copy_n(a2.data(), a2.size(), B.origin());
    array<T, 2, Allocator> out({3, 3}, 0.0, alloc);
    // note transpose to account for fortran ordering
    array_ref<T, 2> Id2(id.data(), {3, 3});
    ma::product(ma::H(B), ma::H(Bi), out);
    verify_approx(out, Id2);
  }
  //SECTION("matrix_inv_batched")
  {
    for (int i = 0; i < 2; i++)
    {
      copy_n(a.data(), a.size(), A[i].origin());
      array<T, 2, Allocator> out({3, 3}, alloc);
      ma::product(ma::H(A[i]({0, 3}, {0, 3})), ma::H(Ai[i]), out);
      array_ref<T, 2> Id2(id.data(), {3, 3});
      verify_approx(out, Id2);
    }
  }
}

TEST_CASE("dense_ma_operations_device_double", "[matrix_operations]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());
  arch::INIT(node);
  {
    using Alloc = device::device_allocator<double>;
    Alloc alloc{};
    test_dense_mat_vec_device<Alloc>(alloc);
    test_dense_mat_mul_device<Alloc>(alloc);
    test_dense_gerf_gqr_device<Alloc>(alloc);
    test_dense_batched_gemm<Alloc>(alloc);
  }
}
TEST_CASE("dense_ma_operations_device_complex", "[matrix_operations]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());
  arch::INIT(node);
  {
    using Alloc = device::device_allocator<std::complex<double>>;
    Alloc alloc{};
    test_dense_mat_vec_device<Alloc>(alloc);
    test_dense_mat_mul_device<Alloc>(alloc);
    test_dense_gerf_gqr_device<Alloc>(alloc);
    test_dense_gerf_gqr_strided_device<Alloc>(alloc);
    test_dense_geqrf_getri_batched_device<Alloc>(alloc);
    test_dense_geqrf_getri_batched_device<Alloc>(alloc);
    test_dense_batched_gemm<Alloc>(alloc);
  }
}
#endif

} // namespace qmcplusplus
