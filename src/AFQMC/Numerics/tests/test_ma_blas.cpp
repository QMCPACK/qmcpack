//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Configuration.h"

#undef APP_ABORT
#define APP_ABORT(x) \
  {                  \
    std::cout << x;  \
    throw;           \
  }

#include <vector>
#include <iostream>


#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Numerics/ma_blas.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

using boost::multi::array;
using boost::multi::array_ref;
using std::vector;
template<std::ptrdiff_t D>
using iextensions = typename boost::multi::iextensions<D>;

namespace qmcplusplus
{
void ma_blas_tests()
{
  vector<double> v = {1., 2., 3.};
  {
    array_ref<double, 1> V(v.data(), iextensions<1u>(v.size()));
    ma::scal(2., V);
    {
      vector<double> v2 = {2., 4., 6.};
      array_ref<double, 1> V2(v2.data(), iextensions<1u>(v2.size()));
      verify_approx(V, V2);
    }
  }

  vector<double> m = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
  array_ref<double, 2> M(m.data(), {3, 3});
  REQUIRE(M.num_elements() == m.size());
  ma::scal(2., M[2]);
  {
    vector<double> m2 = {1., 2., 3., 4., 5., 6., 14., 16., 18.};
    array_ref<double, 2> M2(m2.data(), {3, 3});
    verify_approx(M, M2);
  }

  ma::scal(2., M({0, 3}, 2));
  {
    vector<double> m2 = {1., 2., 6., 4., 5., 12., 14., 16., 36.};
    array_ref<double, 2> M2(m2.data(), {3, 3});
    verify_approx(M, M2);
  }
  ma::scal(2., M({0, 2}, 1));
  {
    vector<double> m2 = {1., 4., 6., 4., 10., 12., 14., 16., 36.};
    array_ref<double, 2> M2(m2.data(), {3, 3});
    verify_approx(M, M2);
  }
  ma::axpy(2., M[1], M[0]); // M[0] += a*M[1]
  {
    vector<double> m2 = {9., 24., 30., 4., 10., 12., 14., 16., 36.};
    array_ref<double, 2> M2(m2.data(), {3, 3});
    verify_approx(M, M2);
  }
  {
    vector<double> m = {9., 24., 30., 9., 4., 10., 12., 7., 14., 16., 36., 1.};
    array_ref<double, 2> M(m.data(), {3, 4});
    REQUIRE(M[2][0] == 14.);
    vector<double> x = {1., 2., 3., 4.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y = {4., 5., 6.};
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));
    ma::gemv<'T'>(1., M, X, 0., Y); // y := M x

    vector<double> y2 = {183., 88., 158.};
    array_ref<double, 1> Y2(y2.data(), iextensions<1u>(y2.size()));
    verify_approx(Y, Y2);
  }
  {
    vector<double> m = {9., 24., 30., 9., 4., 10., 12., 7., 14., 16., 36., 1.};
    array_ref<double, 2> M(m.data(), {3, 4});

    vector<double> x = {1., 2.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y = {4., 5.};
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));

    auto const& mm = M({0, 2}, {0, 2});
    ma::gemv<'T'>(1., M({0, 2}, {0, 2}), X, 0., Y); // y := M x

    vector<double> y2 = {57., 24.};
    array_ref<double, 1> Y2(y2.data(), iextensions<1u>(y2.size()));
    verify_approx(Y, Y2);
  }
  {
    vector<double> m = {9., 24., 30., 4., 10., 12., 14., 16., 36., 4., 9., 1.};
    array_ref<double, 2> M(m.data(), {4, 3});
    REQUIRE(M[2][0] == 14.);
    vector<double> x = {1., 2., 3.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y = {4., 5., 6., 7.};
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));
    ma::gemv<'T'>(1., M, X, 0., Y); // y := M x
    vector<double> y2 = {147., 60., 154., 25.};
    array_ref<double, 1> Y2(y2.data(), iextensions<1u>(y2.size()));
    verify_approx(Y, Y2);
  }
  {
    vector<double> m = {9., 24., 30., 9., 4., 10., 12., 7., 14., 16., 36., 1.};
    array_ref<double, 2> M(m.data(), {3, 4});
    REQUIRE(M[2][0] == 14.);
    vector<double> x = {1., 2., 3.};
    array_ref<double, 1> X(x.data(), iextensions<1u>(x.size()));
    vector<double> y = {4., 5., 6., 7.};
    array_ref<double, 1> Y(y.data(), iextensions<1u>(y.size()));
    ma::gemv<'N'>(1., M, X, 0., Y); // y := M^T x
    vector<double> y2 = {59., 92., 162., 26.};
    array_ref<double, 1> Y2(y2.data(), iextensions<1u>(y2.size()));
    verify_approx(Y, Y2);
  }
  {
    vector<double> a = {9., 24., 30., 2., 4., 10., 12., 9.};
    array_ref<double, 2> A(a.data(), {2, 4});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {9., 24., 6., 8., 4., 10., 2., 5., 14., 16., 9., 0.};
    array_ref<double, 2> B(b.data(), {3, 4});
    REQUIRE(B.num_elements() == b.size());

    vector<double> c(6);
    array_ref<double, 2> C(c.data(), {3, 2});
    REQUIRE(C.num_elements() == c.size());

    ma::gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)

    vector<double> tab = {853., 420., 346., 185., 780., 324.};
    array_ref<double, 2> TAB(tab.data(), {3, 2});
    REQUIRE(TAB.num_elements() == tab.size());

    verify_approx(C, TAB);
  }
  {
    vector<double> a = {9., 24., 30., 4., 10., 12.};
    array_ref<double, 2> A(a.data(), {2, 3});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {
        9., 24., 6., 8., 4., 10., 2., 5.,
    };
    array_ref<double, 2> B(b.data(), {2, 4});
    REQUIRE(B.num_elements() == b.size());

    vector<double> c(12);
    array_ref<double, 2> C(c.data(), {4, 3});
    REQUIRE(C.num_elements() == c.size());


    ma::gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B

    vector<double> tab = {97., 256., 318., 256., 676., 840., 62., 164., 204., 92., 242., 300.};
    array_ref<double, 2> TAB(tab.data(), {4, 3});
    REQUIRE(TAB.num_elements() == tab.size());

    verify_approx(C, TAB);
  }
  {
    vector<double> a = {9., 24., 30., 4., 10., 12., 3., 11., 45., 1., 2., 6.};
    array_ref<double, 2> A(a.data(), {4, 3});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {9., 24., 6., 8., 4., 10., 2., 5., 14., 16., 9., 0.};
    array_ref<double, 2> B(b.data(), {3, 4});
    REQUIRE(B.num_elements() == b.size());

    vector<double> c(9);
    array_ref<double, 2> C(c.data(), {3, 3});
    REQUIRE(C.num_elements() == c.size());

    ma::gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)

    vector<double> tab = {203., 538., 876., 87., 228., 360., 217., 595., 1017.};
    array_ref<double, 2> TAB(tab.data(), {3, 3});
    REQUIRE(TAB.num_elements() == tab.size());
    verify_approx(C, TAB);
  }
  {
    vector<double> a = {9., 24., 30., 4., 10., 12., 14., 16., 36.};
    array_ref<double, 2> A(a.data(), {3, 3});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {9., 24., 4., 4., 10., 1., 14., 16., 3.};
    array_ref<double, 2> B(b.data(), {3, 3});
    REQUIRE(B.num_elements() == b.size());
    vector<double> c(9);
    array_ref<double, 2> C(c.data(), {3, 3});
    REQUIRE(C.num_elements() == c.size());


    ma::gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
    ma::gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)
    ma::gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B
    ma::gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
  }
  {
    vector<double> a = {9., 24., 30., 4., 10., 12.};
    array_ref<double, 2> A(a.data(), {2, 3});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {9., 24., 6., 8., 4., 10., 2., 5., 14., 16., 9., 0.};
    array_ref<double, 2> B(b.data(), {3, 4});
    REQUIRE(B.num_elements() == b.size());

    vector<double> c(8);
    array_ref<double, 2> C(c.data(), {4, 2});
    REQUIRE(C.num_elements() == c.size());

    ma::gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B

    vector<double> tab = {597, 244, 936, 388, 372, 152, 192, 82};
    array_ref<double, 2> TAB(tab.data(), {4, 2});
    REQUIRE(TAB.num_elements() == tab.size());
    verify_approx(C, TAB);
  }
  {
    vector<double> a = {9., 24., 30., 45., 4., 10., 12., 12.};
    array_ref<double, 2> A(a.data(), {2, 4});
    REQUIRE(A.num_elements() == a.size());
    vector<double> b = {9., 24., 56., 4., 10., 78., 14., 16., 90., 6., 9., 18.};
    array_ref<double, 2> B(b.data(), {4, 3});
    REQUIRE(B.num_elements() == b.size());
    vector<double> c = {9., 24., 8., 4., 10., 9.};
    array_ref<double, 2> C(c.data(), {3, 2});
    REQUIRE(C.num_elements() == c.size());
    ma::gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
  }
}

TEST_CASE("test_ma_blas", "[matrix_operations]") { ma_blas_tests(); }

} // namespace qmcplusplus
