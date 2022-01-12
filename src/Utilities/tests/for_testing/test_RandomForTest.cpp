//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "RandomForTest.cpp"

/** \file
 */
namespace qmcplusplus
{

template<typename T = void>
struct RngReferenceVector;

template<>
struct RngReferenceVector<double>
{
  std::vector<double> vec{0.6121701794, 0.1207577838, 0.1690697568, 0.3017894523, 0.4360590181};
};

template<>
struct RngReferenceVector<float>
{
  std::vector<float> vec{0.61217f, 0.12076f,0.16907f,0.30179f,0.43606f};
};

template<>
struct RngReferenceVector<std::complex<double>>
{
std::vector<std::complex<double>> vec{{0.6121701794, 0.1207577838}, {0.1690697568, 0.3017894523}};
};

template<>
struct RngReferenceVector<std::complex<float>>
{
  std::vector<std::complex<float>> vec{{0.61217f, 0.12076f},{0.16907f,0.30179f}};
};

template<typename T>
struct TestFillBufferRngReal : public RngReferenceVector<T>
{
  using ref = RngReferenceVector<T>;
  void operator()()
  {
    testing::RandomForTest<T> rng_for_test;
    std::vector<T> rng_vector(5, 0.0);
    rng_for_test.fillBufferRng(rng_vector.data(), 5);
    auto checkArray = [](auto A, auto B, int n) {
      for (int i = 0; i < n; ++i)
      {
        CHECK(A[i] == Approx(B[i]));
      }
    };
    checkArray(rng_vector.data(), ref::vec.data(), 5);
  }
};

template<typename T>
struct TestFillVecRngReal : public RngReferenceVector<T>
{
  using ref = RngReferenceVector<T>;
  void operator()()
  {
    testing::RandomForTest<T> rng_for_test;
    std::vector<T> rng_vector(5);
    rng_for_test.fillVecRng(rng_vector);
    std::vector<T> reference{0.120758, 0.301789, 0.853906, 0.297716, 0.377862};
    auto checkArray = [](auto A, auto B, int n) {
      for (int i = 0; i < n; ++i)
      {
        CHECK(A[i] == Approx(B[i]));
      }
    };
    checkArray(rng_vector.data(), ref::vec.data(), 5);
  }
};

template<typename T>
struct TestGetVecRngReal : public RngReferenceVector<T>
{
  using ref = RngReferenceVector<T>;
  void operator()()
  {
    testing::RandomForTest<T> rng_for_test;
    std::vector<T> rng_vector = rng_for_test.getRngVec(5);
    std::vector<T> reference{0.120758, 0.301789, 0.853906, 0.297716, 0.377862};
    auto checkArray = [](auto A, auto B, int n) {
      for (int i = 0; i < n; ++i)
      {
        CHECK(A[i] == Approx(B[i]));
      }
    };
    checkArray(rng_vector.data(), ref::vec.data(), 5);
  }
};

template<typename T>
struct TestFillBufferRngComplex : public RngReferenceVector<T>
{
  using ref = RngReferenceVector<T>;
  using VT = typename T::value_type;
  void operator()()
  {
    testing::RandomForTest<VT> rng_for_test;
    std::vector<T> rng_vector(2);
    rng_for_test.fillBufferRng(rng_vector.data(), 2);
    auto checkArray = [](auto A, auto B, int n) {
      for (int i = 0; i < n; ++i)
      {
        CHECK(A[i] == ComplexApprox(B[i]));
      }
    };
    checkArray(rng_vector.data(), ref::vec.data(), 2);
  }
};


template<typename T>
struct TestFillVecRngComplex : public RngReferenceVector<T>
{
  using ref = RngReferenceVector<T>;
  using VT = typename T::value_type;
  void operator()()
  {
    testing::RandomForTest<VT> rng_for_test;
    std::vector<T> rng_vector(2);
    rng_for_test.fillVecRng(rng_vector);
    auto checkArray = [](auto A, auto B, int n) {
      for (int i = 0; i < n; ++i)
      {
        CHECK(A[i] == ComplexApprox(B[i]));
      }
    };
    checkArray(rng_vector.data(), ref::vec.data(), 2);
  }
};

template<typename T>
struct TestGetVecRngComplex : public RngReferenceVector<T>
{
  using ref = RngReferenceVector<T>;
  using VT = typename T::value_type;

  void operator()()
  {
    using ref = RngReferenceVector<T>;
    testing::RandomForTest<VT> rng_for_test;
    std::vector<T> rng_vector = rng_for_test.getRngVecComplex(2);
    auto checkArray = [](auto A, auto B, int n) {
      for (int i = 0; i < n; ++i)
      {
        CHECK(A[i] == ComplexApprox(B[i]));
      }
    };
    checkArray(rng_vector.data(), ref::vec.data(), 2);
  }
};

TEST_CASE("RandomForTest_fillBufferRng_real", "[utilities][for_testing]")
{
  TestFillBufferRngReal<float> tfbrr_f;
  tfbrr_f();
  TestFillBufferRngReal<double> tfbrr_d;
  tfbrr_d();
}

TEST_CASE("RandomForTest_fillVectorRngReal", "[utilities][for_testing]")
{
  TestFillVecRngReal<float> tfvrr_f;
  tfvrr_f();
  TestFillVecRngReal<double> tfvrr_d;
  tfvrr_d();
}

TEST_CASE("RandomForTest_getVecRngReal", "[utilities][for_testing]")
{
  TestGetVecRngReal<float> tgvrr_f;
  tgvrr_f();
  TestGetVecRngReal<double> tgvrr_d;
  tgvrr_d();
}

TEST_CASE("RandomForTest_fillBufferRngComplex", "[utilities][for_testing]")
{
  TestFillBufferRngComplex<std::complex<float>> tfbrc_f;
  tfbrc_f();
  TestFillBufferRngComplex<std::complex<double>> tfbrc_d;
  tfbrc_d();
}

TEST_CASE("RandomForTest_fillVecRngComplex", "[utilities][for_testing]")
{
  TestFillVecRngComplex<std::complex<float>> tfvrc_f;
  tfvrc_f();
  TestFillVecRngComplex<std::complex<double>> tfvrc_d;
  tfvrc_d();
}

TEST_CASE("RandomForTest_getVecRngComplex", "[utilities][for_testing]")
{
  TestGetVecRngComplex<std::complex<float>> tgvrc_f;
  tgvrc_f();
  TestGetVecRngComplex<std::complex<double>> tgvrc_d;
  tgvrc_d();
}


TEST_CASE("RandomForTest_call_operator", "[utilities][for_testing]")
{
  testing::RandomForTest<double> rng_for_test;
  double test = rng_for_test();
  CHECK(test == Approx(0.6121701794));
}


} // namespace qmcplusplus
