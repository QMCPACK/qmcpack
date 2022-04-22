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

#include <type_traits>
#include "Platforms/CPU/BLAS.hpp"
#include "Platforms/CPU/SIMD/simd.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/for_testing/RandomForTest.h"
#include "type_traits/complex_help.hpp"
#include "type_traits/type_mapping.hpp"

namespace qmcplusplus
{

/** Get the right type for rng in case of complex T.
 *  Probably a more elegant way to do this especially in c++17
 */
template<typename T>
using RngValueType = typename std::disjunction<OnTypesEqual<T, float, float>,
                                               OnTypesEqual<T, double, double>,
                                               OnTypesEqual<T, std::complex<float>, float>,
                                               OnTypesEqual<T, std::complex<double>, double>,
                                               default_type<void>>::type;

namespace testing
{

template<typename T, typename = typename std::enable_if<std::is_floating_point<T>::value, void>::type>
void makeRngSpdMatrix(testing::RandomForTest<RngValueType<T>>& rng, Matrix<T>& mat_spd)
{
  int n = mat_spd.rows();
  Matrix<T> mat_a;
  mat_a.resize(n, n);
  rng.fillBufferRng(mat_a.data(), mat_a.size());
  Matrix<T> mat_a_t(mat_a);
  Matrix<T> mat_c;
  mat_c.resize(n, n);
  Matrix<T> mat_u;
  mat_u.resize(n, n);
  Matrix<T> mat_v;
  mat_v.resize(n, n);
  Matrix<T> mat_diag;
  mat_diag.resize(n, n);
  mat_diag = 0.0;
  std::vector<T> singular_values(n);
  std::vector<T> work_vec(5 * n);
  int info;
  const char trans   = 't';
  const char notrans = 'n';
  const char all     = 'a';
  BLAS::gemm(trans, notrans, n, n, n, 1.0, mat_a.data(), n, mat_a_t.data(), n, 0.0, mat_c.data(), n);
  LAPACK::gesvd(all, all, n, n, mat_c.data(), n, singular_values.data(), mat_u.data(), n, mat_v.data(), n,
                work_vec.data(), work_vec.size(), info);
  if (info == 0)
  {
    for (int i = 0; i < n; ++i)
      mat_diag(i, i) = rng();
    mat_diag += 1.0;
    // reuse mat_a_t for this intermediate resutls
    BLAS::gemm(notrans, notrans, n, n, n, 1.0, mat_u.data(), n, mat_diag.data(), n, 0.0, mat_a_t.data(), n);
    // then reuse mat_u
    BLAS::gemm(notrans, notrans, n, n, n, 1.0, mat_a_t.data(), n, mat_v.data(), n, 0.0, mat_u.data(), n);
  }
  simd::transpose(mat_u.data(), n, mat_u.cols(), mat_spd.data(), n, n);
}

template<typename T, IsComplex<T> = true>
void makeRngSpdMatrix(testing::RandomForTest<RngValueType<T>>& rng, Matrix<T>& mat_spd)
{
  int n = mat_spd.rows();
  Matrix<T> mat_a;
  mat_a.resize(n, n);
  rng.fillBufferRng(mat_a.data(), mat_a.size());
  Matrix<T> mat_a_t(mat_a);
  Matrix<T> mat_c;
  mat_c.resize(n, n);
  Matrix<T> mat_u;
  mat_u.resize(n, n);
  Matrix<T> mat_v;
  mat_v.resize(n, n);
  Matrix<T> mat_diag;
  mat_diag.resize(n, n);
  mat_diag = 0.0;
  std::vector<typename T::value_type> singular_values(n);
  std::vector<T> work_vec(5 * n);
  std::vector<typename T::value_type> real_work_vec(5 * n);
  int info;
  const char trans   = 't';
  const char notrans = 'n';
  const char all     = 'a';
  BLAS::gemm(trans, notrans, n, n, n, 1.0, mat_a.data(), n, mat_a_t.data(), n, 0.0, mat_c.data(), n);
  LAPACK::gesvd(all, all, n, n, mat_c.data(), n, singular_values.data(), mat_u.data(), n, mat_v.data(), n,
                work_vec.data(), work_vec.size(), real_work_vec.data(), info);
  if (info == 0)
  {
    for (int i = 0; i < n; ++i)
      mat_diag(i, i) = rng();
    mat_diag += 1.0;
    // reuse mat_a_t for this intermediate resutls
    BLAS::gemm(notrans, notrans, n, n, n, 1.0, mat_u.data(), n, mat_diag.data(), n, 0.0, mat_a_t.data(), n);
    // then reuse mat_u
    BLAS::gemm(notrans, notrans, n, n, n, 1.0, mat_a_t.data(), n, mat_v.data(), n, 0.0, mat_u.data(), n);
  }
  simd::transpose(mat_u.data(), n, mat_u.cols(), mat_spd.data(), n, n);
}

/** Functor to provide scope for rng when making SpdMatrix for testing.
 */
template<typename T>
class MakeRngSpdMatrix
{
public:
  void operator()(Matrix<T>& mat_spd) { makeRngSpdMatrix(rng, mat_spd); }

private:
  testing::RandomForTest<RngValueType<T>> rng;
};

/** make a random Vector
 */
template<typename T, typename = typename std::enable_if<std::is_floating_point<T>::value, void>::type>
void makeRngVector(testing::RandomForTest<RngValueType<T>>& rng, Vector<T>& vec)
{
  int n = vec.size();
  rng.fillBufferRng(vec.data(), vec.size());
}

/** Functor to provide scope for rng when making SpdMatrix for testing.
 */
template<typename T>
class MakeRngVector
{
public:
  void operator()(Vector<T>& vec) { makeRngSpdMatrix(rng, vec); }

private:
  testing::RandomForTest<RngValueType<T>> rng;
};

extern template class MakeRngSpdMatrix<double>;
extern template class MakeRngSpdMatrix<float>;
extern template class MakeRngSpdMatrix<std::complex<double>>;
extern template class MakeRngSpdMatrix<std::complex<float>>;

} // namespace testing
} // namespace qmcplusplus
