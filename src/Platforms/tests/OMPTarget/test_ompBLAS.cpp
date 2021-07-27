//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <memory>
#include <vector>
#include <iostream>
#include "OMPTarget/OMPallocator.hpp"
#include "OMPTarget/ompBLAS.hpp"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <CPU/BLAS.hpp>
#include <checkMatrix.hpp>

namespace qmcplusplus
{

template<typename T>
void test_gemvT(const int N)
{
  using vec_t = Vector<T, OMPallocator<T>>;
  using mat_t = Matrix<T, OMPallocator<T>>;

  ompBLAS::ompBLAS_handle handle;

  vec_t A(N);    // Input vector
  mat_t B(N, N); // Input matrix
  vec_t C(N);    // Result vector (test)
  vec_t D(N);    // Result vector BLAS

  for (int dim = 1; dim <= N; dim++)
  {
    //std::cout << "checking dim " << dim << std::endl;
    handle = dim;
    for (int i = 0; i < dim; i++)
    {
      A[i] = i;
      for (int j = 0; j < dim; j++)
        B.data()[i * dim + j] = i + j;
    }

    A.updateTo();
    B.updateTo();

    // tests omp gemv
    ompBLAS::gemv(handle, 'T', dim, dim, 1.0, B.device_data(), dim, A.device_data(), 1, 0.0, C.device_data(), 1);
    C.updateFrom();

    BLAS::gemv(dim, dim, B.data(), A.data(), D.data());

    int index = 0;
    do
    {
      CHECK(C[index] == D[index]);
      index++;
    } while (index < dim);
  }
}

TEST_CASE("OmpBLAS gemv", "[OMP]")
{
  const int N = 100;
  test_gemvT<float>(N);
  test_gemvT<double>(N);
#if defined(QMC_COMPLEX)
  test_gemvT<std::complex<float>>(N);
  test_gemvT<std::complex<double>>(N);
#endif
}
} // namespace qmcplusplus
