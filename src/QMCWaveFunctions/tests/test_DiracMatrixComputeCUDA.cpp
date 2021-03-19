//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch.hpp>
#include <algorithm>
#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrixComputeCUDA.hpp"
#include "QMCWaveFunctions/tests/CheckMatrix.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"

namespace qmcplusplus
{
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
template<typename T>
using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;
template<typename T>
using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;


TEST_CASE("DiracMatrixComputeCUDA_different_batch_sizes", "[wavefunction][fermion]")
{
  DiracMatrixComputeCUDA<double> dmcc;
  OffloadPinnedMatrix<double> mat_a;
  mat_a.resize(4, 4);
  double A[16] = {2, 5, 8, 7, 5, 2, 2, 8, 7, 5, 6, 6, 5, 4, 4, 8};
  std::copy_n(A, 16, mat_a.data());
  OffloadPinnedVector<std::complex<double>> log_values;
  log_values.resize(1);
  OffloadPinnedMatrix<double> inv_mat_a;
  inv_mat_a.resize(4, 4);
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  dmcc.invert_transpose(*cuda_handles, mat_a, inv_mat_a, log_values);

  // std::cout << mat_a;
  // std::cout << inv_mat_a;
  // std::cout << log_values;

  OffloadPinnedMatrix<double> mat_b;
  mat_b.resize(4, 4);
  double invA[16]{-0.08247423, -0.26804124, 0.26804124, 0.05154639,  0.18556701,  -0.89690722, 0.39690722,  0.13402062,
                  0.24742268,  -0.19587629, 0.19587629, -0.15463918, -0.29896907, 1.27835052,  -0.77835052, 0.06185567};
  std::copy_n(invA, 16, mat_b.data());

  checkMatrix(inv_mat_a, mat_b, std::string(__FILE__), __LINE__);

  OffloadPinnedMatrix<double> mat_a2;
  mat_a2.resize(4, 4);
  std::copy_n(A, 16, mat_a2.data());
  OffloadPinnedMatrix<double> inv_mat_a2;
  inv_mat_a2.resize(4, 4);

  RefVector<const OffloadPinnedMatrix<double>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats{inv_mat_a, inv_mat_a2};

  log_values.resize(2);
  dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values);

  checkMatrix(inv_mat_a, mat_b, std::string(__FILE__), __LINE__);
  checkMatrix(inv_mat_a2, mat_b, std::string(__FILE__), __LINE__);
  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[1] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));

  OffloadPinnedMatrix<double> mat_a3;
  mat_a3.resize(4, 4);
  std::copy_n(A, 16, mat_a3.data());
  OffloadPinnedMatrix<double> inv_mat_a3;
  inv_mat_a3.resize(4, 4);

  a_mats[1] = mat_a3;
  
  RefVector<const OffloadPinnedMatrix<double>> a_mats3{mat_a, mat_a2, mat_a3};
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats3{inv_mat_a, inv_mat_a2, inv_mat_a3};

  log_values.resize(3);
  dmcc.mw_invertTranspose(*cuda_handles, a_mats3, inv_a_mats3, log_values);

  checkMatrix(inv_mat_a, mat_b, std::string(__FILE__), __LINE__);
  checkMatrix(inv_mat_a2, mat_b, std::string(__FILE__), __LINE__);
  checkMatrix(inv_mat_a3, mat_b, std::string(__FILE__), __LINE__);
  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[1] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[2] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));  
}
  
// TEST_CASE("DiracMatrixComputeCUDA_two_walkers", "[wavefunction][fermion]")
// {
//   DiracMatrixComputeCUDA<double> dmcc;
//   OffloadPinnedMatrix<double> mat_a;
//   OffloadPinnedMatrix<double> mat_a2;
//   mat_a.resize(4,4);
//   mat_a2.resize(4,4);

//   double A[16] = {2, 5, 8, 7,
//        5, 2, 2, 8,
//        7, 5, 6, 6,
//     5, 4, 4, 8};
//   std::copy_n(A,16,mat_a.data());
//   std::copy_n(A,16,mat_a2.data());
//   OffloadPinnedVector<std::complex<double>> log_values;
//   log_values.resize(2);
//   OffloadPinnedMatrix<double> inv_mat_a;
//   OffloadPinnedMatrix<double> inv_mat_a2;

//   inv_mat_a.resize(4,4);
//   inv_mat_a2.resize(4,4);

//   auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();

//   RefVector<OffloadPinnedMatrix<double>> a_mats{mat_a, mat_a2};
//   RefVector<OffloadPinnedMatrix<double>> inv_a_mats{inv_mat_a, inv_mat_a2};

//   dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values);

//   OffloadPinnedMatrix<double> mat_b;
//   mat_b.resize(4,4);
//   double invA[16] {-0.08247423, -0.26804124,  0.26804124,  0.05154639,
//         0.18556701, -0.89690722,  0.39690722,  0.13402062,
//         0.24742268, -0.19587629,  0.19587629, -0.15463918,
//     -0.29896907,  1.27835052, -0.77835052,  0.06185567};
//   std::copy_n(invA,16,mat_b.data());

//   checkMatrix(inv_mat_a, mat_b, std::string(__FILE__), __LINE__);
//   checkMatrix(inv_mat_a2, mat_b, std::string(__FILE__), __LINE__);
//   CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328,6.283185307179586}));
//   CHECK(log_values[1] == ComplexApprox(std::complex<double>{5.267858159063328,6.283185307179586}));
// }

} // namespace qmcplusplus
