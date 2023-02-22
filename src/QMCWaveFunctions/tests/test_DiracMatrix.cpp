//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/Fermion/DelayedUpdate.h"
#include "CPU/SIMD/simd.hpp"
#include "checkMatrix.hpp"
#include "createTestMatrix.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;
using LogComplexApprox = Catch::Detail::LogComplexApprox;
  
TEST_CASE("DiracMatrix_identity", "[wavefunction][fermion]")
{
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> m, m_invT;
  LogValueType log_value;
  m.resize(3, 3);
  m_invT.resize(3, 3);

  fillIdentityMatrix(m);

  dm.invert_transpose(m, m_invT, log_value);
  CHECK(log_value == LogComplexApprox(0.0));

  Matrix<ValueType> eye;
  eye.resize(3, 3);
  fillIdentityMatrix(eye);

  CheckMatrixResult check_matrix_result = checkMatrix(m_invT, eye);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("DiracMatrix_inverse", "[wavefunction][fermion]")
{
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a, a_T, a_inv;
  LogValueType log_value;
  a.resize(3, 3);
  a_T.resize(3, 3);
  a_inv.resize(3, 3);

  TestMatrix1::fillInput(a);

  simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
  dm.invert_transpose(a_T, a_inv, log_value);
  CHECK(log_value == LogComplexApprox(TestMatrix1::logDet()));

  Matrix<ValueType> b;
  b.resize(3, 3);

  TestMatrix1::fillInverse(b);

  auto check_matrix_result = checkMatrix(a_inv, b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("DiracMatrix_inverse_matching", "[wavefunction][fermion]")
{
  DiracMatrix<double> dm;

  Matrix<double> a, a_T, a_inv;
  LogValueType log_value;
  a.resize(4, 4);
  a_T.resize(4, 4);
  a_inv.resize(4, 4);

  a(0, 0) = 6;
  a(0, 1) = 5;
  a(0, 2) = 7;
  a(0, 3) = 5;
  a(1, 0) = 2;
  a(1, 1) = 2;
  a(1, 2) = 5;
  a(1, 3) = 4;
  a(2, 0) = 8;
  a(2, 1) = 2;
  a(2, 2) = 6;
  a(2, 3) = 4;
  a(3, 0) = 3;
  a(3, 1) = 8;
  a(3, 2) = 6;
  a(3, 3) = 8;

  simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
  dm.invert_transpose(a_T, a_inv, log_value);
  CHECK(log_value == LogComplexApprox(std::complex<double>{5.50533, 6.28319}));

  Matrix<std::complex<double>> inv_M;
  inv_M.resize(4, 4);
  inv_M(0, 0) = -0.0650406504065041;
  inv_M(0, 1) = -0.2113821138211382;
  inv_M(0, 2) = 0.2113821138211382;
  inv_M(0, 3) = 0.04065040650406502;
  inv_M(1, 0) = 0.3739837398373983;
  inv_M(1, 1) = -0.28455284552845533;
  inv_M(1, 2) = -0.21544715447154467;
  inv_M(1, 3) = 0.016260162601626094;
  inv_M(2, 0) = 0.3902439024390243;
  inv_M(2, 1) = 0.2682926829268292;
  inv_M(2, 2) = -0.2682926829268292;
  inv_M(2, 3) = -0.24390243902439013;
  inv_M(3, 0) = -0.6422764227642275;
  inv_M(3, 1) = 0.1626016260162603;
  inv_M(3, 2) = 0.33739837398373973;
  inv_M(3, 3) = 0.2764227642276421;

  auto check_matrix_result = checkMatrix(inv_M, a_inv);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("DiracMatrix_inverse_matching_2", "[wavefunction][fermion]")
{
  DiracMatrix<double> dm;

  Matrix<double> a, a_T, a_inv;
  LogValueType log_value;
  a.resize(4, 4);
  a_T.resize(4, 4);
  a_inv.resize(4, 4);

  a(0, 0) = 6;
  a(0, 1) = 2;
  a(0, 2) = 8;
  a(0, 3) = 3;
  a(1, 0) = 5;
  a(1, 1) = 2;
  a(1, 2) = 2;
  a(1, 3) = 8;
  a(2, 0) = 7;
  a(2, 1) = 5;
  a(2, 2) = 6;
  a(2, 3) = 6;
  a(3, 0) = 5;
  a(3, 1) = 4;
  a(3, 2) = 4;
  a(3, 3) = 8;

  // lets check Xgetrf against the cuBLAS batched

  simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
  int pivot[4]{-1, -1, -1, -1};
  int status = Xgetrf(a_T.rows(), a_T.cols(), a_T.data(), a_T.cols(), pivot);
  std::vector<double> lu{7.0,
                         0.8571428571428571,
                         0.7142857142857142,
                         0.7142857142857142,
                         5.0,
                         -2.2857142857142856,
                         0.6874999999999998,
                         -0.18750000000000022,
                         6.0,
                         2.8571428571428577,
                         -4.249999999999999,
                         -0.05882352941176502,
                         6.0,
                         -2.1428571428571423,
                         5.1875,
                         3.617647058823531};

  Matrix<double> lu_mat(lu.data(), 4, 4);
  auto check_matrix_result = checkMatrix(lu_mat, a_T);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  dm.invert_transpose(a, a_inv, log_value);
  CHECK(log_value == LogComplexApprox(std::complex<double>{5.50533, 6.28319}));

  Matrix<std::complex<double>> inv_M;
  inv_M.resize(4, 4);
  inv_M(0, 0) = -0.0650406504065041;
  inv_M(0, 1) = -0.2113821138211382;
  inv_M(0, 2) = 0.2113821138211382;
  inv_M(0, 3) = 0.04065040650406502;
  inv_M(1, 0) = 0.3739837398373983;

  inv_M(1, 1) = -0.28455284552845533;
  inv_M(1, 2) = -0.21544715447154467;
  inv_M(1, 3) = 0.016260162601626094;
  inv_M(2, 0) = 0.3902439024390243;
  inv_M(2, 1) = 0.2682926829268292;
  inv_M(2, 2) = -0.2682926829268292;
  inv_M(2, 3) = -0.24390243902439013;
  inv_M(3, 0) = -0.6422764227642275;
  inv_M(3, 1) = 0.1626016260162603;
  inv_M(3, 2) = 0.33739837398373973;
  inv_M(3, 3) = 0.2764227642276421;

  check_matrix_result = checkMatrix(inv_M, a_inv);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

/** This test case is meant to match the cuBLAS_LU::getrf_batched_complex
 *
 */
TEST_CASE("DiracMatrix_inverse_complex", "[wavefunction][fermion]")
{
  DiracMatrix<std::complex<double>> dm;

  Matrix<std::complex<double>> a, a_T, a_inv;
  LogValueType log_value;
  a.resize(4, 4);
  a_T.resize(4, 4);
  a_inv.resize(4, 4);

  a(0, 0) = {2.0, 0.1};
  a(0, 1) = {5.0, 0.1};
  a(0, 2) = {7.0, 0.2};
  a(0, 3) = {5.0, 0.0};
  a(1, 0) = {5.0, 0.1};
  a(1, 1) = {2.0, 0.2};
  a(1, 2) = {5.0, 1.0};
  a(1, 3) = {4.0, -0.1};
  a(2, 0) = {8.0, 0.5};
  a(2, 1) = {2.0, 0.1};
  a(2, 2) = {6.0, -0.2};
  a(2, 3) = {4.0, -0.6};
  a(3, 0) = {7.0, 1.0};
  a(3, 1) = {8.0, 0.5};
  a(3, 2) = {6.0, -0.2};
  a(3, 3) = {8.0, -2.0};

  simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
  int pivot[4]{-1, -1, -1, -1};
  int status = Xgetrf(a_T.rows(), a_T.cols(), a_T.data(), a_T.cols(), pivot);

  double lu[32]{8.0,
                0.5,
                0.8793774319066148,
                0.07003891050583658,
                0.24980544747081712,
                -0.0031128404669260694,
                0.6233463035019455,
                -0.026459143968871595,
                2.0,
                0.1,
                6.248249027237354,
                0.2719844357976654,
                0.7194170575332381,
                -0.01831314754114669,
                0.1212375092639108,
                0.02522449751055713,
                6.0,
                -0.2,
                0.7097276264591441,
                -0.4443579766536965,
                4.999337315778741,
                0.6013141870887196,
                0.26158183940834034,
                0.23245112532996867,
                4.0,
                -0.6,
                4.440466926070039,
                -1.7525291828793774,
                0.840192589866152,
                1.5044529443071093,
                1.0698651110730424,
                -0.10853319738453365};

  Matrix<std::complex<double>> lu_mat(reinterpret_cast<std::complex<double>*>(lu), 4, 4);
  auto check_matrix_result = checkMatrix(lu_mat, a_T);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  dm.invert_transpose(a, a_inv, log_value);
  CHECK(log_value == LogComplexApprox(std::complex<double>{5.603777579195571, 0.12452497406076501}));

  Matrix<std::complex<double>> b;
  b.resize(4, 4);

  b(0, 0)             = {-0.05356228836958328, 0.018778132944668208};
  b(0, 1)             = {0.20217332569060176, 0.10247607676441389};
  b(0, 2)             = {0.1752053760997507, 0.08545803475354152};
  b(0, 3)             = {-0.22019253129809616, -0.23960901438764967};
  b(1, 0)             = {-0.16917709116094917, 0.18307841761769691};
  b(1, 1)             = {-0.6356039515926515, 0.24028909525203923};
  b(1, 2)             = {-0.16030725389326506, -0.2749054789157097};
  b(1, 3)             = {0.9251760746083686, 0.09385511919367814};
  b(2, 0)             = {0.21572431303125872, -0.10633509516999905};
  b(2, 1)             = {0.24869726332502523, -0.11905352270298367};
  b(2, 2)             = {0.16292868571183988, 0.17275739697798137};
  b(2, 3)             = {-0.560684625012445, -0.09607837395378985};
  b(3, 0)             = {0.023435592783503056, -0.030184486280558254};
  b(3, 1)             = {0.09553876547173154, -0.09007339752988484};
  b(3, 2)             = {-0.11510984212629605, -0.020898880528951114};
  b(3, 3)             = {0.05299966349431483, 0.13363053130258684};
  check_matrix_result = checkMatrix(b, a_inv);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}


TEST_CASE("DiracMatrix_update_row", "[wavefunction][fermion]")
{
  DiracMatrix<ValueType> dm;
  DelayedUpdate<ValueType, QMCTraits::QTFull::ValueType> updateEng;
  updateEng.resize(3, 1);

  Matrix<ValueType> a, a_T, a_inv;
  LogValueType log_value;
  a.resize(3, 3);
  a_T.resize(3, 3);
  a_inv.resize(3, 3);

  a(0, 0) = 2.3;
  a(0, 1) = 4.5;
  a(0, 2) = 2.6;
  a(1, 0) = 0.5;
  a(1, 1) = 8.5;
  a(1, 2) = 3.3;
  a(2, 0) = 1.8;
  a(2, 1) = 4.4;
  a(2, 2) = 4.9;

  simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
  dm.invert_transpose(a_T, a_inv, log_value);

  // new row
  Vector<ValueType> v(3), invRow(3);
  v[0] = 1.9;
  v[1] = 2.0;
  v[2] = 3.1;

  updateEng.getInvRow(a_inv, 0, invRow);
  ValueType det_ratio1 = simd::dot(v.data(), invRow.data(), invRow.size());

  ValueType det_ratio = 0.178276269185;
  CHECK(det_ratio1 == ValueApprox(det_ratio));
  updateEng.acceptRow(a_inv, 0, v, det_ratio1);

  Matrix<ValueType> b;
  b.resize(3, 3);

  b(0, 0) = 3.455170657;
  b(0, 1) = -1.35124809;
  b(0, 2) = -0.9233316353;
  b(1, 0) = 0.05476311768;
  b(1, 1) = 0.1591951095;
  b(1, 2) = -0.1362710138;
  b(2, 0) = -2.235099338;
  b(2, 1) = 0.7119205298;
  b(2, 2) = 0.9105960265;

  auto check_matrix_result = checkMatrix(a_inv, b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

} // namespace qmcplusplus
