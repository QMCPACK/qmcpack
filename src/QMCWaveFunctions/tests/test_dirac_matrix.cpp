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

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;

TEST_CASE("DiracMatrix_identity", "[wavefunction][fermion]")
{
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> m, m_invT;
  LogValueType LogValue;
  m.resize(3, 3);
  m_invT.resize(3, 3);

  m(0, 0) = 1.0;
  m(1, 1) = 1.0;
  m(2, 2) = 1.0;

  dm.invert_transpose(m, m_invT, LogValue);
  REQUIRE(LogValue == LogComplexApprox(0.0));

  Matrix<ValueType> eye;
  eye.resize(3, 3);
  eye(0, 0) = 1.0;
  eye(1, 1) = 1.0;
  eye(2, 2) = 1.0;

  auto check_matrix_result = checkMatrix(m_invT, eye);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("DiracMatrix_inverse", "[wavefunction][fermion]")
{
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a, a_T, a_inv;
  LogValueType LogValue;
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
  dm.invert_transpose(a_T, a_inv, LogValue);
  REQUIRE(LogValue == LogComplexApprox(3.78518913425));

  Matrix<ValueType> b;
  b.resize(3, 3);

  b(0, 0) = 0.6159749342;
  b(0, 1) = -0.2408954682;
  b(0, 2) = -0.1646081192;
  b(1, 0) = 0.07923894288;
  b(1, 1) = 0.1496231042;
  b(1, 2) = -0.1428117337;
  b(2, 0) = -0.2974298429;
  b(2, 1) = -0.04586322768;
  b(2, 2) = 0.3927890292;

  auto check_matrix_result = checkMatrix(a_inv, b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("DiracMatrix_inverse_complex", "[wavefunction][fermion]")
{
  DiracMatrix<std::complex<double>> dm;

  Matrix<std::complex<double>> a, a_T, a_inv;
  LogValueType LogValue;
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
  dm.invert_transpose(a_T, a_inv, LogValue);
  CHECK(LogValue == LogComplexApprox(std::complex<double>{5.603777579195571, 0.12452497406076501}));

  Matrix<std::complex<double>> b;
  b.resize(4, 4);

  b(0, 0) = {-0.05356228836958328, 0.018778132944668208};
  b(0, 1) = {-0.16917709116094917, 0.18307841761769691};
  b(0, 2) = {0.21572431303125872, -0.10633509516999905};
  b(0, 3) = {0.023435592783503056, -0.030184486280558254};
  b(1, 0) = {0.20217332569060176, 0.10247607676441389};
  b(1, 1) = {-0.6356039515926515, 0.24028909525203923};
  b(1, 2) = {0.24869726332502523, -0.11905352270298367};
  b(1, 3) = {0.09553876547173154, -0.09007339752988484};
  b(2, 0) = {0.1752053760997507, 0.08545803475354152};
  b(2, 1) = {-0.16030725389326506, -0.2749054789157097};
  b(2, 2) = {0.16292868571183988, 0.17275739697798137};
  b(2, 3) = {-0.11510984212629605, -0.020898880528951114};
  b(3, 0) = {-0.22019253129809616, -0.23960901438764967};
  b(3, 1) = {0.9251760746083686, 0.09385511919367814};
  b(3, 2) = {-0.560684625012445, -0.09607837395378985};
  b(3, 3) = {0.05299966349431483, 0.13363053130258684};

  auto check_matrix_result = checkMatrix(a_inv, b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}


TEST_CASE("DiracMatrix_update_row", "[wavefunction][fermion]")
{
  DiracMatrix<ValueType> dm;
  DelayedUpdate<ValueType, QMCTraits::QTFull::ValueType> updateEng;
  updateEng.resize(3, 1);

  Matrix<ValueType> a, a_T, a_inv;
  LogValueType LogValue;
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
  dm.invert_transpose(a_T, a_inv, LogValue);

  // new row
  Vector<ValueType> v(3), invRow(3);
  v[0] = 1.9;
  v[1] = 2.0;
  v[2] = 3.1;

  updateEng.getInvRow(a_inv, 0, invRow);
  ValueType det_ratio1 = simd::dot(v.data(), invRow.data(), invRow.size());

  ValueType det_ratio = 0.178276269185;
  REQUIRE(det_ratio1 == ValueApprox(det_ratio));
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
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message);}
}

} // namespace qmcplusplus
