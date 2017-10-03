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
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Numerics/OhmmsBlas.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "simd/simd.hpp"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

template <typename T1, typename T2> void check_matrix(Matrix<T1> &a, Matrix<T2> &b)
{
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.rows(); i++) {
    for (int j = 0; j < a.cols(); j++) {
      REQUIRE(a(i,j) == ValueApprox(b(i,j)));
    }
  }
}

TEST_CASE("DiracMatrix_identity", "[wavefunction][fermion]")
{
  DiracMatrix<double> dm;

  Matrix<double> m;
  m.resize(3,3);
  m(0,0) = 1.0;
  m(1,1) = 1.0;
  m(2,2) = 1.0;

  dm.invert(m, true);
  REQUIRE(dm.LogDet == Approx(0.0));

  Matrix<double> eye;
  eye.resize(3,3);
  eye(0,0) = 1.0;
  eye(1,1) = 1.0;
  eye(2,2) = 1.0;

  check_matrix(m, eye);

}

TEST_CASE("DiracMatrix_inverse", "[wavefunction][fermion]")
{
  DiracMatrix<double> dm;

  Matrix<double> a;
  a.resize(3,3);

  a(0,0) = 2.3;
  a(0,1) = 4.5;
  a(0,2) = 2.6;
  a(1,0) = 0.5;
  a(1,1) = 8.5;
  a(1,2) = 3.3;
  a(2,0) = 1.8;
  a(2,1) = 4.4;
  a(2,2) = 4.9;

  dm.invert(a, true);
  REQUIRE(dm.LogDet == Approx(3.78518913425));

  Matrix<double> b;
  b.resize(3,3);

  b(0,0) = 0.6159749342;
  b(0,1) = -0.2408954682;
  b(0,2) = -0.1646081192;
  b(1,0) = 0.07923894288;
  b(1,1) = 0.1496231042;
  b(1,2) = -0.1428117337;
  b(2,0) = -0.2974298429;
  b(2,1) = -0.04586322768;
  b(2,2) = 0.3927890292;

  check_matrix(a, b);
}


TEST_CASE("DiracMatrix_update_row", "[wavefunction][fermion]")
{
  DiracMatrix<double> dm;

  Matrix<double> a;
  a.resize(3,3);

  a(0,0) = 2.3;
  a(0,1) = 4.5;
  a(0,2) = 2.6;
  a(1,0) = 0.5;
  a(1,1) = 8.5;
  a(1,2) = 3.3;
  a(2,0) = 1.8;
  a(2,1) = 4.4;
  a(2,2) = 4.9;

  dm.invert(a,false);

  // new row
  Vector<double> v;
  v.resize(3);
  v(0) = 1.9;
  v(1) = 2.0;
  v(2) = 3.1;

  double det_ratio1 = simd::dot(a[0],v.data(),3);

  double det_ratio = 0.178276269185;
  REQUIRE(det_ratio1 == Approx(det_ratio));
  dm.updateRow(a, v.data(), 0, det_ratio);

  Matrix<double> b;
  b.resize(3,3);

  b(0,0) =  3.455170657;
  b(0,1) =  -1.35124809;
  b(0,2) = -0.9233316353;
  b(1,0) = 0.05476311768;
  b(1,1) = 0.1591951095;
  b(1,2) = -0.1362710138;
  b(2,0) = -2.235099338;
  b(2,1) = 0.7119205298;
  b(2,2) = 0.9105960265;

  check_matrix(a, b);
}


}
