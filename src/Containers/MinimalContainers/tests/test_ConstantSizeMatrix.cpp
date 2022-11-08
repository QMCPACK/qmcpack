//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "MinimalContainers/ConstantSizeMatrix.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"

#include "catch.hpp"


namespace qmcplusplus
{
TEST_CASE("ConstantSizeMatrix basics", "[containers]")
{
  ConstantSizeMatrix<double> cmat(1, 9, 1, 32, 0.0);
  CHECK(cmat.size() == 9);
  CHECK(cmat.capacity() == 32);

  CHECK_NOTHROW(cmat.resize(1, 16));
  CHECK(cmat.size() == 16);
  CHECK(cmat.capacity() == 32);

  CHECK_THROWS(cmat.resize(1, 33));
  CHECK_THROWS(cmat.resize(34, 2));

  CHECK_NOTHROW(cmat.resize(2, 9));
  CHECK(cmat.capacity() == 32);

  CHECK_NOTHROW(cmat.resize(3, 9));
  CHECK(cmat.capacity() == 32);

  //note the odd OhmmsMatrix style access operator semantics
  //If these break you probably breaking other weird QMCPACK code.
  std::vector<double> svec(16, 1.0);
  cmat               = svec;
  double& matrices_d = cmat(0);
  CHECK(matrices_d == Approx(1.0));
  matrices_d = 4.0;
  CHECK(cmat(0) == Approx(4.0));
  CHECK(*(cmat[0]) == Approx(4.0));

  ConstantSizeMatrix<double> cmat2(2, 8, 2, 16, 0.0);
  std::vector<double> svec2(16, 2.0);
  svec.insert(svec.end(), svec2.begin(), svec2.end());
  cmat2 = svec;
  CHECK(*(cmat2[0]) == Approx(1.0));
  CHECK(*(cmat2[0] + 1) == Approx(1.0));
  CHECK(*(cmat2[1]) == Approx(2.0));

  ConstantSizeMatrix<double> cmat3(4, 32, 4, 32, 0.0);
  CHECK_NOTHROW(cmat3.resize(2, 64));
  CHECK(cmat3.size() == 128);
  CHECK(cmat3.capacity() == 128);

  CHECK_NOTHROW(cmat3.resize(8, 16));
  CHECK(cmat3.size() == 128);
  CHECK(cmat3.capacity() == 128);
}

TEST_CASE("ConstantSizeMatrix Ohmms integration", "[containers]")
{
  ConstantSizeMatrix<double> cmat(2, 8, 2, 32, 0.0);
  Matrix<double> omat(2, 9);
  omat = 2.0;
  cmat = omat;

  CHECK(cmat.size() == 18);
  CHECK(*omat[1] == 2.0);
}

} // namespace qmcplusplus
