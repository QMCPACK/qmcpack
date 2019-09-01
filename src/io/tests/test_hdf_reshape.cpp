
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "io/hdf_archive.h"
#include <iostream>
#include <vector>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>

using std::cout;
using std::endl;
using std::vector;

using namespace qmcplusplus;



TEST_CASE("hdf_write_reshape_with_matrix", "[hdf]")
{
  hdf_archive hd;
  bool okay = hd.create("test_write_matrix_reshape.hdf");
  REQUIRE(okay);

  Vector<double> v(6);
  v[0] = 0.0;
  v[1] = 0.1;
  v[2] = 0.2;
  v[3] = 1.0;
  v[4] = 1.1;
  v[5] = 1.2;

  std::array<size_t, 2> shape{2,3};
  hd.writeSlabReshaped(v, shape, "matrix_from_vector");

  hdf_archive hd2;
  hd2.open("test_write_matrix_reshape.hdf");

  Matrix<double> m(2, 3);
  hd2.read(m, "matrix_from_vector");

  REQUIRE(m(0, 0) == Approx(v[0]));
  REQUIRE(m(0, 1) == Approx(v[1]));
  REQUIRE(m(0, 2) == Approx(v[2]));
  REQUIRE(m(1, 0) == Approx(v[3]));
  REQUIRE(m(1, 1) == Approx(v[4]));
  REQUIRE(m(1, 2) == Approx(v[5]));
}
