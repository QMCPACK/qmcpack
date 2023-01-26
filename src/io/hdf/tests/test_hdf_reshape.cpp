
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
#include <iostream>
#include <vector>
#include <complex>
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "type_traits/container_traits_ohmms.h"
#include "hdf/hdf_archive.h"

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

  std::array<int, 2> shape{2, 3};
  hd.writeSlabReshaped(v, shape, "matrix_from_vector");

  std::vector<std::complex<float>> v_cplx(6);
  v_cplx[0] = {0.0, 0.2};
  v_cplx[1] = {0.1, 0.3};
  v_cplx[2] = {0.2, 0.4};
  v_cplx[3] = {1.0, 1.2};
  v_cplx[4] = {1.1, 1.3};
  v_cplx[5] = {1.2, 1.4};
  hd.writeSlabReshaped(v_cplx, shape, "matrix_from_vector_cplx");

  hd.close();

  hdf_archive hd2;
  hd2.open("test_write_matrix_reshape.hdf");

  Matrix<double> m; // resized by read
  hd2.read(m, "matrix_from_vector");

  REQUIRE(m.rows() == 2);
  REQUIRE(m.cols() == 3);
  CHECK(m(0, 0) == Approx(v[0]));
  CHECK(m(0, 1) == Approx(v[1]));
  CHECK(m(0, 2) == Approx(v[2]));
  CHECK(m(1, 0) == Approx(v[3]));
  CHECK(m(1, 1) == Approx(v[4]));
  CHECK(m(1, 2) == Approx(v[5]));

  std::vector<double> vec; // reshaped and resized by read
  hd2.readSlabReshaped(vec, shape, "matrix_from_vector");

  REQUIRE(vec.size() == 6);
  CHECK(vec[0] == Approx(v[0]));
  CHECK(vec[1] == Approx(v[1]));
  CHECK(vec[2] == Approx(v[2]));
  CHECK(vec[3] == Approx(v[3]));
  CHECK(vec[4] == Approx(v[4]));
  CHECK(vec[5] == Approx(v[5]));

  // using hyperslab selection
  Vector<std::complex<float>> vec_cplx; // reshaped and resized by read
  std::array<int, 2> spec{-1, -1};
  hd2.readSlabSelection(vec_cplx, spec, "matrix_from_vector_cplx");

  REQUIRE(vec_cplx.size() == 6);
  CHECK(vec_cplx[0] == ComplexApprox(v_cplx[0]));
  CHECK(vec_cplx[1] == ComplexApprox(v_cplx[1]));
  CHECK(vec_cplx[2] == ComplexApprox(v_cplx[2]));
  CHECK(vec_cplx[3] == ComplexApprox(v_cplx[3]));
  CHECK(vec_cplx[4] == ComplexApprox(v_cplx[4]));
  CHECK(vec_cplx[5] == ComplexApprox(v_cplx[5]));
}
