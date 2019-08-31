
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



TEST_CASE("hdf_write_reorder_with_matrix", "[hdf]")
{
  hdf_archive hd;
  bool okay = hd.create("test_write_matrix_reorder.hdf");
  REQUIRE(okay);

  Vector<double> v(6);
  v[0] = 0.0;
  v[1] = 0.1;
  v[2] = 0.2;
  v[3] = 1.0;
  v[4] = 1.1;
  v[5] = 1.2;

  Matrix<double> matrix_view(v.data(), 2, 3);
  hd.write(matrix_view, "matrix_from_vector");

  hdf_archive hd2;
  hd2.open("test_write_matrix_reorder.hdf");

  Matrix<double> m(2, 3);
  hd2.read(m, "matrix_from_vector");

  REQUIRE(m(0, 0) == Approx(v[0]));
  REQUIRE(m(0, 1) == Approx(v[1]));
  REQUIRE(m(0, 2) == Approx(v[2]));
  REQUIRE(m(1, 0) == Approx(v[3]));
  REQUIRE(m(1, 1) == Approx(v[4]));
  REQUIRE(m(1, 2) == Approx(v[5]));
}


// do this in two different ways.
// 1. directly create a hyperslab_proxy
// 2. use readHyperslab which will create the hyperslab_proxy for you
TEST_CASE("hdf_read_partial", "[hdf]")
{
  hdf_archive hd;
  bool okay = hd.create("test_read_partial.hdf");
  REQUIRE(okay);

  Matrix<double> allData(3, 4);
  Matrix<std::complex<float>> allData_cplx(3, 4);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      allData(i, j) = i + j * 0.1;
      allData_cplx(i, j) = std::complex<float>(i,j * 0.1);
    }
  }

  hd.write(allData, "matrix");
  hd.write(allData_cplx, "matrix_cplx_float");
  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_read_partial.hdf");
  REQUIRE(okay);

  // test getShape
  std::vector<int> datashape;
  okay = hd2.getShape<double>("matrix", datashape);
  REQUIRE(okay);
  REQUIRE(datashape.size() == 2);
  REQUIRE(datashape[0] == 3);
  REQUIRE(datashape[1] == 4);

  okay = hd2.getShape<std::complex<float>>("matrix_cplx_float", datashape);
  REQUIRE(okay);
  REQUIRE(datashape.size() == 2);
  REQUIRE(datashape[0] == 3);
  REQUIRE(datashape[1] == 4);

  //treat std::complex<float> as an array[2]
  okay = hd2.getShape<float>("matrix_cplx_float", datashape);
  REQUIRE(okay);
  REQUIRE(datashape.size() == 3);
  REQUIRE(datashape[0] == 3);
  REQUIRE(datashape[1] == 4);
  REQUIRE(datashape[2] == 2);

  // method 1 (direct utilization of hyperslab_proxy in client code)
  Matrix<double> outbuffer1(1, 4);
  Matrix<double> outbuffer2(3, 1);
  Matrix<double> outbuffer3(1, 1);
  Matrix<double> outbuffer4;

  std::array<int, 2> dims_unused;
  dims_unused[0] = 3;
  dims_unused[1] = 4;
  std::array<int, 2> dims_local;
  dims_local[0] = 1;
  dims_local[1] = 4;
  std::array<int, 2> offsets;
  offsets[0] = 1;
  offsets[1] = 0;

  hyperslab_proxy<Matrix<double>, 2> pxy1(outbuffer1, dims_unused, dims_local, offsets);
  hd2.read(pxy1, "matrix");
  for (int i = 0; i < 4; i++)
  {
    REQUIRE(outbuffer1(0, i) == Approx(allData(1, i)));
  }

  dims_local[0] = 3;
  dims_local[1] = 1;
  offsets[0]    = 0;
  offsets[1]    = 2;
  hyperslab_proxy<Matrix<double>, 2> pxy2(outbuffer2, dims_unused, dims_local, offsets);
  hd2.read(pxy2, "matrix");
  for (int i = 0; i < 3; i++)
    REQUIRE(outbuffer2(i, 0) == Approx(allData(i, 2)));

  // set dims_unused to sero and to be detect
  dims_unused[0] = 0;
  dims_unused[1] = 0;

  dims_local[0] = 1;
  dims_local[1] = 1;
  offsets[0]    = 2;
  offsets[1]    = 0;
  hyperslab_proxy<Matrix<double>, 2> pxy3(outbuffer3, dims_unused, dims_local, offsets);
  hd2.read(pxy3, "matrix");
  REQUIRE(outbuffer3(0, 0) == Approx(allData(2, 0)));

  // mostly the same as outbuffer2 but outbuffer4 resized by hyperslab_proxy
  dims_local[0] = 3;
  dims_local[1] = 1;
  offsets[0]    = 0;
  offsets[1]    = 2;
  hyperslab_proxy<Matrix<double>, 2> pxy4(outbuffer4, dims_unused, dims_local, offsets);
  hd2.read(pxy4, "matrix");
  REQUIRE(outbuffer4.rows() == 3);
  REQUIRE(outbuffer4.cols() == 1);
  for (int i = 0; i < 3; i++)
    REQUIRE(outbuffer2(i, 0) == Approx(allData(i, 2)));

  // method 2 here
  Matrix<double> locob1(1, 4);
  Matrix<double> locob2(3, 1);
  Matrix<double> locob3(1, 1);
  std::array<int, 2> readSpec{1, -1};
  hd2.readHyperslab(locob1, readSpec, "matrix");
  for (int i = 0; i < 4; i++)
  {
    REQUIRE(locob1(0, i) == Approx(allData(1, i)));
  }

  readSpec[0] = -1;
  readSpec[1] = 2;
  hd2.readHyperslab(locob2, readSpec, "matrix");
  for (int i = 0; i < 3; i++)
  {
    REQUIRE(locob2.data()[i] == Approx(allData(i, 2)));
    REQUIRE(locob2(i) == Approx(allData(i,2)));
  }

  readSpec[0] = 2;
  readSpec[1] = 0;
  hd2.readHyperslab(locob3, readSpec, "matrix");
  REQUIRE(locob3.data()[0] == Approx(allData(2, 0)));
  hd2.close();
}
