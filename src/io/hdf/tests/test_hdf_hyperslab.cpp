
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

// do this in two different ways.
// 1. directly create a hyperslab_proxy
// 2. use readSlabSelection which will create the hyperslab_proxy for you
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
      allData(i, j)      = i + j * 0.1;
      allData_cplx(i, j) = std::complex<float>(i, j * 0.1);
    }
  }

  hd.write(allData, "matrix");

  const auto& const_allData_cplx = allData_cplx;
  hd.write(const_allData_cplx, "matrix_cplx_float");
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

  std::array<size_t, 2> dims_unused;
  dims_unused[0] = 3;
  dims_unused[1] = 4;
  std::array<size_t, 2> dims_local;
  dims_local[0] = 1;
  dims_local[1] = 4;
  std::array<size_t, 2> offsets;
  offsets[0] = 1;
  offsets[1] = 0;

  hyperslab_proxy<Matrix<double>, 2> pxy1(outbuffer1, dims_unused, dims_local, offsets);
  hd2.read(pxy1, "matrix");
  for (int i = 0; i < 4; i++)
  {
    CHECK(outbuffer1(0, i) == Approx(allData(1, i)));
  }

  dims_local[0] = 3;
  dims_local[1] = 1;
  offsets[0]    = 0;
  offsets[1]    = 2;
  hyperslab_proxy<Matrix<double>, 2> pxy2(outbuffer2, dims_unused, dims_local, offsets);
  hd2.read(pxy2, "matrix");
  for (int i = 0; i < 3; i++)
    CHECK(outbuffer2(i, 0) == Approx(allData(i, 2)));

  // set dims_unused to sero and to be detect
  dims_unused[0] = 0;
  dims_unused[1] = 0;

  dims_local[0] = 1;
  dims_local[1] = 1;
  offsets[0]    = 2;
  offsets[1]    = 0;
  hyperslab_proxy<Matrix<double>, 2> pxy3(outbuffer3, dims_unused, dims_local, offsets);
  hd2.read(pxy3, "matrix");
  CHECK(outbuffer3(0, 0) == Approx(allData(2, 0)));

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
    CHECK(outbuffer2(i, 0) == Approx(allData(i, 2)));

  // method 2 here
  std::vector<double> locob1;
  Matrix<double> locob2(3, 1);
  Matrix<double> locob3(1, 1);
  std::array<int, 2> readSpec{1, -1};
  hd2.readSlabSelection(locob1, readSpec, "matrix");
  REQUIRE(locob1.size() == 4);
  for (int i = 0; i < 4; i++)
  {
    CHECK(locob1[i] == Approx(allData(1, i)));
  }

  readSpec[0] = -1;
  readSpec[1] = 2;
  hd2.readSlabSelection(locob2, readSpec, "matrix");
  for (int i = 0; i < 3; i++)
  {
    CHECK(locob2.data()[i] == Approx(allData(i, 2)));
    CHECK(locob2(i) == Approx(allData(i, 2)));
  }

  readSpec[0] = 2;
  readSpec[1] = 0;
  hd2.readSlabSelection(locob3, readSpec, "matrix");
  CHECK(locob3.data()[0] == Approx(allData(2, 0)));
  hd2.close();
}
