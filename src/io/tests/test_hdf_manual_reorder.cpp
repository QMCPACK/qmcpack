
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


/*
TEST_CASE("hdf_write_reorder", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_write_reorder.hdf");

  vector<double> v{1.1, -2.2, 3.3, -4.4, 5.5, -6.6, 7.7, -8.8, 9.9};
  vector<hsize_t> dims{3,3};
  bool okay = hd.writeEntry(v, dims, "array_from_vector");
  REQUIRE(okay);  
  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_write_reorder.hdf");
  REQUIRE(okay);

  Tensor<double,3> read_tensor;
  hd2.read(read_tensor, "array_from_vector");
  REQUIRE(read_tensor.size() == v.size());
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      int pos = i*3+j;
      REQUIRE(v[pos] == read_tensor(i,j));
    }
  }
  hd2.close();
}
*/

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


// do this in three ways.
// 1. relying on my overloading in std::vector (do not perfer this, too specific)
// 2. rely on my understanding of hyperslab_proxy (seems to require you to know
//          too much about the dimensions in the dataset and is prone to user error)
// 3. use hyperslab_proxy to handle the back end api, but specify the slice to choose
//          with a vector and query the dataset to figure out the other needed
//          parts.  better than #1 because it does not specially call out vector
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

  // method 1 (relying on vector in hdf_stl)
  vector<int> readSpec{1, -1};
  vector<double> readBuffer(4);

  hd2.read(readBuffer, readSpec, "matrix");
  for (int i = 0; i < 4; i++)
  {
    REQUIRE(readBuffer[i] == Approx(allData(1, i)));
  }

  readSpec[0] = -1;
  readSpec[1] = 2;
  hd2.read(readBuffer, readSpec, "matrix");
  for (int i = 0; i < 3; i++)
  {
    REQUIRE(readBuffer[i] == Approx(allData(i, 2)));
  }

  readSpec[0] = 2;
  readSpec[1] = 0;
  hd2.read(readBuffer, readSpec, "matrix");
  REQUIRE(readBuffer[0] == Approx(allData(2, 0)));

  // method 2 (direct utilization of hyperslab_proxy in client code)
  Matrix<double> outbuffer1(1, 4);
  Matrix<double> outbuffer2(3, 1);
  Matrix<double> outbuffer3(1, 1);

  std::array<int, 2> dims_unused;
  dims_unused[0] = 1;
  dims_unused[1] = 2;
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
  {
    REQUIRE(outbuffer2(i, 0) == Approx(allData(i, 2)));
  }

  dims_local[0] = 1;
  dims_local[1] = 1;
  offsets[0]    = 2;
  offsets[1]    = 0;
  hyperslab_proxy<Matrix<double>, 2> pxy3(outbuffer3, dims_unused, dims_local, offsets);
  hd2.read(pxy3, "matrix");
  REQUIRE(outbuffer3(0, 0) == Approx(allData(2, 0)));

  // method 3 here
  //  buftype sto(4);
/*
  Matrix<double> locob1(1, 4);
  Matrix<double> locob2(3, 1);
  Matrix<double> locob3(1, 1);
  TinyVector<int, 2> locReadSpec{1, -1};
  hd2.read2(locob1, locReadSpec, "matrix");
  for (int i = 0; i < 4; i++)
  {
    REQUIRE(locob1(0, i) == Approx(allData(1, i)));
  }

  locReadSpec[0] = -1;
  locReadSpec[1] = 2;
  hd2.read2(locob2, locReadSpec, "matrix");
  for (int i = 0; i < 3; i++)
  {
    REQUIRE(locob2.data()[i] == Approx(allData(i, 2)));
    REQUIRE(locob2(i) == Approx(allData(i,2)));
  }

  locReadSpec[0] = 2;
  locReadSpec[1] = 0;
  hd2.read2(locob3, locReadSpec, "matrix");
  REQUIRE(locob3.data()[0] == Approx(allData(2, 0)));
*/
  hd2.close();
}
