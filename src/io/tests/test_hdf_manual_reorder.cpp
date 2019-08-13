
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

using std::cout;
using std::endl;
using std::vector;
using namespace qmcplusplus;


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

TEST_CASE("hdf_read_partial", "[hdf]")
{
  hdf_archive hd;
  bool okay = hd.create("test_read_partial.hdf");
  REQUIRE(okay);
  
  Tensor<double,3> read_tensor;
  read_tensor(0,0) = 0.1;
  read_tensor(0,1) = 0.2;
  read_tensor(0,2) = 0.3;
  read_tensor(1,0) = 1.1;
  read_tensor(1,1) = 1.2;
  read_tensor(1,2) = 1.3;
  read_tensor(2,0) = 2.1;
  read_tensor(2,1) = 2.2;
  read_tensor(2,2) = 2.3;

  hd.write(read_tensor, "matrix");
  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_read_partial.hdf");
  REQUIRE(okay);

  vector<int> readSpec{1,-1};
  vector<double> readBuffer(3);

  hd2.read(readBuffer, readSpec, "matrix");
  REQUIRE(readBuffer[0] == read_tensor(1,0));
  REQUIRE(readBuffer[1] == read_tensor(1,1));
  REQUIRE(readBuffer[2] == read_tensor(1,2));

  readSpec[0] = -1;
  readSpec[1] = 2;
  hd2.read(readBuffer, readSpec, "matrix");
  REQUIRE(readBuffer[0] == read_tensor(0,2));
  REQUIRE(readBuffer[1] == read_tensor(1,2));
  REQUIRE(readBuffer[2] == read_tensor(2,2));

  readSpec[0] = 2;
  readSpec[1] = 0;
  hd2.read(readBuffer, readSpec, "matrix");
  REQUIRE(readBuffer[0] == read_tensor(2,0));
  hd2.close();
}
    
