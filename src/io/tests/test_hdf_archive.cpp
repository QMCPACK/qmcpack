
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Message/catch_mpi_main.hpp"

#include "io/hdf_archive.h"
#include <vector>
#include <string>
#include <sstream>


using std::vector;
using std::string;

using namespace qmcplusplus;

TEST_CASE("hdf_archive_empty_file", "[hdf]")
{
  hdf_archive hd;
  bool okay = hd.create("test1.hdf");
  REQUIRE(okay);
  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test1.hdf");
  REQUIRE(okay);
  hd2.close();
}

// Ensure failure if the file we're trying to open does not exist.
TEST_CASE("hdf_archive_failed_to_open", "[hdf]")
{
  hdf_archive hd;
  bool okay = hd.open("should_not_exist.hdf");
  REQUIRE(okay == false);
  hd.close();
}

// Simple scalar data types in hdf_datatype.h
TEST_CASE("hdf_archive_simple_data", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_simple_data.hdf");
  int i = 23;
  bool okay = hd.write(i, "int");
  REQUIRE(okay);

  float f = -2.3;
  okay = hd.write(f, "float");
  REQUIRE(okay);

  double d = 4.5;
  okay = hd.write(d, "double");
  REQUIRE(okay);

  std::complex<float> cf(2.3,3.4);
  okay = hd.write(cf, "complex float");
  REQUIRE(okay);

  hd.close();

  // now read the file and ensure the values are the same

  hdf_archive hd2;
  hd2.open("test_simple_data.hdf");
  int i2;
  okay = hd2.read(i2, "int");
  REQUIRE(okay);
  REQUIRE(i == i2);

  // deliberately out of order
  double d2;
  okay = hd2.read(d2, "double");
  REQUIRE(okay);
  REQUIRE(d == d2);

  double f2;
  okay = hd2.read(f2, "float");
  REQUIRE(okay);
  REQUIRE(f == f2);

  std::complex<float> cf2;
  okay = hd2.read(cf2, "complex float");
  REQUIRE(okay);
  REQUIRE(cf == cf2);

  hd2.close();
}

TEST_CASE("hdf_archive_vector", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_stl.hdf");

  vector<double> v(3);
  v[0] = 2.3;
  v[1] = -100.3;
  v[2] = 135.22;

  bool okay = hd.write(v, "vector_double");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_stl.hdf");
  REQUIRE(okay);

  vector<double> v2;
  okay = hd2.read(v2, "vector_double");
  REQUIRE(v2.size() == 3);
  for (int i = 0; i < v.size(); i++)
  {
    REQUIRE(v[i] == v2[i]);
  }
}

TEST_CASE("hdf_archive_group", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_group.hdf");

  int i = 3;
  bool okay = hd.write(i, "int");
  REQUIRE(okay);

  hd.push("name1");

  int j = 3;
  okay = hd.write(j, "int2");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  hd2.open("test_group.hdf");
  bool int_is_group = hd2.is_group("int");
  REQUIRE(int_is_group == false);

  bool name1_is_group = hd2.is_group("name1");
  REQUIRE(name1_is_group);

  int j2 = 0;
  okay = hd2.read(j2, "name1/int2");
  REQUIRE(okay);
  REQUIRE(j2 == j);

  int j3 = 0;
  hd2.push("name1", false);
  okay = hd2.read(j3, "int2");
  REQUIRE(okay);
  REQUIRE(j3 == j);

  hd2.close();
}

TEST_CASE("hdf_archive_tiny_vector", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_tiny_vector.hdf");

  TinyVector<double, 2> v(2);

  v[0] = 1.2;
  v[1] = 1.3;

  bool okay = hd.write(v, "tiny_vector_double");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  hd2.open("test_tiny_vector.hdf");

  TinyVector<double, 2> v2;
  okay = hd2.read(v2, "tiny_vector_double");
  REQUIRE(okay);
  for (int i = 0; i < v.size(); i++)
  {
    REQUIRE(v[i] == v2[i]);
  }
}

TEST_CASE("hdf_archive_tensor", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_tensor.hdf");

  Tensor<float, 2> v(2);

  v(0,1) = 1.2f;
  v(1,0) = 1.3f;
  v(0,1) = -2.3f;
  v(1,1) = 10.0f;

  bool okay = hd.write(v, "tiny_tensor_float");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  hd2.open("test_tensor.hdf");

  Tensor<float, 2> v2;
  okay = hd2.read(v2, "tiny_tensor_float");
  REQUIRE(okay);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      REQUIRE(v(i,j) == v2(i,j));
    }
  }
}

TEST_CASE("hdf_archive_string", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_string.hdf");

  string s("this is a test");
  bool okay = hd.write(s, "string");
  REQUIRE(okay);


  std::ostringstream o;
  o << "Another test" << std::endl;
  okay = hd.write(o, "ostringstream");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_string.hdf");
  REQUIRE(okay);
  string s2;
  okay = hd2.read(s2, "string");
  REQUIRE(okay);
  REQUIRE(s == s2);

  string o2;
  okay = hd2.read(o2, "ostringstream");
  REQUIRE(okay);
  REQUIRE(o.str() == o2);
}
