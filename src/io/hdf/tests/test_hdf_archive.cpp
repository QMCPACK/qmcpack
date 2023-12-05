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


#include "catch.hpp"

#include "hdf/hdf_archive.h"
#include <vector>
#include <string>
#include <sstream>


using std::string;
using std::vector;

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
  bool b    = true;
  bool okay = hd.writeEntry(b, "bool");
  REQUIRE(okay);

  int i = 23;
  okay  = hd.writeEntry(i, "int");
  REQUIRE(okay);

  float f = -2.3;
  okay    = hd.writeEntry(f, "float");
  REQUIRE(okay);

  const double d = 4.5;
  okay           = hd.writeEntry(d, "double");
  REQUIRE(okay);

  std::complex<float> cf(2.3, 3.4);
  okay = hd.writeEntry(cf, "complex float");
  REQUIRE(okay);

  hd.close();

  // Use the internally checked writes

  hdf_archive hd3;
  hd3.create("test_simple_data.hdf");
  hd3.write(b, "bool");
  hd3.write(i, "int");
  hd3.write(f, "float");
  hd3.write(d, "double");
  hd3.writeEntry(cf, "complex float");
  hd3.close();

  // now read the file and ensure the values are the same

  hdf_archive hd2;
  hd2.open("test_simple_data.hdf");
  bool b2 = false;
  okay    = hd2.readEntry(b2, "bool");
  REQUIRE(okay);
  REQUIRE(b == b2);

  int i2;
  okay = hd2.readEntry(i2, "int");
  REQUIRE(okay);
  REQUIRE(i == i2);

  // deliberately out of order
  double d2;
  okay = hd2.readEntry(d2, "double");
  REQUIRE(okay);
  REQUIRE(d == d2);

  double f2;
  okay = hd2.readEntry(f2, "float");
  REQUIRE(okay);
  REQUIRE(f == f2);

  std::complex<float> cf2;
  okay = hd2.readEntry(cf2, "complex float");
  REQUIRE(okay);
  REQUIRE(cf == cf2);

  // check an error occurs for non-existent entry
  int i666;
  okay = hd2.readEntry(i666, "not an entry");
  REQUIRE(!okay);

  hd2.close();

  // now read the file and ensure the values are the same

  hdf_archive hd4;
  hd4.open("test_simple_data.hdf");
  bool b4 = false;
  hd4.read(b4, "bool");
  REQUIRE(b == b4);

  int i4;
  hd4.read(i4, "int");
  REQUIRE(i == i4);

  // deliberately out of order
  double d4;
  hd4.read(d4, "double");
  REQUIRE(d == d4);

  float f4;
  hd4.read(f4, "float");
  REQUIRE(f == f4);

  std::complex<float> cf4;
  hd4.read(cf4, "complex float");
  REQUIRE(cf == cf4);

  //read scalar
  std::vector<int> datashape;
  okay = hd4.getShape<double>("double", datashape);
  REQUIRE(okay);
  REQUIRE(datashape.size() == 0);

  okay = hd4.getShape<std::complex<float>>("complex float", datashape);
  REQUIRE(okay);
  REQUIRE(datashape.size() == 0);

  okay = hd4.getShape<float>("complex float", datashape);
  REQUIRE(okay);
  REQUIRE(datashape.size() == 1);
  REQUIRE(datashape[0] == 2);

  hd4.close();
}

TEST_CASE("hdf_archive_vector", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_stl.hdf");

  vector<double> v(3);
  v[0] = 2.3;
  v[1] = -100.3;
  v[2] = 135.22;

  bool okay = hd.writeEntry(v, "vector_double");
  REQUIRE(okay);

  const vector<double> v_const(v);
  okay = hd.writeEntry(v_const, "vector_double_const");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_stl.hdf");
  REQUIRE(okay);

  vector<double> v2;
  okay = hd2.readEntry(v2, "vector_double");
  REQUIRE(v2.size() == 3);
  for (int i = 0; i < v.size(); i++)
    CHECK(v[i] == v2[i]);

  vector<double> v2_for_const;
  okay = hd2.readEntry(v2_for_const, "vector_double_const");
  REQUIRE(v2_for_const.size() == 3);
  for (int i = 0; i < v_const.size(); i++)
    CHECK(v_const[i] == v2_for_const[i]);
}

TEST_CASE("hdf_archive_group", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_group.hdf");

  int i     = 3;
  bool okay = hd.writeEntry(i, "int");
  REQUIRE(okay);

  CHECK(hd.group_path_as_string() == "");

  hd.push("name1");

  CHECK(hd.group_path_as_string() == "name1");

  int j = 3;
  okay  = hd.writeEntry(j, "int2");
  REQUIRE(okay);

  hd.push("name2");
  CHECK(hd.group_path_as_string() == "name1/name2");

  hd.close();

  // Check that opening a group on a closed file throws an exception
  REQUIRE_THROWS(hd.push("group"));

  hdf_archive hd2;
  hd2.open("test_group.hdf");
  bool int_is_group = hd2.is_group("int");
  REQUIRE(int_is_group == false);

  bool name1_is_group = hd2.is_group("name1");
  REQUIRE(name1_is_group);

  int j2 = 0;
  okay   = hd2.readEntry(j2, "name1/int2");
  REQUIRE(okay);
  REQUIRE(j2 == j);

  int j3 = 0;
  hd2.push("name1", false);
  okay = hd2.readEntry(j3, "int2");
  REQUIRE(okay);
  REQUIRE(j3 == j);

  REQUIRE_THROWS(hd2.push("nonexistent_group", false));

  hd2.close();
}

TEST_CASE("hdf_archive_scalar_convert", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_scalar_convert.hdf");

  TinyVector<double, 1> v0(1);
  bool okay = hd.writeEntry(v0, "tiny_vector_one");
  REQUIRE(okay);

  double v1(1);
  okay = hd.writeEntry(v1, "tiny_scalar_one");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  hd2.open("test_scalar_convert.hdf");

  TinyVector<double, 1> v2(0);
  okay = hd2.readEntry(v2, "tiny_scalar_one");
  REQUIRE(okay);

  double v3(0);
  okay = hd2.readEntry(v3, "tiny_vector_one");
  REQUIRE(okay);

  REQUIRE(v0[0] == v3);
  REQUIRE(v1 == v2[0]);
}

TEST_CASE("hdf_archive_tiny_vector", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_tiny_vector.hdf");

  TinyVector<double, 2> v(2);
  TinyVector<std::complex<double>, 2> v_cplx(2);

  v_cplx[0] = v[0] = 1.2;
  v_cplx[1] = v[1] = 1.3;

  bool okay = hd.writeEntry(v, "tiny_vector_double");
  REQUIRE(okay);
  okay = hd.writeEntry(v_cplx, "tiny_vector_complex_double");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  hd2.open("test_tiny_vector.hdf");

  TinyVector<double, 2> v2;
  TinyVector<std::complex<double>, 2> v2_cplx(2);
  okay = hd2.readEntry(v2, "tiny_vector_double");
  REQUIRE(okay);
  okay = hd2.readEntry(v2_cplx, "tiny_vector_complex_double");
  REQUIRE(okay);
  for (int i = 0; i < v.size(); i++)
  {
    REQUIRE(v[i] == v2[i]);
    REQUIRE(v_cplx[i] == v2_cplx[i]);
  }
}

TEST_CASE("hdf_archive_tensor", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_tensor.hdf");

  Tensor<float, 2> v(2);

  v(0, 1) = 1.2f;
  v(1, 0) = 1.3f;
  v(0, 1) = -2.3f;
  v(1, 1) = 10.0f;

  bool okay = hd.writeEntry(v, "tiny_tensor_float");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  hd2.open("test_tensor.hdf");

  Tensor<float, 2> v2;
  okay = hd2.readEntry(v2, "tiny_tensor_float");
  REQUIRE(okay);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      REQUIRE(v(i, j) == v2(i, j));
    }
  }
}

TEST_CASE("hdf_archive_string", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_string.hdf");

  string s("this is a test");
  bool okay = hd.writeEntry(s, "string");
  REQUIRE(okay);


  std::ostringstream o;
  o << "Another test" << std::endl;
  okay = hd.writeEntry(o, "ostringstream");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_string.hdf");
  REQUIRE(okay);
  string s2;
  okay = hd2.readEntry(s2, "string");
  REQUIRE(okay);
  REQUIRE(s == s2);

  string o2;
  okay = hd2.readEntry(o2, "ostringstream");
  REQUIRE(okay);
  REQUIRE(o.str() == o2);
}

TEST_CASE("hdf_archive_string_vector", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_string_vector.hdf");

  std::vector<std::string> strings;
  strings.push_back("first");
  // One entry should be longer than 15 characters to avoid the short
  // string optimization and allocate space for the string on the heap
  strings.push_back("really long string");

  bool okay = hd.writeEntry(strings, "string_vector");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_string_vector.hdf");
  REQUIRE(okay);

  std::vector<std::string> strings2;
  okay = hd2.readEntry(strings2, "string_vector");
  REQUIRE(okay);

  REQUIRE(strings2.size() == 2);
  REQUIRE(strings2[0] == "first");
  REQUIRE(strings2[1] == "really long string");
}

TEST_CASE("hdf_archive_dataset_existence_checking", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_dataset_existence_checking.hdf");

  std::vector<uint64_t> numbers;
  numbers.push_back(123456);
  std::string ds_tag = "numbers_vector";

  bool okay = hd.writeEntry(numbers, ds_tag);
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_dataset_existence_checking.hdf");
  REQUIRE(okay);

  REQUIRE(hd2.is_dataset(ds_tag));
  REQUIRE(!hd2.is_dataset("tag_doesnt_exist"));
}

TEST_CASE("hdf_archive_dataset_type_checking", "[hdf]")
{
  hdf_archive hd;
  hd.create("test_dataset_type_checking.hdf");

  std::vector<uint64_t> numbers;
  numbers.push_back(123456);
  std::string ds_tag = "numbers_vector";

  bool okay = hd.writeEntry(numbers, ds_tag);
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_dataset_type_checking.hdf");
  REQUIRE(okay);

  bool is_correct_type = hd2.is_dataset_of_type<uint64_t>(ds_tag);
  REQUIRE(is_correct_type);
  is_correct_type = hd2.is_dataset_of_type<int64_t>(ds_tag);
  REQUIRE(is_correct_type == false);
  REQUIRE_THROWS_AS(hd2.is_dataset_of_type<uint64_t>("tag_doesnt_exist"), std::runtime_error);
}
