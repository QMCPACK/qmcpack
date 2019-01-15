
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "io/hdf_archive.h"
#include <vector>
#include "boost/multi_array.hpp"


using std::vector;

using namespace qmcplusplus;

TEST_CASE("hdf_multiarray_one_dim", "[hdf]")
{

  hdf_archive hd;
  hd.create("test_stl.hdf");

  boost::multi_array<double,1> v(boost::extents[3]);
  v[0] = 2.3;
  v[1] = -100.3;
  v[2] = 135.22;

  bool okay = hd.write(v, "boost_multiarray_one_dim");
  REQUIRE(okay);

  hd.close();

  hdf_archive hd2;
  okay = hd2.open("test_stl.hdf");
  REQUIRE(okay);

  boost::multi_array<double,1> v2(boost::extents[3]);
  okay = hd2.read(v2, "boost_multiarray_one_dim");
  REQUIRE(v2.size() == 3);
  for (int i = 0; i < v.size(); i++)
  {
    REQUIRE(v[i] == v2[i]);
  }

}

