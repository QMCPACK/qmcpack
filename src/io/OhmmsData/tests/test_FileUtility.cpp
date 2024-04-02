//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/FileUtility.h"

#include <string>
#include <string_view>

TEST_CASE("FileUtility", "[xml]")
{
  std::string empty;
  REQUIRE(getExtension(empty) == std::string_view());

  std::string no_extension("this_filename_has_no_extension");
  REQUIRE(getExtension(no_extension) == std::string_view());

  std::string single_extension("filename.h5");
  REQUIRE(getExtension(single_extension) == "h5");

  std::string multiple_extensions("filename.p001.h5");
  REQUIRE(getExtension(multiple_extensions) == "h5");
}
