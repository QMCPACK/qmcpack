//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "hdf/hdf_path.h"

using namespace qmcplusplus;

TEST_CASE("hdf_path_constructor", "[hdf]")
{
  auto gold_string("a/b/c");
  const hdf_path gold_path(gold_string);
  REQUIRE(gold_path.string() == gold_string);

  hdf_path path_sv(std::string_view{gold_string});
  REQUIRE(path_sv.string() == gold_string);
  REQUIRE(path_sv == gold_path);

  hdf_path path_s(std::string{gold_string});
  REQUIRE(path_s.string() == gold_string);
  REQUIRE(path_s == gold_path);

  hdf_path path_path(gold_path);
  REQUIRE(path_path.string() == gold_string);
  REQUIRE(path_path == gold_path);
}

TEST_CASE("hdf_path_append", "[hdf]")
{
  hdf_path path;
  REQUIRE(path.string().empty());
  REQUIRE(path == hdf_path());

  path.append("a");
  REQUIRE(path.string() == "a");
  REQUIRE(path == hdf_path("a"));

  path /= "b";
  REQUIRE(path.string() == "a/b");
  REQUIRE(path == hdf_path("a/b"));

  hdf_path result = path / std::string_view("c");
  REQUIRE(result.string() == "a/b/c");
  REQUIRE(result == hdf_path("a/b/c"));

  result = path / std::string("c");
  REQUIRE(result.string() == "a/b/c");
  REQUIRE(result == hdf_path("a/b/c"));

  result = path / "c/";
  REQUIRE(result.string() == "a/b/c/");
  REQUIRE(result == hdf_path("a/b/c/"));

  result = result / "d";
  REQUIRE(result.string() == "a/b/c/d");
  REQUIRE(result == hdf_path("a/b/c/d"));

  result = path / "/c";
  REQUIRE(result.string() == "/c");
  REQUIRE(result == hdf_path("/c"));

  result /= hdf_path{"d"};
  REQUIRE(result.string() == "/c/d");
  REQUIRE(result == hdf_path("/c/d"));
}

TEST_CASE("hdf_path_concat", "[hdf]")
{
  hdf_path path;
  REQUIRE(path.string().empty());
  REQUIRE(path == hdf_path());

  path.concat("a");
  REQUIRE(path.string() == "a");
  REQUIRE(path == hdf_path("a"));

  path += "b";
  REQUIRE(path.string() == "ab");
  REQUIRE(path == hdf_path("ab"));

  path += hdf_path{"c"};
  REQUIRE(path.string() == "abc");
  REQUIRE(path == hdf_path("abc"));
}

TEST_CASE("hdf_path_remove_path", "[hdf]")
{
  hdf_path path("/a/b/c/d");
  REQUIRE(path == hdf_path("/a/b/c/d"));

  path.remove_subgroup();
  REQUIRE(path == hdf_path("/a/b/c"));

  path.remove_subgroup();
  REQUIRE(path == hdf_path("/a/b"));

  path.remove_subgroup();
  REQUIRE(path == hdf_path("/a"));

  path.remove_subgroup();
  REQUIRE(path == hdf_path("/"));

  path.remove_subgroup();
  REQUIRE(path == hdf_path("/"));

  hdf_path relative_path("a/b/c/d");
  REQUIRE(relative_path == hdf_path("a/b/c/d"));

  relative_path.remove_subgroup();
  REQUIRE(relative_path == hdf_path("a/b/c"));

  relative_path.remove_subgroup();
  REQUIRE(relative_path == hdf_path("a/b"));

  relative_path.remove_subgroup();
  REQUIRE(relative_path == hdf_path("a"));

  relative_path.remove_subgroup();
  REQUIRE(relative_path == hdf_path(""));

  relative_path.remove_subgroup();
  REQUIRE(relative_path == hdf_path(""));
}

TEST_CASE("hdf_path_replace_path", "[hdf]")
{
  hdf_path path("/a/b");
  REQUIRE(path == hdf_path("/a/b"));

  path.replace_subgroup("c");
  REQUIRE(path == hdf_path("/a/c"));

  path.replace_subgroup(std::string{"de"});
  REQUIRE(path == hdf_path("/a/de"));

  path.replace_subgroup(std::string_view{"fs"});
  REQUIRE(path == hdf_path("/a/fs"));
}

TEST_CASE("hdf_path_has_root_path", "[hdf]")
{
  hdf_path root_path("/a/b");
  REQUIRE(root_path.has_root_directory());

  hdf_path relative_path("a/b");
  REQUIRE(!relative_path.has_root_directory());
}
