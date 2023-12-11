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
#include "hdf/hdf_error_suppression.h"

using namespace qmcplusplus;

TEST_CASE("hdf_error_suppression", "[hdf]")
{
  H5E_auto2_t err_func{};
  void* client_data{nullptr};
  // catch main already contains an instance of hdf_error_suppression.
  REQUIRE(hdf_error_suppression::enabled == true);
  H5Eget_auto2(H5E_DEFAULT, &err_func, &client_data);
  REQUIRE(client_data == nullptr);
}

