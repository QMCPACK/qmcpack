//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include <complex>
#include "catch.hpp"
#include "type_traits/template_types.hpp"


namespace qmcplusplus
{

TEST_CASE("convertUPtrToRefvector", "[type_traits]")
{
  struct Dummy {
    double d;
    std::string s;
  };

  UPtrVector<Dummy> uvec;
  for(int i = 0; i < 3; ++i)
    uvec.emplace_back(std::make_unique<Dummy>());

  RefVector<Dummy> rdum(convertUPtrToRefVector(uvec));
}

}
