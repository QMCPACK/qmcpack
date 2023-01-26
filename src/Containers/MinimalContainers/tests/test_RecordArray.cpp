//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory

//////////////////////////////////////////////////////////////////////////////////////

#include "MinimalContainers/RecordArray.hpp"

#include "catch.hpp"


namespace qmcplusplus
{
TEST_CASE("RecordArray basics", "[containers]")
{
  RecordArray<double> records;

  int nentry = 1;
  int nparam = 2;
  records.resize(nentry, nparam);

  records[0][0] = 1.1;
  records[0][1] = 1.2;

  REQUIRE(records.getNumOfParams() == 2);
  REQUIRE(records.getNumOfEntries() == 1);

  CHECK(records[0][0] == Approx(1.1));
  CHECK(records[0][1] == Approx(1.2));
}

} // namespace qmcplusplus
