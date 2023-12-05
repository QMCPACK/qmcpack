//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020  QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "DTModes.h"

namespace qmcplusplus
{
TEST_CASE("DTModes", "[particle]")
{
  DTModes modes;

  modes = DTModes::NEED_FULL_TABLE_ANYTIME | ~DTModes::NEED_TEMP_DATA_ON_HOST;

  CHECK(modes & DTModes::NEED_FULL_TABLE_ANYTIME);
  CHECK(!(modes & DTModes::NEED_TEMP_DATA_ON_HOST));
}

} // namespace qmcplusplus
