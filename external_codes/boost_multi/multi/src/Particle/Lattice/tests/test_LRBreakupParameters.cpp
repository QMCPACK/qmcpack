//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <string>

#include "catch.hpp"

#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/LRBreakupParameters.h"

namespace qmcplusplus
{
/** Lattice is defined but Open BC is also used.
 */
TEST_CASE("LRBreakupParameters", "[lattice]")
{
  LRBreakupParameters<double, 3> myLR;

  TinyVector<TinyVector<double, 3>, 3> R;

  R[0][0] = 1.0;
  R[1][1] = 1.0;
  R[2][2] = 1.0;

  myLR.SetLRCutoffs(R);

  CHECK(myLR.LR_kc == Approx(30.0));

  R[0][0] = 2.0;
  R[1][1] = 2.0;
  R[2][2] = 3.0;

  myLR.SetLRCutoffs(R);

  CHECK(myLR.LR_kc == Approx(15.0));
}

} // namespace qmcplusplus
