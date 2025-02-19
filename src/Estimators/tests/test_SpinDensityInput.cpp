//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "SpinDensityInput.h"
#include "ValidSpinDensityInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
TEST_CASE("SpinDensityInput::readXML", "[estimators]")
{
  using Lattice = CrystalLattice<OHMMS_PRECISION_FULL, OHMMS_DIM>;
  using input = testing::ValidSpinDensityInput;
  for (auto input_xml : input::xml)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    SpinDensityInput sdi(node);

    Lattice lattice;
    if (sdi.get_cell().explicitly_defined == true)
      lattice = sdi.get_cell();
    else
      lattice = testing::makeTestLattice();

    SpinDensityInput::DerivedParameters dev_par = sdi.calculateDerivedParameters(lattice);

    CHECK(dev_par.npoints == 1000);
    TinyVector<int, SpinDensityInput::DIM> grid(10, 10, 10);
    CHECK(dev_par.grid == grid);
    TinyVector<int, SpinDensityInput::DIM> gdims(100, 10, 1);
    CHECK(dev_par.gdims == gdims);
  }
}

} // namespace qmcplusplus
