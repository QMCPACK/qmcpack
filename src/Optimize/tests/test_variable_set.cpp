//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "complex_approx.hpp"

#include "Optimize/VariableSet.h"

#include <stdio.h>
#include <string>

using std::string;

namespace optimize
{
TEST_CASE("VariableSet empty", "[optimize]")
{
  VariableSet vs;

  REQUIRE(vs.is_optimizable() == false);
  REQUIRE(vs.size_of_active() == 0);
  REQUIRE(vs.find("something") == vs.NameAndValue.end());
  REQUIRE(vs.getIndex("something") == -1);
}

TEST_CASE("VariableSet one", "[optimize]")
{
  VariableSet vs;
  VariableSet::value_type first_val(1.123456789);
  vs.insert("first", first_val);
  std::vector<std::string> names{"first"};
  vs.activate(names.begin(), names.end(), true);

  REQUIRE(vs.is_optimizable() == true);
  REQUIRE(vs.size_of_active() == 1);
  REQUIRE(vs.getIndex("first") == 0);
  REQUIRE(vs.name(0) == "first");
  double first_val_real = 1.123456789;
  REQUIRE(std::real(vs[0] ) == Approx(first_val_real));

  std::ostringstream o;
  vs.print(o, 0, false);
  //std::cout << o.str() << std::endl;
  #ifdef QMC_COMPLEX
  REQUIRE(o.str() == "first  (1.123457e+00,0.000000e+00) 0 1  ON 0\n");
  #else
  REQUIRE(o.str() == "first                 1.123457e+00 0 1  ON 0\n");
  #endif

  std::ostringstream o2;
  vs.print(o2, 1, true);
  //std::cout << o2.str() << std::endl;

  #ifdef QMC_COMPLEX
  char formatted_output[] = "  Name                        Value Type Recompute Use Index\n"
                            " ----- ---------------------------- ---- --------- --- -----\n"
                            " first  (1.123457e+00,0.000000e+00)    0         1  ON     0\n";


  REQUIRE(o2.str() == formatted_output);
  #else
  char formatted_output[] = "  Name                        Value Type Recompute Use Index\n"
                            " ----- ---------------------------- ---- --------- --- -----\n"
                            " first                 1.123457e+00    0         1  ON     0\n";


  REQUIRE(o2.str() == formatted_output);
  #endif
}

TEST_CASE("VariableSet output", "[optimize]")
{
  VariableSet vs;
  VariableSet::value_type first_val(11234.56789);
  VariableSet::value_type second_val(0.000256789);
  VariableSet::value_type third_val(-1.2);
  vs.insert("s", first_val);
  vs.insert("second", second_val);
  vs.insert("really_long_name", third_val);
  std::vector<std::string> names{"s", "second", "really_long_name"};
  vs.activate(names.begin(), names.end(), true);

  std::ostringstream o;
  vs.print(o, 0, true);
  //std::cout << o.str() << std::endl;

  #ifdef QMC_COMPLEX
  char formatted_output[] = "            Name                        Value Type Recompute Use Index\n"
                            "---------------- ---------------------------- ---- --------- --- -----\n"
                            "               s  (1.123457e+04,0.000000e+00)    0         1  ON     0\n"
                            "          second  (2.567890e-04,0.000000e+00)    0         1  ON     1\n"
                            "really_long_name (-1.200000e+00,0.000000e+00)    0         1  ON     2\n";
  #else
  char formatted_output[] = "            Name                        Value Type Recompute Use Index\n"
                            "---------------- ---------------------------- ---- --------- --- -----\n"
                            "               s                 1.123457e+04    0         1  ON     0\n"
                            "          second                 2.567890e-04    0         1  ON     1\n"
                            "really_long_name                -1.200000e+00    0         1  ON     2\n";
  #endif

  REQUIRE(o.str() == formatted_output);
}

} // namespace optimize
