//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OptimizableObject.h"
#include "VariableSet.h"
#include "io/hdf/hdf_archive.h"

namespace qmcplusplus
{

class FakeOptimizableObject : public OptimizableObject
{
public:
  FakeOptimizableObject(const std::string& my_name, double var1 = 0.0, double var2 = 0.0) : OptimizableObject(my_name)
  {
    myVars.insert("var1", var1);
    myVars.insert("var2", var2);
  }
  void checkInVariablesExclusive(opt_variables_type& active) override { active.insertFrom(myVars); }

  void resetParametersExclusive(const opt_variables_type& active) override {}

  void writeVariationalParameters(hdf_archive& hout) override
  {
    hout.push("FakeOptimizableObject");
    // Convert to double before writing to HDF.
    // If these get written as complex, but read as double, it can cause memory overwrites.
    // The HDF read call uses the file's type information to determine the amount data read
    //   (2 doubles for complex), but the variable size is just one double.
    double var1 = std::real(myVars["var1"]);
    double var2 = std::real(myVars["var2"]);
    hout.write(var1, "var1");
    hout.write(var2, "var2");
    hout.write(extra_data, "extra_data");
    hout.pop();
  }

  void readVariationalParameters(hdf_archive& hin) override
  {
    hin.push("FakeOptimizableObject", false);

    myVars.clear();
    double tmp_var1;
    double tmp_var2;
    // The HDF read call uses the file's type information to determine the amount data read into the variable.
    hin.read(tmp_var1, "var1");
    hin.read(tmp_var2, "var2");
    hin.read(extra_data, "extra_data");
    myVars.insert("var1", tmp_var1);
    myVars.insert("var2", tmp_var2);
    hin.pop();
  }

  double extra_data;
  opt_variables_type myVars;
};

TEST_CASE("Test OptimizableObject", "[wavefunction]")
{
  FakeOptimizableObject fake_a("functor_a");
  FakeOptimizableObject fake_b("functor_b");

  UniqueOptObjRefs opt_obj_refs;
  opt_obj_refs.push_back(fake_a);
  opt_obj_refs.push_back(fake_b);

  CHECK(opt_obj_refs.size() == 2);

  FakeOptimizableObject fake_c("functor_b");
  opt_obj_refs.push_back(fake_c);
  CHECK(opt_obj_refs.size() == 2);
  CHECK(opt_obj_refs[0].getName() == "functor_a");
  CHECK(opt_obj_refs[1].getName() == "functor_b");
}

TEST_CASE("OptimizableObject HDF output and input", "[wavefunction]")
{
  FakeOptimizableObject fake_a("functor_a", 1.1, 2.3);
  fake_a.extra_data = 3.4;

  hdf_archive hout;
  opt_variables_type opt_vars;
  fake_a.checkInVariablesExclusive(opt_vars);
  opt_vars.writeToHDF("opt_obj.h5", hout);

  fake_a.writeVariationalParameters(hout);

  FakeOptimizableObject fake_2a("functor_a");
  hdf_archive hin;
  opt_variables_type opt_vars2;
  fake_2a.checkInVariablesExclusive(opt_vars2);

  opt_vars2.readFromHDF("opt_obj.h5", hin);
  CHECK(std::real(opt_vars2["var1"]) == Approx(1.1));
  CHECK(std::real(opt_vars2["var2"]) == Approx(2.3));

  fake_2a.readVariationalParameters(hin);

  opt_variables_type opt_vars3;
  fake_2a.checkInVariablesExclusive(opt_vars3);
  CHECK(std::real(opt_vars3["var1"]) == Approx(1.1));
  CHECK(std::real(opt_vars3["var2"]) == Approx(2.3));
  CHECK(fake_2a.extra_data == Approx(3.4));
}
} // namespace qmcplusplus
