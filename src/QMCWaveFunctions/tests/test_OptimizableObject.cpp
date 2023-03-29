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
#include "io/hdf/hdf_archive.h"

namespace qmcplusplus
{

class FakeOptimizableObject : public OptimizableObject
{
public:
  FakeOptimizableObject(const std::string& my_name) : OptimizableObject(my_name) {}
  void checkInVariablesExclusive(opt_variables_type& active) override
  {
    active.insert("var1", var1);
    active.insert("var2", var2);
  }

  void resetParametersExclusive(const opt_variables_type& active) override {}

  double var1 = 1.1;
  double var2 = 2.3;
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
  FakeOptimizableObject fake_a("functor_a");
  hdf_archive hout;
  opt_variables_type opt_vars;

  fake_a.openHDFToSave("opt_obj_vp.h5", hout);
  fake_a.saveVariationalParameters(hout, opt_vars);
  hout.close();

  FakeOptimizableObject fake_a2("functor_a");
  fake_a2.var1 = 0.0;
  fake_a2.var2 = 0.0;
  fake_a2.checkInVariablesExclusive(opt_vars);
  opt_vars.resetIndex();
  hdf_archive hin;
  fake_a2.openHDFToRead("opt_obj_vp.h5", hin);
  fake_a2.readVariationalParameters(hin, opt_vars);
  CHECK(opt_vars.find("var1")->second == ValueApprox(1.1));
  CHECK(opt_vars.find("var2")->second == ValueApprox(2.3));
}

} // namespace qmcplusplus
