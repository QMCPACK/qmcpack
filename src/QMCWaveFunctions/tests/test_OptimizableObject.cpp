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

namespace qmcplusplus
{

class FakeOptimizableObject : public OptimizableObject
{
public:
  FakeOptimizableObject(const std::string& my_name) : OptimizableObject(my_name) {}
  void checkInVariablesExclusive(opt_variables_type& active) override {}
  void resetParametersExclusive(const opt_variables_type& active) override {}
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
} // namespace qmcplusplus
