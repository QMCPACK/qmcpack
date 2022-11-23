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

TEST_CASE("makeRefVector", "[type_traits]")
{
  struct Dummy
  {
    double d;
    std::string s;
  };

  struct DerivedDummy : Dummy
  {
    float f;
  };

  std::vector<DerivedDummy> ddvec;
  for (int i = 0; i < 3; ++i)
    ddvec.push_back(DerivedDummy());

  auto bdum  = makeRefVector<Dummy>(ddvec);
  auto bdum2 = makeRefVector<Dummy>(ddvec);
  CHECK(std::is_same<RefVector<Dummy>, decltype(bdum2)>::value);

  auto bdum3 = makeRefVector<decltype(ddvec)::value_type>(ddvec);
}

TEST_CASE("convertUPtrToRefvector", "[type_traits]")
{
  struct Dummy
  {
    double d;
    std::string s;
  };

  struct DerivedDummy : public Dummy
  {
    DerivedDummy() : e(0) {}
    double e;
  };

  struct OtherDummy
  {
    double f;
  };

  UPtrVector<Dummy> uvec;
  for (int i = 0; i < 3; ++i)
    uvec.emplace_back(std::make_unique<Dummy>());

  RefVector<Dummy> rdum(convertUPtrToRefVector(uvec));
  auto rdum2 = convertUPtrToRefVector(uvec);
  CHECK(std::is_same_v<decltype(rdum), decltype(rdum2)>);
  // Testing to make sure meta programming stops potential ambiguous template resolution.
  UPtrVector<DerivedDummy> ddv;
  ddv.emplace_back(std::make_unique<DerivedDummy>());
  ddv.emplace_back(std::make_unique<DerivedDummy>());

  auto dummy_ref_vec = convertUPtrToRefVector(ddv);
  auto d_ref_vec     = convertUPtrToRefVector<Dummy>(ddv);
  auto dd_ref_vec    = convertUPtrToRefVector(ddv);
  CHECK(!std::is_same_v<decltype(d_ref_vec), decltype(dd_ref_vec)>);

  // This should cause a compilation error. a DerivedDummy cannot be converted to an OtherDummy.
  // RefVector<OtherDummy> od_ref_vec = convertUPtrToRefVector<OtherDummy>(ddv);
}

TEST_CASE("convertPtrToRefvectorSubset", "[type_traits]")
{
  struct Dummy2
  {
    Dummy2(int j) : i(j) {}
    int i;
  };

  std::vector<Dummy2*> pvec;
  for (int i = 0; i < 5; ++i)
    pvec.push_back(new Dummy2(i));

  auto rdum = convertPtrToRefVectorSubset(pvec, 1, 4);

  CHECK(rdum.size() == 4);
  CHECK(rdum[0].get().i == 1);

  for (int i = 0; i < 5; ++i)
    delete pvec[i];
}
} // namespace qmcplusplus
