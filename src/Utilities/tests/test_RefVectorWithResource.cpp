//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "RefVectorWithResource.h"
#include <string>
#include <vector>

namespace qmcplusplus
{
namespace
{
class TestResource : public Resource
{
public:
  TestResource(const std::string& name, int val = 0) : Resource(name), data_(val) {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<TestResource>(*this); }

  int getData() const { return data_; }
  void setData(int val) { data_ = val; }

private:
  int data_;
};

struct ElementBase
{
  virtual ~ElementBase() = default;
  int id = 0;
};

struct ElementDerived : public ElementBase
{
  int extra_value = 0;
};
} // namespace

TEST_CASE("RefVectorWithResource basic operations", "[utilities]")
{
  TestResource res("test_res", 42);
  RefVectorWithResource<int> ref_vec(res);

  CHECK(ref_vec.empty());
  CHECK(ref_vec.size() == 0);
  CHECK(ref_vec.getResource().getName() == "test_res");

  int a = 1;
  int b = 2;
  int c = 3;

  ref_vec.push_back(a);
  ref_vec.push_back(b);
  ref_vec.push_back(c);

  REQUIRE(ref_vec.size() == 3);
  CHECK(ref_vec[0] == 1);
  CHECK(ref_vec[1] == 2);
  CHECK(ref_vec[2] == 3);

  // Modify element through ref_vec
  ref_vec[1] = 20;
  CHECK(b == 20);

  // Check getResource
  auto& captured_res = ref_vec.getResource();
  CHECK(captured_res.getName() == "test_res");
}

TEST_CASE("RefVectorWithResource constructors and CTAD", "[utilities]")
{
  TestResource res("test_res", 100);
  int x = 10, y = 20, z = 30;

  // Constructor with const BaseVec&
  std::vector<std::reference_wrapper<int>> base_vec{x, y, z};
  RefVectorWithResource<int> vec1(res, base_vec);
  REQUIRE(vec1.size() == 3);
  CHECK(vec1[0] == 10);
  CHECK(vec1[1] == 20);
  CHECK(vec1[2] == 30);

  // Constructor with BaseVec&&
  RefVectorWithResource<int> vec2(res, std::vector<std::reference_wrapper<int>>{x, y});
  REQUIRE(vec2.size() == 2);
  CHECK(vec2[0] == 10);
  CHECK(vec2[1] == 20);

  // CTAD
  RefVectorWithResource vec_deduced(res, base_vec);
  static_assert(std::is_same_v<decltype(vec_deduced), RefVectorWithResource<int>>);
  CHECK(vec_deduced.size() == 3);

  // Copy construction and assignment
  RefVectorWithResource<int> vec_copied(vec1);
  REQUIRE(vec_copied.size() == 3);
  CHECK(vec_copied.getResource().getName() == "test_res");

  RefVectorWithResource<int> vec_assigned(res);
  vec_assigned = vec1;
  REQUIRE(vec_assigned.size() == 3);
  CHECK(vec_assigned[2] == 30);

  // Move construction and assignment
  RefVectorWithResource<int> vec_moved(std::move(vec_copied));
  REQUIRE(vec_moved.size() == 3);
  CHECK(vec_moved[0] == 10);

  RefVectorWithResource<int> vec_move_assigned(res);
  vec_move_assigned = std::move(vec2);
  REQUIRE(vec_move_assigned.size() == 2);
  CHECK(vec_move_assigned[1] == 20);
}

TEST_CASE("RefVectorWithResource getCastedResource", "[utilities]")
{
  TestResource res("poly_res", 123);
  int elem = 5;
  RefVectorWithResource<int> ref_vec(res, {elem});

  TestResource& casted_res = ref_vec.getCastedResource<TestResource>();
  CHECK(casted_res.getData() == 123);
  casted_res.setData(456);
  CHECK(res.getData() == 456);

  // Const RefVectorWithResource
  const auto& const_ref_vec = ref_vec;
  CHECK(const_ref_vec.getResource().getName() == "poly_res");
  const TestResource& const_casted_res = const_ref_vec.getCastedResource<TestResource>();
  CHECK(const_casted_res.getData() == 456);
  CHECK(const_ref_vec[0] == 5);
}

TEST_CASE("RefVectorWithResource getCastedElement", "[utilities]")
{
  DummyResource dummy_res("dummy");
  ElementDerived d1;
  ElementDerived d2;
  d1.id          = 1;
  d1.extra_value = 100;
  d2.id          = 2;
  d2.extra_value = 200;

  RefVectorWithResource<ElementBase> vec(dummy_res, {d1, d2});
  REQUIRE(vec.size() == 2);

  CHECK(vec[0].id == 1);
  CHECK(vec[1].id == 2);

  ElementDerived& elem0 = vec.getCastedElement<ElementDerived>(0);
  CHECK(elem0.extra_value == 100);
  elem0.extra_value = 150;
  CHECK(d1.extra_value == 150);

  ElementDerived& elem1 = vec.getCastedElement<ElementDerived>(1);
  CHECK(elem1.extra_value == 200);

  // Test const qualifier matching
  const ElementDerived cd1 = d1;
  const ElementDerived cd2 = d2;
  RefVectorWithResource<const ElementBase> const_vec(dummy_res, {cd1, cd2});
  REQUIRE(const_vec.size() == 2);

  const ElementDerived& celem0 = const_vec.getCastedElement<const ElementDerived>(0);
  CHECK(celem0.extra_value == 150);
  const ElementDerived& celem1 = const_vec.getCastedElement<const ElementDerived>(1);
  CHECK(celem1.extra_value == 200);
}

} // namespace qmcplusplus
