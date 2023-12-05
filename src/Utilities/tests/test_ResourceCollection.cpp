//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{
class MemoryResource : public Resource
{
public:
  MemoryResource(const std::string& name) : Resource(name) {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<MemoryResource>(*this); }

  std::vector<int> data;
};

TEST_CASE("Resource", "[utilities]")
{
  auto mem_res = std::make_unique<MemoryResource>("test_res");
  mem_res->data.resize(5);

  auto res_copy      = mem_res->makeClone();
  auto& res_copy_ref = dynamic_cast<MemoryResource&>(*res_copy);
  REQUIRE(res_copy_ref.data.size() == 5);
}

TEST_CASE("DummyResource", "[utilities]")
{
  DummyResource dummy;
  auto dummy2 = dummy.makeClone();
  REQUIRE(dummy2->getName() == "Dummy");
  DummyResource dummy_alt("dummy_alt_name");
  REQUIRE(dummy_alt.getName() == "dummy_alt_name");
}

class WFCResourceConsumer
{
public:
  void createResource(ResourceCollection& collection)
  {
    auto memory_handle = std::make_unique<MemoryResource>("test_res");
    memory_handle->data.resize(5);
    collection.addResource(std::move(memory_handle));
  }

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<WFCResourceConsumer>& wfcrc_list)
  {
    external_memory_handle = collection.lendResource<MemoryResource>();
  }

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<WFCResourceConsumer>& wfcrc_list)
  {
    collection.takebackResource(external_memory_handle);
  }

  auto& getResourceHandle() { return external_memory_handle; }

private:
  ResourceHandle<MemoryResource> external_memory_handle;
};

TEST_CASE("ResourceCollection", "[utilities]")
{
  ResourceCollection res_collection("abc");
  WFCResourceConsumer wfc, wfc1, wfc2;
  REQUIRE(wfc.getResourceHandle().hasResource() == false);

  wfc.createResource(res_collection);
  REQUIRE(wfc.getResourceHandle().hasResource() == false);

  RefVectorWithLeader wfc_list(wfc, {wfc, wfc1, wfc2});

  {
    ResourceCollectionTeamLock lock(res_collection, wfc_list);
    auto& res_handle = wfc.getResourceHandle();
    REQUIRE(res_handle);

    MemoryResource& mem_res             = res_handle;
    const MemoryResource& const_mem_res = res_handle;
    CHECK(mem_res.data.size() == 5);
    CHECK(const_mem_res.data.size() == 5);
  }

  REQUIRE(wfc.getResourceHandle().hasResource() == false);
}

} // namespace qmcplusplus
