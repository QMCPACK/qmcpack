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

  MemoryResource* makeClone() const override { return new MemoryResource(*this); }

  std::vector<int> data;
};

TEST_CASE("Resource", "[utilities]")
{
  auto mem_res = std::make_unique<MemoryResource>("test_res");
  mem_res->data.resize(5);

  std::unique_ptr<MemoryResource> res_copy(mem_res->makeClone());
  REQUIRE(res_copy->data.size() == 5);
}

class WFCResourceConsumer
{
public:
  void createResource(ResourceCollection& collection)
  {
    external_memory_handle = std::make_unique<MemoryResource>("test_res");
    external_memory_handle->data.resize(5);
    collection.addResource(std::move(external_memory_handle));
  }

  void acquireResource(ResourceCollection& collection)
  {
    auto res_ptr = dynamic_cast<MemoryResource*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error("WFCResourceConsumer::acquireResource dynamic_cast failed");
    external_memory_handle.reset(res_ptr);
  }

  void releaseResource(ResourceCollection& collection)
  {
    collection.takebackResource(std::move(external_memory_handle));
  }

  MemoryResource* getPtr() const { return external_memory_handle.get(); }

private:
  std::unique_ptr<MemoryResource> external_memory_handle;
};

TEST_CASE("ResourceCollection", "[utilities]")
{
  ResourceCollection res_collection("abc");
  WFCResourceConsumer wfc;
  REQUIRE(wfc.getPtr() == nullptr);

  wfc.createResource(res_collection);
  REQUIRE(wfc.getPtr() == nullptr);

  wfc.acquireResource(res_collection);
  REQUIRE(wfc.getPtr() != nullptr);
  REQUIRE(wfc.getPtr()->data.size() == 5);

  res_collection.rewind();
  wfc.releaseResource(res_collection);
  REQUIRE(wfc.getPtr() == nullptr);
}

} // namespace qmcplusplus
