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

#ifndef QMCPLUSPLUS_RESOURCE_H
#define QMCPLUSPLUS_RESOURCE_H

#include <string>
#include <memory>

namespace qmcplusplus
{
class Resource
{
public:
  Resource(const std::string& name) : name_(name) {}
  virtual ~Resource()                                 = default;
  virtual std::unique_ptr<Resource> makeClone() const = 0;
  const std::string& getName() const { return name_; }

private:
  const std::string name_;
  int index_in_collection_ = -1;
  friend class ResourceCollection;
};

/** For the sake of generic code sometimes an dummy resource is needed to pass API.
 */
class DummyResource : public Resource
{
public:
  DummyResource() : Resource("Dummy") {}
  DummyResource(const std::string& name) : Resource(name) {}
  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<DummyResource>(*this); }
};
} // namespace qmcplusplus
#endif
