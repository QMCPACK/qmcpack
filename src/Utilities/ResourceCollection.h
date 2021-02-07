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

#ifndef QMCPLUSPLUS_RESOURCECOLLECTION_H
#define QMCPLUSPLUS_RESOURCECOLLECTION_H

#include <string>
#include <memory>
#include <cstddef>
#include <vector>
#include "Resource.h"

namespace qmcplusplus
{
class ResourceCollection
{
public:
  ResourceCollection(const std::string& name);
  ResourceCollection(const ResourceCollection&);
  const std::string& getName() const { return name_; }
  size_t size() const { return collection_.size(); }
  void printResources() const;

  size_t addResource(std::unique_ptr<Resource>&& res);
  void rewind() { cursor_index_ = 0; }
  std::unique_ptr<Resource> lendResource();
  void takebackResource(std::unique_ptr<Resource>&& res);

private:
  const std::string name_;
  size_t cursor_index_;
  std::vector<std::unique_ptr<Resource>> collection_;
};
} // namespace qmcplusplus
#endif
