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
  const std::string& getName() const { return name_; }
  void printResources();

  size_t addResource(std::unique_ptr<Resource>&& res);
  std::unique_ptr<Resource> lendResource(int index);
  void takebackResource(int index, std::unique_ptr<Resource>&& res);

private:
  const std::string name_;
  std::vector<std::unique_ptr<Resource>> collection;
};
} // namespace qmcplusplus
#endif
