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

#include "ResourceCollection.h"

namespace qmcplusplus
{
ResourceCollection::ResourceCollection(const std::string& name) : name_(name) {}

void ResourceCollection::printResources() {}

size_t ResourceCollection::addResource(std::unique_ptr<Resource>&& res)
{
  size_t id = collection.size();
  collection.emplace_back(std::move(res));
  return id;
}

std::unique_ptr<Resource> ResourceCollection::lendResource(size_t id) { return std::move(collection[id]); }

void ResourceCollection::takebackResource(size_t id, std::unique_ptr<Resource>&& res)
{
  collection[id] = std::move(res);
}

} // namespace qmcplusplus
