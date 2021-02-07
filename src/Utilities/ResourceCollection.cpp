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
#include <iostream>

namespace qmcplusplus
{
ResourceCollection::ResourceCollection(const std::string& name) : name_(name), cursor_index_(0) {}

ResourceCollection::ResourceCollection(const ResourceCollection& ref) : name_(ref.getName()), cursor_index_(0)
{
  for (auto& res : ref.collection_)
    addResource(std::unique_ptr<Resource>(res->makeClone()));
}

void ResourceCollection::printResources() const
{
  std::cout << "list resources in " << getName() << std::endl;
  std::cout << "-------------------------------" << std::endl;
  for (int i = 0; i < collection_.size(); i++)
    std::cout << "resource " << i << "    name: " << collection_[i]->getName()
              << "    address: " << collection_[i].get() << std::endl;
  std::cout << "-------------------------------" << std::endl << std::endl;
}

size_t ResourceCollection::addResource(std::unique_ptr<Resource>&& res)
{
  size_t id = collection_.size();
  collection_.emplace_back(std::move(res));
  return id;
}

std::unique_ptr<Resource> ResourceCollection::lendResource()
{
  if (cursor_index_ >= collection_.size())
    throw std::runtime_error("ResourceCollection::lendResource no more resource to lend. Bug in the code.");
  return std::move(collection_[cursor_index_++]);
}

void ResourceCollection::takebackResource(std::unique_ptr<Resource>&& res)
{
  if (cursor_index_ >= collection_.size())
    throw std::runtime_error("ResourceCollection::takebackResource cannot take back resources more than owned. Bug in the code.");
  collection_[cursor_index_++] = std::move(res);
}

} // namespace qmcplusplus
