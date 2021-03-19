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
#include <Host/OutputManager.h>

namespace qmcplusplus
{
ResourceCollection::ResourceCollection(const std::string& name) : name_(name), cursor_index_(0) {}

ResourceCollection::ResourceCollection(const ResourceCollection& ref) : name_(ref.getName()), cursor_index_(0)
{
  for (auto& res : ref.collection_)
    addResource(std::unique_ptr<Resource>(res->makeClone()), true);
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

size_t ResourceCollection::addResource(std::unique_ptr<Resource>&& res, bool noprint)
{
  size_t index              = collection_.size();
  res->index_in_collection_ = index;
  if (!noprint)
    app_debug_stream() << "Multi walker shared resource \"" << res->getName() << "\" created in resource collection \""
                << name_ << "\" index " << index << std::endl;
  collection_.emplace_back(std::move(res));
  return index;
}

std::unique_ptr<Resource> ResourceCollection::lendResource()
{
  if (cursor_index_ >= collection_.size())
    throw std::runtime_error("ResourceCollection::lendResource BUG no more resource to lend.");
  if (cursor_index_ != collection_[cursor_index_]->index_in_collection_)
    throw std::runtime_error(
        "ResourceCollection::lendResource BUG mismatched cursor index and recorded index in the resource.");
  return std::move(collection_[cursor_index_++]);
}

/** Self explanatory but you must rewind or at least somewhat.
 */  
void ResourceCollection::takebackResource(std::unique_ptr<Resource>&& res)
{
  if (cursor_index_ >= collection_.size())
    throw std::runtime_error("ResourceCollection::takebackResource BUG cannot take back resources more than owned. Don't forget the intuitive call to resource_collection.rewind");
  if (cursor_index_ != res->index_in_collection_)
    throw std::runtime_error(
        "ResourceCollection::takebackResource BUG mismatched cursor index and recorded index in the resource.");
  collection_[cursor_index_++] = std::move(res);
}

} // namespace qmcplusplus
