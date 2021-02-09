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
  std::unique_ptr<Resource> lendResource();
  void takebackResource(std::unique_ptr<Resource>&& res);

  void rewind(size_t cursor = 0) { cursor_index_ = cursor; }

  bool is_lent() const { return is_lent_; }

  void lock() { is_lent_ = true; }
  void unlock() { is_lent_ = false; }

private:
  const std::string name_;
  size_t cursor_index_;
  std::vector<std::unique_ptr<Resource>> collection_;
  bool is_lent_ = false;
};

template <class CONSUMER>
class ResourceCollectionLock
{
public:
  ResourceCollectionLock(ResourceCollection& res_ref, CONSUMER& consumer_ref, size_t cursor = 0)
    : resource(res_ref), consumer(consumer_ref), cursor_begin_(cursor)
  {
    resource.rewind(cursor_begin_);
    consumer.acquireResource(resource);
  }

  ~ResourceCollectionLock()
  {
    resource.rewind(cursor_begin_);
    consumer.releaseResource(resource);
  }

  ResourceCollectionLock(const ResourceCollectionLock&) = delete;
  ResourceCollectionLock(ResourceCollectionLock&&) = delete;

private:
  ResourceCollection& resource;
  CONSUMER& consumer;
  const size_t cursor_begin_;
};

} // namespace qmcplusplus
#endif
