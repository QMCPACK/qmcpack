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
#include "type_traits/RefVectorWithLeader.h"

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

  size_t addResource(std::unique_ptr<Resource>&& res, bool noprint = false);
  std::unique_ptr<Resource> lendResource();
  void takebackResource(std::unique_ptr<Resource>&& res);

  void rewind(size_t cursor = 0) { cursor_index_ = cursor; }

  bool empty() const { return collection_.size() == 0; }

private:
  const std::string name_;
  size_t cursor_index_;
  std::vector<std::unique_ptr<Resource>> collection_;
};

template<class CONSUMER>
class ResourceCollectionLock
{
public:
  ResourceCollectionLock(ResourceCollection& res_ref, CONSUMER& consumer_ref, size_t cursor = 0)
      : resource(res_ref), consumer(consumer_ref), cursor_begin_(cursor), active(!res_ref.empty())
  {
    if (active)
    {
      resource.rewind(cursor_begin_);
      consumer.acquireResource(resource);
    }
  }

  ~ResourceCollectionLock()
  {
    if (active)
    {
      resource.rewind(cursor_begin_);
      consumer.releaseResource(resource);
    }
  }

  ResourceCollectionLock(const ResourceCollectionLock&) = delete;
  ResourceCollectionLock(ResourceCollectionLock&&)      = delete;

private:
  ResourceCollection& resource;
  CONSUMER& consumer;
  const size_t cursor_begin_;
  const bool active;
};


template<class CONSUMER>
class ResourceCollectionTeamLock
{
public:
  ResourceCollectionTeamLock(ResourceCollection& res_ref,
                             const RefVectorWithLeader<CONSUMER>& consumer_ref,
                             size_t cursor = 0)
      : resource(res_ref), consumer(consumer_ref), cursor_begin_(cursor), active(!res_ref.empty())
  {
    if (active)
    {
      resource.rewind(cursor_begin_);
      consumer.getLeader().acquireResource(resource, consumer);
    }
  }

  ~ResourceCollectionTeamLock()
  {
    if (active)
    {
      resource.rewind(cursor_begin_);
      consumer.getLeader().releaseResource(resource, consumer);
    }
  }

  ResourceCollectionTeamLock(const ResourceCollectionTeamLock&) = delete;
  ResourceCollectionTeamLock(ResourceCollectionTeamLock&&)      = delete;

private:
  ResourceCollection& resource;
  const RefVectorWithLeader<CONSUMER>& consumer;
  const size_t cursor_begin_;
  const bool active;
};

} // namespace qmcplusplus
#endif
