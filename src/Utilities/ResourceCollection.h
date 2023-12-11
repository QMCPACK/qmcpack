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
#include "ResourceHandle.h"
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
  bool empty() const { return collection_.size() == 0; }

  size_t addResource(std::unique_ptr<Resource>&& res, bool noprint = false);
  void printResources() const;

  template<class RS>
  ResourceHandle<RS> lendResource() { return dynamic_cast<RS&>(lendResourceImpl()); }

  template<class RS>
  void takebackResource(ResourceHandle<RS>& res_handle) { takebackResourceImpl(res_handle.release()); }

  void rewind(size_t cursor = 0) { cursor_index_ = cursor; }
private:
  Resource& lendResourceImpl();
  void takebackResourceImpl(Resource& res);

  const std::string name_;
  size_t cursor_index_;
  std::vector<std::unique_ptr<Resource>> collection_;
};

/** handles acquire/release resource by the consumer (RefVectorWithLeader type).
 */
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
