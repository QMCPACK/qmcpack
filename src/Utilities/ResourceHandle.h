//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RESOURCHANDLE_H
#define QMCPLUSPLUS_RESOURCHANDLE_H

#include <memory>

namespace qmcplusplus
{
/** ResourceHandle manages the temporary owned resource obtained from a collection
 *
 * Resource copy should be handled at the collection. Individual handle cannot be copied.
 * However, directly using unique_ptr prevents compiler generated copy constructor.
 * For this reason, the ResourceHandle class adds a copy constructor as no-op to facilitate
 * compiler generated copy constructor in classes owning ResourceHandle objects.
 */
template<class RS>
class ResourceHandle : public std::unique_ptr<RS>
{
public:
  ResourceHandle() = default;
  ResourceHandle(const ResourceHandle&) {}
  ResourceHandle(ResourceHandle&&) = default;
  ResourceHandle& operator=(const ResourceHandle&) {};
  ResourceHandle& operator=(ResourceHandle&&) = default;
  using std::unique_ptr<RS>::operator=;
};

} // namespace qmcplusplus
#endif
