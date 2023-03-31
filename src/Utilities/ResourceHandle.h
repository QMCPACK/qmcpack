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

#include <functional>
#include <optional>

namespace qmcplusplus
{
/** ResourceHandle manages the temporary resource referenced from a collection
 */
template<class RS>
class ResourceHandle : private std::optional<std::reference_wrapper<RS>>
{
public:
  ResourceHandle() = default;
  ResourceHandle(RS& res) { this->emplace(res); }
  ResourceHandle& operator=(RS& res)
  {
    this->emplace(res);
    return *this;
  }
  bool hasResource() const { return this->has_value(); }

  RS& getResource() { return this->value(); }
  const RS& getResource() const { return this->value(); }

  operator RS&() { return this->value(); }
  operator const RS&() const { return this->value(); }

  RS& release()
  {
    RS& res = this->value();
    this->reset();
    return res;
  }
};

} // namespace qmcplusplus
#endif
