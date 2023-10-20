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
  using Base = std::optional<std::reference_wrapper<RS>>;

  ResourceHandle() = default;
  ResourceHandle(RS& res) { Base::emplace(res); }

  bool hasResource() const { return Base::has_value(); }
  operator bool() const { return Base::has_value(); }

  RS& getResource() { return Base::value(); }
  const RS& getResource() const { return Base::value(); }

  operator RS&() { return this->value(); }
  operator const RS&() const { return Base::value(); }

  RS& release()
  {
    RS& res = Base::value();
    Base::reset();
    return res;
  }
};

} // namespace qmcplusplus
#endif
