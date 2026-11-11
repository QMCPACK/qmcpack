//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_REFVECTORWITHRESOURCE_H
#define QMCPLUSPLUS_REFVECTORWITHRESOURCE_H

#include <cassert>
#include <vector>
#include <memory>
#include <type_traits>
#include "Resource.h"

namespace qmcplusplus
{
template<typename T>
class RefVectorWithResource : public std::vector<std::reference_wrapper<T>>
{
public:
  using BaseVec = std::vector<std::reference_wrapper<T>>;

  RefVectorWithResource(Resource& res) : res_(res) {}

  RefVectorWithResource(Resource& res, const BaseVec& vec) : BaseVec(vec), res_(res) {}

  RefVectorWithResource(Resource& res, BaseVec&& vec) : BaseVec(std::move(vec)), res_(res) {}

  Resource& getResource() const { return res_; }

  T& operator[](size_t i) const { return BaseVec::operator[](i).get(); }

  template<typename CASTTYPE>
  CASTTYPE& getCastedResource() const
  {
    static_assert(std::is_base_of<Resource, std::remove_cv_t<CASTTYPE>>::value,
                  "CASTTYPE must be a derived type of Resource!");
    assert(dynamic_cast<CASTTYPE*>(&res_.get()) != nullptr);
    return static_cast<CASTTYPE&>(res_.get());
  }

  template<typename CASTTYPE>
  CASTTYPE& getCastedElement(size_t i) const
  {
    static_assert(std::is_const<T>::value == std::is_const<CASTTYPE>::value, "Unmatched const type qualifier!");
    assert(dynamic_cast<CASTTYPE*>(&(*this)[i]) != nullptr);
    return static_cast<CASTTYPE&>((*this)[i]);
  }

private:
  std::reference_wrapper<Resource> res_;
};

template<typename T>
RefVectorWithResource(Resource&, const std::vector<std::reference_wrapper<T>>&) -> RefVectorWithResource<T>;
template<typename T>
RefVectorWithResource(Resource&, std::vector<std::reference_wrapper<T>>&&) -> RefVectorWithResource<T>;

} // namespace qmcplusplus

#endif
