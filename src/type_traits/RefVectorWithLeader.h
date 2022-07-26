//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_REFVECTORWITHLEADER_H
#define QMCPLUSPLUS_REFVECTORWITHLEADER_H

#include <cassert>
#include <vector>
#include <memory>

namespace qmcplusplus
{
template<typename T>
class RefVectorWithLeader : public std::vector<std::reference_wrapper<T>>
{
public:
  RefVectorWithLeader(T& leader) : leader_(leader) {}

  RefVectorWithLeader(T& leader, const std::vector<std::reference_wrapper<T>>& vec) : leader_(leader)
  {
    for (T& element : vec)
      this->push_back(element);
  }

  T& getLeader() const { return leader_; }

  T& operator[](size_t i) const { return std::vector<std::reference_wrapper<T>>::operator[](i).get(); }

  template<typename CASTTYPE>
  CASTTYPE& getCastedLeader() const
  {
    static_assert(std::is_const<T>::value == std::is_const<CASTTYPE>::value, "Unmatched const type qualifier!");
    assert(dynamic_cast<CASTTYPE*>(&leader_.get()) != nullptr);
    return static_cast<CASTTYPE&>(leader_.get());
  }

  template<typename CASTTYPE>
  CASTTYPE& getCastedElement(size_t i) const
  {
    static_assert(std::is_const<T>::value == std::is_const<CASTTYPE>::value, "Unmatched const type qualifier!");
    assert(dynamic_cast<CASTTYPE*>(&(*this)[i]) != nullptr);
    return static_cast<CASTTYPE&>((*this)[i]);
  }

private:
  std::reference_wrapper<T> leader_;
};
} // namespace qmcplusplus

#endif
