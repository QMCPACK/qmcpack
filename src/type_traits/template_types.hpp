//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TEMPLATE_TYPES_HPP
#define QMCPLUSPLUS_TEMPLATE_TYPES_HPP

#include <vector>
#include <functional>
#include <memory>

namespace qmcplusplus
{
/** @name Convenience templated type aliases
 *
 *  when feasible do not use previous aliases in following
 *  that way one line can be scanned to understand a single alias
 *  see: UPtrVector
 *  @{
 */
template<typename T>
using RefVector = std::vector<std::reference_wrapper<T>>;

template<typename T>
using UPtr = std::unique_ptr<T>;

template<typename T>
using UPtrVector = std::vector<std::unique_ptr<T>>;
/** }@ */


// temporary helper function
template<class TR, class T>
static RefVector<TR> makeRefVector(std::vector<T>& vec_list)
{
  RefVector<TR> ref_list;
  ref_list.reserve(vec_list.size());
  for (int i = 0; i < vec_list.size(); ++i)
    ref_list.push_back(vec_list[i]);
  return ref_list;
}

// temporary helper function
template<class T>
static RefVector<T> convertUPtrToRefVector(const UPtrVector<T>& ptr_list)
{
  RefVector<T> ref_list;
  ref_list.reserve(ptr_list.size());
  for (const UPtr<T>& ptr : ptr_list)
    ref_list.push_back(*ptr);
  return ref_list;
}

// temporary helper function
template<class T>
static RefVector<T> convertPtrToRefVector(const std::vector<T*>& ptr_list)
{
  RefVector<T> ref_list;
  ref_list.reserve(ptr_list.size());
  for (auto ptr : ptr_list)
    ref_list.push_back(*ptr);
  return ref_list;
}

// temporary helper function
template<class T>
static std::vector<T*> convert_ref_to_ptr_list(const std::vector<std::reference_wrapper<T>>& ref_list)
{
  std::vector<T*> ptr_list;
  ptr_list.reserve(ref_list.size());
  for (auto& ref : ref_list)
    ptr_list.push_back(&ref.get());
  return ptr_list;
}

// temporary helper function
template<class T>
static std::vector<T*> convertUPtrToPtrVector(const UPtrVector<T>& uptr_list)
{
  std::vector<T*> ptr_list;
  ptr_list.reserve(uptr_list.size());
  for (auto& uptr : uptr_list)
    ptr_list.push_back(uptr.get());
  return ptr_list;
}

} // namespace qmcplusplus
#endif
