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
#include <cassert>
#include "RefVectorWithLeader.h"

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

/** helper function to take vector of class A to refvector of any valid reference type for A
 *
 *  intended usage looks like this
 *  std::vector<DerivedType> vdt
 *  auto refvecbase = makeRefVector<BaseType>(vdt)
 *  or if you just want a refvector of type vdt
 *  auto refvec = makeRefVector<decltype(vdt)::value_type>(vdt)
 *
 *  godbolt.org indicates at least with clang 11&12 we get RVO here.
 *  auto ref_whatevers = makeRefVector<ValidTypeForReference>(whatevers);
 *  makes no extra copy.
 */
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
template<class T2, class T>
static RefVector<T2> convertUPtrToRefVector(const UPtrVector<T>& ptr_list)
{
  RefVector<T2> ref_list;
  ref_list.reserve(ptr_list.size());
  for (const UPtr<T>& ptr : ptr_list)
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

// handling subset
template<class T>
static RefVector<T> convertUPtrToRefVectorSubset(const UPtrVector<T>& ptr_list, int offset, int len)
{
  // check lower and upper bounds
  assert(offset >= 0);
  assert(offset + len <= ptr_list.size());

  RefVector<T> ref_list;
  ref_list.reserve(len);
  for (int i = offset; i < offset + len; i++)
    ref_list.push_back(*ptr_list[i]);

  return ref_list;
}

template<class T>
static RefVector<T> convertPtrToRefVectorSubset(const std::vector<T*>& ptr_list, int offset, int len)
{
  // check lower and upper bounds
  assert(offset >= 0);
  assert(offset + len <= ptr_list.size());

  RefVector<T> ref_list;
  ref_list.reserve(len);
  for (int i = offset; i < offset + len; i++)
    ref_list.push_back(*ptr_list[i]);

  return ref_list;
}

template<class T>
static RefVector<T> convertRefVectorToRefVectorSubset(const RefVector<T>& ref_list, int offset, int len)
{
  // check lower and upper bounds
  assert(offset >= 0);
  assert(offset + len <= ref_list.size());

  RefVector<T> sub_ref_list;
  sub_ref_list.reserve(len);
  for (int i = offset; i < offset + len; i++)
    sub_ref_list.push_back(ref_list[i]);

  return sub_ref_list;
}

} // namespace qmcplusplus
#endif
