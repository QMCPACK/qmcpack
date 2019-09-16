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

namespace qmcplusplus
{

template<typename T>
using RefVector = std::vector<std::reference_wrapper<T>>;

template<typename T>
using UPtrVector = std::vector<std::unique_ptr<T>>;

// temporal helper function
template<class T>
static std::vector<T*> convert_ref_to_ptr_list(const std::vector<std::reference_wrapper<T>>& ref_list)
{
  std::vector<T*> ptr_list;
  ptr_list.reserve(ref_list.size());
  for (auto& ref : ref_list)
    ptr_list.push_back(&ref.get());
  return ptr_list;
}

// temporal helper function
template<class T>
static std::vector<T*> convertUPtrToPtrVector(const UPtrVector<T>& uptr_list)
{
  std::vector<T*> ptr_list;
  ptr_list.reserve(uptr_list.size());
  for (auto& uptr : uptr_list)
    ptr_list.push_back(uptr.get());
  return ptr_list;
}

}
#endif
