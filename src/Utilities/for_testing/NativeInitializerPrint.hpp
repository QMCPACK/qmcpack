//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  The operator<< functions provide output for data structures that can be
 *  used to directly initialize them in test code.
 *
 *  Usage is present in future test code and provides the "coverage" of these
 *  functions.
 */

#ifndef QMCPLUSPLUS_NATIVE_INTITIALIZER_PRINT_HPP
#define QMCPLUSPLUS_NATIVE_INTITIALIZER_PRINT_HPP

#include <string>
#include <iostream>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

/** This wrapper is to allow us to leave the user facing operator<< for classes alone
 */
template<typename OBJECT>
class NativePrint
{
public:
  using type_t = OBJECT;
  NativePrint(const OBJECT& obj) : obj_(obj) {}
  const OBJECT& get_obj() const { return obj_; }

private:
  OBJECT obj_;
};

template<class T, unsigned D>
inline std::ostream& operator<<(std::ostream& out, const NativePrint<TinyVector<T, D>>& np_vec)
{
  out << "{ ";
  auto vec = np_vec.get_obj();
  for (int i = 0; i < D - 1; ++i)
    out << std::setw(12) << std::setprecision(10) << vec[i] << ", ";
  out << std::setw(12) << std::setprecision(10) << vec[D - 1] << " }";
  return out;
}

template<class T>
inline std::ostream& operator<<(std::ostream& out, const NativePrint<std::vector<T>>& np_vec)
{
  out << "{ ";
  auto vec = np_vec.get_obj();
  for (T& t : vec)
    out << std::setprecision(10) << t << ", ";
  out << " }";
  return out;
}

template<>
inline std::ostream& operator<<(std::ostream& out, const NativePrint<std::vector<bool>>& np_vec)
{
  out << "{ ";
  auto vec = np_vec.get_obj();
  for (const bool& t : vec)
  {
    std::string bool_str = t ? "true" : "false";
    out << std::setprecision(10) << bool_str << ", ";
  }
  out << " }";
  return out;
}

template<class T>
inline std::ostream& operator<<(std::ostream& out, const NativePrint<Vector<T>>& np_vec)
{
  out << "{ ";
  auto vec = np_vec.get_obj();
  for (T& t : vec)
    out << std::setprecision(10) << t << ", ";
  out << " }";
  return out;
}


} // namespace qmcplusplus

#endif
