//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TYPE_TESTS_HPP
#define QMCPLUSPLUS_TYPE_TESTS_HPP

namespace qmcplusplus
{
template<typename T>
struct implIsComplex : public std::false_type
{};
template<typename T>
struct implIsComplex<std::complex<T>> : public std::true_type
{};

template<typename T>
using IsComplex = std::enable_if_t<implIsComplex<T>::value, bool>;
template<typename T>
using IsReal = std::enable_if_t<std::is_floating_point<T>::value, bool>;

} // namespace qmcplusplus

#endif
