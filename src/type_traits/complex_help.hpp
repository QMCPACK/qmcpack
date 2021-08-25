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

#ifndef QMCPLUSPLUS_COMPLEX_HELP_HPP
#define QMCPLUSPLUS_COMPLEX_HELP_HPP

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

template <typename T, typename = bool>
struct RealAlias_impl {};

template <typename T>
struct RealAlias_impl<T, IsReal<T>> { using value_type = T; };

template <typename T>
struct RealAlias_impl<T, IsComplex<T>> { using value_type = typename T::value_type; };

/** If you have a function templated on a value that can be real or complex
 *   and you need to get the base Real type if its complex or just the real.
 *
 *  If you try to do this on anything but a fp or a std::complex<fp> you will
 *  get a compilation error.
 */
template <typename T>
using RealAlias = typename RealAlias_impl<T>::value_type;

} // namespace qmcplusplus

#endif
