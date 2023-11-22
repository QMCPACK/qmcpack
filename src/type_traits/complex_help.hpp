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

#include <type_traits>
#include <complex>

namespace qmcplusplus
{
template<typename T>
struct IsComplex_t : public std::false_type
{};
template<typename T>
struct IsComplex_t<std::complex<T>> : public std::true_type
{};

template<typename T>
using IsComplex = std::enable_if_t<IsComplex_t<T>::value, bool>;
template<typename T>
using IsReal = std::enable_if_t<std::is_floating_point<T>::value, bool>;

template<typename T, typename = bool>
struct RealAlias_impl
{};

template<typename T>
struct RealAlias_impl<T, IsReal<T>>
{
  using value_type = T;
};

template<typename T>
struct RealAlias_impl<T, IsComplex<T>>
{
  using value_type = typename T::value_type;
};

/** If you have a function templated on a value that can be real or complex
 *   and you need to get the base Real type if its complex or just the real.
 *
 *  If you try to do this on anything but a fp or a std::complex<fp> you will
 *  get a compilation error.
 */
template<typename T>
using RealAlias = typename RealAlias_impl<T>::value_type;

template<typename TREAL, typename TREF, typename = bool>
struct ValueAlias_impl
{};

template<typename TREAL, typename TREF>
struct ValueAlias_impl<TREAL, TREF, IsReal<TREF>>
{
  using value_type = TREAL;
};

template<typename TREAL, typename TREF>
struct ValueAlias_impl<TREAL, TREF, IsComplex<TREF>>
{
  using value_type = std::complex<TREAL>;
};

/** If you need to make a value type of a given precision based on a reference value type
 *  set the desired POD float point type as TREAL and set the reference type as TREF.
 *  If TREF is real/complex, the generated Value type is real/complex.
 */
template<typename TREAL, typename TREF, typename = std::enable_if_t<std::is_floating_point<TREAL>::value>>
using ValueAlias = typename ValueAlias_impl<TREAL, TREF>::value_type;

///real part of a scalar. Cannot be replaced by std::real due to AFQMC specific needs.
inline float real(const float& c) { return c; }
inline double real(const double& c) { return c; }
inline float real(const std::complex<float>& c) { return c.real(); }
inline double real(const std::complex<double>& c) { return c.real(); }
///imaginary part of a scalar. Cannot be replaced by std::imag due to AFQMC specific needs.
inline float imag(const float& c) { return 0; }
inline double imag(const double& c) { return 0; }
inline float imag(const std::complex<float>& c) { return c.imag(); }
inline double imag(const std::complex<double>& c) { return c.imag(); }
///Workaround to allow conj on scalar to return real instead of complex
inline float conj(const float& c) { return c; }
inline double conj(const double& c) { return c; }
inline std::complex<float> conj(const std::complex<float>& c) { return std::conj(c); }
inline std::complex<double> conj(const std::complex<double>& c) { return std::conj(c); }
//These copy complex->complex, real->real as is.  In event of complex->real, the imaginary part is ignored.
inline void copy_with_complex_cast(const std::complex<double>& source, std::complex<double>& dest) { dest = source; }
inline void copy_with_complex_cast(const std::complex<double>& source, double& dest) { dest = source.real(); }
inline void copy_with_complex_cast(const std::complex<float>& source, std::complex<float>& dest) { dest = source; }
inline void copy_with_complex_cast(const std::complex<float>& source, float& dest) { dest = source.real(); }

} // namespace qmcplusplus

#endif
