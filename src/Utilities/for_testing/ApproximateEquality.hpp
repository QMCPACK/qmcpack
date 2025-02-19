//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_APPROXEQUALITY_HPP
#define QMCPLUSPLUS_APPROXEQUALITY_HPP

#include "catch.hpp"

#include <optional>
#include "type_traits/complex_help.hpp"
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{

template<typename T, IsComplex<T> = true>
bool approxEquality(T val_a, T val_b, std::optional<double> eps)
{
  if (eps)
    return val_a == ComplexApprox(val_b).epsilon(eps.value());
  else
    return val_a == ComplexApprox(val_b);
}

template<typename T, IsReal<T> = true>
bool approxEquality(T val_a, T val_b, std::optional<double> eps)
{
  if (eps)
    return val_a == Approx(val_b).epsilon(eps.value());
  else
    return val_a == Approx(val_b);
}
  
extern template bool approxEquality<float>(float val_a, float val_b, std::optional<double> eps);
extern template bool approxEquality<std::complex<float>>(std::complex<float> val_a,
                                                         std::complex<float> val_b,
                                                         std::optional<double> eps);
extern template bool approxEquality<double>(double val_a, double val_b, std::optional<double> eps);
extern template bool approxEquality<std::complex<double>>(std::complex<double> val_a,
                                                          std::complex<double> val_b,
                                                          std::optional<double> eps);

} // namespace qmcplusplus
#endif
