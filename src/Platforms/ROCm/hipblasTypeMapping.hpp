//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HIPBLAS_TYPE_MAPPING_HPP
#define QMCPLUSPLUS_HIPBLAS_TYPE_MAPPING_HPP

#include <type_traits>
#include "type_traits/type_mapping.hpp"
#include <hipblas/hipblas.h>

namespace qmcplusplus
{

// This saves us writing specific overloads with reinterpret casts for different std::complex to hipblasComplex types.
template<typename T>
using hipblasTypeMap =
    typename std::disjunction<OnTypesEqual<T, float, float>,
                              OnTypesEqual<T, double, double>,
                              OnTypesEqual<T, float*, float*>,
                              OnTypesEqual<T, double*, double*>,
                              OnTypesEqual<T, float**, float**>,
                              OnTypesEqual<T, double**, double**>,
                              OnTypesEqual<T, std::complex<double>, hipblasDoubleComplex>,
                              OnTypesEqual<T, std::complex<float>, hipblasComplex>,
                              OnTypesEqual<T, std::complex<double>*, hipblasDoubleComplex*>,
                              OnTypesEqual<T, std::complex<float>**, hipblasComplex**>,
                              OnTypesEqual<T, std::complex<double>**, hipblasDoubleComplex**>,
                              OnTypesEqual<T, std::complex<float>*, hipblasComplex*>,
                              OnTypesEqual<T, const std::complex<double>*, const hipblasDoubleComplex*>,
                              OnTypesEqual<T, const std::complex<float>*, const hipblasComplex*>,
                              OnTypesEqual<T, const std::complex<float>**, const hipblasComplex**>,
                              OnTypesEqual<T, const std::complex<double>**, const hipblasDoubleComplex**>,
                              OnTypesEqual<T, const std::complex<float>* const*, const hipblasComplex* const*>,
                              OnTypesEqual<T, const std::complex<double>* const*, const hipblasDoubleComplex* const*>,
                              default_type<void>>::type;

template<typename T>
hipblasTypeMap<T> casthipblasType(T var)
{
  return reinterpret_cast<hipblasTypeMap<T>>(var);
}

} // namespace qmcplusplus

#endif // QMCPLUSPLUS_HIPBLAS_TYPE_MAPPING_HPP
