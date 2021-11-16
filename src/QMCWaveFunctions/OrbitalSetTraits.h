//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BasisSetBase.h
 * @brief Declaration of a base class of BasisSet
 */
#ifndef QMCPLUSPLUS_ORBITALSETTRAITS_H
#define QMCPLUSPLUS_ORBITALSETTRAITS_H

#include "Configuration.h"
#include "type_traits/scalar_traits.h"
#include "VariableSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"

namespace qmcplusplus
{
/** dummy class for templated classes
 */
struct DummyGrid
{
  inline void locate(double r) {}
  DummyGrid* makeClone() const { return new DummyGrid; }
};

typedef TinyVector<int, 4> QuantumNumberType;

enum
{
  q_n = 0,
  q_l,
  q_m,
  q_s
};

/** trait class to handel a set of Orbitals
 */
template<typename T>
struct OrbitalSetTraits //: public OrbitalTraits<T>
{
  enum
  {
    DIM = OHMMS_DIM
  };
  typedef typename scalar_traits<T>::real_type RealType;
  typedef typename scalar_traits<T>::value_type ValueType;
  typedef int IndexType;
  typedef TinyVector<RealType, DIM> PosType;
  typedef TinyVector<ValueType, DIM> GradType;
  typedef Tensor<ValueType, DIM> HessType;
  typedef Tensor<ValueType, DIM> TensorType;
  typedef TinyVector<Tensor<ValueType, DIM>, DIM> GradHessType;
  typedef Vector<IndexType> IndexVector_t;
  typedef Vector<ValueType> ValueVector_t;
  typedef Matrix<ValueType> ValueMatrix_t;
  typedef Vector<GradType> GradVector_t;
  typedef Matrix<GradType> GradMatrix_t;
  typedef Vector<HessType> HessVector_t;
  typedef Matrix<HessType> HessMatrix_t;
  typedef Vector<GradHessType> GradHessVector_t;
  typedef Matrix<GradHessType> GradHessMatrix_t;
  typedef VectorSoaContainer<ValueType, DIM + 2> VGLVector_t;
};

///typedef for a set of variables that are varied during an optimization
using opt_variables_type = optimize::VariableSet;
///typedef for a set of variables that can be varied
using variable_map_type = optimize::VariableSet::variable_map_type;

/** evaluate log(psi) as log(|psi|) and phase
 * @param psi real/complex value
 * @return complex<T>(log(|psi|), arg(psi))
 *
 * The return value is always complex regardless of the type of psi.
 * The sign of of a real psi value is represented in the complex part of the return value.
 * The phase of std::log(complex) is in range [-pi, pi] defined in C++
 */
template<typename T>
inline std::complex<T> convertValueToLog(const std::complex<T>& logpsi)
{
  return std::log(logpsi);
}

template<typename T>
inline std::complex<T> convertValueToLog(const T logpsi)
{
  return std::log(std::complex<T>(logpsi));
}

/** evaluate psi based on log(psi)
 * @param logpsi complex value
 * @return exp(log(psi))
 *
 * LogToValue<ValueType>::convert(complex) is the reverse operation of convertValueToLog(ValueType)
 */
template<typename T>
struct LogToValue
{
  template<typename T1>
  inline static T convert(const std::complex<T1>& logpsi)
  {
    return std::real(std::exp(logpsi));
  }
};

template<typename T>
struct LogToValue<std::complex<T>>
{
  template<typename T1, typename = std::enable_if_t<!std::is_same<T, T1>::value>>
  inline static std::complex<T> convert(const std::complex<T1>& logpsi)
  {
    std::complex<T> tmp(std::real(logpsi), std::imag(logpsi));
    return std::exp(tmp);
  }

  inline static std::complex<T> convert(const std::complex<T>& logpsi)
  {
    return std::exp(logpsi);
  }
};

} // namespace qmcplusplus

#endif
