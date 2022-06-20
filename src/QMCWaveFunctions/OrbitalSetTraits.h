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
#include "type_traits/complex_help.hpp"
#include "VariableSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{
/** dummy class for templated classes
 */
struct DummyGrid
{
  inline void locate(double r) {}
  DummyGrid* makeClone() const { return new DummyGrid; }
};

using QuantumNumberType = TinyVector<int, 4>;

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
  using RealType       = RealAlias<T>;
  using ValueType      = T;
  using IndexType      = int;
  using PosType        = TinyVector<RealType, DIM>;
  using GradType       = TinyVector<ValueType, DIM>;
  using HessType       = Tensor<ValueType, DIM>;
  using TensorType     = Tensor<ValueType, DIM>;
  using GradHessType   = TinyVector<Tensor<ValueType, DIM>, DIM>;
  using IndexVector    = Vector<IndexType>;
  using ValueVector    = Vector<ValueType>;
  using ValueMatrix    = Matrix<ValueType>;
  using GradVector     = Vector<GradType>;
  using GradMatrix     = Matrix<GradType>;
  using HessVector     = Vector<HessType>;
  using HessMatrix     = Matrix<HessType>;
  using GradHessVector = Vector<GradHessType>;
  using GradHessMatrix = Matrix<GradHessType>;
  using VGLVector      = VectorSoaContainer<ValueType, DIM + 2>;
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

  inline static std::complex<T> convert(const std::complex<T>& logpsi) { return std::exp(logpsi); }
};

} // namespace qmcplusplus

#endif
