//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "WaveFunctionComponent.h"

namespace qmcplusplus
{

template<class T>
T WaveFunctionComponent::ratio(ParticleSet& P, int iat)
{
  T value;
  do_ratio(value, P, iat);
  return value;
}

// Register supported types
#define declare_type(T) template T WaveFunctionComponent::ratio<T>(ParticleSet & P, int iat);
QMC_FOREACH_TYPE(declare_type)
#undef declare_type


} // namespace qmcplusplus