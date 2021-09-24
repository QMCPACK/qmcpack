//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file NRCOptimizationFunctionWrapper.h
 * @brief Simplifies using NRCOptimization
 */
#ifndef QMCPLUSPLUS_NRCOPTIMIZATIONFUNCTIONWRAPPER_H
#define QMCPLUSPLUS_NRCOPTIMIZATIONFUNCTIONWRAPPER_H

#include "Optimize/NRCOptimization.h"

namespace qmcplusplus
{

// Wrapper class for evaluation of Func.
// The class that implements the objective function no longer needs to inherit NRCOptimization.
// It can instead use this class as a member object.
template<class T>
class NRCOptimizationFunctionWrapper : public NRCOptimization<QMCTraits::RealType>
{
public:
  T& object;
  NRCOptimizationFunctionWrapper(T& o) : object(o) {}
  Return_t Func(Return_t dl) override { return object.costFunc(dl); }
};
} // namespace qmcplusplus
#endif
