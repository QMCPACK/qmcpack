//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OPTIMIZABLEOBJECT_H
#define QMCPLUSPLUS_OPTIMIZABLEOBJECT_H

#include "Configuration.h"
#include "OptimizableObjectT.h"

/**@file OptimizableObject.h
 *@brief Declaration of OptimizableObject
 */
namespace qmcplusplus
{
using opt_variables_type = OptVariablesTypeT<QMCTraits::ValueType>;

using OptVariablesType = OptVariablesTypeT<QMCTraits::ValueType>;

using OptimizableObject = OptimizableObjectT<QMCTraits::ValueType>;

using UniqueOptObjRefs = UniqueOptObjRefsT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
