//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file QuantumVariablesBuilder.h
 */
#ifndef QMCPLUSPLUS_QUANTUM_VARIABLE_BUILDER_H
#define QMCPLUSPLUS_QUANTUM_VARIABLE_BUILDER_H

#include "Particle/QuantumVariables.h"

namespace qmcplusplus
{
/** create QuantumVariables based on kind
 */
std::unique_ptr<QuantumVariables> createQuantumVariables(const QuantumVariableKind kind = QuantumVariableKind::QV_POS);
} // namespace qmcplusplus
#endif
