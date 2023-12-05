//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  Compilation units that construct QMCDriverInput need visibility to the actual input classe
 *  types in the delegation tree.
 */
#ifndef QMCPLUSPLUS_QMCDRIVERINPUTDELEGATES_H
#define QMCPLUSPLUS_QMCDRIVERINPUTDELEGATES_H

#include "EstimatorManagerInput.h"
#include "ScalarEstimatorInputs.h"
#include "MomentumDistributionInput.h"
#include "OneBodyDensityMatricesInput.h"
#include "SpinDensityInput.h"
#include "MagnetizationDensityInput.h"
#include "PerParticleHamiltonianLoggerInput.h"

#endif
