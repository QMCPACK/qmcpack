
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


#ifndef QMCPLUSPLUS_MAKE_ESTIMATOR_MANAGER_H
#define QMCPLUSPLUS_MAKE_ESTIMATOR_MANAGER_H

#include "EstimatorManagerNew.h"
#include <optional>
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
UPtr<EstimatorManagerNew> makeEstimatorManager(const std::optional<EstimatorManagerInput>& global_emi,
                                               const std::optional<EstimatorManagerInput>& driver_emi,
                                               ParticleSet& pset_primary,
                                               TrialWaveFunction& twf_primary,
                                               QMCHamiltonian& ham_primary,
                                               const EstimatorManagerNew::PSPool& pset_pool,
                                               Communicate* comm);
} // namespace qmcplusplus
#endif
