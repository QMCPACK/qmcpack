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

#include "ExampleHeComponent.h"

namespace qmcplusplus
{

template<class T>
T ExampleHeComponent::do_ratioT(ParticleSet& P, int iat)
{
  const auto& ee_table  = P.getDistTableAA(my_table_ee_idx_);
  const auto& ee_dists  = ee_table.getDistances();
  const auto& ee_temp_r = ee_table.getTempDists();

  // only the lower triangle of e-e Distances and Displacements can be used.
  double r12_old = ee_dists[1][0];
  double r12_new = ee_temp_r[iat == 0 ? 1 : 0];

  const auto& ei_table  = P.getDistTableAB(my_table_ei_idx_);
  const auto& ei_dists  = ei_table.getDistances();
  const auto& ei_temp_r = ei_table.getTempDists();

  double r_old = ei_dists[iat][0];
  double r_new = ei_temp_r[0];

  double u_old = A * r12_old / (B * r12_old + 1);
  double u_new = A * r12_new / (B * r12_new + 1);

  double log_v_old = -Z * (r_old)-u_old;
  double log_v_new = -Z * (r_new)-u_new;

  return std::exp(static_cast<T>(log_v_new - log_v_old));
}


} // namespace qmcplusplus