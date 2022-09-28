//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "ExampleHeComponent.h"
#include "OhmmsData/AttributeSet.h"
#include "DistanceTable.h"

/**@file ExampleHeComponent.cpp
 */
namespace qmcplusplus
{
bool ExampleHeComponent::put(xmlNodePtr cur)
{
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));

    if (cname == "var")
    {
      std::string id_in;
      std::string p_name;
      OhmmsAttributeSet rAttrib;
      rAttrib.add(id_in, "id");
      rAttrib.add(p_name, "name");
      rAttrib.put(cur);


      if (p_name == "B")
      {
        ID_B = id_in;
        putContent(B, cur);
        opt_B = true;
      }
    }
    cur = cur->next;
  }

  my_vars_.clear();

  if (opt_B)
    my_vars_.insert(ID_B, B, opt_B, optimize::OTHER_P);


  // Electron-nucleus cusp
  Z = 2.0;

  // Electron-electron cusp
  A = -0.5;
  //A = 0.0;


  return true;
}

ExampleHeComponent::LogValueType ExampleHeComponent::evaluateLog(const ParticleSet& P,
                                                                 ParticleSet::ParticleGradient& G,
                                                                 ParticleSet::ParticleLaplacian& L)
{
  const auto& ee_table  = P.getDistTableAA(my_table_ee_idx_);
  const auto& ee_dists  = ee_table.getDistances();
  const auto& ee_displs = ee_table.getDisplacements();
  // Only the lower triangle is up-to-date after particle-by-particle moves
  double r12  = ee_dists[1][0];
  auto rhat12 = ee_displs[1][0] / r12;

  const auto& ei_table  = P.getDistTableAB(my_table_ei_idx_);
  const auto& ei_dists  = ei_table.getDistances();
  const auto& ei_displs = ei_table.getDisplacements();

  // First index is electron, second index is ion
  double r1 = ei_dists[0][0];
  double r2 = ei_dists[1][0];

  auto rhat1 = ei_displs[0][0] / r1;
  auto rhat2 = ei_displs[1][0] / r2;

  // Normalization for STO is not strictly necessary in this example, but it
  // makes the values directly comparable to existing code
  double norm = 0.5 / std::sqrt(M_PI);

  double du = A / ((B * r12 + 1) * (B * r12 + 1));

  double df1 = -Z;
  G[0]       = -df1 * rhat1 - du * rhat12;
  // Previous line copies all three components and is the same as
  //G[0][0] = -df1*rhat1[0];
  //G[0][1] = -df1*rhat1[1];
  //G[0][2] = -df1*rhat1[2];

  double df2 = -Z;
  G[1]       = -df2 * rhat2 + du * rhat12;

  double del_u = 2 * A / (r12 * (B * r12 + 1) * (B * r12 + 1) * (B * r12 + 1));

  L[0] = -2 * Z / r1 - del_u;
  L[1] = -2 * Z / r2 - del_u;

  double u = A * r12 / (B * r12 + 1) - A / B;

  log_value_ = -Z * (r1 + r2) + std::log(norm * norm) - u;
  return log_value_;
}

ExampleHeComponent::PsiValueType ExampleHeComponent::ratio(ParticleSet& P, int iat)
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

  return std::exp(static_cast<PsiValueType>(log_v_new - log_v_old));
}

ExampleHeComponent::GradType ExampleHeComponent::evalGrad(ParticleSet& P, int iat)
{
  const auto& ei_table  = P.getDistTableAB(my_table_ei_idx_);
  const auto& ei_dists  = ei_table.getDistances();
  const auto& ei_displs = ei_table.getDisplacements();

  double r  = ei_dists[iat][0];
  auto rhat = ei_displs[iat][0] / r;

  const auto& ee_table  = P.getDistTableAA(my_table_ee_idx_);
  const auto& ee_dists  = ee_table.getDistances();
  const auto& ee_displs = ee_table.getDisplacements();

  // only the lower triangle of e-e Distances and Displacements can be used.
  double r12  = ee_dists[1][0];
  auto rhat12 = (iat == 0 ? -ee_displs[1][0] : ee_displs[1][0]) / r12;

  double du = A / ((B * r12 + 1) * (B * r12 + 1));

  return Z * rhat + rhat12 * du;
}

ExampleHeComponent::PsiValueType ExampleHeComponent::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  const auto& ee_table   = P.getDistTableAA(my_table_ee_idx_);
  const auto& ee_dists   = ee_table.getDistances();
  const auto& ee_displs  = ee_table.getDisplacements();
  const auto& ee_temp_r  = ee_table.getTempDists();
  const auto& ee_temp_dr = ee_table.getTempDispls();

  const int jat = (iat == 0 ? 1 : 0);
  // only the lower triangle of e-e Distances and Displacements can be used.
  double r12_old = ee_dists[1][0];
  double r12_new = ee_temp_r[jat];

  auto rhat12 = ee_temp_dr[jat] / r12_new;

  const auto& ei_table   = P.getDistTableAB(my_table_ei_idx_);
  const auto& ei_dists   = ei_table.getDistances();
  const auto& ei_displs  = ei_table.getDisplacements();
  const auto& ei_temp_r  = ei_table.getTempDists();
  const auto& ei_temp_dr = ei_table.getTempDispls();

  double r_old = ei_dists[iat][0];
  double r_new = ei_temp_r[0];

  auto rhat = ei_temp_dr[0] / r_new;

  double du = A / ((B * r12_new + 1) * (B * r12_new + 1));
  double df = -Z;
  grad_iat  = -df * rhat + du * rhat12;

  double u_old = A * r12_old / (B * r12_old + 1);
  double u_new = A * r12_new / (B * r12_new + 1);

  double log_v_old = -Z * (r_old)-u_old;
  double log_v_new = -Z * (r_new)-u_new;

  return std::exp(static_cast<PsiValueType>(log_v_new - log_v_old));
}


ExampleHeComponent::LogValueType ExampleHeComponent::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  return evaluateLog(P, P.G, P.L);
}

std::unique_ptr<WaveFunctionComponent> ExampleHeComponent::makeClone(ParticleSet& tpq) const
{
  return std::make_unique<ExampleHeComponent>(*this);
}

void ExampleHeComponent::resetParametersExclusive(const OptVariablesType& active)
{
  if (my_vars_.size())
  {
    int ia = my_vars_.where(0);
    if (ia > -1)
    {
      int i = 0;
      if (opt_B)
      {
        B = std::real(my_vars_[i++] = active[ia++]);
      }
    }
  }
}

void ExampleHeComponent::evaluateDerivatives(ParticleSet& P,
                                             const OptVariablesType& optvars,
                                             Vector<ValueType>& dlogpsi,
                                             Vector<ValueType>& dhpsioverpsi)
{
  using RealGradType = TinyVector<RealType, 3>;

  double tmpB = std::real(optvars[0]);

  const auto& ee_table   = P.getDistTableAA(my_table_ee_idx_);
  const auto& ee_dists   = ee_table.getDistances();
  const auto& ee_displs  = ee_table.getDisplacements();
  const auto& ee_temp_r  = ee_table.getTempDists();
  const auto& ee_temp_dr = ee_table.getTempDispls();

  double r12  = ee_dists[1][0];
  auto rhat12 = ee_displs[1][0] / r12;

  const auto& ei_table   = P.getDistTableAB(my_table_ei_idx_);
  const auto& ei_dists   = ei_table.getDistances();
  const auto& ei_displs  = ei_table.getDisplacements();
  const auto& ei_temp_r  = ei_table.getTempDists();
  const auto& ei_temp_dr = ei_table.getTempDispls();

  double r1 = ei_dists[0][0];
  double r2 = ei_dists[1][0];

  auto rhat1 = ei_displs[0][0] / r1;
  auto rhat2 = ei_displs[1][0] / r2;

  dlogpsi[0] = A * r12 * r12 / ((tmpB * r12 + 1) * (tmpB * r12 + 1)) - A / (tmpB * tmpB);

  double df       = -Z;
  double du       = A / ((B * r12 + 1) * (B * r12 + 1));
  RealGradType G1 = -df * rhat1 - rhat12 * du;
  RealGradType G2 = -df * rhat2 + rhat12 * du;

  double dudb      = -2 * A * r12 / std::pow(B * r12 + 1, 3);
  RealGradType dG1 = rhat12 * dudb;
  RealGradType dG2 = -1 * rhat12 * dudb;

  double dlap = -6 * A / std::pow((tmpB * r12 + 1), 4);
  dhpsioverpsi[0] =
      dlap + (G1[0] * dG1[0] + G1[1] * dG1[1] + G1[2] * dG1[2]) + (G2[0] * dG2[0] + G2[1] * dG2[1] + G2[2] * dG2[2]);
}


}; // namespace qmcplusplus
