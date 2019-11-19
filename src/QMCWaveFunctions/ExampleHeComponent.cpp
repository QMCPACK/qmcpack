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


#include "QMCWaveFunctions/ExampleHeComponent.h"
#include "OhmmsData/AttributeSet.h"

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

ExampleHeComponent::LogValueType ExampleHeComponent::evaluateLog(ParticleSet& P,
                                                             ParticleSet::ParticleGradient_t& G,
                                                             ParticleSet::ParticleLaplacian_t& L)
{
  const auto& ee_table = P.getDistTable(my_table_ee_idx_);
  // Only the lower triangle is up-to-date after particle-by-particle moves
  double r12  = ee_table.Distances[1][0];
  auto rhat12 = ee_table.Displacements[1][0] / r12;

  const auto& ei_table = P.getDistTable(my_table_ei_idx_);

  // First index is electron, second index is ion
  double r1 = ei_table.Distances[0][0];
  double r2 = ei_table.Distances[1][0];

  auto rhat1 = ei_table.Displacements[0][0] / r1;
  auto rhat2 = ei_table.Displacements[1][0] / r2;

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

  return -Z * (r1 + r2) + std::log(norm * norm) - u;
}

ExampleHeComponent::ValueType ExampleHeComponent::ratio(ParticleSet& P, int iat)
{
  const int jat = (iat == 0 ? 1 : 0);

  const auto& ee_table = P.getDistTable(my_table_ee_idx_);

  // during p-by-p move iat-th row of Distances and Displacements[iat] are up-to-date
  // because setActive is called before ratio
  double r12_old = ee_table.Distances[iat][jat];
  double r12_new = ee_table.Temp_r[jat];

  const auto& ei_table = P.getDistTable(my_table_ei_idx_);

  double r_old = ei_table.Distances[iat][0];
  double r_new = ei_table.Temp_r[0];

  double u_old = A * r12_old / (B * r12_old + 1);
  double u_new = A * r12_new / (B * r12_new + 1);

  double log_v_old = -Z * (r_old)-u_old;
  double log_v_new = -Z * (r_new)-u_new;

  return std::exp((log_v_new - log_v_old));
}

ExampleHeComponent::GradType ExampleHeComponent::evalGrad(ParticleSet& P, int iat)
{
  const auto& ei_table = P.getDistTable(my_table_ei_idx_);

  double r  = ei_table.Distances[iat][0];
  auto rhat = ei_table.Displacements[iat][0] / r;

  const auto& ee_table = P.getDistTable(my_table_ee_idx_);

  const int jat = (iat == 0 ? 1 : 0);

  // during p-by-p move iat-th row of Distances and Displacements[iat] are up-to-date
  // because setActive is called before evalGrad
  double r12  = ee_table.Distances[iat][jat];
  auto rhat12 = ee_table.Displacements[iat][jat] / r12;

  double du = A / ((B * r12 + 1) * (B * r12 + 1));

  return Z * rhat + rhat12 * du;
}

ExampleHeComponent::ValueType ExampleHeComponent::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  const auto& ee_table = P.getDistTable(my_table_ee_idx_);

  const int jat = (iat == 0 ? 1 : 0);

  // during p-by-p move iat-th row of Distances and Displacements[iat] are up-to-date
  // because setActive is called before ratioGrad
  double r12_old = ee_table.Distances[iat][jat];
  double r12_new = ee_table.Temp_r[jat];

  auto rhat12 = ee_table.Temp_dr[jat] / r12_new;

  const auto& ei_table = P.getDistTable(my_table_ei_idx_);

  double r_old = ei_table.Distances[iat][0];
  double r_new = ei_table.Temp_r[0];

  auto rhat = ei_table.Temp_dr[0] / r_new;

  double du = A / ((B * r12_new + 1) * (B * r12_new + 1));
  double df = -Z;
  grad_iat  = -df * rhat + du * rhat12;

  double u_old = A * r12_old / (B * r12_old + 1);
  double u_new = A * r12_new / (B * r12_new + 1);

  double log_v_old = -Z * (r_old)-u_old;
  double log_v_new = -Z * (r_new)-u_new;

  return std::exp((log_v_new - log_v_old));
}


ExampleHeComponent::LogValueType ExampleHeComponent::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  return evaluateLog(P, P.G, P.L);
}

WaveFunctionComponentPtr ExampleHeComponent::makeClone(ParticleSet& tpq) const { return new ExampleHeComponent(*this); }

void ExampleHeComponent::resetParameters(const OptVariablesType& active)
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
                                             std::vector<ValueType>& dlogpsi,
                                             std::vector<ValueType>& dhpsioverpsi)
{
  typedef TinyVector<RealType, 3> RealGradType;


  double tmpB = std::real(optvars[0]);

  const auto& ee_table = P.getDistTable(my_table_ee_idx_);
  double r12           = ee_table.Distances[1][0];
  auto rhat12          = ee_table.Displacements[1][0] / r12;

  const auto& ei_table = P.getDistTable(my_table_ei_idx_);

  double r1 = ei_table.Distances[0][0];
  double r2 = ei_table.Distances[1][0];

  auto rhat1 = ei_table.Displacements[0][0] / r1;
  auto rhat2 = ei_table.Displacements[1][0] / r2;

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
