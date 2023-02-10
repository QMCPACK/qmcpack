//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ForceCeperley.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
ForceCeperley::ForceCeperley(ParticleSet& ions, ParticleSet& elns)
    : ForceBase(ions, elns), d_aa_ID(ions.addTable(ions)), d_ei_ID(elns.addTable(ions))
{
  ReportEngine PRE("ForceCeperley", "ForceCeperley");
  name_  = "Ceperley_Force_Base";
  prefix_ = "HFCep";
  // Defaults
  Rcut    = 0.4;
  m_exp   = 2;
  N_basis = 4;
  forces_ = 0.0;
  ///////////////////////////////////////////////////////////////
  ions.update();
  evaluate_IonIon(forces_ion_ion_);
}

void ForceCeperley::evaluate_IonIon(ParticleSet::ParticlePos& forces) const
{
  forces = 0.0;
  const auto& d_aa(ions_.getDistTableAA(d_aa_ID));
  const ParticleScalar* restrict Zat = ions_.Z.first_address();
  for (size_t ipart = 1; ipart < n_nuc_; ipart++)
  {
    const auto& dist  = d_aa.getDistRow(ipart);
    const auto& displ = d_aa.getDisplRow(ipart);
    for (size_t jpart = 0; jpart < ipart; ++jpart)
    {
      RealType rinv = 1.0 / dist[jpart];
      RealType r3zz = Zat[jpart] * Zat[ipart] * rinv * rinv * rinv;
      forces[jpart] += r3zz * displ[jpart];
      forces[ipart] -= r3zz * displ[jpart];
    }
  }
}

void ForceCeperley::InitMatrix()
{
  Sinv.resize(N_basis, N_basis);
  h.resize(N_basis);
  c.resize(N_basis);
  for (int k = 0; k < N_basis; k++)
  {
    h[k] = std::pow(Rcut, (k + 2)) / static_cast<RealType>(k + 2);
    for (int j = 0; j < N_basis; j++)
    {
      Sinv(k, j) = std::pow(Rcut, (m_exp + k + j + 3)) / static_cast<RealType>(m_exp + k + j + 3);
    }
  }
  // in Numerics/DeterminantOperators.h
  invert_matrix(Sinv, false);
  // in Numerics/MatrixOperators.h
  MatrixOperators::product(Sinv, h, c);
}

ForceCeperley::Return_t ForceCeperley::evaluate(ParticleSet& P)
{
  if (add_ion_ion_ == true)
    forces_ = forces_ion_ion_;
  else
    forces_ = 0.0;
  const auto& d_ab                   = P.getDistTableAB(d_ei_ID);
  const ParticleScalar* restrict Zat = ions_.Z.first_address();
  const ParticleScalar* restrict Qat = P.Z.first_address();
  for (int jat = 0; jat < n_el_; jat++)
  {
    const auto& dist  = d_ab.getDistRow(jat);
    const auto& displ = d_ab.getDisplRow(jat);
    for (int iat = 0; iat < n_nuc_; iat++)
    {
      // electron contribution (special treatment if distance is inside cutoff!)
      RealType r       = dist[iat];
      RealType zoverr3 = Qat[jat] * Zat[iat] / (r * r * r);
      if (r >= Rcut)
      {
        forces_[iat] += zoverr3 * displ[iat];
      }
      else
      {
        RealType g_q = 0.0;
        for (int q = 0; q < N_basis; q++)
        {
          g_q += c[q] * std::pow(r, m_exp + q + 1);
        }
        g_q *= zoverr3;
        // negative sign accounts for definition of target as electrons
        forces_[iat] += g_q * displ[iat];
      }
    }
  }
  return 0.0;
}

bool ForceCeperley::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix_, "name");
  attr.add(ionionforce, "add_ion_ion_");
  attr.put(cur);
  add_ion_ion_ = (ionionforce == "yes") || (ionionforce == "true");
  app_log() << "ionionforce = " << ionionforce << std::endl;
  app_log() << "add_ion_ion_=" << add_ion_ion_ << std::endl;
  app_log() << "first_time_= " << first_time_ << std::endl;
  ParameterSet fcep_param_set;
  fcep_param_set.add(Rcut, "rcut");
  fcep_param_set.add(N_basis, "nbasis");
  fcep_param_set.add(m_exp, "weight_exp");
  fcep_param_set.put(cur);
  app_log() << "    ForceCeperley Parameters" << std::endl;
  app_log() << "        ForceCeperley::Rcut=" << Rcut << std::endl;
  app_log() << "        ForceCeperley::N_basis=" << N_basis << std::endl;
  app_log() << "        ForceCeperley::m_exp=" << m_exp << std::endl;
  InitMatrix();
  return true;
}

std::unique_ptr<OperatorBase> ForceCeperley::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<ForceCeperley>(*this);
}
} // namespace qmcplusplus
