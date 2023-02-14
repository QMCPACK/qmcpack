//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "BareForce.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
BareForce::BareForce(ParticleSet& ions, ParticleSet& elns) : ForceBase(ions, elns), d_ei_id_(elns.addTable(ions))
{
  name_   = "HF_Force_Base";
  prefix_ = "HFBase";
}

std::string BareForce::getClassName() const { return "BareForce"; }

void BareForce::resetTargetParticleSet(ParticleSet& P) {}

std::unique_ptr<OperatorBase> BareForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<BareForce>(*this);
}

void BareForce::registerObservables(std::vector<ObservableHelper>& h5list, hdf_archive& file) const
{
  registerObservablesF(h5list, file);
}

void BareForce::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addObservablesF(plist);
  my_index_ = first_force_index_;
}

void BareForce::setObservables(PropertySetType& plist) { setObservablesF(plist); }

void BareForce::setParticlePropertyList(PropertySetType& plist, int offset) { setParticleSetF(plist, offset); }

BareForce::Return_t BareForce::evaluate(ParticleSet& P)
{
  forces_                                   = forces_ion_ion_;
  const auto& d_ab                          = P.getDistTableAB(d_ei_id_);
  const ParticleSet::Scalar_t* restrict Zat = ions_.Z.first_address();
  const ParticleSet::Scalar_t* restrict Qat = P.Z.first_address();
  //Loop over distinct eln-ion pairs
  for (int jat = 0; jat < d_ab.targets(); jat++)
  {
    const auto& ab_dist  = d_ab.getDistRow(jat);
    const auto& ab_displ = d_ab.getDisplRow(jat);
    for (int iat = 0; iat < d_ab.sources(); iat++)
    {
      Real rinv = 1.0 / ab_dist[iat];
      Real r3zz = Qat[jat] * Zat[iat] * rinv * rinv * rinv;
      forces_[iat] += r3zz * ab_displ[iat];
    }
  }
  tries_++;
  return 0.0;
}

bool BareForce::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix_, "name");
  attr.add(ionionforce, "add_ion_ion_");
  attr.put(cur);
  add_ion_ion_ = (ionionforce == "yes" || ionionforce == "true");
  return true;
}

bool BareForce::get(std::ostream& os) const
{
  os << "Force Base Hamiltonian: " << pair_name_;
  return true;
}

} // namespace qmcplusplus
