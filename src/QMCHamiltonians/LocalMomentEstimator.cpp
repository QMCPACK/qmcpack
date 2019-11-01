//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <QMCHamiltonians/LocalMomentEstimator.h>
#include <Particle/DistanceTableData.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <set>

namespace qmcplusplus
{
LocalMomentEstimator::LocalMomentEstimator(ParticleSet& elns, ParticleSet& srcs)
  : ions(srcs), d_table_ID(elns.addTable(srcs, DT_AOS))
{
  int num_species = elns.groups();
  const SpeciesSet& e_species(elns.getSpeciesSet());
  int ne  = e_species.size();
  int neg = elns.getTotalNum();
  el_id.resize(neg);
  el_nrm.resize(ne);
  //     set up electron identities
  for (int iat = 0; iat < neg; ++iat)
  {
    el_id[iat] = elns.GroupID[iat];
    el_nrm[el_id[iat]] += 1;
  }
  for (int i = 0; i < ne; ++i)
    el_nrm[i] = 1.0 / el_nrm[i];
  int num_srcs = srcs.groups();
  //use the simulation cell radius if any direction is periodic
  if (elns.Lattice.SuperCellEnum)
    Dmax = elns.Lattice.SimulationCellRadius;
  const SpeciesSet& species(srcs.getSpeciesSet());
  int ng = species.size();
  nag    = srcs.getTotalNum();
  ion_id.resize(nag);
  for (int i = 0; i < ng; ++i)
  {
    for (int j(0); j < num_species; j++)
    {
      std::stringstream nm;
      nm << species.speciesName[i] << "_" << e_species.speciesName[j];
      names.push_back(nm.str());
    }
  }
  for (int iat = 0; iat < nag; ++iat)
    ion_id[iat] = srcs.GroupID[iat];
  lm.resize(num_srcs, num_species);
  lm = 0.0;
}

void LocalMomentEstimator::resetTargetParticleSet(ParticleSet& P) { }

LocalMomentEstimator::Return_t LocalMomentEstimator::evaluate(ParticleSet& P)
{
  const auto& d_table = P.getDistTable(d_table_ID);
  lm = 0;
  for (int iat = 0; iat < nag; ++iat)
  {
    int j(0);
#ifndef ENABLE_SOA
    for (int nn = d_table.M[iat]; nn < d_table.M[iat + 1]; ++nn, j++)
    {
      RealType r = d_table.r(nn);
      if (r >= Dmax)
        continue;
      lm(ion_id[iat], el_id[j]) += el_nrm[el_id[j]];
    }
#endif
  }
  return 0.0;
}

void LocalMomentEstimator::registerCollectables(std::vector<observable_helper*>& h5list, hid_t gid) const {}


void LocalMomentEstimator::addObservables(PropertySetType& plist, BufferType& collectables) { addObservables(plist); }


bool LocalMomentEstimator::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(Dmax, "rcut");
  attrib.put(cur);
  return true;
}

bool LocalMomentEstimator::get(std::ostream& os) const
{
  os << myName << " rcut=" << Dmax << std::endl;
  return true;
}

OperatorBase* LocalMomentEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //default constructor is sufficient
  LocalMomentEstimator* myClone = new LocalMomentEstimator(*this);
  myClone->Dmax                 = Dmax;
  myClone->resetTargetParticleSet(qp);
  return myClone;
}

} // namespace qmcplusplus
