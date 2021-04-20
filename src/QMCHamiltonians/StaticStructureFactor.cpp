//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "StaticStructureFactor.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/KContainer.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus
{
StaticStructureFactor::StaticStructureFactor(ParticleSet& P) : Pinit(P)
{
#ifndef USE_REAL_STRUCT_FACTOR
  APP_ABORT("StaticStructureFactor: please recompile with USE_REAL_STRUCT_FACTOR=1");
#endif
  if (P.Lattice.SuperCellEnum == SUPERCELL_OPEN)
    APP_ABORT("StaticStructureFactor is incompatible with open boundary conditions");

  // get particle information
  SpeciesSet& species = P.getSpeciesSet();
  nspecies            = species.size();
  for (int s = 0; s < nspecies; ++s)
    species_name.push_back(species.speciesName[s]);
  reset();
}


void StaticStructureFactor::reset()
{
  myName = "StaticStructureFactor";
  UpdateMode.set(COLLECTABLE, 1);
  ecut     = -1.0;
  nkpoints = -1;
}


std::shared_ptr<OperatorBase> StaticStructureFactor::makeClone(ParticleSet& P, TrialWaveFunction& Psi)
{
  return std::make_shared<StaticStructureFactor>(*this);
}


bool StaticStructureFactor::put(xmlNodePtr cur)
{
  using std::sqrt;
  reset();
  const k2_t& k2_init = Pinit.SK->KLists.ksq;

  std::string write_report = "no";
  OhmmsAttributeSet attrib;
  attrib.add(myName, "name");
  //attrib.add(ecut,"ecut");
  attrib.add(write_report, "report");
  attrib.put(cur);

  if (ecut < 0.0)
    nkpoints = k2_init.size();
  else
  {
    RealType k2cut = 2.0 * ecut;
    nkpoints       = 0;
    for (int i = 0; i < k2_init.size(); ++i)
      if (k2_init[i] < k2cut)
        nkpoints++;
  }

  if (nkpoints == 0)
    APP_ABORT("StaticStructureFactor::put  could not find any kpoints");

  ecut = .5 * k2_init[nkpoints - 1];

  if (write_report == "yes")
    report("  ");

  return true;
}


void StaticStructureFactor::report(const std::string& pad)
{
  app_log() << pad << "StaticStructureFactor report" << std::endl;
  app_log() << pad << "  name     = " << myName << std::endl;
  app_log() << pad << "  ecut     = " << ecut << std::endl;
  app_log() << pad << "  nkpoints = " << nkpoints << std::endl;
  app_log() << pad << "  nspecies = " << nspecies << std::endl;
  for (int s = 0; s < nspecies; ++s)
    app_log() << pad << "    species[" << s << "] = " << species_name[s] << std::endl;
  app_log() << pad << "end StaticStructureFactor report" << std::endl;
}


void StaticStructureFactor::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex = collectables.current();
  std::vector<RealType> tmp(nspecies * 2 * nkpoints); // real & imag parts
  collectables.add(tmp.begin(), tmp.end());
}


void StaticStructureFactor::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const
{
  hid_t sgid            = H5Gcreate(gid, myName.c_str(), 0);
  observable_helper* oh = new observable_helper("kpoints");
  oh->open(sgid); // add to SkAll hdf group
  oh->addProperty(const_cast<std::vector<PosType>&>(Pinit.SK->KLists.kpts_cart), "value");
  h5desc.push_back(oh);
  std::vector<int> ng(2);
  ng[0] = 2;
  ng[1] = nkpoints;
  for (int s = 0; s < nspecies; ++s)
  {
    observable_helper* oh = new observable_helper(species_name[s]);
    oh->set_dimensions(ng, myIndex + s * 2 * nkpoints);
    oh->open(sgid);
    h5desc.push_back(oh);
  }
}


StaticStructureFactor::Return_t StaticStructureFactor::evaluate(ParticleSet& P)
{
  RealType w                     = tWalker->Weight;
  const Matrix<RealType>& rhok_r = P.SK->rhok_r;
  const Matrix<RealType>& rhok_i = P.SK->rhok_i;
  int nkptot                     = rhok_r.cols();
  for (int s = 0; s < nspecies; ++s)
  {
    int kc = myIndex + s * 2 * nkpoints;
    //int kstart  = s*nkptot;
    //for(int k=kstart;k<kstart+nkpoints;++k,++kc)
    //  P.Collectables[kc] += w*rhok_r(k);
    //for(int k=kstart;k<kstart+nkpoints;++k,++kc)
    //  P.Collectables[kc] += w*rhok_i(k);
    for (int k = 0; k < nkpoints; ++k, ++kc)
      P.Collectables[kc] += w * rhok_r(s, k);
    for (int k = 0; k < nkpoints; ++k, ++kc)
      P.Collectables[kc] += w * rhok_i(s, k);
  }
  return 0.0;
}

} // namespace qmcplusplus
