//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SkPot.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
SkPot::SkPot(ParticleSet& source)
{
  sourcePtcl = &source;
  NumSpecies = source.getSpeciesSet().getTotalNum();
  NumK       = source.getSimulationCell().getKLists().numk;
  OneOverN   = 1.0 / static_cast<RealType>(source.getTotalNum());
  Kshell     = source.getSimulationCell().getKLists().kshell;
  MaxKshell  = Kshell.size() - 1;
  RhokTot.resize(NumK);
  Fk.resize(NumK);
  Kmag.resize(MaxKshell);
  OneOverDnk.resize(MaxKshell);
  for (int ks = 0; ks < MaxKshell; ks++)
  {
    Kmag[ks]       = std::sqrt(source.getSimulationCell().getKLists().ksq[Kshell[ks]]);
    OneOverDnk[ks] = 1.0 / static_cast<RealType>(Kshell[ks + 1] - Kshell[ks]);
  }
}

void SkPot::resetTargetParticleSet(ParticleSet& P) { sourcePtcl = &P; }

SkPot::Return_t SkPot::evaluate(ParticleSet& P)
{
  throw std::runtime_error("SkPot::evaluate not implemented. There was an implementation with"
                           " complex-valued storage that may be resurrected using real-valued storage.");
  return value_;
}

bool SkPot::put(xmlNodePtr cur)
{
  OhmmsAttributeSet Tattrib;
  Tattrib.add(K_0, "k0");
  Tattrib.add(V_0, "v0");
  Tattrib.put(cur);
  app_log() << "KSpacePot parameters" << std::endl;
  app_log() << "  k0: " << K_0 << std::endl;
  app_log() << "  v0: " << V_0 << std::endl;
  FillFk();
  return true;
}


bool SkPot::get(std::ostream& os) const { return true; }

std::unique_ptr<OperatorBase> SkPot::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<SkPot> myclone = std::make_unique<SkPot>(*this);
  myclone->FillFk();
  return myclone;
}
} // namespace qmcplusplus
