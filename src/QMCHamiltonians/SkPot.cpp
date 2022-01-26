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
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("SkPot::evaluate(ParticleSet& P)");
#else
  //sum over species
  copy(P.getSK().rhok[0], P.getSK().rhok[0] + NumK, RhokTot.begin());
  for (int i = 1; i < NumSpecies; ++i)
    accumulate_elements(P.getSK().rhok[i], P.getSK().rhok[i] + NumK, RhokTot.begin());
  Vector<ComplexType>::const_iterator iit(RhokTot.begin()), iit_end(RhokTot.end());
  value_ = 0.0;
  for (int i = 0; iit != iit_end; ++iit, ++i)
    value_ += Fk[i] * ((*iit).real() * (*iit).real() + (*iit).imag() * (*iit).imag());
#endif
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
