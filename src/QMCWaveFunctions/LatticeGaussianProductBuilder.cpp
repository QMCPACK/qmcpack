//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "LatticeGaussianProductBuilder.h"
#include "QMCWaveFunctions/LatticeGaussianProduct.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
LatticeGaussianProductBuilder::LatticeGaussianProductBuilder(Communicate* comm, ParticleSet& p, const PtclPoolType& psets)
    : WaveFunctionComponentBuilder(comm, p), ptclPool(psets)
{}

std::unique_ptr<WaveFunctionComponent> LatticeGaussianProductBuilder::buildComponent(xmlNodePtr cur)
{
  ParticleSet& p = targetPtcl;
  // initialize widths to zero; if no user input, then abort
  widthOpt.resize(targetPtcl.getTotalNum(), 0);
  OhmmsAttributeSet oAttrib;
  oAttrib.add(sourceOpt, "source");
  oAttrib.add(nameOpt, "name");
  oAttrib.add(widthOpt, "width");
  oAttrib.put(cur);
  if (nameOpt == "")
  {
    app_warning() << "  LatticeGaussianProductBuilder::put does not have name " << std::endl;
  }

  auto pa_it(ptclPool.find(sourceOpt));
  if (pa_it == ptclPool.end())
  {
    app_error() << "Could not file source ParticleSet " << sourceOpt << " for ion wave function.\n";
  }
  ParticleSet* sourcePtcl = (*pa_it).second;
  auto orb                = std::make_unique<LatticeGaussianProduct>(*sourcePtcl, targetPtcl);
  orb->ParticleAlpha.resize(targetPtcl.getTotalNum());
  orb->ParticleCenter.resize(targetPtcl.getTotalNum());
  int num_nonzero = 0;
  for (int iat = 0; iat < p.getTotalNum(); iat++)
  {
    RealType w = widthOpt[iat];
    if (w > 0.0)
    {
      orb->ParticleCenter[iat] = num_nonzero++;
      orb->ParticleAlpha[iat]  = 0.5 / (w * w);
    }
    else
    {
      orb->ParticleCenter[iat] = -1;
      orb->ParticleAlpha[iat]  = 0.0;
    }
  }
  if (num_nonzero != sourcePtcl->getTotalNum())
  {
    app_error() << "  The number of nonzero widths should be the same as the number of\n"
                << "  centers for the ionwf.\n";
    abort();
  }
  assert(num_nonzero == sourcePtcl->getTotalNum());
  return orb;
}

} // namespace qmcplusplus
