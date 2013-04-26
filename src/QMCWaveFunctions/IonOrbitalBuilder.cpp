//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/IonOrbitalBuilder.h"
#include "QMCWaveFunctions/IonOrbital.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
IonOrbitalBuilder::IonOrbitalBuilder
(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets) :
  OrbitalBuilderBase(p, psi), ptclPool(psets)
{
}

bool
IonOrbitalBuilder:: put(xmlNodePtr cur)
{
  ParticleSet &p = targetPtcl;
  OhmmsAttributeSet oAttrib;
  oAttrib.add(sourceOpt, "source");
  oAttrib.add(nameOpt,   "name"  );
  oAttrib.add(widthOpt,  "width" );
  oAttrib.put(cur);
  if(nameOpt == "")
  {
    app_warning() << "  IonOrbitalBuilder::put does not have name "<< endl;
    return false;
  }
  if (!widthOpt.size())
  {
    app_error() << "  You must specify the \"width\" attribute to "
                << " an ionwf.";
    abort();
  }
  map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if(pa_it == ptclPool.end())
  {
    app_error() << "Could not file source ParticleSet "
                << sourceOpt << " for ion wave function.\n";
    return false;
  }
  ParticleSet* sourcePtcl= (*pa_it).second;
  IonOrbital *orb = new IonOrbital (*sourcePtcl, targetPtcl);
  orb->ParticleAlpha.resize(targetPtcl.getTotalNum());
  orb->ParticleCenter.resize(targetPtcl.getTotalNum());
  int num_nonzero = 0;
  for (int iat=0; iat<p.getTotalNum(); iat++)
  {
    RealType w = widthOpt[iat];
    if (w > 0.0)
    {
      orb->ParticleCenter[iat] = num_nonzero++;
      orb->ParticleAlpha[iat]  = 0.5/(w*w);
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
  assert (num_nonzero == sourcePtcl->getTotalNum());
  targetPsi.addOrbital (orb, nameOpt);
  return true;
}
}
