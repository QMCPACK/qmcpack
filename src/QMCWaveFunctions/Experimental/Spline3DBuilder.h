//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SPLINE3DBUILDER_H
#define QMCPLUSPLUS_SPLINE3DBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"
#include "Numerics/Spline3D/Spline3D.h"
#include "Numerics/Spline3D/Spline3DSet.h"

namespace qmcplusplus
{

class Spline3DBuilder: public OrbitalBuilderBase
{

  //typedef AnalyticOrbitalSet<Spline3D> Spline3DSet_t;
  typedef SingleParticleOrbitalSet<Spline3D> SPOSet_t;
  //static Spline3DSet orbitals;

  Spline3DSet *d_orbitals;
  Grid3D* grid_ref;

public:

  Spline3DBuilder(TrialWaveFunction& a): OrbitalBuilderBase(a),
    d_orbitals(NULL),
    grid_ref(NULL)
  { }

  bool put(xmlNodePtr cur);

  Grid3D* getFullGrid()
  {
    return d_orbitals->getFullGrid();
  }
};
}
#endif
