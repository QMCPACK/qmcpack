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
    
    
#ifndef QMCPLUSPLUS_TRICUBICSPLINEDBUILDER_H
#define QMCPLUSPLUS_TRICUBICSPLINE3DBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/TriCubicSpline.h"
#include "Numerics/Spline3D/TriCubicSplineSet.h"

namespace qmcplusplus
{

class TriCubicSplineBuilder: public OrbitalBuilderBase
{

  //typedef AnalyticOrbitalSet<TriCubicSpline> TriCubicSplineSet_t;
  typedef SingleParticleOrbitalSet<TriCubicSpline> SPOSet_t;
  //static TriCubicSplineSet orbitals;
  //typedef TriCubicSplineSet SPOSet_t;
  TriCubicSplineSet* m_orbitals;
  Grid3D* grid_ref;

public:

  TriCubicSplineBuilder(TrialWaveFunction& a,
                        Grid3D* agrid): OrbitalBuilderBase(a),
    m_orbitals(NULL),
    grid_ref(agrid)
  { }

  //bool put(xmlNodePtr, Grid3D*);
  bool put(xmlNodePtr);

  //    Grid3D* getFullGrid() { return m_orbitals->getFullGrid();}
};
}
#endif
