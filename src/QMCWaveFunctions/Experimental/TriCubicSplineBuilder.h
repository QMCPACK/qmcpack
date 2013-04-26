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
