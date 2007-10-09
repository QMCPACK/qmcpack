//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SPLINE3DBUILDER_H
#define QMCPLUSPLUS_SPLINE3DBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"
#include "Numerics/Spline3D/Spline3D.h"
#include "Numerics/Spline3D/Spline3DSet.h"

namespace qmcplusplus {

  class Spline3DBuilder: public OrbitalBuilderBase {

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

    Grid3D* getFullGrid() { return d_orbitals->getFullGrid();}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
