//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_QMC_MOLECULARORBITALS2GRID3D_H
#define OHMMS_QMC_MOLECULARORBITALS2GRID3D_H

#include "QMC/QMCApps.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Numerics/TriCubicSplineT.h"

namespace ohmmsqmc {

  /** An application to transform Molecular Orbitals on a regular grid
   */
  class MO2Grid3D: public QMCApps {

  public:

    typedef OrbitalBase::RealType RealType;
    typedef OrbitalBase::ValueType ValueType;
    typedef OrbitalBase::PosType PosType;

    ///constructor
    MO2Grid3D(int argc, char** argv);

    ///destructor
    ~MO2Grid3D();

    ///initialization with a file
    bool init();

    ///generate a set of numerical orbitals
    xmlNodePtr generateNumericalOrbitals(xmlNodePtr cur);

    void copyOrbitalSet(map<string,TriCubicSplineT<ValueType>* >& );
  private:

    bool setParticleSets(xmlNodePtr aroot);
    bool setWavefunctions(xmlNodePtr aroot);

    map<string,TriCubicSplineT<ValueType>* > SPOSet;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
