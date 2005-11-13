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
#ifndef QMCPLUSPLUS_MOLECULARORBITALS2GRID3D_H
#define QMCPLUSPLUS_MOLECULARORBITALS2GRID3D_H

#include "QMCApp/QMCAppBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Numerics/TriCubicSplineT.h"

namespace qmcplusplus {

  class ParticleSetPool;

  /** An application to transform Molecular Orbitals on a regular grid
   */
  class MO2Grid3D: public QMCAppBase {

  public:

    typedef OrbitalBase::RealType RealType;
    typedef OrbitalBase::ValueType ValueType;
    typedef OrbitalBase::PosType PosType;

    ///constructor
    MO2Grid3D(int argc, char** argv);

    ///destructor
    ~MO2Grid3D();

    ///validate an input file
    bool validateXML();

    ///do something
    bool execute();

    ///generate a set of numerical orbitals
    xmlNodePtr generateNumericalOrbitals(xmlNodePtr cur);

    /** transfer the wavefunction to other
     * @param other container to which  SPOSet will be copied.
     *
     * This function is introduced for an OrbitalBuilder object
     * which takes a normal XML file for Molecular Orbitals represented
     * by radial orbitals and spherical harmonics and requests
     * mapping to numerical 3D data for QMC. 
     * See QMCWaveFunctions/NumericalOrbitalBuilder.cpp
     */
    void copyOrbitalSet(map<string,TriCubicSplineT<ValueType>* >& other);

  private:

    string InFileRoot;
    ParticleSetPool* ptclPool;
    MCWalkerConfiguration* Electrons;
    ParticleSet* Ions;
    xmlNodePtr dsetPtr;
    xmlNodePtr normalPtr;
    xmlNodePtr corePtr;
    map<string,TriCubicSplineT<ValueType>* > SPOSet;

    bool selectCore(xmlNodePtr cur);
    void getEigVectors(xmlNodePtr cur, const Matrix<RealType>& A);
    xmlNodePtr copyDeterminant(xmlNodePtr cur, bool addg);
    xmlNodePtr copyDeterminantSet(xmlNodePtr cur, xmlNodePtr splinePtr);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
