///////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Kris Delaney
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
/** @file PWOribitalBuilder.h
 * @brief Declaration of a builder class for PWOrbitalSet
 *
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_V0_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_V0_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
namespace qmcplusplus {

  class PWParameterSet;

  /** OrbitalBuilder for Slater determinants in PW basis
  */
  class PWOrbitalBuilder: public OrbitalBuilderBase {

  private:

    typedef PWOrbitalSet::PWBasisPtr PWBasisPtr;

    ///Read routine for HDF wavefunction file version 0.10
    void ReadHDFWavefunction(hid_t hfile);

    ///hdf5 handler to clean up
    hid_t hfileID;
    ///input twist angle
    PosType TwistAngle;
    ///parameter set
    PWParameterSet* myParam;
    //will do something for twist
    PWBasisPtr myBasisSet;
    //Storage for the orbitals and basis is created in PWOSet.
    std::map<std::string,SPOSetBasePtr> PWOSet;

  public:

    ///constructor
    PWOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);
    ~PWOrbitalBuilder();

    ///implement vritual function
    bool put(xmlNodePtr cur);

  private:
    hid_t getH5(xmlNodePtr cur, const char* aname);
    bool putSlaterDet(xmlNodePtr cur);
    bool createPWBasis(xmlNodePtr cur);
    SPOSetBase* createPW(xmlNodePtr cur, int spinIndex);
    void transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
