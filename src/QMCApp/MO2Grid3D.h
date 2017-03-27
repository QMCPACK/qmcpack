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
    
    



#ifndef QMCPLUSPLUS_MOLECULARORBITALS2GRID3D_H
#define QMCPLUSPLUS_MOLECULARORBITALS2GRID3D_H

#include "QMCApp/QMCAppBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Numerics/TriCubicSplineT.h"

namespace qmcplusplus
{

class ParticleSetPool;

/** An application to transform Molecular Orbitals on a regular grid
 */
class MO2Grid3D: public QMCAppBase
{

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
  void copyOrbitalSet(std::map<std::string,TriCubicSplineT<ValueType>* >& other);

private:

  std::string InFileRoot;
  ParticleSetPool* ptclPool;
  MCWalkerConfiguration* Electrons;
  ParticleSet* Ions;
  xmlNodePtr dsetPtr;
  xmlNodePtr normalPtr;
  xmlNodePtr corePtr;
  std::map<std::string,TriCubicSplineT<ValueType>* > SPOSet;

  bool selectCore(xmlNodePtr cur);
  void getEigVectors(xmlNodePtr cur, const Matrix<RealType>& A);
  xmlNodePtr copyDeterminant(xmlNodePtr cur, bool addg);
  xmlNodePtr copyDeterminantSet(xmlNodePtr cur, xmlNodePtr splinePtr);
};
}
#endif
