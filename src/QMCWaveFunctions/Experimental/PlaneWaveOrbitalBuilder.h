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
    
    
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_H

#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/PlaneWaveBasis.h"
#include "QMCWaveFunctions/PlaneWaveOrbitalSet.h"
namespace qmcplusplus
{

/** OrbitalBuilder for Slater determinants in PW basis
*/
class PlaneWaveOrbitalBuilder: public OrbitalBuilderBase
{

private:
  ///Number of up and down particles from ParticleSet
  int nup, ndown, upindx;
  ///Index of spin data from HDF5 file to use for updet and downdet
  int updetspinindex, downdetspinindex;

  ///Read routine for HDF wavefunction file version 0.10
  void ReadHDFWavefunction010(hid_t hfile,double& ecut);

  //Storage for the orbitals and basis is created in PWOSet.
  PlaneWaveOrbitalSet *PWOSet;

public:
  ///constructor
  PlaneWaveOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);
  ~PlaneWaveOrbitalBuilder();

  ///implement vritual function
  bool put(xmlNodePtr cur);
};
}
#endif
