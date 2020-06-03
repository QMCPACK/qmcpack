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


/** @file
 * @brief Declaration of a builder class for PWOrbitalSet
 *
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_V0_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_V0_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
#else
#include "QMCWaveFunctions/PlaneWave/PWRealOrbitalSet.h"
#endif
namespace qmcplusplus
{
class PWParameterSet;
class SlaterDet;

/** OrbitalBuilder for Slater determinants in PW basis
*/
class PWOrbitalBuilder : public WaveFunctionComponentBuilder
{
private:
#if defined(QMC_COMPLEX)
  typedef PWOrbitalSet SPOSetType;
  typedef PWOrbitalSet::PWBasisPtr PWBasisPtr;
#else
  typedef PWRealOrbitalSet SPOSetType;
  typedef PWRealOrbitalSet::PWBasisPtr PWBasisPtr;
#endif

  std::map<std::string, SPOSetPtr> spomap;
  PtclPoolType& ptclPool;

  ///Read routine for HDF wavefunction file version 0.10
  void ReadHDFWavefunction(hid_t hfile);

  ///hdf5 handler to clean up
  hid_t hfileID;
  ///xml node for determinantset
  xmlNodePtr rootNode;
  ///input twist angle
  PosType TwistAngle;
  ///parameter set
  PWParameterSet* myParam;
  //will do something for twist
  PWBasisPtr myBasisSet;
  ////Storage for the orbitals and basis is created in PWOSet.
  //std::map<std::string,SPOSetPtr> PWOSet;
public:
  ///constructor
  PWOrbitalBuilder(Communicate* comm, ParticleSet& els, PtclPoolType& psets);
  ~PWOrbitalBuilder();

  ///implement vritual function
  WaveFunctionComponent* buildComponent(xmlNodePtr cur) override;

private:
  hid_t getH5(xmlNodePtr cur, const char* aname);
  WaveFunctionComponent* putSlaterDet(xmlNodePtr cur);
  bool createPWBasis(xmlNodePtr cur);
  SPOSet* createPW(xmlNodePtr cur, int spinIndex);
#if defined(QMC_COMPLEX)
  void transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc);
#endif
};
} // namespace qmcplusplus
#endif
