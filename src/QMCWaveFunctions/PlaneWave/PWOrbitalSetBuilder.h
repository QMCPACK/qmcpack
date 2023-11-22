//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 * @brief Declaration of a builder class for PWOrbitalSet
 *
 */
#ifndef QMCPLUSPLUS_PWORBITAL_BUILDER_H
#define QMCPLUSPLUS_PWORBITAL_BUILDER_H

#include "SPOSetBuilder.h"
#include "hdf/hdf_archive.h"
#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
#else
#include "QMCWaveFunctions/PlaneWave/PWRealOrbitalSet.h"
#endif

namespace qmcplusplus
{
struct PWParameterSet;
class SlaterDet;

/** OrbitalBuilder for Slater determinants in PW basis
*/
class PWOrbitalSetBuilder : public SPOSetBuilder
{
private:
#if defined(QMC_COMPLEX)
  using SPOSetType = PWOrbitalSet;
#else
  using SPOSetType = PWRealOrbitalSet;
#endif

  /// target particle set
  const ParticleSet& targetPtcl;
  ///xml node for determinantset
  xmlNodePtr rootNode{nullptr};
  ///input twist angle
  PosType TwistAngle;
  ///parameter set
  std::unique_ptr<PWParameterSet> myParam;
  //will do something for twist
  std::unique_ptr<PWBasis> myBasisSet;
  ///hdf5 handler to clean up
  hdf_archive hfile;

public:
  ///constructor
  PWOrbitalSetBuilder(const ParticleSet& p, Communicate* comm, xmlNodePtr cur);
  ~PWOrbitalSetBuilder() override;

  /// create an sposet from xml and save the resulting SPOSet
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

private:
  bool getH5(xmlNodePtr cur, const char* aname);
  bool createPWBasis();
  std::unique_ptr<SPOSet> createPW(xmlNodePtr cur, const std::string& objname, int spinIndex);
#if defined(QMC_COMPLEX)
  void transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc);
#endif
};
} // namespace qmcplusplus
#endif
