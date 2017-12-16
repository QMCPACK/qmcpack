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


#ifndef QMCPLUSPLUS_GTO_BUILDER_H
#define QMCPLUSPLUS_GTO_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"
#include "Numerics/GaussianBasisSet.h"
#include "io/hdf_archive.h"
namespace qmcplusplus
{

class GTOBuilder: public QMCTraits
{

public:

  typedef GaussianCombo<RealType>                    RadialOrbitalType;
  typedef SphericalBasisSet<RadialOrbitalType>     CenteredOrbitalType;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;
  ///the species
  std::string m_species;
  ///constructor
  GTOBuilder(xmlNodePtr cur=NULL);

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter)
  {
    m_orbitals = oset;
    m_species = acenter;
  }

  bool addGrid(xmlNodePtr cur)
  {
    return true;
  }

  bool addGridH5(hdf_archive &hin)
  {
    return true;
  }

  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);
  bool addRadialOrbitalH5(hdf_archive &hin, const QuantumNumberType& nlms)
  {
    return true;
  }

  bool putCommon(xmlNodePtr cur);

};
}
#endif
