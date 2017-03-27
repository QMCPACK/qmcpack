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
    
    
#ifndef QMCPLUSPLUS_GAUSSIANORBITAL_MOLECULARORBITALS_H
#define QMCPLUSPLUS_GAUSSIANORBITAL_MOLECULARORBITALS_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBasis.h"
#include "Numerics/GaussianBasisSet.h"

namespace qmcplusplus
{

/**Class to add a set of Gaussian Type atomic orbital basis functions
 *to the collection of basis functions.
 *
 @brief Example of a ROT (Radial Orbital Type)
*/
class GTOMolecularOrbitals: public OrbitalBuilderBase
{
public:

  typedef GaussianCombo<RealType>                    RadialOrbitalType;
  typedef SphericalOrbitalSet<RadialOrbitalType>     CenteredOrbitalType;
  typedef MolecularOrbitalBasis<CenteredOrbitalType> BasisSetType;

  ///constructor
  GTOMolecularOrbitals(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions);

  ///implement vritual function
  bool put(xmlNodePtr cur);

  ///returns a BasisSet
  BasisSetType* addBasisSet(xmlNodePtr cur);

private:

  enum {DONOT_EXPAND=0, GAUSSIAN_EXPAND=1, NATURAL_EXPAND};

  bool Normalized;
  ParticleSet& IonSys;
  BasisSetType*      BasisSet;
  DistanceTableData* d_table;
  std::map<std::string,int>    RnlID;
  std::map<std::string,int>    CenterID;
  ///map for (n,l,m,s) to its quantum number index
  std::map<std::string,int> nlms_id;

  int expandYlm(const std::string& rnl, const QuantumNumberType& nlms,
                int num, CenteredOrbitalType* aos, xmlNodePtr cur1,
                int expandlm);

};
}
#endif
