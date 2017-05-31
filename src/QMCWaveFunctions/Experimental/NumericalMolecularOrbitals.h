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
    
    
#ifndef QMCPLUSPLUS_NUMERICAL_MOLECULARORBITALS_H
#define QMCPLUSPLUS_NUMERICAL_MOLECULARORBITALS_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

class GridMolecularOrbitals;

class NumericalMolecularOrbitals: public OrbitalBuilderBase
{

  GridMolecularOrbitals *Original;

public:

  /** constructor
   * \param wfs reference to the wavefunction
   * \param ions reference to the ions
   * \param els reference to the electrons
   */
  NumericalMolecularOrbitals(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   *
   */
  bool put(xmlNodePtr cur);

private:

};
}
#endif
