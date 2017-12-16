//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GENERIC_MOLECULARORBITALBUILDER_H
#define QMCPLUSPLUS_GENERIC_MOLECULARORBITALBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

class ParticleSet;

/** Builder class to create Fermi wavefunctions using MO Basis.
 *
 * Any molecular orbital has to use the distance table between the ions and electrons.
 * The traits of MO Basis is that each basis function is
 * \f$ \phi_{nlm} =  \frac{R_{nl}}{r^l} r^{l}\Re (Y_{lm}) \f$
 * where \f$ \frac{R_{nl}}{r^l} \f$ is any radial grid function and
 * \f$ \Re (Y_{lm}) \f$ is a real spherical harmonic and \f$ r^l \Re (Y_{lm})\f$ handled by
 * SphericalTensor object.
 * Depending upon the input, one can convert the functions on the radial
 * grids or can carry on the calculations using the input functions.
 */
struct MolecularOrbitalBuilder: public OrbitalBuilderBase
{

  typedef std::map<std::string,ParticleSet*> PtclPoolType;

  MolecularOrbitalBuilder(ParticleSet& p, TrialWaveFunction& psi,
                          PtclPoolType& psets):OrbitalBuilderBase(p,psi),ptclPool(psets) { }

  bool put(xmlNodePtr cur);
  bool putSpecial(xmlNodePtr cur);
  bool putOpen(const std::string& fname_in);

  ///need ParticleSetPool
  PtclPoolType& ptclPool;
};


}
#endif
