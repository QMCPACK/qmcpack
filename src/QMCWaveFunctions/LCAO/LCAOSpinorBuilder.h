//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOA_LCAO_SPINOR_BUILDER_H
#define QMCPLUSPLUS_SOA_LCAO_SPINOR_BUILDER_H

#include "QMCWaveFunctions/LCAO/LCAOrbitalBuilder.h"

namespace qmcplusplus
{
/** @file LCAOSpinorBuidler.h
   *
   * Derives from LCAOrbitalBuilder.h. Overrides createSPOSetFromXML method to read up and
   * down channel from HDF5 and construct SpinorSet
   * 
   */
class LCAOSpinorBuilder : public LCAOrbitalBuilder
{
public:
  /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
  LCAOSpinorBuilder(ParticleSet& els, ParticleSet& ions, Communicate* comm, xmlNodePtr cur);
  SPOSet* createSPOSetFromXML(xmlNodePtr cur);

private:
  bool loadMO(LCAOrbitalSet& up, LCAOrbitalSet& dn, xmlNodePtr cur);
  bool putFromH5(LCAOrbitalSet& up, LCAOrbitalSet& dn);
};
} // namespace qmcplusplus
#endif
