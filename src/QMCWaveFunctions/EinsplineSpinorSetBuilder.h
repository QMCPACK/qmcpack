//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file EinsplineSetBuilder.h
 *
 * Builder class for einspline-based SPOSet objects.
 */
#ifndef QMCPLUSPLUS_EINSPLINE_SPINORSET_BUILDER_H
#define QMCPLUSPLUS_EINSPLINE_SPINORSET_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"

class Communicate;

namespace qmcplusplus
{

/** EinsplineSpinorSet builder
 */
class EinsplineSpinorSetBuilder : public SPOSetBuilder
{
  typedef std::map<std::string, ParticleSet*> PtclPoolType;
public:
  ///constructor
  EinsplineSpinorSetBuilder(ParticleSet& p, PtclPoolType& psets, Communicate* comm, xmlNodePtr cur):SPOSetBuilder(comm){};

  ///destructor
  ~EinsplineSpinorSetBuilder(){};

  /** initialize the Antisymmetric wave function for electrons
   * @param cur the current xml node
   */
  SPOSet* createSPOSetFromXML(xmlNodePtr cur){return nullptr;};

};

} // namespace qmcplusplus


#endif
