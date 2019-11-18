//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPINORBASISSETFACTORY_H
#define QMCPLUSPLUS_SPINORBASISSETFACTORY_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

/** @file SpinorSetBuilderFactory.h
 *  This class derives from SPOSetBuilderFactory.  It's purpose is to override createSPOSetBuilder so that 
 *   spinor appropriate factories can be parsed and called.  Example, <spinorset_builder> knows to route einspline
 *   construction through EinsplineSpinorSetBuilder as oppposed to EinsplineSetBuilder. 
 */

namespace qmcplusplus
{

class SpinorSetBuilderFactory  : public SPOSetBuilderFactory
{
public:
  typedef std::map<std::string, ParticleSet*> PtclPoolType;


  SpinorSetBuilderFactory(Communicate* comm, ParticleSet& els, PtclPoolType& psets);

  ~SpinorSetBuilderFactory(){};
  SPOSetBuilder* createSPOSetBuilder(xmlNodePtr rootNode);

};
} // namespace qmcplusplus
#endif
