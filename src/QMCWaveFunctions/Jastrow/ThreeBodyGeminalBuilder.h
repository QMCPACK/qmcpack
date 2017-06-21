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
    
    
/**@file ThreeBodyGeminalBuilder.h
 *@brief declaration of three-body jastrow of Geminal functions
 */
#ifndef QMCPLUSPLUS_THREEBODY_GEMINALBUILDER_H
#define QMCPLUSPLUS_THREEBODY_GEMINALBUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
namespace qmcplusplus
{

class ThreeBodyGeminal;
class BasisSetBuilder;

/**@ingroup WFSBuilder
 * @brief An abstract class for wave function builders
 */
class ThreeBodyGeminalBuilder: public OrbitalBuilderBase
{

public:

  ThreeBodyGeminalBuilder(ParticleSet& els, TrialWaveFunction& wfs,
                          ParticleSet& ions);

  /// process a xml node at cur
  bool put(xmlNodePtr cur);

protected:

  ParticleSet& sourcePtcl;
  BasisSetBuilder* basisBuilder;
};
}
#endif
