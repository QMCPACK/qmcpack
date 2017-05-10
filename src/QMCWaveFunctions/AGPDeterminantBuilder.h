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
    
    
/**@file AGPDeterminantBuilder.h
 *@brief declaration of a builder class for AGPDeterminant
 */
#ifndef QMCPLUSPLUS_AGPDETERMINANT_GEMINALBUILDER_H
#define QMCPLUSPLUS_AGPDETERMINANT_GEMINALBUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetFactory.h"
namespace qmcplusplus
{

class AGPDeterminant;

/**@ingroup WFSBuilder
 * @brief An abstract class for wave function builders
 */
class AGPDeterminantBuilder: public OrbitalBuilderBase
{

public:

  AGPDeterminantBuilder(ParticleSet& els, TrialWaveFunction& wfs, PtclPoolType& pset);

  /// process a xml node at cur
  bool put(xmlNodePtr cur);

protected:

  ///reference to a PtclPoolType
  PtclPoolType& ptclPool;
  ///basiset Factory
  BasisSetFactory* myBasisSetFactory;
  ///AGPDeterminant
  AGPDeterminant* agpDet;
  std::string funcOpt;
  std::string transformOpt;

  template <typename BasisBuilderT>
  bool createAGP(BasisBuilderT* abuilder, xmlNodePtr cur);

};
}
#endif
