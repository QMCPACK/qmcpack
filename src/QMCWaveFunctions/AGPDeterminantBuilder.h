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

#include "QMCWaveFunctions/AGPDeterminant.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{

/**@ingroup WFSBuilder
 * @brief An abstract class for wave function builders
 */
class AGPDeterminantBuilder : public WaveFunctionComponentBuilder
{
public:
  AGPDeterminantBuilder(Communicate* comm, ParticleSet& els, const PtclPoolType& pset);

  /// process a xml node at cur
  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

protected:
  ///reference to a PtclPoolType
  const PtclPoolType& ptclPool;
  ///basiset Factory
  std::unique_ptr<SPOSetBuilderFactory> mySPOSetBuilderFactory;
  ///AGPDeterminant
  std::unique_ptr<AGPDeterminant> agpDet;
  std::string funcOpt;
  std::string transformOpt;

  template<typename BasisBuilderT>
  bool createAGP(BasisBuilderT* abuilder, xmlNodePtr cur);
};
} // namespace qmcplusplus
#endif
