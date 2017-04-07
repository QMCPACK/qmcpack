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
    
    
#ifndef QMCPLUSPLUS_GENERIC_JASTROW_BUILDER_H
#define QMCPLUSPLUS_GENERIC_JASTROW_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

//forward declaration
class ParticleSet;

/**@ingroup WFSBuilder
 *A builder class to add a one- or two-body Jastrow function to a TrialWaveFunction
 *
 *A xml node with OrbtialBuilderBase::jastrow_tag is parsed recursively.
 */
struct JastrowBuilder: public OrbitalBuilderBase
{

  JastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets);

  /**@param cur the current xmlNodePtr to be processed by JastrowBuilder
   *@return true if succesful
   */
  bool put(xmlNodePtr cur);

  template<class JeeType>
  bool createTwoBodySpin(xmlNodePtr cur, JeeType* j2);


  template<class JeeType>
  bool createTwoBodyNoSpin(xmlNodePtr cur, JeeType* j2);


  template<class JneType>
  bool createOneBody(xmlNodePtr cur, JneType* j1);

  ///need ParticleSetPool
  PtclPoolType& ptclPool;

  /// the element name for Correlation, "correlation"
  std::string corr_tag;
};
}
#endif
