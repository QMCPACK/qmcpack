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
    
    
#ifndef QMCPLUSPLUS_NUMERICAL_JASTROW_AA_BUILDER_H
#define QMCPLUSPLUS_NUMERICAL_JASTROW_AA_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/NumericalJastrowFunctor.h"

namespace qmcplusplus
{

//forward declaration
class ParticleSet;

/**@ingroup WFSBuilder
 *A builder class to add a numerical two-body Jastrow function to a TrialWaveFunction
 *
 *A xml node with OrbtialBuilderBase::jastrow_tag is parsed recursively.
 */
struct NJAABuilder: public OrbitalBuilderBase
{

  typedef NumericalJastrow<RealType> FuncType;
  typedef FuncType::FNIN              InFuncType;
  typedef FuncType::FNOUT             OutFuncType;

  ///true for spin-indepedent
  bool IgnoreSpin;
  ///pointer to <grid/>
  xmlNodePtr gridPtr;
  ///unique analytic functions
  std::vector<InFuncType*> InFunc;

  NJAABuilder(ParticleSet& p, TrialWaveFunction& psi);

  /**@param cur the current xmlNodePtr to be processed by NumericalJastrowBuilder
   *@return true if succesful
   */
  bool put(xmlNodePtr cur);

  /** build UniqueIn and InFunc
   *
   * Initialize the input analytic functions
   */
  bool putInFunc(xmlNodePtr cur);

  /** create an analytic function for speciesA-speciesB pair
   * @param jastfunction name of JastrowFunctor
   * @param speciesA species index
   * @param speciesB species index
   *
   * Check IgnoreSpin flag and speciesA-speciesB to instantiate a proper
   * function, e.g., PadeJastrow
   */
  InFuncType* createInFunc(const std::string& jastfunction, int speciesA, int speciesB);

};
}
#endif
