//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file WaveFunctionFactory.h
 *@brief Declaration of a WaveFunctionFactory
 */
#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_FACTORY_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_FACTORY_H

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetFactory.h"
#include "Message/MPIObjectBase.h"
namespace qmcplusplus
{

/** Factory class to build a many-body wavefunction
 */
struct WaveFunctionFactory: public MPIObjectBase
{

  typedef std::map<std::string,ParticleSet*> PtclPoolType;
  ///target ParticleSet
  ParticleSet* targetPtcl;
  ///many-body wavefunction object
  TrialWaveFunction* targetPsi;
  ///reference to the PtclPoolType
  PtclPoolType&  ptclPool;
  ///input node for a many-body wavefunction
  xmlNodePtr myNode;
  ///builder tree
  std::vector<OrbitalBuilderBase*> psiBuilder;

  /** constructor
   * @param qp quantum particleset
   * @param pset pool of particlesets
   * @param c  communicator
   */
  WaveFunctionFactory(ParticleSet* qp, PtclPoolType& pset, Communicate* c);

  ~WaveFunctionFactory();

  void setPsi(TrialWaveFunction* psi);

  ///read from xmlNode
  bool put(xmlNodePtr cur);

  ///reset member data
  void reset();

  /** process xmlNode to populate targetPsi
   */
  bool build(xmlNodePtr cur, bool buildtree=true);

  /** add Fermion wavefunction term */
  bool addFermionTerm(xmlNodePtr cur);

  /** add finite-difference linear response wavefunction term */
  bool addFDLRTerm(xmlNodePtr cur);

  /** add an OrbitalBuilder and the matching xml node
   * @param b OrbitalBuilderBase*
   * @oaram cur xmlNode for b
   * @return true if successful
   */
  bool addNode(OrbitalBuilderBase* b, xmlNodePtr cur);

  void setCloneSize(int np);

  WaveFunctionFactory* clone(ParticleSet* qp, int ip, const std::string& aname);

  std::vector<WaveFunctionFactory*> myClones;
};

}
#endif
