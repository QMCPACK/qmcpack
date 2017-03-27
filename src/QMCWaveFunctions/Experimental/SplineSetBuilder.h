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
    
    
#ifndef QMCPLUSPLUS_SPLINESETBUILDER_H
#define QMCPLUSPLUS_SPLINESETBUILDER_H

#include "Numerics/TriCubicSplineT.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"

namespace qmcplusplus
{

/**@ingroup WFSBuilder
 * A builder class for a set of Spline functions
 */
class SplineSetBuilder: public BasisSetBuilder
{

public:

  typedef TriCubicSplineT<ValueType,RealType>           SPOType;
  typedef TriCubicSplineT<ValueType,RealType>::GridType GridType;
  typedef SingleParticleOrbitalSet<SPOType>             SPOSetType;
  typedef std::map<std::string,ParticleSet*> PtclPoolType;

  /** constructor
   * @param p target ParticleSet
   * @param psets a set of ParticleSet objects
   */
  SplineSetBuilder(ParticleSet& p, PtclPoolType& psets);

  bool put(xmlNodePtr cur);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   */
  SPOSetBase* createSPOSet(xmlNodePtr cur);

private:
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///reference to a ParticleSetPool
  PtclPoolType& ptclPool;
  ///global GridType*, should be generalized for any number of grids
  GridType* GridXYZ;
  ///set of SPOType*
  std::map<std::string,SPOType* > NumericalOrbitals;
  ///set of SPOSetType*
  std::map<std::string,SPOSetType*> SPOSet;
};
}
#endif
