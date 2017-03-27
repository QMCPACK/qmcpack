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
    
    
/** @file SplineBasisSet.h
 * @brief Define SplineBasisSet
 */
#ifndef QMCPLUSPLUS_SPLINE3DBASISSET_H
#define QMCPLUSPLUS_SPLINE3DBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/Spline3D/SplineGridHandler.h"

namespace qmcplusplus
{

/** a BasisSetBase to handle real-space wavefunctions
 *
 * template parameters
 * - GT grid type, e.g., TricubicBsplineGrid.
 * - ST storage type, e.g., Array<T,D>
 */
template<typename GT, typename ST>
class SplineBasisSet:
  public BasisSetBase<QMCTraits::ValueType>,
  public SplineGridHandler<QMCTraits::RealType,OHMMS_DIM>
{
public:
  typedef GT   GridType;
  typedef ST   StorageType;

  ///spline grid
  GridType bKnots;
  ///grid hanlder
  //SplineGridHandler<RealType,3> gridHander;
  ///bspline data
  std::vector<const StorageType*> P;

  SplineBasisSet() {}

  //// IMPLEMENT virtual functions of BasisSetBase
  ///do nothing
  void resetParameters(OptimizableSetType& optVariables) {}
  ///resize the basis set
  void setBasisSetSize(int nbs)
  {
    BasisSetSize=nbs;
    nbs-=P.size();
    while(nbs>0)
    {
      P.push_back(0);
      nbs--;
    }
  }
  ///reset the target particle set
  void resetTargetParticleSet(ParticleSet& e) {}

  void evaluateForPtclMove(const ParticleSet& e, int iat)
  {
    if(Orthorhombic)
    {
      bKnots.Find(e.R[iat][0],e.R[iat][1],e.R[iat][2]);
    }
    else
    {
      PosType ru(Lattice.toUnit(e.R[iat]));
      bKnots.Find(ru[0],ru[1],ru[2]);
    }
    for(int j=0; j<BasisSetSize; j++)
      Phi[j]=bKnots.evaluate(*P[j]);
  }

  void evaluateAllForPtclMove(const ParticleSet& e, int iat)
  {
    if(Orthorhombic)
    {
      bKnots.FindAll(e.R[iat][0],e.R[iat][1],e.R[iat][2]);
      for(int j=0; j<BasisSetSize; j++)
        Phi[j]=bKnots.evaluate(*P[j],dPhi[j],d2Phi[j]);
    }
    else
    {
      PosType ru(Lattice.toUnit(e.R[iat]));
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      for(int j=0; j<BasisSetSize; j++)
      {
        Phi[j]=bKnots.evaluate(*P[j],gu,hess);
        dPhi[j]=dot(Lattice.G,gu);
      }
    }
  }

  void evaluateForWalkerMove(const ParticleSet& e)
  {
    std::cout << "WHO IS CALLING SplineBasisSet:: evaluateForWalkerMove(P)  " << std::endl;
  }

  void evaluateForWalkerMove(const ParticleSet& e, int iat)
  {
    evaluateAllForPtclMove(e,iat);
  }

  /** add basis orbital
   * @param i basisset orbital index
   * @param data raw data on a grid
   * @param curP splined data
   */
  void add(int i, const StorageType& data, StorageType* curP)
  {
    if(i<BasisSetSize)
    {
      bKnots.Init(data,*curP);
      P[i]=curP;
    }
  }

  /** add basis orbital
   * @param i basisset orbital index
   * @param curP splined data
   */
  void add(int i, StorageType* curP)
  {
    if(i<BasisSetSize)
      P[i]=curP;
  }

};

}
#endif
