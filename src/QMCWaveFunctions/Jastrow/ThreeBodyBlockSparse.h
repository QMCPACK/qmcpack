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
    
    
#ifndef QMCPLUSPLUS_THREEBODY_BLOCKSPARSE_H
#define QMCPLUSPLUS_THREEBODY_BLOCKSPARSE_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 * @brief ThreeBodyBlockSparse functions
 */
class ThreeBodyBlockSparse: public OrbitalBase
{

public:

  typedef BasisSetBase<RealType> BasisSetType;

  ///constructor
  ThreeBodyBlockSparse(const ParticleSet& ions, ParticleSet& els);

  ~ThreeBodyBlockSparse();

  //implement virtual functions for optimizations
  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
  void resetParameters(const opt_variables_type& active);
  void reportStatus(std::ostream& os);
  //evaluate the distance table with els
  void resetTargetParticleSet(ParticleSet& P);

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L);

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  ValueType ratio(ParticleSet& P, int iat);

  void restore(int iat);

  void acceptMove(ParticleSet& P, int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  void setBasisSet(BasisSetType* abasis)
  {
    GeminalBasis=abasis;
  }

  bool put(xmlNodePtr cur);

  //set blocks
  void setBlocks(const std::vector<int>& blockspergroup);

  ///reference to the center
  const ParticleSet& CenterRef;
  ///assign same blocks for the group
  bool SameBlocksForGroup;
  ///index of the table for source-target
  int myTableIndex;
  ///size of the localized basis set
  int BasisSize;
  ///number of particles
  int NumPtcls;
  ///offset of the index
  int IndexOffSet;
  /** temporary value for update */
  RealType diffVal;
  ///root name for Lambda compoenents
  std::string ID_Lambda;
  /** Y(iat,ibasis) value of the iat-th ortbial, the basis index ibasis
   */
  Matrix<RealType> Y;
  /** dY(iat,ibasis) value of the iat-th ortbial, the basis index ibasis
   */
  Matrix<PosType>  dY;
  /** d2Y(iat,ibasis) value of the iat-th ortbial, the basis index ibasis
   */
  Matrix<RealType> d2Y;
  /** V(i,j) = Lambda(k,kk) U(i,kk)
   */
  Matrix<RealType> V;

  /** Symmetric matrix connecting Geminal Basis functions */
  Matrix<RealType> Lambda;
  /** boolean to enable/disable optmization of Lambda(i,j) component */
  Matrix<int> FreeLambda;
  std::vector<IndexType> BlocksPerGroup;
  std::vector<IndexType> Blocks;
  std::vector<IndexType> BlockOffset;
  std::vector<IndexType> BlockID;
  std::vector<Matrix<RealType>* > LambdaBlocks;

  /** Uk[i] = \sum_j dot(U[i],V[j]) */
  Vector<RealType> Uk;

  /** Gradient for update mode */
  Matrix<PosType> dUk;

  /** Laplacian for update mode */
  Matrix<RealType> d2Uk;

  /** temporary Laplacin for update */
  Vector<RealType> curLap, tLap;
  /** temporary Gradient for update */
  Vector<PosType> curGrad, tGrad;
  /** tempory Lambda*newY for update */
  Vector<RealType> curV;
  /** tempory Lambda*(newY-Y(iat)) for update */
  Vector<RealType> delV;
  /** tempory Lambda*(newY-Y(iat)) for update */
  Vector<RealType> curVal;

  RealType *FirstAddressOfdY;
  RealType *LastAddressOfdY;
  RealType *FirstAddressOfgU;
  RealType *LastAddressOfgU;

  /** Geminal basis function */
  BasisSetType *GeminalBasis;

  /** evaluateLog and store data for particle-by-particle update */
  void evaluateLogAndStore(ParticleSet& P);

  void checkLambda();
};
}
#endif

