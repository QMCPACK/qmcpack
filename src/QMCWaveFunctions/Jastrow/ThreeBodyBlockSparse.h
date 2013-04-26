//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  void reportStatus(ostream& os);
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

  /** later merge the loop */
  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL);

  /** later merge the loop */
  ValueType logRatio(ParticleSet& P, int iat,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL);

  void restore(int iat);

  void acceptMove(ParticleSet& P, int iat);

  inline void update(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL,
                     int iat);

  RealType registerData(ParticleSet& P, PooledData<RealType>& buf);

  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

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
  string ID_Lambda;
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
  vector<IndexType> BlocksPerGroup;
  vector<IndexType> Blocks;
  vector<IndexType> BlockOffset;
  vector<IndexType> BlockID;
  vector<Matrix<RealType>* > LambdaBlocks;

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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1796 $   $Date: 2007-02-22 11:40:21 -0600 (Thu, 22 Feb 2007) $
 * $Id: ThreeBodyBlockSparse.h 1796 2007-02-22 17:40:21Z jnkim $
 ***************************************************************************/

