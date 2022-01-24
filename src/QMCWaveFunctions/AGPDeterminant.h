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


/** @file AGPDeterminant.h
 * @brief Declaration of AGPDeterminant for pairing orbitals.
 */
#ifndef QMCPLUSPLUS_AGP_DIRACDETERMINANT_H
#define QMCPLUSPLUS_AGP_DIRACDETERMINANT_H
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus
{
class AGPDeterminant : public WaveFunctionComponent
{
public:
  ///define BasisSetType with RealType
  using BasisSetType = BasisSetBase<RealType>;
  using IndexVector  = BasisSetType::IndexVector;
  using ValueVector  = BasisSetType::ValueVector;
  using ValueMatrix  = BasisSetType::ValueMatrix;
  using GradVector   = BasisSetType::GradVector;
  using GradMatrix   = BasisSetType::GradMatrix;

  BasisSetType* GeminalBasis;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  AGPDeterminant(BasisSetType* bs = nullptr);

  ///default destructor
  ~AGPDeterminant() override;

  void checkInVariables(opt_variables_type& active) override;
  void checkOutVariables(const opt_variables_type& active) override;
  void resetParameters(const opt_variables_type& active) override;
  void reportStatus(std::ostream& os) override;

  ///reset the size: with the number of particles and number of orbtials
  void resize(int nup, int ndown);

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  PsiValueType ratio(ParticleSet& P, int iat) override;

  void ratioUp(ParticleSet& P, int iat);

  void ratioDown(ParticleSet& P, int iat);

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  void resizeByWalkers(int nwalkers);

  /** Calculate the log value of the Dirac determinant for particles
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@return the value of the determinant
   *
   *\f$ (first,first+nel). \f$  Add the gradient and laplacian
   *contribution of the determinant to G(radient) and L(aplacian)
   *for local energy calculations.
   */
  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  ///Total number of particles
  int NumPtcls;
  ///number of major spins
  int Nup;
  ///number of minor spins
  int Ndown;
  ///size of the basis set
  int BasisSize;

  ///index of the particle (or row)
  int WorkingIndex;

  /////Current determinant value
  //ValueType CurrentDet;

  ///coefficient of the up/down block
  ValueMatrix Lambda;

  ///coefficient of the major block
  ValueMatrix LambdaUP;

  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix psiM, psiM_temp;


  /**  Transient data for gradient and laplacian evaluation
   *
   * \f$phiD(j,k) = \sum_{j^{'}} \lambda_{j^{'},j} \phi_k(r_j) \f$
   * j runs over the particle index index
   */
  ValueMatrix phiT;

  /// temporary container for testing
  ValueMatrix psiMinv;
  /// store gradients
  GradMatrix dY;
  /// store laplacians
  ValueMatrix d2Y;
  /// temporary determinant-related matrix for gradients
  GradMatrix dpsiU, dpsiD;
  /// temporary determinant-related matrix for laplacians
  ValueMatrix d2psiU, d2psiD;

  /// value of single-particle orbital for particle-by-particle update
  /** temporary vector for a particle-by-particle move
   *
   * phiTv = Lambda Y(iat)
   */
  ValueVector phiTv;
  ValueVector psiU, psiD;
  GradVector dpsiUv, dpsiDv;
  ValueVector d2psiUv, d2psiDv;
  ValueVector workV1, workV2;
  ValueVector WorkSpace;
  IndexVector Pivot;

  ///current ratio
  PsiValueType curRatio;
  ///cummulate ratio for particle-by-particle update
  RealType cumRatio;
  ///address of  dpsiU[0][0]
  BasisSetType::ValueType* FirstAddressOfdVU;
  ///address of FirstAddressOfdVU+OHMMS_DIM*Nup*Nup
  BasisSetType::ValueType* LastAddressOfdVU;
  ///address of  dpsiD[0][0]
  BasisSetType::ValueType* FirstAddressOfdVD;
  ///address of FirstAddressOfdVD+OHMMS_DIM*Ndown*Nup
  BasisSetType::ValueType* LastAddressOfdVD;
  ///address of myG[0][0]
  ParticleSet::SingleParticleValue* FirstAddressOfG;
  ///address of FirstAddressOfG+OHMMS_DIM*NumPtcls
  ParticleSet::SingleParticleValue* LastAddressOfG;
  ///address of dY[0][0]
  BasisSetType::ValueType* FirstAddressOfdY;
  ///address of FirstAddressOfdY+NumPtcls*BasisSize
  BasisSetType::ValueType* LastAddressOfdY;

  ParticleSet::ParticleGradient myG, myG_temp;
  ParticleSet::ParticleLaplacian myL, myL_temp;

  void evaluateLogAndStore(const ParticleSet& P);
};
} // namespace qmcplusplus
#endif
