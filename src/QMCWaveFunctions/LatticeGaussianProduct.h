//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file LatticeGaussianProduct.h
 * @brief Simple gaussian functions used for orbitals for ions
 */
#ifndef QMCPLUSPLUS_LATTICE_GAUSSIAN_PRODUCT
#define QMCPLUSPLUS_LATTICE_GAUSSIAN_PRODUCT
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
/** A composite Orbital
 */
struct LatticeGaussianProduct : public WaveFunctionComponent
{
private:
  ParticleAttrib<RealType> U, d2U;
  ParticleAttrib<PosType> dU;
  RealType *FirstAddressOfdU, *LastAddressOfdU;
  ///table index
  int myTableID;
  ///orbital centers
  ParticleSet& CenterRef;
  int NumTargetPtcls, NumCenters;
  RealType curVal, curLap;
  PosType curGrad;

public:
  std::vector<RealType> ParticleAlpha;
  std::vector<int> ParticleCenter;

  LatticeGaussianProduct(ParticleSet& centers, ParticleSet& ptcls);

  ~LatticeGaussianProduct();

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& o);

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& o);

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os);

  /** reset the parameters during optimizations
   */
  void resetParameters(const opt_variables_type& active);

  void resetTargetParticleSet(ParticleSet& P);

  RealType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  ValueType ratio(ParticleSet& P, int iat);

  void acceptMove(ParticleSet& P, int iat);

  void restore(int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  GradType evalGrad(ParticleSet& P, int iat);

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);


  WaveFunctionComponent* makeClone(ParticleSet& tqp) const;

  void evaluateLogAndStore(ParticleSet& P, ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL);
};
} // namespace qmcplusplus
#endif
