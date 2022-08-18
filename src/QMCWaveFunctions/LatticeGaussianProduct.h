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
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{
/** A composite Orbital
 */
class LatticeGaussianProduct : public WaveFunctionComponent
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

  ~LatticeGaussianProduct() override;

  std::string getClassName() const override { return "LatticeGaussianProduct"; }

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  PsiValueType ratio(ParticleSet& P, int iat) override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  void restore(int iat) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override
  {}

  void evaluateLogAndStore(const ParticleSet& P, ParticleSet::ParticleGradient& dG, ParticleSet::ParticleLaplacian& dL);
};
} // namespace qmcplusplus
#endif
