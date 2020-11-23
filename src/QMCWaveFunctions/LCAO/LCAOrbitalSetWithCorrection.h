//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_WITH_CORRECTION_TEMP_H
#define QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_WITH_CORRECTION_TEMP_H

#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/LCAO/LCAOrbitalSet.h"
#include "QMCWaveFunctions/LCAO/SoaCuspCorrectionBasisSet.h"


namespace qmcplusplus
{
/** class to add cusp correction to LCAOrbitalSet.
   *
   */
struct LCAOrbitalSetWithCorrection : public LCAOrbitalSet
{
  SoaCuspCorrection cusp;

  /** constructor
     * @param ions
     * @param els
     * @param bs pointer to the BasisSet
     * @param rl report level
     */
  LCAOrbitalSetWithCorrection(ParticleSet& ions, ParticleSet& els, std::unique_ptr<basis_type>&& bs, bool optimize);

  LCAOrbitalSetWithCorrection(const LCAOrbitalSetWithCorrection& in) = default;

  SPOSet* makeClone() const override;

  void setOrbitalSetSize(int norbs) override;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi) override;

  void evaluateVGL(const ParticleSet& P,
                   int iat,
                   ValueVector_t& psi,
                   GradVector_t& dpsi,
                   ValueVector_t& d2psi) override;

  void evaluateVGH(const ParticleSet& P,
                   int iat,
                   ValueVector_t& psi,
                   GradVector_t& dpsi,
                   HessVector_t& grad_grad_psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet,
                            GGGMatrix_t& grad_grad_grad_logdet) override;

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet) override;
};
} // namespace qmcplusplus
#endif
