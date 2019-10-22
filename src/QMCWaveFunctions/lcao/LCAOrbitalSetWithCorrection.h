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
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include "QMCWaveFunctions/lcao/SoaCuspCorrectionBasisSet.h"


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
  LCAOrbitalSetWithCorrection(ParticleSet& ions, ParticleSet& els, basis_type* bs, bool optimize);

  LCAOrbitalSetWithCorrection(const LCAOrbitalSetWithCorrection& in) = default;

  SPOSet* makeClone() const;

  void setOrbitalSetSize(int norbs);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet,
                            GGGMatrix_t& grad_grad_grad_logdet);

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet);
};
} // namespace qmcplusplus
#endif
