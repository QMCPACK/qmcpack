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
#include "LCAOrbitalSet.h"
#include "SoaCuspCorrection.h"


namespace qmcplusplus
{
/** class to add cusp correction to LCAOrbitalSet.
   *
   */
class LCAOrbitalSetWithCorrection : public SPOSet
{
public:
  using basis_type = LCAOrbitalSet::basis_type;
  /** constructor
     * @param ions
     * @param els
     * @param bs pointer to the BasisSet
     * @param rl report level
     */
  LCAOrbitalSetWithCorrection(const std::string& my_name,
                              ParticleSet& ions,
                              ParticleSet& els,
                              std::unique_ptr<basis_type>&& bs);

  LCAOrbitalSetWithCorrection(const LCAOrbitalSetWithCorrection& in) = default;

  std::string getClassName() const final { return "LCAOrbitalSetWithCorrection"; }

  std::unique_ptr<SPOSet> makeClone() const final;

  void setOrbitalSetSize(int norbs) final;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) final;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) final;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) final;

  friend class LCAOrbitalBuilder;

private:
  LCAOrbitalSet lcao;

  SoaCuspCorrection cusp;
};
} // namespace qmcplusplus
#endif
