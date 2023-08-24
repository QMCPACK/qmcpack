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

#include "QMCWaveFunctions/SPOSetT.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "LCAOrbitalSetT.h"
#include "SoaCuspCorrection.h"


namespace qmcplusplus
{
/** class to add cusp correction to LCAOrbitalSet.
   *
   */

template<typename T>
class LCAOrbitalSetWithCorrectionT : public SPOSetT<T>
{
public:
  using basis_type  = typename LCAOrbitalSetT<T>::basis_type;
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
  using GradMatrix  = typename SPOSetT<T>::GradMatrix;
  /** constructor
     * @param ions
     * @param els
     * @param bs pointer to the BasisSet
     * @param rl report level
     */
  LCAOrbitalSetWithCorrectionT(const std::string& my_name,
                               ParticleSet& ions,
                               ParticleSet& els,
                               std::unique_ptr<basis_type>&& bs);

  LCAOrbitalSetWithCorrectionT(const LCAOrbitalSetWithCorrectionT& in) = default;

  std::string getClassName() const final { return "LCAOrbitalSetWithCorrectionT"; }

  std::unique_ptr<SPOSetT<T>> makeClone() const final;

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
  LCAOrbitalSetT<T> lcao;

  SoaCuspCorrection cusp;
};
} // namespace qmcplusplus
#endif
