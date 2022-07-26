//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONSTANTSPOSET_H
#define QMCPLUSPLUS_CONSTANTSPOSET_H

#include <memory>
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"

#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{
/** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   * SoA verson of LCOrtbitalSet
   * Localized basis set is always real 
   */
struct ConstantSPOSet : public SPOSet
{
public:

  ConstantSPOSet(){};
  ConstantSPOSet(ValueMatrix& vals)
  {
    psi_=vals; 
    OrbitalSetSize=psi_.cols();
  }
  std::unique_ptr<SPOSet> makeClone() const override{ return std::make_unique<ConstantSPOSet>();};


  void checkInVariables(opt_variables_type& active) override
  {
    APP_ABORT("ConstantSPOSet should not call checkInVariables");
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    APP_ABORT("ConstantSPOSet should not call checkOutVariables");
  }

  ///reset
  void resetParameters(const opt_variables_type& active) override
  {
    APP_ABORT("ConstantSPOSet should not call resetParameters");
  }

  /** set the OrbitalSetSize and Identity=false and initialize internal storages
    */
  void setOrbitalSetSize(int norbs) override 
  {
    APP_ABORT("ConstantSPOSet should not call setOrbitalSetSize()");
  }


  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override 
  {
    if (psi.size() != OrbitalSetSize)
      APP_ABORT("Borked");
    for(int iorb=0; iorb<OrbitalSetSize; iorb++)
      psi[iorb]=psi_(iat,iorb);
  };

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override 
  {
    for(int iorb=0; iorb<OrbitalSetSize; iorb++)
      psi[iorb]=psi_(iat,iorb);
    dpsi=0.0;
    d2psi=0.0;
  };
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override 
  {
    logdet=psi_ ;
  };



protected:
  ///number of Single-particle orbitals


private:
  ValueMatrix psi_;
};
} // namespace qmcplusplus
#endif
