//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 Raymond Clay and QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONSTANTSPOSET_H
#define QMCPLUSPLUS_CONSTANTSPOSET_H

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{
/** Constant SPOSet for testing purposes.  Fixed N_elec x N_orb matrices storing value, gradients, and laplacians, etc.,
   *  These values are accessed through standard SPOSet calls like evaluateValue, evaluateVGL, etc.
   *  Exists to provide deterministic and known output to objects requiring SPOSet evaluations.      
   *
   */
class ConstantSPOSet : public SPOSet
{
public:
  ConstantSPOSet(const std::string& my_name) = delete;

  //Constructor needs number of particles and number of orbitals.  This is the minimum
  //amount of information needed to sanely construct all data members and perform size
  //checks later.
  ConstantSPOSet(const std::string& my_name, const int nparticles, const int norbitals);

  std::unique_ptr<SPOSet> makeClone() const override;

  std::string getClassName() const override;

  void checkOutVariables(const opt_variables_type& active) override;

  void setOrbitalSetSize(int norbs) override;

  /**
  * @brief Setter method to set \phi_j(r_i). Stores input matrix in ref_psi_.
  * @param Nelec x Nion ValueType matrix of \phi_j(r_i)
  * @return void
  */ 
  void setRefVals(const ValueMatrix& vals);
  /**
  * @brief Setter method to set \nabla_i \phi_j(r_i). Stores input matrix in ref_egrad_.
  * @param Nelec x Nion GradType matrix of \grad_i \phi_j(r_i)
  * @return void
  */ 
  void setRefEGrads(const GradMatrix& grads);
  /**
  * @brief Setter method to set \nabla^2_i \phi_j(r_i). Stores input matrix in ref_elapl_.
  * @param Nelec x Nion GradType matrix of \grad^2_i \phi_j(r_i)
  * @return void
  */ 
  void setRefELapls(const ValueMatrix& lapls);

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

private:
  const int numparticles_; /// evaluate_notranspose arrays are nparticle x norb matrices.
                           /// To ensure consistent array sizing and enforcement,
                           /// we agree at construction how large these matrices will be.
                           /// norb is stored in SPOSet::OrbitalSetSize.  

  //Value, electron gradient, and electron laplacian at "reference configuration".
  //i.e. before any attempted moves.

  ValueMatrix ref_psi_;
  GradMatrix ref_egrad_;
  ValueMatrix ref_elapl_;
};
} // namespace qmcplusplus
#endif
