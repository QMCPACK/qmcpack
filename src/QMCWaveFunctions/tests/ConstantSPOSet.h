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
/** Constant SPOSet for testing purposes.  Fixed N_elec x N_orb matrices storing value, gradients, and laplacians.
   *
   */
struct ConstantSPOSet : public SPOSet
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

  void setRefVals(const ValueMatrix& vals);
  void setRefEGrads(const GradMatrix& grads);
  void setRefELapls(const ValueMatrix& lapls);

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

protected:
private:
  int numparticles_;
  
  //Value, electron gradient, and electron laplacian at "reference configuration".  
  //i.e. before any attempted moves.
  
  ValueMatrix ref_psi_;
  GradMatrix ref_egrad_;
  ValueMatrix ref_elapl_;
};
} // namespace qmcplusplus
#endif
