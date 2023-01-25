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
  ConstantSPOSet(const std::string& my_name, const int nparticles, const int norbitals): SPOSet(my_name), numparticles_(nparticles)
  {
    OrbitalSetSize = norbitals; 
    ref_psi_.resize(numparticles_,OrbitalSetSize);
    ref_egrad_.resize(numparticles_,OrbitalSetSize);
    ref_elapl_.resize(numparticles_,OrbitalSetSize);
    
    ref_psi_ = 0.0;
    ref_egrad_ = 0.0;
    ref_elapl_ = 0.0;
    
  }

  std::unique_ptr<SPOSet> makeClone() const override
  {
    auto myclone = std::make_unique<ConstantSPOSet>(my_name_, numparticles_, OrbitalSetSize);
    myclone->setRefVals(ref_psi_);
    myclone->setRefEGrads(ref_egrad_);
    myclone->setRefELapls(ref_elapl_);
    return myclone;
  };

  std::string getClassName() const override { return "ConstantSPOSet"; }

  void checkOutVariables(const opt_variables_type& active) override
  {
    APP_ABORT("ConstantSPOSet should not call checkOutVariables");
  }

  void setOrbitalSetSize(int norbs) override { APP_ABORT("ConstantSPOSet should not call setOrbitalSetSize()"); }

  void setRefVals(const ValueMatrix& vals)
  { 
    assert(vals.cols()==OrbitalSetSize);
    assert(vals.rows()==numparticles_);
    ref_psi_ = vals; 
  };
  void setRefEGrads(const GradMatrix& grads)
  { 
    assert(grads.cols()==OrbitalSetSize);
    assert(grads.rows()==numparticles_);
    ref_egrad_ = grads; 
  };
  void setRefELapls(const ValueMatrix& lapls)
  { 
    assert(lapls.cols()==OrbitalSetSize);
    assert(lapls.rows()==numparticles_);
    ref_elapl_ = lapls; 
  };

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override
  {
    assert(psi.size() == OrbitalSetSize);
    for (int iorb = 0; iorb < OrbitalSetSize; iorb++)
      psi[iorb] = ref_psi_(iat, iorb);
  };

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override
  {
    for (int iorb = 0; iorb < OrbitalSetSize; iorb++)
    {
      psi[iorb]   = ref_psi_(iat, iorb);
      dpsi[iorb]  = ref_egrad_(iat, iorb);
      d2psi[iorb] = ref_elapl_(iat, iorb);
    }
  };

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override
  {
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector v(logdet[i], logdet.cols());
      GradVector g(dlogdet[i], dlogdet.cols());
      ValueVector l(d2logdet[i], d2logdet.cols());
      evaluateVGL(P, iat, v, g, l);
    }
  }


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
