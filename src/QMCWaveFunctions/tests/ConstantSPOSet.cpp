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

#include "QMCWaveFunctions/tests/ConstantSPOSet.h"

namespace qmcplusplus
{
ConstantSPOSet::ConstantSPOSet(const std::string& my_name, const int nparticles, const int norbitals)
    : SPOSet(my_name), numparticles_(nparticles)
{
  OrbitalSetSize = norbitals;
  ref_psi_.resize(numparticles_, OrbitalSetSize);
  ref_egrad_.resize(numparticles_, OrbitalSetSize);
  ref_elapl_.resize(numparticles_, OrbitalSetSize);

  ref_psi_   = 0.0;
  ref_egrad_ = 0.0;
  ref_elapl_ = 0.0;
};

std::unique_ptr<SPOSet> ConstantSPOSet::makeClone() const
{
  auto myclone = std::make_unique<ConstantSPOSet>(my_name_, numparticles_, OrbitalSetSize);
  myclone->setRefVals(ref_psi_);
  myclone->setRefEGrads(ref_egrad_);
  myclone->setRefELapls(ref_elapl_);
  return myclone;
};

std::string ConstantSPOSet::getClassName() const { return "ConstantSPOSet"; };

void ConstantSPOSet::checkOutVariables(const opt_variables_type& active)
{
  APP_ABORT("ConstantSPOSet should not call checkOutVariables");
};

void ConstantSPOSet::setOrbitalSetSize(int norbs) { APP_ABORT("ConstantSPOSet should not call setOrbitalSetSize()"); }

void ConstantSPOSet::setRefVals(const ValueMatrix& vals)
{
  assert(vals.cols() == OrbitalSetSize);
  assert(vals.rows() == numparticles_);
  ref_psi_ = vals;
};
void ConstantSPOSet::setRefEGrads(const GradMatrix& grads)
{
  assert(grads.cols() == OrbitalSetSize);
  assert(grads.rows() == numparticles_);
  ref_egrad_ = grads;
};
void ConstantSPOSet::setRefELapls(const ValueMatrix& lapls)
{
  assert(lapls.cols() == OrbitalSetSize);
  assert(lapls.rows() == numparticles_);
  ref_elapl_ = lapls;
};

void ConstantSPOSet::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  const auto* vp = dynamic_cast<const VirtualParticleSet*>(&P);
  int ptcl = vp ? vp->refPtcl : iat;
  assert(psi.size() == OrbitalSetSize);
  for (int iorb = 0; iorb < OrbitalSetSize; iorb++)
    psi[iorb] = ref_psi_(ptcl, iorb);
};

void ConstantSPOSet::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  for (int iorb = 0; iorb < OrbitalSetSize; iorb++)
  {
    psi[iorb]   = ref_psi_(iat, iorb);
    dpsi[iorb]  = ref_egrad_(iat, iorb);
    d2psi[iorb] = ref_elapl_(iat, iorb);
  }
};

void ConstantSPOSet::evaluate_notranspose(const ParticleSet& P,
                                          int first,
                                          int last,
                                          ValueMatrix& logdet,
                                          GradMatrix& dlogdet,
                                          ValueMatrix& d2logdet)
{
  for (int iat = first, i = 0; iat < last; ++iat, ++i)
  {
    ValueVector v(logdet[i], logdet.cols());
    GradVector g(dlogdet[i], dlogdet.cols());
    ValueVector l(d2logdet[i], d2logdet.cols());
    evaluateVGL(P, iat, v, g, l);
  }
}
} //namespace qmcplusplus
