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

template<typename T>
ConstantSPOSet<T>::ConstantSPOSet(const std::string& my_name, const int nparticles, const int norbitals)
    : SPOSetT<T>(my_name), numparticles_(nparticles)
{
  SPOSet::OrbitalSetSize = norbitals;
  ref_psi_.resize(numparticles_, SPOSet::OrbitalSetSize);
  ref_egrad_.resize(numparticles_, SPOSet::OrbitalSetSize);
  ref_elapl_.resize(numparticles_, SPOSet::OrbitalSetSize);

  ref_psi_   = 0.0;
  ref_egrad_ = 0.0;
  ref_elapl_ = 0.0;
};

template<typename T>
std::unique_ptr<SPOSetT<T>> ConstantSPOSet<T>::makeClone() const
{
  auto myclone = std::make_unique<ConstantSPOSet>(SPOSet::my_name_, numparticles_, SPOSet::OrbitalSetSize);
  myclone->setRefVals(ref_psi_);
  myclone->setRefEGrads(ref_egrad_);
  myclone->setRefELapls(ref_elapl_);
  return myclone;
};

template<typename T>
std::string ConstantSPOSet<T>::getClassName() const
{
  return "ConstantSPOSet";
};

template<typename T>
void ConstantSPOSet<T>::checkOutVariables(const OptVariables& active)
{
  APP_ABORT("ConstantSPOSet should not call checkOutVariables");
};

template<typename T>
void ConstantSPOSet<T>::setOrbitalSetSize(int norbs)
{
  APP_ABORT("ConstantSPOSet should not call setOrbitalSetSize()");
}

template<typename T>
void ConstantSPOSet<T>::setRefVals(const ValueMatrix& vals)
{
  assert(vals.cols() == SPOSet::OrbitalSetSize);
  assert(vals.rows() == numparticles_);
  ref_psi_ = vals;
};

template<typename T>
void ConstantSPOSet<T>::setRefEGrads(const GradMatrix& grads)
{
  assert(grads.cols() == SPOSet::OrbitalSetSize);
  assert(grads.rows() == numparticles_);
  ref_egrad_ = grads;
};

template<typename T>
void ConstantSPOSet<T>::setRefELapls(const ValueMatrix& lapls)
{
  assert(lapls.cols() == SPOSet::OrbitalSetSize);
  assert(lapls.rows() == numparticles_);
  ref_elapl_ = lapls;
};

template<typename T>
void ConstantSPOSet<T>::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  const auto* vp = dynamic_cast<const VirtualParticleSet*>(&P);
  int ptcl       = vp ? vp->refPtcl : iat;
  assert(psi.size() == SPOSet::OrbitalSetSize);
  for (int iorb = 0; iorb < SPOSet::OrbitalSetSize; iorb++)
    psi[iorb] = ref_psi_(ptcl, iorb);
};

template<typename T>
void ConstantSPOSet<T>::evaluateVGL(const ParticleSet& P,
                                    int iat,
                                    ValueVector& psi,
                                    GradVector& dpsi,
                                    ValueVector& d2psi)
{
  for (int iorb = 0; iorb < SPOSet::OrbitalSetSize; iorb++)
  {
    psi[iorb]   = ref_psi_(iat, iorb);
    dpsi[iorb]  = ref_egrad_(iat, iorb);
    d2psi[iorb] = ref_elapl_(iat, iorb);
  }
};

template<typename T>
void ConstantSPOSet<T>::evaluate_notranspose(const ParticleSet& P,
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

#if !defined(MIXED_PRECISION)
template class ConstantSPOSet<double>;
template class ConstantSPOSet<std::complex<double>>;
#endif
template class ConstantSPOSet<float>;
template class ConstantSPOSet<std::complex<float>>;


} //namespace qmcplusplus
