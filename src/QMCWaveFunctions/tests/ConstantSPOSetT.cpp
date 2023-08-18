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

#include "ConstantSPOSetT.h"

namespace qmcplusplus
{

template<class T>
ConstantSPOSetT<T>::ConstantSPOSetT(const std::string& my_name, const int nparticles, const int norbitals)
    : SPOSetT<T>(my_name), numparticles_(nparticles)
{
  this->OrbitalSetSize = norbitals;
  ref_psi_.resize(numparticles_, this->OrbitalSetSize);
  ref_egrad_.resize(numparticles_, this->OrbitalSetSize);
  ref_elapl_.resize(numparticles_, this->OrbitalSetSize);

  ref_psi_   = 0.0;
  ref_egrad_ = 0.0;
  ref_elapl_ = 0.0;
}

template<class T>
std::unique_ptr<SPOSetT<T>> ConstantSPOSetT<T>::makeClone() const
{
  auto myclone = std::make_unique<ConstantSPOSetT<T>>(this->my_name_, numparticles_, this->OrbitalSetSize);
  myclone->setRefVals(ref_psi_);
  myclone->setRefEGrads(ref_egrad_);
  myclone->setRefELapls(ref_elapl_);
  return myclone;
}

template<class T>
void ConstantSPOSetT<T>::checkOutVariables(const opt_variables_type& active)
{
  APP_ABORT("ConstantSPOSet should not call checkOutVariables");
};

template<class T>
void ConstantSPOSetT<T>::setOrbitalSetSize(int norbs)
{
  APP_ABORT("ConstantSPOSet should not call setOrbitalSetSize()");
}

template<class T>
void ConstantSPOSetT<T>::setRefVals(const ValueMatrix& vals)
{
  assert(vals.cols() == this->OrbitalSetSize);
  assert(vals.rows() == numparticles_);
  ref_psi_ = vals;
}

template<class T>
void ConstantSPOSetT<T>::setRefEGrads(const GradMatrix& grads)
{
  assert(grads.cols() == this->OrbitalSetSize);
  assert(grads.rows() == numparticles_);
  ref_egrad_ = grads;
}

template<class T>
void ConstantSPOSetT<T>::setRefELapls(const ValueMatrix& lapls)
{
  assert(lapls.cols() == this->OrbitalSetSize);
  assert(lapls.rows() == numparticles_);
  ref_elapl_ = lapls;
}

template<class T>
void ConstantSPOSetT<T>::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  const auto* vp = dynamic_cast<const VirtualParticleSet*>(&P);
  int ptcl       = vp ? vp->refPtcl : iat;
  assert(psi.size() == this->OrbitalSetSize);
  for (int iorb = 0; iorb < this->OrbitalSetSize; iorb++)
    psi[iorb] = ref_psi_(ptcl, iorb);
}

template<class T>
void ConstantSPOSetT<T>::evaluateVGL(const ParticleSet& P,
                                     int iat,
                                     ValueVector& psi,
                                     GradVector& dpsi,
                                     ValueVector& d2psi)
{
  for (int iorb = 0; iorb < this->OrbitalSetSize; iorb++)
  {
    psi[iorb]   = ref_psi_(iat, iorb);
    dpsi[iorb]  = ref_egrad_(iat, iorb);
    d2psi[iorb] = ref_elapl_(iat, iorb);
  }
}

template<class T>
void ConstantSPOSetT<T>::evaluate_notranspose(const ParticleSet& P,
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

template class ConstantSPOSetT<float>;
template class ConstantSPOSetT<double>;
template class ConstantSPOSetT<std::complex<float>>;
template class ConstantSPOSetT<std::complex<double>>;

} //namespace qmcplusplus
