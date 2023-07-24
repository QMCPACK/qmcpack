//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FAKESPOTT_H
#define QMCPLUSPLUS_FAKESPOTT_H

#include "QMCWaveFunctions/SPOSetT.h"

namespace qmcplusplus
{
template<class T>
class FakeSPOT : public SPOSetT<T>
{
public:
  Matrix<T> a;
  Matrix<T> a2;
  Vector<T> v;
  Matrix<T> v2;

  using ValueVector = typename SPOSetT<T>::ValueVector;
  using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using GradMatrix  = typename SPOSetT<T>::GradMatrix;
  using GradType    = typename SPOSetT<T>::GradType;

  typename SPOSetT<T>::GradVector gv;

  FakeSPOT();

  ~FakeSPOT() override = default;

  std::string getClassName() const override { return "FakeSPO"; }

  std::unique_ptr<SPOSetT<T>> makeClone() const override;

  virtual void report() {}

  void setOrbitalSetSize(int norbs) override;

  void evaluateValue(const ParticleSetT<T>& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSetT<T>& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSetT<T>& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;
};

} // namespace qmcplusplus
#endif
