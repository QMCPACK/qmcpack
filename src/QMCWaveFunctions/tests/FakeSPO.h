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


#ifndef QMCPLUSPLUS_FAKESPO_H
#define QMCPLUSPLUS_FAKESPO_H

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{

template<typename VALUE>
class FakeSPO : public SPOSetT<VALUE>
{
public:
  enum
  {
    DIM     = OHMMS_DIM,
  };
  using SPOSet             = SPOSetT<VALUE>;
  using ValueType          = typename SPOSet::ValueType;
  using GradType           = typename SPOSet::GradType;
  using ValueVector        = typename SPOSet::ValueVector;
  using GradVector         = typename SPOSet::GradVector;
  using ValueMatrix        = typename SPOSet::ValueMatrix;
  using GradMatrix         = typename SPOSet::GradMatrix;

  Matrix<ValueType> a;
  Matrix<ValueType> a2;
  Vector<ValueType> v;
  Matrix<ValueType> v2;

  GradVector gv;

  FakeSPO();
  ~FakeSPO() override {}

  std::string getClassName() const override { return "FakeSPO"; }

  std::unique_ptr<SPOSet> makeClone() const override;
  virtual void report() {}
  void setOrbitalSetSize(int norbs) override;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

private:
  using SPOSet::OrbitalSetSize;
};

extern template class FakeSPO<QMCTraits::QTBase::RealType>;
extern template class FakeSPO<QMCTraits::QTBase::ComplexType>;

} // namespace qmcplusplus
#endif
