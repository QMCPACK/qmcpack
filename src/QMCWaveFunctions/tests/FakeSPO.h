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

class FakeSPO : public SPOSet
{
public:
  Matrix<ValueType> a;
  Matrix<ValueType> a2;
  Vector<ValueType> v;
  Matrix<ValueType> v2;

  SPOSet::GradVector_t gv;
  
  FakeSPO();
  virtual ~FakeSPO() {}
  virtual SPOSet* makeClone() const{ return new FakeSPO; }
  virtual void report() {}
  virtual void resetParameters(const opt_variables_type& optVariables) {}
  virtual void setOrbitalSetSize(int norbs);

  virtual void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi);

  virtual void evaluateVGL(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    ValueMatrix_t& d2logdet);
};

} // namespace qmcplusplus
#endif
