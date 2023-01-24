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

  ConstantSPOSet(const std::string& my_name, const ValueMatrix& vals) : SPOSet(my_name)
  {
    const int nrows = vals.rows();
    const int ncols = vals.cols();
    OrbitalSetSize  = ncols;

    grad_.resize(nrows, ncols);
    lapl_.resize(nrows, ncols);

    psi_  = vals;
    grad_ = 0.0;
    lapl_ = 0.0;
  };

  ConstantSPOSet(const std::string& my_name, const ValueMatrix& vals, const GradMatrix& grads) : SPOSet(my_name)
  {
    const int nrows = vals.rows();
    const int ncols = vals.cols();
    OrbitalSetSize  = ncols;

    assert(grads.rows() == nrows);
    assert(grads.cols() == ncols);

    lapl_.resize(nrows, ncols);

    psi_  = vals;
    grad_ = grads;
    lapl_ = 0.0;
  }

  ConstantSPOSet(const std::string& my_name, const ValueMatrix& vals, const GradMatrix& grads, const ValueMatrix& lapls)
      : SPOSet(my_name)
  {
    const int nrows = vals.rows();
    const int ncols = vals.cols();
    OrbitalSetSize  = ncols;

    assert(grads.rows() == nrows);
    assert(grads.cols() == ncols);
    assert(lapls.rows() == nrows);
    assert(lapls.cols() == ncols);

    psi_  = vals;
    grad_ = grads;
    lapl_ = lapls;
  }


  std::unique_ptr<SPOSet> makeClone() const override
  {
    auto myclone = std::make_unique<ConstantSPOSet>(my_name_, psi_, grad_, lapl_);
    return myclone;
  };

  std::string getClassName() const override { return "ConstantSPOSet"; }

  void checkOutVariables(const opt_variables_type& active) override
  {
    APP_ABORT("ConstantSPOSet should not call checkOutVariables");
  }

  void setOrbitalSetSize(int norbs) override { APP_ABORT("ConstantSPOSet should not call setOrbitalSetSize()"); }


  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override
  {
    assert(psi.size() == OrbitalSetSize);
    for (int iorb = 0; iorb < OrbitalSetSize; iorb++)
      psi[iorb] = psi_(iat, iorb);
  };

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override
  {
    for (int iorb = 0; iorb < OrbitalSetSize; iorb++)
    {
      psi[iorb]   = psi_(iat, iorb);
      dpsi[iorb]  = grad_(iat, iorb);
      d2psi[iorb] = lapl_(iat, iorb);
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
  ValueMatrix psi_;
  GradMatrix grad_;
  ValueMatrix lapl_;
};
} // namespace qmcplusplus
#endif
