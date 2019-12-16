//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineSet.h
 *
 * BsplineSet<SplineAdoptor> is a SPOSet class to work with determinant classes
 */
#ifndef QMCPLUSPLUS_EINSPLINE_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_ADOPTOR_H

#include <Lattice/CrystalLattice.h>
#include <spline/einspline_engine.hpp>
#include <spline/einspline_util.hpp>
#include <simd/allocator.hpp>

namespace qmcplusplus
{
/** BsplineSet<SplineAdoptor>, a SPOSet
 * @tparam SplineAdoptor implements evaluation functions that matched the storage requirements.
 *
 * Equivalent to EinsplineSetExtended<Storage>
 * Storage is now handled by SplineAdoptor class that is specialized for precision, storage etc.
 * @todo Make SplineAdoptor be a member not the base class. This is needed
 * to make MultiBsplineSet (TBD) which has multiple SplineAdoptors for distributed
 * cases.
 */
struct BsplineSet : public SPOSet
{
  // propagate SPOSet virtual functions
  using SPOSet::finalizeConstruction;
  using SPOSet::evaluateValue;
  using SPOSet::evaluateDetRatios;
  using SPOSet::evaluateVGL;
  using SPOSet::mw_evaluateVGL;
  using SPOSet::evaluateVGH;
  using SPOSet::evaluateVGHGH;

  virtual SPOSet* makeClone() const override = 0;

  void resetParameters(const opt_variables_type& active) override {}

  void resetTargetParticleSet(ParticleSet& e) override {}

  void setOrbitalSetSize(int norbs) override
  {
    OrbitalSetSize = norbs;
    //SplineAdoptor::first_spo=0;
    //SplineAdoptor::last_spo=norbs;
  }

  virtual void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet) override
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector_t v(logdet[i], OrbitalSetSize);
      GradVector_t g(dlogdet[i], OrbitalSetSize);
      ValueVector_t l(d2logdet[i], OrbitalSetSize);
      evaluateVGL(P, iat, v, g, l);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    HessMatrix_t& grad_grad_logdet) override
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    typedef HessMatrix_t::value_type hess_type;
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector_t v(logdet[i], OrbitalSetSize);
      GradVector_t g(dlogdet[i], OrbitalSetSize);
      HessVector_t h(grad_grad_logdet[i], OrbitalSetSize);
      evaluateVGH(P, iat, v, g, h);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    HessMatrix_t& grad_grad_logdet,
                                    GGGMatrix_t& grad_grad_grad_logdet) override
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    typedef HessMatrix_t::value_type hess_type;
    typedef GGGMatrix_t::value_type ghess_type;
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector_t v(logdet[i], OrbitalSetSize);
      GradVector_t g(dlogdet[i], OrbitalSetSize);
      HessVector_t h(grad_grad_logdet[i], OrbitalSetSize);
      GGGVector_t gh(grad_grad_grad_logdet[i], OrbitalSetSize);
      evaluateVGHGH(P, iat, v, g, h, gh);
    }
  }

  virtual void evaluateGradSource(const ParticleSet& P,
                                  int first,
                                  int last,
                                  const ParticleSet& source,
                                  int iat_src,
                                  GradMatrix_t& gradphi) override
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }

  virtual void evaluateGradSource(const ParticleSet& P,
                                  int first,
                                  int last,
                                  const ParticleSet& source,
                                  int iat_src,
                                  GradMatrix_t& grad_phi,
                                  HessMatrix_t& grad_grad_phi,
                                  GradMatrix_t& grad_lapl_phi) override
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }
};

} // namespace qmcplusplus
#endif
