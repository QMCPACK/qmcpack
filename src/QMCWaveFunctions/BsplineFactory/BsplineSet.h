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
template<typename SplineAdoptor>
struct BsplineSet: public SPOSet, public SplineAdoptor
{
  typedef typename SplineAdoptor::SplineType SplineType;
  typedef typename SplineAdoptor::PointType  PointType;
  typedef typename SplineAdoptor::DataType  DataType;

  ///** default constructor */
  //BsplineSet() { }

  SPOSet* makeClone() const
  {
    return new BsplineSet<SplineAdoptor>(*this);
  }

  /** set_spline to the big table
   * @param psi_r starting address of real part of psi(ispline)
   * @param psi_i starting address of imaginary part of psi(ispline)
   * @param twist twist id, reserved to sorted adoptor, ignored
   * @param ispline index of this spline function
   * @param level refinement level
   *
   * Each adoptor handles the map with respect to the twist, state index and refinement level
   */
  template<typename CT>
  void set_spline(CT* spline_r, CT* spline_i, int twist, int ispline, int level)
  {
    SplineAdoptor::set_spline(spline_r,spline_i,twist,ispline,level);
  }

  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    SplineAdoptor::evaluate_v(P,iat,psi);
  }

  inline void evaluateDetRatios(const VirtualParticleSet& VP, ValueVector_t& psi, const ValueVector_t& psiinv, std::vector<ValueType>& ratios)
  {
    assert(psi.size() == psiinv.size());
    SplineAdoptor::evaluateDetRatios(VP, psi, psiinv, ratios);
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    SplineAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
  {
    SplineAdoptor::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
  }

  void resetParameters(const opt_variables_type& active)
  { }

  void resetTargetParticleSet(ParticleSet& e)
  { }

  void setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
    //SplineAdoptor::first_spo=0;
    //SplineAdoptor::last_spo=norbs;
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      ValueVector_t v(logdet[i],OrbitalSetSize);
      GradVector_t  g(dlogdet[i],OrbitalSetSize);
      ValueVector_t l(d2logdet[i],OrbitalSetSize);
      SplineAdoptor::evaluate_vgl(P,iat,v,g,l);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    typedef HessMatrix_t::value_type hess_type;
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      ValueVector_t v(logdet[i],OrbitalSetSize);
      GradVector_t  g(dlogdet[i],OrbitalSetSize);
      HessVector_t  h(grad_grad_logdet[i],OrbitalSetSize);
      SplineAdoptor::evaluate_vgh(P,iat,v,g,h);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    typedef HessMatrix_t::value_type hess_type;
    typedef GGGMatrix_t::value_type  ghess_type;
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      ValueVector_t v(logdet[i],OrbitalSetSize);
      GradVector_t  g(dlogdet[i],OrbitalSetSize);
      HessVector_t  h(grad_grad_logdet[i],OrbitalSetSize);
      GGGVector_t  gh(grad_grad_grad_logdet[i],OrbitalSetSize);
      SplineAdoptor::evaluate_vghgh(P,iat,v,g,h,gh);
    }
  }
  virtual void evaluateGradSource (const ParticleSet &P
                                     , int first, int last, const ParticleSet &source
                                     , int iat_src, GradMatrix_t &gradphi)
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last,
                                     const ParticleSet &source, int iat_src,
                                     GradMatrix_t &grad_phi,
                                     HessMatrix_t &grad_grad_phi,
                                     GradMatrix_t &grad_lapl_phi)
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }
};

}
#endif
