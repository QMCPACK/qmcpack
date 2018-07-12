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
 * BsplineSet<SplineAdoptor> is a SPOSetBase class to work with determinant classes
 */
#ifndef QMCPLUSPLUS_EINSPLINE_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_ADOPTOR_H

#include <Lattice/CrystalLattice.h>
#include <spline/einspline_engine.hpp>
#include <spline/einspline_util.hpp>
#include <simd/allocator.hpp>

namespace qmcplusplus
{

/** BsplineSet<SplineAdoptor>, a SPOSetBase
 * @tparam SplineAdoptor implements evaluation functions that matched the storage requirements.
 *
 * Equivalent to EinsplineSetExtended<Storage>
 * Storage is now handled by SplineAdoptor class that is specialized for precision, storage etc.
 * @todo Make SplineAdoptor be a member not the base class. This is needed
 * to make MultiBsplineSet (TBD) which has multiple SplineAdoptors for distributed
 * cases.
 */
template<typename SplineAdoptor>
struct BsplineSet: public SPOSetBase, public SplineAdoptor
{
  typedef typename SplineAdoptor::SplineType SplineType;
  typedef typename SplineAdoptor::PointType  PointType;
  typedef typename SplineAdoptor::DataType  DataType;

  ///** default constructor */
  //BsplineSet() { }

  SPOSetBase* makeClone() const
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

  inline ValueType RATIO(const ParticleSet& P, int iat, const ValueType* restrict arow)
  {
    //this is just an example how to resuse t_logpsi
    int ip=omp_get_thread_num()*2;
    // YYYY: need to fix
    //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
    return ValueType();
  }

  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    SplineAdoptor::evaluate_v(P,iat,psi);
  }

  inline void evaluateValues(VirtualParticleSet& VP, ValueMatrix_t& psiM)
  {
    SplineAdoptor::evaluateValues(VP, psiM);
  }

  inline size_t estimateMemory(const int nP)
  {
    return SplineAdoptor::estimateMemory(nP);
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    SplineAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);

#if 0
    //debug GL combo
    CONSTEXPR double eps=std::numeric_limits<float>::epsilon();
    ValueVector_t psi_copy(psi);
    GLVector_t gl(psi.size());
    SplineAdoptor::evaluate_vgl_combo(P,iat,psi_copy,gl);
    auto gradX=gl.data(0);
    auto gradY=gl.data(1);
    auto gradZ=gl.data(2);
    auto lap=gl.data(3);
    double v_err=0, g_err=0, l_err=0;
    for(size_t i=0; i<psi.size(); ++i)
    {
      v_err+=std::abs(psi[i]-psi_copy[i]);
      double dx=std::abs(dpsi[i][0]-gradX[i]);
      double dy=std::abs(dpsi[i][1]-gradY[i]);
      double dz=std::abs(dpsi[i][2]-gradZ[i]);
      g_err+=std::sqrt(dx*dx+dy*dy+dz*dz);
      l_err+=std::abs(d2psi[i]-lap[i]);
    }
    if(v_err>eps || g_err > eps || l_err>eps)
      std::cout << "ERROR " << v_err << " " << g_err << " " << l_err << std::endl;
#endif
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

  /** einspline does not need any other state data */
  void evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl, bool newp)
  {
    SplineAdoptor::evaluate_vgl_combo(P,iat,vgl);
  }

};

}
#endif
