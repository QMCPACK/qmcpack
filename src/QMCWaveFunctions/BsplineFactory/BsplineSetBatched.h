//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file BsplineSetBatched.h
 *
 * BsplineSet<SplineAdoptor, Batching> is a SPOSet class to work with determinant classes
 */
#ifndef QMCPLUSPLUS_BSPLINE_SET_BATCHED_H
#define QMCPLUSPLUS_BSPLINE_SET_BATCHED_H

#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include <Lattice/CrystalLattice.h>
#include <spline/einspline_engine.hpp>
#include <spline/einspline_util.hpp>
#include <simd/allocator.hpp>
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "QMCWaveFunctions/Batching.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSetBatched.h"

namespace qmcplusplus
{

template<typename SA>
struct BsplineSet<SA, Batching::BATCHED>: public SA,
                                          public SPOSet<Batching::BATCHED>
{
  typedef typename SA::SplineType SplineType;
  typedef typename SA::PointType  PointType;
  typedef typename SA::DataType  DataType;

  Batching B_ = Batching::BATCHED;

  BsplineSet() : SPOSet<Batching::BATCHED>::WhatAmI("BsplineSet<SA, Batching::BATCHED>") {}
  
  SPOSet* makeClone() const
  {
    return new BsplineSet<SA, Batching::BATCHED>(*this);
  }

  // SPOSet has these pure virtal members. Not making another wrapper for them.
  // Should be possible to template SPOSet and put in common class.
  virtual void setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
    //SplineAdoptor::first_spo=0;
    //SplineAdoptor::last_spo=norbs;
  }
  
  virtual void resetParameters(const opt_variables_type& active)
  { }
  virtual void resetTargetParticleSet(ParticleSet& e)
  { }

  ///** default constructor */
  //BsplineSet() { }

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
    SA::set_spline(spline_r,spline_i,twist,ispline,level);
  }

  QMCTraits::ValueType RATIO(int iat, const QMCTraits::ValueType* restrict arow)
  {
  //this is just an example how to resuse t_logpsi
    int ip=omp_get_thread_num()*2;
    // YYYY: need to fix
    //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
    return QMCTraits::ValueType();
  }

  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    APP_ABORT("Single evaluations not supported by BsplineSet<BATCHED>");
    //SA::evaluate_v(P,iat,psi);
  }
  
  inline void evaluate(const std::vector<ParticleSet*>& P, int iat, ValueVector_t& psi)
  {
    SA::evaluate_v(P,iat,psi);
  }

  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::ValueVector_t& d2psi)
  {
    APP_ABORT("Single evaluations not supported by BsplineSet<BATCHED>");
  }

  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::HessVector_t& grad_grad_psi)
  {
    APP_ABORT("Single evaluations not supported by BsplineSet<BATCHED>");
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::ValueMatrix_t& d2logdet)
  {
    APP_ABORT("Single evaluations not supported by BsplineSet<BATCHED>");
  }

  
  inline void evaluateValues(const std::vector<VirtualParticleSet*>& VP, ValueMatrix_t& psiM, ValueAlignedVector_t& SPOMem)
  {
    SA::evaluateValues(VP, psiM, SPOMem);
  }

  inline size_t estimateMemory(const int nP)
  {
    return SA::estimateMemory(nP);
  }

  inline void evaluate(const std::vector<ParticleSet*>& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    SA::evaluate_vgl(P,iat,psi,dpsi,d2psi);

  }

  inline void evaluate(const std::vector<ParticleSet*>& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
  {
    SA::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
  }

  void evaluate_notranspose(const std::vector<ParticleSet*>& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      ValueVector_t v(logdet[i],OrbitalSetSize);
      GradVector_t  g(dlogdet[i],OrbitalSetSize);
      ValueVector_t l(d2logdet[i],OrbitalSetSize);
      SA::evaluate_vgl(P,iat,v,g,l);
    }
  }

  virtual void evaluate_notranspose(const std::vector<ParticleSet*>& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Single evaluations not supported by BsplineSet<BATCHED>");
  }

  /** einspline does not need any other state data */
  void evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
  {
    APP_ABORT("Single evaluations not supported by BsplineSet<BATCHED>");
  }

};

}

#endif
