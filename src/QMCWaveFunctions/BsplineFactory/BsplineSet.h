//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
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
#ifndef QMCPLUSPLUS_BSPLINE_SET_H
#define QMCPLUSPLUS_BSPLINE_SET_H

#include <Lattice/CrystalLattice.h>
#include <spline/einspline_engine.hpp>
#include <spline/einspline_util.hpp>
#include <simd/allocator.hpp>
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Batching.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "WhatAmI.h"

#ifdef QMC_CUDA
#include "QMCWaveFunctions/SPOSetBatched.h"
#endif
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
 * SA SplineAdoptor type
 * PST ParticleSet type
 */
  
// Unspecialized version of BsplineSet,

template<typename SA, Batching batching = Batching::SINGLE>
struct BsplineSet; //: public SA {};
// {
//   ///** default constructor */
//   //BsplineSet() { }

//   virtual void resetParameters(const opt_variables_type& active)
//   { }
//   virtual void resetTargetParticleSet(ParticleSet& e)
//   { }


// };

template<typename SA>
struct BsplineSet<SA, Batching::SINGLE>: public SA,
					 public SPOSet<Batching::SINGLE>
{
  typedef typename SA::SplineType SplineType;
  typedef typename SA::PointType  PointType;
  typedef typename SA::DataType  DataType;

  SPOSet* makeClone() const
  {
    return new BsplineSet<SA, Batching::SINGLE>(*this);
  }
  
  Batching B_ = Batching::SINGLE;

  virtual void setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
    //SplineAdoptor::first_spo=0;
    //SplineAdoptor::last_spo=norbs;
  }
  // SPOSet has these pure virtal members. Not making another wrapper for them.
  virtual void resetParameters(const opt_variables_type& active)
  { }
  virtual void resetTargetParticleSet(ParticleSet& e)
  { }

  ///** default constructor */
  BsplineSet() : SPOSet<Batching::SINGLE>::WhatAmI("BsplineSet<SA, Batching::SINGLE>") {}

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

  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    SA::evaluate_v(P,iat,psi);
  }

  inline void evaluateValues(const VirtualParticleSet& VP, ValueMatrix_t& psiM, ValueAlignedVector_t& SPOMem)
  {
    SA::evaluateValues(VP, psiM, SPOMem);
  }

  inline size_t estimateMemory(const int nP)
  {
    return SA::estimateMemory(nP);
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    SA::evaluate_vgl(P,iat,psi,dpsi,d2psi);
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
  {
    SA::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
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
      SA::evaluate_vgl(P,iat,v,g,l);
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
      SA::evaluate_vgh(P,iat,v,g,h);
    }
  }

};

  
// template<typename SA, class batching>
// inline typename QMCTraits::ValueType BsplineSet<SA, batching>::RATIO(const batching& P, int iat, const QMCTraits::ValueType* restrict arow)
// {
//     //this is just an example how to resuse t_logpsi
//     int ip=omp_get_thread_num()*2;
//     // YYYY: need to fix
//     //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
//     return QMCTraits::ValueType();
// }

// template<typename SA>
// class BsplineSet<SA, ParticleSet>
// BsplineSet<SA, ParticleSet>:: QMCTraits::ValueType RATIO(const ParticleSet& P, int iat, const QMCTraits::ValueType* restrict arow)
//   {
//     //this is just an example how to resuse t_logpsi
//     int ip=omp_get_thread_num()*2;
//     // YYYY: need to fix
//     //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
//     return QMCTraits::ValueType();
//   }


  
// template<typename SA, class batching>
// inline typename QMCTraits::ValueType BsplineSet<SA, batching>::RATIO(const batching& P, int iat, const QMCTraits::ValueType* restrict arow)
// {
//     //this is just an example how to resuse t_logpsi
//     int ip=omp_get_thread_num()*2;
//     // YYYY: need to fix
//     //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
//     return QMCTraits::ValueType();
// }

// template<typename SA>
// class BsplineSet<SA, ParticleSet>
// BsplineSet<SA, ParticleSet>:: QMCTraits::ValueType RATIO(const ParticleSet& P, int iat, const QMCTraits::ValueType* restrict arow)
//   {
//     //this is just an example how to resuse t_logpsi
//     int ip=omp_get_thread_num()*2;
//     // YYYY: need to fix
//     //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
//     return QMCTraits::ValueType();
//   }  
}
#endif
