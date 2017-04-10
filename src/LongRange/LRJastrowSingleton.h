//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRJastrowSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRJASTROWSINGLETON_H
#define QMCPLUSPLUS_LRJASTROWSINGLETON_H

#include "LongRange/LRHandlerTemp.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include <limits>

namespace qmcplusplus
{
/** JastrowFunctor
 *
 * A Func for LRHandlerTemp.  Four member functions have to be provided
 *
 * - reset(T volume) : reset the normalization factor
 * - operator() (T r, T rinv) : return a value of the original function e.g., 1.0/r
 * - Fk(T k, T rc)
 * - Xk(T k, T rc)
 *
 */

template<class T=double>
struct JastrowFunctor
{
  T Rs;
  T SqrtRs;
  T OneOverSqrtRs;
  T NormFactor;
  inline JastrowFunctor() {}

  void reset(ParticleSet& ref)
  {
    reset(ref.getTotalNum(),ref.Lattice.Volume);
    //NormFactor=4.0*M_PI/ref.Lattice.Volume;
    //T Density=ref.getTotalNum()/ref.Lattice.Volume;
    //Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //SqrtRs=std::sqrt(Rs);
    //OneOverSqrtRs = 1.0 / SqrtRs;
  }

  void reset(ParticleSet& ref, T rs)
  {
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    Rs = rs;
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 /SqrtRs;
  }

  /** reset by the number of particles and the volume
   * @param n  number of particles
   * @param vol volume
   */
  void reset(int n, T vol)
  {
    NormFactor=4.0*M_PI/vol;
    T Density=static_cast<T>(n)/vol;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
  }

  inline T operator()(T r, T rinv)
  {
    if(r< std::numeric_limits<T>::epsilon())
      return SqrtRs-0.5*r;
    else
      return Rs*rinv*(1.0-std::exp(-r*OneOverSqrtRs));
    //if (r > 1e-10) return Rs*rinv*(1.0 - std::exp(-r*OneOverSqrtRs));
    //return 1.0 / OneOverSqrtRs - 0.5 * r;
  }

  inline T df(T r, T rinv)
  {
    if(r< std::numeric_limits<T>::epsilon())
      return -0.5+r*OneOverSqrtRs/3.0;
    else
    {
      T exponential = std::exp(-r*OneOverSqrtRs);
      return -Rs*rinv*rinv*(1.0 - exponential) + exponential*rinv*SqrtRs;
    }
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T coskr = std::cos(k*rc);
    T sinkr = std::sin(k*rc);
    T oneOverK = 1.0/k;
    return -NormFactor * Rs *
           (coskr*oneOverK*oneOverK
            - std::exp(-rc*OneOverSqrtRs)*(coskr - OneOverSqrtRs * sinkr * oneOverK)/(k*k+1.0/Rs));
  }
};


struct LRJastrowSingleton
{

  typedef LRHandlerTemp<JastrowFunctor<double>,LPQHIBasis> LRHandlerType;

  static LRHandlerType* JastrowHandler;

  static LRHandlerType* getHandler(ParticleSet& ref, double kc);
};
}
#endif
