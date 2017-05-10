//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_EWALD_HANLDER_H
#define QMCPLUSPLUS_EWALD_HANLDER_H

#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{

/* LR breakup for the standard Ewald method
 *
 * Quasi-2D Ewald method : J. Phys.: Condens. Matter 16, 891 (2004)
 * http://iopscience.iop.org/0953-8984/16/6/017/
 * Note that \f$ \simga \rightarrow 1/\sigma\f$
 * It is possible to use 3D Ewald but for the bulk system, the optimal breakup method
 * is used.
 */
class TwoDEwaldHandler: public LRHandlerBase
{

public:
  /// Related to the Gaussian width: \f$ v_l = v(r)erf(\sigma r)\f$
  RealType Sigma;
  RealType Volume;
  ///store |k|
  std::vector<RealType> kMag;
  /// Constructor
  TwoDEwaldHandler(ParticleSet& ref, RealType kc_in=-1.0)
    : LRHandlerBase(kc_in)
  {
    Sigma=LR_kc=ref.Lattice.LR_kc;
  }

  /** "copy" constructor
   * @param aLR LRHandlerTemp
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  TwoDEwaldHandler(const TwoDEwaldHandler& aLR, ParticleSet& ref);

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new TwoDEwaldHandler(*this,ref);
  }

  void initBreakup(ParticleSet& ref);

  void Breakup(ParticleSet& ref, RealType rs_in)
  {
    initBreakup(ref);
  }

  void resetTargetParticleSet(ParticleSet& ref) { }

  inline RealType evaluate(RealType r, RealType rinv)
  {
    return erfc(r*Sigma)*rinv;
  }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline RealType evaluateLR(RealType r)
  {
    return -erf(r*Sigma)/r;
  }

  inline RealType evaluateSR_k0()
  {
    return 2.0*std::sqrt(M_PI)/(Sigma*Volume);
  }

  inline RealType evaluateLR_r0()
  {
    return 2.0*Sigma/std::sqrt(M_PI);
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline RealType srDf(RealType r, RealType rinv)
  {
    return 0.0;
  }

  void fillFk(KContainer& KList);

};
}
#endif
