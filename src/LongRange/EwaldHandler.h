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
class EwaldHandler: public LRHandlerBase
{

public:
  ///type of supercell
  int SuperCellEnum;
  /// Related to the Gaussian width: \f$ v_l = v(r)erf(\sigma r)\f$
  mRealType Sigma;
  ///Volume of the supercell
  mRealType Volume;
  ///Area of the supercell: always z is the slab direction
  mRealType Area;
  /** Define prefactors for the mixed boundary conditions
   *
   * For quasi-2D (see Appendix A of JPC)
   * PreFactors[0] = \f$ \frac{2\pi}{A}\f$
   * PreFactors[1] = \f$ \frac{2\pi}{A}\frac{1}{\sigma\pi}\f$
   * PreFactors[2] = \f$ \frac{2\pi}{A}\frac{1}{\sigma\pi}\f$
   * PreFactors[3] = \f$ 2\frac{\sqrt{\pi}}{A*\sigma}-\frac{2\pi}{k*A} erfc(\frac{k}{2\sigma}\f$
   */
  TinyVector<mRealType,4> PreFactors;
  ///store |k|
  std::vector<mRealType> kMag;
  /// Constructor
  EwaldHandler(ParticleSet& ref, mRealType kc_in=-1.0)
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
  EwaldHandler(const EwaldHandler& aLR, ParticleSet& ref);

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new EwaldHandler(*this,ref);
  }

  void initBreakup(ParticleSet& ref);

  void Breakup(ParticleSet& ref, mRealType rs_in)
  {
    initBreakup(ref);
  }

  void resetTargetParticleSet(ParticleSet& ref) { }

  inline mRealType evaluate(mRealType r, mRealType rinv)
  {
    return erfc(r*Sigma)*rinv;
  }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline mRealType evaluateLR(mRealType r)
  {
    return -erf(r*Sigma)/r;
  }

  inline mRealType evaluateSR_k0()
  {
    return 0.0;
  }

  inline mRealType evaluateLR_r0()
  {
    return 2.0*Sigma/std::sqrt(M_PI)+PreFactors[3];
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline mRealType srDf(mRealType r, mRealType rinv)
  {
    return 0.0;
  }

  void fillFk(KContainer& KList);

  /** evaluate k-dependent
   */
  mRealType evaluate_slab(mRealType z, const std::vector<int>& kshell
                         , const pComplexType* restrict rk1, const pComplexType* restrict rk2);

  /** evaluate k=0 term at z
   * @param z distance in the slab direction
   * @param zp z*Sigma
   * @return \f$X*z*erf(z*Sgima) + Y*exp^{-z^2*Simga^2}]\f$
   *
   * Here \f$ X=\frac{2\pi}{A}\f$ and \f$ Y=\frac{2\sqrt{\pi}}{A*Simga}\f$
   */
  inline mRealType SlabFunc0(mRealType z, mRealType zp)
  {
    return PreFactors[0]*z*erf(zp)+PreFactors[1]*std::exp(-zp*zp);
  }

  /** evaluate k!=0 term at z
   * @param ks index of this kshell
   * @param z distance in the slab direction
   * @param zp z*Sigma
   */
  inline mRealType SlabFuncK(int ks, mRealType z, mRealType zp)
  {
    mRealType expkz=std::exp(kMag[ks]*z);
    mRealType kz0=PreFactors[2]*kMag[ks];//could save this
    return erfc(kz0-zp)/expkz+erfc(kz0+zp)*expkz;
  }
};
}
#endif
