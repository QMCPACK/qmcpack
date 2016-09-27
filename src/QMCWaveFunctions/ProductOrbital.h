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
    
    
/** @file ProductOrbital.h
 * @brief Declaration of \f$\Psi^c=\prod_i \psi_i({\bf R})\f$
 */
#ifndef QMCPLUSPLUS_GENERIC_PRODUCT_WITHCONSTRAINTS_H
#define QMCPLUSPLUS_GENERIC_PRODUCT_WITHCONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"

namespace qmcplusplus
{

/** A composite Orbital
 */
struct ProductOrbital: public OrbitalBase
{

  ///A list of OrbitalBase*
  std::vector<OrbitalBase*> Psi;
  /** Contraints on Psi
   *
   * Constraints reset optimizable variables that may be shared by Psi.
   */
  OrbitalConstraintsBase* Constraints;

  ProductOrbital(OrbitalConstraintsBase* control):
    Constraints(control)
  {
    Optimizable=true;
    OrbitalName="ProductOrbital";
  }

  ~ProductOrbital();

  void setContraints(OrbitalConstraintsBase* control)
  {
    Constraints=control;
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& o);

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& o);

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os);

  /** reset the parameters during optimizations
   */
  void resetParameters(const opt_variables_type& active);

  void resetTargetParticleSet(ParticleSet& P);

  ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL);

  ValueType ratio(ParticleSet& P, int iat);

  void acceptMove(ParticleSet& P, int iat);

  void restore(int iat);

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat);

  RealType
  registerData(ParticleSet& P, BufferType& buf);

  RealType
  updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);

  void
  copyFromBuffer(ParticleSet& P, BufferType& buf);

  RealType
  evaluateLog(ParticleSet& P,BufferType& buf);

  OrbitalBase* makeClone(ParticleSet& tqp) const;

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
