//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ProductOrbital.h
 * @brief Declaration of \f$\Psi^c=\prod_i \psi_i({\bf R})\f$
 */
#ifndef QMCPLUSPLUS_GENERIC_PRODUCT_WITHCONSTRAINTS_H
#define QMCPLUSPLUS_GENERIC_PRODUCT_WITHCONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"

namespace qmcplusplus {

  /** A composite Orbital
   */
  struct ProductOrbital: public OrbitalBase 
  {

    ///A list of OrbitalBase* 
    vector<OrbitalBase*> Psi;
    /** Contraints on Psi
     *
     * Constraints reset optimizable variables that may be shared by Psi. 
     */
    OrbitalConstraintsBase* Constraints;

    ProductOrbital(OrbitalConstraintsBase* control):
      Constraints(control) {
        Optimizable=true;
        OrbitalName="ProductOrbital";
      }

    ~ProductOrbital();

    void setContraints(OrbitalConstraintsBase* control) {
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
    void reportStatus(ostream& os);

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
