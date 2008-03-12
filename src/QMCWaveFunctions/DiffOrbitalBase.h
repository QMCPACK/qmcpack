//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DIFFERENTIAL_ORBITAL_BASE_H
#define QMCPLUSPLUS_DIFFERENTIAL_ORBITAL_BASE_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"

/**@file DiffOrbitalBase.h
 *@brief Declaration of DiffOrbitalBase, NumericalDiffOrbital and AnalyticDiffOrbital
 */
namespace qmcplusplus {


  /**@defgroup OrbitalComponent Orbital group
   * @brief Base base class for derivatives of OrbitalBase
   *
   * Derived classes implement the differentiate function which evaluates
   * - \f$\fraction{\partial \log\Psi}{\partial \alpha}\$
   * - \f$\nabla \fraction{\partial \log\Psi}{\partial \alpha}\$
   * - \f$\nabla^2 \fraction{\partial \log\Psi}{\partial \alpha}\$
   */
  struct DiffOrbitalBase
  {
    //@{typedefs inherited from OrbitalBase
    typedef OrbitalBase::RealType           RealType;
    typedef OrbitalBase::ValueType          ValueType;
    typedef OrbitalBase::GradVectorType     GradVectorType;
    typedef OrbitalBase::ValueVectorType   ValueVectorType;
    typedef OrbitalBase::OptimizableSetType OptimizableSetType;
    //@}
    //
    ///\f$\frac{\partial_{\alpha}\Psi}{\Psi} = \partial_{\alpha}\log\Psi\f$
    RealType dLogPsi;
    ///\f$\frac{H\partial_{\alpha}\Psi}{\Psi} = \partial_{\alpha}\nabla^2 \log\Psi\f$
    RealType dHPsi;
    ///\f$\nabla \partial_{\alpha} log\Psi\f$
    GradVectorType gradLogPsi;
    ///\f$\nabla^2 \partial_{\alpha} log\Psi\f$
    ValueVectorType lapLogPsi;
    /** target parameter
     *
     * OptimizableSetType is derived from a map and data_type is pair<string,T>
     */
    pair<string,RealType> targetParam;
    //OptimizableSetType::data_type targetParam;

    /// default constructor
    inline DiffOrbitalBase(){ }

    ///default destructor
    virtual ~DiffOrbitalBase(){ }

    /** resize internal container
     * @param nptcls the number of particles
     */
    void resize(int nptcls);

    /**  set the name of the optimizable variable 
     * @param a name of the variable
     * @param v initial value of the variable
     */
    void setParameter(const string& a, RealType v);

    /** true, if optVariables has this object as active variables
     */
    bool isOptimizable(OptimizableSetType& optVariables);

    /** evaluate derivatives at \f$\{R\}\f$
     * @param P current configuration
     * @param ke0 current kinetic energy
     */
    void evaluateDerivatives(ParticleSet& P, RealType ke0);

    /** reset the parameters during optimizations
     */
    virtual void resetParameters(OptimizableSetType& optVariables)=0;

    /** reset properties, e.g., distance tables, for a new target ParticleSet
     * @param P ParticleSet
     */
    virtual void resetTargetParticleSet(ParticleSet& P)=0;

    /** evaluate the value of the orbital
     * @param P active ParticleSet
     * @return \f$\partial_{\alpha}\log\Psi\f$
     */
    virtual RealType differentiate(ParticleSet& P) = 0;

  };

  /** a generic DiffOrbitalBase using a finite-difference method
   *
   * This class handles any orbital whose analytic derivatives are not implemented nor easily available.
   */
  struct NumericalDiffOrbital: public DiffOrbitalBase
  {

    NumericalDiffOrbital(OrbitalBase* orb=0): refOrbital(orb) {}

    void resetParameters(OptimizableSetType& optVariables);
    void resetTargetParticleSet(ParticleSet& P);
    RealType differentiate(ParticleSet& P);

    ///pointer to the reference orbital
    OrbitalBase* refOrbital;
  };

  /** a generic DiffOrbitalBase using analytic derivates
   */
  struct AnalyticDiffOrbital: public DiffOrbitalBase
  {

    AnalyticDiffOrbital(OrbitalBase* orb=0): refOrbital(orb) {}

    void resetParameters(OptimizableSetType& optVariables);
    void resetTargetParticleSet(ParticleSet& P);
    RealType differentiate(ParticleSet& P);

    ///pointer to the reference orbital
    OrbitalBase* refOrbital;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/

