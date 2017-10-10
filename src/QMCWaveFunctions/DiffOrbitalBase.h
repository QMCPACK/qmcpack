//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DIFFERENTIAL_ORBITAL_BASE_H
#define QMCPLUSPLUS_DIFFERENTIAL_ORBITAL_BASE_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

/**@file DiffOrbitalBase.h
 *@brief Declaration of DiffOrbitalBase, NumericalDiffOrbital and AnalyticDiffOrbital
 */
namespace qmcplusplus
{



/**@defgroup OrbitalComponent Orbital group
 * @brief Base base class for derivatives of OrbitalBase
 *
 * Derived classes implement the differentiate function which evaluates
 * - \f$\fraction{\partial \log\Psi}{\partial \alpha}\$
 * - \f$\nabla \fraction{\partial \log\Psi}{\partial \alpha}\$
 * - \f$\nabla^2 \fraction{\partial \log\Psi}{\partial \alpha}\$
 * Each object handles one or more parameters during the optimiation.
 * The data type of refOrbital, std::vector<OrbitalBase*> is intended for the cases
 * when a common variable is used by several OrbitalBase class, e.g.,
 * TwoBodyJastrowOrbital<PadeFunctor> and OneBodyJastrowOrbital<PadeFunctor>
 * with one scaling parameter.
 */
struct DiffOrbitalBase
{
  enum {SourceIndex  = OrbitalBase::SourceIndex,
        VisitorIndex = OrbitalBase::VisitorIndex,
        WalkerIndex  = OrbitalBase::WalkerIndex
       };

  //@{typedefs inherited from OrbitalBase
  typedef OrbitalBase::RealType             RealType;
  typedef OrbitalBase::ValueType            ValueType;
  typedef OrbitalBase::PosType              PosType;
  typedef ParticleSet::ParticleGradient_t   GradVectorType;
  typedef ParticleSet::ParticleLaplacian_t  ValueVectorType;
  //@}
  /** list of reference orbitals which contribute to the derivatives
   */
  std::vector<OrbitalBase*> refOrbital;

  /// default constructor
  DiffOrbitalBase(OrbitalBase* orb=0);

  ///default destructor
  virtual ~DiffOrbitalBase() { }

  /** add an OrbitalBase*
   */
  inline void addOrbital(OrbitalBase* psi)
  {
    refOrbital.push_back(psi);
  }

  /** initialize the object
   *
   * Only those that need extra checks and resizing implement this function
   */
  virtual void initialize() {}

  ///prepare internal data for a new particle sets
  virtual void resetTargetParticleSet(ParticleSet& P)=0;

  /** evaluate derivatives at \f$\{R\}\f$
   * @param P current configuration
   * @param ke0 current kinetic energy
   */
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   std::vector<RealType>& dlogpsi,
                                   std::vector<RealType>& dhpsioverpsi)=0;

  virtual void evaluateDerivRatios(ParticleSet& VP, const opt_variables_type& optvars, Matrix<ValueType>& dratios);

  virtual void multiplyDerivsByOrbR(std::vector<RealType>& dlogpsi)
  {
    for (int i=0; i<refOrbital.size(); ++i)
    {
      RealType myrat = std::exp(refOrbital[i]->LogValue)*std::cos(refOrbital[i]->PhaseValue);
      for(int j=0; j<refOrbital[i]->myVars.size(); j++)
      {
        int loc=refOrbital[j]->myVars.where(j);
        dlogpsi[loc] *= myrat;
      }
    }
  };

  /** check out optimizable variables
   */
  virtual void checkOutVariables(const opt_variables_type& optvars)=0;

  /** reset the parameters during optimizations
   */
  virtual void resetParameters(const opt_variables_type& optvars)=0;

  /** make clone
    * @param tqp target Quantum ParticleSet
    */
  virtual DiffOrbitalBasePtr makeClone(ParticleSet& tqp) const;


};

/** a generic DiffOrbitalBase using a finite-difference method for a single optimizable parameter.
 *
 * This class handles any orbital whose analytic derivatives are not implemented nor easily available.
 */
struct NumericalDiffOrbital: public DiffOrbitalBase
{

  NumericalDiffOrbital(OrbitalBase* orb=0): DiffOrbitalBase(orb) {}

  void resetTargetParticleSet(ParticleSet& P);
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi);
  void checkOutVariables(const opt_variables_type& optvars);
  void resetParameters(const opt_variables_type& optvars);

  ///\f$\nabla \partial_{\alpha} log\Psi\f$
  GradVectorType gradLogPsi, dg_p, dg_m;
  ///\f$\nabla^2 \partial_{\alpha} log\Psi\f$
  ValueVectorType lapLogPsi, dl_p, dl_m;
};

/** a generic DiffOrbitalBase using analytic derivates for a single optimizable parameter
 */
struct AnalyticDiffOrbital: public DiffOrbitalBase
{

  AnalyticDiffOrbital(OrbitalBase* orb=0): DiffOrbitalBase(orb)  { }

  void resetTargetParticleSet(ParticleSet& P);
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi);
  void checkOutVariables(const opt_variables_type& optvars);
  void resetParameters(const opt_variables_type& optvars);

  ///\f$\nabla \partial_{\alpha} log\Psi\f$
  GradVectorType gradLogPsi;
  ///\f$\nabla^2 \partial_{\alpha} log\Psi\f$
  ValueVectorType lapLogPsi;
  ///get the index in the variable list
  int MyIndex;
};
}
#endif

