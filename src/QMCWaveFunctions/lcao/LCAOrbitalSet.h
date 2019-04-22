//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H
#define QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H

#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"

#include <Numerics/MatrixOperators.h>
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{
/** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   * SoA verson of LCOrtbitalSet
   * Localized basis set is always real 
   */
struct LCAOrbitalSet : public SPOSet
{
public:
  typedef SoaBasisSetBase<ValueType> basis_type;
  typedef basis_type::vgl_type vgl_type;
  typedef basis_type::vgh_type vgh_type;
  typedef basis_type::vghgh_type vghgh_type;
  ///pointer to the basis set
  basis_type* myBasisSet;
  ///number of Single-particle orbitals
  IndexType BasisSetSize;
  /** pointer to matrix containing the coefficients
     *
     * makeClone makes a shallow copy
     */
  ValueMatrix_t* C;
  /// Scratch space for the initial coefficents before the rotation is applied
  ValueMatrix_t m_init_B;
  /// true if SPO parameters (orbital rotation parameters) have been supplied by input
  bool params_supplied;
  /// list of supplied orbital rotation parameters
  std::vector<RealType> params;

  ///true if C is an identity matrix
  bool Identity;
  ///if true, do not clean up
  bool IsCloned;
  ///Temp(BasisSetSize) : Row index=V,Gx,Gy,Gz,L
  vgl_type Temp;
  ///Tempv(OrbitalSetSize) Tempv=C*Temp
  vgl_type Tempv;

  //These are temporary VectorSoAContainers to hold value, gradient, and hessian for
  //all basis or SPO functions evaluated at a given point.
  //Nbasis x [1(value)+3(gradient)+6(hessian)]
  vgh_type Temph;
  //Norbitals x [1(value)+3(gradient)+6(hessian)]
  vgh_type Temphv;

  //These are temporary VectorSoAContainers to hold value, gradient, hessian, and
  // gradient hessian for all basis or SPO functions evaluated at a given point.
  //Nbasis x [1(value)+3(gradient)+6(hessian)+10(grad_hessian)]
  vghgh_type Tempgh;
  //Nbasis x [1(value)+3(gradient)+6(hessian)+10(grad_hessian)]
  vghgh_type Tempghv;

  //vector that contains active orbital rotation parameter indices
  std::vector<std::pair<int, int>> m_act_rot_inds;
  /** constructor
     * @param bs pointer to the BasisSet
     */
  LCAOrbitalSet(basis_type* bs = nullptr);

  LCAOrbitalSet(const LCAOrbitalSet& in) = default;

  virtual ~LCAOrbitalSet();

  SPOSet* makeClone() const;

  /// create optimizable orbital rotation parameters
  void buildOptVariables(const std::vector<std::pair<int, int>>& rotations);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi,
                           const ValueType& psiCurrent,
                           const std::vector<ValueType>& Coeff,
                           const std::vector<size_t>& C2node_up,
                           const std::vector<size_t>& C2node_dn,
                           const ValueVector_t& detValues_up,
                           const ValueVector_t& detValues_dn,
                           const GradMatrix_t& grads_up,
                           const GradMatrix_t& grads_dn,
                           const ValueMatrix_t& lapls_up,
                           const ValueMatrix_t& lapls_dn,
                           const ValueMatrix_t& M_up,
                           const ValueMatrix_t& M_dn,
                           const ValueMatrix_t& Minv_up,
                           const ValueMatrix_t& Minv_dn,
                           const GradMatrix_t& B_grad,
                           const ValueMatrix_t& B_lapl,
                           const std::vector<int>& detData_up,
                           const size_t N1,
                           const size_t N2,
                           const size_t NP1,
                           const size_t NP2,
                           const std::vector<std::vector<int>>& lookup_tbl);


  void checkInVariables(opt_variables_type& active)
  {
    if (Optimizable && !IsCloned)
    {
      if (myVars.size())
        active.insertFrom(myVars);
      else
        Optimizable = false;
    }
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    if (Optimizable && !IsCloned)
      myVars.getIndex(active);
  }


  ///reset
  void resetParameters(const opt_variables_type& active)
  {
#if !defined(QMC_COMPLEX)
    if (Optimizable)
    {
      std::vector<RealType> param(m_act_rot_inds.size());
      for (int i = 0; i < m_act_rot_inds.size(); i++)
      {
        int loc  = myVars.where(i);
        param[i] = myVars[i] = active[loc];
      }
      apply_rotation(param);
    }
#endif
  }

  ///reset the target particleset
  void resetTargetParticleSet(ParticleSet& P)
  {
    //myBasisSet->resetTargetParticleSet(P);
  }

  /** set the OrbitalSetSize
    */
  virtual void setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
    Tempv.resize(OrbitalSetSize);
    Temphv.resize(OrbitalSetSize);
    Tempghv.resize(OrbitalSetSize);
  }

  /** set the basis set
    */
  void setBasisSet(basis_type* bs);

  /** return the size of the basis set
    */
  int getBasisSetSize() const { return (myBasisSet == nullptr) ? 0 : myBasisSet->getBasisSetSize(); }

  bool setIdentity(bool useIdentity);

  void checkObject() const
  {
    if (!(OrbitalSetSize == C->rows() && BasisSetSize == C->cols()))
      APP_ABORT("   LCAOrbitalSet::checkObject Linear coeffient for LCAOrbitalSet is not consistent with the input.");
  }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector_t& psi,
                         const ValueVector_t& psiinv,
                         std::vector<ValueType>& ratios);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi);

  void evaluate(const ParticleSet& P, 
               int iat, 
               ValueVector_t& psi, 
               GradVector_t& dpsi, 
               HessVector_t& grad_grad_psi,
               GGGVector_t& grad_grad_grad_psi);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet,
                            GGGMatrix_t& grad_grad_grad_logdet);

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix_t& gradphi);

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix_t& grad_phi,
                          HessMatrix_t& grad_grad_phi,
                          GradMatrix_t& grad_lapl_phi);

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet);

private:
  //helper functions to handl Identity
  void evaluate_vgl_impl(const vgl_type& temp, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) const;

  void evaluate_vgl_impl(const vgl_type& temp,
                         int i,
                         ValueMatrix_t& logdet,
                         GradMatrix_t& dlogdet,
                         ValueMatrix_t& d2logdet) const;
  //These two functions unpack the data in vgh_type temp object into wavefunction friendly data structures.
  //This unpacks temp into vectors psi, dpsi, and d2psi.
  void evaluate_vgh_impl(const vgh_type& temp, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& d2psi) const;

  //Unpacks temp into the ith row (or electron index) of logdet, dlogdet, dhlogdet.
  void evaluate_vgh_impl(const vgh_type& temp,
                         int i,
                         ValueMatrix_t& logdet,
                         GradMatrix_t& dlogdet,
                         HessMatrix_t& dhlogdet) const;
  //Unpacks data in vghgh_type temp object into wavefunction friendly data structures for value, gradient, hessian
  //and gradient hessian.  
  void evaluate_vghgh_impl(const vghgh_type& temp, 
                           ValueVector_t& psi, 
                           GradVector_t& dpsi, 
                           HessVector_t& d2psi,
                           GGGVector_t& dghpsi) const;

  void evaluate_vghgh_impl(const vghgh_type& temp,
                         int i,
                         ValueMatrix_t& logdet,
                         GradMatrix_t& dlogdet,
                         HessMatrix_t& dhlogdet,
                         GGGMatrix_t& dghlogdet) const;
  
  void evaluate_ionderiv_vgl_impl(const vghgh_type& temp,
                         int i,
                         GradMatrix_t& dlogdet,
                         HessMatrix_t& dglogdet,
                         GradMatrix_t& dllogdet) const;
  

#if !defined(QMC_COMPLEX)
  //function to perform orbital rotations
  void apply_rotation(const std::vector<RealType>& param);

  //helper function to apply_rotation
  void exponentiate_antisym_matrix(ValueMatrix_t& mat);


  //helper function to evaluatederivative; evaluate orbital rotation parameter derivative using table method
  void table_method_eval(std::vector<RealType>& dlogpsi,
                         std::vector<RealType>& dhpsioverpsi,
                         const ParticleSet::ParticleLaplacian_t& myL_J,
                         const ParticleSet::ParticleGradient_t& myG_J,
                         const size_t nel,
                         const size_t nmo,
                         const ValueType& psiCurrent,
                         const std::vector<RealType>& Coeff,
                         const std::vector<size_t>& C2node_up,
                         const std::vector<size_t>& C2node_dn,
                         const ValueVector_t& detValues_up,
                         const ValueVector_t& detValues_dn,
                         const GradMatrix_t& grads_up,
                         const GradMatrix_t& grads_dn,
                         const ValueMatrix_t& lapls_up,
                         const ValueMatrix_t& lapls_dn,
                         const ValueMatrix_t& M_up,
                         const ValueMatrix_t& M_dn,
                         const ValueMatrix_t& Minv_up,
                         const ValueMatrix_t& Minv_dn,
                         const GradMatrix_t& B_grad,
                         const ValueMatrix_t& B_lapl,
                         const std::vector<int>& detData_up,
                         const size_t N1,
                         const size_t N2,
                         const size_t NP1,
                         const size_t NP2,
                         const std::vector<std::vector<int>>& lookup_tbl);
#endif
};
} // namespace qmcplusplus
#endif
