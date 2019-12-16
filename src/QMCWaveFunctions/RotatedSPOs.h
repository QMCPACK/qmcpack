//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) QMCPACK developers.
////
//// File developed by: Sergio D. Pineda Flores, sergio_pinedaflores@berkeley.edu, University of California, Berkeley
////                    Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
////                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////
//// File created by: Sergio D. Pineda Flores, sergio_pinedaflores@berkeley.edu, University of California, Berkeley
////////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_ROTATION_HELPER_H
#define QMCPLUSPLUS_ROTATION_HELPER_H

#include "QMCWaveFunctions/SPOSet.h"


namespace qmcplusplus
{
class RotatedSPOs : public SPOSet
{
public:
  //constructor
  RotatedSPOs(SPOSet* spos);
  //destructor
  ~RotatedSPOs();

  //vector that contains active orbital rotation parameter indices
  std::vector<std::pair<int, int>> m_act_rot_inds;

  //function to perform orbital rotations
  void apply_rotation(const std::vector<RealType>& param, bool use_stored_copy);

  //helper function to apply_rotation
  void exponentiate_antisym_matrix(ValueMatrix_t& mat);

  //A particualr SPOSet used for Orbitals
  SPOSet* Phi;

  /// true if SPO parameters (orbital rotation parameters) have been supplied by input
  bool params_supplied;
  /// list of supplied orbital rotation parameters
  std::vector<RealType> params;

  bool IsCloned;

  /// the number of electrons of the majority spin
  size_t nel_major_;

  SPOSet* makeClone() const override;

  // myG_temp (myL_temp) is the Gradient (Laplacian) value of of the Determinant part of the wfn
  // myG_J is the Gradient of the all other parts of the wavefunction (typically just the Jastrow).
  //       It represents \frac{\nabla\psi_{J}}{\psi_{J}}
  // myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}} . The Laplacian portion
  // IMPORTANT NOTE:  The value of P.L holds \nabla^2 ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J will hold
  ParticleSet::ParticleGradient_t myG_temp, myG_J;
  ParticleSet::ParticleLaplacian_t myL_temp, myL_J;

  ValueMatrix_t Bbar;
  ValueMatrix_t psiM_inv;
  ValueMatrix_t psiM_all;
  GradMatrix_t dpsiM_all;
  ValueMatrix_t d2psiM_all;


  // Single Slater creation
  void buildOptVariables(const size_t nel) override;
  // For the MSD case rotations must be created in MultiSlaterFast class
  void buildOptVariables(const std::vector<std::pair<int, int>>& rotations) override;


  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi,
                           const int& FirstIndex,
                           const int& LastIndex) override;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi,
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
                           const std::vector<std::vector<int>>& lookup_tbl) override;

  //helper function to evaluatederivative; evaluate orbital rotation parameter derivative using table method
  void table_method_eval(std::vector<ValueType>& dlogpsi,
                         std::vector<ValueType>& dhpsioverpsi,
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

  void checkInVariables(opt_variables_type& active) override
  {
    //reset parameters to zero after coefficient matrix has been updated
    for (int k = 0; k < myVars.size(); ++k)
      myVars[k] = 0.0;

    if (Optimizable && !IsCloned)
    {
      if (myVars.size())
        active.insertFrom(myVars);
      Phi->storeParamsBeforeRotation();
    }
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    if (Optimizable && !IsCloned)
    {
      myVars.getIndex(active);
    }
  }

  ///reset
  void resetParameters(const opt_variables_type& active) override
  {
    if (Optimizable && !IsCloned)
    {
      std::vector<RealType> param(m_act_rot_inds.size());
      for (int i = 0; i < m_act_rot_inds.size(); i++)
      {
        int loc  = myVars.where(i);
        param[i] = myVars[i] = active[loc];
      }
      apply_rotation(param, true);
    }
  }
  //*********************************************************************************
  //the following functions simply call Phi's corresponding functions
  void resetTargetParticleSet(ParticleSet& P) override { Phi->resetTargetParticleSet(P); }

  void setOrbitalSetSize(int norbs) override { Phi->setOrbitalSetSize(norbs); }

  //  void setBasisSet(basis_type* bs);

  int getBasisSetSize() { return Phi->getBasisSetSize(); }

  //  bool setIdentity(bool useIdentity)
  //  {return Phi->setIdentity(useIdentity); }

  void checkObject() { Phi->checkObject(); }

  void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi) override
  {
    assert(psi.size() <= OrbitalSetSize);
    Phi->evaluateValue(P, iat, psi);
  }


  void evaluateVGL(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) override
  {
    assert(psi.size() <= OrbitalSetSize);
    Phi->evaluateVGL(P, iat, psi, dpsi, d2psi);
  }

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector_t& psi,
                         const ValueVector_t& psiinv,
                         std::vector<ValueType>& ratios) override
  {
    Phi->evaluateDetRatios(VP, psi, psiinv, ratios);
  }

  void evaluateVGH(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi) override
  {
    assert(psi.size() <= OrbitalSetSize);
    Phi->evaluateVGH(P, iat, psi, dpsi, grad_grad_psi);
  }


  void evaluateVGHGH(const ParticleSet& P,
                     int iat,
                     ValueVector_t& psi,
                     GradVector_t& dpsi,
                     HessVector_t& grad_grad_psi,
                     GGGVector_t& grad_grad_grad_psi) override
  {
    Phi->evaluateVGHGH(P, iat, psi, dpsi, grad_grad_psi, grad_grad_grad_psi);
  }


  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet) override
  {
    Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, d2logdet);
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet) override
  {
    Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, grad_grad_logdet);
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet,
                            GGGMatrix_t& grad_grad_grad_logdet) override
  {
    Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, grad_grad_logdet, grad_grad_grad_logdet);
  }

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix_t& grad_phi) override
  {
    Phi->evaluateGradSource(P, first, last, source, iat_src, grad_phi);
  }

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix_t& grad_phi,
                          HessMatrix_t& grad_grad_phi,
                          GradMatrix_t& grad_lapl_phi) override
  {
    Phi->evaluateGradSource(P, first, last, source, iat_src, grad_phi, grad_grad_phi, grad_lapl_phi);
  }

  //  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet)
  //  {Phi->evaluateThridDeriv(P, first, last, grad_grad_grad_logdet); }
};

} //namespace qmcplusplus

#endif
