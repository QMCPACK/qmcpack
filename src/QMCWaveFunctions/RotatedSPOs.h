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
  RotatedSPOs(std::unique_ptr<SPOSet>&& spos);
  //destructor
  ~RotatedSPOs() override;

  // Vector of rotation matrix indices
  using RotationIndices = std::vector<std::pair<int, int>>;

  // Active orbital rotation parameter indices
  RotationIndices m_act_rot_inds;
  RotationIndices m_full_rot_inds;

  // Construct a list of the matrix indices for non-zero rotation parameters.
  // (The structure for a sparse representation of the matrix)
  // Only core->active rotations are created.
  static void createRotationIndices(int nel, int nmo, RotationIndices& rot_indices);

  // Construct a list of the matrix indices for all possible rotation parameters.
  // They are orderd such that the first part of the list overlaps with the rotations
  // created by createRotationindices
  static void createRotationIndicesFull(int nel, int nmo, RotationIndices& rot_indices);

  // Fill in antisymmetric matrix from the list of rotation parameter indices
  // and a list of parameter values.
  // This function assumes rot_mat is properly sized upon input and is set to zero.
  static void constructAntiSymmetricMatrix(const RotationIndices& rot_indices,
                                           const std::vector<ValueType>& param,
                                           ValueMatrix& rot_mat);

  // Extract the list of rotation parameters from the entries in an antisymmetric matrix
  // This function expects rot_indices and param are the same length.
  static void extractParamsFromAntiSymmetricMatrix(const RotationIndices& rot_indices,
                                                   const ValueMatrix& rot_mat,
                                                   std::vector<ValueType>& param);

  // Apply rotation in delta_param to rotation in old_param.
  // Apply that rotation to MO coefficients and return the new rotation parameters in new_param
  // The size of delta_param is expected to be the size of the active rotations
  // The size of old_param and new_param are expected to the size of all rotations
  void apply_delta_rotation(const std::vector<RealType>& delta_param,
                            const std::vector<RealType>& old_param,
                            std::vector<RealType>& new_param);


  //function to perform orbital rotations
  void apply_rotation(const std::vector<RealType>& param, bool use_stored_copy);

  void apply_full_rotation(const std::vector<RealType>& full_param, bool use_stored_copy);

  // Apply the list of rotations in the history.  Used for initializing from a file.
  void apply_rotation_history();

  // Compute matrix exponential of an antisymmetric matrix (result is rotation matrix)
  static void exponentiate_antisym_matrix(ValueMatrix& mat);

  // Compute matrix log of rotation matrix to produce antisymmetric matrix
  static void log_antisym_matrix(ValueMatrix& mat);

  //A particular SPOSet used for Orbitals
  std::unique_ptr<SPOSet> Phi;

  /// Set the rotation parameters (usually from input file)
  void setRotationParameters(const std::vector<RealType>& param_list);

  /// the number of electrons of the majority spin
  size_t nel_major_;

  std::unique_ptr<SPOSet> makeClone() const override;

  // myG_temp (myL_temp) is the Gradient (Laplacian) value of of the Determinant part of the wfn
  // myG_J is the Gradient of the all other parts of the wavefunction (typically just the Jastrow).
  //       It represents \frac{\nabla\psi_{J}}{\psi_{J}}
  // myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}} . The Laplacian portion
  // IMPORTANT NOTE:  The value of P.L holds \nabla^2 ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J will hold
  ParticleSet::ParticleGradient myG_temp, myG_J;
  ParticleSet::ParticleLaplacian myL_temp, myL_J;

  ValueMatrix Bbar;
  ValueMatrix psiM_inv;
  ValueMatrix psiM_all;
  GradMatrix dpsiM_all;
  ValueMatrix d2psiM_all;


  // Single Slater creation
  void buildOptVariables(size_t nel) override;

  // For the MSD case rotations must be created in MultiSlaterDetTableMethod class
  void buildOptVariables(const RotationIndices& rotations) override;


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
                           const ValueVector& detValues_up,
                           const ValueVector& detValues_dn,
                           const GradMatrix& grads_up,
                           const GradMatrix& grads_dn,
                           const ValueMatrix& lapls_up,
                           const ValueMatrix& lapls_dn,
                           const ValueMatrix& M_up,
                           const ValueMatrix& M_dn,
                           const ValueMatrix& Minv_up,
                           const ValueMatrix& Minv_dn,
                           const GradMatrix& B_grad,
                           const ValueMatrix& B_lapl,
                           const std::vector<int>& detData_up,
                           const size_t N1,
                           const size_t N2,
                           const size_t NP1,
                           const size_t NP2,
                           const std::vector<std::vector<int>>& lookup_tbl) override;

  void evaluateDerivativesWF(ParticleSet& P,
                             const opt_variables_type& optvars,
                             std::vector<ValueType>& dlogpsi,
                             const QTFull::ValueType& psiCurrent,
                             const std::vector<ValueType>& Coeff,
                             const std::vector<size_t>& C2node_up,
                             const std::vector<size_t>& C2node_dn,
                             const ValueVector& detValues_up,
                             const ValueVector& detValues_dn,
                             const ValueMatrix& M_up,
                             const ValueMatrix& M_dn,
                             const ValueMatrix& Minv_up,
                             const ValueMatrix& Minv_dn,
                             const std::vector<int>& detData_up,
                             const std::vector<std::vector<int>>& lookup_tbl) override;

  //helper function to evaluatederivative; evaluate orbital rotation parameter derivative using table method
  void table_method_eval(std::vector<ValueType>& dlogpsi,
                         std::vector<ValueType>& dhpsioverpsi,
                         const ParticleSet::ParticleLaplacian& myL_J,
                         const ParticleSet::ParticleGradient& myG_J,
                         const size_t nel,
                         const size_t nmo,
                         const ValueType& psiCurrent,
                         const std::vector<RealType>& Coeff,
                         const std::vector<size_t>& C2node_up,
                         const std::vector<size_t>& C2node_dn,
                         const ValueVector& detValues_up,
                         const ValueVector& detValues_dn,
                         const GradMatrix& grads_up,
                         const GradMatrix& grads_dn,
                         const ValueMatrix& lapls_up,
                         const ValueMatrix& lapls_dn,
                         const ValueMatrix& M_up,
                         const ValueMatrix& M_dn,
                         const ValueMatrix& Minv_up,
                         const ValueMatrix& Minv_dn,
                         const GradMatrix& B_grad,
                         const ValueMatrix& B_lapl,
                         const std::vector<int>& detData_up,
                         const size_t N1,
                         const size_t N2,
                         const size_t NP1,
                         const size_t NP2,
                         const std::vector<std::vector<int>>& lookup_tbl);

  void table_method_evalWF(std::vector<ValueType>& dlogpsi,
                           const size_t nel,
                           const size_t nmo,
                           const ValueType& psiCurrent,
                           const std::vector<RealType>& Coeff,
                           const std::vector<size_t>& C2node_up,
                           const std::vector<size_t>& C2node_dn,
                           const ValueVector& detValues_up,
                           const ValueVector& detValues_dn,
                           const ValueMatrix& M_up,
                           const ValueMatrix& M_dn,
                           const ValueMatrix& Minv_up,
                           const ValueMatrix& Minv_dn,
                           const std::vector<int>& detData_up,
                           const std::vector<std::vector<int>>& lookup_tbl);

  void checkInVariables(opt_variables_type& active) override
  {
    if (Optimizable)
    {
      if (myVars.size())
        active.insertFrom(myVars);
    }
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    if (Optimizable)
    {
      myVars.getIndex(active);
    }
  }

  ///reset
  void resetParameters(const opt_variables_type& active) override;

  void saveExtraParameters(hdf_archive& hout, int id) override;
  void readExtraParameters(hdf_archive& hin, int id) override;


  //*********************************************************************************
  //the following functions simply call Phi's corresponding functions
  void setOrbitalSetSize(int norbs) override { Phi->setOrbitalSetSize(norbs); }

  void checkObject() const override { Phi->checkObject(); }

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override
  {
    assert(psi.size() <= OrbitalSetSize);
    Phi->evaluateValue(P, iat, psi);
  }


  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override
  {
    assert(psi.size() <= OrbitalSetSize);
    Phi->evaluateVGL(P, iat, psi, dpsi, d2psi);
  }

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<ValueType>& ratios) override
  {
    Phi->evaluateDetRatios(VP, psi, psiinv, ratios);
  }

  void evaluateVGH(const ParticleSet& P,
                   int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override
  {
    assert(psi.size() <= OrbitalSetSize);
    Phi->evaluateVGH(P, iat, psi, dpsi, grad_grad_psi);
  }


  void evaluateVGHGH(const ParticleSet& P,
                     int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override
  {
    Phi->evaluateVGHGH(P, iat, psi, dpsi, grad_grad_psi, grad_grad_grad_psi);
  }


  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override
  {
    Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, d2logdet);
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) override
  {
    Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, grad_grad_logdet);
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) override
  {
    Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, grad_grad_logdet, grad_grad_grad_logdet);
  }

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi) override
  {
    Phi->evaluateGradSource(P, first, last, source, iat_src, grad_phi);
  }

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi,
                          HessMatrix& grad_grad_phi,
                          GradMatrix& grad_lapl_phi) override
  {
    Phi->evaluateGradSource(P, first, last, source, iat_src, grad_phi, grad_grad_phi, grad_lapl_phi);
  }

  //  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix& grad_grad_grad_logdet)
  //  {Phi->evaluateThridDeriv(P, first, last, grad_grad_grad_logdet); }

private:
  /// true if SPO parameters (orbital rotation parameters) have been supplied by input
  bool params_supplied;
  /// list of supplied orbital rotation parameters
  std::vector<RealType> params;

  // List of all rotation parameters
  opt_variables_type myVarsFull;
};

} //namespace qmcplusplus

#endif
