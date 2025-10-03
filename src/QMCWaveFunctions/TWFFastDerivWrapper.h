//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//                      Anouar Benali, abenali.sci@gmail.com, Qubit Pharmaceuticals
//
// File created by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TWFFASTDERIVWRAPPER_H
#define QMCPLUSPLUS_TWFFASTDERIVWRAPPER_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "ResourceCollection.h"
#include "Utilities/ResourceCollection.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{
/**
 *  TWFFastDerivWrapper is a wrapper class for TrialWavefunction that provides separate and low level access to the Jastrow and 
 *  SPOSet objects.  This is so that observables can be recast in matrix form and their derivatives taken efficiently. 
 *  Currently this is hard coded for ground state slater jastrow wave functions, but generalization to
 *  arbitrary occupations and multideterminants are straightforward and will come online as the code is tested and validated. 
 *
 *  Please see : J. Chem. Phys. 144, 194105 (2016) https://doi.org/10.1063/1.4948778 for implementation details and formalism.
 */
class TWFFastDerivWrapper
{
public:
  using ValueMatrix = SPOSet::ValueMatrix;
  using GradMatrix  = SPOSet::GradMatrix;
  using HessMatrix  = SPOSet::HessMatrix;
  using IndexType   = QMCTraits::IndexType;
  using RealType    = QMCTraits::RealType;
  using ValueType   = QMCTraits::ValueType;
  using GradType    = QMCTraits::GradType;
  using ValueVector = SPOSet::ValueVector;
  using GradVector  = SPOSet::GradVector;

  // used when we get const reference to *MultiDiracDeterminants::detData
  template<typename DT>
  using OffloadVector = Vector<DT, OffloadPinnedAllocator<DT>>;

  TWFFastDerivWrapper() = default;
  /** @brief Add a particle group.
   *
   *  Here, a "group" corresponds to a subset of particles which are antisymmetric with 
   *  respect to eachother.  TWFFastDerivWrapper ensures that there is a binding between the groupid
   *  in ParticleSet and the SPOSet associated with that particle group.  This function stores
   *  the ParticleSet groupid and SPOSet in a vector for lookup and communication with QMCPACK conventions, 
   *  but is agnostic to the order of group registration or evaluation.  
   *
   *  @param[in] P.  ParticleSet
   *  @param[in] groupid.  ParticleSet groupid to which SPOSet corresponds.
   *  @param[in] spo.  Pointer to SPOSet.  
   *  @return void.
   */
  void addGroup(const ParticleSet& P, const IndexType groupid, SPOSet* spo);
  inline void addJastrow(WaveFunctionComponent* j) { jastrow_list_.push_back(j); };

  /**
   * @brief add a MultiSlaterDet
   * 
   * this is called when registering a MSD
   * the caller should also iterate over MultiDiracDeterminants in the MSD and register the associated SPOSets
   * 
   * @param P. ParticleSet 
   * @param msd. Pointer to MultiSlaterDetTableMethod
   */
  void addMultiSlaterDet(const ParticleSet& P, const WaveFunctionComponent* msd);
  inline bool hasMultiSlaterDet() const { return multislaterdet_; }
  inline const WaveFunctionComponent& getMultiSlaterDet() const
  {
    if (!multislaterdet_)
    {
      throw std::runtime_error("Error: No MultiSlaterDet registered for this TWFFastDerivWrapper");
    }
    return *multislaterdet_;
  }


  /** @brief Takes particle set groupID and returns the TWF internal index for it.  
   *
   *  ParticleSet groups can be registered in whichever order.  However, the internal indexing 
   *  of TWFFastDerivWrapper keeps track on a first-come, first serve basis.  That is, if I register 
   *  particle groups 3, 1, and 0 in that order, then the internal indexing goes like
   *  0->3, 1->1, 2->0.  Thus, this function looks up where in the internal indexing scheme
   *  ParticleSet gid is located.  This is necessary, since the matrix list data structures follow
   *  the internal TWF indexing.  
   *
   *  @param[in] gid. ParticleSet group ID to look up.  
   *  @return The corresponding internal groupID.  
   */
  IndexType getTWFGroupIndex(const IndexType gid) const;

  /** @ingroup Query functions
   * @{
   * @brief These are query functions to the internal state of various groups and SPOSet info.
   *   All group indices will refer to the internal group indices.  Convert from ParticleSet to 
   *   TWF group first!
   *
   *   Source of truth for orbital sizes will be the individual SPOSets.  Particle group sizes
   *   will be ParticleSet in conjunction with groupID maps.  
   */
  inline IndexType numGroups() const { return spos_.size(); };
  SPOSet* getSPOSet(const IndexType sid) const { return spos_[sid]; };
  inline IndexType numOrbitals(const IndexType sid) const { return spos_[sid]->size(); };
  /** @} */

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //These are convenience functions/wrappers to SPOSet calls.  Idea being that observables just need      //
  //to make calls to this object to build the auxiliary matrices required for fast derivative computation.//
  //On the other hand, it wouldn't be unreasonable to make the observables do the SPOSet calls themselves.//
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** @brief Returns the non-rectangular slater matrices M_ij=phi_j(r_i) (N_particle x N_orbs) for each species group. 
   *
   *  @param[in] P particle set.
   *  @param[in,out] mmat Output vector of slater matrices.  Each vector entry corresponds to a different particle group.
   *  @return Void
   */
  void getM(const ParticleSet& P, std::vector<ValueMatrix>& mmat) const;

  /** @brief Returns value of all orbitals (relevant to given species/group) at a particular particle coordinate.
   *
   *  @param[in] P particle set.
   *  @param[in] iel particle ID.
   *  @param[in,out] val Vector of phi_i(r_iel) for all i=0,Norbs.
   *  @return Void
   */
  IndexType getRowM(const ParticleSet& P, const IndexType iel, ValueVector& val) const;

  /** @brief Returns value, gradient, and laplacian matrices for all orbitals and all particles, species by species. 
   *
   *  @param[in] P particle set.
   *  @param[in,out] mvec Slater matrix M_ij=phi_j(r_i) for each species group.
   *  @param[in,out] gmat electron gradient of slater matrix [G_ij]_a = d/dr_a,i phi_j(r_i).  a=x,y,z  
   *  @param[in,out] lmat electron laplacian of slater matrix [L_ij] = \nabla^2_i phi_j(r_i).
   *  @return Void
   */
  void getEGradELaplM(const ParticleSet& P,
                      std::vector<ValueMatrix>& mvec,
                      std::vector<GradMatrix>& gmat,
                      std::vector<ValueMatrix>& lmat) const;

  /** @brief Returns x,y,z components of ion gradient of slater matrices.
   *
   *  @param[in] P particle set.
   *  @param[in] source ion particle set.
   *  @param[in]*** iat ion ID w.r.t. which to take derivative.
   *  @param[in,out] dmvec Slater matrix d/dR_{iat,a} M_ij=d/dR_{iat,a} phi_j(r_i) for each species group.  First index is a=x,y,z.
   *  @return Void
   */
  void getIonGradM(const ParticleSet& P,
                   const ParticleSet& source,
                   const int iat,
                   std::vector<std::vector<ValueMatrix>>& dmvec) const;

  /** @brief Returns x,y,z components of ion gradient of slater matrices and their laplacians..
   *
   *  @param[in] P particle set.
   *  @param[in] source ion particle set.
   *  @param[in] iat ion ID w.r.t. which to take derivative.
   *  @param[in,out] dmvec Slater matrix d/dR_{iat,a} M_ij=d/dR_{iat,a} phi_j(r_i) for each species group.  
   *                 First index is a=x,y,z.
   *  @param[in,out] dlmat Slater matrix d/dR_{iat,a} L_ij=d/dR_{iat,a} \nabla^2_i phi_j(r_i) for each species group.  
   *                 First index is a=x,y,z.
   *  @return Void
   */
  void getIonGradIonGradELaplM(const ParticleSet& P,
                               const ParticleSet& source,
                               int iat,
                               std::vector<std::vector<ValueMatrix>>& dmvec,
                               std::vector<std::vector<GradMatrix>>& dgmat,
                               std::vector<std::vector<ValueMatrix>>& dlmat) const;


  /** @brief Evaluates value, gradient, and laplacian of the total jastrow.  So of J in exp(J).
    *
    * @param[in] Particle set.
    * @param[in,out] Electron gradients.
    * @param[in,out] Electron laplacians.
    * @return Jastrow value
    */
  RealType evaluateJastrowVGL(const ParticleSet& P,
                              ParticleSet::ParticleGradient& G,
                              ParticleSet::ParticleLaplacian& L) const;

  /** @brief Given a proposed move of electron iel from r->r', computes exp(J')/exp(J)
    *
    * @param[in] P electron particle set.
    * @param[in] iel electron being moved.
    * @return Jastrow ratio.
    */
  RealType evaluateJastrowRatio(ParticleSet& P, const int iel) const;

  /** @brief Given a proposed move of electron iel from r->r', computes exp(J(r'))/exp(J(r))
    *   and grad_iel(J(r')).)
    *
    * @param[in] P electron particle set.
    * @param[in] iel electron being moved.
    * @param[in,out] grad grad_iel(J(r')). So iel electron gradient at new position.
    * @return Jastrow ratio.
    */
  RealType calcJastrowRatioGrad(ParticleSet& P, const int iel, GradType& grad) const;

  /** @brief Return ionic gradient of J(r).
    *
    * @param[in] P electron particle set.
    * @param[in] source ion particle set.
    * @param[in] iat Ion to take derivative w.r.t.
    * @return Gradient of J(r) w.r.t. ion iat.
    */
  GradType evaluateJastrowGradSource(ParticleSet& P, ParticleSet& source, const int iat) const;

  /** @brief Return ionic gradients of J(r), grad_iel(J(r)), and lapl_iel(J(r)).
    *
    * @param[in] P electron particle set.
    * @param[in] source ion particle set.
    * @param[in] iat Ion to take derivative w.r.t.
    * @param[in,out] grad_grad iat ion gradient of electron gradients of J(r).  
    * @param[in,out] grad_lapl iat ion gradient of electron laplacians of J(r).
    * @return Gradient of J(r) w.r.t. ion iat.
    */
  GradType evaluateJastrowGradSource(ParticleSet& P,
                                     ParticleSet& source,
                                     const int iat,
                                     TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                                     TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) const;

  /** @brief Takes sub matrices of full SPOSet quantities (evaluated on all particles and all orbitals), consistent with ground
   *   state occupations.
   *
   *  @param[in] A non-rectangular matrices of SPOSet derived quantities, evaluated on all orbitals and all particles.
   *  @param[in,out] Aslice square matrices consistent with a ground state occupation.
   *  @return Void
   */
  void getGSMatrices(const std::vector<ValueMatrix>& A, std::vector<ValueMatrix>& Aslice) const;

  /** @brief Calculates derivative of observable via Tr[M^{-1} dB - X * dM ].  Consistent with ground state occupation.
   *
   *  @param[in] Minv. inverse of slater matrices for ground state occupations. 
   *  @param[in] X.  X=M^-1 B M^-1.  B observable matrix, and is consistent with ground state occupation.
   *  @param[in] dM. Target derivative of M, and is consistent with ground state occupation.
   *  @param[in] dB. Target derivative of B, and is consistent with ground state occupation.
   *  @return Derivative of O psi/psi = Tr[M^{-1} dB - X * dM ]
   */
  ValueType computeGSDerivative(const std::vector<ValueMatrix>& Minv,
                                const std::vector<ValueMatrix>& X,
                                const std::vector<ValueMatrix>& dM,
                                const std::vector<ValueMatrix>& dB) const;

  /** @brief Calculates derivative of observable and wavefunction for several MultiDiracDeterminant objects
   *         calculated quantities are relative to refdet (first element is 0 in output arrays)
   *
   *  @param[in] Minv_Mv. M^-1.M(virt)
   *  @param[in] Minv_B.  M^-1.B
   *  @param[in] Minv_dM.  M^-1.dM
   *  @param[in] Minv_dB.  M^-1.dB
   *  @param[in] mdd_spo_ids. mapping from index in mdds into SPOSets (i.e. into first dim of vec of matrices Minv, X, dM, dB, etc.)
   *  @param[in] mdds. vector of MultiDiracDeterminants; evaluate derivative of observable for all determinants
   *  @param[out] dvals_dmu_O. d/dmu(O D[i][j]/D[i][j]) for MultiDiracDet i, excited det j (also includes GS det at j==0)
   *  @param[out] dvals_dmu. d/dmu(log(D[i][j]) for MultiDiracDet i, excited det j (also includes GS det at j==0)
   *  @return 
   */
  void computeMDDerivatives_dmu(const std::vector<ValueMatrix>& Minv_Mv,
                                const std::vector<ValueMatrix>& Minv_B,
                                const std::vector<ValueMatrix>& Minv_dM,
                                const std::vector<ValueMatrix>& Minv_dB,
                                const std::vector<IndexType>& mdd_spo_ids,
                                const std::vector<const WaveFunctionComponent*>& mdds,
                                std::vector<ValueVector>& dvals_dmu_O,
                                std::vector<ValueVector>& dvals_dmu) const;

  /** @brief Calculates observable for several MultiDiracDeterminant objects
   *         calculated quantities are relative to refdet (first element is 0 in output arrays)
   *
   *  @param[in] Minv_Mv. M^-1.M(virt)
   *  @param[in] Minv_B.  M^-1.B
   *  @param[in] mdd_spo_ids. mapping from index in mdds into SPOSets (i.e. into first dim of vec of matrices Minv, X, dM, dB, etc.)
   *  @param[in] mdds. vector of MultiDiracDeterminants; evaluate derivative of observable for all determinants
   *  @param[out] dvals_O. (O D[i][j]/D[i][j]) for MultiDiracDet i, excited det j (also includes GS det at j==0)
   *  @return 
   */
  void computeMDDerivatives_Obs(const std::vector<ValueMatrix>& Minv_Mv,
                                const std::vector<ValueMatrix>& Minv_B,
                                const std::vector<IndexType>& mdd_spo_ids,
                                const std::vector<const WaveFunctionComponent*>& mdds,
                                std::vector<ValueVector>& dvals_O) const;

  /** @brief Uses per-det lists of derivatives to calculate total excited det contributions
   *
   *  @param[in] mdds. vector of MultiDiracDeterminants; evaluate derivative of observable for all determinants
   *  @param[in] dvals_dmu_O. d/dmu(O D[i][j]/D[i][j]) for MultiDiracDet i, excited det j (also includes GS det at j==0)
   *  @param[in] dvals_O. (O D[i][j]/D[i][j]) for MultiDiracDet i, excited det j (also includes GS det at j==0)
   *  @param[in] dvals_dmu. d/dmu(log(D[i][j]) for MultiDiracDet i, excited det j (also includes GS det at j==0)
   *  @return  {d_mu(OPsi/Psi), d_mu(log(Psi)), OPsi/Psi}
   */
  std::tuple<TWFFastDerivWrapper::ValueType, TWFFastDerivWrapper::ValueType, TWFFastDerivWrapper::ValueType>
      computeMDDerivatives_total(const std::vector<const WaveFunctionComponent*>& mdds,
                                 const std::vector<ValueVector>& dvals_dmu_O,
                                 const std::vector<ValueVector>& dvals_O,
                                 const std::vector<ValueVector>& dvals_dmu) const;

  //////////////////////////////////////////////////////////////////////////////////////////////
  //And now we just have some helper functions for doing useful math on our lists of matrices.//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /** @brief Helper function that inverts all slater matrices in our species list.
   *
   *  @param[in] M. List of slater matrices for each species.  These are square and consistent with some occupation.
   *  @param[in,out] Minv. The species by species list of inverted matrices from M.
   *  @return Void.
   */
  void invertMatrices(const std::vector<ValueMatrix>& M, std::vector<ValueMatrix>& Minv);

  /** @brief Build the auxiliary X=M^-1 B M^-1 matrix.
   *
   *  @param[in] Minv. List of slater matrix inverses M^-1 for a given occupation. 
   *  @param[in] B. Observable auxiliary matrix for a given occupation.
   *  @param[in,out] X. M^-1*B*M^-1 is stored in this list of matrices.
   *  @return Void. 
   */
  void buildX(const std::vector<ValueMatrix>& Minv, const std::vector<ValueMatrix>& B, std::vector<ValueMatrix>& X);

  /** @brief Build auxiliary matrices.
   *
   *  @param[in] Minv. List of slater matrix inverses M^-1 for a given occupation. 
   *  @param[in] B. Observable auxiliary matrix for a given occupation.
   *  @param[in] M. List of slater matrices (rectangular, including virtuals)
   *  @param[out] X. M^-1*B*M^-1 is stored in this list of matrices.
   *  @param[out] Minv_B. M^-1*B is stored in this list of matrices.
   *  @param[out] Minv_Mv. M^-1*Mvirt is stored in this list of matrices.
   *  @return Void. 
   */
  void buildIntermediates(const std::vector<ValueMatrix>& Minv,
                          const std::vector<ValueMatrix>& B,
                          const std::vector<ValueMatrix>& M,
                          std::vector<ValueMatrix>& X,
                          std::vector<ValueMatrix>& Minv_B,
                          std::vector<ValueMatrix>& Minv_Mv);

  /** @brief Build auxiliary matrices for ion derivs
   * first dim is x/y/z for dB/dM/Minv_dB/Minv_dM
   *
   *  @param[in] Minv. List of slater matrix inverses M^-1 for a given occupation. 
   *  @param[in] dB. list of derivative of observable auxiliary matrix for a given occupation.
   *  @param[in] dM. List of derivative of slater matrices (rectangular, including virtuals)
   *  @param[out] Minv_dB. M^-1*dB is stored in this list of matrices.
   *  @param[out] Minv_dM. M^-1*dM is stored in this list of matrices.
   *  @return Void. 
   */
  void buildIntermediates_dmu(const std::vector<ValueMatrix>& Minv,
                              const std::vector<std::vector<ValueMatrix>>& dB,
                              const std::vector<std::vector<ValueMatrix>>& dM,
                              std::vector<std::vector<ValueMatrix>>& Minv_dB,
                              std::vector<std::vector<ValueMatrix>>& Minv_dM);

  /** @brief Goes through a list of matrices and zeros them out.  Does this in place.
   *
   *  @param[in,out] A. The list of matrices to be zeroed out.  After call, A is all zeros.
   *  @return Void.
   */
  void wipeMatrices(std::vector<ValueMatrix>& A);

  /** @brief Goes through a list of vectors and zeros them out.  Does this in place.
   *
   *  @param[in,out] A. The list of vectors to be zeroed out.  After call, A is all zeros.
   *  @return Void.
   */
  void wipeVectors(std::vector<ValueVector>& A);

  /** @brief Returns Tr(A*B).  Actually, we assume a block diagonal structure, so this is 
   *    really Sum_i Tr(A_i*B_i), where i is the species index.
   * 
   *  @param[in] A.  Vector of matrices A.  
   *  @param[in] B.  Vector of matrices B.
   *  @return Value of Sum_i Tr(A_i*B_i).
   */
  ValueType trAB(const std::vector<ValueMatrix>& A, const std::vector<ValueMatrix>& B);

  /**
   * @brief transform X = A[o,v] - A[o,o].B[o,v]
   * 
   * @param[in] A matrix with dims [o,o+v]
   * @param[in] B matrix with dims [o,v]
   * @param[out] X matrix with dims [o,v]
   */
  void transform_Av_AoBv(const ValueMatrix& A, const ValueMatrix& B, ValueMatrix& X) const;

  // Minimal resource interface so TWF can call it like other components
  static void createResource(ResourceCollection& collection);

  static void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<TWFFastDerivWrapper>& wrappers);

  static void releaseResource(ResourceCollection&, const RefVectorWithLeader<TWFFastDerivWrapper>& wrappers);


private:
  std::vector<SPOSet*> spos_;
  std::vector<IndexType> groups_;
  std::vector<ValueMatrix> psi_M_;
  std::vector<ValueMatrix> psi_M_inv_;
  std::vector<WaveFunctionComponent*> jastrow_list_;
  // pointer to MultiSlaterDetTableMethod if one has been registered, otherwise nullptr
  // access constituent MultiDiracDets and SPOsets through this slaterdet (associate with spos_ via group ID)
  const WaveFunctionComponent* multislaterdet_ = nullptr;

  struct TWFFastDerivWrapperMultiWalkerMem;
  ResourceHandle<TWFFastDerivWrapperMultiWalkerMem> mw_mem_handle_;
};

/**@}*/
} // namespace qmcplusplus
#endif
