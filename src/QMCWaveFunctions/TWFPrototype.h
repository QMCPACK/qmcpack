//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


/**@file TWFPrototype.h
 *@brief Declaration of TWFPrototype
 */
#ifndef QMCPLUSPLUS_TWFPROTOTYPE_H
#define QMCPLUSPLUS_TWFPROTOTYPE_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"
namespace qmcplusplus
{
/**
 *  TWFPrototype is a wrapper class for TrialWavefunction that provides separate and low level access to the Jastrow and 
 *  SPOSet objects.  This is so that observables can be recast in matrix form and their derivatives taken efficiently. 
 *  Currently this is hard coded for ground state slater jastrow wave functions, but generalization to
 *  arbitrary occupations and multideterminants are straightforward and will come online as the code is tested and validated. 
 *
 *  Please see : J. Chem. Phys. 144, 194105 (2016) https://doi.org/10.1063/1.4948778 for implementation details and formalism.
 */
class TWFPrototype
{
public:
  using ValueMatrix_t = SPOSet::ValueMatrix_t;
  using GradMatrix_t  = SPOSet::GradMatrix_t;
  using HessMatrix_t  = SPOSet::HessMatrix_t;
  using IndexType     = QMCTraits::IndexType;
  using RealType      = QMCTraits::RealType;
  using ValueType     = QMCTraits::ValueType;

  using ValueVector_t = SPOSet::ValueVector_t;
  using GradVector_t  = SPOSet::GradVector_t;

  TWFPrototype();

  
  /** @brief Add a particle group.
   *
   *  Here, a "group" corresponds to a subset of particles which are antisymmetric with 
   *  respect to eachother.  TWFPrototype ensures that there is a binding between the groupid
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
  inline void add_jastrow(WaveFunctionComponent* j) { J.push_back(j); };

  /** @brief Takes particle set groupID and returns the TWF internal index for it.  
   *
   *  ParticleSet groups can be registered in whichever order.  However, the internal indexing 
   *  of TWFPrototype keeps track on a first-come, first serve basis.  That is, if I register 
   *  particle groups 3, 1, and 0 in that order, then the internal indexing goes like
   *  0->3, 1->1, 2->0.  Thus, this function looks up where in the internal indexing scheme
   *  ParticleSet gid is located.  This is necessary, since the matrix list data structures follow
   *  the internal TWF indexing.  
   *
   *  @param[in] gid. ParticleSet group ID to look up.  
   *  @return The corresponding internal groupID.  
   */  
  IndexType getTWFGroupIndex(const IndexType gid);

  /** @ingroup Query functions
   * @{
   * @brief These are query functions to the internal state of various groups and SPOSet info.
   *   All group indices will refer to the internal group indices.  Convert from ParticleSet to 
   *   TWF group first!
   *
   *   Source of truth for orbital sizes will be the individual SPOSets.  Particle group sizes
   *   will be ParticleSet in conjunction with groupID maps.  
   */
  inline IndexType numGroups() { return spos.size(); };
  SPOSet* getSPOSet(const IndexType sid) { return spos[sid]; };
  inline IndexType numOrbitals(const IndexType sid) { return spos[sid]->size(); };
  inline IndexType numParticles(const IndexType sid) { return num_ptcls[sid]; }; 
  /** @} */

  /** @brief Returns log(Psi).  Should be consistent with QMCPACK unwrapped TrialWavefunction.
   *
   *  @param[in] P. Particle set.
   *  @return log(Psi).
   */
  RealType evaluateLog(ParticleSet& P);

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
  void get_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mmat);

  /** @brief Returns value of all orbitals (relevant to given species/group) at a particular particle coordinate.
   *
   *  @param[in] P particle set.
   *  @param[in]**** iel particle ID.
   *  @param[in,out] val Vector of phi_i(r_iel) for all i=0,Norbs.
   *  @return Void
   */
  IndexType get_M_row(const ParticleSet& P, IndexType iel, ValueVector_t& val);

  /** @brief Returns the ion derivative of all orbitals (relevant to given species/group) at a particular particle coordinate.
   *
   *  @param[in] P particle set.
   *  @param[in] source ion particle set.
   *  @param[in] iel particle ID.
   *  @param[in]m iat_source ion index w.r.t. which the ion derivative is taken
   *  @param[in,out]  dval Vector of d/dR_iat_source phi_i(r_iel) for all i=0,Norbs. First index is X derivative, second Y derivative, etc.
   *  @return Void
   */
  IndexType get_igrad_row(const ParticleSet& P,
                          const ParticleSet& source,
                          IndexType iel,
                          IndexType iat_source,
                          std::vector<ValueVector_t>& dval);

  /** @brief Returns value, gradient, and laplacian matrices for all orbitals and all particles, species by species. 
   *
   *  @param[in] P particle set.
   *  @param[in,out] mvec Slater matrix M_ij=phi_j(r_i) for each species group.
   *  @param[in,out] gmat electron gradient of slater matrix [G_ij]_a = d/dr_a,i phi_j(r_i).  a=x,y,z  
   *  @param[in,out] lmat electron laplacian of slater matrix [L_ij] = \nabla^2_i phi_j(r_i).
   *  @return Void
   */
  void get_egrad_elapl_M(const ParticleSet& P,
                         std::vector<ValueMatrix_t>& mvec,
                         std::vector<GradMatrix_t>& gmat,
                         std::vector<ValueMatrix_t>& lmat);

  /** @brief Returns x,y,z components of ion gradient of slater matrices.
   *
   *  @param[in] P particle set.
   *  @param[in] source ion particle set.
   *  @param[in]*** iat ion ID w.r.t. which to take derivative.
   *  @param[in,out] dmvec Slater matrix d/dR_{iat,a} M_ij=d/dR_{iat,a} phi_j(r_i) for each species group.  First index is a=x,y,z.
   *  @return Void
   */
  void get_igrad_M(const ParticleSet& P,
                   const ParticleSet& source,
                   int iat,
                   std::vector<std::vector<ValueMatrix_t>>& dmvec);

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
  void get_igrad_igradelapl_M(const ParticleSet& P,
                              const ParticleSet& source,
                              int iat,
                              std::vector<std::vector<ValueMatrix_t>>& dmvec,
                              std::vector<std::vector<ValueMatrix_t>>& dlmat);


  /** @brief Takes sub matrices of full SPOSet quantities (evaluated on all particles and all orbitals), consistent with ground
   *   state occupations.
   *
   *  @param[in] A non-rectangular matrices of SPOSet derived quantities, evaluated on all orbitals and all particles.
   *  @param[in,out] Aslice square matrices consistent with a ground state occupation.
   *  @return Void
   */
  void get_gs_matrix(const std::vector<ValueMatrix_t>& A, std::vector<ValueMatrix_t>& Aslice);

  /** @brief Calculates derivative of observable via Tr[M^{-1} dB - X * dM ].  Consistent with ground state occupation.
   *
   *  @param[in] Minv. inverse of slater matrices for ground state occupations. 
   *  @param[in] X.  X=M^-1 B M^-1.  B observable matrix, and is consistent with ground state occupation.
   *  @param[in] dM. Target derivative of M, and is consistent with ground state occupation.
   *  @param[in] dB. Target derivative of B, and is consistent with ground state occupation.
   *  @return Derivative of O psi/psi = Tr[M^{-1} dB - X * dM ]
   */
  ValueType compute_gs_derivative(const std::vector<ValueMatrix_t>& Minv,
                                  const std::vector<ValueMatrix_t>& X,
                                  const std::vector<ValueMatrix_t>& dM,
                                  const std::vector<ValueMatrix_t>& dB);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //And now we just have some helper functions for doing useful math on our lists of matrices.//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /** @brief Helper function that inverts all slater matrices in our species list.
   *
   *  @param[in] M. List of slater matrices for each species.  These are square and consistent with some occupation.
   *  @param[in,out] Minv. The species by species list of inverted matrices from M.
   *  @return Void.
   */
  void invert_M(const std::vector<ValueMatrix_t>& M, std::vector<ValueMatrix_t>& Minv);

  /** @brief Helper function that inverts all slater matrices in our species list.
   *
   *  @param[in] Minv. List of slater matrix inverses M^-1 for a given occupation. 
   *  @param[in] B. Observable auxiliary matrix for a given occupation.
   *  @param[in,out] X. M^-1*B*M^-1 is stored in this list of matrices.
   *  @return Void. 
   */
  void build_X(const std::vector<ValueMatrix_t>& Minv,
               const std::vector<ValueMatrix_t>& B,
               std::vector<ValueMatrix_t>& X);

  /** @brief Goes through a list of matrices and zeros them out.  Does this in place.
   *
   *  @param[in,out] A. The list of matrices to be zeroed out.  After call, A is all zeros.
   *  @return Void.
   */
  void wipe_matrix(std::vector<ValueMatrix_t>& A);

  /** @brief Returns Tr(A*B).  Actually, we assume a block diagonal structure, so this is 
   *    really Sum_i Tr(A_i*B_i), where i is the species index.
   * 
   *  @param[in] A.  Vector of matrices A.  
   *  @param[in] B.  Vector of matrices B.
   *  @return Value of Sum_i Tr(A_i*B_i).
   */
  ValueType trAB(const std::vector<ValueMatrix_t>& A, const std::vector<ValueMatrix_t>& B);


private:
  std::vector<IndexType> num_ptcls;
  std::vector<IndexType> num_orbs;
  std::vector<SPOSet*> spos;
  std::vector<IndexType> groups;
  std::vector<ValueMatrix_t> psiM;
  std::vector<ValueMatrix_t> psiMinv;
  std::vector<WaveFunctionComponent*> J;

  bool initialized;
};

/**@}*/
} // namespace qmcplusplus
#endif
