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

#include <memory>
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"

#include "Numerics/MatrixOperators.h"
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
  using basis_type = SoaBasisSetBase<ValueType>;
  using vgl_type   = basis_type::vgl_type;
  using vgh_type   = basis_type::vgh_type;
  using vghgh_type = basis_type::vghgh_type;

  ///pointer to the basis set
  std::unique_ptr<basis_type> myBasisSet;
  /// pointer to matrix containing the coefficients
  std::shared_ptr<ValueMatrix> C;

  /** constructor
     * @param bs pointer to the BasisSet
     */
  LCAOrbitalSet(std::unique_ptr<basis_type>&& bs, bool optimize);

  LCAOrbitalSet(const LCAOrbitalSet& in);

  std::unique_ptr<SPOSet> makeClone() const override;

  void storeParamsBeforeRotation() override { C_copy = *C; }

  void applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy) override;

  void checkInVariables(opt_variables_type& active) override
  {
    APP_ABORT("LCAOrbitalSet should not call checkInVariables");
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    APP_ABORT("LCAOrbitalSet should not call checkOutVariables");
  }

  ///reset
  void resetParameters(const opt_variables_type& active) override
  {
    APP_ABORT("LCAOrbitalSet should not call resetParameters");
  }

  /** set the OrbitalSetSize and Identity=false and initialize internal storages
    */
  void setOrbitalSetSize(int norbs) override;

  /** return the size of the basis set
    */
  int getBasisSetSize() const override { return (myBasisSet == nullptr) ? 0 : myBasisSet->getBasisSetSize(); }

  bool isIdentity() const { return Identity; };

  /** check consistency between Identity and C
    *
    */
  void checkObject() const override;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<ValueType>& ratios) override;

  void evaluateVGH(const ParticleSet& P,
                   int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override;

  void evaluateVGHGH(const ParticleSet& P,
                     int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) override;

  //NOTE:  The data types get complicated here, so here's an overview of the
  //       data types associated with ionic derivatives, and how to get their data.
  //
  //NOTE:  These data structures hold the data for one particular ion, and so the ID is implicit.
  //       It's up to the user to keep track of which ion these derivatives refer to.
  //
  // 1.) GradMatrix grad_phi:  Holds the ionic derivatives of each SPO for each electron.
  //            Example:  grad_phi[iel][iorb][idim].  iel  -- electron index.
  //                                                iorb -- orbital index.
  //                                                idim  -- cartesian index of ionic derivative.
  //                                                        X=0, Y=1, Z=2.
  //
  // 2.) HessMatrix grad_grad_phi:  Holds the ionic derivatives of the electron gradient components
  //                                   for each SPO and each electron.
  //            Example:  grad_grad_phi[iel][iorb](idim,edim)  iel  -- electron index.
  //                                                           iorb -- orbital index.
  //                                                           idim -- ionic derivative's cartesian index.
  //                                                              X=0, Y=1, Z=2
  //                                                           edim -- electron derivative's cartesian index.
  //                                                              x=0, y=1, z=2.
  //
  // 3.) GradMatrix grad_lapl_phi:  Holds the ionic derivatives of the electron laplacian for each SPO and each electron.
  //            Example:  grad_lapl_phi[iel][iorb][idim].  iel  -- electron index.
  //                                                       iorb -- orbital index.
  //                                                       idim -- cartesian index of ionic derivative.
  //                                                           X=0, Y=1, Z=2.

  /**
 * \brief Calculate ion derivatives of SPO's.
 *  
 *  @param P Electron particle set.
 *  @param first index of first electron 
 *  @@param last index of last electron
 *  @param source Ion particle set.
 *  @param iat_src  Index of ion.
 *  @param gradphi Container storing ion gradients for all particles and all orbitals.
 */
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi) override;

  /**
 * \brief Calculate ion derivatives of SPO's, their gradients, and their laplacians.
 *  
 *  @param P Electron particle set.
 *  @param first index of first electron.
 *  @@param last index of last electron
 *  @param source Ion particle set.
 *  @param iat_src  Index of ion.
 *  @param grad_phi Container storing ion gradients for all particles and all orbitals.
 *  @param grad_grad_phi Container storing ion gradients of electron gradients for all particles and all orbitals.
 *  @param grad_lapl_phi Container storing ion gradients of SPO laplacians for all particles and all orbitals.
 */
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi,
                          HessMatrix& grad_grad_phi,
                          GradMatrix& grad_lapl_phi) override;

  void evaluateGradSourceRow(const ParticleSet& P,
                             int iel,
                             const ParticleSet& source,
                             int iat_src,
                             GradVector& grad_phi) override;

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix& grad_grad_grad_logdet) override;

protected:
  ///number of Single-particle orbitals
  const IndexType BasisSetSize;
  /// a copy of the original C before orbital rotation is applied;
  ValueMatrix C_copy;

  ///true if C is an identity matrix
  bool Identity;
  ///Temp(BasisSetSize) : Row index=V,Gx,Gy,Gz,L
  vgl_type Temp;
  ///Tempv(OrbitalSetSize) Tempv=C*Temp
  vgl_type Tempv;

  ///These are temporary VectorSoAContainers to hold value, gradient, and hessian for
  ///all basis or SPO functions evaluated at a given point.
  ///Nbasis x [1(value)+3(gradient)+6(hessian)]
  vgh_type Temph;
  ///Norbitals x [1(value)+3(gradient)+6(hessian)]
  vgh_type Temphv;

  ///These are temporary VectorSoAContainers to hold value, gradient, hessian, and
  /// gradient hessian for all basis or SPO functions evaluated at a given point.
  ///Nbasis x [1(value)+3(gradient)+6(hessian)+10(grad_hessian)]
  vghgh_type Tempgh;
  ///Nbasis x [1(value)+3(gradient)+6(hessian)+10(grad_hessian)]
  vghgh_type Tempghv;

private:
  ///helper functions to handl Identity
  void evaluate_vgl_impl(const vgl_type& temp, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) const;

  void evaluate_vgl_impl(const vgl_type& temp,
                         int i,
                         ValueMatrix& logdet,
                         GradMatrix& dlogdet,
                         ValueMatrix& d2logdet) const;
  ///These two functions unpack the data in vgh_type temp object into wavefunction friendly data structures.
  ///This unpacks temp into vectors psi, dpsi, and d2psi.
  void evaluate_vgh_impl(const vgh_type& temp, ValueVector& psi, GradVector& dpsi, HessVector& d2psi) const;

  ///Unpacks temp into the ith row (or electron index) of logdet, dlogdet, dhlogdet.
  void evaluate_vgh_impl(const vgh_type& temp,
                         int i,
                         ValueMatrix& logdet,
                         GradMatrix& dlogdet,
                         HessMatrix& dhlogdet) const;
  ///Unpacks data in vghgh_type temp object into wavefunction friendly data structures for value, gradient, hessian
  ///and gradient hessian.
  void evaluate_vghgh_impl(const vghgh_type& temp,
                           ValueVector& psi,
                           GradVector& dpsi,
                           HessVector& d2psi,
                           GGGVector& dghpsi) const;

  void evaluate_vghgh_impl(const vghgh_type& temp,
                           int i,
                           ValueMatrix& logdet,
                           GradMatrix& dlogdet,
                           HessMatrix& dhlogdet,
                           GGGMatrix& dghlogdet) const;


  ///Unpacks data in vgl object and calculates/places ionic gradient result into dlogdet.
  void evaluate_ionderiv_v_impl(const vgl_type& temp, int i, GradMatrix& dlogdet) const;

  ///Unpacks data in vgl object and calculates/places ionic gradient of value,
  ///  electron gradient, and electron laplacian result into dlogdet, dglogdet, and dllogdet respectively.
  void evaluate_ionderiv_vgl_impl(const vghgh_type& temp,
                                  int i,
                                  GradMatrix& dlogdet,
                                  HessMatrix& dglogdet,
                                  GradMatrix& dllogdet) const;

  ///Unpacks data in vgl object and calculates/places ionic gradient of a single row (phi_j(r)) into dlogdet.
  void evaluate_ionderiv_v_row_impl(const vgl_type& temp, GradVector& dlogdet) const;
};
} // namespace qmcplusplus
#endif
