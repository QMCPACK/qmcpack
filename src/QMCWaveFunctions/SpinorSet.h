//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by:  Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPINORSET_H
#define QMCPLUSPLUS_SPINORSET_H

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{
/** Class for Melton & Mitas style Spinors.
 *
 */
class SpinorSet : public SPOSet
{
public:
  ///name of the class
  std::string className;

  /** constructor */
  SpinorSet();
  ~SpinorSet() override = default;

  //This class is initialized by separately building the up and down channels of the spinor set and
  //then registering them.
  void set_spos(std::unique_ptr<SPOSet>&& up, std::unique_ptr<SPOSet>&& dn);
  /// reset parameters to the values from optimizer
  void resetParameters(const opt_variables_type& optVariables) override;

  /** set the OrbitalSetSize
   * @param norbs number of single-particle orbitals
   */
  void setOrbitalSetSize(int norbs) override;

  //gets the BasisSetSize from the underlying SPOSet that make up the spinor
  int getBasisSetSize() const override;


  /** evaluate the values of this spinor set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   * @param dpsi gradients of the SPO
   * @param d2psi laplacians of the SPO
   */
  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   * @param dpsi gradients of the SPO
   * @param d2psi laplacians of the SPO
   * @param dspin spin gradient of the SPO
   */
  void evaluateVGL_spin(const ParticleSet& P,
                        int iat,
                        ValueVector& psi,
                        GradVector& dpsi,
                        ValueVector& d2psi,
                        ValueVector& dspin) override;

  /** evaluate the values, gradients and laplacians and spin gradient of this single-particle orbital sets of multiple walkers
   * @param spo_list the list of SPOSet pointers in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat active particle
   * @param psi_v_list the list of value vector pointers in a walker batch
   * @param dpsi_v_list the list of gradient vector pointers in a walker batch
   * @param d2psi_v_list the list of laplacian vector pointers in a walker batch
   * @param dspin_v_list the list of spin gradients vector pointers in a walker batch
   */
  void mw_evaluateVGLWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector>& psi_v_list,
                              const RefVector<GradVector>& dpsi_v_list,
                              const RefVector<ValueVector>& d2psi_v_list,
                              const RefVector<ValueVector>& dspin_v_list) const override;

  /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last) particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param d2logdet laplacians
   *
   */
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

  void mw_evaluate_notranspose(const RefVectorWithLeader<SPOSet>& spo_list,
                               const RefVectorWithLeader<ParticleSet>& P_list,
                               int first,
                               int last,
                               const RefVector<ValueMatrix>& logdet_list,
                               const RefVector<GradMatrix>& dlogdet_list,
                               const RefVector<ValueMatrix>& d2logdet_list) const override;

  void evaluate_notranspose_spin(const ParticleSet& P,
                                 int first,
                                 int last,
                                 ValueMatrix& logdet,
                                 GradMatrix& dlogdet,
                                 ValueMatrix& d2logdet,
                                 ValueMatrix& dspinlogdet) override;
  /** Evaluate the values, spin gradients, and spin laplacians of single particle spinors corresponding to electron iat.
   *  @param P current particle set.
   *  @param iat electron index.
   *  @param spinor values.
   *  @param spin gradient values. d/ds phi(r,s).
   *
   */
  void evaluate_spin(const ParticleSet& P, int iat, ValueVector& psi, ValueVector& dpsi) override;

  std::unique_ptr<SPOSet> makeClone() const override;

private:
  //Sposet for the up and down channels of our spinors.
  std::unique_ptr<SPOSet> spo_up;
  std::unique_ptr<SPOSet> spo_dn;

  //temporary arrays for holding the values of the up and down channels respectively.
  ValueVector psi_work_up;
  ValueVector psi_work_down;

  //temporary arrays for holding the gradients of the up and down channels respectively.
  GradVector dpsi_work_up;
  GradVector dpsi_work_down;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  ValueVector d2psi_work_up;
  ValueVector d2psi_work_down;

  //Same as above, but these are the full matrices containing all spinor/particle combinations.
  ValueMatrix logpsi_work_up;
  ValueMatrix logpsi_work_down;

  GradMatrix dlogpsi_work_up;
  GradMatrix dlogpsi_work_down;

  ValueMatrix d2logpsi_work_up;
  ValueMatrix d2logpsi_work_down;
};

} // namespace qmcplusplus
#endif
