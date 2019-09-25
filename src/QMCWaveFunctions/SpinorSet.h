//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
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
  ~SpinorSet();
 
  void set_spos(std::shared_ptr<SPOSet> up, std::shared_ptr<SPOSet> dn);
  /// reset parameters to the values from optimizer
  void resetParameters(const opt_variables_type& optVariables) override;

  /** reset the target particleset
   *  this is used to reset the pointer to ion-electron distance table needed by LCAO basis set.
   *  Ye: Only AoS needs it, SoA LCAO doesn't need this. Reseting pointers is a state machine very hard to maintain.
   *  This interface should be removed with AOS.
   */
  void resetTargetParticleSet(ParticleSet& P) override;

  /** set the OrbitalSetSize
   * @param norbs number of single-particle orbitals
   * Ye: I prefer to remove this interface in the future. SPOSet builders need to handle the size correctly.
   * It doesn't make sense allowing to set the value at any place in the code.
   */
  void setOrbitalSetSize(int norbs) override;


  /** evaluate the values of this spinor set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi) override;

  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   * @param dpsi gradients of the SPO
   * @param d2psi laplacians of the SPO
   */
  void evaluate(const ParticleSet& P,
                        int iat,
                        ValueVector_t& psi,
                        GradVector_t& dpsi,
                        ValueVector_t& d2psi) override;

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
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    ValueMatrix_t& d2logdet) override;

  void evaluate_spin(const ParticleSet& P,
                        int iat,
                        ValueVector_t& psi,
                        ValueVector_t& dpsi,
                        ValueVector_t& d2psi);
                   

private:
  std::shared_ptr<SPOSet> spo_up;
  std::shared_ptr<SPOSet> spo_dn;

  ValueVector_t psi_work_up;
  ValueVector_t psi_work_down;
  
  GradVector_t dpsi_work_up;
  GradVector_t dpsi_work_down;

  ValueVector_t d2psi_work_up;
  ValueVector_t d2psi_work_down;

  ValueMatrix_t logpsi_work_up;
  ValueMatrix_t logpsi_work_down;
  
  GradMatrix_t dlogpsi_work_up;
  GradMatrix_t dlogpsi_work_down;
  
  ValueMatrix_t d2logpsi_work_up;
  ValueMatrix_t d2logpsi_work_down;
  
};

}
#endif
