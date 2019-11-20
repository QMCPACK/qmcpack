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
  ~SpinorSet() =  default;
 
  //This class is initialized by separately building the up and down channels of the spinor set and
  //then registering them.  
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
  /** Evaluate the values, spin gradients, and spin laplacians of single particle spinors corresponding to electron iat.
   *  @param P current particle set.
   *  @param iat electron index.
   *  @param spinor values.
   *  @param spin gradient values. d/ds phi(r,s).
   *  @param spin laplacian values. d^2/ds^2 phi(r,s). 
   *
   */
  void evaluate_spin(const ParticleSet& P,
                        int iat,
                        ValueVector_t& psi,
                        ValueVector_t& dpsi,
                        ValueVector_t& d2psi);
                   

private:
  //Sposet for the up and down channels of our spinors.  
  std::shared_ptr<SPOSet> spo_up;
  std::shared_ptr<SPOSet> spo_dn;

  //temporary arrays for holding the values of the up and down channels respectively.
  ValueVector_t psi_work_up;
  ValueVector_t psi_work_down;
  
  //temporary arrays for holding the gradients of the up and down channels respectively.
  GradVector_t dpsi_work_up;
  GradVector_t dpsi_work_down;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  ValueVector_t d2psi_work_up;
  ValueVector_t d2psi_work_down;

  //Same as above, but these are the full matrices containing all spinor/particle combinations.
  ValueMatrix_t logpsi_work_up;
  ValueMatrix_t logpsi_work_down;
  
  GradMatrix_t dlogpsi_work_up;
  GradMatrix_t dlogpsi_work_down;
  
  ValueMatrix_t d2logpsi_work_up;
  ValueMatrix_t d2logpsi_work_down;
  
};

}
#endif
