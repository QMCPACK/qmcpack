//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SINGLEPARTICLEORBITALSETBASE_H
#define QMCPLUSPLUS_SINGLEPARTICLEORBITALSETBASE_H

#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"
#include "Particle/VirtualParticleSet.h"

#include "io/hdf_archive.h"
#if !defined(ENABLE_SOA)
#include "Message/CommOperators.h"
#endif
#include "QMCWaveFunctions/SPOSetTypeAliases.h"

#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

#include "Batching.h"
#include "WhatAmI.h"
namespace qmcplusplus
{

/** base class for Single-particle orbital sets
 *
 * SPOSet stands for S(ingle)P(article)O(rbital)SetBase which contains
 * a number of single-particle orbitals with capabilities of evaluating \f$ \psi_j({\bf r}_i)\f$
 */

template<Batching batching = Batching::SINGLE>
class SPOSet;

template<>
class SPOSet<Batching::SINGLE> : public virtual WhatAmI
{
public:
    // Because so many SPOSet<Batching::SINGLE> children exist and
  //expect to inherit typedefs
  using SSTA = SPOSetTypeAliases;
  using QMCT = QMCTraits;
  using ValueType = QMCTraits::ValueType;
  using IndexVector_t = SSTA::IndexVector_t;
  using ValueVector_t = SSTA::ValueVector_t;
  using ValueAlignedVector_t = SSTA::ValueAlignedVector_t;
  using ValueMatrix_t = SSTA::ValueMatrix_t;
  using GradVector_t = SSTA::GradVector_t;
  using GradMatrix_t = SSTA::GradMatrix_t;
  using HessVector_t = SSTA::HessVector_t;
  using HessMatrix_t = SSTA::HessMatrix_t;
  using HessType = SSTA::HessType;
  using HessArray_t = SSTA::HessArray_t;
  using GradHessType = SSTA::GGGType;
  //using GradHessVector_t = SSTA::GGGVector_t;
  using GGGVector_t = SSTA::GGGVector_t;
  //using GradHessMatrix_t = SSTA::GGGMatrix_t;
  using GGGMatrix_t = SSTA::GGGMatrix_t;
  using VGLVector_t = SSTA::VGLVector_t;
  using Walker_t = SSTA::Walker_t;
  using GGGType = SSTA::GGGType;
  typedef std::map<std::string,SPOSet<>*> SPOPool_t;

  static constexpr int DIM = 3;
  
  ///index in the builder list of sposets
  int builder_index;
  ///true if SPO is optimizable
  bool Optimizable;

  ///flag to calculate ionic derivatives
  bool ionDerivs;
  ///number of Single-particle orbitals
  QMCT::IndexType OrbitalSetSize;
  /// Optimizable variables
  opt_variables_type myVars;
  ///name of the basis set
  std::string className;
  /** name of the object
   *
   * Several user classes can own SPOSet and use objectName as counter
   */
  std::string objectName;
#if !defined(ENABLE_SOA)
  ///true if C is an identity matrix
  bool Identity;
  ///if true, do not clean up
  bool IsCloned;
  ///number of Single-particle orbitals
  IndexType BasisSetSize;
  /** pointer matrix containing the coefficients
   *
   * makeClone makes a shallow copy
   */
  SSTA::ValueMatrix_t* C;
  ///occupation number
  Vector<RealType> Occ;
  ///Pass Communicator
  Communicate *myComm;
#endif
  
  /** constructor */
  SPOSet();

  /** destructor */
  virtual ~SPOSet()
  {
#if !defined(ENABLE_SOA)
    if(!IsCloned && C!= nullptr) delete C;
#endif
  }

  /** return the size of the orbital set
   */
  inline int size() const
  {
    return OrbitalSetSize;
  }

  /** print basic SPOSet information
   */
  void basic_report(const std::string& pad="");
  
  /** print SPOSet information
   */
  virtual void report(const std::string& pad="")
  {
    basic_report(pad);
  }


  /** return the size of the orbital set
   */
  inline int getOrbitalSetSize() const
  {
    return OrbitalSetSize;
  }


#if !defined(ENABLE_SOA)
  int getBasisSetSize() const { return BasisSetSize; }

  bool setIdentity(bool useIdentity);

  void checkObject();

  ///get C and Occ
  bool put(xmlNodePtr cur);
#else
  virtual int getBasisSetSize() const { return 0; }

  virtual void checkObject() const {}
#endif

  ///reset
  virtual void resetParameters(const opt_variables_type& optVariables)=0;

  virtual void checkInVariables(opt_variables_type& active) {}
  virtual void checkOutVariables(const opt_variables_type& active) {}

  // Evaluate the derivative of the optimized orbitals with
  // respect to the parameters
  virtual void evaluateDerivatives
  (ParticleSet& P, int iat, const opt_variables_type& active,
   SSTA::ValueMatrix_t& d_phi, SSTA::ValueMatrix_t& d_lapl_phi) {}


  ///reset the target particleset
  virtual void resetTargetParticleSet(ParticleSet& P)=0;
  /** set the OrbitalSetSize
   * @param norbs number of single-particle orbitals
   */
  virtual void setOrbitalSetSize(int norbs)=0;

  /** estimate the memory needs for evaluating SPOs of particles in the size of ValueType
   * @param nP, number of particles.
   */
  virtual size_t
  estimateMemory(const int nP) { return 0; }


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  returns whether this is an LCOrbitalSetOpt object
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual bool is_of_type_LCOrbitalSetOpt() const { return false; }

  virtual QMCT::PosType get_k(int orb)
  {
    return QMCT::PosType();
  }

  /** make a clone of itself
   */
  virtual SPOSet<>* makeClone() const;

  virtual bool transformSPOSet()
  {
    return true;
  }

  // Routine to set up data for the LCOrbitalSetOpt child class specifically
  // Should be left empty for other derived classes
  virtual void init_LCOrbitalSetOpt(const double mix_factor=0.0) { };

  // Routine to update internal data for the LCOrbitalSetOpt child class specifically
  // Should be left empty for other derived classes
  virtual void rotate_B(const std::vector<QMCT::RealType> &rot_mat) { };

#if !defined(ENABLE_SOA)
protected:
  bool putOccupation(xmlNodePtr occ_ptr);
  bool putFromXML(xmlNodePtr coeff_ptr);
  bool putFromH5(const char* fname, xmlNodePtr coeff_ptr);
#endif

  virtual void
  evaluate (const ParticleSet& P, QMCT::PosType &r, SSTA::ValueVector_t &psi)
  {
    app_error() << "Need specialization for SPOSet::evaluate "
                << "(const ParticleSet& P, PosType &r).\n";
    abort();
  }

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat, SSTA::ValueVector_t& psi)=0;

  /** evaluate VGL of SPOs using SoA container for gl
   */
  virtual void
  evaluateVGL(const ParticleSet& P, int iat, SSTA::VGLVector_t& vgl);

  /** evaluate values for the virtual moves, e.g., sphere move for nonlocalPP
   * @param VP virtual particle set
   * @param psiM single-particle orbitals psiM(i,j) for the i-th particle and the j-th orbital
   * @param SPOMem scratch space for SPO value evaluation, alignment is required.
   */
  virtual void
  evaluateValues(const VirtualParticleSet& VP, SSTA::ValueMatrix_t& psiM, SSTA::ValueAlignedVector_t& SPOMem);


  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::ValueVector_t& d2psi)=0;

  /** evaluate the values, gradients and hessians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::HessVector_t& grad_grad_psi)=0;

  virtual void
  evaluateThirdDeriv(const ParticleSet& P, int first, int last
                     , SSTA::GGGMatrix_t& grad_grad_grad_logdet);

  
  /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last) particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param d2logdet laplacians
   *
   */
  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::ValueMatrix_t& d2logdet)=0;

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet);

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet, SSTA::GGGMatrix_t& grad_grad_grad_logdet);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src, SSTA::GradMatrix_t &gradphi);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src
                                   , SSTA::GradMatrix_t &grad_phi, SSTA::HessMatrix_t &grad_grad_phi, SSTA::GradMatrix_t &grad_lapl_phi);

};

#if defined(ENABLE_SMARTPOINTER)
typedef boost::shared_ptr<SPOSet> SPOSetPtr;
#else
  //typedef SPOSet* SPOSetPtr;
#endif

}
#endif


