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
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Message/CommOperators.h"

#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

namespace qmcplusplus
{

/** base class for Single-particle orbital sets
 *
 * SPOSetBase stands for S(ingle)P(article)O(rbital)SetBase which contains
 * a number of single-particle orbitals with capabilities of evaluating \f$ \psi_j({\bf r}_i)\f$
 */
class SPOSetBase: public QMCTraits
{
public:
  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
  typedef OrbitalSetTraits<ValueType>::GradMatrix_t  GradMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessVector_t  HessVector_t;
  typedef OrbitalSetTraits<ValueType>::HessMatrix_t  HessMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType      HessType;
  typedef Array<HessType,OHMMS_DIM>                  HessArray_t;
  typedef OrbitalSetTraits<ValueType>::GradHessType  GGGType;
  typedef OrbitalSetTraits<ValueType>::GradHessVector_t GGGVector_t;
  typedef OrbitalSetTraits<ValueType>::GradHessMatrix_t GGGMatrix_t;
  typedef OrbitalSetTraits<ValueType>::RefVector_t      RefVector_t;
  typedef OrbitalSetTraits<ValueType>::VGLVector_t      VGLVector_t;
  typedef ParticleSet::Walker_t                      Walker_t;
  typedef std::map<std::string,SPOSetBase*> SPOPool_t;
  
  ///index in the builder list of sposets
  int builder_index;
  ///true if C is an identity matrix
  bool Identity;
  ///true if SPO is optimizable
  bool Optimizable;
  ///flag to calculate ionic derivatives
  bool ionDerivs;
  ///if true, can use GL type, default=false
  bool CanUseGLCombo;
  ///if true, need distance tables
  bool NeedDistanceTables;
  ///if true, do not clean up
  bool IsCloned;
  ///number of Single-particle orbitals
  IndexType OrbitalSetSize;
  ///number of Single-particle orbitals
  IndexType BasisSetSize;
  ///index of the particle
  IndexType ActivePtcl;
  /** pointer matrix containing the coefficients
   *
   * makeClone makes a shallow copy
   */
  ValueMatrix_t* C;
  ///occupation number
  Vector<RealType> Occ;
  /// Optimizable variables
  opt_variables_type myVars;
  ///name of the basis set
  std::string className;
  /** name of the object
   *
   * Several user classes can own SPOSetBase and use objectName as counter
   */
  std::string objectName;
  
  ///Pass Communicator
  Communicate *myComm;
  
  /** constructor */
  SPOSetBase();

  /** destructor */
  virtual ~SPOSetBase();

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

  inline int getBasisSetSize() const
  {
    return BasisSetSize;
  }

  bool setIdentity(bool useIdentity);

  void checkObject();

  ///get C and Occ
  bool put(xmlNodePtr cur);

  virtual bool put(xmlNodePtr cur, SPOPool_t &spo_pool)
  {
    return put(cur);
  }

  ///reset
  virtual void resetParameters(const opt_variables_type& optVariables)=0;

  virtual void checkInVariables(opt_variables_type& active) {}
  virtual void checkOutVariables(const opt_variables_type& active) {}

  // Evaluate the derivative of the optimized orbitals with
  // respect to the parameters
  virtual void evaluateDerivatives
  (ParticleSet& P, int iat, const opt_variables_type& active,
   ValueMatrix_t& d_phi, ValueMatrix_t& d_lapl_phi) {}


  ///reset the target particleset
  virtual void resetTargetParticleSet(ParticleSet& P)=0;
  /** set the OrbitalSetSize
   * @param norbs number of single-particle orbitals
   */
  virtual void setOrbitalSetSize(int norbs)=0;

  virtual void
  evaluate (const ParticleSet& P, PosType &r, ValueVector_t &psi)
  {
    app_error() << "Need specialization for SPOSetBase::evaluate "
                << "(const ParticleSet& P, PosType &r).\n";
    abort();
  }

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)=0;

  /** compute dot_product of new row and old row */
  virtual ValueType RATIO(const ParticleSet& P, int iat, const ValueType*
      restrict arow);

  /** evaluate VGL using SoA container for gl
   *
   * If newp is true, use particle set data for the proposed move
   */
  virtual void
    evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl, bool newp);

  /** evaluate values for the virtual moves, e.g., sphere move for nonlocalPP
   * @param VP virtual particle set
   * @param psiM single-particle orbitals psiM(i,j) for the i-th particle and the j-th orbital
   */
  virtual void
  evaluateValues(const ParticleSet& VP, ValueMatrix_t& psiM);

  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)=0;

  /** evaluate the values, gradients and hessians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)=0;

  /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last)particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param d2logdet laplacians
   *
   * Call evaluate_notranspose to build logdet
   */
#if 0
  virtual void
  evaluate(const ParticleSet& P, int first, int last, ValueMatrix_t &t_logpsi
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

  virtual void
  evaluate(const ParticleSet& P, int first, int last, ValueMatrix_t &t_logpsi
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  virtual void
  evaluate(const ParticleSet& P, int first, int last, ValueMatrix_t &t_logpsi
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);
#endif

  virtual void
  evaluateThirdDeriv(const ParticleSet& P, int first, int last
                     , GGGMatrix_t& grad_grad_grad_logdet);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  returns whether this is an LCOrbitalSetOpt object
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual bool is_of_type_LCOrbitalSetOpt() const { return false; }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)=0;

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src, GradMatrix_t &gradphi);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src
                                   , GradMatrix_t &grad_phi, HessMatrix_t &grad_grad_phi, GradMatrix_t &grad_lapl_phi);

  virtual void evaluateBasis (const ParticleSet &P, int first, int last
                              , ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl);

  virtual void evaluateForDeriv (const ParticleSet &P, int first, int last
                                 , ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl);

  virtual inline void setpm(int x) {};

  virtual void copyParamsFromMatrix (const opt_variables_type& active
                                     , const ValueMatrix_t &mat, std::vector<RealType> &destVec);

  virtual PosType get_k(int orb)
  {
    return PosType();
  }

  /** make a clone of itself
   */
  virtual SPOSetBase* makeClone() const;

  virtual bool transformSPOSet()
  {
    return true;
  }

  // Routine to set up data for the LCOrbitalSetOpt child class specifically
  // Should be left empty for other derived classes
  virtual void init_LCOrbitalSetOpt(const double mix_factor=0.0) { };

  // Routine to update internal data for the LCOrbitalSetOpt child class specifically
  // Should be left empty for other derived classes
  virtual void rotate_B(const std::vector<RealType> &rot_mat) { };

#ifdef QMC_CUDA

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param r is the position of the particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate (const ParticleSet& P, const PosType& r, std::vector<RealType> &psi);

  virtual void initGPU() {  }

  //////////////////////////////////////////
  // Walker-parallel vectorized functions //
  //////////////////////////////////////////
  virtual void
  reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool) { }

  virtual void
  evaluate (std::vector<Walker_t*> &walkers, int iat, gpu::device_vector<CudaValueType*> &phi);

  virtual void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &new_pos
                         , gpu::device_vector<CudaValueType*> &phi);

  virtual void
  evaluate (std::vector<Walker_t*> &walkers,
            std::vector<PosType> &new_pos,
            gpu::device_vector<CudaValueType*> &phi,
            gpu::device_vector<CudaValueType*> &grad_lapl_list,
            int row_stride);

  virtual void
  evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi);
  virtual void
  evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi);
#endif


protected:
  bool putOccupation(xmlNodePtr occ_ptr);
  bool putFromXML(xmlNodePtr coeff_ptr);
  bool putFromH5(const char* fname, xmlNodePtr coeff_ptr);
};

#if defined(ENABLE_SMARTPOINTER)
typedef boost::shared_ptr<SPOSetBase> SPOSetBasePtr;
#else
typedef SPOSetBase*                   SPOSetBasePtr;
#endif

}
#endif


