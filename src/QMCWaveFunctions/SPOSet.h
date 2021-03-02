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
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#ifdef QMC_CUDA
#include "type_traits/CUDATypes.h"
#endif

namespace qmcplusplus
{

class ResourceCollection;

/** base class for Single-particle orbital sets
 *
 * SPOSet stands for S(ingle)P(article)O(rbital)Set which contains
 * a number of single-particle orbitals with capabilities of evaluating \f$ \psi_j({\bf r}_i)\f$
 */
class SPOSet : public QMCTraits
{
public:
  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t GradVector_t;
  typedef OrbitalSetTraits<ValueType>::GradMatrix_t GradMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessVector_t HessVector_t;
  typedef OrbitalSetTraits<ValueType>::HessMatrix_t HessMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType HessType;
  typedef Array<HessType, OHMMS_DIM> HessArray_t;
  typedef OrbitalSetTraits<ValueType>::GradHessType GGGType;
  typedef OrbitalSetTraits<ValueType>::GradHessVector_t GGGVector_t;
  typedef OrbitalSetTraits<ValueType>::GradHessMatrix_t GGGMatrix_t;
  typedef OrbitalSetTraits<ValueType>::VGLVector_t VGLVector_t;
  typedef ParticleSet::Walker_t Walker_t;
  typedef std::map<std::string, SPOSet*> SPOPool_t;

  /** constructor */
  SPOSet(bool use_OMP_offload = false, bool ion_deriv = false, bool optimizable = false);

  /** destructor
   *
   * Derived class destructor needs to pay extra attention to freeing memory shared among clones of SPOSet.
   */
  virtual ~SPOSet() {}

  // accessor function to Optimizable
  inline bool isOptimizable() const { return Optimizable; }

  /** return the size of the orbital set
   * Ye: this needs to be replaced by getOrbitalSetSize();
   */
  inline int size() const { return OrbitalSetSize; }

  inline const std::string& getClassName() const { return className; }

  /** print basic SPOSet information
   */
  void basic_report(const std::string& pad = "") const;

  /** print SPOSet information
   */
  virtual void report(const std::string& pad = "") const { basic_report(pad); }


  /** return the size of the orbitals
   */
  inline int getOrbitalSetSize() const { return OrbitalSetSize; }

  /** Query if this SPOSet uses OpenMP offload
  */
  inline bool isOMPoffload() const { return useOMPoffload; }

  /** Query if this SPOSet has an explicit ion dependence. returns true if it does.
  */
  inline bool hasIonDerivs() const { return ionDerivs; }

  /// return the size of the basis set if there is any
  virtual int getBasisSetSize() const { return 0; }

  /// check a few key parameters before putting the SPO into a determinant
  virtual void checkObject() const {}

  /// create optimizable orbital rotation parameters
  // Single Slater creation
  virtual void buildOptVariables(const size_t nel) {}
  // For the MSD case rotations must be created in MultiSlaterFast class
  virtual void buildOptVariables(const std::vector<std::pair<int, int>>& rotations) {}
  // store parameters before getting destroyed by rotation.
  virtual void storeParamsBeforeRotation() {}
  // apply rotation to all the orbitals
  virtual void applyRotation(const ValueMatrix_t& rot_mat, bool use_stored_copy = false)
  {
    std::ostringstream o;
    o << "SPOSet::applyRotation is not implemented by " << className << std::endl;
    APP_ABORT(o.str());
  }
  /// reset parameters to the values from optimizer
  virtual void resetParameters(const opt_variables_type& optVariables) = 0;

  /// check in/out parameters to the global list of parameters used by the optimizer
  virtual void checkInVariables(opt_variables_type& active) {}
  virtual void checkOutVariables(const opt_variables_type& active) {}

  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   std::vector<ValueType>& dlogpsi,
                                   std::vector<ValueType>& dhpsioverpsi,
                                   const int& FirstIndex,
                                   const int& LastIndex)
  {}
  /** Evaluate the derivative of the optimized orbitals with respect to the parameters
   *  this is used only for MSD, to be refined for better serving both single and multi SD
   */
  virtual void evaluateDerivatives(ParticleSet& P,
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
                                   const std::vector<std::vector<int>>& lookup_tbl)
  {}

  /** Evaluate the derivative of the optimized orbitals with respect to the parameters
   *  this is used only for MSD, to be refined for better serving both single and multi SD
   */
  virtual void evaluateDerivativesWF(ParticleSet& P,
                                     const opt_variables_type& optvars,
                                     std::vector<ValueType>& dlogpsi,
                                     const QTFull::ValueType& psiCurrent,
                                     const std::vector<ValueType>& Coeff,
                                     const std::vector<size_t>& C2node_up,
                                     const std::vector<size_t>& C2node_dn,
                                     const ValueVector_t& detValues_up,
                                     const ValueVector_t& detValues_dn,
                                     const ValueMatrix_t& M_up,
                                     const ValueMatrix_t& M_dn,
                                     const ValueMatrix_t& Minv_up,
                                     const ValueMatrix_t& Minv_dn,
                                     const std::vector<int>& detData_up,
                                     const std::vector<std::vector<int>>& lookup_tbl)
  {}

  /** set the OrbitalSetSize
   * @param norbs number of single-particle orbitals
   * Ye: I prefer to remove this interface in the future. SPOSet builders need to handle the size correctly.
   * It doesn't make sense allowing to set the value at any place in the code.
   */
  virtual void setOrbitalSetSize(int norbs) = 0;

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi) = 0;

  /** evaluate determinant ratios for virtual moves, e.g., sphere move for nonlocalPP
   * @param VP virtual particle set
   * @param psi values of the SPO, used as a scratch space if needed
   * @param psiinv the row of inverse slater matrix corresponding to the particle moved virtually
   * @param ratios return determinant ratios
   */
  virtual void evaluateDetRatios(const VirtualParticleSet& VP,
                                 ValueVector_t& psi,
                                 const ValueVector_t& psiinv,
                                 std::vector<ValueType>& ratios);

  /** evaluate determinant ratios for virtual moves, e.g., sphere move for nonlocalPP, of multiple walkers
   * @param spo_list the list of SPOSet pointers in a walker batch
   * @param vp_list a list of virtual particle sets in a walker batch
   * @param psi_list a list of values of the SPO, used as a scratch space if needed
   * @param invRow_ptr_list a list of pointers to the rows of inverse slater matrix corresponding to the particles moved virtually
   * @param ratios_list a list of returning determinant ratios
   */
  virtual void mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                    const RefVector<const VirtualParticleSet>& vp_list,
                                    const RefVector<ValueVector_t>& psi_list,
                                    const std::vector<const ValueType*>& invRow_ptr_list,
                                    std::vector<std::vector<ValueType>>& ratios_list) const;

  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   * @param dpsi gradients of the SPO
   * @param d2psi laplacians of the SPO
   */
  virtual void evaluateVGL(const ParticleSet& P,
                           int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           ValueVector_t& d2psi) = 0;

  /** evaluate the values, gradients and laplacians of this single-particle orbital sets of multiple walkers
   * @param spo_list the list of SPOSet pointers in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat active particle
   * @param psi_v_list the list of value vector pointers in a walker batch
   * @param dpsi_v_list the list of gradient vector pointers in a walker batch
   * @param d2psi_v_list the list of laplacian vector pointers in a walker batch
   */
  virtual void mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& spo_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector_t>& psi_v_list,
                              const RefVector<GradVector_t>& dpsi_v_list,
                              const RefVector<ValueVector_t>& d2psi_v_list) const;

  /** evaluate the values, gradients and laplacians of this single-particle orbital sets
   *  and determinant ratio and grads of multiple walkers
   * @param spo_list the list of SPOSet pointers in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat active particle
   * @param phi_vgl_v orbital values, gradients and laplacians of all the walkers
   * @param psi_ratio_grads_v determinant ratio and grads of all the walkers
   */
  virtual void mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                              const RefVectorWithLeader<ParticleSet>& P_list,
                                              int iat,
                                              const std::vector<const ValueType*>& invRow_ptr_list,
                                              VGLVector_t& phi_vgl_v,
                                              std::vector<ValueType>& ratios,
                                              std::vector<GradType>& grads) const;

  /** evaluate the values, gradients and hessians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   * @param dpsi gradients of the SPO
   * @param grad_grad_psi hessians of the SPO
   */
  virtual void evaluateVGH(const ParticleSet& P,
                           int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           HessVector_t& grad_grad_psi);

  /** evaluate the values, gradients, hessians, and grad hessians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   * @param dpsi gradients of the SPO
   * @param grad_grad_psi hessians of the SPO
   * @param grad_grad_grad_psi grad hessians of the SPO
   */
  virtual void evaluateVGHGH(const ParticleSet& P,
                             int iat,
                             ValueVector_t& psi,
                             GradVector_t& dpsi,
                             HessVector_t& grad_grad_psi,
                             GGGVector_t& grad_grad_grad_psi);

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void evaluate_spin(const ParticleSet& P, int iat, ValueVector_t& psi, ValueVector_t& dpsi);

  /** evaluate the third derivatives of this single-particle orbital set
   * @param P current ParticleSet
   * @param first first particle
   * @param last last particle
   * @param grad_grad_grad_logdet third derivatives of the SPO
   */
  virtual void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet);

  /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last) particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param d2logdet laplacians
   *
   */
  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    ValueMatrix_t& d2logdet) = 0;

  virtual void mw_evaluate_notranspose(const RefVectorWithLeader<SPOSet>& spo_list,
                                       const RefVectorWithLeader<ParticleSet>& P_list,
                                       int first,
                                       int last,
                                       const RefVector<ValueMatrix_t>& logdet_list,
                                       const RefVector<GradMatrix_t>& dlogdet_list,
                                       const RefVector<ValueMatrix_t>& d2logdet_list) const;

  /** evaluate the values, gradients and hessians of this single-particle orbital for [first,last) particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param grad_grad_logdet hessians
   *
   */
  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    HessMatrix_t& grad_grad_logdet);

  /** evaluate the values, gradients, hessians and third derivatives of this single-particle orbital for [first,last) particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param grad_grad_logdet hessians
   * @param grad_grad_grad_logdet third derivatives
   *
   */
  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    HessMatrix_t& grad_grad_logdet,
                                    GGGMatrix_t& grad_grad_grad_logdet);

  /** evaluate the gradients of this single-particle orbital
   *  for [first,last) target particles with respect to the given source particle
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param iat_src source particle index
   * @param gradphi gradients
   *
   */
  virtual void evaluateGradSource(const ParticleSet& P,
                                  int first,
                                  int last,
                                  const ParticleSet& source,
                                  int iat_src,
                                  GradMatrix_t& gradphi);

  /** evaluate the gradients of values, gradients, laplacians of this single-particle orbital
   *  for [first,last) target particles with respect to the given source particle
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param iat_src source particle index
   * @param gradphi gradients of values
   * @param grad_grad_phi gradients of gradients
   * @param grad_lapl_phi gradients of laplacians
   *
   */
  virtual void evaluateGradSource(const ParticleSet& P,
                                  int first,
                                  int last,
                                  const ParticleSet& source,
                                  int iat_src,
                                  GradMatrix_t& grad_phi,
                                  HessMatrix_t& grad_grad_phi,
                                  GradMatrix_t& grad_lapl_phi);

  /** access the k point related to the given orbital */
  virtual PosType get_k(int orb) { return PosType(); }

  /** initialize a shared resource and hand it to collection
   */
  virtual void createResource(ResourceCollection& collection) const {}

  /** acquire a shared resource from collection
   */
  virtual void acquireResource(ResourceCollection& collection) {}

  /** return a shared resource to collection
   */
  virtual void releaseResource(ResourceCollection& collection) {}

  /** make a clone of itself
   * every derived class must implement this to have threading working correctly.
   */
  virtual SPOSet* makeClone() const;

  /** Used only by cusp correction in AOS LCAO.
   * Ye: the SoA LCAO moves all this responsibility to the builder.
   * This interface should be removed with AoS.
   */
  virtual bool transformSPOSet() { return true; }

  /** finalize the construction of SPOSet
   *
   * for example, classes serving accelerators may need to transfer data from host to device
   * after the host side objects are built.
   */
  virtual void finalizeConstruction() {}

  /// set object name
  void setName(const std::string& name) { myName = name; }
  /// return object name
  const std::string& getName() const { return myName; }

#ifdef QMC_CUDA
  /** Evaluate the SPO value at an explicit position.
   * Ye: This is used only for debugging the CUDA code and should be removed.
   */
  virtual void evaluate(const ParticleSet& P, PosType& r, ValueVector_t& psi);

  using CTS = CUDAGlobalTypes;

  //////////////////////////////////////////
  // Walker-parallel vectorized functions //
  //////////////////////////////////////////
  virtual void reserve(PointerPool<gpu::device_vector<CTS::ValueType>>& pool) {}

  virtual void evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::ValueType*>& phi);

  virtual void evaluate(std::vector<Walker_t*>& walkers,
                        std::vector<PosType>& new_pos,
                        gpu::device_vector<CTS::ValueType*>& phi);

  virtual void evaluate(std::vector<Walker_t*>& walkers,
                        std::vector<PosType>& new_pos,
                        gpu::device_vector<CTS::ValueType*>& phi,
                        gpu::device_vector<CTS::ValueType*>& grad_lapl_list,
                        int row_stride);

  virtual void evaluate(std::vector<Walker_t*>& walkers,
                        std::vector<PosType>& new_pos,
                        gpu::device_vector<CTS::ValueType*>& phi,
                        gpu::device_vector<CTS::ValueType*>& grad_lapl_list,
                        int row_stride,
                        int k,
                        bool klinear);

  virtual void evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi);
  virtual void evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi);
#endif

protected:
  ///true, if the derived class uses OpenMP offload and statisfies a few assumptions
  const bool useOMPoffload;
  ///true, if the derived class has non-zero ionic derivatives.
  const bool ionDerivs;
  ///true if SPO is optimizable
  const bool Optimizable;
  ///number of Single-particle orbitals
  IndexType OrbitalSetSize;
  /// Optimizable variables
  opt_variables_type myVars;
  ///name of the class
  std::string className;
  /// name of the object, unique identifier
  std::string myName;
};

typedef SPOSet* SPOSetPtr;

} // namespace qmcplusplus
#endif
