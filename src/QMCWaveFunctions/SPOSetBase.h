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

#include <stdexcept>
#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

namespace qmcplusplus
{

class OrbitalBase;

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
  typedef ParticleSet::Walker_t                      Walker_t;
  typedef std::map<std::string,SPOSetBase*> SPOPool_t;

  ///index in the builder list of sposets
  int builder_index;
  ///true if C is an identity matrix
  bool Identity;
  ///true if SPO is optimizable
  bool Optimizable;
  ///true if precomputed distance tables are needed
  bool NeedsDistanceTable;
  ///flag to calculate ionic derivatives
  bool ionDerivs;
  ///total number of orbitals
  IndexType TotalOrbitalSize;
  ///number of Single-particle orbitals
  IndexType OrbitalSetSize;
  ///number of Single-particle orbitals
  IndexType BasisSetSize;
  ///index of the particle
  IndexType ActivePtcl;
  ///matrix to store temporary value before transpose
  ValueMatrix_t t_logpsi;
  ///matrix containing the coefficients
  ValueMatrix_t C;
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

  /** constructor
   */
  SPOSetBase()
    :Identity(false),TotalOrbitalSize(0),OrbitalSetSize(0),BasisSetSize(0),
    NeedsDistanceTable(false),
    ActivePtcl(-1),Optimizable(false),ionDerivs(false),builder_index(-1)
  {
    className="invalid";
  }

  /** destructor
   */
  virtual ~SPOSetBase() {}

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


  bool setIdentity(bool useIdentity)
  {
    Identity = useIdentity;

    if ( (OrbitalSetSize > 0) && (BasisSetSize > 0) )
      C.resize(OrbitalSetSize,BasisSetSize);
    else {
      app_error() << "either OrbitalSetSize or BasisSetSize has an invalid value !!\n";
      app_error() << "OrbitalSetSize = " << OrbitalSetSize << std::endl;
      app_error() << "BasisSetSize = " << BasisSetSize << std::endl;
      abort();
    }

    if (OrbitalSetSize <= BasisSetSize) {
      for (int i=0; i<OrbitalSetSize; i++)
        C(i,i) = 1.0;
    }
    else {
      for (int i=0; i<BasisSetSize; i++)
        C(i,i) = 1.0;
    }

    return true;
  }

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
  virtual void
  evaluate(const ParticleSet& P, int first, int last
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

  virtual void
  evaluate(const ParticleSet& P, int first, int last
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  virtual void
  evaluate(const ParticleSet& P, int first, int last
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

  virtual void
  evaluateThirdDeriv(const ParticleSet& P, int first, int last
                     , GGGMatrix_t& grad_grad_grad_logdet);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  returns whether this is an LCOrbitalSetOpt object
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual bool is_of_type_LCOrbitalSetOpt() const { return false; }

  virtual OrbitalBase * tf_component();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Evaluates the values, x,y,z derivatives, and x-y-z-summed second derivatives of the
  ///         specified linear combinations of basis set orbitals at the specified particles'
  ///         positions.
  ///
  /// \param[in]      P            Object containing information on particle positions.
  /// \param[in]      mt           the move type: 'p' for particle move, 'w' for walker move
  /// \param[in]      ostart       Iterator for the start of the index range specifying which linear combinations of orbitals to evaluate.
  /// \param[in]      oend         Iterator for the end   of the index range specifying which linear combinations of orbitals to evaluate.
  /// \param[in]      pstart       Iterator for the start of the index range specifying which particles' positions to use.
  /// \param[in]      pend         Iterator for the end   of the index range specifying which particles' positions to use.
  /// \param[in,out]  vmat         On input, points to an array of length (# of linear combinations) * (# of particle).
  ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix
  ///                              of the values of the specified linear combinations for the specified particles' positions.
  /// \param[in,out]  gmat         On input, points to an array of length (# of linear combinations) * (# of particle).
  ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix,
  ///                              each element of which is a length 3 vector containing the x,y,z gradients of the values in vmat.
  /// \param[in,out]  lmat         On input, points to an array of length (# of linear combinations) * (# of particle).
  ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix,
  ///                              each element of which is the sum of x^2, y^2, and z^2 second derivatives of the values in vmat.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class IntIter> void evaluate_notranspose_iterators(const ParticleSet& P,
                                                              const char mt,
                                                              IntIter ostart,
                                                              IntIter oend,
                                                              IntIter pstart,
                                                              IntIter pend,
                                                              ValueType * const vmat,
                                                              GradType * const gmat,
                                                              ValueType * const lmat) {

    // prepare the vector of orbital indices
    std::vector<int>::const_iterator oe = SPOSetBase::prepare_index_vector(ostart, oend, this->m_oidx);

    // prepare the vector of particle indices
    std::vector<int>::const_iterator pe = SPOSetBase::prepare_index_vector(pstart, pend, this->m_pidx);

    // call the child class's evaluate function
    this->evaluate_notranspose_general(P, mt, m_oidx.begin(), oe, m_pidx.begin(), pe, vmat, gmat, lmat);

  }

  virtual void evaluate_notranspose_general(const ParticleSet& P,
                                            const char mt,
                                            std::vector<int>::const_iterator ostart,
                                            std::vector<int>::const_iterator oend,
                                            std::vector<int>::const_iterator pstart,
                                            std::vector<int>::const_iterator pend,
                                            ValueType * const vmat,
                                            GradType * const gmat,
                                            ValueType * const lmat) {

    throw std::runtime_error("SPOSetBase::evaluate_notranspose_general not implemented");

  }

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

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Enlarges the supplied vector if it is not big enough
  ///
  /// \param[in,out]  v              the vector
  /// \param[in]      n              the minimum length we want the vector to have
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class T> static void ensure_vector_is_big_enough(T & v, const size_t n) {
    if ( v.size() < n )
      v.resize(n);
  }


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

  virtual void set_orbital_mixing_factor(const double factor);

  /// \brief  vector to hold orbital indices
  std::vector<int> m_oidx;

  /// \brief  vector to hold particle indices
  std::vector<int> m_pidx;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Place the indices from a range of indices specified by iterators into a vector.
  ///         Note that the vector may be longer than the index range (it is not shrunk to fit it)
  ///         but that an iterator to the end of the range is returned.
  ///
  /// \param[in]      start       iterator for the start of the range (should dereference to int)
  /// \param[in]      end         iterator for the end   of the range (should dereference to int)
  /// \param[in]      vec         vector to store the range in
  ///
  /// \return  iterator to the end of the entered range
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class IntIter> static std::vector<int>::iterator prepare_index_vector(IntIter start, IntIter end, std::vector<int> & vec) {

    // get the length
    int length = 0;
    for (IntIter s = start; s != end; s++)
      length++;

    // expand the vector if necessary
    SPOSetBase::ensure_vector_is_big_enough(vec, length);

    // put the values in the vector
    std::copy(start, end, vec.begin());

    // return an iterator to the end of the range inside the vector
    return ( vec.begin() + length );

  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Place a range of contiguous indices [start,end) into a vector.
  ///         Note that the vector may be longer than the index range (it is not shrunk to fit it)
  ///         but an iterator to the end of the range returned.
  ///
  /// \param[in]      start       start of the range
  /// \param[in]      end         end of the range
  /// \param[in]      vec         vector to store the range in
  ///
  /// \return  iterator to the end of the entered range
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  static std::vector<int>::iterator prepare_index_vector_contiguous(const int start, const int end, std::vector<int> & vec) {

    // check sanity
    if ( end < start )
      throw std::runtime_error("end is less than start in SPOSetBase::prepare_index_vector_contiguous");

    // expand the vector if necessary
    SPOSetBase::ensure_vector_is_big_enough(vec, end - start);

    // put the range into the vector
    std::vector<int>::iterator it = vec.begin();
    for(int i = start; i < end; i++, it++)
      *it = i;

    // return the end of the range
    return it;

  }

};

#if defined(ENABLE_SMARTPOINTER)
typedef boost::shared_ptr<SPOSetBase> SPOSetBasePtr;
#else
typedef SPOSetBase*                   SPOSetBasePtr;
#endif

}
#endif


