//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file BasisSetBase.h
 * @brief Declaration of a base class of BasisSet
 */
#ifndef QMCPLUSPLUS_BASISSETBASE_H
#define QMCPLUSPLUS_BASISSETBASE_H

#include "Particle/ParticleSet.h"
#include "Message/MPIObjectBase.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/SPOSetInfo.h"
#include <QMCWaveFunctions/SPOSetInputInfo.h>
#include "QMCWaveFunctions/SPOSetBase.h"


namespace qmcplusplus
{

/** base class for a basis set
 *
 * Define a common storage for the derived classes and
 * provides  a minimal set of interfaces to get/set BasisSetSize.
 */
template<typename T>
struct BasisSetBase: public OrbitalSetTraits<T>
{

  enum {MAXINDEX=2+OHMMS_DIM};
  typedef typename OrbitalSetTraits<T>::RealType      RealType;
  typedef typename OrbitalSetTraits<T>::ValueType     ValueType;
  typedef typename OrbitalSetTraits<T>::IndexType     IndexType;
  typedef typename OrbitalSetTraits<T>::HessType      HessType;
  typedef typename OrbitalSetTraits<T>::IndexVector_t IndexVector_t;
  typedef typename OrbitalSetTraits<T>::ValueVector_t ValueVector_t;
  typedef typename OrbitalSetTraits<T>::ValueMatrix_t ValueMatrix_t;
  typedef typename OrbitalSetTraits<T>::GradVector_t  GradVector_t;
  typedef typename OrbitalSetTraits<T>::GradMatrix_t  GradMatrix_t;
  typedef typename OrbitalSetTraits<T>::HessVector_t  HessVector_t;
  typedef typename OrbitalSetTraits<T>::HessMatrix_t  HessMatrix_t;
  typedef TinyVector<HessType, OHMMS_DIM>                     GGGType;
  typedef Vector<GGGType>                             GGGVector_t;
  typedef Matrix<GGGType>                             GGGMatrix_t;


  ///size of the basis set
  IndexType BasisSetSize;
  ///index of the particle
  IndexType ActivePtcl;
  ///counter to keep track
  unsigned long Counter;
  ///phi[i] the value of the i-th basis set
  ValueVector_t Phi;
  ///dphi[i] the gradient of the i-th basis set
  GradVector_t  dPhi;
  ///d2phi[i] the laplacian of the i-th basis set
  ValueVector_t d2Phi;
  ///grad_grad_Phi[i] the full hessian of the i-th basis set
  HessVector_t  grad_grad_Phi;
  ///grad_grad_grad_Phi the full hessian of the i-th basis set
  GGGVector_t  grad_grad_grad_Phi;
  ///container to store value, laplacian and gradient
  ValueMatrix_t Temp;

  ValueMatrix_t Y;
  GradMatrix_t dY;
  ValueMatrix_t d2Y;

  ///default constructor
  BasisSetBase():BasisSetSize(0), ActivePtcl(-1), Counter(0) { }
  ///virtual destructor
  virtual ~BasisSetBase() { }
  /** resize the container */
  void resize(int ntargets)
  {
    if(BasisSetSize)
    {
      Phi.resize(BasisSetSize);
      dPhi.resize(BasisSetSize);
      d2Phi.resize(BasisSetSize);
      grad_grad_Phi.resize(BasisSetSize);
      grad_grad_grad_Phi.resize(BasisSetSize);
      Temp.resize(BasisSetSize,MAXINDEX);
      Y.resize(ntargets,BasisSetSize);
      dY.resize(ntargets,BasisSetSize);
      d2Y.resize(ntargets,BasisSetSize);
    }
    else
    {
      app_error() << "  BasisSetBase::BasisSetSize == 0" << std::endl;
    }
  }

  ///clone the basis set
  virtual BasisSetBase* makeClone() const=0;
  /** return the basis set size */
  inline IndexType getBasisSetSize() const
  {
    return BasisSetSize;
  }

  /**@{ functions to perform optimizations  */
  /** checkIn optimizable variables */
  virtual void checkInVariables(opt_variables_type& active)
  { }
  /** checkOut optimizable variables */
  virtual void checkOutVariables(const opt_variables_type& active)
  { }
  /** reset parameters */
  virtual void resetParameters(const opt_variables_type& active)
  {}
  /**@}*/
  ///resize the basis set
  virtual void setBasisSetSize(int nbs) = 0;
  ///reset the target particle set
  virtual void resetTargetParticleSet(ParticleSet& P)=0;

  virtual void evaluateWithHessian(const ParticleSet& P, int iat)=0;
  virtual void evaluateWithThirdDeriv(const ParticleSet& P, int iat)=0;
  virtual void evaluateThirdDerivOnly(const ParticleSet& P, int iat)=0;
  virtual void evaluateForWalkerMove(const ParticleSet& P)=0;
  virtual void evaluateForWalkerMove(const ParticleSet& P, int iat) =0;
  virtual void evaluateForPtclMove(const ParticleSet& P, int iat) =0;
  virtual void evaluateAllForPtclMove(const ParticleSet& P, int iat) =0;
  virtual void evaluateForPtclMoveWithHessian(const ParticleSet& P, int iat)=0;
};

/** Base for real basis set
 *
 * Equivalent to BasisSetBase with minimum requirements
 * Used by lcao
 */
template<typename T>
struct RealBasisSetBase
{
  typedef T value_type;
  typedef VectorSoaContainer<T,OHMMS_DIM+2> vgl_type;
  ///size of the basis set
  int BasisSetSize;

  inline int getBasisSetSize()
  {
    return BasisSetSize;
  }

  virtual RealBasisSetBase<T>* makeClone() const = 0;
  virtual void setBasisSetSize(int nbs)=0;
  virtual void evaluateVGL(const ParticleSet& P, int iat, vgl_type& vgl)=0;
  virtual void evaluateV(const ParticleSet& P, int iat, value_type* restrict vals)=0;
};


/** base class for the real BasisSet builder
 *
 * \warning {
 * We have not quite figured out how to use real/complex efficiently.
 * There are three cases we have to deal with
 * - real basis functions and real coefficients
 * - real basis functions and complex coefficients
 * - complex basis functions and complex coefficients
 * For now, we decide to keep both real and complex basis sets and expect
 * the user classes {\bf KNOW} what they need to use.
 * }
 */
struct BasisSetBuilder: public QMCTraits, public MPIObjectBase
{
  typedef std::map<std::string,SPOSetBase*> SPOPool_t;
  typedef std::vector<int> indices_t;
  typedef std::vector<RealType> energies_t;


  /// whether implementation conforms only to legacy standard
  bool legacy;

  /// state info of all possible states available in the basis
  std::vector<SPOSetInfo*> states;

  /// list of all sposets created by this builder
  std::vector<SPOSetBase*> sposets;

  BasisSetBuilder();
  virtual ~BasisSetBuilder() {}
  virtual bool put(xmlNodePtr cur)=0;

  /// reserve space for states (usually only one set, multiple for e.g. spin dependent einspline)
  void reserve_states(int nsets=1);

  /// allow modification of state information
  inline void modify_states(int index=0)
  {
    states[index]->modify();
  }

  /// clear state information
  inline void clear_states(int index=0)
  {
    states[index]->clear();
  }

  /// create an sposet from xml (legacy)
  virtual SPOSetBase* createSPOSetFromXML(xmlNodePtr cur)=0;

  /// create an sposet from a general xml request
  virtual SPOSetBase* createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input_info);

  /// create an sposet from xml and save the resulting SPOSet
  SPOSetBase* createSPOSet(xmlNodePtr cur);


};

}
#endif
