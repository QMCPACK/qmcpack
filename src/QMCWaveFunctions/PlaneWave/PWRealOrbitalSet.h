//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file PWRealOrbitalSet.h
 * @brief Define PWRealOrbitalSet derived from SPOSetBase
 *
 * This is a specialized single-particle orbital set for real trial
 * wavefunctions and enabled with QMC_COMPLEX=0
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_REALORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_REALORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{

class PWRealOrbitalSet: public SPOSetBase
{

public:

  typedef PWBasis                    BasisSet_t;
#if defined(ENABLE_SMARTPOINTER)
  typedef boost::shared_ptr<PWBasis> PWBasisPtr;
#else
  typedef PWBasis*                   PWBasisPtr;
#endif

  /** inherit the enum of BasisSet_t */
  enum {PW_VALUE=   BasisSet_t::PW_VALUE,
        PW_LAP=     BasisSet_t::PW_LAP,
        PW_GRADX=   BasisSet_t::PW_GRADX,
        PW_GRADY=   BasisSet_t::PW_GRADY,
        PW_GRADZ=   BasisSet_t::PW_GRADZ,
        PW_MAXINDEX=BasisSet_t::PW_MAXINDEX
       };

  /** default constructor
  */
  PWRealOrbitalSet(): OwnBasisSet(false)
  {
  }

  /** delete BasisSet only it owns this
   *
   * Builder takes care of who owns what
   */
  ~PWRealOrbitalSet();

  /** resize  the orbital base
   * @param bset PWBasis
   * @param nbands number of bands
   * @param cleaup if true, owns PWBasis. Will clean up.
   */
  void resize(PWBasisPtr bset, int nbands, bool cleanup=false);

  /** add eigenstate for jorb-th orbital
   * @param coefs real input data
   * @param jorb orbital index
   */
  void addVector(const std::vector<RealType>& coefs,int jorb);

  /** add eigenstate for jorb-th orbital
   * @param coefs complex input data
   * @param jorb orbital index
   */
  void addVector(const std::vector<ComplexType>& coefs,int jorb);

  void resetParameters(const opt_variables_type& optVariables);

  void setOrbitalSetSize(int norbs);

  void resetTargetParticleSet(ParticleSet& P);

  inline ValueType evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return real(BLAS::dot(BasisSetSize,CC[ib],myBasisSet->Zv.data()));
  }

  void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }


  /** boolean
   *
   * If true, this has to delete the BasisSet
   */
  bool OwnBasisSet;
  ///TwistAngle of this PWRealOrbitalSet
  PosType TwistAngle;
  ///My basis set
  PWBasisPtr myBasisSet;
  ///Plane-wave coefficients of complex: (iband,g-vector)
  Matrix<ComplexType> CC;
  /// temporary array to perform gemm operation
  Matrix<ComplexType> Temp;
  ///temporary complex vector before assigning to a real psi
  Vector<ComplexType> tempPsi;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1769 $   $Date: 2007-02-17 17:38:34 -0600 (Sat, 17 Feb 2007) $
 * $Id: PWRealOrbitalSet.h 1769 2007-02-17 23:38:34Z jnkim $
 ***************************************************************************/
