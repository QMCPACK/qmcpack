//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file PWRealOrbitalSetT.h
 * @brief Define PWRealOrbitalSetT derived from SPOSetT
 *
 * This is a specialized single-particle orbital set for real trial
 * wavefunctions and enabled with QMC_COMPLEX=0
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_REALORBITALSETT_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_REALORBITALSETT_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSetT.h"
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
template<class T>
class PWRealOrbitalSetT : public SPOSetT<T>
{
public:
  using BasisSet_t = PWBasis;
  using PWBasisPtr = PWBasis*;

  using IndexType   = typename SPOSetT<T>::IndexType;
  using RealType    = typename SPOSetT<T>::RealType;
  using ComplexType = typename SPOSetT<T>::ComplexType;
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using GradMatrix  = typename SPOSetT<T>::GradMatrix;
  using HessMatrix  = typename SPOSetT<T>::HessMatrix;
  using PosType     = typename SPOSetT<T>::PosType;

  /** inherit the enum of BasisSet_t */
  enum
  {
    PW_VALUE    = BasisSet_t::PW_VALUE,
    PW_LAP      = BasisSet_t::PW_LAP,
    PW_GRADX    = BasisSet_t::PW_GRADX,
    PW_GRADY    = BasisSet_t::PW_GRADY,
    PW_GRADZ    = BasisSet_t::PW_GRADZ,
    PW_MAXINDEX = BasisSet_t::PW_MAXINDEX
  };

  /** default constructor
  */
  PWRealOrbitalSetT(const std::string& my_name)
      : SPOSetT<T>(my_name), OwnBasisSet(false), myBasisSet(nullptr), BasisSetSize(0)
  {}

  std::string getClassName() const override { return "PWRealOrbitalSetT"; }

  /** delete BasisSet only it owns this
   *
   * Builder takes care of who owns what
   */
  ~PWRealOrbitalSetT() override;

  std::unique_ptr<SPOSetT<T>> makeClone() const override;

  /** resize  the orbital base
   * @param bset PWBasis
   * @param nbands number of bands
   * @param cleaup if true, owns PWBasis. Will clean up.
   */
  void resize(PWBasisPtr bset, int nbands, bool cleanup = false);

  /** add eigenstate for jorb-th orbital
   * @param coefs real input data
   * @param jorb orbital index
   */
  void addVector(const std::vector<RealType>& coefs, int jorb);

  /** add eigenstate for jorb-th orbital
   * @param coefs complex input data
   * @param jorb orbital index
   */
  void addVector(const std::vector<ComplexType>& coefs, int jorb);

  void setOrbitalSetSize(int norbs) override;

  inline T evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return real(BLAS::dot(BasisSetSize, CC[ib], myBasisSet->Zv.data()));
  }

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

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
                            HessMatrix& grad_grad_logdet) override
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
  ///number of basis
  IndexType BasisSetSize;
  ///Plane-wave coefficients of complex: (iband,g-vector)
  Matrix<ComplexType> CC;
  /// temporary array to perform gemm operation
  Matrix<ComplexType> Temp;
  ///temporary complex vector before assigning to a real psi
  Vector<ComplexType> tempPsi;
};
} // namespace qmcplusplus
#endif
