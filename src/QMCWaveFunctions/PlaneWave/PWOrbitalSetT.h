//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

/** @file PWOrbitalSetT.h
 * @brief Definition of member functions of Plane-wave basis set
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSETT_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSETT_BLAS_H

#include "CPU/BLAS.hpp"
#include "QMCWaveFunctions/PlaneWave/PWBasisT.h"
#include "QMCWaveFunctions/SPOSetT.h"
#include "type_traits/complex_help.hpp"

namespace qmcplusplus
{

template<class T>
class PWOrbitalSetT : public SPOSetT<T>
{
public:
  using RealType    = typename SPOSetT<T>::RealType;
  using ComplexType = T;
  using PosType     = typename SPOSetT<T>::PosType;
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
  using GradMatrix  = typename SPOSetT<T>::GradMatrix;
  using GradType    = typename SPOSetT<T>::GradType;
  using IndexType   = typename SPOSetT<T>::IndexType;

  using BasisSet_t = PWBasisT<T>;
  using PWBasisPtr = PWBasisT<T>*;

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
  PWOrbitalSetT<T>(const std::string& my_name)
      : SPOSetT<T>(my_name), OwnBasisSet(false), myBasisSet(nullptr), BasisSetSize(0), C(nullptr), IsCloned(false)
  {}

  std::string getClassName() const override { return "PWOrbitalSetT"; }

  /** delete BasisSet only it owns this
     *
     * Builder takes care of who owns what
     */
  ~PWOrbitalSetT<T>() override;

  std::unique_ptr<SPOSetT<T>> makeClone() const override;
  /** resize  the orbital base
     * @param bset PWBasis
     * @param nbands number of bands
     * @param cleaup if true, owns PWBasis. Will clean up.
     */
  void resize(PWBasisPtr bset, int nbands, bool cleanup = false);

  /** Builder class takes care of the assertion
     */
  void addVector(const std::vector<ComplexType>& coefs, int jorb);
  void addVector(const std::vector<RealType>& coefs, int jorb);

  void setOrbitalSetSize(int norbs) override;

  inline T evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return BLAS::dot(BasisSetSize, (*C)[ib], myBasisSet->Zv.data());
  }

  void evaluateValue(const ParticleSetT<T>& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSetT<T>& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSetT<T>& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;

  /** boolean
     *
     * If true, this has to delete the BasisSet
     */
  bool OwnBasisSet;
  /// TwistAngle of this PWOrbitalSetT
  PosType TwistAngle;
  /// My basis set
  PWBasisPtr myBasisSet;
  /// number of basis
  IndexType BasisSetSize;
  /** pointer to matrix containing the coefficients
     *
     * makeClone makes a shallow copy and flag IsCloned
     */
  ValueMatrix* C;
  /// if true, do not clean up
  bool IsCloned;

  /** temporary array to perform gemm operation */
  Matrix<T> Temp;
};
} // namespace qmcplusplus
#endif
