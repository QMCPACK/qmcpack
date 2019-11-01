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


/** @file PWOrbitalSet.h
 * @brief Definition of member functions of Plane-wave basis set
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{
class PWOrbitalSet : public SPOSet
{
public:
  typedef PWBasis BasisSet_t;
  typedef PWBasis* PWBasisPtr;

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
  PWOrbitalSet() : OwnBasisSet(false), myBasisSet(nullptr), BasisSetSize(0), C(nullptr), IsCloned(false)
  {
    className = "PWOrbitalSet";
  }

  /** delete BasisSet only it owns this
   *
   * Builder takes care of who owns what
   */
  ~PWOrbitalSet();

  SPOSet* makeClone() const;
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

  void resetParameters(const opt_variables_type& optVariables);

  void setOrbitalSetSize(int norbs);

  void resetTargetParticleSet(ParticleSet& P);

  inline ValueType evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return BLAS::dot(BasisSetSize, (*C)[ib], myBasisSet->Zv.data());
  }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet);

  /** boolean
   *
   * If true, this has to delete the BasisSet
   */
  bool OwnBasisSet;
  ///TwistAngle of this PWOrbitalSet
  PosType TwistAngle;
  ///My basis set
  PWBasisPtr myBasisSet;
  ///number of basis
  IndexType BasisSetSize;
  /** pointer to matrix containing the coefficients
   *
   * makeClone makes a shallow copy and flag IsCloned
   */
  ValueMatrix_t* C;
  ///if true, do not clean up
  bool IsCloned;
  /////Plane-wave coefficients: (iband,g-vector)
  //Matrix<ValueType> Coefs;
  /** temporary array to perform gemm operation */
  Matrix<ValueType> Temp;
};
} // namespace qmcplusplus
#endif
