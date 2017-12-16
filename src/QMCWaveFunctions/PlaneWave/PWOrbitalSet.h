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
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{

class PWOrbitalSet: public SPOSetBase
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
  PWOrbitalSet(): OwnBasisSet(false)
  {
  }

  /** delete BasisSet only it owns this
   *
   * Builder takes care of who owns what
   */
  ~PWOrbitalSet();

  SPOSetBase* makeClone() const;
  /** resize  the orbital base
   * @param bset PWBasis
   * @param nbands number of bands
   * @param cleaup if true, owns PWBasis. Will clean up.
   */
  void resize(PWBasisPtr bset, int nbands, bool cleanup=false);

  /** Builder class takes care of the assertion
  */
  void addVector(const std::vector<ComplexType>& coefs,int jorb);
  void addVector(const std::vector<RealType>& coefs,int jorb);

  void resetParameters(const opt_variables_type& optVariables);

  void setOrbitalSetSize(int norbs);

  void resetTargetParticleSet(ParticleSet& P);

  inline ValueType evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return BLAS::dot(BasisSetSize,(*C)[ib],myBasisSet->Zv.data());
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

  /** boolean
   *
   * If true, this has to delete the BasisSet
   */
  bool OwnBasisSet;
  ///TwistAngle of this PWOrbitalSet
  PosType TwistAngle;
  ///My basis set
  PWBasisPtr myBasisSet;
  /////Plane-wave coefficients: (iband,g-vector)
  //Matrix<ValueType> Coefs;
  /** temporary array to perform gemm operation */
  Matrix<ValueType> Temp;
};
}
#endif
