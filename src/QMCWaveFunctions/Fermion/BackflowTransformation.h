//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BACKFLOW_TRANSFORMATION_H
#define QMCPLUSPLUS_BACKFLOW_TRANSFORMATION_H

#include "Configuration.h"
#include <map>
#include <cmath>
#include "Particle/ParticleSet.h"
#include "DistanceTable.h"
#include "Particle/ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{
class BackflowTransformation
{
public:
  typedef BackflowFunctionBase::WFBufferType WFBufferType;

  // All BF quantities should be real, so eliminating complex (ValueType) possibility
  enum
  {
    DIM = OHMMS_DIM
  };
  typedef OHMMS_PRECISION RealType;
  typedef int IndexType;
  typedef TinyVector<RealType, DIM> PosType;
  typedef TinyVector<RealType, DIM> GradType;
  typedef Tensor<RealType, DIM> HessType;
  typedef Vector<IndexType> IndexVector_t;
  typedef Vector<GradType> GradVector_t;
  typedef Matrix<GradType> GradMatrix_t;
  typedef Vector<HessType> HessVector_t;
  typedef Matrix<HessType> HessMatrix_t;

  typedef Array<HessType, 3> HessArray_t;

  typedef std::map<std::string, ParticleSet*> PtclPoolType;
  //typedef Array<GradType,3>       GradArray_t;
  //typedef Array<PosType,3>        PosArray_t;

  ///number of quantum particles
  int NumTargets;

  /// active particle in pbyp moves
  int activeParticle;

  /// quasiparticle coordinates
  ParticleSet QP;

  // number of variational parameters
  int numParams;

  /** current update mode */
  int UpdateMode;

  /** enum for a update mode */
  enum
  {
    ORB_PBYP_RATIO,   /*!< particle-by-particle ratio only */
    ORB_PBYP_ALL,     /*!< particle-by-particle, update Value-Gradient-Laplacian */
    ORB_PBYP_PARTIAL, /*!< particle-by-particle, update Value and Grdient */
    ORB_WALKER,       /*!< walker update */
    ORB_ALLWALKER     /*!< all walkers update */
  };

  // map index of variables from local arrays to outside world
  std::map<int, int> optIndexMap;

  // cutoff of radial funtions
  RealType cutOff;

  // pos of first optimizable variable in global array
  int numVarBefore;

  /// Distance Table
  const int myTableIndex_;

  // matrix of laplacians
  // /vec{B(i)} = sum_{k} /grad_{k}^2 /vec{x_i}
  GradVector_t Bmat;

  GradMatrix_t Bmat_full, Bmat_temp;

  // matrix of first derivatives
  // A(i,j)[a,b] = (Grad_i)_a (x_j)_b
  //               i,j:particle index
  //               a,b=(x,y,z)
  // notice that A(i,j) is a symmetric matrix, improve later
  HessMatrix_t Amat, Amat_temp;

  // \nabla_a A_{i,j}^{\alpha,\beta}
  // derivative of A matrix with respect to var. prms.
  HessArray_t Xmat;

  // \sum_i \nabla_a B_{i,j}^{\alpha}
  GradMatrix_t Ymat;

  // \nabla_a x_i^{\alpha}
  GradMatrix_t Cmat;

  RealType *FirstOfP, *LastOfP;
  RealType *FirstOfA, *LastOfA;
  RealType *FirstOfB, *LastOfB;
  RealType *FirstOfA_temp, *LastOfA_temp;
  RealType *FirstOfB_temp, *LastOfB_temp;

  // Identity
  HessType HESS_ID;
  HessType DummyHess;

  std::vector<std::unique_ptr<BackflowFunctionBase>> bfFuns;

  std::map<std::string, int> sources;
  std::vector<std::string> names;

  /// new qp coordinates for pbyp moves.
  ParticleSet::ParticlePos_t newQP;
  ParticleSet::ParticlePos_t oldQP;

  //Vector<PosType> storeQP;
  Vector<PosType> storeQP;

  /// store index of qp coordinates that changed during pbyp move
  std::vector<int> indexQP, index;

  opt_variables_type myVars;

  BackflowTransformation(ParticleSet& els);

  void copyFrom(const BackflowTransformation& tr, ParticleSet& targetPtcl);

  std::unique_ptr<BackflowTransformation> makeClone(ParticleSet& tqp) const;

  ~BackflowTransformation();

  bool put(xmlNodePtr cur) { return true; }

  void acceptMove(const ParticleSet& P, int iat);

  void restore(int iat = 0);

  void checkInVariables(opt_variables_type& active);

  void reportStatus(std::ostream& os);

  void checkOutVariables(const opt_variables_type& active);

  bool isOptimizable();

  void resetParameters(const opt_variables_type& active);

  void registerData(ParticleSet& P, WFBufferType& buf);

  void updateBuffer(ParticleSet& P, WFBufferType& buf, bool redo);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /** calculate quasi-particle coordinates only
   */
  void transformOnly(const ParticleSet& P);

  /** calculate new quasi-particle coordinates after pbyp move
   */
  void evaluatePbyP(const ParticleSet& P, int iat);

  /** calculate new quasi-particle coordinates after pbyp move
   */
  void evaluatePbyPWithGrad(const ParticleSet& P, int iat);

  /** calculate new quasi-particle coordinates after pbyp move
   */
  void evaluatePbyPAll(const ParticleSet& P, int iat);

  /** calculate only Bmat. Assume that QP and Amat are current
   *  This is used in pbyp moves, in updateBuffer()
   */
  void evaluateBmatOnly(const ParticleSet& P, int iat);

  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  void evaluate(const ParticleSet& P);

  /** calculate quasi-particle coordinates and store in Pnew
   */
  void evaluate(const ParticleSet& P, ParticleSet& Pnew);

  void evaluateDerivatives(const ParticleSet& P);

  void testDeriv(const ParticleSet& P);

  void testPbyP(ParticleSet& P);
};

} // namespace qmcplusplus

#endif
