//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BACKFLOW_FUNCTIONBASE_H
#define QMCPLUSPLUS_BACKFLOW_FUNCTIONBASE_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Configuration.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{
/**  Base class for backflow transformations.
 *  FT is an optimizable functor class that implements the radial function
 *  Any class used for Jastrow functions should work
 */
class BackflowFunctionBase //: public OrbitalSetTraits<QMCTraits::ValueType>
{
public:
  using Walker_t     = ParticleSet::Walker_t;
  using WFBufferType = Walker_t::WFBuffer_t;

  // All BF quantities should be real, so eliminating complex (ValueType) possibility
  enum
  {
    DIM = OHMMS_DIM
  };
  using RealType    = OHMMS_PRECISION;
  using IndexType   = int;
  using PosType     = TinyVector<RealType, DIM>;
  using GradType    = TinyVector<RealType, DIM>;
  using HessType    = Tensor<RealType, DIM>;
  using IndexVector = Vector<IndexType>;
  using GradVector  = Vector<GradType>;
  using GradMatrix  = Matrix<GradType>;
  using HessVector  = Vector<HessType>;
  using HessMatrix  = Matrix<HessType>;

  using HessArray = Array<HessType, 3>;
  //using GradArray_t = Array<GradType,3>      ;
  //using PosArray_t = Array<PosType,3>       ;

  /** enum for a update mode */
  enum
  {
    ORB_PBYP_RATIO,   /*!< particle-by-particle ratio only */
    ORB_PBYP_ALL,     /*!< particle-by-particle, update Value-Gradient-Laplacian */
    ORB_PBYP_PARTIAL, /*!< particle-by-particle, update Value and Grdient */
    ORB_WALKER,       /*!< walker update */
    ORB_ALLWALKER     /*!< all walkers update */
  };

  ///Reference to the center
  ParticleSet& CenterSys;
  ///number of centers, e.g., ions
  int NumCenters;
  ///number of quantum particles
  int NumTargets;
  // number of variational parameters own by the radial function
  int numParams;
  // index of first parameter in derivative array
  int indexOfFirstParam;
  // temporary storage for derivatives
  std::vector<TinyVector<RealType, 3>> derivs;

  // mmorales: all quantities produced by BF transformations
  //           should be real, so change everything here to ???<RealType>
  Matrix<PosType> UIJ;
  Vector<PosType> UIJ_temp;

  HessMatrix AIJ;
  HessVector AIJ_temp;

  GradMatrix BIJ;
  GradVector BIJ_temp;

  RealType *FirstOfU, *LastOfU;
  RealType *FirstOfA, *LastOfA;
  RealType *FirstOfB, *LastOfB;

  bool uniqueFunctions;
  opt_variables_type myVars;

  BackflowFunctionBase(ParticleSet& ions, ParticleSet& els)
      : CenterSys(ions), numParams(0), indexOfFirstParam(-1), uniqueFunctions(false)
  {
    NumCenters = CenterSys.getTotalNum(); // in case
    NumTargets = els.getTotalNum();
  }

  void resize(int NT, int NC)
  {
    NumTargets = NT;
    NumCenters = NC;
    UIJ.resize(NumTargets, NumCenters);
    UIJ = 0;
    AIJ.resize(NumTargets, NumCenters);
    AIJ = 0;
    BIJ.resize(NumTargets, NumCenters);
    BIJ = 0;
    UIJ_temp.resize(NumCenters);
    UIJ_temp = 0;
    AIJ_temp.resize(NumCenters);
    AIJ_temp = 0;
    BIJ_temp.resize(NumCenters);
    BIJ_temp = 0;
  }

  virtual std::unique_ptr<BackflowFunctionBase> makeClone(ParticleSet& tqp) const = 0;

  virtual ~BackflowFunctionBase(){};

  virtual void acceptMove(int iat, int UpdateType) = 0;

  virtual void restore(int iat, int UpdateType) = 0;

  virtual void reportStatus(std::ostream& os) = 0;

  virtual void resetParameters(const opt_variables_type& active) = 0;

  virtual void checkInVariables(opt_variables_type& active) = 0;

  virtual void checkOutVariables(const opt_variables_type& active) = 0;

  virtual bool isOptimizable() = 0;

  virtual int indexOffset() = 0;

  // Note: numParams should be set in Builder class, so it is known here
  inline int setParamIndex(int n)
  {
    indexOfFirstParam = n;
    return numParams;
  }

  virtual void registerData(WFBufferType& buf) = 0;

  void updateBuffer(WFBufferType& buf)
  {
    buf.put(FirstOfU, LastOfU);
    buf.put(FirstOfA, LastOfA);
    buf.put(FirstOfB, LastOfB);
  }

  void copyFromBuffer(WFBufferType& buf)
  {
    buf.get(FirstOfU, LastOfU);
    buf.get(FirstOfA, LastOfA);
    buf.get(FirstOfB, LastOfB);
  }

  /** calculate quasi-particle coordinates only
   */
  virtual void evaluate(const ParticleSet& P, ParticleSet& QP) = 0;

  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  virtual void evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix& Bmat, HessMatrix& Amat) = 0;

  /** calculate quasi-particle coordinates after pbyp move
   */
  virtual void evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos& newQP, const std::vector<int>& index) = 0;

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  virtual void evaluatePbyP(const ParticleSet& P,
                            ParticleSet::ParticlePos& newQP,
                            const std::vector<int>& index,
                            HessMatrix& Amat) = 0;

  /** calculate quasi-particle coordinates, Bmat and Amat after pbyp move
   */
  virtual void evaluatePbyP(const ParticleSet& P,
                            ParticleSet::ParticlePos& newQP,
                            const std::vector<int>& index,
                            GradMatrix& Bmat,
                            HessMatrix& Amat) = 0;

  /** calculate quasi-particle coordinates after pbyp move
   */
  virtual void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos& newQP) = 0;

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  virtual void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos& newQP, HessMatrix& Amat) = 0;

  /** calculate quasi-particle coordinates, Bmat and Amat after pbyp move
   */
  virtual void evaluatePbyP(const ParticleSet& P,
                            int iat,
                            ParticleSet::ParticlePos& newQP,
                            GradMatrix& Bmat,
                            HessMatrix& Amat) = 0;

  /** calculate only Bmat
   *  This is used in pbyp moves, in updateBuffer()
   */
  virtual void evaluateBmatOnly(const ParticleSet& P, GradMatrix& Bmat_full) = 0;

  /** calculate quasi-particle coordinates, Bmat and Amat
   *  calculate derivatives wrt to variational parameters
   */
  virtual void evaluateWithDerivatives(const ParticleSet& P,
                                       ParticleSet& QP,
                                       GradMatrix& Bmat,
                                       HessMatrix& Amat,
                                       GradMatrix& Cmat,
                                       GradMatrix& Ymat,
                                       HessArray& Xmat) = 0;
};

} // namespace qmcplusplus

#endif
