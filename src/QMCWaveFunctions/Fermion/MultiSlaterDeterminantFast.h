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
    
    
#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANTFAST_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANTFAST_ORBITAL_H
#include <Configuration.h>
#include <QMCWaveFunctions/OrbitalBase.h>
#include <QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h>
#include <QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h>
#include <QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h>
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief An AntiSymmetric OrbitalBase composed of a linear combination of SlaterDeterminants.
 *
 *\f[
 *MS({\bf R}) = \sum_n c_n S_n({\bf R})
 *\f].
 *
 *The (S)ingle(P)article(O)rbitalSet template parameter is an
 *engine which fills in single-particle orbital terms.
 *
 \f[
 \frac{\nabla_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n \nabla_i D_n}
 {\sum_{n=1}^M c_n D_n}
 \f]
 \f[
 \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N(\nabla_i
 S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
 \f]
 The Laplacian
 \f[
 \frac{\nabla^2_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N
 (\nabla_i^2S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
 \f]
 */
class MultiSlaterDeterminantFast: public OrbitalBase, public FermionBase
{

public:

  void registerTimers();
  NewTimer RatioTimer,RatioGradTimer,RatioAllTimer,UpdateTimer,EvaluateTimer;
  NewTimer Ratio1Timer,Ratio1GradTimer,Ratio1AllTimer, AccRejTimer;

  typedef MultiDiracDeterminantBase*    DiracDeterminantPtr;
  typedef SPOSetBase*                   SPOSetBasePtr;
  typedef SPOSetProxyForMSD*            SPOSetProxyPtr;
  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
  typedef OrbitalSetTraits<ValueType>::HessMatrix_t  HessMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType      HessType;
  typedef Array<HessType,3>                          HessArray_t;
  typedef TinyVector<HessType, OHMMS_DIM>                    GGGType;
  typedef Vector<GGGType>                            GGGVector_t;
  typedef Matrix<GGGType>                            GGGMatrix_t;
  typedef ParticleSet::Walker_t                      Walker_t;


  ///constructor
  MultiSlaterDeterminantFast(ParticleSet& targetPtcl,MultiDiracDeterminantBase* up, MultiDiracDeterminantBase* dn);

  ///destructor
  ~MultiSlaterDeterminantFast();

  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
  void resetParameters(const opt_variables_type& active);
  void reportStatus(std::ostream& os);

  void resetTargetParticleSet(ParticleSet& P);

  ///set BF pointers
  void setBF(BackflowTransformation* bf)
  {
    usingBF=true;
    BFTrans=bf;
    Dets[0]->setBF(bf);
    Dets[1]->setBF(bf);
  }

  ValueType
  evaluate_vgl_impl(ParticleSet& P
           ,ParticleSet::ParticleGradient_t& g_tmp
           ,ParticleSet::ParticleLaplacian_t& l_tmp);

  ValueType
  evaluate(ParticleSet& P
           ,ParticleSet::ParticleGradient_t& G
           ,ParticleSet::ParticleLaplacian_t& L);

  RealType
  evaluateLog(ParticleSet& P //const DistanceTableData* dtable,
              , ParticleSet::ParticleGradient_t& G
              , ParticleSet::ParticleLaplacian_t& L);

  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  ValueType evalGrad_impl(ParticleSet& P, int iat, bool newpos, GradType& g_at);

  ValueType ratio(ParticleSet& P, int iat);
  ValueType ratio_impl(ParticleSet& P, int iat);
  void evaluateRatiosAlltoOne(ParticleSet& P, int iat)
  {
    // the base class routine may probably work, just never tested.
    // it can also be highly optimized with a specialized implementation.
    APP_ABORT(" Need to implement MultiSlaterDeterminantFast::evaluateRatiosAlltoOne. \n");
  }

  void acceptMove(ParticleSet& P, int iat);
  void restore(int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi);

  void resize(int,int);
  void initialize();

  void testMSD(ParticleSet& P, int iat);

  size_t NP;
  size_t nels_up,nels_dn;
  size_t FirstIndex_up;
  size_t FirstIndex_dn;
  size_t ActiveSpin;
  bool usingCSF;
  bool IsCloned;
  ValueType curRatio;
  ValueType psiCurrent;

  // assume Dets[0]: up, Dets[1]:down
  std::vector<MultiDiracDeterminantBase*> Dets;
  std::map<std::string,size_t> SPOSetID;

  // map determinant in linear combination to unique det list
  std::vector<size_t>* C2node_up;
  std::vector<size_t>* C2node_dn;
  std::vector<RealType>* C;

  ParticleSet::ParticleGradient_t myG,myG_temp;
  ParticleSet::ParticleLaplacian_t myL,myL_temp;
  ValueVector_t laplSum_up;
  ValueVector_t laplSum_dn;

  //optimizable variable is shared with the clones
  opt_variables_type* myVars;
  // coefficients of csfs, these are only used during optm
  std::vector<RealType>* CSFcoeff;
  // number of dets per csf
  std::vector<size_t>* DetsPerCSF;
  // coefficient of csf expansion (smaller dimension)
  std::vector<RealType>* CSFexpansion;

  // transformation
  BackflowTransformation *BFTrans;
  bool usingBF;

  // temporary storage for evaluateDerivatives
  ParticleSet::ParticleGradient_t gmPG;
  Matrix<RealType> dpsia_up, dLa_up;
  Matrix<RealType> dpsia_dn, dLa_dn;
  Array<GradType,OHMMS_DIM> dGa_up, dGa_dn;

// debug, erase later
//      MultiSlaterDeterminant *msd;

};

}
#endif
