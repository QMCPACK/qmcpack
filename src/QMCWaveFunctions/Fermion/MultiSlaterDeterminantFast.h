//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  typedef SPOSetBase*              SPOSetBasePtr;
  typedef SPOSetProxyForMSD*             SPOSetProxyPtr;
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
  void reportStatus(ostream& os);

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
  evaluate(ParticleSet& P
           ,ParticleSet::ParticleGradient_t& G
           ,ParticleSet::ParticleLaplacian_t& L);

  RealType
  evaluateLog(ParticleSet& P //const DistanceTableData* dtable,
              , ParticleSet::ParticleGradient_t& G
              , ParticleSet::ParticleLaplacian_t& L);

  RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L,
              PooledData<RealType>& buf,
              bool fillBuffer );

  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  ValueType ratio(ParticleSet& P, int iat
                  , ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL);

  ValueType ratio(ParticleSet& P, int iat);
  void acceptMove(ParticleSet& P, int iat);
  void restore(int iat);

  void update(ParticleSet& P
              , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL
              , int iat);
  RealType evaluateLog(ParticleSet& P,BufferType& buf);
  RealType registerData(ParticleSet& P, BufferType& buf);
  void registerDataForDerivatives(ParticleSet& P, BufferType& buf, int storageType=0);
  virtual void memoryUsage_DataForDerivatives(ParticleSet& P,long& orbs_only,long& orbs, long& invs, long& dets)
  {
    Dets[0]->memoryUsage_DataForDerivatives(P,orbs_only,orbs,invs,dets);
    Dets[1]->memoryUsage_DataForDerivatives(P,orbs_only,orbs,invs,dets);
  }
  RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, BufferType& buf);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           vector<RealType>& dlogpsi,
                           vector<RealType>& dhpsioverpsi);

//      void evaluateDerivatives(ParticleSet& P,
//                                       const opt_variables_type& optvars,
//                                       vector<RealType>& dlogpsi,
//                                       vector<RealType>& dhpsioverpsi,
//                                       PooledData<RealType>& buf);

  void resize(int,int);

  void testMSD(ParticleSet& P, int iat);

  int NP;
  int nels_up,nels_dn;
  int FirstIndex_up;
  int FirstIndex_dn;
  bool usingCSF;

  // assume Dets[0]: up, Dets[1]:down
  vector<MultiDiracDeterminantBase*> Dets;

  vector<int> DetID;

  map<string,int> SPOSetID;

  // map determinant in linear combination to unique det list
  vector<int> C2node_up;
  vector<int> C2node_dn;

  vector<RealType> C;

  ValueType curRatio;
  ValueType psiCurrent;

  ValueType* FirstAddressOfG;
  ValueType* LastAddressOfG;

  ParticleSet::ParticleGradient_t myG,myG_temp;
  ParticleSet::ParticleLaplacian_t myL,myL_temp;
  ValueVector_t laplSum_up;
  ValueVector_t laplSum_dn;

  opt_variables_type myVars;

// CSFs
  // coefficients of csfs, these are only used during optm
  vector<RealType> CSFcoeff;
  // number of dets per csf
  vector<int> DetsPerCSF;
  // coefficient of csf expansion (smaller dimension)
  vector<RealType> CSFexpansion;

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
/***************************************************************************
 * $RCSfile$   $Author: miguel.mmorales $
 * $Revision: 4791 $   $Date: 2010-05-12 12:08:35 -0500 (Wed, 12 May 2010) $
 * $Id: MultiSlaterDeterminantFast.h 4791 2010-05-12 17:08:35Z miguel.mmorales $
 ***************************************************************************/
