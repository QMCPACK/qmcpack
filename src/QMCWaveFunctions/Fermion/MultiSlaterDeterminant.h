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


#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANT_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANT_ORBITAL_H
#include <Configuration.h>
#include <QMCWaveFunctions/Fermion/DiracDeterminant.h>
#include <QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h>
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"

namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief An AntiSymmetric WaveFunctionComponent composed of a linear combination of SlaterDeterminants.
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
class MultiSlaterDeterminant : public WaveFunctionComponent
{
public:
  void registerTimers();
  NewTimer &RatioTimer, &RatioGradTimer, &RatioAllTimer, &UpdateTimer, &EvaluateTimer;
  NewTimer &Ratio1Timer, &Ratio1GradTimer, &Ratio1AllTimer, &AccRejTimer, &evalOrbTimer;

  typedef DiracDeterminantBase* DiracDeterminantBasePtr;
  typedef SPOSet* SPOSetPtr;
  typedef SPOSetProxyForMSD* SPOSetProxyPtr;
  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t GradVector_t;
  typedef OrbitalSetTraits<ValueType>::HessMatrix_t HessMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType HessType;
  typedef Array<HessType, 3> HessArray_t;
  typedef TinyVector<HessType, 3> GGGType;
  typedef Vector<GGGType> GGGVector_t;
  typedef Matrix<GGGType> GGGMatrix_t;
  typedef ParticleSet::Walker_t Walker_t;


  ///constructor
  MultiSlaterDeterminant(ParticleSet& targetPtcl, SPOSetProxyPtr upspo, SPOSetProxyPtr dnspo);

  ///destructor
  ~MultiSlaterDeterminant();

  virtual void checkInVariables(opt_variables_type& active);
  virtual void checkOutVariables(const opt_variables_type& active);
  virtual void resetParameters(const opt_variables_type& active);
  virtual void reportStatus(std::ostream& os);

  void resetTargetParticleSet(ParticleSet& P);

  ///set BF pointers
  virtual void setBF(BackflowTransformation* BFTrans) {}

  virtual ValueType evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  virtual LogValueType evaluateLog(ParticleSet& P //const DistanceTableData* dtable,
                               ,
                               ParticleSet::ParticleGradient_t& G,
                               ParticleSet::ParticleLaplacian_t& L);

  virtual GradType evalGrad(ParticleSet& P, int iat);
  virtual ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  virtual ValueType ratio(ParticleSet& P, int iat);
  virtual void acceptMove(ParticleSet& P, int iat);
  virtual void restore(int iat);

  virtual void registerData(ParticleSet& P, WFBufferType& buf);
  virtual LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false);
  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  virtual WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const;
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   std::vector<RealType>& dlogpsi,
                                   std::vector<RealType>& dhpsioverpsi);

  virtual void resize(int, int);

  /**
    add a new SlaterDeterminant with coefficient c to the
    list of determinants
    */
  //int NumOrbitals_ground,NumOrbitals_total;
  int nels_up, nels_dn;
  int NumUniqueDets_up;
  int NumUniqueDets_dn;
  std::vector<int> DetID;

  int FirstIndex_up, LastIndex_up;
  int FirstIndex_dn, LastIndex_dn;

  std::map<std::string, int> SPOSetID;

  SPOSetProxyPtr spo_up;
  SPOSetProxyPtr spo_dn;

  std::vector<DiracDeterminantBasePtr> dets_up;
  std::vector<DiracDeterminantBasePtr> dets_dn;

  // map determinant in linear combination to unique det list
  std::vector<size_t> C2node_up;
  std::vector<size_t> C2node_dn;

  std::vector<ValueType> C;

  // lap(#uniqueDet,part#)
  ValueVector_t detValues_up;
  ValueVector_t detValues_dn;

  // UGLY, how do I get around this? I want to use GradMatrix instead...
  // grads(#uniqueDet,part#)
  std::vector<ParticleSet::ParticleGradient_t> grads_up;
  std::vector<ParticleSet::ParticleGradient_t> grads_dn;

  // lap(#uniqueDet,part#)
  std::vector<ParticleSet::ParticleLaplacian_t> lapls_up;
  std::vector<ParticleSet::ParticleLaplacian_t> lapls_dn;

  // grads(#uniqueDet,part#)
  std::vector<ParticleSet::ParticleGradient_t> tempgrad;

  // lap(#uniqueDet,part#)
  std::vector<ParticleSet::ParticleLaplacian_t> templapl;

  ValueType curRatio;
  ValueType psiCurrent;
  ValueVector_t detsRatios;
  ValueVector_t tempstorage_up;
  ValueVector_t tempstorage_dn;
  GradVector_t grad_temp;

  ParticleSet::ParticleGradient_t myG;
  ParticleSet::ParticleLaplacian_t myL;

  opt_variables_type myVars;

  // CSFs
  bool usingCSF;
  // coefficients of csfs, these are only used during optm
  std::vector<ValueType> CSFcoeff;
  // number of dets per csf
  std::vector<size_t> DetsPerCSF;
  // coefficiesize_tof csf expansion (smaller dimension)
  std::vector<RealType> CSFexpansion;
};

} // namespace qmcplusplus
#endif
