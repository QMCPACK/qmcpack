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


#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANTWITHBACKFLOW_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANTWITHBACKFLOW_ORBITAL_H
#include <Configuration.h>
#include <QMCWaveFunctions/Fermion/DiracDeterminant.h>
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include <QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h>
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief MultiSlaterDeterminantWithBackflow
 */
class MultiSlaterDeterminantWithBackflow : public MultiSlaterDeterminant
{
public:
  ///constructor
  MultiSlaterDeterminantWithBackflow(ParticleSet& targetPtcl,
                                     SPOSetProxyPtr upspo,
                                     SPOSetProxyPtr dnspo,
                                     BackflowTransformation* tr);

  ///destructor
  ~MultiSlaterDeterminantWithBackflow();

  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
  void resetParameters(const opt_variables_type& active);
  void reportStatus(std::ostream& os);

  ///set BF pointers
  void setBF(BackflowTransformation* bf)
  {
    BFTrans = bf;
    for (int i = 0; i < dets_up.size(); i++)
      dets_up[i]->setBF(bf);
    for (int i = 0; i < dets_dn.size(); i++)
      dets_dn[i]->setBF(bf);
  }

  ValueType evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  LogValueType evaluateLog(ParticleSet& P //const DistanceTableData* dtable,
                           ,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L);

  GradType evalGrad(ParticleSet& P, int iat);
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  PsiValueType ratio(ParticleSet& P, int iat);
  void acceptMove(ParticleSet& P, int iat);
  void restore(int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);
  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false);
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const;
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi);

  void resize(int, int);

  // transformation
  BackflowTransformation* BFTrans;

  // temporary storage for evaluateDerivatives
  Matrix<RealType> dpsia_up, dLa_up;
  Matrix<RealType> dpsia_dn, dLa_dn;
  Array<GradType, 3> dGa_up, dGa_dn;
};

} // namespace qmcplusplus
#endif
