//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RPA_JASTROW_H
#define QMCPLUSPLUS_RPA_JASTROW_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"

namespace qmcplusplus
{
/** JastrowBuilder using RPA functor
 *  Modification of RPAJastrow
 *
 */
class RPAJastrow : public WaveFunctionComponent
{
  using HandlerType = LRHandlerBase;
  using FuncType    = BsplineFunctor<RealType>;
  using GridType    = LinearGrid<RealType>;

public:
  RPAJastrow(ParticleSet& target);

  ~RPAJastrow() override;

  bool put(xmlNodePtr cur);

  void buildOrbital(const std::string& name,
                    const std::string& UL,
                    const std::string& US,
                    const std::string& RF,
                    RealType R,
                    RealType K);

  void makeShortRange();
  void makeLongRange();

  void setHandler(std::unique_ptr<HandlerType> Handler) { myHandler = std::move(Handler); };

  /** check out optimizable variables
    */
  void checkOutVariables(const opt_variables_type& o) override;

  /** check in an optimizable parameter
        * @param o a super set of optimizable variables
    */
  void checkInVariables(opt_variables_type& o) override;

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os) override;

  /** reset the parameters during optimizations
    */
  void resetParameters(const opt_variables_type& active) override;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  PsiValueType ratio(ParticleSet& P, int iat) override;
  GradType evalGrad(ParticleSet& P, int iat) override;
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  void restore(int iat) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

private:
  bool IgnoreSpin;
  bool DropLongRange;
  bool DropShortRange;
  RealType Rs;
  RealType Kc;
  RealType Rcut;
  std::string ID_Rs;
  std::string rpafunc;
  std::string MyName;

  ///@}
  /** main handler
   */
  std::unique_ptr<HandlerType> myHandler;
  ///object to handle the long-range part
  kSpaceJastrow* LongRangeRPA;
  ///@{objects to handle the short-range part
  ///two-body Jastrow function
  WaveFunctionComponent* ShortRangeRPA;
  ///numerical function owned by ShortRangeRPA
  FuncType* nfunc;
  GridType* myGrid;
  ParticleSet& targetPtcl;
  ///A list of WaveFunctionComponent*
  UPtrVector<WaveFunctionComponent> Psi;
};
} // namespace qmcplusplus
#endif
