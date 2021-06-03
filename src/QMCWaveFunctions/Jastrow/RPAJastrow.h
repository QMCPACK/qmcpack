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
struct RPAJastrow : public WaveFunctionComponent
{
  typedef LRHandlerBase HandlerType;
  typedef BsplineFunctor<RealType> FuncType;
  typedef LinearGrid<RealType> GridType;

  RPAJastrow(ParticleSet& target);

  ~RPAJastrow();

  bool put(xmlNodePtr cur);

  void buildOrbital(const std::string& name,
                    const std::string& UL,
                    const std::string& US,
                    const std::string& RF,
                    RealType R,
                    RealType K);

  void makeShortRange();
  void makeLongRange();

  void setHandler(HandlerType* Handler) { myHandler = Handler; };

  /** check out optimizable variables
    */
  void checkOutVariables(const opt_variables_type& o);

  /** check in an optimizable parameter
        * @param o a super set of optimizable variables
    */
  void checkInVariables(opt_variables_type& o);

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os);

  /** reset the parameters during optimizations
    */
  void resetParameters(const opt_variables_type& active);

  LogValueType evaluateLog(const ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  PsiValueType ratio(ParticleSet& P, int iat);
  GradType evalGrad(ParticleSet& P, int iat);
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false);

  void restore(int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  WaveFunctionComponent* makeClone(ParticleSet& tqp) const;

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
  HandlerType* myHandler;
  ///object to handle the long-range part
  kSpaceJastrow* LongRangeRPA;
  ///@{objects to handle the short-range part
  ///two-body Jastrow function
  WaveFunctionComponent* ShortRangeRPA;
  ///numerical function owned by ShortRangeRPA
  FuncType* nfunc;
  GridType* myGrid;
  ///adaptor function to initialize nfunc
  ShortRangePartAdapter<RealType>* SRA;
  ParticleSet& targetPtcl;
  ///A list of WaveFunctionComponent*
  std::vector<WaveFunctionComponent*> Psi;
};
} // namespace qmcplusplus
#endif
