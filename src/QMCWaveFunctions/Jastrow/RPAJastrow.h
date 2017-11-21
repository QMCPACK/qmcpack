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

#include "QMCWaveFunctions/OrbitalBase.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"


namespace qmcplusplus
{

/** JastrowBuilder using RPA functor
 *  Modification of RPAJastrow
 *
 */
struct RPAJastrow: public OrbitalBase
{
  typedef LRHandlerBase HandlerType;
  typedef BsplineFunctor<RealType> FuncType;
  typedef LinearGrid<RealType> GridType;

  RPAJastrow(ParticleSet& target, bool is_manager);

  ~RPAJastrow();

  bool put(xmlNodePtr cur);

  void buildOrbital(const std::string& name, const std::string& UL
                    , const std::string& US, const std::string& RF, RealType R, RealType K);

  void makeShortRange();
  void makeLongRange();

  void setHandler(HandlerType* Handler)
  {
    myHandler=Handler;
  };

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

  void resetTargetParticleSet(ParticleSet& P);

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  ValueType ratio(ParticleSet& P, int iat);
  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  void acceptMove(ParticleSet& P, int iat);

  void restore(int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  OrbitalBase* makeClone(ParticleSet& tqp) const;

private:

  bool IsManager;
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
  LRTwoBodyJastrow* LongRangeRPA;
  ///@{objects to handle the short-range part
  ///two-body Jastrow function
  OrbitalBase* ShortRangeRPA;
  ///numerical function owned by ShortRangeRPA
  FuncType* nfunc;
  GridType* myGrid;
  ///adaptor function to initialize nfunc
  ShortRangePartAdapter<RealType>* SRA;
  ParticleSet& targetPtcl;
  ///A list of OrbitalBase*
  std::vector<OrbitalBase*> Psi;
};
}
#endif
