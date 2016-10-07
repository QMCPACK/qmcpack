//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_BACKFLOW_BUILDER_H
#define QMCPLUSPLUS_BACKFLOW_BUILDER_H
//#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee_kSpace.h"
#include "QMCWaveFunctions/Fermion/Backflow_eI.h"
#include "QMCWaveFunctions/Fermion/GaussianFunctor.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRRPABFeeHandlerTemp.h"
#include "Particle/ParticleSet.h"
#include "Configuration.h"
#include <map>
#include <cmath>

namespace qmcplusplus
{

class BackflowBuilder: public OrbitalBuilderBase
{

  typedef LRHandlerBase HandlerType;
  typedef LinearGrid<RealType> GridType;
  typedef std::map<std::string,ParticleSet*>   PtclPoolType;

public:

  BackflowBuilder(ParticleSet& p, PtclPoolType& pool, TrialWaveFunction& psi);

  ~BackflowBuilder();

  bool put(xmlNodePtr cur);

  BackflowTransformation* getBFTrans()
  {
    return BFTrans;
  }

  RealType cutOff;

private:

  PtclPoolType& ptclPool;
  BackflowTransformation* BFTrans;
  bool IgnoreSpin;
  RealType Rs;
  RealType Kc;
  RealType Rcut;
  bool OneBody;
  bool TwoBody;

  HandlerType* myHandler;

  void addOneBody(xmlNodePtr cur);

  void addTwoBody(xmlNodePtr cur);

  void addRPA(xmlNodePtr cur);

  void makeShortRange_oneBody();

  void makeLongRange_oneBody();

  void makeShortRange_twoBody(xmlNodePtr cur, Backflow_ee<BsplineFunctor<RealType> > *tbf, std::vector<int>& offsets);

  void makeLongRange_twoBody(xmlNodePtr cur, Backflow_ee_kSpace *tbf, std::vector<int>& offsets);


};

}

#endif
