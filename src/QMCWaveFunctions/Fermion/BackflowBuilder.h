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
  typedef map<string,ParticleSet*>   PtclPoolType;

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

  void makeShortRange_twoBody(xmlNodePtr cur, Backflow_ee<BsplineFunctor<double> > *tbf, vector<int>& offsets);

  void makeLongRange_twoBody(xmlNodePtr cur, Backflow_ee_kSpace *tbf, vector<int>& offsets);


};

}

#endif
