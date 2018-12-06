//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COUNTING_JASTROW_BUILDER_H
#define QMCPLUSPLUS_COUNTING_JASTROW_BUILDER_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

//#include "Utilities/IteratorUtility.h"
//#include "Utilities/ProgressReportEngine.h"
//#include "OhmmsData/AttributeSet.h"
//#include "OhmmsData/ParameterSet.h"
//
//#include <string>


namespace qmcplusplus
{

struct CountingJastrowBuilder: public WaveFunctionComponentBuilder
{
  CountingJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur);
  
 private:
  typedef WaveFunctionComponent::RealType RT;
  
  ///jastrow/@name 
  std::string NameOpt;
  ///jastrow/@type
  std::string TypeOpt;
  ///jastrow/@function
  std::string Jastfunction;

  void createJC(xmlNodePtr cur);

};

}

#endif
