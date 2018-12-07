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
  ///jastrow/@source
  std::string SourceOpt;
  ///jastrow/@region
  std::string RegionOpt;

  template<template<class> class RegionType>
  bool createCJ(xmlNodePtr cur);
};

}

#endif
