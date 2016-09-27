//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by:  Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////



#ifndef QMCPLUSPLUS_REPTILEESTIMATOR_H
#define QMCPLUSPLUS_REPTILEESTIMATOR_H

#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/RecordProperty.h"
#include "Estimators/accumulators.h"
#include "QMCDrivers/MultiChain.h"

namespace qmcplusplus
{

struct ReptileEstimator: public QMCTraits
{
  typedef accumulator_set<RealType> accumulator_type;

  ReptileEstimator() {}
  ReptileEstimator(MultiChain* polymer)
  {
    Reptile=polymer;
  }
  virtual ~ReptileEstimator() {}
//     virtual void put(xmlNodePtr cur) {}
//     virtual void put(xmlNodePtr cur, MCWalkerConfiguration& refWalker) {put(cur);}
  virtual void evaluate(MultiChain::iterator first, MultiChain::iterator last, int ipsi)=0;
  virtual void addNames(std::vector<std::string>& names)=0;
  virtual void setValues(std::vector<accumulator_type>& data, RealType weight)
  {
    std::vector<RealType>::iterator Vit(Values.begin());
    std::vector<accumulator_type>::iterator Dit(data.begin());
    Dit+=myIndex;
    //these must be the same length.
    for(; Vit!=Values.end(); Vit++,Dit++)
      (*Dit)((*Vit),weight);
  };



  int myIndex;
  std::vector<RealType> Values;
  MultiChain* Reptile;
};
}
#endif
