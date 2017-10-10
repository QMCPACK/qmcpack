//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_POLYMERESTIMATOR_H
#include "QMCDrivers/PolymerChain.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

class PolymerEstimator
{

  PolymerChain& Reptile;
  int Middle;
  int Counter;

  int ReptileLength;
  int nPsi;
  int EpotLength;
  int EpotSize;
  int ElocSize;
  std::ofstream* fout;
  std::ofstream* OutLocEn;
  std::ofstream* OutPotEn;

  Vector<double> AvgLocalEnergy;
  Vector<double> AvgPotentialEnergy;
  Vector<double> TotalWeight;

  std::vector<double> PEavg;
  std::vector<double> PE2;

public:

  PolymerEstimator(PolymerChain& reptile, int npsi=1);

  ~PolymerEstimator()
  {
    clean();
  }

  void clean();

  void resetReportSettings(const std::string& aname);

  inline void reset()
  {
    Counter = 0;
    for(int i=0; i<PEavg.size(); i++)
      PEavg[i]=0.0;
    for(int i=0; i<PE2.size(); i++)
      PE2[i]=0.0;
    AvgLocalEnergy=0.0;
    AvgPotentialEnergy=0.0;
    TotalWeight=0.0;
  }

  void report(int iter);

  void accumulate();

};
}
#endif
