//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  ofstream* fout;
  ofstream* OutLocEn;
  ofstream* OutPotEn;

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

  void resetReportSettings(const string& aname);

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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 759 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: PolymerEstimator.h 759 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
