//////////////////////////////////////////////////////////////////
// (c) Copyright 2009- by Jeremy McMinis and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef QMCPLUSPLUS_FROMP_ANALYSIS_H
#define QMCPLUSPLUS_FROMP_ANALYSIS_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/ForwardWalking/HDF5_FW.h"
#include "QMCApp/ParticleSetPool.h"

namespace qmcplusplus
{

class QMCUpdateBase;
class ParticleSetPool;
/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class FRSingleOMP: public QMCDriver , public CloneManager
{
  typedef map<string,ParticleSet*> PtclPoolType;
public:
  /// Constructor.
  FRSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
              HamiltonianPool& hpool, ParticleSetPool& pset, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);

private:
  FRSingleOMP(const FRSingleOMP& a): QMCDriver(a), CloneManager(a), ptclPool(a.ptclPool) { }
  /// Copy operator (disabled).
  FRSingleOMP& operator=(const FRSingleOMP&)
  {
    return *this;
  }

  //main drivers
  void FWOneStep();
  void transferParentsOneGeneration();

  //helper functions
  void fillIDMatrix();
  void resetWeights();
  int getNumberOfSamples(int omittedSteps);
  //debugging functions
  void printIDs(vector<long> vi);
  void printInts(vector<int> vi);

  string xmlrootName;
  stringstream fname;
  int doWeights ;
  int weightFreq, weightLength, numSteps;
  int gridDivs;
  double overL;
  ParticleSetPool&  ptclPool;
  string ionroot;

  vector<int> walkersPerBlock;
  vector<int> pointsToCalculate;
  vector<vector<int> > Weights;
  vector<vector<long> > IDs, PIDs, realPIDs, realIDs;
//       HDF5_FW_float hdf_float_data;
  HDF5_FW_long hdf_ID_data,hdf_PID_data;

//       HDF5_FW_weights hdf_WGT_data;
  HDF5_FW_density hdf_Den_data;
  hid_t c_file;
  string WIDstring,PIDstring;
  int verbose,startStep;
  int gensTransferred;
};
}

#endif

/***************************************************************************
 * $RCSfile: VMCSingle.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCSingle.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
