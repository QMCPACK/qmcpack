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
  typedef std::map<std::string,ParticleSet*> PtclPoolType;
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
  void printIDs(std::vector<long> vi);
  void printInts(std::vector<int> vi);

  std::string xmlrootName;
  std::stringstream fname;
  int doWeights ;
  int weightFreq, weightLength, numSteps;
  int gridDivs;
  double overL;
  ParticleSetPool&  ptclPool;
  std::string ionroot;

  std::vector<int> walkersPerBlock;
  std::vector<int> pointsToCalculate;
  std::vector<std::vector<int> > Weights;
  std::vector<std::vector<long> > IDs, PIDs, realPIDs, realIDs;
//       HDF5_FW_float hdf_float_data;
  HDF5_FW_long hdf_ID_data,hdf_PID_data;

//       HDF5_FW_weights hdf_WGT_data;
  HDF5_FW_density hdf_Den_data;
  hid_t c_file;
  std::string WIDstring,PIDstring;
  int verbose,startStep;
  int gensTransferred;
};
}

#endif

