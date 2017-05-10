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


#ifndef QMCPLUSPLUS_FWMPI_ANALYSIS_H
#define QMCPLUSPLUS_FWMPI_ANALYSIS_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/ForwardWalking/HDF5_FW.h"

namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class FWSingleMPI: public QMCDriver , public CloneManager
{
public:
  /// Constructor.
  FWSingleMPI(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
              HamiltonianPool& hpool, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);

private:
  FWSingleMPI(const FWSingleMPI& a): QMCDriver(a), CloneManager(a) { }
  /// Copy operator (disabled).
  FWSingleMPI& operator=(const FWSingleMPI&)
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
  int doWeights, doObservables, doDat;
  int weightFreq, weightLength, numSteps;

  std::vector<int> walkersPerBlock;
  std::vector<int> pointsToCalculate;
  std::vector<std::vector<int> > Weights;
  std::vector<std::vector<long> > IDs, PIDs, realPIDs, realIDs;
  HDF5_FW_float hdf_float_data;
  HDF5_FW_long hdf_ID_data,hdf_PID_data;
  HDF5_FW_observables hdf_OBS_data;
  HDF5_FW_weights hdf_WGT_data;
  hid_t c_file;
  std::string WIDstring,PIDstring;
  int verbose,startStep;
  int gensTransferred;
};
}

#endif

