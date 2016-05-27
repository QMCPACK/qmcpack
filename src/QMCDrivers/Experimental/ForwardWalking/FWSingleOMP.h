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
#ifndef QMCPLUSPLUS_FWOMP_ANALYSIS_H
#define QMCPLUSPLUS_FWOMP_ANALYSIS_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/ForwardWalking/HDF5_FW.h"

namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class FWSingleOMP: public QMCDriver , public CloneManager
{
public:
  /// Constructor.
  FWSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
              HamiltonianPool& hpool,WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);

private:
  FWSingleOMP(const FWSingleOMP& a): QMCDriver(a), CloneManager(a) { }
  /// Copy operator (disabled).
  FWSingleOMP& operator=(const FWSingleOMP&)
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

/***************************************************************************
 * $RCSfile: VMCSingle.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCSingle.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
