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
#ifndef QMCPLUSPLUS_FW_ANALYSIS_H
#define QMCPLUSPLUS_FW_ANALYSIS_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class FWSingle: public QMCDriver
{
public:
  /// Constructor.
  FWSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);

private:
  FWSingle(const FWSingle& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  FWSingle& operator=(const FWSingle&)
  {
    return *this;
  }

  //main drivers
  void FWOneStep();
  void transferParentsOneGeneration();

  //helper functions
  void fillIDMatrix();
  void fillWalkerPositionsandWeights(int nstep);
  void readInLong(int step, std::string IDstring, std::vector<long>& data_out);
  void readInFloat(int step, std::vector<float>& data_out);
  void resetWeights();
  int getNumberOfSamples(int omittedSteps);
  //debugging functions
  void printIDs(std::vector<long> vi);
  void printInts(std::vector<int> vi);

  std::string xmlrootName;
  std::stringstream fname;
  int weightFreq, weightLength, numSteps;

  std::vector<int> walkersPerBlock;
  std::vector<int> pointsToCalculate;
  std::vector<std::vector<int> > Weights;
  std::vector<std::vector<long> > IDs, PIDs, realPIDs, realIDs;
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
