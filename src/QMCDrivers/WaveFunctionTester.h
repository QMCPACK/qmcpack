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
#ifndef QMCPLUSPLUS_WAVEFUNCTIONTEST_H
#define QMCPLUSPLUS_WAVEFUNCTIONTEST_H

#include "QMCDrivers/QMCDriver.h"
#include "QMCApp/ParticleSetPool.h"
namespace qmcplusplus
{


/** Information for output of relative error in wavefunction derivatives
  vs. finite difference delta.
*/
class FiniteDiffErrData : public QMCTraits
{
public:
  FiniteDiffErrData();
  bool put(xmlNodePtr q);

  int particleIndex;
  int gradientComponentIndex;
  std::string outputFile;
};

/** Test the correctness of TrialWaveFunction for the values,
    gradients and laplacians
*/
class WaveFunctionTester: public QMCDriver
{
public:
  /// Constructor.
  WaveFunctionTester(MCWalkerConfiguration& w,
                     TrialWaveFunction& psi,
                     QMCHamiltonian& h,
                     ParticleSetPool& ptclPool,
                     WaveFunctionPool& ppool);

  ~WaveFunctionTester();

  bool run();
  bool put(xmlNodePtr q);

private:
  ParticleSetPool &PtclPool;
  ParticleSet::ParticlePos_t deltaR;
  std::string checkRatio, checkClone, checkHamPbyP, sourceName, wftricks, checkEloc;
  std::string checkBasic, checkRatioV;
  xmlNodePtr myNode;
  double deltaParam;
  double toleranceParam;
  bool outputDeltaVsError;
  FiniteDiffErrData DeltaVsError;
 
  /// Copy Constructor (disabled)
  WaveFunctionTester(const WaveFunctionTester& a):
    QMCDriver(a), PtclPool(a.PtclPool) { }
  /// Copy Operator (disabled)
  WaveFunctionTester& operator=(const WaveFunctionTester&)
  {
    return *this;
  }
  /** basic tests for G and L */
  void runBasicTest();
  /** the basic ratios check */
  void runRatioTest();
  void runRatioTest2();
  /** test ratios with virtual moves */
  void runRatioV();
  /** test clone implementations of new wavefunctions and operators */
  void runCloneTest();
  void runDerivTest();
  void runDerivNLPPTest();
  void runDerivCloneTest();
  void runGradSourceTest();
  void runZeroVarianceTest();
  void runwftricks();
  void runNodePlot();
  void printEloc();

  // compute numerical gradient and laplacian
  void computeNumericalGrad(RealType delta,
                            ParticleSet::ParticleGradient_t &G_fd,
                            ParticleSet::ParticleLaplacian_t &L_fd);

  bool checkGradients(int lower_iat, int upper_iat,
                      ParticleSet::ParticleGradient_t &G,
                      ParticleSet::ParticleLaplacian_t &L,
                      ParticleSet::ParticleGradient_t &G_fd,
                      ParticleSet::ParticleLaplacian_t &L_fd,
                      std::stringstream &log,
                      int indent=0);

  bool checkGradientAtConfiguration(MCWalkerConfiguration::Walker_t* W1,
                                    std::stringstream &fail_log,
                                    bool &ignore);

  //vector<RealType> Mv3(std::vector<std::vector<RealType> >& M, std::vector<RealType>& v);

  std::ofstream fout;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
