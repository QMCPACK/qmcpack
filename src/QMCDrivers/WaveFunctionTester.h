//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  bool checkSlaterDet; // flag to perform determinant-resolved test of SlaterDet
  std::string checkSlaterDetOption;
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
