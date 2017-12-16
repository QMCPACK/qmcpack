//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_VMC_CUDA_H
#define QMCPLUSPLUS_VMC_CUDA_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCcuda: public QMCDriver
{
public:
  /// Constructor.
  VMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool);

  bool run();
  bool runWithDrift();

  /// advance walkers without drift
  void advanceWalkers();
  /// advance walkers with drift
  void advanceWalkersWithDrift();

  bool put(xmlNodePtr cur);
  RealType fillOverlapHamiltonianMatrices(Matrix<RealType>& LeftM, Matrix<RealType>& RightM);
  inline void setOpt(bool o)
  {
    forOpt=o;
  };

private:
  ///previous steps
  int prevSteps;
  ///previous stepsbetweensamples
  int prevStepsBetweenSamples;
  /// tau/mass
  RealType m_tauovermass;
  /// Whether to use drift or not
  std::string UseDrift;
  ///period for walker dump
  int myPeriod4WalkerDump;
  /// Copy Constructor (disabled)
  VMCcuda(const VMCcuda& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCcuda& operator=(const VMCcuda&)
  {
    return *this;
  }
  ///hide initialization from the main function
  bool checkBounds (std::vector<PosType> &newpos, std::vector<bool> &valid);

  void resetRun();

  opt_variables_type dummy;
  int numParams;
  Matrix<RealType> d_logpsi_dalpha, d_hpsioverpsi_dalpha;
  RealType w_beta,w_alpha;
  RealType E_avg, V_avg;
  std::string GEVtype;
  bool forOpt;

  ///These are the values we collect to build the Matrices GLOBAL
  Matrix<RealType> Olp, Ham, Ham2;
  std::vector<RealType> D_E, HD2, HD, D;
  RealType sE,sE2,sE4,sW,sN;

  inline void clearComponentMatrices()
  {
    Olp=0.0;
    Ham=0.0;
    Ham2=0.0;
    for(int i=0; i<D_E.size(); i++)
    {
      D_E[i]=0.0;
      HD[i]=0.0;
      HD2[i]=0.0;
      D[i]=0.0;
    }
    sE=0;
    sE2=0;
    sE4=0;
    sW=0;
    sN=0;
  }

  void resizeForOpt(int n)
  {
    Olp.resize(n,n);
    Ham.resize(n,n);
    Ham2.resize(n,n);
    D_E.resize(n);
    HD.resize(n);
    HD2.resize(n);
    D.resize(n);
    clearComponentMatrices();
  }
};
}

#endif
