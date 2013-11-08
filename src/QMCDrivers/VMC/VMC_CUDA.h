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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  string UseDrift;
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
  bool checkBounds (vector<PosType> &newpos, vector<bool> &valid);

  void resetRun();

  opt_variables_type dummy;
  int numParams;
  Matrix<ValueType> d_logpsi_dalpha, d_hpsioverpsi_dalpha;
  RealType w_beta,w_alpha;
  RealType E_avg, V_avg;
  string GEVtype;
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
/***************************************************************************
 * $RCSfile: VMCcuda.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCcuda.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
