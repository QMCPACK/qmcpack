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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_VMCMULTIRENYI_OMP_H
#define QMCPLUSPLUS_VMCMULTIRENYI_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

struct EEholder
{
  vector<vector<int> > g_w_num;
  vector<vector<vector<int> > > g_w_pindex;

  int NG, NW; //groups walkers

  void resize(int nw, int ng)
  {
    NG=ng;
    NW=nw;
    g_w_num.resize(ng);
    g_w_pindex.resize(ng);
    for(int i(0); i<ng; i++)
    {
      g_w_num[i].resize(nw,0);
      g_w_pindex[i].resize(nw);
    }
  }
};

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCMultiRenyiOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  VMCMultiRenyiOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                   HamiltonianPool& hpool, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
  //inline vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  ///number of warmup steps
  int myWarmupSteps;
  ///number of RN warmup steps
  int myRNWarmupSteps;
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///option to enable/disable drift equation or RN for VMC
  string UseDrift;
  ///compute EE
  string computeEE, csEE;
  int EEN, EENR;
  RealType EEdR;
  RealType csoffset;
  RealType vsize, rsize;
  stringstream N_dat, s_dat, ee_dat;
  void godeeper(int lvl,vector<Matrix<RealType> >& ratios, vector<RealType>& estimateA,vector<RealType>& estimateV,vector<RealType>& estimateS,vector<RealType>& estimateN,vector<RealType>& tv,vector<int>& itr);
  fstream file_out;
  ///Ways to set rn constant
  RealType logoffset,logepsilon;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  VMCMultiRenyiOMP(const VMCMultiRenyiOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  VMCMultiRenyiOMP& operator=(const VMCMultiRenyiOMP&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile: VMCMultiRenyiOMP.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCMultiRenyiOMP.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
