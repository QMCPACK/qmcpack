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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file DMCOMPOPT.h
 * @brief Define OpenMP-able DMC Driver.
 */
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/CommOperators.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers
 * @brief DMC driver using OpenMP paragra
 *
 * This is the main DMC driver with MPI/OpenMP loops over the walkers.
 */
class DMCOMPOPT: public QMCDriver, public CloneManager
{
public:

  /// Constructor.
  DMCOMPOPT(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
            HamiltonianPool& hpool,WaveFunctionPool& ppool);

  bool run();
  bool put(xmlNodePtr cur)
  {
    return true;
  };
  void setTau(RealType i);
//     void resetComponents(xmlNodePtr cur);
  void fillVectors(std::vector<RealType>& d, std::vector<RealType>& hd, RealType& e, Matrix<RealType>& olp)
  {
    myComm->allreduce(D);
    myComm->allreduce(HD);
    myComm->allreduce(DT);
    myComm->allreduce(D2);
    myComm->allreduce(D_E);
    myComm->allreduce(Overlap);
    std::vector<RealType> g_stats(6,0);
    g_stats[0]=sE;
    g_stats[1]=sE2;
    g_stats[2]=ssE;
    g_stats[3]=sW;
    g_stats[4]=smW;
    g_stats[5]=sN;
    myComm->allreduce(g_stats);
    RealType nrm = 1.0/g_stats[5];
    RealType wgtNrm = std::exp(t*g_stats[2]*nrm);
    RealType avgInv =std::sqrt(g_stats[3]/g_stats[5]);
    RealType avgInvWgt = std::sqrt((g_stats[3]/g_stats[5])/wgtNrm);
    RealType avgWgt = std::sqrt(g_stats[4]/g_stats[5]);
    if ((printderivs=="yes")&&(myComm->rank()==0))
    {
      stringstream fn;
      fn<<RootName.c_str()<<".derivs";
      ofstream d_out(fn.str().c_str());
      d_out.precision(6);
      d_out<<"#CSF    D        HD"<<endl;
      for (int i=0; i<NumOptimizables; i++)
        d_out<<i<<"  "<<nrm*D[i]<<"  "<<nrm*HD[i]<<endl;
    }
//       app_log()<<"Ebar: "<<g_stats[2]/g_stats[5]<<endl;
//       app_log()<<"Weights: "<<avgWgt<<" "<<avgInvWgt<<" "<<(g_stats[4]/g_stats[5])*wgtNrm<<" "<<(g_stats[3]/g_stats[5])/wgtNrm<<" "<<wgtNrm<<endl;
    std::vector<RealType> DT2(DT.size(),0);
    for(int i=0; i<NumOptimizables; i++)
      DT2[i] = (Overlap(i,i)==0)?0:1.0/std::sqrt(nrm*Overlap(i,i));
    for(int i=0; i<NumOptimizables; i++)
      for(int j=0; j<NumOptimizables; j++)
      {
        Overlap(i,j) *= nrm;//DT2[i]*DT2[j];
        olp(i+1,j+1) = Overlap(i,j);
      }
    for (int i=0; i<NumOptimizables; i++)
      olp(i+1,0) = olp(0,i+1)= 0.0;
    olp(0,0)=1.0;
//       RealType dpNrm = std::sqrt(4.0/(avgInvWgt + avgWgt + 2.0))*nrm;
    for(int i=0; i<D.size(); i++)
      d[i] = DT2[i]*(D[i]/avgWgt - DT[i]/avgInv)*nrm;
    for(int i=0; i<D.size(); i++)
      hd[i] = DT2[i]*(D[i]/avgWgt)*nrm;
//       std::vector<RealType> d2(d.size(),0);
//       for(int i=0;i<D.size();i++) d2[i] = DT2[i]*(D[i] - DT[i]/avgInvWgt);
    e=E_avg;
//       RealType Det= invert_matrix(Overlap,true);
//       app_log()<<Det<<endl;
//
//       for (int i=0; i<NumOptimizables; i++) app_log()<< DT2[i] <<" ";
//       app_log()<<endl;
    for (int i=0; i<NumOptimizables; i++)
      app_log()<< DT2[i]*(D[i]/avgWgt)*nrm<<" ";
    app_log()<<endl;
    for (int i=0; i<NumOptimizables; i++)
      app_log()<< DT2[i]*(DT[i]/avgInv)*nrm<<" ";
    app_log()<<endl;
//
    for (int i=0; i<NumOptimizables; i++)
      app_log()<<d[i]<<" ";
    app_log()<<endl;
//
//       for (int i=0; i<NumOptimizables; i++)
//       {
//         d[i]=0.0;
//         for (int j=0; j<NumOptimizables; j++)
//           d[i] += d0[j]*Overlap(i,j);
//       }
//
//       for (int i=0; i<NumOptimizables; i++) app_log()<<d[i]<<" ";
//       app_log()<<endl;
  }

private:
  ///Index to determine what to do when node crossing is detected
  IndexType KillNodeCrossing;
  ///Interval between branching
  IndexType BranchInterval;
  ///hdf5 file name for Branch conditions
  string BranchInfo;
  ///input string to determine kill walkers or not
  string KillWalker;
  ///input string to determine swap walkers among mpi processors
  string SwapWalkers;
  ///input string to determine to use reconfiguration
  string Reconfiguration;
  ///input string to determine to use nonlocal move
  string NonLocalMove;
  ///input string to benchmark OMP performance
  string BenchMarkRun;
  ///input string to use fast gradient
  string UseFastGrad;
  ///input to control maximum age allowed for walkers.
  IndexType mover_MaxAge;

  int NumOptimizables;
  RealType E_avg, V_avg, t;
  std::vector<RealType> D_E, D2, D, DT, HD;
  Matrix<RealType> Overlap;
  RealType sE,sE2,ssE,sN;
  RealType sW,smW;
  int myPeriod4WalkerDump, wlen, Eindx;
  int samples_this_node;
  bool firsttime;
  string printderivs;
  std::vector<opt_variables_type> dummyOptVars;


  void resetUpdateEngines();
  void benchMark();
  /// Copy Constructor (disabled)
  DMCOMPOPT(const DMCOMPOPT& a): QMCDriver(a), CloneManager(a) { }
  /// Copy operator (disabled).
  DMCOMPOPT& operator=(const DMCOMPOPT&)
  {
    return *this;
  }

  inline void clearComponentMatrices()
  {
    Overlap=0.0;
    for(int i=0; i<D_E.size(); i++)
    {
      D_E[i]=0.0;
      D[i]=0.0;
      HD[i]=0.0;
      DT[i]=0.0;
      D2[i]=0.0;
    }
    sE=0;
    sE2=0;
    ssE=0;
    sW=0;
    smW=0;
    sN=0;
  }

  void resizeForOpt(int n)
  {
    Overlap.resize(n,n);
    D_E.resize(n);
    D.resize(n);
    HD.resize(n);
    DT.resize(n);
    D2.resize(n);
    clearComponentMatrices();
  }

  void fillComponentMatrices(std::vector<RealType>& d, std::vector<RealType>& hd, Walker_t& w)
  {
    std::vector<RealType> g_stats(5,0);
    RealType E_L = w.getPropertyBase()[LOCALENERGY];
    RealType E_L2= E_L*E_L;
    sE +=E_L;
    sE2+=E_L2;
    sN+=1;
    ssE+=w.getPropertyHistorySum(Eindx,wlen);
    RealType invfnw = std::exp(t*w.getPropertyHistorySum(Eindx,wlen));
    sW +=invfnw;
    smW+=1.0/invfnw;
    for (int i=0; i<NumOptimizables; i++)
    {
      RealType di  = d[i];
      D_E[i]+=   invfnw*di*E_L;
      DT[i]  +=   di*invfnw;
      D[i]  +=   di;
      HD[i]  +=   hd[i];
      D2[i] +=   di*di;
      for (int j=0; j<NumOptimizables; j++)
        Overlap(i,j) += di*d[j]*invfnw;
    }
  }
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCOMPOPT.h,v $   $Author: jnkim $
 * $Revision: 1604 $   $Date: 2007-01-06 11:00:24 -0600 (Sat, 06 Jan 2007) $
 * $Id: DMCOMPOPT.h 1604 2007-01-06 17:00:24Z jnkim $
 ***************************************************************************/
