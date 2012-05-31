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
#ifndef QMCPLUSPLUS_QMCRENYIUPDATEBASE_H
#define QMCPLUSPLUS_QMCRENYIUPDATEBASE_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{
  class QMCRenyiUpdateBase: public QMCUpdateBase {
  public:
    /// Constructor.
    QMCRenyiUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg, int order=2);
        
    ~QMCRenyiUpdateBase();
    
    void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);
    void initWalkers(WalkerIter_t it, WalkerIter_t it_end);
    
    virtual void check_region(WalkerIter_t it, WalkerIter_t it_end, RealType v, string shape, ParticleSet::ParticlePos_t& ed, ParticleSet::ParticlePos_t& Center, int maxN, int minN);
    virtual int get_region(ParticleSet::ParticlePos_t& Pos,int iat);
    
    virtual void double_check_region(WalkerIter_t it, WalkerIter_t it_end);
    
    void put_in_box(PosType& Pos);
    RealType get_stats(std::vector<RealType>& n)
    {
      RealType nrm=1.0/cnt;
      for (int i(0);i<n.size();i++) n[i]+=nrm*n_region[i];
      return nrm*regions[NumPtcl+3];
    }

    void clear_stats()
    {
      std::fill(n_region.begin(),n_region.end(),0);
      cnt=0;
      regions[NumPtcl+3]=0;
    }
    
    void print_all();


//     RealType advanceWalkerForEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios);

    vector<NewTimer*> myTimers;
    vector<MCWalkerConfiguration*> W_vec;
    vector<TrialWaveFunction*> Psi_vec;
    vector<ParticleSet::ParticlePos_t*> deltaR_vec;
    
    
    int RenyiOrder;
    int mxN,mnN;
    int cnt;
    ///patch geometry (sphere,...
    string computeEE;
    ///patch size (sphere:r^2,...
    RealType vsize;
    ParticleSet::ParticlePos_t C,Edge;
    std::vector<RealType> regions, n_region;
      
  };
  
}

#endif
/***************************************************************************
 * $RCSfile: VMCRenyiOMP.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCRenyiOMP.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
