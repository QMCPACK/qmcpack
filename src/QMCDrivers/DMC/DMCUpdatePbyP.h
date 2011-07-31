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
#ifndef QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus {

  class DMCUpdatePbyPWithRejection: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdatePbyPWithRejection(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdatePbyPWithRejection();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:
    vector<NewTimer*> myTimers;
  };

  class DMCUpdatePbyPWithRejectionFast: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdatePbyPWithRejectionFast();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:
    vector<NewTimer*> myTimers;
  };


  class DMCUpdatePbyPWithKill: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdatePbyPWithKill(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdatePbyPWithKill();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:

    /// Copy Constructor (disabled)
    DMCUpdatePbyPWithKill(const DMCUpdatePbyPWithKill& a): QMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCUpdatePbyPWithKill& operator=(const DMCUpdatePbyPWithKill&) { return *this;}
    vector<NewTimer*> myTimers;

  };
  
  class RNDMCUpdatePbyPFast: public QMCUpdateBase {

  public:

    RNDMCUpdatePbyPFast(MCWalkerConfiguration& w, MCWalkerConfiguration& wg, TrialWaveFunction& psi, TrialWaveFunction& guide, 
        QMCHamiltonian& h, RandomGenerator_t& rg);

    ~RNDMCUpdatePbyPFast();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
    
    void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

//     void estimateNormWalkers(vector<TrialWaveFunction*>& pclone
//     , vector<MCWalkerConfiguration*>& wclone
//     , vector<QMCHamiltonian*>& hclone
//     , vector<RandomGenerator_t*>& rng
//     , vector<RealType>& ratio_i_0);

  private:
    MCWalkerConfiguration W_G;
    vector<NewTimer*> myTimers;
    int maxS;
    RealType efn;
    int estimateCrossings, maxcopy;
    
  };
  
  
  class RNDMCUpdatePbyPCeperley: public QMCUpdateBase {
    
    public:
      
      /// Constructor.
      RNDMCUpdatePbyPCeperley(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
                              QMCHamiltonian& h, RandomGenerator_t& rg);
                                       ///destructor
      ~RNDMCUpdatePbyPCeperley();
                                       
      void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
      void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);
      
    private:
      vector<NewTimer*> myTimers;
      int maxS;
      RealType efn;
      int estimateCrossings;
  };  
}

#endif
/***************************************************************************
 * $RCSfile: DMCUpdatePbyP.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
