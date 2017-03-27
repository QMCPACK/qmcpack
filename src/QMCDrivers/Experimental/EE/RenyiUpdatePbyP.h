//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_RENYI_PARTICLEBYPARTICLE_UPDATE_H
#define QMCPLUSPLUS_RENYI_PARTICLEBYPARTICLE_UPDATE_H
#include "QMCDrivers/EE/QMCRenyiUpdateBase.h"

namespace qmcplusplus
{

class RenyiUpdatePbyP: public QMCRenyiUpdateBase
{
public:
  /// Constructor.
  RenyiUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                  QMCHamiltonian& h, RandomGenerator_t& rg, int order);

  ~RenyiUpdatePbyP();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
  void plotSwapAmplitude(WalkerIter_t it, WalkerIter_t it_end,Matrix<RealType>& averageSwaps);

private:
};

//   class RenyiUpdateAll: public QMCRenyiUpdateBase {
//   public:
//     /// Constructor.
//     RenyiUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//         QMCHamiltonian& h, RandomGenerator_t& rg, int order);
//
//     ~RenyiUpdateAll();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//
// //     RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);
//
//   private:
//   };
//
//   class RenyiUpdateAllWithDrift: public QMCRenyiUpdateBase {
//   public:
//     /// Constructor.
//     RenyiUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//         QMCHamiltonian& h, RandomGenerator_t& rg, int order);
//
//     ~RenyiUpdateAllWithDrift();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//
// //     RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);
//
//   private:
//     std::vector<ParticleSet::ParticlePos_t*> drifts;
//   };

//   class VMCUpdatePbyPWithDrift: public QMCUpdateBase {
//   public:
//     /// Constructor.
//     VMCUpdatePbyPWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//         QMCHamiltonian& h, RandomGenerator_t& rg);
//
//     ~VMCUpdatePbyPWithDrift();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//
//   private:
//     std::vector<NewTimer*> myTimers;
//   };

//   class VMCUpdatePbyPWithDriftFast: public QMCUpdateBase {
//   public:
//     /// Constructor.
//     VMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//         QMCHamiltonian& h, RandomGenerator_t& rg);
//
//     ~VMCUpdatePbyPWithDriftFast();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
// //     for linear opt CS
//     void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i);
//     RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);
//     RealType advanceWalkerForCSEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios, std::vector<RealType>& weights, std::vector<RealType>& logs );
//   private:
//     std::vector<NewTimer*> myTimers;
//   };
//
//     /** @ingroup QMCDrivers  ParticleByParticle
//    *@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
//    */
//   class VMCUpdateRenyiWithDriftFast: public QMCUpdateBase {
//   public:
//     /// Constructor.
//     VMCUpdateRenyiWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//         QMCHamiltonian& h, RandomGenerator_t& rg);
//
//     ~VMCUpdateRenyiWithDriftFast();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//
//   private:
//     std::vector<NewTimer*> myTimers;
//   };


//   class VMCUpdatePbyPSampleRN: public QMCUpdateBase {
//   public:
//     /// Constructor.
//     VMCUpdatePbyPSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide,
//         QMCHamiltonian& h, RandomGenerator_t& rg);
//
//     ~VMCUpdatePbyPSampleRN();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//     void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i);
//     void setLogEpsilon(RealType eps)
//     {
//       logEpsilon=eps;
// //       app_log()<<eps<< std::endl;
//     }
//     void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);
//
//     void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//     , std::vector<MCWalkerConfiguration*>& wclone
//     , std::vector<QMCHamiltonian*>& hclone
//     , std::vector<RandomGenerator_t*>& rng
//     , std::vector<RealType>& ratio_i_0);
//
//   private:
//     std::vector<NewTimer*> myTimers;
//     //prefactor multiplying the guiding function
//     RealType logEpsilon;
//   };
}

#endif
