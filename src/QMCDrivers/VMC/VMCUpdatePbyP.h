//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_VMC_PARTICLEBYPARTICLE_UPDATE_H
#define QMCPLUSPLUS_VMC_PARTICLEBYPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{


/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCUpdatePbyP: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdatePbyP();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//     void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i);

//     void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//     , std::vector<MCWalkerConfiguration*>& wclone
//     , std::vector<QMCHamiltonian*>& hclone
//     , std::vector<RandomGenerator_t*>& rng
//     , std::vector<RealType>& ratio_i_0);

  RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);
//     RealType advanceWalkerForCSEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios, std::vector<RealType>& weights, std::vector<RealType>& logs );

private:
  std::vector<NewTimer*> myTimers;
};

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
 */
class VMCUpdatePbyPWithDrift: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdatePbyPWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                         QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdatePbyPWithDrift();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  std::vector<NewTimer*> myTimers;
};

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
 */
class VMCUpdatePbyPWithDriftFast: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                             QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdatePbyPWithDriftFast();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
//     for linear opt CS
//     void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i);
  RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);
//     RealType advanceWalkerForCSEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios, std::vector<RealType>& weights, std::vector<RealType>& logs );
private:
  std::vector<NewTimer*> myTimers;
};

/** @ingroup QMCDrivers  ParticleByParticle
*@brief Implements the VMC algorithm using particle-by-particle move with the drift equation.
*/
class VMCUpdateRenyiWithDriftFast: public QMCUpdateBase
{
public:
  /// Constructor.
  VMCUpdateRenyiWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                              QMCHamiltonian& h, RandomGenerator_t& rg);

  ~VMCUpdateRenyiWithDriftFast();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
  std::vector<NewTimer*> myTimers;
};


/** @ingroup QMCDrivers  ParticleByParticle
*@brief Implements the VMC algorithm using particle-by-particle move. Samples |Psi| to increase number of walkers near nodes.
*/
//   class VMCUpdatePbyPSampleRN: public QMCUpdateBase {
//   public:
//     /// Constructor.
//     VMCUpdatePbyPSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide,
//         QMCHamiltonian& h, RandomGenerator_t& rg);
//
//     ~VMCUpdatePbyPSampleRN();
//
//     void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);
// //     void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i);
//     void setLogEpsilon(RealType eps)
//     {
//       logEpsilon=eps;
// //       app_log()<<eps<< std::endl;
//     }
//     void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);
//
// //     void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
// //     , std::vector<MCWalkerConfiguration*>& wclone
// //     , std::vector<QMCHamiltonian*>& hclone
// //     , std::vector<RandomGenerator_t*>& rng
// //     , std::vector<RealType>& ratio_i_0);
//
//   private:
//     std::vector<NewTimer*> myTimers;
//     //prefactor multiplying the guiding function
//     RealType logEpsilon;
//   };
}

#endif
