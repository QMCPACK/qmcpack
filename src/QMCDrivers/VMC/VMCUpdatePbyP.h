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

class Communicate;

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
 *  @brief   Implements the VMC algorithm using particle-by-particle moves with a nodeless guiding function
 */
class VMCUpdatePbyPNodeless : public QMCUpdateBase {

  public:

    VMCUpdatePbyPNodeless(MCWalkerConfiguration & w,
                          TrialWaveFunction & psi,
                          QMCHamiltonian & h,
                          RandomGenerator_t & rg,
                          const ParticleSet & ips,
                          const RealType eps);

    ~VMCUpdatePbyPNodeless();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

    RealType advanceWalkerForEE(Walker_t& w1,
                                std::vector<PosType>& dR,
                                std::vector<int>& iats,
                                std::vector<int>& rs,
                                std::vector<RealType>& ratios);

    bool put(xmlNodePtr cur);

    static void reset_history(Communicate * const comm, const int nthread);
    static void process_history(Communicate * const comm, const int nthread, const bool record);
    static void reset_tfl();

    template <typename T>
    static void check_all_positive(const T & container, const std::string & name) {
      for (const auto & d : container)
        if ( d <= 0 )
          APP_ABORT( ( "ERROR: an element of " + name + " was not positive" ) );
    }

  private:

    RealType get_nodeless_gf_sqn(const RealType logpsi,
                                 const ParticleSet::ParticleGradient_t & grad,
                                 const ParticleSet::ParticleLaplacian_t & lap) const;

    RealType get_nodeless_gf_sqn_new(const RealType logpsi) const;

    RealType init_nodeless(const ParticleSet & P, const RealType tfl);

    RealType update_nodeless(const ParticleSet & P, const int iat, const RealType tfl);

    void nodeless_accept(const ParticleSet & P, const int iat, const RealType tfl);

    /// \brief  whether we need to re-initialize the nodeless guiding function
    bool nodelessInitialized;

    /// \brief  coefficient for the laplacian in the guiding function formula
    RealType NodelessEpsilon;

    /// \brief  coefficient for the exponential damping in the guiding function formula
    RealType NodelessAlpha;
    RealType NodelessBeta;

    /// \brief  saved value of the nodeless guiding function
    RealType savedGF;
    RealType savedGFTmp;

    /// \brief  overall counting group penalty's exponent to the nodeless adjustment
    RealType cgPenaltyExponent;
    RealType cgPenaltyExponentTmp;

    /// \brief  overall min distance penalty to the nodeless adjustment
    RealType mdPenalty;
    RealType mdPenaltyTmp;

    /// \brief  min distance penalties for each particle
    std::vector<RealType> mdPenalties;
    RealType mdPenaltiesTmp;

    /// \brief  particle number count for each counting region
    std::vector<RealType> cgCounts;
    std::vector<RealType> cgCountsTmp;

    /// \brief  matrix of un-normalized counting function values at each particle position
    std::vector<RealType> cgUnormalized;
    std::vector<RealType> cgUnormalizedTmp;

    /// \brief  vector of counting normalizations
    std::vector<RealType> cgNorms;
    RealType cgNormTmp;

    /// \brief  standard deviations for the counting groups' count penalty
    std::vector<RealType> cgCountSigmas;

    /// \brief  target number of electrons for each counting group
    std::vector<RealType> cgCountNelecs;

    /// \brief  starting index for looking up each counting group's guassians
    std::vector<int> cgGaussStarts;

    /// \brief  ending index for looking up each counting group's guassians
    std::vector<int> cgGaussEnds;

    /// \brief  coefficients for counting groups' linear combinations of gaussians
    std::vector<RealType> cgGaussAlphas;

    /// \brief  standard deviations for the counting groups' gaussians
    std::vector<RealType> cgGaussSigmas;

    /// \brief  switching speeds for the min distance cutoffs
    std::vector<RealType> mdBetas;

    /// \brief  switching distances for the min distance cutoffs
    std::vector<RealType> mdDists;

    /// \brief  coordinates of the counting groups' gaussians' centers
    std::vector<ParticleSet::SingleParticlePos_t> cgGaussCenters;

    /// \brief  coordinates of the min distance cutoff centers
    std::vector<ParticleSet::SingleParticlePos_t> mdCenters;

    /// \brief  vector containing the positions of the ions (i.e. the atomic nuclei)
    std::vector<ParticleSet::SingleParticlePos_t> IonPositions;

    /// \brief  history of the sampled configurations' trial function logarithms
    static std::vector<std::vector<RealType> > tfl_history;

    /// \brief  average of trial function logarithms
    static RealType tfl_avg;

    /// \brief  standard deviation of trial function logarithms
    static RealType tfl_sdv;

//    /// \brief  whether to record tfl_avg and tfl_sdv from the next sample
//    static bool setNGMag;

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
