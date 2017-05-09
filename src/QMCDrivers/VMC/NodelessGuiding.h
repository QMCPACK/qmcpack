//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 QMCPACK developers.
//
// File developed by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California Berkeley
//
// File created by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VMC_NODELESSGUIDING_H
#define QMCPLUSPLUS_VMC_NODELESSGUIDING_H

#include <vector>
#include <string>

#include "QMCDrivers/QMCUpdateBase.h"

class Communicate;

namespace qmcplusplus
{

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Class that controls sampling via a nodeless guiding function
///
///////////////////////////////////////////////////////////////////////////////////////////////////
class VMCUpdatePbyPNodeless : public QMCUpdateBase
{

  // public member functions
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

  // public member functions (static)
  public:

    static void reset_history(Communicate * const comm, const int nthread);

    static void process_history(Communicate * const comm, const int nthread, const bool record);

    static void reset_tfl();

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Checks that all elements of a container are positive.
    ///
    /// \param[in]      container   the container to check
    /// \param[in]      name        name of the container for use in error printing
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    static void check_all_positive(const T & container, const std::string & name)
    {
      for (const auto & d : container)
        if ( d <= 0 )
          APP_ABORT( ( "ERROR: an element of " + name + " was not positive" ) );
    }


  // private member functions
  private:

    RealType get_nodeless_gf_sqn(const RealType logpsi,
                                 const ParticleSet::ParticleGradient_t & grad,
                                 const ParticleSet::ParticleLaplacian_t & lap) const;

    RealType get_nodeless_gf_sqn_new(const RealType logpsi) const;

    RealType init_nodeless(const ParticleSet & P, const RealType tfl);

    RealType update_nodeless(const ParticleSet & P, const int iat, const RealType tfl);

    void nodeless_accept(const ParticleSet & P, const int iat, const RealType tfl);

  // private member data
  private:

    /// \brief  whether we need to re-initialize the nodeless guiding function
    bool nodelessInitialized;

    /// \brief  coefficient determining the strength of the nodeless adjustment
    RealType NodelessEpsilon;

    /// \brief  coefficient for the old guiding function formula
    RealType NodelessAlpha;

    /// \brief  coefficient for the old guiding function formula
    RealType NodelessBeta;

    /// \brief  saved value of the nodeless guiding function
    RealType savedGF;

    /// \brief  temporary value of the nodeless guiding function
    RealType savedGFTmp;

    /// \brief  overall counting group penalty's exponent to the nodeless adjustment
    RealType cgPenaltyExponent;

    /// \brief  temporary overall counting group penalty's exponent to the nodeless adjustment
    RealType cgPenaltyExponentTmp;

    /// \brief  overall min distance penalty to the nodeless adjustment
    RealType mdPenalty;

    /// \brief  temporary overall min distance penalty to the nodeless adjustment
    RealType mdPenaltyTmp;

    /// \brief  min distance penalties for each particle
    std::vector<RealType> mdPenalties;

    /// \brief  temporary min distance penalty
    RealType mdPenaltiesTmp;

    /// \brief  particle number count for each counting region
    std::vector<RealType> cgCounts;

    /// \brief  temporary particle number count for each counting region
    std::vector<RealType> cgCountsTmp;

    /// \brief  matrix of un-normalized counting function values at each particle position
    std::vector<RealType> cgUnormalized;

    /// \brief  temporary matrix of un-normalized counting function values at each particle position
    std::vector<RealType> cgUnormalizedTmp;

    /// \brief  vector of counting normalizations
    std::vector<RealType> cgNorms;

    /// \brief  a temporary counting group normalization
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

  // private member data (static)
  private:

    /// \brief  running totals for count averages
    static std::vector<std::vector<RealType> > cga_hist;

    /// \brief  running totals for count variances
    static std::vector<std::vector<RealType> > cgv_hist;

    /// \brief  history of the sampled configurations' trial function logarithms
    static std::vector<RealType> tla_hist;

    /// \brief  history of the square of the sampled configurations' trial function logarithms
    static std::vector<RealType> tlv_hist;

    /// \brief  running totals for number of samples;
    static std::vector<RealType> nsp_hist;

    /// \brief  average of trial function logarithms
    static RealType tfl_avg;

    /// \brief  standard deviation of trial function logarithms
    static RealType tfl_sdv;

    /// \brief  whether or not nodeless guiding is being used
    static bool usingNodelessGuiding;

};

}

#endif
