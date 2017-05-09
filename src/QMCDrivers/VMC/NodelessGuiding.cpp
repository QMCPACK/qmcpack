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

#include <cmath>
#include <algorithm>
#include <boost/format.hpp>

#include "QMCDrivers/VMC/NodelessGuiding.h"
#include "Numerics/Blasf.h"
#include "OhmmsData/XMLcxx11Helper.h"
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperatorsMPI.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{

// initialize static data members
std::vector<std::vector<VMCUpdatePbyPNodeless::RealType> > VMCUpdatePbyPNodeless::cga_hist;
std::vector<std::vector<VMCUpdatePbyPNodeless::RealType> > VMCUpdatePbyPNodeless::cgv_hist;
std::vector<VMCUpdatePbyPNodeless::RealType>               VMCUpdatePbyPNodeless::tla_hist;
std::vector<VMCUpdatePbyPNodeless::RealType>               VMCUpdatePbyPNodeless::tlv_hist;
std::vector<VMCUpdatePbyPNodeless::RealType>               VMCUpdatePbyPNodeless::nsp_hist;
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::tfl_avg = 0.0;
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::tfl_sdv = -1.0;
bool VMCUpdatePbyPNodeless::usingNodelessGuiding = false;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Constructor for the VMCUpdate class that uses particle-by-particle moves and a
///         nodeless guiding function
///
/// \param[in,out]  w      the set of walkers to use
/// \param[in,out]  psi    the trial function, which in general has nodes
/// \param[in,out]  h      the Hamiltonian
/// \param[in,out]  rg     random number generator
/// \param[in]      ips    ion particle set (tells us where the ions are positioned)
/// \param[in]      eps    coefficient for nodeless guiding adjustment
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless(MCWalkerConfiguration & w,
                                             TrialWaveFunction & psi,
                                             QMCHamiltonian & h,
                                             RandomGenerator_t & rg,
                                             const ParticleSet & ips,
                                             const RealType eps)
: QMCUpdateBase(w,psi,h,rg)
, NodelessEpsilon(eps)
, IonPositions(ips.R.begin(), ips.R.end())
, nodelessInitialized(false)
{

  // set the flag to ensure we are using nodeless guiding
  #pragma omp single
  {
    usingNodelessGuiding = true;
  }

  // ensure some sanity
  if ( NodelessEpsilon < RealType(0) )
    APP_ABORT("VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless was given a negative epsilon");
  //if ( NodelessAlpha < RealType(0) )
  //  APP_ABORT("VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless was given a negative alpha");
  //if ( NodelessBeta < RealType(0) )
  //  APP_ABORT("VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless was given a negative beta");

  //// print the ion positions
  //std::cout << std::endl;
  //std::cout << "printing ion positions in VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless" << std::endl;
  //for (auto it = IonPositions.begin(); it != IonPositions.end(); it++)
  //  std::cout << *it << std::endl;
  //std::cout << std::endl;
  //APP_ABORT("Stopping after printing ion particle positions for Nodeless guiding");

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  destructor
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::~VMCUpdatePbyPNodeless()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Read parameters of the nodeless guiding function from xml
///
/// \param[in,out]  cur      pointer to the xml node to read from
///
///////////////////////////////////////////////////////////////////////////////////////////////////
bool VMCUpdatePbyPNodeless::put(xmlNodePtr cur)
{

  // do the usual base class put
  const bool base_success = QMCUpdateBase::put(cur);

  // we will no longer be initialized after this
  nodelessInitialized = false;

  // find and process nodeless guiding node
  bool my_success = false;
  for ( auto kur : getXMLChildList(cur, "nodelessGuiding") ) {

    // remember that we found something
    my_success = true;

    // reset guiding function parameters
    cgCountSigmas.clear();
    cgCountNelecs.clear();
    cgGaussStarts.clear();
    cgGaussEnds.clear();
    cgGaussAlphas.clear();
    cgGaussSigmas.clear();
    cgGaussCenters.clear();
    mdBetas.clear();
    mdDists.clear();
    mdCenters.clear();

    // read counting groups' info
    for ( auto chl : getXMLChildList(kur, "countGroup") ) {

      // read counting group standard deviation and target number of electrons
      cgCountSigmas.push_back(-1.0);
      cgCountNelecs.push_back(0.0);
      getXMLAttributes(chl, *cgCountSigmas.rbegin(), "sigma", *cgCountNelecs.rbegin(), "nelec");

      // record the start of this counting group's gaussians
      cgGaussStarts.push_back(cgGaussSigmas.size());

      // read the standard deviations and centers of the gaussians
      for ( auto khl : getXMLChildList(chl, "gaussian") ) {

        // <gaussian alpha="1.0" sigma="1.3" type="Array"> 0.1 0.2 0.3 </gaussian>

        // read alpha and sigma values for this gaussian
        cgGaussAlphas.push_back(-1.0);
        cgGaussSigmas.push_back(-1.0);
        getXMLAttributes(khl, *cgGaussAlphas.rbegin(), "alpha", *cgGaussSigmas.rbegin(), "sigma");

        // read the coordinates of the gaussian's center
        cgGaussCenters.push_back(ParticleSet::SingleParticlePos_t());
        if ( ! putContent(cgGaussCenters.rbegin()->begin(), cgGaussCenters.rbegin()->end(), khl) )
          APP_ABORT("ERROR: problem reading coordinates of a countGroup gaussian center");

      }

      // record the end of this counting group's gaussians
      cgGaussEnds.push_back(cgGaussSigmas.size());

    }

    // read minimum distance information
    for ( auto chl : getXMLChildList(kur, "minDistCenter") ) {

      // <minDistCenter beta="1.0" dist="1.3"> -3.1 3.3 4.7 </minDistCenter>

      // read switching speed and distance parameters
      mdBetas.push_back(-1.0);
      mdDists.push_back(-1.0);
      getXMLAttributes(chl, *mdBetas.rbegin(), "beta", *mdDists.rbegin(), "dist");

      // read the coordinates of the center to measure distance from
      mdCenters.push_back(ParticleSet::SingleParticlePos_t());
      if ( ! putContent(mdCenters.rbegin()->begin(), mdCenters.rbegin()->end(), chl) )
        APP_ABORT("ERROR: problem reading coordinates of a min distance center");

    }

    // check that things that should be positive are
    check_all_positive(cgCountSigmas, "cgCountSigmas");
    check_all_positive(cgGaussAlphas, "cgGaussAlphas");
    check_all_positive(cgGaussSigmas, "cgGaussSigmas");
    check_all_positive(mdBetas, "mdBetas");
    check_all_positive(mdDists, "mdDists");

    // check that there is at least one minimum distance center
    if ( mdBetas.empty() )
      APP_ABORT("ERROR: Found no minimum distance centers while reading nodelessGuiding entry in xml input");

    // print what we have
    const bool to_print = false;
    if ( to_print && OHMMS::Controller->rank() == 0 && omp_get_thread_num() == 0 ) {
      app_log() << std::endl;
      app_log() << omp_get_num_threads() << " threads are running through VMCUpdatePbyPNodeless::put" << std::endl;
      app_log() << std::endl;
      app_log() << "NodelessEpsilon = " << NodelessEpsilon << std::endl;
      app_log() << std::endl;
      for (int i = 0; i < cgCountSigmas.size(); i++) {
        app_log() << "Count group:" << std::endl;
        app_log() << "  count sigma = " << cgCountSigmas.at(i) << std::endl;
        app_log() << "  count nelec = " << cgCountNelecs.at(i) << std::endl;
        app_log() << "  Gaussians" << std::endl;
        for (int j = cgGaussStarts.at(i); j < cgGaussEnds.at(i); j++) {
          app_log() << "    alpha = " << cgGaussAlphas.at(j) << std::endl;
          app_log() << "    sigma = " << cgGaussSigmas.at(j) << std::endl;
          app_log() << "   center = " << cgGaussCenters.at(j) << std::endl;
        }
        app_log() << std::endl;
      }
      for (int i = 0; i < mdBetas.size(); i++) {
        app_log() << "Min Dist Center:" << std::endl;
        app_log() << "    beta = " << mdBetas.at(i) << std::endl;
        app_log() << "    dist = " << mdDists.at(i) << std::endl;
        app_log() << "  center = " << mdCenters.at(i) << std::endl;
        app_log() << std::endl;
      }
    }

  }

  // make sure we found something
  if (!my_success)
    APP_ABORT("ERROR: failed to find nodelessGuiding entry in xml input");

  // return whether success was had
  return base_success && my_success;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Empties the history vectors and ensures they are the correct size.
///
/// \param[in,out]  comm     communicator for averaging over all processes
/// \param[in]      nthread  number of threads per process
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::reset_history(Communicate * const comm, const int nthread)
{

  // nothing to do if we are not using nodeless guiding
  if ( !usingNodelessGuiding )
    return;

  // reset running totals
  #pragma omp single
  {
    cga_hist.assign(nthread, std::vector<RealType>());
    cgv_hist.assign(nthread, std::vector<RealType>());
    tla_hist.assign(nthread, 0.0);
    tlv_hist.assign(nthread, 0.0);
    nsp_hist.assign(nthread, 0.0);
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes and stores the average and standard deviations of the history of trial
///         funciton logarithms and then clears the history.
///
/// \param[in,out]  comm     communicator for averaging over all processes
/// \param[in]      nthread  number of threads per process
/// \param[in]      record   whether the function should record tfl data
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::process_history(Communicate * const comm, const int nthread, const bool record)
{

  // nothing to do if we are not using nodeless guiding
  if ( !usingNodelessGuiding )
    return;

  // collect totals
  const int nc = ( cga_hist.at(0).empty() ? 0 : cga_hist.at(0).size() );
  std::vector<RealType> sums(2*nc+3, 0.0);
  RealType & tla = sums.at(2*nc+0);
  RealType & tlv = sums.at(2*nc+1);
  RealType & nsp = *sums.rbegin();
  for (int ip = 0; ip < cga_hist.size(); ip++) // loop over different threads' data
  {
    if ( nc ) daxpy(nc, 1.0, &cga_hist.at(ip).at(0), 1, &sums.at(0*nc), 1); // averages
    if ( nc ) daxpy(nc, 1.0, &cgv_hist.at(ip).at(0), 1, &sums.at(1*nc), 1); // averages of squares
    tla += tla_hist.at(ip); // trial logarithm average
    tlv += tlv_hist.at(ip); // trial logarithm average of squares
    nsp += nsp_hist.at(ip); // number of samples
  }
  comm->allreduce(sums);
  dscal(sums.size()-1, 1.0 / nsp, &sums.at(0), 1);

  if ( nsp > 0.1 ) // do nothing if we have no previous sample data to work with
  {

    app_log() << std::endl
              << " === Printing Statistics Related to Nodeless Guiding ===" << std::endl
              << std::endl;

    if ( record && tfl_sdv < 0.0 )
      app_log() << std::endl
                << " Note: these statistics are warmup statistics" << std::endl
                << std::endl;

    // compute variance of trial function logarithm
    tlv = std::sqrt( std::abs( tlv - tla * tla ) );

    // print trial function logarithm stats
    app_log() << std::endl;
    app_log() << "  Statistics of |log(psi)|: " << std::endl;
    app_log() << std::endl;
    app_log() << boost::format("               samples = %11.0f") % nsp << std::endl;
    app_log() << boost::format("               average = %18.6f") % tla << std::endl;
    app_log() << boost::format("    standard deviation = %18.6f") % tlv << std::endl;
    app_log() << std::endl;

    // print counting group stats
    if ( nc ) app_log() << "  Statistics of counting groups: " << std::endl;
    if ( nc ) app_log() << std::endl;
    for (int i = 0; i < nc; i++)
    {
      sums.at(nc+i) = std::sqrt( std::abs( sums.at(nc+i) - sums.at(i) * sums.at(i) ) );
      app_log() << boost::format("  count group %3i:  avg = %8.4f     sdev = %8.4f") % i % sums.at(i) % sums.at(nc+i) << std::endl;
    }
    if ( nc ) app_log() << std::endl;

    // only record the values if requested and if we don't have values already
    if ( record && tfl_sdv < 0.0 ) {
      tfl_avg = tla;
      tfl_sdv = tlv;
      app_log() << boost::format("Note: after warmup, nodeless guiding tfl_avg was set to %.6f") % tfl_avg << std::endl;
      app_log() << std::endl;
    }

    app_log() << " === Done Printing Statistics Related to Nodeless Guiding ===" << std::endl
              << std::endl;

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Resets info stored for trial function logarithm.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::reset_tfl()
{
  tfl_avg = 0.0;
  tfl_sdv = -1.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Initialize internal nodeless guiding info for the provided configuration.
///
/// \param[in]      P        holds the configuration information
/// \param[in]      tfl      trial function logarithm
///
/// \return  the value of the overall guiding function after initialization
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::init_nodeless(const ParticleSet & P, const RealType tfl)
{

  // if we don't need to initialize, just return the already saved guiding function value
  if ( nodelessInitialized ) return savedGF;

  // get dimensions
  const int np = P.R.size(); // number of particles
  const int nc = cgCountSigmas.size(); // number of counting groups
  const int nd = mdBetas.size(); // number of min distance centers

  // get a temporary particle position that will be useful
  ParticleSet::SingleParticlePos_t tpp;

  // initialize un-normalized counting functions and their normalizing constants
  cgUnormalized.assign(np * nc, 0);
  cgUnormalizedTmp.assign(nc, 0);
  cgNorms.assign(np, 0);
  for (int i = 0; i < np; i++) {
    for (int k = 0; k < nc; k++) {
      for (int l = cgGaussStarts.at(k); l < cgGaussEnds.at(k); l++) {
        tpp = cgGaussCenters.at(l) - P.R[i];
        cgUnormalized.at(k+i*nc) += cgGaussAlphas.at(l) * std::exp( -dot(tpp,tpp) / ( 2.0 * cgGaussSigmas.at(l) * cgGaussSigmas.at(l) ) );
      }
      cgNorms.at(i) += cgUnormalized.at(k+i*nc);
    }
  }

  // compute each counting group's total count
  cgCounts.assign(nc, 0);
  cgCountsTmp.assign(nc, 0);
  for (int i = 0; i < np; i++)
    for (int k = 0; k < nc; k++)
      cgCounts.at(k) += cgUnormalized.at(k+i*nc) / cgNorms.at(i);

  // get sum of counting group penalty exponents
  cgPenaltyExponent = 0;
  for (int k = 0; k < nc; k++)
    cgPenaltyExponent -=   ( cgCountNelecs.at(k) - cgCounts.at(k) ) * ( cgCountNelecs.at(k) - cgCounts.at(k) )
                         / ( 2.0 * cgCountSigmas.at(k) * cgCountSigmas.at(k) );

  // get product of min distance penalties
  mdPenalty = 1;
  mdPenalties.assign(np, 1);
  for (int i = 0; i < np; i++) {
    RealType max_val = 0;
    for (int k = 0; k < nd; k++) {
      tpp = mdCenters.at(k) - P.R[i];
      max_val = std::max(max_val, 1.0 / ( 1.0 + std::exp( mdBetas.at(k) * ( std::sqrt( std::abs( dot(tpp,tpp) ) ) - mdDists.at(k) ) ) ) );
    }
    mdPenalties.at(i) = max_val;
    mdPenalty *= max_val;
  }

  // initialize the nodeless adjustment as the epsilon-scaled "average" trial function value
  RealType nodelessAdj = NodelessEpsilon * std::exp( 2.0 * tfl_avg );

  // penalize the nodeless adjustment by the min distance and counting group penalties
  nodelessAdj *= mdPenalty * std::exp(cgPenaltyExponent);

  // if we're not using nodeless guiding yet, ues the usual |Psi|^2 guiding function
  if ( tfl_sdv < 0.0 )
    savedGF = std::exp( 2.0 * tfl );

  // otherwise, use the trial function square norm plus penalized nodeless adjustment
  else
    savedGF = std::exp( 2.0 * tfl ) + nodelessAdj;

  // record that we are now initialized
  nodelessInitialized = true;

  // return the guiding function value
  return savedGF;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Update internal nodeless guiding info after a single particle move
///
/// \param[in]      P        holds the configuration information
/// \param[in]      iat      index of the moved particle
/// \param[in]      tfl      trial function logarithm after the move
///
/// \return  the value of the overall guiding function after the update
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::update_nodeless(const ParticleSet & P, const int iat, const RealType tfl)
{

  // problem if we are not initialized
  if ( !nodelessInitialized )
    APP_ABORT("VMCUpdatePbyPNodeless not initialized when we called update_nodeless");

  // get dimensions
  const int np = P.R.size(); // number of particles
  const int nc = cgCountSigmas.size(); // number of counting groups
  const int nd = mdBetas.size(); // number of min distance centers

  // get index component that will be used many times
  const int iat_nc = iat * nc;

  // get a temporary particle position that will be useful
  ParticleSet::SingleParticlePos_t tpp;

  // remove the moved particle's old contribution to each count
  for (int k = 0; k < nc; k++)
    cgCountsTmp[k] = cgCounts[k] - cgUnormalized[k+iat_nc] / cgNorms[iat];

  // initialize un-normalized counting function and its norm for the moved particle
  cgNormTmp = 0;
  for (int k = 0; k < nc; k++)
  {
    cgUnormalizedTmp[k] = 0;
    for (int l = cgGaussStarts[k]; l < cgGaussEnds[k]; l++)
    {
      tpp = cgGaussCenters[l] - P.R[iat];
      cgUnormalizedTmp[k] += cgGaussAlphas[l] * std::exp( -dot(tpp,tpp) / ( 2.0 * cgGaussSigmas[l] * cgGaussSigmas[l] ) );
    }
    cgNormTmp += cgUnormalizedTmp[k];
  }

  // add in the moved particle's new contribution to each count
  for (int k = 0; k < nc; k++)
    cgCountsTmp[k] += cgUnormalizedTmp[k] / cgNormTmp;

  // get sum of counting group penalty exponents
  cgPenaltyExponentTmp = 0;
  for (int k = 0; k < nc; k++)
    cgPenaltyExponentTmp -=   ( cgCountNelecs[k] - cgCountsTmp[k] ) * ( cgCountNelecs[k] - cgCountsTmp[k] )
                            / ( 2.0 * cgCountSigmas[k] * cgCountSigmas[k] );

  // update product of min distance penalties
  {
    RealType max_val = 0;
    for (int k = 0; k < nd; k++) {
      tpp = mdCenters[k] - P.R[iat];
      max_val = std::max(max_val, 1.0 / ( 1.0 + std::exp( mdBetas[k] * ( std::sqrt( std::abs( dot(tpp,tpp) ) ) - mdDists[k] ) ) ) );
    }
    mdPenaltyTmp = mdPenalty * max_val / mdPenalties[iat];
    mdPenaltiesTmp = max_val;
  }

  // initialize the nodeless adjustment as the epsilon-scaled "average" trial function value
  RealType nodelessAdj = NodelessEpsilon * std::exp( 2.0 * tfl_avg );

  // penalize the nodeless adjustment by the min distance and counting group penalties
  nodelessAdj *= mdPenaltyTmp * std::exp(cgPenaltyExponentTmp);

  // if we're not using nodeless guiding yet, ues the usual |Psi|^2 guiding function
  if ( tfl_sdv < 0.0 )
    savedGFTmp = std::exp( 2.0 * tfl );

  // remember the new the trial function square norm plus penalized nodeless adjustment
  else
    savedGFTmp = std::exp( 2.0 * tfl ) + nodelessAdj;

  // return the guiding function value
  return savedGFTmp;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Upon accepting a move, permanantly record data for the new configuration
///
/// \param[in]      P        holds the configuration information
/// \param[in]      iat      index of the moved particle
/// \param[in]      tfl      trial function logarithm after the move
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::nodeless_accept(const ParticleSet & P, const int iat, const RealType tfl)
{

  // problem if we are not initialized
  if ( !nodelessInitialized )
    APP_ABORT("VMCUpdatePbyPNodeless not initialized when we called nodeless_accept");

  // get dimensions
  const int nc = cgCountSigmas.size(); // number of counting groups

  // record the new counting region counts
  std::copy(cgCountsTmp.begin(), cgCountsTmp.end(), cgCounts.begin());

  // record new un-normalized counting function data
  std::copy(cgUnormalizedTmp.begin(), cgUnormalizedTmp.end(), &cgUnormalized.at(iat*nc));

  // record new norm
  cgNorms[iat] = cgNormTmp;

  // record new sum of counting group penalty exponents
  cgPenaltyExponent = cgPenaltyExponentTmp;

  // record new min distance penalty information
  mdPenalties[iat] = mdPenaltiesTmp;
  mdPenalty = mdPenaltyTmp;

  // record new guiding function value
  savedGF = savedGFTmp;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that moves the supplied group of walkers according to a Metropolis-Hastings
///         random walk over the nodeless guiding function.
///
/// \param[in,out]  it       iterator to the first walker to be moved
/// \param[in,out]  it_end   iterator to one past the last walker to be moved
/// \param[in]      measure  ??? (not sure what this is for, it does not get used here)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{

  // initialize a counter for the number of walkers
  int nw = 0;

  // loop over the walkers
  for (; it != it_end; it++) {

    // get a reference to the current walker
    Walker_t& thisWalker(**it);

    // load the walker information into our particle set
    W.loadWalker(thisWalker,true);

    // get a reference to where this walker stores its trial function state
    Walker_t::Buffer_t & w_buffer(thisWalker.DataSet);

    // load trial wave function state for this walker into our trial wave function object
    Psi.copyFromBuffer(W, w_buffer);

    // initialize the old configuration's trial funciton logpsi and nodeless guiding function square norm
    Psi.evaluateLog(W);
    RealType old_tfl = Psi.getLogPsi();
    //RealType old_sqn = this->get_nodeless_gf_sqn(old_tfl, W.G, W.L);
    //RealType old_sqn = this->get_nodeless_gf_sqn_new(old_tfl);
    RealType old_sqn = this->init_nodeless(W, old_tfl);

    // loop over the sub-steps in the sampling walk
    for (int si = 0; si < nSubSteps; si++) {

      // create this sub-step's random displacements for each particle 
      makeGaussRandomWithEngine(deltaR, RandomGen);

      // initialize flag to tell if any particles moved during this substep
      bool stuck = true;

      // loop over types of particles
      for(int ig = 0; ig < W.groups(); ig++) {

        // get the mass-modified time step
        RealType sqrttau = std::sqrt(Tau*MassInvS[ig]);

        // loop over particles of this species
        for (int iat = W.first(ig); iat < W.last(ig); iat++) {

          // get the proposed displacement for this particle
          const mPosType dr = sqrttau*deltaR[iat];

          // set up the move in our particle set
          const bool move_is_legal = W.makeMoveAndCheck(iat, dr);

          // reject illegal moves
          if ( ! move_is_legal ) {
            nReject++;
            continue;
          }

//          // get the trial function ratio and the change in the gradient and laplacian for the proposed move
//          const RealType ratio = Psi.ratio(W, iat, dG, dL);
//
//          // get the gradient and laplacian at the proposed configuration
//          G = W.G + dG;
//          L = W.L + dL;

          // get the trial function ratio
          const RealType ratio = Psi.ratio(W,iat);

          // get the log of the trial function at the new configuration
          const RealType new_tfl = std::log(std::abs(ratio)) + old_tfl;

          // get the square norm of the nodeless guiding function
          //const RealType new_sqn = this->get_nodeless_gf_sqn(new_tfl, G, L);
          //const RealType new_sqn = this->get_nodeless_gf_sqn_new(new_tfl);

          // compute what would happen to the nodeless guiding object if we accepted the move
          const RealType new_sqn = this->update_nodeless(W, iat, new_tfl);

          // if the ratio of square norms satisfies the Metropolis condition, accept the move
          if ( RandomGen() < new_sqn / old_sqn ) {

            nAccept++;
            stuck = false;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            //W.G = G;
            //W.L = L;
            old_tfl = new_tfl;
            old_sqn = new_sqn;

            // record what happens to the nodeless guiding object now that the move is accepted
            this->nodeless_accept(W, iat, old_tfl);

//            {
//              // get the sum of the laplacian and square gradient of log(Psi)
//              RealType d = 0.0;
//              for (int i = 0; i < W.G.size(); i++) {
//                d += dot(W.G[i], W.G[i]);
//                d += W.L[i];
//              }
//
//              // get |Psi|^2 and |laplacian of Psi|^2
//              const RealType psi2 = std::exp( 2.0 * new_tfl );
//              const RealType lap2 = psi2 * d * d;
//            }

          // otherwise, reject the move
          } else {

            nReject++;
            W.rejectMove(iat);
            Psi.rejectMove(iat);

          }

        } // end iat loop over particles

      } // end ig loop over types of particles

      // if nothing moved, increment the all-rejected counter
      nAllRejected += ( stuck ? 1 : 0 );

      //// record our new position, gradient, and laplacian
      //thisWalker.R = W.R;
      //thisWalker.G = W.G;
      //thisWalker.L = W.L;

      // record our new position
      thisWalker.R = W.R;

    } // end si loop over sub-steps

    // for now, we are assuming there is only one walker to advance
    if ( ++nw > 1 )
      APP_ABORT("VMCUpdatePbyPNodeless::advanceWalkers encountered more than one walker");

    // get thread number
    const int ip = omp_get_thread_num();

    // save running totals for trial funciton logarithm
    tla_hist.at(ip) += old_tfl;
    tlv_hist.at(ip) += old_tfl * old_tfl;
    nsp_hist.at(ip) += 1.0;

    // save running totals for counting groups
    if ( cga_hist.at(ip).size() != cgCounts.size() ) cga_hist.at(ip).assign(cgCounts.size(), 0.0);
    if ( cgv_hist.at(ip).size() != cgCounts.size() ) cgv_hist.at(ip).assign(cgCounts.size(), 0.0);
    for (int i = 0; i < cgCounts.size(); i++)
    {
      cga_hist[ip][i] += cgCounts[i];
      cgv_hist[ip][i] += cgCounts[i] * cgCounts[i];
    }

    // save the trial wave function state
    const RealType logpsi = Psi.updateBuffer(W, w_buffer, false);

    // save the particle set information
    W.saveWalker(thisWalker);

    // evaluate the local energy
    EstimatorRealType eloc = H.evaluate(W);

    //std::printf("psi^2 = %10.2e    g = %10.2e    eloc = %10.2e    toAvg = %10.2e\n", std::exp(2.0*logpsi), old_sqn, eloc, std::exp(2.0*logpsi) * eloc / old_sqn);

    // save some basic info about the underlying trial function
    thisWalker.resetProperty(logpsi, Psi.getPhase(), eloc);

    // save the logarithm of the guiding function
    thisWalker.LogGuiding = 0.5 * std::log(old_sqn);

    // not sure what these do
    H.auxHevaluate(W, thisWalker);
    H.saveProperty(thisWalker.getPropertyBase());

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not implemented...
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::advanceWalkerForEE(Walker_t& w1,
                                                                          std::vector<PosType>& dR,
                                                                          std::vector<int>& iats,
                                                                          std::vector<int>& rs,
                                                                          std::vector<RealType>& ratios)
{
  APP_ABORT("VMCUpdatePbyPNodeless::advanceWalkerForEE not implemented");
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the square norm of the nodeless guiding function.
///         This guiding function is defined to be
///
///         g = sqrt( psi^2 + epsilon * (laplacian of psi)^2 / ( 1 + exp( alpha * ( logpsi - mu + sigma ) / sigma ) ) )
///
///         where epsilon and alpha are user-chosen parameters and mu and sigma are the average and standard deviation
///         of a previously-sampled set of trial function logarithm values.
///
/// \param[in]      logpsi   the log of the trial function
/// \param[in]      grad     the gradient of the log of the trial function
/// \param[in]      lap      the laplacian of the log of the trial function
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::get_nodeless_gf_sqn(const RealType logpsi,
                                                                           const ParticleSet::ParticleGradient_t & grad,
                                                                           const ParticleSet::ParticleLaplacian_t & lap) const
{

  // get the sum of the laplacian and square gradient of log(Psi)
  RealType d = 0.0;
  for (int i = 0; i < grad.size(); i++) {
    d += dot(grad[i], grad[i]);
    d += lap[i];
  }

  // get |Psi|^2 and |laplacian of Psi|^2
  const RealType psi2 = std::exp( 2.0 * logpsi );
  const RealType lap2 = psi2 * d * d;

  // If we don't have an average and standard deviation for a previously-taken set of
  // trial function logarithms, use a simpler nodeless guiding function instead.
  // The idea is to use this for the very first warmup when we have no history
  // to work with and then to switch to the general formula after that.
  if ( tfl_sdv < 0.0 )
    return psi2 + NodelessEpsilon * lap2 / 100.0 ;

  // return the square norm of the nodeless guiding function
  return psi2 + NodelessEpsilon * lap2 / ( 1.0 + std::exp( NodelessAlpha * ( logpsi - tfl_avg + tfl_sdv ) / tfl_sdv ) );

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the new nodeless guiding function.
///         This guiding function is defined to be
///
///                   psi^2 + epsilon * <psi^2> * prod_i ( 1 / ( 1 + exp( beta * ( rm_i - alpha ) ) ) )
///
///         where epsilon, beta, and alpha are user-chosen parameters rm_i is the distance from
///         the ith electron to the nearest nucleus
///
/// \param[in]      logpsi   the log of the trial function
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::get_nodeless_gf_sqn_new(const RealType logpsi) const
{

  // If we don't have an average and standard deviation for a previously-taken set of
  // trial function logarithms, use a simpler nodeless guiding function instead.
  // The idea is to use this for the very first warmup when we have no history
  // to work with and then to switch to the general formula after that.
  if ( tfl_sdv < 0.0 )
    return std::exp( 2.0 * logpsi );

  // !!! THIS IS NOT QUITE RIGHT.  YOU NEED TO ACTUALLY AVERAGE PSI^2 !!!
  // initialize return based on average psi2
  RealType retval = NodelessEpsilon * std::exp( 2.0 * tfl_avg );

  ParticleSet::SingleParticlePos_t temp_vec;

  // loop over types of particles
  for(int ig = 0; ig < W.groups(); ig++) {

    //std::cout << "ig = " << ig << std::endl;

    // loop over particles of this species
    for (int iat = W.first(ig); iat < W.last(ig); iat++) {

      //auto temp_vec = 1.0 * W.R[iat];

      //std::cout << "  iat = " << iat;
      //std::cout << "    W.R[iat] = " << W.R[iat] << std::endl;

      // get the smallest electron-ion distance for this electron
      RealType min_dist = 1.0e100;
      for (auto it = IonPositions.begin(); it != IonPositions.end(); it++) {
        temp_vec = *it - W.R[iat];
        min_dist = std::min(min_dist, std::abs(std::sqrt(dot(temp_vec, temp_vec))));
      }

      //std::printf("  %10.2e", min_dist);

      // apply the distance penalty for this electron
      retval *= 1.0 / ( 1.0 + std::exp( NodelessBeta * ( min_dist - NodelessAlpha ) ) );

    }

  }

  // print what happened
  //std::printf("    psi^2 = %10.2e    adjust by %10.2e", std::exp( 2.0 * logpsi ), retval);
  //std::cout << std::endl;

  // add psi2
  retval += std::exp( 2.0 * logpsi );

  //APP_ABORT("VMCUpdatePbyPNodeless::get_nodeless_gf_sqn_new stopping here");

  // return nodeless guiding function value
  return retval;

}

} // end namespace qmcplusplus
