//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at
// Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at
//                    Urbana-Champaign Jeremy McMinnis, jmcminis@gmail.com,
//                    University of Illinois at Urbana-Champaign Jeongnim Kim,
//                    jeongnim.kim@gmail.com, University of Illinois at
//                    Urbana-Champaign Cynthia Gu, zg1@ornl.gov, Oak Ridge
//                    National Laboratory Ye Luo, yeluo@anl.gov, Argonne
//                    National Laboratory Mark A. Berrill, berrillma@ornl.gov,
//                    Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "MCWalkerConfigurationT.h"

#include "LongRange/StructFact.h"
#include "Message/CommOperators.h"
#include "Message/Communicate.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/MCSample.h"
#include "Particle/ReptileT.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/IteratorUtility.h"
#include "hdf/HDFVersion.h"
#include "hdf/hdf_hyperslab.h"

#include <map>

namespace qmcplusplus
{
template <typename T>
MCWalkerConfigurationT<T>::MCWalkerConfigurationT(
    const SimulationCellT<T>& simulation_cell,
    const DynamicCoordinateKind kind) :
    ParticleSetT<T>(simulation_cell, kind),
    ReadyForPbyP(false),
    UpdateMode(Update_Walker),
    reptile(0),
    Polymer(0)
{
}

template <typename T>
MCWalkerConfigurationT<T>::MCWalkerConfigurationT(
    const MCWalkerConfigurationT& mcw) :
    ParticleSetT<T>(mcw),
    ReadyForPbyP(false),
    UpdateMode(Update_Walker),
    Polymer(0)
{
    samples.clearEnsemble();
    samples.setMaxSamples(mcw.getMaxSamples());
    this->setWalkerOffsets(mcw.getWalkerOffsets());
    this->Properties = mcw.Properties;
}

template <typename T>
MCWalkerConfigurationT<T>::~MCWalkerConfigurationT() = default;

template <typename T>
void
MCWalkerConfigurationT<T>::createWalkers(int n)
{
    const int old_nw = this->getActiveWalkers();
    WalkerConfigurationsT<T>::createWalkers(n, this->TotalNum);
    // no pre-existing walkers, need to initialized based on particleset.
    if (old_nw == 0)
        for (auto& awalker : this->walker_list_) {
            awalker->R = this->R;
            awalker->spins = this->spins;
        }
    resizeWalkerHistories();
}

template <typename T>
void
MCWalkerConfigurationT<T>::resize(int numWalkers, int numPtcls)
{
    if (this->TotalNum && this->walker_list_.size())
        app_warning()
            << "MCWalkerConfiguration::resize cleans up the walker list."
            << std::endl;
    const int old_nw = this->getActiveWalkers();
    ParticleSetT<T>::resize(unsigned(numPtcls));
    WalkerConfigurationsT<T>::resize(numWalkers, this->TotalNum);
    // no pre-existing walkers, need to initialized based on particleset.
    if (old_nw == 0)
        for (auto& awalker : this->walker_list_) {
            awalker->R = this->R;
            awalker->spins = this->spins;
        }
}

/** Make Metropolis move to the walkers and save in a temporary array.
 * @param it the iterator of the first walker to work on
 * @param tauinv  inverse of the time step
 *
 * R + D + X
 */
template <typename T>
void
MCWalkerConfigurationT<T>::sample(iterator it, RealType tauinv)
{
    throw std::runtime_error("MCWalkerConfiguration::sample obsolete");
    //  makeGaussRandom(R);
    //  R *= tauinv;
    //  R += (*it)->R + (*it)->Drift;
}

/** reset the Property container of all the walkers
 */
template <typename T>
void
MCWalkerConfigurationT<T>::resetWalkerProperty(int ncopy)
{
    int m(this->PropertyList.size());
    app_log() << "  Resetting Properties of the walkers " << ncopy << " x " << m
              << std::endl;
    try {
        this->Properties.resize(ncopy, m);
    }
    catch (std::domain_error& de) {
        app_error() << de.what() << '\n'
                    << "This is likely because some object has attempted to "
                       "add walker properties\n"
                    << " in excess of WALKER_MAX_PROPERTIES.\n"
                    << "build with cmake ... "
                       "-DWALKER_MAX_PROPERTIES=at_least_properties_required"
                    << std::endl;
        APP_ABORT("Fatal Exception");
    }

    for (auto& walker : this->walker_list_) {
        walker->resizeProperty(ncopy, m);
        walker->Weight = 1.0;
    }
    resizeWalkerHistories();
}

template <typename T>
void
MCWalkerConfigurationT<T>::resizeWalkerHistories()
{
    // using std::vector<std::vector<RealType> > is too costly.
    int np = this->PropertyHistory.size();
    if (np)
        for (int iw = 0; iw < this->walker_list_.size(); ++iw)
            this->walker_list_[iw]->PropertyHistory = this->PropertyHistory;
    np = this->PHindex.size();
    if (np)
        for (int iw = 0; iw < this->walker_list_.size(); ++iw)
            this->walker_list_[iw]->PHindex = this->PHindex;
    ;
}

/** allocate the SampleStack
 * @param n number of samples per thread
 */
template <typename T>
void
MCWalkerConfigurationT<T>::setNumSamples(int n)
{
    samples.clearEnsemble();
    samples.setMaxSamples(n);
}

/** save the current walkers to SampleStack
 */
template <typename T>
void
MCWalkerConfigurationT<T>::saveEnsemble()
{
    saveEnsemble(this->walker_list_.begin(), this->walker_list_.end());
}

/** save the [first,last) walkers to SampleStack
 */
template <typename T>
void
MCWalkerConfigurationT<T>::saveEnsemble(iterator first, iterator last)
{
    for (; first != last; first++) {
        samples.appendSample(MCSample(**first));
    }
}
/** load a single sample from SampleStack
 */
template <typename T>
void
MCWalkerConfigurationT<T>::loadSample(ParticleSetT<T>& pset, size_t iw) const
{
    samples.loadSample(pset, iw);
}

/** load SampleStack to walker_list_
 */
template <typename T>
void
MCWalkerConfigurationT<T>::loadEnsemble()
{
    using WP = WalkerProperties::Indexes;
    int nsamples = std::min(samples.getMaxSamples(), samples.getNumSamples());
    if (samples.empty() || nsamples == 0)
        return;
    typename Walker_t::PropertyContainer_t prop(
        1, this->PropertyList.size(), 1, WP::MAXPROPERTIES);
    this->walker_list_.resize(nsamples);
    for (int i = 0; i < nsamples; ++i) {
        auto awalker = std::make_unique<Walker_t>(this->TotalNum);
        awalker->Properties.copy(prop);
        samples.getSample(i).convertToWalker(*awalker);
        this->walker_list_[i] = std::move(awalker);
    }
    resizeWalkerHistories();
    samples.clearEnsemble();
}

template <typename T>
bool
MCWalkerConfigurationT<T>::dumpEnsemble(
    std::vector<MCWalkerConfigurationT<T>*>& others, HDFWalkerOutput& out,
    int np, int nBlock)
{
    WalkerConfigurationsT<T> wctemp;
    for (auto* mcwc : others) {
        const auto& astack(mcwc->getSampleStack());
        const size_t sample_size =
            std::min(mcwc->getMaxSamples(), mcwc->numSamples());
        for (int j = 0; j < sample_size; ++j) {
            const auto& sample = astack.getSample(j);
            const size_t num_ptcls = sample.getNumPtcls();
            auto awalker = std::make_unique<Walker_t>(num_ptcls);
            sample.convertToWalker(*awalker);
            wctemp.push_back(std::move(awalker));
        }
    }
    const int w = wctemp.getActiveWalkers();
    if (w == 0)
        return false;

    // The following code assumes the same amount of active walkers on all the
    // MPI ranks
    std::vector<int> nwoff(np + 1, 0);
    for (int ip = 0; ip < np; ++ip)
        nwoff[ip + 1] = nwoff[ip] + w;
    wctemp.setWalkerOffsets(nwoff);
    out.dump(wctemp, nBlock);
    return true;
}

template <typename T>
int
MCWalkerConfigurationT<T>::getMaxSamples() const
{
    return samples.getMaxSamples();
}

template <typename T>
void
MCWalkerConfigurationT<T>::loadEnsemble(
    std::vector<MCWalkerConfigurationT<T>*>& others, bool doclean)
{
    using WP = WalkerProperties::Indexes;
    std::vector<int> off(others.size() + 1, 0);
    for (int i = 0; i < others.size(); ++i) {
        off[i + 1] = off[i] +
            std::min(others[i]->getMaxSamples(), others[i]->numSamples());
    }
    int nw_tot = off.back();
    if (nw_tot) {
        typename Walker_t::PropertyContainer_t prop(
            1, this->PropertyList.size(), 1, WP::MAXPROPERTIES);
        while (this->walker_list_.size())
            this->pop_back();
        this->walker_list_.resize(nw_tot);
        for (int i = 0; i < others.size(); ++i) {
            SampleStackT<T>& astack(others[i]->getSampleStack());
            for (int j = 0, iw = off[i]; iw < off[i + 1]; ++j, ++iw) {
                auto awalker = std::make_unique<Walker_t>(this->TotalNum);
                awalker->Properties.copy(prop);
                astack.getSample(j).convertToWalker(*awalker);
                this->walker_list_[iw] = std::move(awalker);
            }
            if (doclean)
                others[i]->clearEnsemble();
        }
    }
    if (doclean)
        resizeWalkerHistories();
}

template <typename T>
void
MCWalkerConfigurationT<T>::clearEnsemble()
{
    samples.clearEnsemble();
}

#ifndef QMC_COMPLEX
#ifndef MIXED_PRECISION
template class MCWalkerConfigurationT<double>;
#else
template class MCWalkerConfigurationT<float>;
#endif
#else
#ifndef MIXED_PRECISION
template class MCWalkerConfigurationT<std::complex<double>>;
#else
template class MCWalkerConfigurationT<std::complex<float>>;
#endif
#endif

} // namespace qmcplusplus
