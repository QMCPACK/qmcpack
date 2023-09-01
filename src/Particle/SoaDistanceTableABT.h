//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DTDIMPL_ABT_H
#define QMCPLUSPLUS_DTDIMPL_ABT_H

#include "Concurrency/OpenMP.h"
#include "Lattice/ParticleBConds3DSoa.h"
#include "Particle/DistanceTableT.h"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a
 * transposed form
 */
template <typename T, unsigned D, int SC>
struct SoaDistanceTableABT :
    public DTD_BConds<typename ParticleSetTraits<T>::RealType, D, SC>,
    public DistanceTableABT<T>
{
    using RealType = typename DistanceTableABT<T>::RealType;
    using PosType = typename DistanceTableABT<T>::PosType;
    using IndexType = typename DistanceTableABT<T>::IndexType;

    SoaDistanceTableABT(
        const ParticleSetT<T>& source, ParticleSetT<T>& target) :
        DTD_BConds<RealType, D, SC>(source.getLattice()),
        DistanceTableABT<T>(source, target, DTModes::ALL_OFF),
        evaluate_timer_(createGlobalTimer(std::string("DTAB::evaluate_") +
                target.getName() + "_" + source.getName(),
            timer_level_fine)),
        move_timer_(createGlobalTimer(std::string("DTAB::move_") +
                target.getName() + "_" + source.getName(),
            timer_level_fine)),
        update_timer_(createGlobalTimer(std::string("DTAB::update_") +
                target.getName() + "_" + source.getName(),
            timer_level_fine))
    {
        resize();
    }

    void
    resize()
    {
        if (this->num_sources_ * this->num_targets_ == 0)
            return;

        // initialize memory containers and views
        const int num_sources_padded = getAlignedSize<RealType>(this->num_sources_);
        this->distances_.resize(this->num_targets_);
        this->displacements_.resize(this->num_targets_);
        for (int i = 0; i < this->num_targets_; ++i) {
            this->distances_[i].resize(num_sources_padded);
            this->displacements_[i].resize(num_sources_padded);
        }

        // The padding of temp_r_ and temp_dr_ is necessary for the memory copy
        // in the update function temp_r_ is padded explicitly while temp_dr_ is
        // padded internally
        this->temp_r_.resize(num_sources_padded);
        this->temp_dr_.resize(this->num_sources_);
    }

    SoaDistanceTableABT() = delete;
    SoaDistanceTableABT(const SoaDistanceTableABT&) = delete;

    /** evaluate the full table */
    inline void
    evaluate(ParticleSetT<T>& P) override
    {
        ScopedTimer local_timer(evaluate_timer_);
#pragma omp parallel
        {
            int first, last;
            FairDivideAligned(this->num_sources_, getAlignment<RealType>(),
                omp_get_num_threads(), omp_get_thread_num(), first, last);

            // be aware of the sign of Displacement
            for (int iat = 0; iat < this->num_targets_; ++iat)
                DTD_BConds<RealType, D, SC>::computeDistances(P.R[iat],
                    this->origin_.getCoordinates().getAllParticlePos(),
                    this->distances_[iat].data(), this->displacements_[iat],
                    first, last);
        }
    }

    /// evaluate the temporary pair relations
    inline void
    move(const ParticleSetT<T>& P, const PosType& rnew, const IndexType iat,
        bool prepare_old) override
    {
        ScopedTimer local_timer(move_timer_);
        DTD_BConds<RealType, D, SC>::computeDistances(rnew,
            this->origin_.getCoordinates().getAllParticlePos(), this->temp_r_.data(),
            this->temp_dr_, 0, this->num_sources_);
        // If the full table is not ready all the time, overwrite the current
        // value. If this step is missing, DT values can be undefined in case a
        // move is rejected.
        if (!(this->modes_ & DTModes::NEED_FULL_TABLE_ANYTIME) && prepare_old)
            DTD_BConds<RealType, D, SC>::computeDistances(P.R[iat],
                this->origin_.getCoordinates().getAllParticlePos(),
                this->distances_[iat].data(), this->displacements_[iat], 0,
                this->num_sources_);
    }

    /// update the stripe for jat-th particle
    inline void
    update(IndexType iat) override
    {
        ScopedTimer local_timer(update_timer_);
        std::copy_n(this->temp_r_.data(), this->num_sources_,
            this->distances_[iat].data());
        for (int idim = 0; idim < D; ++idim)
            std::copy_n(this->temp_dr_.data(idim), this->num_sources_,
                this->displacements_[iat].data(idim));
    }

    int
    get_first_neighbor(
        IndexType iat, RealType& r, PosType& dr, bool newpos) const override
    {
        RealType min_dist = std::numeric_limits<RealType>::max();
        int index = -1;
        if (newpos) {
            for (int jat = 0; jat < this->num_sources_; ++jat)
                if (this->temp_r_[jat] < min_dist) {
                    min_dist = this->temp_r_[jat];
                    index = jat;
                }
            if (index >= 0) {
                r = min_dist;
                dr = this->temp_dr_[index];
            }
        }
        else {
            for (int jat = 0; jat < this->num_sources_; ++jat)
                if (this->distances_[iat][jat] < min_dist) {
                    min_dist = this->distances_[iat][jat];
                    index = jat;
                }
            if (index >= 0) {
                r = min_dist;
                dr = this->displacements_[iat][index];
            }
        }
        assert(index >= 0 && index < this->num_sources_);
        return index;
    }

private:
    /// timer for evaluate()
    NewTimer& evaluate_timer_;
    /// timer for move()
    NewTimer& move_timer_;
    /// timer for update()
    NewTimer& update_timer_;
};
} // namespace qmcplusplus
#endif
