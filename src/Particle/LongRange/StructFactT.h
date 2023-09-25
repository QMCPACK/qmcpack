//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of
// Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_STRUCTFACTT_H
#define QMCPLUSPLUS_STRUCTFACTT_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/ParticleSetTraits.h"
#include <NewTimer.h>
#include <OMPTarget/OffloadAlignedAllocators.hpp>
#include <Resource.h>
#include <type_traits/template_types.hpp>

namespace qmcplusplus
{
template <typename T>
class ParticleSetT;
template <typename T>
class KContainerT;
template <typename T>
struct SKMultiWalkerMemT;

/** @ingroup longrange
 *\brief Calculates the structure-factor for a particle set
 *
 * Structure factor per species
 *   Rhok[alpha][k] \f$ \equiv \rho_{k}^{\alpha} = \sum_{i} e^{i{\bf k}\cdot{\bf
 *r_i}}\f$ Structure factor per particle eikr[i][k]
 */
template <typename T>
class StructFactT
{
public:
    // Typedef for the lattice-type
    using ParticleLayout = typename LatticeParticleTraits<T>::ParticleLayout;
    using RealType = typename ParticleSetTraits<T>::RealType;

    static constexpr auto DIM = ParticleSetTraits<T>::DIM;

    /** enumeration for the methods to handle mixed bconds
     *
     * Allow overwriting lattice::SuperCellEnum to use D-dim k-point sets with
     * mixed BC
     */
    int SuperCellEnum;
    /// 2-D container for the phase
    Matrix<RealType> rhok_r, rhok_i;
    Matrix<RealType> eikr_r, eikr_i;
    /** Constructor - copy ParticleSet and init. k-shells
     * @param lattice long range box
     * @param kc cutoff for k
     *
     * At least in the batched version Structure factor is _NOT_ valid
     * after construction.
     */
    StructFactT(const ParticleLayout& lattice, const KContainerT<T>& k_lists);
    /// desructor
    ~StructFactT();

    /**  Update Rhok if all particles moved
     */
    void
    updateAllPart(const ParticleSetT<T>& P);

    /** Update RhoK for all particles for multiple walkers particles.
     *
     *  In batched context until this is called StructFact is invalid and will
     * cause a crash if any Hamiltonian using StructFact indirectly through
     * ParticleSet is evaluated.
     */
    static void
    mw_updateAllPart(const RefVectorWithLeader<StructFactT>& sk_list,
        const RefVectorWithLeader<ParticleSetT<T>>& p_list,
        SKMultiWalkerMemT<T>& mw_mem);

    /** @brief switch on the storage per particle
     * if StorePerParticle was false, this function allocates memory and
     * precompute data if StorePerParticle was true, this function is no-op
     */
    void
    turnOnStorePerParticle(const ParticleSetT<T>& P);

    /// accessor of StorePerParticle
    bool
    isStorePerParticle() const
    {
        return StorePerParticle;
    }

    /// accessor of k_lists_
    const KContainerT<T>&
    getKLists() const
    {
        return k_lists_;
    }

private:
    /// Compute all rhok elements from the start
    void
    computeRhok(const ParticleSetT<T>& P);
    /** resize the internal data
     * @param nkpts
     * @param num_species number of species
     * @param num_ptcls number of particles
     */
    void
    resize(int nkpts, int num_species, int num_ptcls);

    /// K-Vector List.
    const KContainerT<T>& k_lists_;
    /** Whether intermediate data is stored per particle. default false
     * storing data per particle needs significant amount of memory but some
     * calculation may request it. storing data per particle specie is more
     * cost-effective
     */
    bool StorePerParticle;
    /// timer for updateAllPart
    NewTimer& update_all_timer_;
};

/// multi walker shared memory buffer
template <typename T>
struct SKMultiWalkerMemT : public Resource
{
    using RealType = typename StructFactT<T>::RealType;

    /// dist displ for temporary and old pairs
    Matrix<RealType, OffloadPinnedAllocator<RealType>> nw_rhok;

    SKMultiWalkerMemT() : Resource("SKMultiWalkerMem")
    {
    }

    SKMultiWalkerMemT(const SKMultiWalkerMemT&) : SKMultiWalkerMemT()
    {
    }

    std::unique_ptr<Resource>
    makeClone() const override
    {
        return std::make_unique<SKMultiWalkerMemT>(*this);
    }
};

} // namespace qmcplusplus

#endif
