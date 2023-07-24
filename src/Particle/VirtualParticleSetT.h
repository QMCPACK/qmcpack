//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VIRTUAL_PARTICLESETT_H
#define QMCPLUSPLUS_VIRTUAL_PARTICLESETT_H

#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "Particle/ParticleSetT.h"
#include <ResourceHandle.h>

namespace qmcplusplus
{
// forward declaration.
class NonLocalECPComponent;
template <typename T>
struct NLPPJob;
struct VPMultiWalkerMem;

/** A ParticleSet that handles virtual moves of a selected particle of a given
 * physical ParticleSet Virtual moves are defined as moves being proposed but
 * will never be accepted. VirtualParticleSet is introduced to avoid changing
 * any internal states of the physical ParticleSet. For this reason, the
 * physical ParticleSet is always marked const. It is heavily used by non-local
 * PP evaluations.
 */
template <typename T>
class VirtualParticleSetT : public ParticleSetT<T>
{
public:
    using RealType = typename ParticleSetT<T>::RealType;
    using PosType = typename ParticleSetT<T>::PosType;

private:
    /// true, if virtual particles are on a sphere for NLPP
    bool onSphere;
    /// multi walker resource
    ResourceHandle<VPMultiWalkerMem> mw_mem_handle_;

    Vector<int, OffloadPinnedAllocator<int>>&
    getMultiWalkerRefPctls();

    /// ParticleSet this object refers to after makeMoves
    std::optional<std::reference_wrapper<const ParticleSetT<T>>> refPS;

public:
    /// Reference particle
    int refPtcl;
    /// Reference source particle, used when onSphere=true
    int refSourcePtcl;

    /// ParticleSet this object refers to
    const ParticleSetT<T>&
    getRefPS() const
    {
        return refPS.value();
    }

    inline bool
    isOnSphere() const
    {
        return onSphere;
    }

    const Vector<int, OffloadPinnedAllocator<int>>&
    getMultiWalkerRefPctls() const;

    /** constructor
     * @param p ParticleSet whose virtual moves are handled by this object
     * @param nptcl number of virtual particles
     * @param dt_count_limit distance tables corresepond to [0, dt_count_limit)
     * of the reference particle set are created
     */
    VirtualParticleSetT(
        const ParticleSetT<T>& p, int nptcl, size_t dt_count_limit = 0);

    ~VirtualParticleSetT();

    /// initialize a shared resource and hand it to a collection
    void
    createResource(ResourceCollection& collection) const;
    /** acquire external resource and assocaite it with the list of ParticleSet
     * Note: use RAII ResourceCollectionTeamLock whenever possible
     */
    static void
    acquireResource(ResourceCollection& collection,
        const RefVectorWithLeader<VirtualParticleSetT>& vp_list);
    /** release external resource
     * Note: use RAII ResourceCollectionTeamLock whenever possible
     */
    static void
    releaseResource(ResourceCollection& collection,
        const RefVectorWithLeader<VirtualParticleSetT>& vp_list);

    /** move virtual particles to new postions and update distance tables
     * @param refp reference particle set
     * @param jel reference particle that all the VP moves from
     * @param deltaV Position delta for virtual moves.
     * @param sphere set true if VP are on a sphere around the reference source
     * particle
     * @param iat reference source particle
     */
    void
    makeMoves(const ParticleSetT<T>& refp, int jel,
        const std::vector<PosType>& deltaV, bool sphere = false, int iat = -1);

    /** move virtual particles to new postions and update distance tables
     * @param refp reference particle set
     * @param jel reference particle that all the VP moves from
     * @param deltaV Position delta for virtual moves.
     * @param deltaS Spin delta for virtual moves.
     * @param sphere set true if VP are on a sphere around the reference source
     * particle
     * @param iat reference source particle
     */
    void
    makeMovesWithSpin(const ParticleSetT<T>& refp, int jel,
        const std::vector<PosType>& deltaV, const std::vector<RealType>& deltaS,
        bool sphere = false, int iat = -1);

    static void
    mw_makeMoves(const RefVectorWithLeader<VirtualParticleSetT>& vp_list,
        const RefVectorWithLeader<ParticleSetT<T>>& p_list,
        const RefVector<const std::vector<PosType>>& deltaV_list,
        const RefVector<const NLPPJob<RealType>>& joblist, bool sphere);

    static void
    mw_makeMovesWithSpin(
        const RefVectorWithLeader<VirtualParticleSetT>& vp_list,
        const RefVectorWithLeader<ParticleSetT<T>>& p_list,
        const RefVector<const std::vector<PosType>>& deltaV_list,
        const RefVector<const std::vector<RealType>>& deltaS_list,
        const RefVector<const NLPPJob<RealType>>& joblist, bool sphere);

    static RefVectorWithLeader<ParticleSetT<T>>
    RefVectorWithLeaderParticleSet(
        const RefVectorWithLeader<VirtualParticleSetT>& vp_list)
    {
        RefVectorWithLeader<ParticleSetT<T>> ref_list(vp_list.getLeader());
        ref_list.reserve(ref_list.size());
        for (VirtualParticleSetT& vp : vp_list)
            ref_list.push_back(vp);
        return ref_list;
    }

    static size_t
    countVPs(const RefVectorWithLeader<const VirtualParticleSetT>& vp_list)
    {
        size_t nVPs = 0;
        for (const VirtualParticleSetT& vp : vp_list)
            nVPs += vp.getTotalNum();
        return nVPs;
    }

    static size_t
    countVPs(const RefVectorWithLeader<VirtualParticleSetT>& vp_list)
    {
        size_t nVPs = 0;
        for (const VirtualParticleSetT& vp : vp_list)
            nVPs += vp.getTotalNum();
        return nVPs;
    }
};
} // namespace qmcplusplus
#endif