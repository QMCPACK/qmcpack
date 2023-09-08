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

/** @file VirtualParticleSet.cpp
 * A proxy class to the quantum ParticleSet
 */

#include "VirtualParticleSetT.h"

#include "Particle/DistanceTableT.h"
#include "Particle/createDistanceTableT.h"
#include "QMCHamiltonians/NLPPJob.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{

struct VPMultiWalkerMem : public Resource
{
    /// multi walker reference particle
    Vector<int, OffloadPinnedAllocator<int>> mw_refPctls;

    VPMultiWalkerMem() : Resource("VPMultiWalkerMem")
    {
    }

    VPMultiWalkerMem(const VPMultiWalkerMem&) : VPMultiWalkerMem()
    {
    }

    std::unique_ptr<Resource>
    makeClone() const override
    {
        return std::make_unique<VPMultiWalkerMem>(*this);
    }
};

template <typename T>
VirtualParticleSetT<T>::VirtualParticleSetT(
    const ParticleSetT<T>& p, int nptcl, size_t dt_count_limit) :
    ParticleSetT<T>(p.getSimulationCell())
{
    this->setName("virtual");

    // initialize local data structure
    this->setSpinor(p.isSpinor());
    this->TotalNum = nptcl;
    this->R.resize(nptcl);
    if (this->isSpinor())
        this->spins.resize(nptcl);
    this->coordinates_->resize(nptcl);

    // create distancetables
    assert(dt_count_limit <= p.getNumDistTables());
    if (dt_count_limit == 0)
        dt_count_limit = p.getNumDistTables();
    for (int i = 0; i < dt_count_limit; ++i)
        if (p.getDistTable(i).getModes() & DTModes::NEED_VP_FULL_TABLE_ON_HOST)
            this->addTable(p.getDistTable(i).get_origin());
        else
            this->addTable(p.getDistTable(i).get_origin(),
                DTModes::MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST);
}

template <typename T>
VirtualParticleSetT<T>::~VirtualParticleSetT() = default;

template <typename T>
Vector<int, OffloadPinnedAllocator<int>>&
VirtualParticleSetT<T>::getMultiWalkerRefPctls()
{
    return mw_mem_handle_.getResource().mw_refPctls;
}

template <typename T>
const Vector<int, OffloadPinnedAllocator<int>>&
VirtualParticleSetT<T>::getMultiWalkerRefPctls() const
{
    return mw_mem_handle_.getResource().mw_refPctls;
}

template <typename T>
void
VirtualParticleSetT<T>::createResource(ResourceCollection& collection) const
{
    collection.addResource(std::make_unique<VPMultiWalkerMem>());
    ParticleSetT<T>::createResource(collection);
}

template <typename T>
void
VirtualParticleSetT<T>::acquireResource(ResourceCollection& collection,
    const RefVectorWithLeader<VirtualParticleSetT>& vp_list)
{
    auto& vp_leader = vp_list.getLeader();
    vp_leader.mw_mem_handle_ = collection.lendResource<VPMultiWalkerMem>();

    auto p_list = RefVectorWithLeaderParticleSet(vp_list);
    ParticleSetT<T>::acquireResource(collection, p_list);
}

template <typename T>
void
VirtualParticleSetT<T>::releaseResource(ResourceCollection& collection,
    const RefVectorWithLeader<VirtualParticleSetT>& vp_list)
{
    collection.takebackResource(vp_list.getLeader().mw_mem_handle_);
    auto p_list = RefVectorWithLeaderParticleSet(vp_list);
    ParticleSetT<T>::releaseResource(collection, p_list);
}

/// move virtual particles to new postions and update distance tables
template <typename T>
void
VirtualParticleSetT<T>::makeMoves(const ParticleSetT<T>& refp, int jel,
    const std::vector<PosType>& deltaV, bool sphere, int iat)
{
    if (sphere && iat < 0)
        throw std::runtime_error(
            "VirtualParticleSet::makeMoves is invoked incorrectly, the flag "
            "sphere=true requires iat specified!");
    onSphere = sphere;
    refPS = refp;
    refPtcl = jel;
    refSourcePtcl = iat;
    assert(this->R.size() == deltaV.size());
    for (size_t ivp = 0; ivp < this->R.size(); ivp++)
        this->R[ivp] = refp.R[jel] + deltaV[ivp];
    if (refp.isSpinor())
        for (size_t ivp = 0; ivp < this->R.size(); ivp++)
            this->spins[ivp] = refp.spins[jel]; // no spin deltas in this API
    this->update();
}

/// move virtual particles to new postions and update distance tables
template <typename T>
void
VirtualParticleSetT<T>::makeMovesWithSpin(const ParticleSetT<T>& refp, int jel,
    const std::vector<PosType>& deltaV, const std::vector<RealType>& deltaS,
    bool sphere, int iat)
{
    assert(refp.isSpinor());
    if (sphere && iat < 0)
        throw std::runtime_error(
            "VirtualParticleSet::makeMovesWithSpin is invoked incorrectly, the "
            "flag sphere=true requires iat specified!");
    onSphere = sphere;
    refPS = refp;
    refPtcl = jel;
    refSourcePtcl = iat;
    assert(this->R.size() == deltaV.size());
    assert(this->spins.size() == deltaS.size());
    for (size_t ivp = 0; ivp < this->R.size(); ivp++) {
        this->R[ivp] = refp.R[jel] + deltaV[ivp];
        this->spins[ivp] = refp.spins[jel] + deltaS[ivp];
    }
    this->update();
}

template <typename T>
void
VirtualParticleSetT<T>::mw_makeMoves(
    const RefVectorWithLeader<VirtualParticleSetT>& vp_list,
    const RefVectorWithLeader<ParticleSetT<T>>& refp_list,
    const RefVector<const std::vector<PosType>>& deltaV_list,
    const RefVector<const NLPPJob<RealType>>& joblist, bool sphere)
{
    auto& vp_leader = vp_list.getLeader();
    vp_leader.onSphere = sphere;
    vp_leader.refPS = refp_list.getLeader();

    const size_t nVPs = countVPs(vp_list);
    auto& mw_refPctls = vp_leader.getMultiWalkerRefPctls();
    mw_refPctls.resize(nVPs);

    RefVectorWithLeader<ParticleSetT<T>> p_list(vp_leader);
    p_list.reserve(vp_list.size());

    size_t ivp = 0;
    for (int iw = 0; iw < vp_list.size(); iw++) {
        VirtualParticleSetT& vp(vp_list[iw]);
        const std::vector<PosType>& deltaV(deltaV_list[iw]);
        const NLPPJob<RealType>& job(joblist[iw]);

        vp.onSphere = sphere;
        vp.refPS = refp_list[iw];
        vp.refPtcl = job.electron_id;
        vp.refSourcePtcl = job.ion_id;
        assert(vp.R.size() == deltaV.size());
        for (size_t k = 0; k < vp.R.size(); k++, ivp++) {
            vp.R[k] = refp_list[iw].R[vp.refPtcl] + deltaV[k];
            if (vp_leader.isSpinor())
                vp.spins[k] =
                    refp_list[iw]
                        .spins[vp.refPtcl]; // no spin deltas in this API
            mw_refPctls[ivp] = vp.refPtcl;
        }
        p_list.push_back(vp);
    }
    assert(ivp == nVPs);

    mw_refPctls.updateTo();
    ParticleSetT<T>::mw_update(p_list);
}

template <typename T>
void
VirtualParticleSetT<T>::mw_makeMovesWithSpin(
    const RefVectorWithLeader<VirtualParticleSetT>& vp_list,
    const RefVectorWithLeader<ParticleSetT<T>>& refp_list,
    const RefVector<const std::vector<PosType>>& deltaV_list,
    const RefVector<const std::vector<RealType>>& deltaS_list,
    const RefVector<const NLPPJob<RealType>>& joblist, bool sphere)
{
    auto& vp_leader = vp_list.getLeader();
    if (!vp_leader.isSpinor())
        throw std::runtime_error(
            "VirtualParticleSet::mw_makeMovesWithSpin should not be called if "
            "particle sets aren't spionor types");
    vp_leader.onSphere = sphere;
    vp_leader.refPS = refp_list.getLeader();

    const size_t nVPs = countVPs(vp_list);
    auto& mw_refPctls = vp_leader.getMultiWalkerRefPctls();
    mw_refPctls.resize(nVPs);

    RefVectorWithLeader<ParticleSetT<T>> p_list(vp_leader);
    p_list.reserve(vp_list.size());

    size_t ivp = 0;
    for (int iw = 0; iw < vp_list.size(); iw++) {
        VirtualParticleSetT& vp(vp_list[iw]);
        const std::vector<PosType>& deltaV(deltaV_list[iw]);
        const std::vector<RealType>& deltaS(deltaS_list[iw]);
        const NLPPJob<RealType>& job(joblist[iw]);

        vp.onSphere = sphere;
        vp.refPS = refp_list[iw];
        vp.refPtcl = job.electron_id;
        vp.refSourcePtcl = job.ion_id;
        assert(vp.R.size() == deltaV.size());
        assert(vp.spins.size() == deltaS.size());
        assert(vp.R.size() == vp.spins.size());
        for (size_t k = 0; k < vp.R.size(); k++, ivp++) {
            vp.R[k] = refp_list[iw].R[vp.refPtcl] + deltaV[k];
            vp.spins[k] = refp_list[iw].spins[vp.refPtcl] + deltaS[k];
            mw_refPctls[ivp] = vp.refPtcl;
        }
        p_list.push_back(vp);
    }
    assert(ivp == nVPs);

    mw_refPctls.updateTo();
    ParticleSetT<T>::mw_update(p_list);
}

#ifndef QMC_COMPLEX
template class VirtualParticleSetT<double>;
template class VirtualParticleSetT<float>;
#else
template class VirtualParticleSetT<std::complex<double>>;
template class VirtualParticleSetT<std::complex<float>>;
#endif
} // namespace qmcplusplus
