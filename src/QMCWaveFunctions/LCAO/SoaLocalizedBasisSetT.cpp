//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#include "SoaLocalizedBasisSetT.h"

#include "MultiFunctorAdapter.h"
#include "MultiQuinticSpline1D.h"
#include "Numerics/SoaCartesianTensor.h"
#include "Numerics/SoaSphericalTensor.h"
#include "Particle/DistanceTableT.h"
#include "SoaAtomicBasisSetT.h"

#include <memory>

namespace qmcplusplus
{
template <class COT, typename ORBT>
SoaLocalizedBasisSetT<COT, ORBT>::SoaLocalizedBasisSetT(
    ParticleSetT<ORBT>& ions, ParticleSetT<ORBT>& els) :
    ions_(ions),
    myTableIndex(els.addTable(ions,
        DTModes::NEED_FULL_TABLE_ANYTIME |
            DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
    SuperTwist(0.0)
{
    NumCenters = ions.getTotalNum();
    NumTargets = els.getTotalNum();
    LOBasisSet.resize(ions.getSpeciesSet().getTotalNum());
    BasisOffset.resize(NumCenters + 1);
    BasisSetSize = 0;
}

template <class COT, typename ORBT>
SoaLocalizedBasisSetT<COT, ORBT>::SoaLocalizedBasisSetT(
    const SoaLocalizedBasisSetT& a) :
    SoaBasisSetBaseT<ORBT>(a),
    NumCenters(a.NumCenters),
    NumTargets(a.NumTargets),
    ions_(a.ions_),
    myTableIndex(a.myTableIndex),
    SuperTwist(a.SuperTwist),
    BasisOffset(a.BasisOffset)
{
    LOBasisSet.reserve(a.LOBasisSet.size());
    for (auto& elem : a.LOBasisSet)
        LOBasisSet.push_back(std::make_unique<COT>(*elem));
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::setPBCParams(
    const TinyVector<int, 3>& PBCImages, const TinyVector<double, 3> Sup_Twist,
    const std::vector<ORBT>& phase_factor)
{
    for (int i = 0; i < LOBasisSet.size(); ++i)
        LOBasisSet[i]->setPBCParams(PBCImages, Sup_Twist, phase_factor);

    SuperTwist = Sup_Twist;
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::setBasisSetSize(int nbs)
{
    const auto& IonID(ions_.GroupID);
    if (BasisSetSize > 0 && nbs == BasisSetSize)
        return;

    if (auto& mapping = ions_.get_map_storage_to_input(); mapping.empty()) {
        // evaluate the total basis dimension and offset for each center
        BasisOffset[0] = 0;
        for (int c = 0; c < NumCenters; c++)
            BasisOffset[c + 1] =
                BasisOffset[c] + LOBasisSet[IonID[c]]->getBasisSetSize();
        BasisSetSize = BasisOffset[NumCenters];
    }
    else {
        // when particles are reordered due to grouping, AOs need to restore the
        // input order to match MOs.
        std::vector<int> map_input_to_storage(mapping.size());
        for (int c = 0; c < NumCenters; c++)
            map_input_to_storage[mapping[c]] = c;

        std::vector<size_t> basis_offset_input_order(BasisOffset.size(), 0);
        for (int c = 0; c < NumCenters; c++)
            basis_offset_input_order[c + 1] = basis_offset_input_order[c] +
                LOBasisSet[IonID[map_input_to_storage[c]]]->getBasisSetSize();

        for (int c = 0; c < NumCenters; c++)
            BasisOffset[c] = basis_offset_input_order[mapping[c]];

        BasisSetSize = basis_offset_input_order[NumCenters];
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::queryOrbitalsForSType(
    const std::vector<bool>& corrCenter, std::vector<bool>& is_s_orbital) const
{
    const auto& IonID(ions_.GroupID);
    for (int c = 0; c < NumCenters; c++) {
        int idx = BasisOffset[c];
        int bss = LOBasisSet[IonID[c]]->BasisSetSize;
        std::vector<bool> local_is_s_orbital(bss);
        LOBasisSet[IonID[c]]->queryOrbitalsForSType(local_is_s_orbital);
        for (int k = 0; k < bss; k++) {
            if (corrCenter[c]) {
                is_s_orbital[idx++] = local_is_s_orbital[k];
            }
            else {
                is_s_orbital[idx++] = false;
            }
        }
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::evaluateVGL(
    const ParticleSetT<ORBT>& P, int iat, vgl_type& vgl)
{
    const auto& IonID(ions_.GroupID);
    const auto& coordR = P.activeR(iat);
    const auto& d_table = P.getDistTableAB(myTableIndex);
    const auto& dist = (P.getActivePtcl() == iat) ? d_table.getTempDists() :
                                                    d_table.getDistRow(iat);
    const auto& displ = (P.getActivePtcl() == iat) ? d_table.getTempDispls() :
                                                     d_table.getDisplRow(iat);

    PosType Tv;
    for (int c = 0; c < NumCenters; c++) {
        Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
        Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
        Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
        LOBasisSet[IonID[c]]->evaluateVGL(
            P.getLattice(), dist[c], displ[c], BasisOffset[c], vgl, Tv);
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::mw_evaluateVGL(
    const RefVectorWithLeader<ParticleSetT<ORBT>>& P_list, int iat,
    OffloadMWVGLArray& vgl_v)
{
    for (size_t iw = 0; iw < P_list.size(); iw++) {
        const auto& IonID(ions_.GroupID);
        const auto& coordR = P_list[iw].activeR(iat);
        const auto& d_table = P_list[iw].getDistTableAB(myTableIndex);
        const auto& dist = (P_list[iw].getActivePtcl() == iat) ?
            d_table.getTempDists() :
            d_table.getDistRow(iat);
        const auto& displ = (P_list[iw].getActivePtcl() == iat) ?
            d_table.getTempDispls() :
            d_table.getDisplRow(iat);

        PosType Tv;

        // number of walkers * BasisSetSize
        auto stride = vgl_v.size(1) * BasisSetSize;
        assert(BasisSetSize == vgl_v.size(2));
        vgl_type vgl_iw(vgl_v.data_at(0, iw, 0), BasisSetSize, stride);

        for (int c = 0; c < NumCenters; c++) {
            Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
            Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
            Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
            LOBasisSet[IonID[c]]->evaluateVGL(P_list[iw].getLattice(), dist[c],
                displ[c], BasisOffset[c], vgl_iw, Tv);
        }
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::evaluateVGH(
    const ParticleSetT<ORBT>& P, int iat, vgh_type& vgh)
{
    const auto& IonID(ions_.GroupID);
    const auto& d_table = P.getDistTableAB(myTableIndex);
    const auto& dist = (P.getActivePtcl() == iat) ? d_table.getTempDists() :
                                                    d_table.getDistRow(iat);
    const auto& displ = (P.getActivePtcl() == iat) ? d_table.getTempDispls() :
                                                     d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++) {
        LOBasisSet[IonID[c]]->evaluateVGH(
            P.getLattice(), dist[c], displ[c], BasisOffset[c], vgh);
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::evaluateVGHGH(
    const ParticleSetT<ORBT>& P, int iat, vghgh_type& vghgh)
{
    // APP_ABORT("SoaLocalizedBasisSetT::evaluateVGH() not implemented\n");

    const auto& IonID(ions_.GroupID);
    const auto& d_table = P.getDistTableAB(myTableIndex);
    const auto& dist = (P.getActivePtcl() == iat) ? d_table.getTempDists() :
                                                    d_table.getDistRow(iat);
    const auto& displ = (P.getActivePtcl() == iat) ? d_table.getTempDispls() :
                                                     d_table.getDisplRow(iat);
    for (int c = 0; c < NumCenters; c++) {
        LOBasisSet[IonID[c]]->evaluateVGHGH(
            P.getLattice(), dist[c], displ[c], BasisOffset[c], vghgh);
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::evaluateV(
    const ParticleSetT<ORBT>& P, int iat, ORBT* restrict vals)
{
    const auto& IonID(ions_.GroupID);
    const auto& coordR = P.activeR(iat);
    const auto& d_table = P.getDistTableAB(myTableIndex);
    const auto& dist = (P.getActivePtcl() == iat) ? d_table.getTempDists() :
                                                    d_table.getDistRow(iat);
    const auto& displ = (P.getActivePtcl() == iat) ? d_table.getTempDispls() :
                                                     d_table.getDisplRow(iat);

    PosType Tv;
    for (int c = 0; c < NumCenters; c++) {
        Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
        Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
        Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
        LOBasisSet[IonID[c]]->evaluateV(
            P.getLattice(), dist[c], displ[c], vals + BasisOffset[c], Tv);
    }
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::mw_evaluateValue(
    const RefVectorWithLeader<ParticleSetT<ORBT>>& P_list, int iat,
    OffloadMWVArray& v)
{
    for (size_t iw = 0; iw < P_list.size(); iw++)
        evaluateV(P_list[iw], iat, v.data_at(iw, 0));
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::evaluateGradSourceV(
    const ParticleSetT<ORBT>& P, int iat, const ParticleSetT<ORBT>& ions,
    int jion, vgl_type& vgl)
{
    // We need to zero out the temporary array vgl.
    auto* restrict gx = vgl.data(1);
    auto* restrict gy = vgl.data(2);
    auto* restrict gz = vgl.data(3);

    for (int ib = 0; ib < BasisSetSize; ib++) {
        gx[ib] = 0;
        gy[ib] = 0;
        gz[ib] = 0;
    }

    const auto& IonID(ions_.GroupID);
    const auto& d_table = P.getDistTableAB(myTableIndex);
    const auto& dist = (P.getActivePtcl() == iat) ? d_table.getTempDists() :
                                                    d_table.getDistRow(iat);
    const auto& displ = (P.getActivePtcl() == iat) ? d_table.getTempDispls() :
                                                     d_table.getDisplRow(iat);

    PosType Tv;
    Tv[0] = Tv[1] = Tv[2] = 0;
    // Since LCAO's are written only in terms of (r-R), ionic derivatives only
    // exist for the atomic center that we wish to take derivatives of.
    // Moreover, we can obtain an ion derivative by multiplying an electron
    // derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For
    // now, just note this is the electron VGL function.
    LOBasisSet[IonID[jion]]->evaluateVGL(
        P.getLattice(), dist[jion], displ[jion], BasisOffset[jion], vgl, Tv);
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::evaluateGradSourceVGL(
    const ParticleSetT<ORBT>& P, int iat, const ParticleSetT<ORBT>& ions,
    int jion, vghgh_type& vghgh)
{
    // We need to zero out the temporary array vghgh.
    auto* restrict gx = vghgh.data(1);
    auto* restrict gy = vghgh.data(2);
    auto* restrict gz = vghgh.data(3);

    auto* restrict hxx = vghgh.data(4);
    auto* restrict hxy = vghgh.data(5);
    auto* restrict hxz = vghgh.data(6);
    auto* restrict hyy = vghgh.data(7);
    auto* restrict hyz = vghgh.data(8);
    auto* restrict hzz = vghgh.data(9);

    auto* restrict gxxx = vghgh.data(10);
    auto* restrict gxxy = vghgh.data(11);
    auto* restrict gxxz = vghgh.data(12);
    auto* restrict gxyy = vghgh.data(13);
    auto* restrict gxyz = vghgh.data(14);
    auto* restrict gxzz = vghgh.data(15);
    auto* restrict gyyy = vghgh.data(16);
    auto* restrict gyyz = vghgh.data(17);
    auto* restrict gyzz = vghgh.data(18);
    auto* restrict gzzz = vghgh.data(19);

    for (int ib = 0; ib < BasisSetSize; ib++) {
        gx[ib] = 0;
        gy[ib] = 0;
        gz[ib] = 0;

        hxx[ib] = 0;
        hxy[ib] = 0;
        hxz[ib] = 0;
        hyy[ib] = 0;
        hyz[ib] = 0;
        hzz[ib] = 0;

        gxxx[ib] = 0;
        gxxy[ib] = 0;
        gxxz[ib] = 0;
        gxyy[ib] = 0;
        gxyz[ib] = 0;
        gxzz[ib] = 0;
        gyyy[ib] = 0;
        gyyz[ib] = 0;
        gyzz[ib] = 0;
        gzzz[ib] = 0;
    }

    // Since jion is indexed on the source ions not the ions_ the distinction
    // between ions_ and ions is extremely important.
    const auto& IonID(ions.GroupID);
    const auto& d_table = P.getDistTableAB(myTableIndex);
    const auto& dist = (P.getActivePtcl() == iat) ? d_table.getTempDists() :
                                                    d_table.getDistRow(iat);
    const auto& displ = (P.getActivePtcl() == iat) ? d_table.getTempDispls() :
                                                     d_table.getDisplRow(iat);

    // Since LCAO's are written only in terms of (r-R), ionic derivatives only
    // exist for the atomic center that we wish to take derivatives of.
    // Moreover, we can obtain an ion derivative by multiplying an electron
    // derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For
    // now, just note this is the electron VGL function.

    LOBasisSet[IonID[jion]]->evaluateVGHGH(
        P.getLattice(), dist[jion], displ[jion], BasisOffset[jion], vghgh);
}

template <class COT, typename ORBT>
void
SoaLocalizedBasisSetT<COT, ORBT>::add(int icenter, std::unique_ptr<COT> aos)
{
    LOBasisSet[icenter] = std::move(aos);
}

// TODO: this should be redone with template template parameters

template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<double>, SoaCartesianTensor<double>,
        double>,
    double>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<double>, SoaCartesianTensor<double>,
        std::complex<double>>,
    std::complex<double>>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<float>, SoaCartesianTensor<float>,
        float>,
    float>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<float>, SoaCartesianTensor<float>,
        std::complex<float>>,
    std::complex<float>>;

template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<double>, SoaSphericalTensor<double>,
        double>,
    double>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<double>, SoaSphericalTensor<double>,
        std::complex<double>>,
    std::complex<double>>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<float>, SoaSphericalTensor<float>,
        float>,
    float>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiQuinticSpline1D<float>, SoaSphericalTensor<float>,
        std::complex<float>>,
    std::complex<float>>;

template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<double>>,
        SoaCartesianTensor<double>, double>,
    double>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<double>>,
        SoaCartesianTensor<double>, std::complex<double>>,
    std::complex<double>>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<float>>,
        SoaCartesianTensor<float>, float>,
    float>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<float>>,
        SoaCartesianTensor<float>, std::complex<float>>,
    std::complex<float>>;

template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<double>>,
        SoaSphericalTensor<double>, double>,
    double>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<double>>,
        SoaSphericalTensor<double>, std::complex<double>>,
    std::complex<double>>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<float>>,
        SoaSphericalTensor<float>, float>,
    float>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<GaussianCombo<float>>,
        SoaSphericalTensor<float>, std::complex<float>>,
    std::complex<float>>;

template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<double>>,
        SoaCartesianTensor<double>, double>,
    double>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<double>>,
        SoaCartesianTensor<double>, std::complex<double>>,
    std::complex<double>>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<float>>,
        SoaCartesianTensor<float>, float>,
    float>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<float>>,
        SoaCartesianTensor<float>, std::complex<float>>,
    std::complex<float>>;

template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<double>>,
        SoaSphericalTensor<double>, double>,
    double>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<double>>,
        SoaSphericalTensor<double>, std::complex<double>>,
    std::complex<double>>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<float>>,
        SoaSphericalTensor<float>, float>,
    float>;
template class SoaLocalizedBasisSetT<
    SoaAtomicBasisSetT<MultiFunctorAdapter<SlaterCombo<float>>,
        SoaSphericalTensor<float>, std::complex<float>>,
    std::complex<float>>;
} // namespace qmcplusplus