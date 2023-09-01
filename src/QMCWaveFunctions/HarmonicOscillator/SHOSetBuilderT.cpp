//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National
// Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National
//                    Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National
// Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "SHOSetBuilderT.h"

#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/SPOSetInputInfo.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/string_utils.h"

namespace qmcplusplus
{
template <class T>
SHOSetBuilderT<T>::SHOSetBuilderT(ParticleSetT<T>& P, Communicate* comm) :
    SPOSetBuilderT<T>("SHO", comm),
    Ps(P)
{
    this->ClassName = "SHOSetBuilderT";
    this->legacy = false;
    app_log() << "Constructing SHOSetBuilderT" << std::endl;
    reset();
}

template <class T>
SHOSetBuilderT<T>::~SHOSetBuilderT() = default;

template <class T>
void
SHOSetBuilderT<T>::reset()
{
    nstates = 0;
    mass = -1.0;
    energy = -1.0;
    length = -1.0;
    center = 0.0;
}

template <class T>
std::unique_ptr<SPOSetT<T>>
SHOSetBuilderT<T>::createSPOSetFromXML(xmlNodePtr cur)
{
    APP_ABORT("SHOSetBuilderT::createSPOSetFromXML  SHOSetBuilder should not "
              "use legacy interface");

    app_log() << "SHOSetBuilderT::createSHOSet(xml) " << std::endl;

    SPOSetInputInfo input(cur);

    return createSPOSet(cur, input);
}

template <class T>
std::unique_ptr<SPOSetT<T>>
SHOSetBuilderT<T>::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input)
{
    app_log() << "SHOSetBuilderT::createSHOSet(indices) " << std::endl;
    reset();

    // read parameters
    std::string spo_name = "sho";
    OhmmsAttributeSet attrib;
    attrib.add(spo_name, "name");
    attrib.add(spo_name, "id");
    attrib.add(mass, "mass");
    attrib.add(energy, "energy");
    attrib.add(energy, "frequency");
    attrib.add(length, "length");
    attrib.add(center, "center");
    attrib.add(nstates, "size");
    attrib.put(cur);

    if (energy < 0.0)
        energy = 1.0;
    if (mass < 0.0 && length < 0.0)
        length = 1.0;
    if (mass < 0.0)
        mass = 1.0 / (energy * length * length);
    else if (length < 0.0)
        length = 1.0 / std::sqrt(mass * energy);

    // initialize states and/or adjust basis
    int smax = -1;
    if (input.has_index_info)
        smax = std::max(smax, input.max_index());
    if (input.has_energy_info) {
        smax = std::max(smax, (int)std::ceil(input.max_energy() / energy));
    }
    if (smax < 0)
        APP_ABORT("SHOSetBuilderT::Initialize\n  invalid basis size");
    update_basis_states(smax);

    // create sho state request
    indices_t& indices = input.get_indices(this->states);
    std::vector<SHOState*> sho_states;
    for (int i = 0; i < indices.size(); ++i)
        sho_states.push_back(basis_states[indices[i]]);

    // make the sposet
    auto sho =
        std::make_unique<SHOSetT<T>>(spo_name, length, center, sho_states);

    sho->report("  ");
    return sho;
}

template <class T>
void
SHOSetBuilderT<T>::update_basis_states(int smax)
{
    int states_required = smax - basis_states.size() + 1;
    if (states_required > 0) {
        RealType N = smax + 1;
        if (QMCTraits::DIM == 1)
            nmax = smax;
        else if (QMCTraits::DIM == 2)
            nmax = std::ceil(.5 * std::sqrt(8. * N + 1.) - 1.5);
        else if (QMCTraits::DIM == 3) {
            RealType f = std::exp(1.0 / 3.0 *
                std::log(81. * N + 3. * std::sqrt(729. * N * N - 3.)));
            nmax = std::ceil(f / 3. + 1. / f - 2.);
        }
        else
            APP_ABORT("SHOSetBuilderT::update_basis_states  dimensions other "
                      "than 1, 2, or 3 are not supported");
        int ndim = nmax + 1;
        ind_dims[QMCTraits::DIM - 1] = 1;
        for (int d = QMCTraits::DIM - 2; d > -1; --d)
            ind_dims[d] = ind_dims[d + 1] * ndim;
        int s = 0;
        int ntot = pow(ndim, QMCTraits::DIM);
        TinyVector<int, QMCTraits::DIM> qnumber;
        for (int m = 0; m < ntot; ++m) {
            int n = 0; // principal quantum number
            int nrem = m;
            for (int d = 0; d < QMCTraits::DIM; ++d) {
                int i = nrem / ind_dims[d];
                nrem -= i * ind_dims[d];
                qnumber[d] = i;
                n += i;
            }
            if (n <= nmax) {
                SHOState* st;
                if (s < basis_states.size())
                    st = basis_states[s];
                else {
                    st = new SHOState();
                    basis_states.add(st);
                }
                RealType e = energy * (n + .5 * QMCTraits::DIM);
                st->set(qnumber, e);
                s++;
            }
        }
        basis_states.energy_sort(1e-6, true);
    }

    // reset energy scale even if no states need to be added
    for (int i = 0; i < basis_states.size(); ++i) {
        SHOState& state = *basis_states[i];
        const TinyVector<int, QMCTraits::DIM>& qnumber = state.quantum_number;
        int n = 0;
        for (int d = 0; d < QMCTraits::DIM; ++d)
            n += qnumber[d];
        state.energy = energy * (n + .5 * QMCTraits::DIM);
    }

    // somewhat redundant, but necessary
    this->clear_states(0);
    this->states[0]->finish(basis_states.states);

    if (basis_states.size() <= smax)
        APP_ABORT("SHOSetBuilderT::update_basis_states  failed to make enough "
                  "states");
}

template <class T>
void
SHOSetBuilderT<T>::report(const std::string& pad)
{
    app_log() << pad << "SHOSetBuilderT report" << std::endl;
    app_log() << pad << "  dimension = " << QMCTraits::DIM << std::endl;
    app_log() << pad << "  mass      = " << mass << std::endl;
    app_log() << pad << "  frequency = " << energy << std::endl;
    app_log() << pad << "  energy    = " << energy << std::endl;
    app_log() << pad << "  length    = " << length << std::endl;
    app_log() << pad << "  center    = " << center << std::endl;
    app_log() << pad << "  nstates   = " << nstates << std::endl;
    app_log() << pad << "  nmax      = " << nmax << std::endl;
    app_log() << pad << "  ind_dims  = " << ind_dims << std::endl;
    app_log() << pad << "  # basis states = " << basis_states.size()
              << std::endl;
    app_log() << pad << "  basis_states" << std::endl;
    for (int s = 0; s < basis_states.size(); ++s)
        basis_states[s]->report(pad + "  " + int2string(s) + " ");
    app_log() << pad << "end SHOSetBuilderT report" << std::endl;
    app_log().flush();
}

template class SHOSetBuilderT<double>;
template class SHOSetBuilderT<float>;
template class SHOSetBuilderT<std::complex<double>>;
template class SHOSetBuilderT<std::complex<float>>;

} // namespace qmcplusplus
