//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SIMULATIONCELLT_H
#define QMCPLUSPLUS_SIMULATIONCELLT_H

#include "LongRange/KContainerT.h"
#include "ParticleSetTraits.h"

namespace qmcplusplus
{
template <typename T>
class ParticleSetPoolT;

template <typename T>
class SimulationCellT
{
public:
    using Lattice = typename LatticeParticleTraits<T>::ParticleLayout;

    SimulationCellT();
    SimulationCellT(const Lattice& lattice);

    const Lattice&
    getLattice() const
    {
        return lattice_;
    }
    const Lattice&
    getPrimLattice() const
    {
        return primative_lattice_;
    }
    const Lattice&
    getLRBox() const
    {
        return LRBox_;
    }

    void
    resetLRBox();

    /// access k_lists_ read only
    const KContainerT<T>&
    getKLists() const
    {
        return k_lists_;
    }

private:
    /// simulation cell lattice
    Lattice lattice_;
    /// Primative cell lattice
    Lattice primative_lattice_;
    /// long-range box
    Lattice LRBox_;

    /// K-Vector List.
    KContainerT<T> k_lists_;

    friend class ParticleSetPoolT<T>;
};
} // namespace qmcplusplus
#endif
