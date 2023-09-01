//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National
// Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National
// Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SOA_LCAO_SPINOR_BUILDERT_H
#define QMCPLUSPLUS_SOA_LCAO_SPINOR_BUILDERT_H

#include "QMCWaveFunctions/LCAO/LCAOrbitalBuilderT.h"

namespace qmcplusplus
{
/** @file LCAOSpinorBuidler.h
 *
 * Derives from LCAOrbitalBuilder.h. Overrides createSPOSetFromXML method to
 * read up and down channel from HDF5 and construct SpinorSet
 *
 */
template <class T>
class LCAOSpinorBuilderT : public LCAOrbitalBuilderT<T>
{
public:
    using BasisSet_t = typename LCAOrbitalBuilderT<T>::BasisSet_t;
    using RealType = typename LCAOrbitalBuilderT<T>::RealType;
    using ValueType = typename LCAOrbitalBuilderT<T>::ValueType;

    /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     *
     * Derives from LCAOrbitalBuilder, but will require an h5_path to be set
     */
    LCAOSpinorBuilderT(ParticleSetT<T>& els, ParticleSetT<T>& ions,
        Communicate* comm, xmlNodePtr cur);

    /** creates and returns SpinorSet
     *
     * Creates an up and down LCAOrbitalSet
     * calls LCAOSpinorBuilder::loadMO to build up and down from the H5 file
     * registers up and down into a SpinorSet and returns
     */
    std::unique_ptr<SPOSetT<T>>
    createSPOSetFromXML(xmlNodePtr cur) override;

private:
    /** load the up and down MO sets
     *
     * checks to make sure not PBC and initialize the Occ vector.
     * call putFromH5 to parse the up and down MO coefficients
     */
    bool
    loadMO(LCAOrbitalSetT<T>& up, LCAOrbitalSetT<T>& dn, xmlNodePtr cur);

    /** parse h5 file for spinor info
     *
     * assumes the h5 file has KPTS_0/eigenset_0(_imag) for the real/imag part
     * of up component of spinor assumes the h5 file as KPTS_0/eigenset_1(_imag)
     * for the real/imag part of dn component of spinor reads the various
     * coefficient matricies and broadcast after this, we have up/dn
     * LCAOrbitalSet that can be registered to the SpinorSet
     */
    bool
    putFromH5(LCAOrbitalSetT<T>& up, LCAOrbitalSetT<T>& dn, xmlNodePtr);
};
} // namespace qmcplusplus
#endif
