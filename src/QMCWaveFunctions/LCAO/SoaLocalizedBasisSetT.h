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

/** @file SoaLocalizedBasisSetT.h
 * @brief A derived class from BasisSetBase
 *
 * This is intended as a replacement for MolecularWaveFunctionComponent and
 * any other localized basis set.
 */
#ifndef QMCPLUSPLUS_SOA_LOCALIZEDBASISSETT_H
#define QMCPLUSPLUS_SOA_LOCALIZEDBASISSETT_H

#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "QMCWaveFunctions/BasisSetBaseT.h"

#include <memory>

namespace qmcplusplus
{
/** A localized basis set derived from SoaBasisSetBase<ORBT>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 * The template parameter ORBT denotes the orbital value return type
 */
template <class COT, typename ORBT>
class SoaLocalizedBasisSetT : public SoaBasisSetBaseT<ORBT>
{
public:
    using RealType = typename COT::RealType;
    using BaseType = SoaBasisSetBaseT<ORBT>;
    using ValueType = ORBT;

    using vgl_type = typename BaseType::vgl_type;
    using vgh_type = typename BaseType::vgh_type;
    using vghgh_type = typename BaseType::vghgh_type;
    using PosType = typename ParticleSetT<ORBT>::PosType;
    using OffloadMWVGLArray = typename BaseType::OffloadMWVGLArray;
    using OffloadMWVArray = typename BaseType::OffloadMWVArray;

    using BaseType::BasisSetSize;

    /// number of centers, e.g., ions
    size_t NumCenters;
    /// number of quantum particles
    size_t NumTargets;
    /// ion particle set
    const ParticleSetT<ORBT>& ions_;
    /// number of quantum particles
    const int myTableIndex;
    /// Global Coordinate of Supertwist read from HDF5
    PosType SuperTwist;

    /** container to store the offsets of the basis functions for each center
     * Due to potential reordering of ions, offsets can be in any order.
     */
    std::vector<size_t> BasisOffset;

    /** container of the unique pointers to the Atomic Orbitals
     *
     * size of LOBasisSet = number  of unique centers
     */
    std::vector<std::unique_ptr<COT>> LOBasisSet;

    /** constructor
     * @param ions ionic system
     * @param els electronic system
     */
    SoaLocalizedBasisSetT(ParticleSetT<ORBT>& ions, ParticleSetT<ORBT>& els);

    /** copy constructor */
    SoaLocalizedBasisSetT(const SoaLocalizedBasisSetT& a);

    /** makeClone */
    BaseType*
    makeClone() const override
    {
        return new SoaLocalizedBasisSetT<COT, ORBT>(*this);
    }

    /** set Number of periodic Images to evaluate the orbitals.
        Set to 0 for non-PBC, and set manually in the input.
        Passes the pre-computed phase factor for evaluation of complex
       wavefunction. If WF is real Phase_factor is real and equals 1 if gamma or
       -1 if non-Gamma.
    */
    void
    setPBCParams(const TinyVector<int, 3>& PBCImages,
        const TinyVector<double, 3> Sup_Twist,
        const std::vector<ORBT>& phase_factor);

    /** set BasisSetSize and allocate mVGL container
     */
    void
    setBasisSetSize(int nbs) override;

    /**  Determine which orbitals are S-type.  Used by cusp correction.
     */
    void
    queryOrbitalsForSType(const std::vector<bool>& corrCenter,
        std::vector<bool>& is_s_orbital) const override;

    /** compute VGL
     * @param P quantum particleset
     * @param iat active particle
     * @param vgl Matrix(5,BasisSetSize)
     * @param trialMove if true, use getTempDists()/getTempDispls()
     */
    void
    evaluateVGL(const ParticleSetT<ORBT>& P, int iat, vgl_type& vgl) override;

    /** compute V using packed array with all walkers
     * @param P_list list of quantum particleset (one for each walker)
     * @param iat active particle
     * @param v   Array(n_walkers, BasisSetSize)
     */
    void
    mw_evaluateValue(const RefVectorWithLeader<ParticleSetT<ORBT>>& P_list,
        int iat, OffloadMWVArray& v) override;

    /** compute VGL using packed array with all walkers
     * @param P_list list of quantum particleset (one for each walker)
     * @param iat active particle
     * @param vgl   Array(n_walkers, 5, BasisSetSize)
     */
    void
    mw_evaluateVGL(const RefVectorWithLeader<ParticleSetT<ORBT>>& P_list,
        int iat, OffloadMWVGLArray& vgl) override;

    /** compute VGH
     * @param P quantum particleset
     * @param iat active particle
     * @param vgl Matrix(10,BasisSetSize)
     * @param trialMove if true, use getTempDists()/getTempDispls()
     */
    void
    evaluateVGH(const ParticleSetT<ORBT>& P, int iat, vgh_type& vgh) override;

    /** compute VGHGH
     * @param P quantum particleset
     * @param iat active particle
     * @param vghgh Matrix(20,BasisSetSize)
     * @param trialMove if true, use getTempDists()/getTempDispls()
     */
    void
    evaluateVGHGH(
        const ParticleSetT<ORBT>& P, int iat, vghgh_type& vghgh) override;

    /** compute values for the iat-paricle move
     *
     * Always uses getTempDists() and getTempDispls()
     * Tv is a translation vector; In PBC, in order to reduce the number
     * of images that need to be summed over when generating the AO the
     * nearest image displacement, dr, is used. Tv corresponds to the
     * translation that takes the 'general displacement' (displacement
     * between ion position and electron position) to the nearest image
     * displacement. We need to keep track of Tv because it must be add
     * as a phase factor, i.e., exp(i*k*Tv).
     */
    void
    evaluateV(
        const ParticleSetT<ORBT>& P, int iat, ORBT* restrict vals) override;

    void
    evaluateGradSourceV(const ParticleSetT<ORBT>& P, int iat,
        const ParticleSetT<ORBT>& ions, int jion, vgl_type& vgl) override;

    void
    evaluateGradSourceVGL(const ParticleSetT<ORBT>& P, int iat,
        const ParticleSetT<ORBT>& ions, int jion, vghgh_type& vghgh) override;

    /** add a new set of Centered Atomic Orbitals
     * @param icenter the index of the center
     * @param aos a set of Centered Atomic Orbitals
     */
    void
    add(int icenter, std::unique_ptr<COT> aos);
};
} // namespace qmcplusplus
#endif