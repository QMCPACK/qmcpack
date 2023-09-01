//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore
// National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of
//                    Illinois at Urbana-Champaign Jaron T. Krogel,
//                    krogeljt@ornl.gov, Oak Ridge National Laboratory Jeongnim
//                    Kim, jeongnim.kim@gmail.com, University of Illinois at
//                    Urbana-Champaign Ye Luo, yeluo@anl.gov, Argonne National
//                    Laboratory Mark A. Berrill, berrillma@ornl.gov, Oak Ridge
//                    National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SOA_LCAO_ORBITAL_BUILDERT_H
#define QMCPLUSPLUS_SOA_LCAO_ORBITAL_BUILDERT_H

#include "QMCWaveFunctions/LCAO/LCAOrbitalSetT.h"
#include "QMCWaveFunctions/SPOSetBuilderT.h"

#include <map>

namespace qmcplusplus
{
/** SPOSetBuilder using new LCAOrbitalSet and Soa versions
 *
 * Reimplement MolecularSPOSetBuilder
 * - support both CartesianTensor and SphericalTensor
 */
template <typename T>
class LCAOrbitalBuilderT : public SPOSetBuilderT<T>
{
public:
    using BasisSet_t = typename LCAOrbitalSetT<T>::basis_type;
    using RealType = typename LCAOrbitalSetT<T>::RealType;
    using ValueType = typename LCAOrbitalSetT<T>::ValueType;
    using PosType = typename LCAOrbitalSetT<T>::PosType;

    /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
    LCAOrbitalBuilderT(ParticleSetT<T>& els, ParticleSetT<T>& ions,
        Communicate* comm, xmlNodePtr cur);
    ~LCAOrbitalBuilderT() override;
    std::unique_ptr<SPOSetT<T>>
    createSPOSetFromXML(xmlNodePtr cur) override;

protected:
    /// target ParticleSet
    ParticleSetT<T>& targetPtcl;
    /// source ParticleSet
    ParticleSetT<T>& sourcePtcl;
    /// localized basis set map
    std::map<std::string, std::unique_ptr<BasisSet_t>> basisset_map_;
    /// if true, add cusp correction to orbitals
    bool cuspCorr;
    /// Path to HDF5 Wavefunction
    std::string h5_path;
    /// Number of periodic Images for Orbital evaluation
    TinyVector<int, 3> PBCImages;
    /// Coordinates Super Twist
    PosType SuperTwist;
    /// Periodic Image Phase Factors. Correspond to the phase from the
    /// PBCImages. Computed only once.
    std::vector<T> PeriodicImagePhaseFactors;
    /// Store Lattice parameters from HDF5 to use in PeriodicImagePhaseFactors
    Tensor<double, 3> Lattice;

    /// Enable cusp correction
    bool doCuspCorrection;

    /** create basis set
     *
     * Use ao_traits<T,I,J> to match (ROT)x(SH) combo
     */
    template <int I, int J>
    BasisSet_t*
    createBasisSet(xmlNodePtr cur);
    template <int I, int J>
    BasisSet_t*
    createBasisSetH5();

    // The following items were previously in SPOSet
    /// occupation number
    Vector<RealType> Occ;
    bool
    loadMO(LCAOrbitalSetT<T>& spo, xmlNodePtr cur);
    bool
    putOccupation(LCAOrbitalSetT<T>& spo, xmlNodePtr occ_ptr);
    bool
    putFromXML(LCAOrbitalSetT<T>& spo, xmlNodePtr coeff_ptr);
    bool
    putFromH5(LCAOrbitalSetT<T>& spo, xmlNodePtr coeff_ptr);
    bool
    putPBCFromH5(LCAOrbitalSetT<T>& spo, xmlNodePtr coeff_ptr);
    // the dimensions of Ctemp are determined by the dataset on file
    void
    LoadFullCoefsFromH5(hdf_archive& hin, int setVal, PosType& SuperTwist,
        Matrix<std::complex<RealType>>& Ctemp, bool MultiDet);
    // the dimensions of Creal are determined by the dataset on file
    void
    LoadFullCoefsFromH5(hdf_archive& hin, int setVal, PosType& SuperTwist,
        Matrix<RealType>& Creal, bool Multidet);
    void
    EvalPeriodicImagePhaseFactors(PosType SuperTwist,
        std::vector<RealType>& LocPeriodicImagePhaseFactors);
    void
    EvalPeriodicImagePhaseFactors(PosType SuperTwist,
        std::vector<std::complex<RealType>>& LocPeriodicImagePhaseFactors);
    /** read matrix from h5 file
     * \param[in] hin: hdf5 arhive to be read from
     * \param setname: where to read from in hdf5 archive
     * \param[out] Creal: matrix read from h5
     *
     * added in header to allow use from derived class LCAOSpinorBuilder as well
     */
    void
    readRealMatrixFromH5(hdf_archive& hin, const std::string& setname,
        Matrix<RealType>& Creal) const;

private:
    /// enable cusp correction
    std::unique_ptr<SPOSetT<T>>
    createWithCuspCorrection(xmlNodePtr cur, const std::string& spo_name,
        std::string cusp_file, std::unique_ptr<BasisSet_t>&& myBasisSet);
    /// load a basis set from XML input
    std::unique_ptr<BasisSet_t>
    loadBasisSetFromXML(xmlNodePtr cur, xmlNodePtr parent);
    /// load a basis set from h5 file
    std::unique_ptr<BasisSet_t>
    loadBasisSetFromH5(xmlNodePtr parent);
    /// determine radial orbital type based on "keyword" and "transform"
    /// attributes
    int
    determineRadialOrbType(xmlNodePtr cur) const;
};

} // namespace qmcplusplus
#endif
