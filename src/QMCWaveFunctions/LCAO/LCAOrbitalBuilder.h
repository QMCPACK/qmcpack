//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOA_LCAO_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_SOA_LCAO_ORBITAL_BUILDER_H

#include <map>
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/LCAO/LCAOrbitalSet.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"

namespace qmcplusplus
{
/** SPOSetBuilder using new LCAOrbitalSet and Soa versions
   *
   * Reimplement MolecularSPOSetBuilder
   * - support both CartesianTensor and SphericalTensor
   */
class LCAOrbitalBuilder : public SPOSetBuilder
{
public:
  using BasisSet_t = LCAOrbitalSet::basis_type;
  /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
  LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, Communicate* comm, xmlNodePtr cur);
  ~LCAOrbitalBuilder() override;
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

  // For testing
  using BasissetMap = std::map<std::string, std::unique_ptr<BasisSet_t>>;
  const BasissetMap& getBasissetMap() const { return basisset_map_; }

protected:
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///source ParticleSet
  ParticleSet& sourcePtcl;
  /// localized basis set map
  std::map<std::string, std::unique_ptr<BasisSet_t>> basisset_map_;
  /// if true, add cusp correction to orbitals
  bool cuspCorr;
  ///Path to HDF5 Wavefunction
  std::string h5_path;
  ///Number of periodic Images for Orbital evaluation
  TinyVector<int, 3> PBCImages;
  ///Coordinates Super Twist
  PosType SuperTwist;
  ///Periodic Image Phase Factors. Correspond to the phase from the PBCImages. Computed only once.
  Vector<ValueType, OffloadPinnedAllocator<ValueType>> PeriodicImagePhaseFactors;
  Array<RealType, 2, OffloadPinnedAllocator<RealType>> PeriodicImageDisplacements;
  ///Store Lattice parameters from HDF5 to use in PeriodicImagePhaseFactors
  Tensor<double, 3> Lattice;

  /// Enable cusp correction
  bool doCuspCorrection;
  /// Captured gpu input string
  std::string useGPU;

  /** Create a localized basis set from an XML input node
   *
   * Processes the atomicBasisSet elements per ion species and builds
   * the complete localized basis set. Uses ao_traits<T,I,J> to match
   * the appropriate (Radial Orbital Type) x (Spherical Harmonics) combinations.
   *
   * @param cur pointer to the XML node containing the basis set definitions
   * @return a unique pointer to the newly constructed BasisSet_t
   */
  template<int I, int J>
  std::unique_ptr<BasisSet_t> createBasisSet(xmlNodePtr cur) const;

  /** Create a localized basis set from an HDF5 file
   *
   * Reads atomic basis set parameters and definitions from an HDF5 archive
   * (specified by h5_path) and builds the localized basis set. Uses ao_traits<T,I,J>
   * to match the appropriate combinations.
   *
   * @return a unique pointer to the newly constructed BasisSet_t
   */
  template<int I, int J>
  std::unique_ptr<BasisSet_t> createBasisSetH5() const;

  // The following items were previously in SPOSet
  ///occupation number
  Vector<RealType> Occ;
  bool loadMO(LCAOrbitalSet& spo, xmlNodePtr cur);
  bool putOccupation(LCAOrbitalSet& spo, xmlNodePtr occ_ptr);
  bool putFromXML(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr);
  bool putFromH5(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr);
  bool putPBCFromH5(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr);
  // the dimensions of Ctemp are determined by the dataset on file
  void LoadFullCoefsFromH5(hdf_archive& hin,
                           int setVal,
                           PosType& SuperTwist,
                           Matrix<std::complex<RealType>>& Ctemp,
                           bool MultiDet);
  // the dimensions of Creal are determined by the dataset on file
  void LoadFullCoefsFromH5(hdf_archive& hin, int setVal, PosType& SuperTwist, Matrix<RealType>& Creal, bool Multidet);
  void EvalPeriodicImagePhaseFactors(
      PosType SuperTwist,
      Vector<RealType, OffloadPinnedAllocator<RealType>>& LocPeriodicImagePhaseFactors,
      Array<RealType, 2, OffloadPinnedAllocator<RealType>>& LocPeriodicImageDisplacements);
  void EvalPeriodicImagePhaseFactors(
      PosType SuperTwist,
      Vector<std::complex<RealType>, OffloadPinnedAllocator<std::complex<RealType>>>& LocPeriodicImagePhaseFactors,
      Array<RealType, 2, OffloadPinnedAllocator<RealType>>& LocPeriodicImageDisplacements);
  /** read matrix from h5 file
   * \param[in] hin: hdf5 arhive to be read from
   * \param setname: where to read from in hdf5 archive
   * \param[out] Creal: matrix read from h5
   *
   * added in header to allow use from derived class LCAOSpinorBuilder as well
   */
  void readRealMatrixFromH5(hdf_archive& hin,
                            const std::string& setname,
                            Matrix<LCAOrbitalBuilder::RealType>& Creal) const;

private:
  /** Load and construct a complete localized basis set entirely from XML
   *
   * Determines the radial orbital type from the XML attributes, dispatches
   * the creation to the appropriate template instantiation of createBasisSet,
   * and returns the resulting BasisSet_t.
   *
   * @param cur pointer to the current XML node being parsed
   * @param parent pointer to the parent XML node
   * @return unique pointer to the constructed basis set
   */
  std::unique_ptr<BasisSet_t> loadBasisSetFromXML(xmlNodePtr cur, xmlNodePtr parent) const;

  /** Load and construct a complete localized basis set partially/fully from an HDF5 file
   *
   * Determines the radial orbital type from the XML attributes, dispatches
   * the creation to the appropriate template instantiation of createBasisSetH5
   * (which reads the bulk of the data from the HDF5 archive), and returns it.
   *
   * @param parent pointer to the parent XML node containing configuration attributes
   * @return unique pointer to the constructed basis set
   */
  std::unique_ptr<BasisSet_t> loadBasisSetFromH5(xmlNodePtr parent) const;
  ///determine radial orbital type based on "keyword" and "transform" attributes
  int determineRadialOrbType(xmlNodePtr cur) const;
};


} // namespace qmcplusplus
#endif
