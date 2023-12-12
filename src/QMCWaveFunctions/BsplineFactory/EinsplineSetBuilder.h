//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file EinsplineSetBuilder.h
 *
 * Builder class for einspline-based SPOSet objects.
 */
#ifndef QMCPLUSPLUS_EINSPLINE_SET_BUILDER_H
#define QMCPLUSPLUS_EINSPLINE_SET_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/BandInfo.h"
#include <filesystem>
#include <map>

#define PW_COEFF_NORM_TOLERANCE 1e-6

class Communicate;

namespace qmcplusplus
{
///forward declaration of BsplineReader
struct BsplineReader;

// Helper needed for TwistMap
struct Int3less
{
  bool operator()(const TinyVector<int, 3>& a, const TinyVector<int, 3>& b) const
  {
    if (a[0] > b[0])
      return false;
    if (a[0] < b[0])
      return true;
    if (a[1] > b[1])
      return false;
    if (a[1] < b[1])
      return true;
    if (a[2] > b[2])
      return false;
    if (a[2] < b[2])
      return true;
    return false;
  }
};
struct Int4less
{
  bool operator()(const TinyVector<int, 4>& a, const TinyVector<int, 4>& b) const
  {
    for (int i = 0; i < 4; i++)
    {
      if (a[i] > b[i])
        return false;
      if (a[i] < b[i])
        return true;
    }
    return false;
  }
};


/** construct a name for spline SPO set
 */
struct H5OrbSet
{
  ///index for the spin set
  int SpinSet;
  ///number of orbitals that belong to this set
  int NumOrbs;
  ///name of the HDF5 file
  std::filesystem::path FileName;
  /** true if a < b
   *
   * The ordering
   * - name
   * - spin set
   * - number of orbitals
   */
  bool operator()(const H5OrbSet& a, const H5OrbSet& b) const
  {
    if (a.FileName == b.FileName)
    {
      if (a.SpinSet == b.SpinSet)
        return a.NumOrbs < b.NumOrbs;
      else
        return a.SpinSet < b.SpinSet;
    }
    else
      return a.FileName < b.FileName;
  }

  H5OrbSet(std::filesystem::path name, int spinSet, int numOrbs)
      : SpinSet(spinSet), NumOrbs(numOrbs), FileName(std::move(name))
  {}
  H5OrbSet() = default;
};

/** EinsplineSet builder
 */
class EinsplineSetBuilder : public SPOSetBuilder
{
public:
  using PSetMap      = std::map<std::string, const std::unique_ptr<ParticleSet>>;
  using UnitCellType = CrystalLattice<ParticleSet::Scalar_t, DIM>;

  ///reference to the particleset pool
  const PSetMap& ParticleSets;
  ///quantum particle set
  ParticleSet& TargetPtcl;
  ///ionic system
  ParticleSet* SourcePtcl;

  /**  Helper vector for sorting bands
   */
  std::vector<std::unique_ptr<std::vector<BandInfo>>> FullBands;

  /// reader to use BsplineReader
  std::unique_ptr<BsplineReader> MixedSplineReader;

  ///This is true if we have the orbital derivatives w.r.t. the ion positions
  bool HaveOrbDerivs;
  ///root XML node with href, sort, tilematrix, twistnum, source, precision,truncate,version
  xmlNodePtr XMLRoot;

  std::map<H5OrbSet, SPOSet*, H5OrbSet> SPOSetMap;

  ///constructor
  EinsplineSetBuilder(ParticleSet& p, const PSetMap& psets, Communicate* comm, xmlNodePtr cur);

  ///destructor
  ~EinsplineSetBuilder() override;

  /** initialize the Antisymmetric wave function for electrons
   * @param cur the current xml node
   */
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

  /** initialize with the existing SPOSet */
  std::unique_ptr<SPOSet> createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info) override;

  //////////////////////////////////////
  // HDF5-related data  and functions //
  //////////////////////////////////////
  hdf_archive H5File;
  std::filesystem::path H5FileName;
  // HDF5 orbital file version
  typedef enum
  {
    QMCPACK,
    ESHDF
  } FormatType;
  FormatType Format;
  TinyVector<int, 3> Version;
  std::string parameterGroup, ionsGroup, eigenstatesGroup;
  std::vector<int> Occ;
  bool ReadOrbitalInfo(bool skipChecks = false);
  bool ReadOrbitalInfo_ESHDF(bool skipChecks = false);
  void BroadcastOrbitalInfo();
  bool CheckLattice();

  /** read gvectors for each twist
   * @return true, if psi_g is found
   */
  bool ReadGvectors_ESHDF();

  Tensor<double, OHMMS_DIM> Lattice, RecipLattice, LatticeInv, SuperLattice, GGt;
  UnitCellType SuperCell, PrimCell, PrimCellInv;
  int NumBands, NumElectrons, NumSpins, NumTwists;
  int MaxNumGvecs;
  double MeshFactor;
  RealType MatchingTol;
  TinyVector<int, 3> MeshSize;
  std::vector<std::vector<TinyVector<int, 3>>> Gvecs;

  Vector<int> IonTypes;
  Vector<TinyVector<double, OHMMS_DIM>> IonPos;
  // mapping the ions in the supercell to the primitive cell
  std::vector<int> Super2Prim;

  /////////////////////////////
  // Twist angle information //
  /////////////////////////////
  // The "true" twist number after analyzing twistnum, twist XML input and h5
  int twist_num_;
  // primitive cell k-points from DFT calculations
  std::vector<TinyVector<double, OHMMS_DIM>> primcell_kpoints;
  // primitive cell to supercell tiling matrix
  Tensor<int, OHMMS_DIM> TileMatrix;
  // This vector stores which twist indices will be used by this clone
  std::vector<TinyVector<int, OHMMS_DIM>> UseTwists;
  std::vector<int> IncludeTwists, DistinctTwists;
  /// if false, splines are conceptually complex valued
  bool use_real_splines_;
  int NumDistinctOrbitals;
  // This is true if the corresponding twist in DistinctTwists should
  // should be used to generate two distinct orbitals from the real and
  // imaginary parts.
  std::vector<bool> MakeTwoCopies;
  // This maps a 3-integer twist index into the twist number in the file
  std::map<TinyVector<int, OHMMS_DIM>, int, Int3less> TwistMap;

  bool TwistPair(PosType a, PosType b) const;
  void TileIons();
  void OccupyBands(int spin, int sortBands, int numOrbs, bool skipChecks = false);
  void OccupyBands_ESHDF(int spin, int sortBands, int numOrbs);

  ////////////////////////////////
  // Atomic orbital information //
  ////////////////////////////////
  struct CenterInfo
  {
    std::vector<int> lmax, spline_npoints, GroupID;
    std::vector<double> spline_radius, cutoff, inner_cutoff, non_overlapping_radius;
    std::vector<TinyVector<double, OHMMS_DIM>> ion_pos;
    int Ncenters;

    CenterInfo() : Ncenters(0){};

    void resize(int ncenters)
    {
      Ncenters = ncenters;
      lmax.resize(ncenters, -1);
      spline_npoints.resize(ncenters, -1);
      GroupID.resize(ncenters, 0);
      spline_radius.resize(ncenters, -1.0);
      inner_cutoff.resize(ncenters, -1.0);
      non_overlapping_radius.resize(ncenters, -1.0);
      cutoff.resize(ncenters, -1.0);
      ion_pos.resize(ncenters);
    }
  } AtomicCentersInfo;

  // This returns the path in the HDF5 file to the group for orbital
  // with twist ti and band bi
  std::string OrbitalPath(int ti, int bi);

  /////////////////////////////////////////////////////////////
  // Information to avoid storing the same orbitals twice in //
  // spin-restricted calculations.                           //
  /////////////////////////////////////////////////////////////
  int LastSpinSet, NumOrbitalsRead;

  std::string occ_format;
  int particle_hole_pairs;
  bool makeRotations;

protected:
  /** broadcast SortBands
   * @param N number of state
   * @param root true if it is the i/o node
   */
  void bcastSortBands(int splin, int N, bool root);

  /** a specific but clean code path in createSPOSetFromXML, for PBC, double, ESHDF
   * @param cur the current xml node
   */
  void set_metadata(int numOrbs,
                    int twist_num_inp,
                    const TinyVector<double, OHMMS_DIM>& twist_inp,
                    bool skipChecks = false);

  /** analyze twists of orbitals in h5 and determinine twist_num_
   * @param twist_num_inp twistnum XML input
   * @param twist_inp twst XML input
   */
  void AnalyzeTwists2(const int twist_num_inp, const TinyVector<double, OHMMS_DIM>& twist_inp);

  /// twistnum_inp == -9999 to indicate no given input after parsing XML
  static constexpr int TWISTNUM_NO_INPUT = -9999;
  /// twist_inp[i] <= -9999 to indicate no given input after parsing XML
  static constexpr double TWIST_NO_INPUT = -9999;
};

} // namespace qmcplusplus


#endif
