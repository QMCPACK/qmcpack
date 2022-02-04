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
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "Numerics/HDFNumericAttrib.h"
#include <map>

#define PW_COEFF_NORM_TOLERANCE 1e-6

class Communicate;

namespace qmcplusplus
{
///forward declaration of BsplineReaderBase
struct BsplineReaderBase;

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
  ///type of orbitals defined
  int OrbitalType;
  ///index for the spin set
  int SpinSet;
  ///number of orbitals that belong to this set
  int NumOrbs;
  ///name of the HDF5 file
  std::string FileName;
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

  H5OrbSet(const H5OrbSet& a) : SpinSet(a.SpinSet), NumOrbs(a.NumOrbs), FileName(a.FileName) {}
  H5OrbSet(std::string name, int spinSet, int numOrbs) : SpinSet(spinSet), NumOrbs(numOrbs), FileName(name) {}
  H5OrbSet() {}
};

/** EinsplineSet builder
 */
class EinsplineSetBuilder : public SPOSetBuilder
{
public:
  using PtclPoolType = std::map<std::string, ParticleSet*>;
  using UnitCellType = CrystalLattice<ParticleSet::Scalar_t, DIM>;

  ///reference to the particleset pool
  const PtclPoolType& ParticleSets;
  ///quantum particle set
  ParticleSet& TargetPtcl;
  ///ionic system
  ParticleSet* SourcePtcl;

  /**  Helper vector for sorting bands
   */
  std::vector<std::unique_ptr<std::vector<BandInfo>>> FullBands;

  /// reader to use BsplineReaderBase
  std::unique_ptr<BsplineReaderBase> MixedSplineReader;

  ///This is true if we have the orbital derivatives w.r.t. the ion positions
  bool HaveOrbDerivs;
  ///root XML node with href, sort, tilematrix, twistnum, source, precision,truncate,version
  xmlNodePtr XMLRoot;

  std::map<H5OrbSet, SPOSet*, H5OrbSet> SPOSetMap;

  ///constructor
  EinsplineSetBuilder(ParticleSet& p, const PtclPoolType& psets, Communicate* comm, xmlNodePtr cur);

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
  hid_t H5FileID;
  std::string H5FileName;
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
  bool HasCoreOrbs;
  bool ReadOrbitalInfo(bool skipChecks = false);
  bool ReadOrbitalInfo_ESHDF(bool skipChecks = false);
  void BroadcastOrbitalInfo();
  bool CheckLattice();

  /** read gvectors for each twist
   * @return true, if psi_g is found
   */
  bool ReadGvectors_ESHDF();

  /** set tiling properties of oset
   * @param oset spline-orbital engine to be initialized
   * @param numOrbs number of orbitals that belong to oset
   */
  template<typename SPE>
  inline void setTiling(SPE* oset, int numOrbs)
  {
    oset->TileFactor = TileFactor;
    oset->Tiling     = (TileFactor[0] * TileFactor[1] * TileFactor[2] != 1);
    oset->PrimLattice.set(Lattice);
    oset->SuperLattice.set(SuperLattice);
    oset->GGt = GGt;
    oset->setOrbitalSetSize(numOrbs);
  }


  Tensor<double, OHMMS_DIM> Lattice, RecipLattice, LatticeInv, SuperLattice, GGt;
  UnitCellType SuperCell, PrimCell, PrimCellInv;
  int NumBands, NumElectrons, NumSpins, NumTwists, NumCoreStates;
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
  std::vector<TinyVector<double, OHMMS_DIM>> TwistAngles;
  //     integer index of sym operation from the irreducible brillion zone
  std::vector<int> TwistSymmetry;
  //     number of twists equivalent to this one in the big DFT grid
  std::vector<int> TwistWeight;

  TinyVector<int, OHMMS_DIM> TileFactor;
  Tensor<int, OHMMS_DIM> TileMatrix;
  TinyVector<int, OHMMS_DIM> TwistMesh;
  // This vector stores which twist indices will be used by this
  // clone
  std::vector<TinyVector<int, OHMMS_DIM>> UseTwists;
  std::vector<int> IncludeTwists, DistinctTwists;
  /// if false, splines are conceptually complex valued
  bool use_real_splines_;
  int NumDistinctOrbitals, NumCoreOrbs, NumValenceOrbs;
  // This is true if the corresponding twist in DistinctTwists should
  // should be used to generate two distinct orbitals from the real and
  // imaginary parts.
  std::vector<bool> MakeTwoCopies;
  inline bool TwistPair(PosType a, PosType b);
  // This maps a 3-integer twist index into the twist number in the file
  std::map<TinyVector<int, OHMMS_DIM>, int, Int3less> TwistMap;
  void TileIons();
  void OccupyBands(int spin, int sortBands, int numOrbs, bool skipChecks = false);
  void OccupyBands_ESHDF(int spin, int sortBands, int numOrbs);

  void CopyBands(int numOrbs);

  /////////////////////////////
  // Muffin-tin information  //
  /////////////////////////////
  int NumMuffinTins;
  std::vector<double> MT_APW_radii;
  std::vector<Vector<double>> MT_APW_rgrids;
  std::vector<int> MT_APW_lmax;
  std::vector<int> MT_APW_num_radial_points;
  std::vector<TinyVector<double, OHMMS_DIM>> MT_centers;

  ////////////////////////////////
  // Atomic orbital information //
  ////////////////////////////////
  std::vector<AtomicOrbital<std::complex<double>>> AtomicOrbitals;

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
  std::string CoreStatePath(int ti, int bi);
  std::string MuffinTinPath(int ti, int bi, int tin);

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
   * @return true, if core is found
   */
  bool bcastSortBands(int splin, int N, bool root);

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
