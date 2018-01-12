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

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/BandInfo.h"
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "QMCWaveFunctions/EinsplineSet.h"
#include "Numerics/HDFNumericAttrib.h"
#include <map>

class Communicate;

namespace qmcplusplus
{

///forward declaraton of BsplineReaderBase
class BsplineReaderBase;

// Helper needed for TwistMap
struct Int3less
{
  bool operator()(const TinyVector<int,3>& a, const TinyVector<int,3> &b) const
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
  bool operator()(const TinyVector<int,4>& a, const TinyVector<int,4>&b)
  const
  {
    for (int i=0; i<4; i++)
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
  bool operator()(const H5OrbSet &a, const H5OrbSet &b) const
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

  H5OrbSet (const H5OrbSet &a) :
    FileName(a.FileName), SpinSet(a.SpinSet), NumOrbs(a.NumOrbs)
  { }
  H5OrbSet ( std::string name, int spinSet, int numOrbs) :
    FileName(name), SpinSet(spinSet), NumOrbs(numOrbs)
  { }
  H5OrbSet()
  { }
};

/** EinsplineSet builder
 */
class EinsplineSetBuilder : public BasisSetBuilder
{
public:

  typedef std::map<std::string,ParticleSet*> PtclPoolType;
  typedef CrystalLattice<ParticleSet::Scalar_t,DIM> UnitCellType;

  ///reference to the particleset pool
  PtclPoolType &ParticleSets;
  ///quantum particle set
  ParticleSet &TargetPtcl;
  ///ionic system
  ParticleSet *SourcePtcl;
  ///index for the ion-el distance table
  int myTableIndex;

  /**  Helper vector for sorting bands
   */
  //std::vector<BandInfo> SortBands;
  std::vector<std::vector<BandInfo>*> FullBands;

  /// reader to use BsplineReaderBase
  BsplineReaderBase *MixedSplineReader;

  ///This is true if we have the orbital derivatives w.r.t. the ion positions
  bool HaveOrbDerivs;
  ///root XML node with href, sort, tilematrix, twistnum, source, precision,truncate,version
  xmlNodePtr XMLRoot;

  ////static std::map<H5OrbSet,multi_UBspline_3d_d*,H5OrbSet> ExtendedMap_d;
  ////static std::map<H5OrbSet,multi_UBspline_3d_z*,H5OrbSet> ExtendedMap_z;
  ////static std::map<H5OrbSet,EinsplineSetExtended<double>*,H5OrbSet> ExtendedSetMap_d;
  //static std::map<H5OrbSet,SPOSetBase*,H5OrbSet> SPOSetMap;
  std::map<H5OrbSet,SPOSetBase*,H5OrbSet> SPOSetMap;

  ///constructor
  EinsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur);

  ///destructor
  ~EinsplineSetBuilder();

  /** process xml node to initialize the builder */
  bool put (xmlNodePtr cur);

  /** initialize the Antisymmetric wave function for electrons
   * @param cur the current xml node
   */
  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);

  /** a specific but clean code path in createSPOSetFromXML, for PBC, double, ESHDF
   * @param cur the current xml node
   */
  void set_metadata(int numOrbs, int TwistNum_inp);

  /** initialize with the existing SPOSet */
  SPOSetBase* createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input_info);

  //////////////////////////////////////
  // HDF5-related data  and functions //
  //////////////////////////////////////
  hid_t H5FileID;
  std::string H5FileName;
  // HDF5 orbital file version
  typedef enum {QMCPACK, ESHDF} FormatType;
  FormatType Format;
  TinyVector<int,3> Version;
  std::string parameterGroup, ionsGroup, eigenstatesGroup;
  std::vector<int> Occ;
  bool HasCoreOrbs;
  bool ReadOrbitalInfo ();
  bool ReadOrbitalInfo_ESHDF ();
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
    oset->Tiling = (TileFactor[0]*TileFactor[1]*TileFactor[2] != 1);
    oset->PrimLattice  = Lattice;
    oset->SuperLattice = SuperLattice;
    //oset->GGt=dot(transpose(oset->PrimLattice.G), oset->PrimLattice.G);
    oset->GGt=GGt;
    oset->setOrbitalSetSize (numOrbs);
    oset->BasisSetSize   = numOrbs;
  }


  Tensor<double,OHMMS_DIM> Lattice, RecipLattice, LatticeInv, SuperLattice, GGt;
  UnitCellType SuperCell, PrimCell, PrimCellInv;
  int NumBands, NumElectrons, NumSpins, NumTwists, NumCoreStates;
  int MaxNumGvecs;
  double MeshFactor;
  RealType BufferLayer;
  RealType MatchingTol;
  TinyVector<int,3> MeshSize;
  std::vector<std::vector<TinyVector<int,3> > > Gvecs;

  //fftw_plan FFTplan;
  //Array<std::complex<double>,3> FFTbox;

  Vector<int> IonTypes;
  Vector<TinyVector<double,OHMMS_DIM> > IonPos;
  // mapping the ions in the supercell to the primitive cell
  std::vector<int> Super2Prim;

  /////////////////////////////
  // Twist angle information //
  /////////////////////////////
  // This stores which "true" twist number I am using
  int TwistNum;
  TinyVector<double,OHMMS_DIM> givenTwist;
  std::vector<TinyVector<double,OHMMS_DIM> > TwistAngles;
//     integer index of sym operation from the irreducible brillion zone
  std::vector<int> TwistSymmetry;
//     number of twists equivalent to this one in the big DFT grid
  std::vector<int> TwistWeight;

  TinyVector<int,OHMMS_DIM> TileFactor;
  Tensor<int,OHMMS_DIM> TileMatrix;
  TinyVector<int,OHMMS_DIM> TwistMesh;
  // This vector stores which twist indices will be used by this
  // clone
  std::vector<TinyVector<int,OHMMS_DIM> > UseTwists;
  std::vector<int> IncludeTwists, DistinctTwists;
  bool UseRealOrbitals;
  int NumDistinctOrbitals, NumCoreOrbs, NumValenceOrbs;
  // This is true if the corresponding twist in DistinctTwists should
  // should be used to generate two distinct orbitals from the real and
  // imaginary parts.
  std::vector<bool> MakeTwoCopies;
  inline bool TwistPair (PosType a, PosType b);
  // This maps a 3-integer twist index into the twist number in the file
  std::map<TinyVector<int,OHMMS_DIM>,int,Int3less> TwistMap;
  //void AnalyzeTwists();
  void AnalyzeTwists2();
  void TileIons();
  void OccupyBands(int spin, int sortBands, int numOrbs);
  void OccupyBands_ESHDF(int spin, int sortBands, int numOrbs);

#if 0
  void ReadBands      (int spin, EinsplineSetExtended<std::complex<double> >* orbitalSet);
  void ReadBands_ESHDF(int spin, EinsplineSetExtended<std::complex<double> >* orbitalSet);
  void ReadBands      (int spin, EinsplineSetExtended<        double  >* orbitalSet);
  void ReadBands_ESHDF(int spin, EinsplineSetExtended<        double  >* orbitalSet);
#endif

  void CopyBands(int numOrbs);

  /////////////////////////////
  // Muffin-tin information  //
  /////////////////////////////
  int NumMuffinTins;
  std::vector<double> MT_APW_radii;
  std::vector<Vector<double> > MT_APW_rgrids;
  std::vector<int> MT_APW_lmax;
  std::vector<int> MT_APW_num_radial_points;
  std::vector<TinyVector<double, OHMMS_DIM> > MT_centers;

  ////////////////////////////////
  // Atomic orbital information //
  ////////////////////////////////
  std::vector<AtomicOrbital<std::complex<double> > > AtomicOrbitals;

  struct CenterInfo
  {
    std::vector<int> lmax, spline_npoints, GroupID;
    std::vector<double> spline_radius, cutoff;
    std::vector<TinyVector<double,OHMMS_DIM> > ion_pos;
    int Ncenters;

    CenterInfo(): Ncenters(0) {};

    void resize(int ncenters)
    {
      Ncenters=ncenters;
      lmax.resize(ncenters, -1);
      spline_npoints.resize(ncenters, -1);
      GroupID.resize(ncenters);
      spline_radius.resize(ncenters, -1.0);
      cutoff.resize(ncenters, -1.0);
      ion_pos.resize(ncenters);
    }
  } AtomicCentersInfo;

  // This returns the path in the HDF5 file to the group for orbital
  // with twist ti and band bi
  std::string OrbitalPath   (int ti, int bi);
  std::string CoreStatePath (int ti, int bi);
  std::string MuffinTinPath (int ti, int bi, int tin);

  /////////////////////////////////////////////////////////////
  // Information to avoid storing the same orbitals twice in //
  // spin-restricted calculations.                           //
  /////////////////////////////////////////////////////////////
  int LastSpinSet, NumOrbitalsRead;

  std::string occ_format;
  RealType qafm;
  int particle_hole_pairs;
  bool makeRotations;
#if 0
  std::vector<RealType> rotationMatrix;
  std::vector<int> rotatedOrbitals;
  void RotateBands_ESHDF(int spin, EinsplineSetExtended<std::complex<double > >* orbitalSet);
  void RotateBands_ESHDF(int spin, EinsplineSetExtended<double>* orbitalSet);
#endif

  /** broadcast SortBands
   * @param N number of state
   * @param root true if it is the i/o node
   * @return true, if core is found
   */
  bool bcastSortBands(int splin, int N, bool root);

  int MyToken;
  inline void update_token(const char* f, int l, const char* msg) 
  {
    app_debug() << "TOKEN=" << MyToken << " " << msg << " " << f << " " << l << std::endl; MyToken++;
  }
  //inline void update_token(const char* f, int l, const char* msg) 
  //{}
};

}


#endif
