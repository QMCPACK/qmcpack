//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_EINSPLINE_SET_BUILDER_H
#define QMCPLUSPLUS_EINSPLINE_SET_BUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/EinsplineSet.h"
#include "Numerics/HDFNumericAttrib.h"
#include <map>

class Communicate;

namespace qmcplusplus {

  // Helper needed for TwistMap
  struct Int3less
  {
    bool operator()(const TinyVector<int,3>& a, const TinyVector<int,3> &b) const
    {
      if (a[0] > b[0]) return false;
      if (a[0] < b[0]) return true;
      if (a[1] > b[1]) return false;
      if (a[1] < b[1]) return true;
      if (a[2] > b[2]) return false;
      if (a[2] < b[2]) return true;	
      return false;
    }
  };
  struct Int4less
  { 
    bool operator()(const TinyVector<int,4>& a, const TinyVector<int,4>&b) 
      const
    {
      for (int i=0; i<4; i++) {
	if (a[i] > b[i]) return false;
	if (a[i] < b[i]) return true;
      }
      return false;
    }
  };
  

  struct H5OrbSet {
    string FileName;
    int SpinSet;
    int NumOrbs;
    bool operator()(const H5OrbSet &a, const H5OrbSet &b)
    {
      if (a.FileName == b.FileName) {
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
    H5OrbSet (string name, int spinSet, int numOrbs) :
      FileName(name), SpinSet(spinSet), NumOrbs(numOrbs)
    { }
    H5OrbSet() 
    { }
  };

  struct BandInfo {
    int TwistIndex, BandIndex, Spin;
    double Energy;
    // This is true if we should make distinct copies 
    // represeninting a +k, -k pair
    bool MakeTwoCopies;
    // True if this state is a core state
    bool IsCoreState;
    inline bool operator<(BandInfo other) const
    { return Energy < other.Energy; }
  };

  class EinsplineSetBuilder : public BasisSetBuilder {
  public:
    //////////////////////
    // Type definitions //
    //////////////////////
    typedef map<string,ParticleSet*> PtclPoolType;

    EinsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur);
    
    ~EinsplineSetBuilder();
    
    bool put (xmlNodePtr cur);

        /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     */
    SPOSetBase* createSPOSet(xmlNodePtr cur);
    
  protected:
    // Type definitions
    typedef CrystalLattice<RealType,OHMMS_DIM> UnitCellType;

    // Helper vector for sorting bands
    std::vector<BandInfo> SortBands;

    // The actual orbital set we're building
    EinsplineSet *OrbitalSet, *LastOrbitalSet;
    typedef EinsplineOrb<complex<double>,OHMMS_DIM> OrbType;
    // The map key is (spin, twist, band, center)
    static std::map<TinyVector<int,4>,OrbType*,Int4less> OrbitalMap;
    
    static std::map<H5OrbSet,multi_UBspline_3d_d*,H5OrbSet> ExtendedMap_d;
    static std::map<H5OrbSet,multi_UBspline_3d_z*,H5OrbSet> ExtendedMap_z;

    xmlNodePtr XMLRoot;

    //////////////////////////////////////
    // HDF5-related data  and functions //
    //////////////////////////////////////
    hid_t H5FileID;
    string H5FileName;
    // HDF5 orbital file version
    TinyVector<int,2> Version;
    string parameterGroup, ionsGroup, eigenstatesGroup;
    bool HaveLocalizedOrbs;
    bool ReadOrbitalInfo ();
    void BroadcastOrbitalInfo();


    Tensor<double,OHMMS_DIM> Lattice, RecipLattice, LatticeInv, SuperLattice;
    UnitCellType SuperCell, PrimCell, PrimCellInv;
    int NumBands, NumElectrons, NumSpins, NumTwists, 
      NumCoreStates;
    Vector<int> IonTypes;
    Vector<TinyVector<double,OHMMS_DIM> > IonPos;
    /////////////////////////////
    // Twist angle information //
    /////////////////////////////
    // This stores which "true" twist number I am using
    int TwistNum;
    std::vector<PosType> TwistAngles;
    TinyVector<int,OHMMS_DIM> TileFactor;
    Tensor<int,OHMMS_DIM> TileMatrix;
    TinyVector<int,OHMMS_DIM> TwistMesh;
    // This vector stores which twist indices will be used by this
    // clone 
    vector<TinyVector<int,OHMMS_DIM> > UseTwists;
    vector<int> IncludeTwists, DistinctTwists;
    int NumDistinctOrbitals, NumCoreOrbs, NumValenceOrbs;
    // This is true if the corresponding twist in DistinctTwists should
    // should be used to generate two distinct orbitals from the real and 
    // imaginary parts.
    vector<bool> MakeTwoCopies;
    inline bool TwistPair (PosType a, PosType b);
    // This maps a 3-integer twist index into the twist number in the file
    map <TinyVector<int,OHMMS_DIM>,int,Int3less> TwistMap;
    void AnalyzeTwists();
    void AnalyzeTwists2();
    void TileIons();
    void OccupyBands(int spin, bool sortBands);
    void ReadBands(int spin, EinsplineSetLocal* orbitalSet);
    void ReadBands(int spin, EinsplineSetExtended<complex<double> >* orbitalSet);
    void ReadBands(int spin, EinsplineSetExtended<        double  >* orbitalSet);
    void CopyBands(int numOrbs);
    
    /////////////////////////////
    // Muffin-tin information  //
    /////////////////////////////
    int NumMuffinTins;
    std::vector<RealType> MT_APW_radii;
    std::vector<Vector<double> > MT_APW_rgrids;
    std::vector<int> MT_APW_lmax;
    std::vector<int> MT_APW_num_radial_points;
    std::vector<PosType> MT_centers;

    // This returns the path in the HDF5 file to the group for orbital
    // with twist ti and band bi
    string OrbitalPath   (int ti, int bi);
    string CoreStatePath (int ti, int bi);
    string MuffinTinPath (int ti, int bi, int tin);

    /////////////////////////////////////////////////////////////
    // Information to avoid storing the same orbitals twice in //
    // spin-restricted calculations.                           //
    /////////////////////////////////////////////////////////////
    int LastSpinSet, NumOrbitalsRead;
  }; 
}


#endif
