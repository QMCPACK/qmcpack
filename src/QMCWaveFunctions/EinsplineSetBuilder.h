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

    // The actual orbital set we're building
    EinsplineSet *OrbitalSet, *LastOrbitalSet;
    typedef EinsplineOrb<complex<double>,OHMMS_DIM> OrbType;
    // The map key is (spin, twist, band, center)
    static std::map<TinyVector<int,4>,OrbType*,Int4less> OrbitalMap;

    xmlNodePtr XMLRoot;
    hid_t H5FileID;
    string H5FileName;
    // HDF5 orbital file version
    TinyVector<int,2> Version;

    Tensor<double,OHMMS_DIM> Lattice, RecipLattice, LatticeInv, SuperLattice;
    UnitCellType SuperCell, PrimCell, PrimCellInv;
    int NumBands, NumElectrons, NumSpins, NumTwists;
    Vector<int> IonTypes;
    Vector<TinyVector<double,OHMMS_DIM> > IonPos;
    /////////////////////////////
    // Twist angle information //
    /////////////////////////////
    // This stores which "true" twist number I am using
    int TwistNum;
    Vector<PosType> TwistAngles;
    TinyVector<int,OHMMS_DIM> TileFactor;
    Tensor<int,OHMMS_DIM> TileMatrix;
    TinyVector<int,OHMMS_DIM> TwistMesh;
    // This vector stores which twist indices will be used by this
    // clone 
    vector<TinyVector<int,OHMMS_DIM> > UseTwists;
    vector<int> IncludeTwists;
    // This maps a 3-integer twist index into the twist number in the file
    map <TinyVector<int,OHMMS_DIM>,int,Int3less> TwistMap;
    void AnalyzeTwists();
    void AnalyzeTwists2();
    void OccupyAndReadBands(int spin);
    void OccupyAndReadBands2(int spin, bool sortBands);
    void CopyBands(int numOrbs);

    /////////////////////////////////////////////////////////////
    // Information to avoid storing the same orbitals twice in //
    // spin-restricted calculations.                           //
    /////////////////////////////////////////////////////////////
    int LastSpinSet, NumOrbitalsRead;
  };

}


#endif
