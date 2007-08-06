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

namespace qmcplusplus {
  
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
    // Helper needed for TwistMap
    struct Int3less
    {
      bool operator()(TinyVector<int,3> a, TinyVector<int,3> b) const
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

    // The actual orbital set we're building
    EinsplineSetBase *OrbitalSet;

    xmlNodePtr XMLRoot;
    hid_t H5FileID;
    string H5FileName;
    Tensor<double,OHMMS_DIM> Lattice, RecipLattice, LatticeInv;
    int NumBands, NumElectrons, NumSpins, NumTwists;
    Vector<int> IonTypes;
    Vector<TinyVector<double,OHMMS_DIM> > IonPos;
    // Twist angle information
    Vector<PosType> TwistAngles;
    TinyVector<int,OHMMS_DIM> Tile;
    TinyVector<int,OHMMS_DIM> TwistMesh;
    // This maps a 3-integer twist index into the twist number in the file
    map <TinyVector<int,OHMMS_DIM>,int,Int3less> TwistMap;
    void AnalyzeTwists();

  };

}


#endif
