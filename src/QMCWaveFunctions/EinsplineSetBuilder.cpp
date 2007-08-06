//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim
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

#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"

namespace qmcplusplus {
  EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur) 
    : XMLRoot(cur), Tile(1,1,1)
  {
  }

  EinsplineSetBuilder::~EinsplineSetBuilder()
  {

  }


  bool 
  EinsplineSetBuilder::put(xmlNodePtr cur) {
    string hdfName;
    OhmmsAttributeSet attribs;
    attribs.add (hdfName, "href");
    attribs.put(XMLRoot);
    cerr << "In EinsplineSetBuilder::put.  href=\"" << hdfName << "\".\n";
  }

  SPOSetBase*
  EinsplineSetBuilder::createSPOSet(xmlNodePtr cur) {
    OhmmsAttributeSet attribs;
    attribs.add (H5FileName, "href");
    attribs.add (Tile, "tile");
    attribs.put (XMLRoot);
    cerr << "In EinsplineSetBuilder::createSPOSet().  href=\"" 
	 << H5FileName << "\".\n";
    H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    if (H5FileID < 0) {
      app_error() << "Could not open HDF5 file \"" << H5FileName 
		  << "\" in EinsplineSetBuilder::createSPOSet.  Aborting.\n";
      abort();
    }

    //////////////////////////////////////////////////
    // Read basic parameters from the orbital file. //
    //////////////////////////////////////////////////
    HDFAttribIO<Tensor<double,3> > h_Lattice(Lattice), h_RecipLattice(RecipLattice);
    h_Lattice.read      (H5FileID, "/parameters_0/lattice");
    h_RecipLattice.read (H5FileID, "/parameters_0/reciprocal_lattice");
    fprintf (stderr, 
	     "Lattice = \n   [ %8.5f %8.5f %8.5f\n"
	     "     %8.5f %8.5f %8.5f\n"
	     "     %8.5f %8.5f %8.5f ]\n", Lattice(0,0), Lattice(0,1), Lattice(0,2),
	     Lattice(1,0), Lattice(1,1), Lattice(1,2), Lattice(2,0), Lattice(2,1), Lattice(2,2));
    HDFAttribIO<int> h_NumBands(NumBands), h_NumElectrons(NumElectrons), 
      h_NumSpins(NumSpins), h_NumTwists(NumTwists);
    h_NumBands.read     (H5FileID, "/parameters_0/num_bands");
    h_NumElectrons.read (H5FileID, "/parameters_0/num_electrons");
    h_NumSpins.read     (H5FileID, "/parameters_0/num_spins");
    h_NumTwists.read    (H5FileID, "/parameters_0/num_twists");
    fprintf (stderr, "bands = %d, elecs = %d, spins = %d, twists = %d\n",
	     NumBands, NumElectrons, NumSpins, NumTwists);
    if (Tile[0]!=1 || Tile[1]!=1 || Tile[2]!=1)
      fprintf (stderr, "  Using a %dx%dx%d tiling factor.\n", 
	       Tile[0], Tile[1], Tile[2]);

    //////////////////////////////////
    // Read ion types and locations //
    //////////////////////////////////
    HDFAttribIO<Vector<int> >                 h_IonTypes(IonTypes);
    HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
    h_IonTypes.read (H5FileID, "/ions_2/atom_types");
    h_IonPos.read   (H5FileID, "/ions_2/pos");

    
  }
}
