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
  }

  SPOSetBase*
  EinsplineSetBuilder::createSPOSet(xmlNodePtr cur) {
    OhmmsAttributeSet attribs;
    int numOrbs = 0;
    attribs.add (H5FileName, "href");
    attribs.put (XMLRoot);
    attribs.add (Tile,       "tile");
    attribs.add (numOrbs,    "size");
    attribs.add (numOrbs,    "norbs");
    attribs.put (cur);
    if (numOrbs == 0) {
      app_error() << "You must specify the number of orbitals in the input file.\n";
      abort();
    }
    else 
      cerr << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";

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
	     "  Lattice = \n    [ %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f ]\n", 
	     Lattice(0,0), Lattice(0,1), Lattice(0,2), 
	     Lattice(1,0), Lattice(1,1), Lattice(1,2), 
	     Lattice(2,0), Lattice(2,1), Lattice(2,2));
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	LatticeInv(i,j) = RecipLattice(i,j)/(2.0*M_PI);
    HDFAttribIO<int> h_NumBands(NumBands), h_NumElectrons(NumElectrons), 
      h_NumSpins(NumSpins), h_NumTwists(NumTwists);
    h_NumBands.read     (H5FileID, "/parameters_0/num_bands");
    h_NumElectrons.read (H5FileID, "/parameters_0/num_electrons");
    h_NumSpins.read     (H5FileID, "/parameters_0/num_spins");
    h_NumTwists.read    (H5FileID, "/parameters_0/num_twists");
    fprintf (stderr, "  bands = %d, elecs = %d, spins = %d, twists = %d\n",
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

    ///////////////////////////
    // Read the twist angles //
    ///////////////////////////
    TwistAngles.resize(NumTwists);
    for (int ti=0; ti<NumTwists; ti++) {
      ostringstream path;
      path << "eigenstates_3/twist_" << ti << "/twist_angle";
      HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
      h_Twist.read (H5FileID, path.str().c_str());
      fprintf (stderr, "  Found twist angle (%5.3f, %5.3f, %5.3f)\n", 
	       TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    }
    AnalyzeTwists();
    fprintf (stderr, "  Found a valid %dx%dx%x twist mesh.\n", 
	     TwistMesh[0], TwistMesh[1], TwistMesh[2]);

    //////////////////////////////////
    // Create the OrbitalSet object //
    //////////////////////////////////
    // HACK HACK HACK, for now
    bool localized = false;
    if (localized)
      OrbitalSet = new EinsplineSetLocalized;
    else
      OrbitalSet = new EinsplineSetExtended;

    /////////////////////////
    // Setup internal data //
    /////////////////////////

    // Lattice information
    OrbitalSet->TileFactor = Tile;
    OrbitalSet->Tiling = Tile[0]!=1 || Tile[1]!=1 || Tile[2]!=1;
    Tensor<double,3> superLattice;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	superLattice(i,j) = (double)Tile[i]*Lattice(i,j);
    OrbitalSet->PrimLattice  = Lattice;
    OrbitalSet->SuperLattice = superLattice;
        
    // For now, just try reading k=0 data
    OrbitalSet->Orbitals.resize(NumBands);
    for (int bi=0; bi<NumBands; bi++) {
      ostringstream groupPath;
      groupPath << "/eigenstates_3/twist_0/band_" << bi << "/";
      OrbitalSet->Orbitals[bi].read(H5FileID, groupPath.str());
    }
    OrbitalSet->BasisSetSize   = NumBands;
    OrbitalSet->setOrbitalSetSize (numOrbs);

    return OrbitalSet;
  }
  
  
  // This function analyzes the twist vectors to see if they lay on a
  // valid k-point mesh.  It flags errors an aborts if they do not.
  // As a side-effect, it sets up TwistMap, which maps [ix,iy,iz] into
  // a single integer twist index.
  void
  EinsplineSetBuilder::AnalyzeTwists()
  {
    PosType minTwist(TwistAngles[0]), maxTwist(TwistAngles[0]);
    for (int ti=0; ti<NumTwists; ti++) 
      for (int i=0; i<3; i++) {
	minTwist[i] = min (TwistAngles[ti][i], minTwist[i]);
	maxTwist[i] = max (TwistAngles[ti][i], maxTwist[i]);
      }
    // The difference between maxTwist and minTwist should be of the
    // form (n-1)/n.  Therefore, we determine n by
    PosType nf;
    nf[0] = -1.0/((maxTwist[0]-minTwist[0]) -1.0);
    nf[1] = -1.0/((maxTwist[1]-minTwist[1]) -1.0);
    nf[2] = -1.0/((maxTwist[2]-minTwist[2]) -1.0);
    
    bool meshOK = true;
    // Make sure they are close to integers
    meshOK = meshOK && (std::fabs(nf[0] - round(nf[0]))<1.0e-6);
    meshOK = meshOK && (std::fabs(nf[1] - round(nf[1]))<1.0e-6);
    meshOK = meshOK && (std::fabs(nf[2] - round(nf[2]))<1.0e-6);
    if (!meshOK) {
      app_error() << "It appears that the twist angles in file " 
		  << H5FileName << " do not form a valid mesh.  Aborting.\n";
      abort();
    }
    
    TinyVector<int,3> n((int)round(nf[0]), (int)round(nf[1]), (int)round(nf[2]));
    TwistMesh = n;
    
    // Now, make sure we have all the k-points in the lattice
    PosType twist;
    for (int ix=0; ix<n[0]; ix++) 
      for (int iy=0; iy<n[1]; iy++)
	for (int iz=0; iz<n[2]; iz++) {
	  twist[0] = 
	    minTwist[0] + (double)ix/(double)(n[0]-1)*(maxTwist[0]-minTwist[0]);
	  twist[1] = 
	    minTwist[1] + (double)iy/(double)(n[1]-1)*(maxTwist[1]-minTwist[1]);
	  twist[2] = 
	    minTwist[2] + (double)iz/(double)(n[2]-1)*(maxTwist[2]-minTwist[2]);
	  bool twistFound = false;
	  for (int ti=0; ti<NumTwists; ti++) {
	    PosType diff = TwistAngles[ti]-twist;
	    if (dot(diff,diff)<1.0e-8) {
	      twistFound = true;
	      TinyVector<int,3> tindex (ix, iy, iz);	      
	      TwistMap[tindex] = ti;
	    }
	  }
	  if (!twistFound) {
	    fprintf (stderr, "Missing twist vector (%8.4f, %8.4f, %8.4f) "
		     "in CheckkPointMesh.\n", twist[0], twist[1], twist[2]);
	    abort();
	  }
	}
    
    // If we got this far, we have a valid twist mesh.  Now check to
    // see if the mesh is commensurate with the titling factor
    if (((TwistMesh[0] % Tile[0]) != 0) || ((TwistMesh[1] % Tile[1]) != 0) ||
	((TwistMesh[2] % Tile[2]) != 0)) {
      app_error() << "The tiling factor, " << Tile[0] << "x" << Tile[1] << "x" << Tile[2] 
		  << " is not commensurate with the k-point mesh, " 
		  << TwistMesh[0] << "x" << TwistMesh[1] << "x" << TwistMesh[2] << ".\n";
      abort();
    }
  }
}


