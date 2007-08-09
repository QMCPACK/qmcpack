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
#include <vector>

namespace qmcplusplus {
  EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur) 
    : XMLRoot(cur), TileFactor(1,1,1), TwistNum(0)
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
    attribs.add (TileFactor, "tile");
    attribs.put (XMLRoot);
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
    // Check the version
    HDFAttribIO<TinyVector<int,2> > h_Version(Version);
    h_Version.read (H5FileID, "/version");
    string parameterGroup, ionsGroup, eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) {
      parameterGroup  = "/parameters_0";
      ionsGroup       = "/ions_2";
      eigenstatesGroup = "/eigenstates_3";
    }
    else if (Version[0]==0 && Version[1]==20) {
      parameterGroup  = "/parameters";
      ionsGroup       = "/ions";
      eigenstatesGroup = "/eigenstates";
    }
    else {
      cerr << "Unknown HDF5 orbital file version " 
	   << Version[0] << "." << Version[1] << "\n";
      abort();
    }
    fprintf (stderr, "  HDF5 orbital file version %d.%d.\n", Version[0], Version[1]);
    HDFAttribIO<Tensor<double,3> > h_Lattice(Lattice), h_RecipLattice(RecipLattice);
    h_Lattice.read      (H5FileID, (parameterGroup+"/lattice").c_str());
    
    h_RecipLattice.read (H5FileID, (parameterGroup+"/reciprocal_lattice").c_str());
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
    h_NumBands.read     (H5FileID, (parameterGroup+"/num_bands").c_str());
    h_NumElectrons.read (H5FileID, (parameterGroup+"/num_electrons").c_str());
    h_NumSpins.read     (H5FileID, (parameterGroup+"/num_spins").c_str());
    h_NumTwists.read    (H5FileID, (parameterGroup+"/num_twists").c_str());
    fprintf (stderr, "  bands = %d, elecs = %d, spins = %d, twists = %d\n",
	     NumBands, NumElectrons, NumSpins, NumTwists);
    if (TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1)
      fprintf (stderr, "  Using a %dx%dx%d tiling factor.\n", 
	       TileFactor[0], TileFactor[1], TileFactor[2]);

    //////////////////////////////////
    // Read ion types and locations //
    //////////////////////////////////
    HDFAttribIO<Vector<int> >                 h_IonTypes(IonTypes);
    HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
    h_IonTypes.read (H5FileID, (ionsGroup+"/atom_types").c_str());
    h_IonPos.read   (H5FileID, (ionsGroup+"/pos").c_str());

    ///////////////////////////
    // Read the twist angles //
    ///////////////////////////
    TwistAngles.resize(NumTwists);
    for (int ti=0; ti<NumTwists; ti++) {
      ostringstream path;
      if ((Version[0]==0 && Version[1]==11) || NumTwists > 0)
	path << eigenstatesGroup << "/twist_" << ti << "/twist_angle";
      else
	path << eigenstatesGroup << "/twist/twist_angle";
      HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
      h_Twist.read (H5FileID, path.str().c_str());
      fprintf (stderr, "  Found twist angle (%5.3f, %5.3f, %5.3f)\n", 
	       TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    }
    AnalyzeTwists();
    fprintf (stderr, "  Found a valid %dx%dx%d twist mesh.\n", 
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
    OrbitalSet->TileFactor = TileFactor;
    OrbitalSet->Tiling = TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1;
    Tensor<double,3> superLattice;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	superLattice(i,j) = (double)TileFactor[i]*Lattice(i,j);
    OrbitalSet->PrimLattice  = Lattice;
    OrbitalSet->SuperLattice = superLattice;
    OrbitalSet->GGt=dot(OrbitalSet->PrimLattice.G,
			transpose(OrbitalSet->PrimLattice.G));
        
//     OrbitalSet->Orbitals.resize(UseTwists.size()*NumBands);
//     for (int ti=0; ti<UseTwists.size(); ti++) {
//       int tindex = TwistMap[UseTwists[ti]];
//       for (int bi=0; bi<NumBands; bi++) {
// 	ostringstream groupPath;
// 	if ((Version[0]==0 && Version[1]==11) || NumTwists > 0)
// 	  groupPath << eigenstatesGroup << "/twist_" 
// 		    << tindex << "/band_" << bi << "/";
// 	else
// 	  groupPath << eigenstatesGroup << "/twist/band_" << bi << "/";
// 	OrbitalSet->Orbitals[bi].read(H5FileID, groupPath.str());
//       }
//     }
    
    OrbitalSet->setOrbitalSetSize (numOrbs);
    // Now, figure out occupation for the bands and read them
    OccupyAndReadBands();
    OrbitalSet->BasisSetSize   = numOrbs;

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
    if (((TwistMesh[0] % TileFactor[0]) != 0) || 
	((TwistMesh[1] % TileFactor[1]) != 0) ||
	((TwistMesh[2] % TileFactor[2]) != 0)) {
      app_error() << "The tiling factor, " 
		  << TileFactor[0] << "x" << TileFactor[1] << "x" << TileFactor[2] 
		  << " is not commensurate with the k-point mesh, " 
		  << TwistMesh[0] << "x" << TwistMesh[1] << "x" << TwistMesh[2] << ".\n";
      abort();
    }

    TinyVector<int,3> untiledMesh (TwistMesh[0]/TileFactor[0], 
				   TwistMesh[1]/TileFactor[1], 
				   TwistMesh[2]/TileFactor[2]);
    // Finally, let's decide which twist vectors we're supposed to
    // read
    fprintf (stderr, "  After untiling by %dx%dx%d, we are left with a %dx%dx%d k-point mesh.\n",
	     TileFactor[0], TileFactor[1], TileFactor[2],
	     untiledMesh[0], untiledMesh[1], untiledMesh[2]);
    TinyVector<int,3> offset;
    offset[0] = TwistNum/(untiledMesh[2]*untiledMesh[1]);
    offset[1] = (TwistNum%(untiledMesh[2]*untiledMesh[1])) / untiledMesh[1];
    offset[2] = (TwistNum%(untiledMesh[2]*untiledMesh[1])) % untiledMesh[1];
//     map<TinyVector<int,3>, int>::iterator iter;
//     for (iter = TwistMap.begin(); iter!=TwistMap.end(); iter++)
//       cerr << "TwistMap = " << (*iter).first 
// 	   << ", " << (*iter).second << endl;
    
    fprintf (stderr, "  Including twist vectors:\n");
    UseTwists.clear();
    for (int tx=0; tx<TileFactor[0]; tx++)
      for (int ty=0; ty<TileFactor[1]; ty++)
	for (int tz=0; tz<TileFactor[2]; tz++) {
	  TinyVector<int,3> tIndex;
	  tIndex = offset;
	  tIndex[0] += tx*untiledMesh[0];
	  tIndex[1] += ty*untiledMesh[1];
	  tIndex[2] += tz*untiledMesh[2];
	  UseTwists.push_back(tIndex);
	  int ti = TwistMap[tIndex];
 	  fprintf (stderr, "tIndex = (%d, %d, %d)  ti = %d\n", 
 		   tIndex[0], tIndex[1], tIndex[2], ti);
	    
	  fprintf (stderr, "    (%6.3f %6.3f %6.3f)\n", 
		   TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
	}
  }
  

  class BandInfo {
  public:
    int TwistIndex, BandIndex;
    double Energy;
    inline bool operator<(BandInfo other) const
    { return Energy < other.Energy; }
  };

  void
  EinsplineSetBuilder::OccupyAndReadBands()
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    std::vector<BandInfo> SortBands;
    for (int ti=0; ti<UseTwists.size(); ti++) {
      int tindex = TwistMap[UseTwists[ti]];
      for (int bi=0; bi<NumBands; bi++) {
	BandInfo band;
	band.TwistIndex = tindex;
	band.BandIndex  = bi;
	
	// Read eigenenergy from file
	ostringstream ePath;
	if ((Version[0]==0 && Version[1]==11) || NumTwists > 0)
	  ePath << eigenstatesGroup << "/twist_" 
		    << tindex << "/band_" << bi << "/eigenvalue";
	else
	  ePath << eigenstatesGroup << "/twist/band_" << bi << "/eigenvalue";
	
	HDFAttribIO<double> h_energy(band.Energy);
	h_energy.read (H5FileID, ePath.str().c_str());
	SortBands.push_back(band);
      }
    }
    // Now sort the bands by energy
    sort (SortBands.begin(), SortBands.end());
    // Read in the occupied bands
    OrbitalSet->Orbitals.resize(OrbitalSet->getOrbitalSetSize());
    cerr << "Orbitals size = " << OrbitalSet->Orbitals.size() << endl;
    for (int i=0; i<OrbitalSet->getOrbitalSetSize(); i++) {
      int ti   = SortBands[i].TwistIndex;
      int bi   = SortBands[i].BandIndex;
      double e = SortBands[i].Energy;
      ostringstream groupPath;
      if ((Version[0]==0 && Version[1]==11) || NumTwists > 0)
	groupPath << eigenstatesGroup << "/twist_" 
		  << ti << "/band_" << bi << "/";
      else
	groupPath << eigenstatesGroup << "/twist/band_" << bi << "/";
      
      PosType twist, k;
      twist = TwistAngles[ti];
      Tensor<double,3> G = OrbitalSet->PrimLattice.G;
      k = 2.0*M_PI*(twist[0]*OrbitalSet->PrimLattice.Gv[0] +
		    twist[1]*OrbitalSet->PrimLattice.Gv[1] +
		    twist[2]*OrbitalSet->PrimLattice.Gv[2]);
	
//       k[0] = twist[0]*G(0,0) + twist[1]*G(0,1) + twist[2]*G(0,2);
//       k[1] = twist[0]*G(1,0) + twist[1]*G(1,1) + twist[2]*G(1,2);
//       k[2] = twist[0]*G(2,0) + twist[1]*G(2,1) + twist[2]*G(2,2);
      fprintf (stderr, "  ti=%d  bi=%d energy=%8.5f k=(%6.4f, %6.4f, %6.4f)\n", 
	       ti, bi, e, k[0], k[1], k[2]);
      
      OrbitalSet->Orbitals[i].kVec = k;
      OrbitalSet->Orbitals[i].read(H5FileID, groupPath.str());
    }
  }
}


