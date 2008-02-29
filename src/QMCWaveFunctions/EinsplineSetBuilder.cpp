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
  std::map<TinyVector<int,4>,EinsplineSetBuilder::OrbType*,Int4less> 
  EinsplineSetBuilder::OrbitalMap;

  EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, 
      PtclPoolType& psets, xmlNodePtr cur) 
    : XMLRoot(cur), TileFactor(1,1,1), TwistNum(0), LastSpinSet(-1), NumOrbitalsRead(-1)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	TileMatrix(i,j) = 0;
  }

  EinsplineSetBuilder::~EinsplineSetBuilder()
  {

  }


  bool 
  EinsplineSetBuilder::put(xmlNodePtr cur) {
    string hdfName;
    OhmmsAttributeSet attribs;
    attribs.add (hdfName, "href");
    return attribs.put(XMLRoot);
  }

  SPOSetBase*
  EinsplineSetBuilder::createSPOSet(xmlNodePtr cur) {
    OhmmsAttributeSet attribs;
    int numOrbs = 0;
    bool sortBands = true;
    attribs.add (H5FileName, "href");
    attribs.add (TileFactor, "tile");
    attribs.add (sortBands, "sort");
    attribs.add (TileMatrix, "tilematrix");
    attribs.put (XMLRoot);
    attribs.add (numOrbs,    "size");
    attribs.add (numOrbs,    "norbs");
    attribs.put (cur);

    // The tiling can be set by a simple vector, (e.g. 2x2x2), or by a
    // full 3x3 matrix of integers.  If the tilematrix was not set in
    // the input file... 
    bool matrixNotSet = true;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	matrixNotSet = matrixNotSet && (TileMatrix(i,j) == 0);
    // then set the matrix to what may have been specified in the
    // tiling vector
    if (matrixNotSet) 
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	  TileMatrix(i,j) = (i==j) ? TileFactor(i) : 0;

    fprintf (stderr, " [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d ]\n",
	     TileMatrix(0,0), TileMatrix(0,1), TileMatrix(0,2),
	     TileMatrix(1,0), TileMatrix(1,1), TileMatrix(1,2),
	     TileMatrix(2,0), TileMatrix(2,1), TileMatrix(2,2));

    if (numOrbs == 0) {
      app_error() << "You must specify the number of orbitals in the input file.\n";
      abort();
    }
    else 
      cerr << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";

    ///////////////////////////////////////////////
    // Read occupation information from XML file //
    ///////////////////////////////////////////////
    cur = cur->children;
    int spinSet = -1;
    while (cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "occupation") {
	string occ_mode("ground");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(occ_mode,"mode");
        oAttrib.add(spinSet,"spindataset");
        oAttrib.put(cur);
	if(occ_mode != "ground") {
	  app_error() << "Only ground state occupation currently supported "
		      << "in EinsplineSetBuilder.\n";
	  abort();
	}
      }
      cur = cur->next;
    }


    
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
    fprintf (stderr, "  HDF5 orbital file version %d.%d\n", Version[0], Version[1]);
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
//     for (int i=0; i<3; i++)
//       for (int j=0; j<3; j++)
// 	superLattice(i,j) = (double)TileFactor[i]*Lattice(i,j);
    SuperLattice = dot(TileMatrix, Lattice);

    fprintf (stderr, 
	     "  Lattice = \n    [ %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f ]\n", 
	     Lattice(0,0), Lattice(0,1), Lattice(0,2), 
	     Lattice(1,0), Lattice(1,1), Lattice(1,2), 
	     Lattice(2,0), Lattice(2,1), Lattice(2,2));
    fprintf (stderr, 
	     "  SuperLattice = \n    [ %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f ]\n", 
	     SuperLattice(0,0), SuperLattice(0,1), SuperLattice(0,2), 
	     SuperLattice(1,0), SuperLattice(1,1), SuperLattice(1,2), 
	     SuperLattice(2,0), SuperLattice(2,1), SuperLattice(2,2));
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
      if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
	path << eigenstatesGroup << "/twist_" << ti << "/twist_angle";
      else
	path << eigenstatesGroup << "/twist/twist_angle";
      HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
      h_Twist.read (H5FileID, path.str().c_str());
      fprintf (stderr, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n", 
	       TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    }

    ////////////////////////////////////////////////////////////////
    // Determine whether we need comp //
    ////////////////////////////////////////////////////////////////
    bool HaveLocalizedOrbs = false;
    
    for (int ti=0; ti<NumTwists; ti++) {
      for (int bi=0; bi<NumBands; bi++) {
	double radius = 0.0;
	ostringstream path;
	if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
	  path << eigenstatesGroup << "/twist_" << ti << "/band_"
	       << bi << "/radius";
	else
	  path << eigenstatesGroup << "/twist/band_" << bi << "/radius";
	HDFAttribIO<double>  h_radius(radius);
	h_radius.read(H5FileID, path.str().c_str());
	HaveLocalizedOrbs = HaveLocalizedOrbs || (radius > 0.0);
      }
    }

    //////////////////////////////////////////////////////////////////////////////
    // Now, analyze the k-point mesh to figure out the what k-points are needed //
    //////////////////////////////////////////////////////////////////////////////
    AnalyzeTwists2();

    //////////////////////////////////
    // Create the OrbitalSet object //
    //////////////////////////////////
    if (HaveLocalizedOrbs) 
      OrbitalSet = new EinsplineSetLocal;
    else
      OrbitalSet = new EinsplineSetExtended<complex<double> >;

    /////////////////////////
    // Setup internal data //
    /////////////////////////
    
    // Lattice information
    OrbitalSet->TileFactor = TileFactor;
    OrbitalSet->Tiling = 
      TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1;
    PrimCell = Lattice;
    SuperCell = SuperLattice;
    OrbitalSet->PrimLattice  = Lattice;
    OrbitalSet->SuperLattice = SuperLattice;
    OrbitalSet->GGt=dot(OrbitalSet->PrimLattice.G,
			transpose(OrbitalSet->PrimLattice.G));
    OrbitalSet->setOrbitalSetSize (numOrbs);
    OrbitalSet->BasisSetSize   = numOrbs;
    
    
    if (HaveLocalizedOrbs) {
      EinsplineSetLocal *restrict orbitalSet = 
	dynamic_cast<EinsplineSetLocal*>(OrbitalSet);

#pragma omp critical(read_einspline_orbs)
      {
	if ((spinSet == LastSpinSet) && (numOrbs <= NumOrbitalsRead))
	  CopyBands(numOrbs);
	else {
	  // Now, figure out occupation for the bands and read them
	  OccupyBands(spinSet, sortBands);
	  ReadBands (spinSet, orbitalSet);
	}
      }
      // Now, store what we have read
      LastOrbitalSet = OrbitalSet;
      LastSpinSet = spinSet;
      NumOrbitalsRead = numOrbs;
    }
    // Otherwise, use EinsplineSetExtended
    else {
      EinsplineSetExtended<complex<double> > *restrict orbitalSet = 
	dynamic_cast<EinsplineSetExtended<complex<double> >*>(OrbitalSet);
      OccupyBands(spinSet, sortBands);
      ReadBands(spinSet, orbitalSet);
    }
    return OrbitalSet;
  }
  

  inline TinyVector<double,3>
  IntPart (TinyVector<double,3> twist)
  {
    return TinyVector<double,3> (round(twist[0]), round(twist[1]), round(twist[2]));
  }
  
  inline TinyVector<double,3>
  FracPart (TinyVector<double,3> twist)
  {
    return twist - IntPart (twist);
  }
  
  void
  EinsplineSetBuilder::AnalyzeTwists2()
  {
    Tensor<double,3> S;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	S(i,j) = (double)TileMatrix(i,j);
    
    vector<PosType> superFracs;
    // This holds to which supercell kpoint each primitive k-point belongs
    vector<int> superIndex;

    int numPrimTwists = TwistAngles.size();

    for (int ki=0; ki<numPrimTwists; ki++) {
      PosType primTwist = TwistAngles[ki];
      PosType superTwist = dot (S, primTwist);
      PosType kp = PrimCell.k_cart(primTwist);
      PosType ks = SuperCell.k_cart(superTwist);
      if (dot(ks-kp, ks-kp) > 1.0e-12) {
	app_error() << "Primitive and super k-points do not agree.  Error in coding.\n";
	abort();
      }
      PosType frac = FracPart (superTwist);
      bool found = false;
      for (int j=0; j<superFracs.size(); j++) {
	PosType diff = frac - superFracs[j];
	if (dot(diff,diff)<1.0e-12) {
	  found = true;
	  superIndex.push_back(j);
	}
      }
      if (!found) {
	superIndex.push_back(superFracs.size());
	superFracs.push_back(frac);
      }
    }
    int numSuperTwists = superFracs.size();
    cerr << "Found " << numSuperTwists << " distinct supercell twists.\n";

    // For each supercell twist, create a list of primitive twists which
    // belong to it.
    vector<vector<int> > superSets;
    superSets.resize(numSuperTwists);
    for (int ki=0; ki<numPrimTwists; ki++)
      superSets[superIndex[ki]].push_back(ki);

    for (int si=0; si<numSuperTwists; si++) {
      fprintf (stderr, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n",
	       si, superFracs[si][0], superFracs[si][1], superFracs[si][2]);
      fprintf (stderr, "  Using k-points: ");
      for (int i=0; i<superSets[si].size(); i++) 
	fprintf (stderr, " %d", superSets[si][i]);
      fprintf (stderr, "\n");
    }

    // Now check to see that each supercell twist has the right twists
    // to tile the primitive cell orbitals.
    int numTwistsNeeded = abs(TileMatrix.det());
    for (int si=0; si<numSuperTwists; si++) {
      // First make sure we have enough points
      if (superSets[si].size() != numTwistsNeeded) {
	fprintf (stderr, "Super twist %d should own %d k-points, but owns %d.\n",
		 si, numTwistsNeeded, superSets[si].size());
	abort();
      }
      // Now, make sure they are all distinct
      int N = superSets[si].size();
      for (int i=0; i<N; i++) {
	PosType twistPrim_i  = TwistAngles[superSets[si][i]];
	PosType twistSuper_i = dot (S, twistPrim_i);
	PosType superInt_i   = IntPart (twistSuper_i);
	for (int j=i+1; j<N; j++) {
	  PosType twistPrim_j  = TwistAngles[superSets[si][j]];
	  PosType twistSuper_j = dot (S, twistPrim_j);
	  PosType superInt_j   = IntPart (twistSuper_j);
	  if (dot(superInt_i-superInt_j, superInt_i-superInt_j) < 1.0e-6) {
	    cerr << "Identical k-points detected in super twist set "
		 << si << endl;
	    abort();
	  }
	}
      }
    }
    if (TwistNum >= superSets.size()) {
      app_error() << "Trying to use supercell twist " << TwistNum
		  << " when only " << superSets.size() << " sets exist.\n"
		  << "Please select a twist number between 0 and "
		  << superSets.size()-1 << ".\n";
      abort();
    }

    // Finally, record which k-points to include on this group of
    // processors, which have been assigned supercell twist TwistNum
    IncludeTwists.clear();
    for (int i=0; i<superSets[TwistNum].size(); i++)
      IncludeTwists.push_back(superSets[TwistNum][i]);
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
    // see if the mesh is commensurate with the tiling factor
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
  

  void
  EinsplineSetBuilder::OccupyBands(int spin, bool sortBands)
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    SortBands.clear();
    for (int ti=0; ti<IncludeTwists.size(); ti++) {
      int tindex = IncludeTwists[ti];
      for (int bi=0; bi<NumBands; bi++) {
	BandInfo band;
	band.TwistIndex = tindex;
	band.BandIndex  = bi;
	
	// Read eigenenergy from file
	ostringstream ePath, sPath;
	if ((Version[0]==0 && Version[1]==11) || NumTwists > 1) {
	  ePath << eigenstatesGroup << "/twist_" 
		    << tindex << "/band_" << bi << "/eigenvalue";
	  sPath << eigenstatesGroup << "/twist_" 
		    << tindex << "/band_" << bi << "/spin";
	}
	else if (NumBands > 1) {
	  ePath << eigenstatesGroup << "/twist/band_" << bi << "/eigenvalue";
	  sPath << eigenstatesGroup << "/twist/band_" << bi << "/spin";
	}
	else {
	  ePath << eigenstatesGroup << "/twist/band/eigenvalue";
	  sPath << eigenstatesGroup << "/twist/band/spin";
	}
	
	HDFAttribIO<double> h_energy(band.Energy);
	HDFAttribIO<int> h_spin(band.Spin);
	band.Energy = -1.01e100;
	h_energy.read (H5FileID, ePath.str().c_str());
	if (band.Energy > -1.0e100) {
	  h_spin.read   (H5FileID, sPath.str().c_str());
	  if (band.Spin == spin)
	    SortBands.push_back(band);
	}
      }
    }
    // Now sort the bands by energy
    if (sortBands) {
      cerr << "Sorting the bands now:\n";
      sort (SortBands.begin(), SortBands.end());
    }
  }


  void
  EinsplineSetBuilder::ReadBands (int spin, EinsplineSetLocal* orbitalSet)
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    // Read in the occupied bands
    orbitalSet->Orbitals.resize(orbitalSet->getOrbitalSetSize());
    cerr << "Orbitals size = " << orbitalSet->Orbitals.size() << endl;
    int iorb  = 0;
    int iband = 0;
    while (iorb < orbitalSet->getOrbitalSetSize()) {
      int ti   = SortBands[iband].TwistIndex;
      int bi   = SortBands[iband].BandIndex;
      double e = SortBands[iband].Energy;
      PosType twist, k;
      twist = TwistAngles[ti];
      Tensor<double,3> G = orbitalSet->PrimLattice.G;
      k = orbitalSet->PrimLattice.k_cart(twist);
      fprintf (stderr, "  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f)\n", 
	       ti, bi, e, k[0], k[1], k[2]);
      
      EinsplineOrb<complex<double>,OHMMS_DIM> *orb;
      // Check to see if we have already read this orbital, perhaps on
      // another processor in this OpenMP node.
      std::map<TinyVector<int,4>,OrbType*,Int4less>::iterator iter =
	OrbitalMap.find(TinyVector<int,4>(spin,ti,bi,0));      
      if (iter != OrbitalMap.end()) 
	orb = iter->second;
      else { // The orbital has not yet been read, so read it now
	orb = new EinsplineOrb<complex<double>,OHMMS_DIM>;
	OrbitalMap[TinyVector<int,4>(spin, ti, bi, 0)] = orb;
	orb->kVec = k;
	orb->Lattice = SuperLattice;
	ostringstream groupPath;
	
	if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
	  groupPath << eigenstatesGroup << "/twist_" 
		    << ti << "/band_" << bi << "/";
	else if (NumBands > 1)
	  groupPath << eigenstatesGroup << "/twist/band_" << bi << "/";
	else 
	  groupPath << eigenstatesGroup << "/twist/band/";
	
	orb->read(H5FileID, groupPath.str());
      }
      orbitalSet->Orbitals[iorb] = orb;
      iorb++;
      if (orb->uCenters.size() > 1) 
	cerr << "Making " << orb->uCenters.size() << " copies of band "
	     << iband << endl;
      // If the orbital has more than one center associated with it,
      // make copies of the orbital, changing only the center
      // associated with it.
      for (int icopy=1; icopy<orb->uCenters.size(); icopy++) {
	iter = OrbitalMap.find(TinyVector<int,4>(spin,ti,bi,icopy));
	EinsplineOrb<complex<double>,OHMMS_DIM> *orbCopy;
	if (iter != OrbitalMap.end()) 
	  orbCopy = iter->second;
	else {
	  orbCopy = new EinsplineOrb<complex<double>,OHMMS_DIM>(*orb);
	  OrbitalMap[TinyVector<int,4>(spin, ti, bi, icopy)] = orbCopy;
	  orbCopy->uCenter = orbCopy->uCenters[icopy];
	  if (orb->Reflections.size() > icopy)
	    orbCopy->Reflection = orbCopy->Reflections[icopy];
	  orbCopy->Center  = orb->Lattice.toCart(orbCopy->uCenter);
	}
	orbitalSet->Orbitals[iorb] = orbCopy;
	iorb++;
      }
      iband++;      
    }
  }
  

  void
  EinsplineSetBuilder::ReadBands 
  (int spin, EinsplineSetExtended<double>* orbitalSet)
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";
    
  }

  string
  EinsplineSetBuilder::OrbitalPath(int ti, int bi)
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    ostringstream groupPath;
    
    if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
      groupPath << eigenstatesGroup << "/twist_" 
		<< ti << "/band_" << bi << "/";
    else if (NumBands > 1)
      groupPath << eigenstatesGroup << "/twist/band_" << bi << "/";
    else 
      groupPath << eigenstatesGroup << "/twist/band/";

    return groupPath.str();
  }

  void
  EinsplineSetBuilder::ReadBands 
  (int spin, EinsplineSetExtended<complex<double> >* orbitalSet)
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    // Find the orbital mesh size
    int ti   = SortBands[0].TwistIndex;
    int bi   = SortBands[0].BandIndex;
    Array<complex<double>,3> rawData, splineData;
    string vectorName = OrbitalPath (ti, bi) + "eigenvector";
    HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
    h_rawData.read(H5FileID, vectorName.c_str());

    Ugrid x_grid, y_grid, z_grid;
    BCtype_z xBC, yBC, zBC;

    int nx, ny, nz;
    nx = rawData.size(0); 
    ny = rawData.size(1);
    nz = rawData.size(2);
    xBC.lCode = PERIODIC;    xBC.rCode = PERIODIC;
    yBC.lCode = PERIODIC;    yBC.rCode = PERIODIC;
    zBC.lCode = PERIODIC;    zBC.rCode = PERIODIC;
    x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = nx-1;
    y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = ny-1;
    z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = nz-1;

    // Create the multiUBspline object
    cerr << "OrbitalSetSize = " << orbitalSet->getOrbitalSetSize() << endl;
    orbitalSet->MultiSpline = create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC,
							  orbitalSet->getOrbitalSetSize());
    splineData.resize(nx-1, ny-1, nz-1);
    for (int ix=0; ix<(nx-1); ix++)
      for (int iy=0; iy<(ny-1); iy++)
	for (int iz=0; iz<(nz-1); iz++)
	  splineData(ix,iy,iz) = rawData(ix,iy,iz);

    set_multi_UBspline_3d_z (orbitalSet->MultiSpline, 0, splineData.data());
    orbitalSet->kPoints.resize(orbitalSet->getOrbitalSetSize());
        
    int iorb  = 1;    
    while (iorb < orbitalSet->getOrbitalSetSize()) {
      int ti   = SortBands[iorb].TwistIndex;
      int bi   = SortBands[iorb].BandIndex;
      double e = SortBands[iorb].Energy;
      PosType twist, k;
      twist = TwistAngles[ti];
      k = orbitalSet->PrimLattice.k_cart(twist);
      orbitalSet->kPoints[iorb] = k;
      fprintf (stderr, "  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f)\n", 
	       ti, bi, e, k[0], k[1], k[2]);
      
      vectorName = OrbitalPath (ti, bi) + "eigenvector";
      HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
      h_rawData.read(H5FileID, vectorName.c_str());
      if ((rawData.size(0) != nx) ||
	  (rawData.size(1) != ny) ||
	  (rawData.size(2) != nz)) {
	fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
	fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
	abort();
      }

      for (int ix=0; ix<(nx-1); ix++)
	for (int iy=0; iy<(ny-1); iy++)
	  for (int iz=0; iz<(nz-1); iz++)
	    splineData(ix,iy,iz) = rawData(ix,iy,iz);
      set_multi_UBspline_3d_z (orbitalSet->MultiSpline, iorb, splineData.data());
     
      iorb++;
    }
  }
  


  void
  EinsplineSetBuilder::CopyBands(int numOrbs) 
  {
    if (dynamic_cast<EinsplineSetLocal*>(OrbitalSet)) {
      EinsplineSetLocal *restrict locOrbitalSet = 
	dynamic_cast<EinsplineSetLocal*>(OrbitalSet);
      EinsplineSetLocal *restrict lastOrbitalSet =
	dynamic_cast<EinsplineSetLocal*>(LastOrbitalSet);
      locOrbitalSet->Orbitals.resize(numOrbs);
      for (int i=0; i<numOrbs; i++) 
	locOrbitalSet->Orbitals[i] = lastOrbitalSet->Orbitals[i];
    }
  }
}
