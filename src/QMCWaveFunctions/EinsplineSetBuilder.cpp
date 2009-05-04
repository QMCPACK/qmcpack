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
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <vector>
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "ParticleIO/ESHDFParticleParser.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus {
  std::map<TinyVector<int,4>,EinsplineSetBuilder::OrbType*,Int4less> 
  EinsplineSetBuilder::OrbitalMap;
  std::map<H5OrbSet,multi_UBspline_3d_z*,H5OrbSet>
  EinsplineSetBuilder::ExtendedMap_z;
  std::map<H5OrbSet,multi_UBspline_3d_d*,H5OrbSet>
  EinsplineSetBuilder::ExtendedMap_d;

  std::map<H5OrbSet,SPOSetBase*,H5OrbSet>  EinsplineSetBuilder::SPOSetMap;


  EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, 
      PtclPoolType& psets, xmlNodePtr cur) 
    : XMLRoot(cur), TileFactor(1,1,1), TwistNum(0), LastSpinSet(-1), 
      NumOrbitalsRead(-1), NumMuffinTins(0), NumCoreStates(0),
      ParticleSets(psets), TargetPtcl(p), H5FileID(-1),
      Format(QMCPACK)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	TileMatrix(i,j) = 0;
  }

  EinsplineSetBuilder::~EinsplineSetBuilder()
  {
    DEBUG_MEMORY("EinsplineSetBuilder::~EinsplineSetBuilder");
    if(H5FileID>=0) H5Fclose(H5FileID);
  }


  bool 
  EinsplineSetBuilder::put(xmlNodePtr cur) {
    string hdfName;
    OhmmsAttributeSet attribs;
    attribs.add (hdfName, "href");
    return attribs.put(XMLRoot);
  }

  bool
  EinsplineSetBuilder::ReadOrbitalInfo_ESHDF()
  {
    app_log() << "  Reading orbital file in ESHDF format.\n";
    TinyVector<int,3> version;
    HDFAttribIO<TinyVector<int,3> > h_version(version);
    h_version.read (H5FileID, "/version");
    app_log() << "  HDF5 orbital file version " 
	      << version[0] << "." << version[1] << "." << version[2] << endl;

    HDFAttribIO<Tensor<double,3> > h_Lattice(Lattice);
    h_Lattice.read      (H5FileID, "/supercell/primitive_vectors");
    RecipLattice = 2.0*M_PI*inverse(Lattice);
    SuperLattice = dot(TileMatrix, Lattice);

    char buff[1000];

    snprintf (buff, 1000, 
	     "  Lattice = \n    [ %9.6f %9.6f %9.6f\n"
	     "      %9.6f %9.6f %9.6f\n"
	     "      %9.6f %9.6f %9.6f ]\n", 
	     Lattice(0,0), Lattice(0,1), Lattice(0,2), 
	     Lattice(1,0), Lattice(1,1), Lattice(1,2), 
	     Lattice(2,0), Lattice(2,1), Lattice(2,2));
    app_log() << buff;
    snprintf (buff, 1000, 
              "  SuperLattice = \n    [ %9.6f %9.6f %9.6f\n"
                  "      %9.6f %9.6f %9.6f\n"
                  "      %9.6f %9.6f %9.6f ]\n", 
	      SuperLattice(0,0), SuperLattice(0,1), SuperLattice(0,2), 
	      SuperLattice(1,0), SuperLattice(1,1), SuperLattice(1,2), 
	      SuperLattice(2,0), SuperLattice(2,1), SuperLattice(2,2));
    app_log() << buff;
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	LatticeInv(i,j) = RecipLattice(i,j)/(2.0*M_PI);

    int have_dpsi = 0;
    int NumAtomicOrbitals = 0;
    HDFAttribIO<int> h_NumBands(NumBands), h_NumElectrons(NumElectrons), 
      h_NumSpins(NumSpins), h_NumTwists(NumTwists), h_NumCore(NumCoreStates),
      h_NumMuffinTins(NumMuffinTins), h_have_dpsi(have_dpsi), 
      h_NumAtomicOrbitals(NumAtomicOrbitals);
    NumCoreStates = NumMuffinTins;

    h_NumBands.read      (H5FileID, "/electrons/kpoint_0/spin_0/number_of_states");
    h_NumCore.read       (H5FileID, "/electrons/kpoint_0/spin_0/number_of_core_states");
    h_NumElectrons.read  (H5FileID, "/electrons/number_of_electrons");
    h_NumSpins.read      (H5FileID, "/electrons/number_of_spins");
    h_NumTwists.read     (H5FileID, "/electrons/number_of_kpoints");
    h_NumMuffinTins.read (H5FileID, "/muffin_tins/number_of_tins");
    h_have_dpsi.read     (H5FileID, "/electrons/have_dpsi");
    h_NumAtomicOrbitals.read (H5FileID, "/electrons/number_of_atomic_orbitals");
    HaveOrbDerivs = have_dpsi;
    app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons 
	      << ", spins=" << NumSpins << ", twists=" << NumTwists 
	      << ", muffin tins=" << NumMuffinTins << endl;
    app_log() << "atomic orbital=" << NumAtomicOrbitals << endl;
    if (TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1)
      app_log() << "  Using a " << TileFactor[0] << "x" << TileFactor[1] 
		<< "x" << TileFactor[2] << " tiling factor.\n";

    //////////////////////////////////
    // Read ion types and locations //
    //////////////////////////////////
    Vector<int> species_ids;
    HDFAttribIO<Vector<int> > h_species_ids(species_ids);
    h_species_ids.read (H5FileID, "/atoms/species_ids");
    int num_species;
    HDFAttribIO<int> h_num_species(num_species);
    h_num_species.read (H5FileID, "/atoms/number_of_species");
    vector<int> atomic_numbers(num_species);
    for (int isp=0; isp<num_species; isp++) {
      ostringstream name;
      name << "/atoms/species_" << isp << "/atomic_number";
      HDFAttribIO<int> h_atomic_number (atomic_numbers[isp]);
      h_atomic_number.read(H5FileID, name.str().c_str());
    }
    IonTypes.resize(species_ids.size());
    for (int i=0; i<species_ids.size(); i++) 
      IonTypes[i] = atomic_numbers[species_ids[i]];

    HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
    h_IonPos.read   (H5FileID, "/atoms/positions");
    for (int i=0; i<IonTypes.size(); i++)
      app_log() << "Atom type(" << i << ") = " << IonTypes(i) << endl;

    /////////////////////////////////////
    // Read atomic orbital information //
    /////////////////////////////////////
    AtomicOrbitals.resize(NumAtomicOrbitals);
    for (int iat=0; iat<NumAtomicOrbitals; iat++) {
      AtomicOrbital<complex<double> > &orb = AtomicOrbitals[iat];
      int lmax, polynomial_order, spline_points;
      RealType cutoff_radius, polynomial_radius, spline_radius;
      PosType position;
      HDFAttribIO<int> h_lmax(lmax), h_polynomial_order(polynomial_order),
	h_spline_points(spline_points);
      HDFAttribIO<RealType> h_cutoff_radius(cutoff_radius), 
	h_polynomial_radius(polynomial_radius),
	h_spline_radius(spline_radius);
      HDFAttribIO<PosType> h_position(position);
      ostringstream groupstream;
      groupstream << "/electrons/atomic_orbital_" << iat << "/";
      string groupname = groupstream.str();
      h_lmax.read              (H5FileID, (groupname + "lmax"             ).c_str());
      h_polynomial_order.read  (H5FileID, (groupname + "polynomial_order" ).c_str());
      h_spline_points.read     (H5FileID, (groupname + "spline_points"    ).c_str());
      h_cutoff_radius.read     (H5FileID, (groupname + "cutoff_radius"    ).c_str());
      h_polynomial_radius.read (H5FileID, (groupname + "polynomial_radius").c_str());
      h_spline_radius.read     (H5FileID, (groupname + "spline_radius"    ).c_str());
      h_position.read          (H5FileID, (groupname + "position"         ).c_str());
      orb.set_pos (position);
      orb.set_lmax (lmax);
      orb.set_cutoff (cutoff_radius);
      orb.set_spline (spline_radius, spline_points);
      orb.set_polynomial (polynomial_radius, polynomial_order);
    }
      

    ///////////////////////////
    // Read the twist angles //
    ///////////////////////////
    TwistAngles.resize(NumTwists);
    for (int ti=0; ti<NumTwists; ti++) {
      ostringstream path;
      path << "/electrons/kpoint_" << ti << "/reduced_k";
      HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
      h_Twist.read (H5FileID, path.str().c_str());
      snprintf (buff, 1000, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n", 
	       TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
      app_log() << buff;
    }

    //////////////////////////////////////////////////////////
    // If the density has not been set in TargetPtcl, and   //
    // the density is available, read it in and save it     //
    // in TargetPtcl.                                       //
    //////////////////////////////////////////////////////////

    // FIXME:  add support for more than one spin density
    if (!TargetPtcl.Density_G.size()) {
      HDFAttribIO<vector<TinyVector<int,OHMMS_DIM> > > 
	h_reduced_gvecs(TargetPtcl.DensityReducedGvecs);
      HDFAttribIO<Array<RealType,OHMMS_DIM> > 
	h_density_r (TargetPtcl.Density_r);
      TinyVector<int,3> mesh;
      h_reduced_gvecs.read (H5FileID, "/electrons/density/gvectors");

      int numG = TargetPtcl.DensityReducedGvecs.size();
      // Convert primitive G-vectors to supercell G-vectors
      // Also, flip sign since ESHDF format uses opposite sign convention
      for (int iG=0; iG < numG; iG++) 
	TargetPtcl.DensityReducedGvecs[iG] = 
	  -1 * dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
      app_log() << "  Read " << numG << " density G-vectors.\n";

      for (int ispin=0; ispin<NumSpins; ispin++) {
	ostringstream density_r_path, density_g_path;
	density_r_path << "/electrons/density/spin_" << ispin << "/density_r";
	density_g_path << "/electrons/density/spin_" << ispin << "/density_g";
	h_density_r.read (H5FileID, density_r_path.str().c_str());

	if (TargetPtcl.DensityReducedGvecs.size()) {
	  app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
	  vector<ComplexType> density_G;
	  HDFAttribIO<vector<ComplexType > > h_density_G (density_G);
	  h_density_G.read (H5FileID, density_g_path.str().c_str());
	  if (!density_G.size()) {
	    app_error() << "  Density reduced G-vectors defined, but not the"
			<< " density.\n";
	    abort();
	  }
	  else {
	    if (ispin == 0)
	      TargetPtcl.Density_G = density_G;
	    else
	      for (int iG=0; iG<density_G.size(); iG++)
		TargetPtcl.Density_G[iG] += density_G[iG];
	  }
	}
      }
    }
    HaveLocalizedOrbs = false;
    return true;


  }







  bool
  EinsplineSetBuilder::ReadOrbitalInfo()
  {
    //H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    if (H5FileID < 0) {
      app_error() << "Could not open HDF5 file \"" << H5FileName 
		  << "\" in EinsplineSetBuilder::createSPOSet.  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::ReadOrbitalInfo");
    }
    
    // Read format
    std::string format;
    HDFAttribIO<string> h_format(format);
    h_format.read(H5FileID, "/format");
    if (format.find("ES")<format.size())
    {
      Format = ESHDF;
      return ReadOrbitalInfo_ESHDF();
    }


    //////////////////////////////////////////////////
    // Read basic parameters from the orbital file. //
    //////////////////////////////////////////////////
    // Check the version
    HDFAttribIO<TinyVector<int,2> > h_Version(Version);
    h_Version.read (H5FileID, "/version");
    app_log() << "  HDF5 orbital file version " 
	      << Version[0] << "." << Version[1] << "\n";
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
      ostringstream o;
      o << "Unknown HDF5 orbital file version " << Version[0] << "." << Version[1] << "\n";
      APP_ABORT(o.str());
      //abort();
    }
    HDFAttribIO<Tensor<double,3> > h_Lattice(Lattice), h_RecipLattice(RecipLattice);
    h_Lattice.read      (H5FileID, (parameterGroup+"/lattice").c_str());
    
    h_RecipLattice.read (H5FileID, (parameterGroup+"/reciprocal_lattice").c_str());
    SuperLattice = dot(TileMatrix, Lattice);

    char buff[1000];

    snprintf (buff, 1000, 
	     "  Lattice = \n    [ %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f\n"
	     "      %8.5f %8.5f %8.5f ]\n", 
	     Lattice(0,0), Lattice(0,1), Lattice(0,2), 
	     Lattice(1,0), Lattice(1,1), Lattice(1,2), 
	     Lattice(2,0), Lattice(2,1), Lattice(2,2));
    app_log() << buff;
    snprintf (buff, 1000, 
              "  SuperLattice = \n    [ %13.12f %13.12f %13.12f\n"
                  "      %13.12f %13.12f %13.12f\n"
                  "      %13.12f %13.12f %13.12f ]\n", 
	      SuperLattice(0,0), SuperLattice(0,1), SuperLattice(0,2), 
	      SuperLattice(1,0), SuperLattice(1,1), SuperLattice(1,2), 
	      SuperLattice(2,0), SuperLattice(2,1), SuperLattice(2,2));
    app_log() << buff;
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	LatticeInv(i,j) = RecipLattice(i,j)/(2.0*M_PI);
    HDFAttribIO<int> h_NumBands(NumBands), h_NumElectrons(NumElectrons), 
      h_NumSpins(NumSpins), h_NumTwists(NumTwists), h_NumCore(NumCoreStates),
      h_NumMuffinTins(NumMuffinTins);
    NumCoreStates = NumMuffinTins = 0;
    h_NumBands.read      (H5FileID, (parameterGroup+"/num_bands").c_str());
    h_NumCore.read       (H5FileID, (parameterGroup+"/num_core_states").c_str());
    h_NumElectrons.read  (H5FileID, (parameterGroup+"/num_electrons").c_str());
    h_NumSpins.read      (H5FileID, (parameterGroup+"/num_spins").c_str());
    h_NumTwists.read     (H5FileID, (parameterGroup+"/num_twists").c_str());
    h_NumMuffinTins.read (H5FileID, (parameterGroup+"/muffin_tins/num_tins").c_str());
    app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons 
	      << ", spins=" << NumSpins << ", twists=" << NumTwists 
	      << ", muffin tins=" << NumMuffinTins << endl;
    // fprintf (stderr, "  bands = %d, elecs = %d, spins = %d, twists = %d\n",
    // 	     NumBands, NumElectrons, NumSpins, NumTwists);
    if (TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1)
      app_log() << "  Using a " << TileFactor[0] << "x" << TileFactor[1] 
		<< "x" << TileFactor[2] << " tiling factor.\n";

    /////////////////////////////////
    // Read muffin tin information //
    /////////////////////////////////
    MT_APW_radii.resize(NumMuffinTins);
    MT_APW_rgrids.resize(NumMuffinTins);
    MT_APW_lmax.resize(NumMuffinTins);
    MT_APW_num_radial_points.resize(NumMuffinTins);
    MT_centers.resize(NumMuffinTins);
    for (int tin=0; tin<NumMuffinTins; tin++) {
      ostringstream MTstream;
      if (NumMuffinTins > 1)
	MTstream << parameterGroup << "/muffin_tins/muffin_tin_" << tin;
      else
	MTstream << parameterGroup << "/muffin_tins/muffin_tin";
      string MTgroup = MTstream.str();
      HDFAttribIO<int> h_lmax(MT_APW_lmax[tin]), 
	h_num_radial_points(MT_APW_num_radial_points[tin]);
      HDFAttribIO<RealType> h_radius (MT_APW_radii[tin]);
      HDFAttribIO<PosType> h_center (MT_centers[tin]);
      HDFAttribIO<Vector<double> > h_rgrid (MT_APW_rgrids[tin]);
      h_lmax.read              (H5FileID, (MTgroup+"/lmax").c_str());
      h_num_radial_points.read (H5FileID, (MTgroup+"/num_radial_points").c_str());
      h_radius.read            (H5FileID, (MTgroup+"/radius").c_str());
      h_center.read            (H5FileID, (MTgroup+"/center").c_str());
      h_rgrid.read             (H5FileID, (MTgroup+"/r"     ).c_str());
    }

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
      snprintf (buff, 1000, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n", 
	       TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
      app_log() << buff;
    }

    /////////////////////////////////////////////////////////
    // Determine whether we need to use localized orbitals //
    /////////////////////////////////////////////////////////
    HaveLocalizedOrbs = false;
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

    //////////////////////////////////////////////////////////
    // If the density has not been set in TargetPtcl, and   //
    // the density is available, read it in and save it     //
    // in TargetPtcl.                                       //
    //////////////////////////////////////////////////////////
    if (!TargetPtcl.Density_G.size()) {
      HDFAttribIO<vector<TinyVector<int,OHMMS_DIM> > > 
	h_reduced_gvecs(TargetPtcl.DensityReducedGvecs);
      HDFAttribIO<Array<RealType,OHMMS_DIM> > 
	h_density_r (TargetPtcl.Density_r);
      h_reduced_gvecs.read (H5FileID, "/density/reduced_gvecs");
      h_density_r.read (H5FileID,     "/density/rho_r");

      int numG = TargetPtcl.DensityReducedGvecs.size();
      // Convert primitive G-vectors to supercell G-vectors
      for (int iG=0; iG < numG; iG++) 
	TargetPtcl.DensityReducedGvecs[iG] = 
	  dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);

      app_log() << "  Read " << numG << " density G-vectors.\n";
      if (TargetPtcl.DensityReducedGvecs.size()) {
	app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
	HDFAttribIO<vector<ComplexType > > h_density_G (TargetPtcl.Density_G);
	h_density_G.read (H5FileID, "/density/rho_G");
	if (!TargetPtcl.Density_G.size()) {
	  app_error() << "  Density reduced G-vectors defined, but not the"
		      << " density.\n";
	  abort();
	}
      }
    }

    return true;
  }

  void
  EinsplineSetBuilder::BroadcastOrbitalInfo()
  {
    if(myComm->size() == 1) return;

    int numIons = IonTypes.size();
    int numAtomicOrbitals = AtomicOrbitals.size();
    int numDensityGvecs = TargetPtcl.DensityReducedGvecs.size();
    PooledData<RealType> abuffer;
    PooledData<int>       aibuffer;
    aibuffer.add(Version.begin(),Version.end()); //myComm->bcast(Version);
    aibuffer.add(Format);
    abuffer.add(Lattice.begin(),Lattice.end());//myComm->bcast(Lattice);
    abuffer.add(RecipLattice.begin(),RecipLattice.end()); //myComm->bcast(RecipLattice);
    abuffer.add(SuperLattice.begin(),SuperLattice.end()); //myComm->bcast(SuperLattice);
    abuffer.add(LatticeInv.begin(),LatticeInv.end()); //myComm->bcast(LatticeInv);
    aibuffer.add(NumBands); //myComm->bcast(NumBands);
    aibuffer.add(NumElectrons); //myComm->bcast(NumElectrons);
    aibuffer.add(NumSpins); //myComm->bcast(NumSpins);
    aibuffer.add(NumTwists); //myComm->bcast(NumTwists);
    aibuffer.add(numIons); //myComm->bcast(numIons);
    aibuffer.add(NumMuffinTins);
    aibuffer.add(numAtomicOrbitals);
    aibuffer.add(numDensityGvecs);
    aibuffer.add((int)HaveOrbDerivs);

    myComm->bcast(abuffer);
    myComm->bcast(aibuffer);

    if(myComm->rank())
    {
      abuffer.rewind();
      aibuffer.rewind();
      aibuffer.get(Version.begin(),Version.end());
      aibuffer.get(Format);
      abuffer.get(Lattice.begin(),Lattice.end());
      abuffer.get(RecipLattice.begin(),RecipLattice.end());
      abuffer.get(SuperLattice.begin(),SuperLattice.end());
      abuffer.get(LatticeInv.begin(),LatticeInv.end());
      aibuffer.get(NumBands);
      aibuffer.get(NumElectrons);
      aibuffer.get(NumSpins);
      aibuffer.get(NumTwists);
      aibuffer.get(numIons);
      aibuffer.get(NumMuffinTins);
      aibuffer.get(numAtomicOrbitals);
      aibuffer.get(numDensityGvecs);
      aibuffer.get(HaveOrbDerivs);
      MT_APW_radii.resize(NumMuffinTins);
      MT_APW_lmax.resize(NumMuffinTins);
      MT_APW_rgrids.resize(NumMuffinTins);
      MT_APW_num_radial_points.resize(NumMuffinTins);
      MT_centers.resize(NumMuffinTins);
      TargetPtcl.DensityReducedGvecs.resize(numDensityGvecs);
      TargetPtcl.Density_G.resize(numDensityGvecs);
      AtomicOrbitals.resize(numAtomicOrbitals);
    }

    vector<int> rgrids_sizes(NumMuffinTins);
    for (int tin=0; tin<NumMuffinTins; tin++) 
      rgrids_sizes[tin] = MT_APW_rgrids[tin].size();
    
    myComm->bcast(rgrids_sizes);
    if (myComm->rank())
      for (int tin=0; tin<NumMuffinTins; tin++)
	MT_APW_rgrids[tin].resize(rgrids_sizes[tin]);
    

    if (IonTypes.size() != numIons) {
      IonTypes.resize(numIons);
      IonPos.resize(numIons);
    }

    //new buffer
    PooledData<RealType> bbuffer;
    PooledData<int> bibuffer;
    for(int i=0; i<numIons; ++i) bibuffer.add(IonTypes[i]);
    //myComm->bcast(IonTypes);
    
    bbuffer.add(&IonPos[0][0],&IonPos[0][0]+OHMMS_DIM*numIons);
    //myComm->bcast(IonPos);

    if (TwistAngles.size() != NumTwists) TwistAngles.resize(NumTwists);
    bbuffer.add(&TwistAngles[0][0],&TwistAngles[0][0]+OHMMS_DIM*NumTwists);
    //myComm->bcast(TwistAngles);
    
    bibuffer.add(HaveLocalizedOrbs);
    //myComm->bcast(HaveLocalizedOrbs);

    bbuffer.add(MT_APW_radii.begin(), MT_APW_radii.end());
    bibuffer.add(MT_APW_lmax.begin(),  MT_APW_lmax.end());
    bibuffer.add(MT_APW_num_radial_points.begin(), 
		 MT_APW_num_radial_points.end());
      

    bbuffer.add(&(MT_centers[0][0]), &(MT_centers[0][0])+OHMMS_DIM*NumMuffinTins);
    for (int i=0; i<NumMuffinTins; i++) 
      bbuffer.add(MT_APW_rgrids[i].begin(), MT_APW_rgrids[i].end());
    bibuffer.add(&(TargetPtcl.DensityReducedGvecs[0][0]),
		&(TargetPtcl.DensityReducedGvecs[0][0])+numDensityGvecs*OHMMS_DIM);
    bbuffer.add(&(TargetPtcl.Density_G[0]),
		&(TargetPtcl.Density_G[0]) + numDensityGvecs);
    for (int iat=0; iat<numAtomicOrbitals; iat++) {
      AtomicOrbital<complex<double> > &orb = AtomicOrbitals[iat];
      bibuffer.add (orb.SplinePoints);
      bibuffer.add (orb.PolyOrder);
      bibuffer.add (orb.lMax);
      bibuffer.add (orb.Numlm);
      bbuffer.add  (&orb.Pos[0], &orb.Pos[0]+OHMMS_DIM);
      bbuffer.add  (orb.CutoffRadius);
      bbuffer.add  (orb.SplineRadius);
      bbuffer.add  (orb.PolyRadius);
    }

    myComm->bcast(bbuffer);
    myComm->bcast(bibuffer);

    if(myComm->rank())
    {
      bbuffer.rewind();
      bibuffer.rewind();
      for(int i=0; i<numIons; ++i) bibuffer.get(IonTypes[i]);
      bbuffer.get(&IonPos[0][0],&IonPos[0][0]+OHMMS_DIM*numIons);
      bbuffer.get(&TwistAngles[0][0],&TwistAngles[0][0]+OHMMS_DIM*NumTwists);
      bibuffer.get(HaveLocalizedOrbs);
      bbuffer.get(MT_APW_radii.begin(), MT_APW_radii.end());
      bibuffer.get(MT_APW_lmax.begin(),  MT_APW_lmax.end());
      bibuffer.get(MT_APW_num_radial_points.begin(), 
		   MT_APW_num_radial_points.end());
      bbuffer.get(&(MT_centers[0][0]), 
		  &(MT_centers[0][0])+OHMMS_DIM*NumMuffinTins);
      for (int i=0; i<NumMuffinTins; i++) 
	bbuffer.get(MT_APW_rgrids[i].begin(), MT_APW_rgrids[i].end());
      bibuffer.get(&(TargetPtcl.DensityReducedGvecs[0][0]),
		  &(TargetPtcl.DensityReducedGvecs[0][0])+
		   numDensityGvecs*OHMMS_DIM);
      bbuffer.get(&(TargetPtcl.Density_G[0]),
		  &(TargetPtcl.Density_G[0]) + numDensityGvecs);
      for (int iat=0; iat<numAtomicOrbitals; iat++) {
	AtomicOrbital<complex<double> > &orb = AtomicOrbitals[iat];
	bibuffer.get (orb.SplinePoints);
	bibuffer.get (orb.PolyOrder);
	bibuffer.get (orb.lMax);
	bibuffer.get (orb.Numlm);
	bbuffer.get  (&orb.Pos[0], &orb.Pos[0]+OHMMS_DIM);
	bbuffer.get  (orb.CutoffRadius);
	bbuffer.get  (orb.SplineRadius);
	bbuffer.get  (orb.PolyRadius);
      }

    }
  }

  SPOSetBase*
  EinsplineSetBuilder::createSPOSet(xmlNodePtr cur) {
    OhmmsAttributeSet attribs;
    int numOrbs = 0;
    bool sortBands = true;
    string sourceName;
    attribs.add (H5FileName, "href");
    attribs.add (TileFactor, "tile");
    attribs.add (sortBands,  "sort");
    attribs.add (TileMatrix, "tilematrix");
    attribs.add (TwistNum,   "twistnum");
    attribs.add (sourceName, "source");
    attribs.put (XMLRoot);
    attribs.add (numOrbs,    "size");
    attribs.add (numOrbs,    "norbs");
    attribs.put (cur);

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
          APP_ABORT("EinsplineSetBuilder::createSPOSet");
	}
      }
      cur = cur->next;
    }


    H5OrbSet set(H5FileName, spinSet, numOrbs);
    std::map<H5OrbSet,SPOSetBase*,H5OrbSet>::iterator iter;
    iter = SPOSetMap.find (set);
    if (iter != SPOSetMap.end()) {
      app_log() << "SPOSet parameters match in EinsplineSetBuilder:  "
		<< "cloning EinsplineSet object.\n";
      return iter->second->makeClone();
    }


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

    if (myComm->rank() == 0) 
      fprintf (stderr, " [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d ]\n",
	       TileMatrix(0,0), TileMatrix(0,1), TileMatrix(0,2),
	       TileMatrix(1,0), TileMatrix(1,1), TileMatrix(1,2),
	       TileMatrix(2,0), TileMatrix(2,1), TileMatrix(2,2));
    
    if (numOrbs == 0) {
      app_error() << "You must specify the number of orbitals in the input file.\n";
      APP_ABORT("EinsplineSetBuilder::createSPOSet");
    }
    else 
      app_log() << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";


    /////////////////////////////////////////////////////////////////
    // Read the basic orbital information, without reading all the //
    // orbitals themselves.                                        //
    /////////////////////////////////////////////////////////////////
    if (myComm->rank() == 0) 
      if (!ReadOrbitalInfo()) {
	app_error() << "Error reading orbital info from HDF5 file.  Aborting.\n";
        APP_ABORT("EinsplineSetBuilder::createSPOSet");
      }
    
    BroadcastOrbitalInfo();

    ///////////////////////////////////////////////////////////////////
    // Now, analyze the k-point mesh to figure out the what k-points //
    // are needed                                                    //
    ///////////////////////////////////////////////////////////////////
    PrimCell.set(Lattice);
    SuperCell.set(SuperLattice);
    for (int iat=0; iat<AtomicOrbitals.size(); iat++)
      AtomicOrbitals[iat].Lattice = Lattice;

    // Copy supercell into the ParticleSets
//     app_log() << "Overwriting XML lattice with that from the ESHDF file.\n";
//     PtclPoolType::iterator piter;
//     for(piter = ParticleSets.begin(); piter != ParticleSets.end(); piter++)
//       piter->second->Lattice.copy(SuperCell);

    AnalyzeTwists2();

    //////////////////////////////////
    // Create the OrbitalSet object //
    //////////////////////////////////
    if (HaveLocalizedOrbs) 
      OrbitalSet = new EinsplineSetLocal;
    else if (UseRealOrbitals)
      OrbitalSet = new EinsplineSetExtended<double>;
    else
      OrbitalSet = new EinsplineSetExtended<complex<double> >;


    ///JEREMY    OrbitalSet->resetSourceParticleSet(*ParticleSets["i"]);
    /////////////////////////
    // Setup internal data //
    /////////////////////////
    
    // Lattice information
    OrbitalSet->TileFactor = TileFactor;
    OrbitalSet->Tiling = (TileFactor[0]*TileFactor[1]*TileFactor[2] != 1);
    //OrbitalSet->Tiling = 
    //  TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1;

    OrbitalSet->PrimLattice  = Lattice;
    OrbitalSet->SuperLattice = SuperLattice;
    OrbitalSet->GGt=dot(OrbitalSet->PrimLattice.G,
			transpose(OrbitalSet->PrimLattice.G));
    OrbitalSet->setOrbitalSetSize (numOrbs);
    OrbitalSet->BasisSetSize   = numOrbs;
    TileIons();
    
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
      if (UseRealOrbitals) { 
	EinsplineSetExtended<double> *restrict orbitalSet =
	  dynamic_cast<EinsplineSetExtended<double>* > (OrbitalSet);    
        OccupyBands(spinSet, sortBands);
#pragma omp critical(read_extended_orbs)
	{ 
	  if (Format == ESHDF)
	    ReadBands_ESHDF(spinSet,orbitalSet);
	  else
	    ReadBands(spinSet, orbitalSet); 
	}
      }
      else {
	EinsplineSetExtended<complex<double> > *restrict orbitalSet = 
	  dynamic_cast<EinsplineSetExtended<complex<double> >*>(OrbitalSet);
	OccupyBands(spinSet, sortBands);
#pragma omp critical(read_extended_orbs)
	{ 
	  if (Format == ESHDF)
	    ReadBands_ESHDF(spinSet,orbitalSet);
	  else
	    ReadBands(spinSet, orbitalSet); 
	}
      }
    }

#ifndef QMC_COMPLEX
    if (myComm->rank()==0 && OrbitalSet->MuffinTins.size() > 0) {
      FILE *fout  = fopen ("TestMuffins.dat", "w");
      Vector<double> phi(numOrbs), lapl(numOrbs);
      Vector<PosType> grad(numOrbs);
      ParticleSet P;
      P.R.resize(6);
      for (int i=0; i<P.R.size(); i++)
	P.R[i] = PosType (0.0, 0.0, 0.0);
      PosType N = 0.25*PrimCell.a(0) + 0.25*PrimCell.a(1) + 0.25*PrimCell.a(2);
      for (double x=-1.0; x<=1.0; x+=0.0000500113412) {
	// for (double x=-0.003; x<=0.003; x+=0.0000011329343481381) {
	P.R[0] = x * (PrimCell.a(0) + 0.914*PrimCell.a(1) + 
		      0.781413*PrimCell.a(2));
	double r = std::sqrt(dot(P.R[0], P.R[0]));
	double rN = std::sqrt(dot(P.R[0]-N, P.R[0]-N));
	OrbitalSet->evaluate(P, 0, phi, grad, lapl);
	// OrbitalSet->evaluate(P, 0, phi);
	fprintf (fout, "%1.12e ", r*x/std::fabs(x));
	for (int j=0; j<numOrbs; j++) {
	  double gmag = std::sqrt(dot(grad[j],grad[j]));
	  fprintf (fout, "%16.12e ", 
		   /*phi[j]*phi[j]**/(-5.0/r  -0.5*lapl[j]/phi[j]));
	  // double E = -5.0/r -0.5*lapl[j]/phi[j];
	  fprintf (fout, "%16.12e ", phi[j]);
	  fprintf (fout, "%16.12e ", gmag);
	}
	fprintf (fout, "\n");
      }
      fclose(fout);
    }
#endif

    SPOSetMap[set] = OrbitalSet;
    
    if (sourceName.size() && (ParticleSets.find(sourceName) == ParticleSets.end()))
    {
      app_log() << "  EinsplineSetBuilder creates a ParticleSet " << sourceName << endl;
      ParticleSet* ions=new ParticleSet;
      ions->Lattice=TargetPtcl.Lattice;
      ESHDFIonsParser ap(*ions,H5FileID,myComm);
      ap.put(XMLRoot);
      ap.expand(TileMatrix);
      ions->setName(sourceName);
      ParticleSets[sourceName]=ions;

      //overwrite the lattice and assign random
      if(TargetPtcl.Lattice.SuperCellEnum)
      { 
        TargetPtcl.Lattice=ions->Lattice;
        makeUniformRandom(TargetPtcl.R);
        TargetPtcl.R.setUnit(PosUnit::LatticeUnit);
        TargetPtcl.convert2Cart(TargetPtcl.R);
	TargetPtcl.createSK();
      }
    }

    return OrbitalSet;
  }
  
  //////////////////////////////////////////////////////////////
  // Create the ion ParticleSet from the data in the HDF file //
  //////////////////////////////////////////////////////////////
  void
  EinsplineSetBuilder::CreateIonParticleSet(string sourceName)
  {
    //    ParticleSet &pTemp = *(new MCWalkerConfiguration);
    ParticleSet &pTemp = *(new ParticleSet);
    pTemp.setName (sourceName);
    SpeciesSet& tspecies(pTemp.getSpeciesSet());
    

    ParticleSets[sourceName] = &pTemp;
  }

  ////////////////////////////////////////////////////
  // Tile the ion positions according to TileMatrix //
  ////////////////////////////////////////////////////
  void 
  EinsplineSetBuilder::TileIons()
  {
    Vector<TinyVector<double, OHMMS_DIM> > primPos   = IonPos;
    Vector<int>                            primTypes = IonTypes;
    int numCopies = abs(TileMatrix.det());
    IonTypes.resize(primPos.size()*numCopies);
    IonPos.resize  (primPos.size()*numCopies);

    int maxCopies = 10;

    typedef TinyVector<double,3> Vec3;

    int index=0;

    for (int i0=-maxCopies; i0<=maxCopies; i0++)    
      for (int i1=-maxCopies; i1<=maxCopies; i1++)
	for (int i2=-maxCopies; i2<=maxCopies; i2++) 
	  for (int iat=0; iat < primPos.size(); iat++) {
	    Vec3 r     = primPos[iat];
	    Vec3 uPrim = PrimCell.toUnit(r);
	    for (int i=0; i<3; i++)   uPrim[i] -= std::floor(uPrim[i]);
	    r = PrimCell.toCart(uPrim) + (double)i0*PrimCell.a(0) + 
	      (double)i1*PrimCell.a(1) + (double)i2*PrimCell.a(2);
	    Vec3 uSuper = SuperCell.toUnit(r);
	    if ((uSuper[0] >= -1.0e-6) && (uSuper[0] < 0.9999) &&
		(uSuper[1] >= -1.0e-6) && (uSuper[1] < 0.9999) &&
		(uSuper[2] >= -1.0e-6) && (uSuper[2] < 0.9999)) {
	      IonPos[index]= r;
	      IonTypes[index]= primTypes[iat];
	      index++;
	    }
	  }
    if (index != primPos.size()*numCopies) {
      app_error() << "The number of tiled ions, " << IonPos.size() 
		  << ", does not match the expected number of "
		  << primPos.size()*numCopies << ".  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::TileIons()");
    }
    if (myComm->rank() == 0) {
      fprintf (stderr, "Supercell ion positions = \n");
      for (int i=0; i<IonPos.size(); i++)
	fprintf (stderr, "   [%14.12f %14.12f %14.12f]\n",
		 IonPos[i][0], IonPos[i][1], IonPos[i][2]);
    }
  }


  //inline TinyVector<double,3>
  //IntPart (TinyVector<double,3> twist)
  //{
  //  return TinyVector<double,3> (round(twist[0]-1.0e-10), round(twist[1]-1.0e-10), round(twist[2]-1.0e-10));
  //}
  //
  //inline TinyVector<double,3>
  //FracPart (TinyVector<double,3> twist)
  //{
  //  return twist - IntPart (twist);
  //}
  
  template<typename T>
  inline TinyVector<T,3>
  IntPart (const TinyVector<T,3>& twist)
  {
    return TinyVector<T,3> (round(twist[0]-1.0e-6),
                            round(twist[1]-1.0e-6), 
                            round(twist[2]-1.0e-6));
  }
  
  template<typename T>
  inline TinyVector<T,3>
  FracPart (const TinyVector<T,3>& twist)
  {
    return twist - IntPart (twist);
  }
  
  bool EinsplineSetBuilder::TwistPair (PosType a, PosType b)
  {
    bool pair = true;
    for (int n=0; n<OHMMS_DIM; n++) {
      double d = a[n] + b[n];
      if (std::fabs(d - round(d)) > 1.0e-8)
	pair = false;
    }
    return pair;
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
    app_log() << "Found " << numSuperTwists << " distinct supercell twists.\n";

    // For each supercell twist, create a list of primitive twists which
    // belong to it.
    vector<vector<int> > superSets;
    superSets.resize(numSuperTwists);
    for (int ki=0; ki<numPrimTwists; ki++)
      superSets[superIndex[ki]].push_back(ki);

    if (myComm->rank() == 0) 
      for (int si=0; si<numSuperTwists; si++) {
	fprintf (stderr, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n",
		 si, superFracs[si][0], superFracs[si][1], superFracs[si][2]);
	fprintf (stderr, "  Using k-points: ");
	for (int i=0; i<superSets[si].size(); i++) 
	  fprintf (stderr, " %d", superSets[si][i]);
	fprintf (stderr, "\n");
      }

    // Check supertwist for this node
    if (!myComm->rank()) 
      fprintf (stderr, "  Using supercell twist %d:  [ %9.5f %9.5f %9.5f]\n",
	       TwistNum, superFracs[TwistNum][0], superFracs[TwistNum][1],
	       superFracs[TwistNum][2]);
    
#ifndef QMC_COMPLEX
    // Check to see if supercell twist is okay to use with real wave
    // functions 
    for (int dim=0; dim<OHMMS_DIM; dim++) {
      double t = 2.0*superFracs[TwistNum][dim];
      if (std::fabs(t - round(t)) > 1.0e-10) {
	app_error() << "Cannot use this super twist with real wavefunctions.\n" 
		    << "Please recompile with QMC_COMPLEX=1.\n";
	abort();
      }
    }
#endif

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
	    app_error() << "Identical k-points detected in super twist set "
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

    // Now, find out which twists are distinct
    DistinctTwists.clear();
#ifndef QMC_COMPLEX
    vector<int> copyTwists;
    for (int i=0; i<IncludeTwists.size(); i++) {
      int ti        = IncludeTwists[i];
      PosType twist_i = TwistAngles[ti];
      bool distinct=true;
      for (int j=i+1; j<IncludeTwists.size(); j++) {
	int tj = IncludeTwists[j];
	PosType twist_j = TwistAngles[tj];
	PosType sum  = twist_i + twist_j;
	PosType diff = twist_i - twist_j;
	if (TwistPair (twist_i, twist_j)) 
	  distinct = false;
      }
      if (distinct)
	DistinctTwists.push_back(ti);
      else
	copyTwists.push_back(ti);
    }
    // Now determine which distinct twists require two copies
    MakeTwoCopies.resize(DistinctTwists.size());
    for (int i=0; i<DistinctTwists.size(); i++) {
      MakeTwoCopies[i] = false;
      int ti = DistinctTwists[i];
      PosType twist_i = TwistAngles[ti];
      for (int j=0; j<copyTwists.size(); j++) {
	int tj = copyTwists[j];
	PosType twist_j = TwistAngles[tj];
	if (TwistPair(twist_i, twist_j))
	  MakeTwoCopies[i] = true;
      }
      if (myComm->rank() == 0)
	fprintf (stderr, "Using %d copies of twist angle [%6.3f, %6.3f, %6.3f]\n",
		 MakeTwoCopies[i] ? 2 : 1, twist_i[0], twist_i[1], twist_i[2]);
    }
    
    // Find out if we can make real orbitals
    UseRealOrbitals = true;
    for (int i=0; i < DistinctTwists.size(); i++) {
      int ti = DistinctTwists[i];
      PosType twist = TwistAngles[ti];
      for (int j=0; j<OHMMS_DIM; j++)
	if (std::fabs(twist[j]-0.0) > 1.0e-8 &&
	    std::fabs(twist[j]-0.5) > 1.0e-8 &&
	    std::fabs(twist[j]+0.5) > 1.0e-8)
	  UseRealOrbitals = false;
    }

    if (UseRealOrbitals && (DistinctTwists.size() > 1)) {
      app_log() << "***** Use of real orbitals is possible, but not currently implemented\n"
		<< "      with more than one twist angle.\n";
      UseRealOrbitals = false;
    }
    
    if (UseRealOrbitals) 
      app_log() << "Using real orbitals.\n";
    else
      app_log() << "Using complex orbitals.\n";
    
#else
    DistinctTwists.resize(IncludeTwists.size());
    MakeTwoCopies.resize(IncludeTwists.size());
    for (int i=0; i<IncludeTwists.size(); i++) {
      DistinctTwists[i] = IncludeTwists[i];
      MakeTwoCopies[i] = false;
    }

#endif
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
    
    app_log() << "  Including twist vectors:\n";
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
	  app_log() << "tIndex = (" << tIndex[0] << ", " << tIndex[1] << ", "
		    << tIndex[2] << ")\n";
 	  // fprintf (stderr, "tIndex = (%d, %d, %d)  ti = %d\n", 
 	  // 	   tIndex[0], tIndex[1], tIndex[2], ti);

	  char buff[100];
	  snprintf (buff, 100, "    (%6.3f %6.3f %6.3f)\n", 
		    TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
	  app_log() << buff;
	}
  }


  void
  EinsplineSetBuilder::OccupyBands_ESHDF(int spin, bool sortBands)
  {
    if (myComm->rank() != 0) 
      return;

    SortBands.clear();
    for (int ti=0; ti<DistinctTwists.size(); ti++) {
      int tindex = DistinctTwists[ti];
      // First, read valence states
      ostringstream ePath;
      ePath << "/electrons/kpoint_" << tindex << "/spin_" 
	    << spin << "/eigenvalues";
      vector<double> eigvals;
      HDFAttribIO<vector<double> > h_eigvals(eigvals);
      h_eigvals.read(H5FileID, ePath.str().c_str());
      for (int bi=0; bi<NumBands; bi++) {
	BandInfo band;
	band.IsCoreState = false;
	band.TwistIndex = tindex;
	band.BandIndex  = bi;
	band.MakeTwoCopies = MakeTwoCopies[ti];
	band.Energy = eigvals[bi];
	if (band.Energy > -1.0e100) 
	  SortBands.push_back(band);
      }
      // Now, read core states
      for (int cs=0; cs<NumCoreStates; cs++) {
	BandInfo band;
	band.IsCoreState = true;
	band.TwistIndex = tindex;
	band.BandIndex  = cs;
	band.MakeTwoCopies = MakeTwoCopies[ti];
	HDFAttribIO<double> h_energy(band.Energy);
	
	h_energy.read   (H5FileID, (CoreStatePath(ti,cs)+"eigenvalue").c_str());
	if (band.Energy > -1.0e100) 
	  SortBands.push_back(band);
      }
    }
    // Now sort the bands by energy
    if (sortBands) {
      app_log() << "Sorting the bands now:\n";
      sort (SortBands.begin(), SortBands.end());
    }
    
    int orbIndex = 0;
    int numOrbs = 0;
    NumValenceOrbs=0;
    NumCoreOrbs=0;
    while (numOrbs < OrbitalSet->getOrbitalSetSize()) {
      if (SortBands[orbIndex].MakeTwoCopies)
	numOrbs += 2;
      else
	numOrbs++;
      if (SortBands[orbIndex].IsCoreState)
	NumCoreOrbs++;
      else
	NumValenceOrbs++;
      orbIndex++;
    }
    NumDistinctOrbitals = orbIndex;
    app_log() << "We will read " << NumDistinctOrbitals << " distinct orbitals.\n";
    app_log() << "There are " << NumCoreOrbs << " core states and " 
	      << NumValenceOrbs << " valence states.\n";
  }

  

  void
  EinsplineSetBuilder::OccupyBands(int spin, bool sortBands)
  {
    if (myComm->rank() != 0) 
      return;

    if (Format == ESHDF) {
      OccupyBands_ESHDF (spin, sortBands);
      return;
    }

    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    SortBands.clear();
    for (int ti=0; ti<DistinctTwists.size(); ti++) {
      int tindex = DistinctTwists[ti];
      // First, read valence states
      for (int bi=0; bi<NumBands; bi++) {
	BandInfo band;
	band.IsCoreState = false;
	band.TwistIndex = tindex;
	band.BandIndex  = bi;
	band.MakeTwoCopies = MakeTwoCopies[ti];
	
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
      // Now, read core states
      for (int cs=0; cs<NumCoreStates; cs++) {
	BandInfo band;
	band.IsCoreState = true;
	band.TwistIndex = tindex;
	band.BandIndex  = cs;
	band.MakeTwoCopies = MakeTwoCopies[ti];
	HDFAttribIO<double> h_energy(band.Energy);
	
	h_energy.read   (H5FileID, (CoreStatePath(ti,cs)+"eigenvalue").c_str());
	if (band.Energy > -1.0e100) 
	  SortBands.push_back(band);
      }
    }
    // Now sort the bands by energy
    if (sortBands) {
      app_log() << "Sorting the bands now:\n";
      sort (SortBands.begin(), SortBands.end());
    }
    
    int orbIndex = 0;
    int numOrbs = 0;
    NumValenceOrbs=0;
    NumCoreOrbs=0;
    while (numOrbs < OrbitalSet->getOrbitalSetSize()) {
      if (SortBands[orbIndex].MakeTwoCopies)
	numOrbs += 2;
      else
	numOrbs++;
      if (SortBands[orbIndex].IsCoreState)
	NumCoreOrbs++;
      else
	NumValenceOrbs++;
      orbIndex++;
    }
    NumDistinctOrbitals = orbIndex;
    app_log() << "We will read " << NumDistinctOrbitals << " distinct orbitals.\n";
    app_log() << "There are " << NumCoreOrbs << " core states and " 
	      << NumValenceOrbs << " valence states.\n";
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
	app_log() << "Making " << orb->uCenters.size() << " copies of band "
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
    bool root = myComm->rank()==0;
    // bcast other stuff
    myComm->bcast (NumDistinctOrbitals);
    myComm->bcast (NumValenceOrbs);
    myComm->bcast (NumCoreOrbs);
    int N = NumDistinctOrbitals;

    orbitalSet->kPoints.resize(N);
    orbitalSet->MakeTwoCopies.resize(N);
    orbitalSet->StorageValueVector.resize(N); orbitalSet->BlendValueVector.resize(N);
    orbitalSet->StorageLaplVector.resize(N);  orbitalSet->BlendLaplVector.resize(N);
    orbitalSet->StorageGradVector.resize(N);  orbitalSet->BlendGradVector.resize(N);
    orbitalSet->StorageHessVector.resize(N);
    orbitalSet->phase.resize(N);
    orbitalSet->eikr.resize(N);
    orbitalSet->NumValenceOrbs = NumValenceOrbs;
    orbitalSet->NumCoreOrbs    = NumCoreOrbs;
    // Read in k-points
    int numOrbs = orbitalSet->getOrbitalSetSize();
    int num = 0;
    if (root) {
      for (int iorb=0; iorb<N; iorb++) {
	int ti = SortBands[iorb].TwistIndex;
	PosType twist  = TwistAngles[ti];
	orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
	orbitalSet->MakeTwoCopies[iorb] = 
	  (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
	num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
      }
      PosType twist0 = TwistAngles[SortBands[0].TwistIndex];
      for (int i=0; i<OHMMS_DIM; i++)
	if (std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
	  orbitalSet->HalfG[i] = 1;
	else
	  orbitalSet->HalfG[i] = 0;
    }
    myComm->bcast(orbitalSet->kPoints);
    myComm->bcast(orbitalSet->MakeTwoCopies);
    myComm->bcast(orbitalSet->HalfG);

    // First, check to see if we have already read this in
    H5OrbSet set(H5FileName, spin, N);
    
    int nx, ny, nz, bi, ti;
    Array<complex<double>,3> rawData;
    Array<double,3>         splineData;

    if (root) {
      // Find the orbital mesh size
      int i=0;
      while (SortBands[i].IsCoreState) i++;
      ti = SortBands[i].TwistIndex;
      bi = SortBands[i].BandIndex;
      string vectorName = OrbitalPath (ti, bi) + "eigenvector";
      HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
      h_rawData.read(H5FileID, vectorName.c_str());
      nx = rawData.size(0); 
      ny = rawData.size(1);
      nz = rawData.size(2);
      splineData.resize(nx-1, ny-1, nz-1);
      PosType ru;
      for (int ix=0; ix<(nx-1); ix++) {
	ru[0] = (RealType)ix / (RealType)(nx-1);
	for (int iy=0; iy<(ny-1); iy++) {
	  ru[1] = (RealType)iy / (RealType)(ny-1);
	  for (int iz=0; iz<(nz-1); iz++) {
	    ru[2] = (RealType)iz / (RealType)(nz-1);
	    double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
	    double s, c;
	    sincos(phi, &s, &c);
	    complex<double> phase(c,s);
	    complex<double> z = phase*rawData(ix,iy,iz);
	    splineData(ix,iy,iz) = z.imag();
	  }
	}
      }

      PosType twist, k;
      twist = TwistAngles[ti];
      k = orbitalSet->PrimLattice.k_cart(twist);
      double e = SortBands[i].Energy;
    }

    TinyVector<int,3> nxyz(nx,ny,nz);
    myComm->bcast(nxyz);
    nx=nxyz[0]; ny=nxyz[1]; nz=nxyz[2];

    if (!root) splineData.resize(nx-1,ny-1,nz-1);

    myComm->bcast(splineData);
    Ugrid x_grid, y_grid, z_grid;
    BCtype_d xBC, yBC, zBC;

    if (orbitalSet->HalfG[0]) 
      { xBC.lCode = ANTIPERIODIC;    xBC.rCode = ANTIPERIODIC; }
    else
      { xBC.lCode = PERIODIC;        xBC.rCode = PERIODIC; }

    if (orbitalSet->HalfG[1]) 
      { yBC.lCode = ANTIPERIODIC;    yBC.rCode = ANTIPERIODIC; }
    else
      { yBC.lCode = PERIODIC;        yBC.rCode = PERIODIC; }

    if (orbitalSet->HalfG[2]) 
      { zBC.lCode = ANTIPERIODIC;    zBC.rCode = ANTIPERIODIC; }
    else
      { zBC.lCode = PERIODIC;        zBC.rCode = PERIODIC; }

    x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = nx-1;
    y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = ny-1;
    z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = nz-1;

    // Create the multiUBspline object
    orbitalSet->MultiSpline = 
      create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
    
   set_multi_UBspline_3d_d (orbitalSet->MultiSpline, 0, splineData.data());

    //////////////////////////////////////
    // Create the MuffinTin APW splines //
    //////////////////////////////////////
    orbitalSet->MuffinTins.resize(NumMuffinTins);
    for (int tin=0; tin<NumMuffinTins; tin++) {
      orbitalSet->MuffinTins[tin].Atom = tin;
      orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
      orbitalSet->MuffinTins[tin].set_lattice(Lattice);
      orbitalSet->MuffinTins[tin].init_APW 
	(MT_APW_rgrids[tin], MT_APW_lmax[tin], 
	 NumValenceOrbs);
    }
           
    int iorb  = 0;
    int icore = 0;
    int ival = 0;
    while (iorb < N) {
      bool isCore;
      if (root)  isCore = SortBands[iorb].IsCoreState;
      myComm->bcast (isCore);
      if (isCore) {
	int atom, l, m=0;
	double rmax;
	Vector<double> g, r;
	PosType twist, k;
	if (root) {
	  ti   = SortBands[iorb].TwistIndex;
	  bi   = SortBands[iorb].BandIndex;
	  double e = SortBands[iorb].Energy;
	  twist = TwistAngles[ti];
	  k = orbitalSet->PrimLattice.k_cart(twist);
	  string atomName = CoreStatePath (ti, bi) + "atom";
	  string gName    = CoreStatePath (ti, bi) + "g";
	  string rMaxName = CoreStatePath (ti, bi) + "rmax";
	  string lName    = CoreStatePath (ti, bi) + "l";
	  string kName    = CoreStatePath (ti, bi) + "k";
	  string rName    = CoreStatePath (ti, bi) + "r";
	  HDFAttribIO<int> h_atom(atom), h_l(l);
	  HDFAttribIO<double> h_rmax(rmax);
	  HDFAttribIO<Vector<double> > h_g(g);
	  HDFAttribIO<Vector<double> > h_r(r);
	  h_atom.read (H5FileID, atomName.c_str());
	  h_l.read    (H5FileID,    lName.c_str());
	  h_rmax.read (H5FileID, rMaxName.c_str());
	  h_g.read    (H5FileID,   gName.c_str());
	  h_r.read    (H5FileID,   rName.c_str());
	  
	  fprintf (stderr, "  Core state:     ti=%3d  bi=%3d energy=%8.5f "
		   "k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
		   ti, bi, e, k[0], k[1], k[2], myComm->rank());
	}
	myComm->bcast (atom);  myComm->bcast(rmax);	
	myComm->bcast (l);     myComm->bcast (k);
	int ng = g.size();     myComm->bcast(ng);
	if (g.size() != ng) {
	  g.resize(ng);
	  r.resize(ng);
	}
	myComm->bcast (g);
	myComm->bcast (r);

	double Z = (double)IonTypes(atom);
	OrbitalSet->MuffinTins[atom].addCore (l, m, r, g, k, Z);
	icore++;
      }
      else {
	if (root) {
	  int ti   = SortBands[iorb].TwistIndex;
	  int bi   = SortBands[iorb].BandIndex;
	  double e = SortBands[iorb].Energy;
	  PosType twist, k;
	  twist = TwistAngles[ti];
	  k = orbitalSet->PrimLattice.k_cart(twist);
	  fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
		   ti, bi, e, k[0], k[1], k[2], myComm->rank());
	  
	  string vectorName = OrbitalPath (ti, bi) + "eigenvector";
	  HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
	  h_rawData.read(H5FileID, vectorName.c_str());
	  if ((rawData.size(0) != nx) ||
	      (rawData.size(1) != ny) ||
	      (rawData.size(2) != nz)) {
	    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
	    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
	    abort();
	  }
	  PosType ru;
	  for (int ix=0; ix<(nx-1); ix++) {
	    ru[0] = (RealType)ix / (RealType)(nx-1);
	    for (int iy=0; iy<(ny-1); iy++) {
	      ru[1] = (RealType)iy / (RealType)(ny-1);
	      for (int iz=0; iz<(nz-1); iz++) {
		ru[2] = (RealType)iz / (RealType)(nz-1);
		double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
		double s, c;
		sincos(phi, &s, &c);
		complex<double> phase(c,s);
		complex<double> z = phase*rawData(ix,iy,iz);
		splineData(ix,iy,iz) = z.real();
	      }
	    }
	  }
	}
	myComm->bcast(splineData);
	set_multi_UBspline_3d_d 
	  (orbitalSet->MultiSpline, ival, splineData.data());
	
	// Now read muffin tin data
	for (int tin=0; tin<NumMuffinTins; tin++) {
	  // app_log() << "Reading data for muffin tin " << tin << endl;
	  PosType twist, k;
	  int lmax = MT_APW_lmax[tin];
	  int numYlm = (lmax+1)*(lmax+1);
	  Array<complex<double>,2> 
	    u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
	  Array<complex<double>,1> du_lm_dr (numYlm);
	  if (root) {
	    int ti   = SortBands[iorb].TwistIndex;
	    int bi   = SortBands[iorb].BandIndex;
	    twist = TwistAngles[ti];
	    k = orbitalSet->PrimLattice.k_cart(twist);
	    string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
	    string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
	    HDFAttribIO<Array<complex<double>,2> > h_u_lm_r(u_lm_r);
	    HDFAttribIO<Array<complex<double>,1> > h_du_lm_dr(du_lm_dr);
	    h_u_lm_r.read(H5FileID, uName.c_str());
	    h_du_lm_dr.read(H5FileID, duName.c_str());
	  }
	  myComm->bcast(u_lm_r);
	  myComm->bcast(du_lm_dr);
	  myComm->bcast(k);
	  double Z = (double)IonTypes(tin);
	  OrbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
	}
	ival++;
      } // valence state
      iorb++;
    }
    ExtendedMap_d[set] = orbitalSet->MultiSpline;
  }


  void
  EinsplineSetBuilder::ReadBands_ESHDF
  (int spin, EinsplineSetExtended<complex<double > >* orbitalSet)
  {
    bool root = myComm->rank()==0;
    // bcast other stuff
    myComm->bcast (NumDistinctOrbitals);
    myComm->bcast (NumValenceOrbs);
    myComm->bcast (NumCoreOrbs);
    int N = NumDistinctOrbitals;

    orbitalSet->kPoints.resize(N);
    orbitalSet->MakeTwoCopies.resize(N);
    orbitalSet->StorageValueVector.resize(N); orbitalSet->BlendValueVector.resize(N);
    orbitalSet->StorageLaplVector.resize(N);  orbitalSet->BlendLaplVector.resize(N);
    orbitalSet->StorageGradVector.resize(N);  orbitalSet->BlendGradVector.resize(N);
    orbitalSet->StorageHessVector.resize(N);
    orbitalSet->phase.resize(N);
    orbitalSet->eikr.resize(N);
    orbitalSet->NumValenceOrbs = NumValenceOrbs;
    orbitalSet->NumCoreOrbs    = NumCoreOrbs;
    // Read in k-points
    int numOrbs = orbitalSet->getOrbitalSetSize();
    int num = 0;
    if (root) {
      for (int iorb=0; iorb<N; iorb++) {
	int ti = SortBands[iorb].TwistIndex;
	PosType twist  = TwistAngles[ti];
	orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
	orbitalSet->MakeTwoCopies[iorb] = 
	  (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
	num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
      }
    }
    myComm->bcast(orbitalSet->kPoints);
    myComm->bcast(orbitalSet->MakeTwoCopies);

    // First, check to see if we have already read this in
    H5OrbSet set(H5FileName, spin, N);
    
    int nx, ny, nz, bi, ti;
    Array<complex<double>,3> splineData;
    TinyVector<int,3> mesh;

    // Find the orbital mesh size
    if (root) {
      HDFAttribIO<TinyVector<int,3> > h_mesh(mesh);
      h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
      h_mesh.read (H5FileID, "/electrons/mesh");
    }
    myComm->bcast(mesh);
    nx=mesh[0]; ny=mesh[1]; nz=mesh[2];
    splineData.resize(nx,ny,nz);

    Ugrid x_grid, y_grid, z_grid;
    BCtype_z xBC, yBC, zBC;

    xBC.lCode = PERIODIC;        xBC.rCode = PERIODIC;
    yBC.lCode = PERIODIC;        yBC.rCode = PERIODIC;
    zBC.lCode = PERIODIC;        zBC.rCode = PERIODIC;
    x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = nx;
    y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = ny;
    z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = nz;

    // Create the multiUBspline object
    orbitalSet->MultiSpline = 
      create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
    
    //////////////////////////////////////
    // Create the MuffinTin APW splines //
    //////////////////////////////////////
    orbitalSet->MuffinTins.resize(NumMuffinTins);
    for (int tin=0; tin<NumMuffinTins; tin++) {
      orbitalSet->MuffinTins[tin].Atom = tin;
      orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
      orbitalSet->MuffinTins[tin].set_lattice(Lattice);
      orbitalSet->MuffinTins[tin].init_APW 
	(MT_APW_rgrids[tin], MT_APW_lmax[tin], 
	 NumValenceOrbs);
    }
    for (int iat=0; iat<AtomicOrbitals.size(); iat++) {
      AtomicOrbitals[iat].set_num_bands(NumValenceOrbs);
      AtomicOrbitals[iat].allocate();
    }
      
           
    int iorb  = 0;
    int icore = 0;
    int ival = 0;

    int isComplex;
    if (root) {
      HDFAttribIO<int> h_isComplex(isComplex);
      h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");

      if (!isComplex) {
	app_error() << "Expected complex orbitals in ES-HDF file, but found real ones.\n";
	abort();
      }
    }

    while (iorb < N) {
      bool isCore;
      if (root)  isCore = SortBands[iorb].IsCoreState;
      myComm->bcast (isCore);
      if (isCore) {
	app_error() << "Core states not supported by ES-HDF yet.\n";
	abort();
      }
      else {
	PosType twist;
	if (root) {
	  int ti   = SortBands[iorb].TwistIndex;
	  int bi   = SortBands[iorb].BandIndex;
	  double e = SortBands[iorb].Energy;
	  PosType k;
	  twist = TwistAngles[ti];
	  k = orbitalSet->PrimLattice.k_cart(twist);
	  fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
		   ti, bi, e, k[0], k[1], k[2], myComm->rank());

	  ostringstream path;
	  path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
	  string psiName = path.str() + "psi_r";

	  HDFAttribIO<Array<complex<double>,3> >  h_splineData(splineData);
	  h_splineData.read(H5FileID, psiName.c_str());
	  if ((splineData.size(0) != nx) ||
	      (splineData.size(1) != ny) ||
	      (splineData.size(2) != nz)) {
	    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
	    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
	    abort();
	  }
	}
	myComm->bcast(splineData);
	myComm->bcast(twist);
	set_multi_UBspline_3d_z 
	  (orbitalSet->MultiSpline, ival, splineData.data());
      
	// Read atomic orbital information
	for (int iat=0; iat<AtomicOrbitals.size(); iat++) {
	  app_log() << "Reading orbital " << iat << " for band " << ival << endl;
	  AtomicOrbital<complex<double> > &orb = AtomicOrbitals[iat];
	  Array<complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm), 
	    poly_coefs(orb.PolyOrder+1,orb.Numlm);
	  if (root) { 
	    int ti   = SortBands[iorb].TwistIndex;
	    int bi   = SortBands[iorb].BandIndex;
	    ostringstream path;
	    path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
	    AtomicOrbital<complex<double> > &orb = AtomicOrbitals[iat];
	    ostringstream spline_path, poly_path;
	    spline_path << path.str() << "radial_spline_" << iat;
	    poly_path   << path.str() << "poly_coefs_"    << iat;
	    HDFAttribIO<Array<complex<double>,2> > h_radial_spline(radial_spline);
	    HDFAttribIO<Array<complex<double>,2> > h_poly_coefs(poly_coefs);
	    h_radial_spline.read(H5FileID, spline_path.str().c_str());
	    h_poly_coefs.read   (H5FileID, poly_path.str().c_str());
	    // cerr << "radial_spline.size = (" << radial_spline.size(0) 
	    // 	 << ", " << radial_spline.size(1) << ")\n";
	    // cerr << "poly_coefs.size = (" << poly_coefs.size(0) 
	    // 	 << ", " << poly_coefs.size(1) << ")\n";
	  }
	  myComm->bcast(radial_spline);
	  myComm->bcast(poly_coefs);
	  AtomicOrbitals[iat].set_band (ival, radial_spline, poly_coefs, twist);
	}

	// Now read muffin tin data
	for (int tin=0; tin<NumMuffinTins; tin++) {
	  // app_log() << "Reading data for muffin tin " << tin << endl;
	  PosType twist, k;
	  int lmax = MT_APW_lmax[tin];
	  int numYlm = (lmax+1)*(lmax+1);
	  Array<complex<double>,2> 
	    u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
	  Array<complex<double>,1> du_lm_dr (numYlm);
	  if (root) {
	    int ti   = SortBands[iorb].TwistIndex;
	    int bi   = SortBands[iorb].BandIndex;
	    twist = TwistAngles[ti];
	    k = orbitalSet->PrimLattice.k_cart(twist);
	    string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
	    string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
	    HDFAttribIO<Array<complex<double>,2> > h_u_lm_r(u_lm_r);
	    HDFAttribIO<Array<complex<double>,1> > h_du_lm_dr(du_lm_dr);
	    h_u_lm_r.read(H5FileID, uName.c_str());
	    h_du_lm_dr.read(H5FileID, duName.c_str());
	  }
	  myComm->bcast(u_lm_r);
	  myComm->bcast(du_lm_dr);
	  myComm->bcast(k);
	  double Z = (double)IonTypes(tin);
	  OrbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
	}
	ival++;
      } // valence state
      iorb++;
    }
    orbitalSet->AtomicOrbitals = AtomicOrbitals;
    for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
      orbitalSet->AtomicOrbitals[i].registerTimers();

    ExtendedMap_z[set] = orbitalSet->MultiSpline;
  }



  void
  EinsplineSetBuilder::ReadBands_ESHDF
  (int spin, EinsplineSetExtended<double>* orbitalSet)
  {
    vector<AtomicOrbital<double> > realOrbs(AtomicOrbitals.size());
    for (int iat=0; iat<realOrbs.size(); iat++) {
      AtomicOrbital<complex<double> > &corb (AtomicOrbitals[iat]);
      realOrbs[iat].set_pos  (corb.Pos);
      realOrbs[iat].set_lmax (corb.lMax);
      realOrbs[iat].set_cutoff (corb.CutoffRadius);
      realOrbs[iat].set_spline (corb.SplineRadius, corb.SplinePoints);
      realOrbs[iat].set_polynomial (corb.PolyRadius, corb.PolyOrder);
      realOrbs[iat].Lattice = corb.Lattice;
    }

    bool root = myComm->rank()==0;
    // bcast other stuff
    myComm->bcast (NumDistinctOrbitals);
    myComm->bcast (NumValenceOrbs);
    myComm->bcast (NumCoreOrbs);
    int N = NumDistinctOrbitals;

    orbitalSet->kPoints.resize(N);
    orbitalSet->MakeTwoCopies.resize(N);
    orbitalSet->StorageValueVector.resize(N); orbitalSet->BlendValueVector.resize(N);
    orbitalSet->StorageLaplVector.resize(N);  orbitalSet->BlendLaplVector.resize(N);
    orbitalSet->StorageGradVector.resize(N);  orbitalSet->BlendGradVector.resize(N);
    orbitalSet->StorageHessVector.resize(N);
    orbitalSet->phase.resize(N);
    orbitalSet->eikr.resize(N);
    orbitalSet->NumValenceOrbs = NumValenceOrbs;
    orbitalSet->NumCoreOrbs    = NumCoreOrbs;
    orbitalSet->FirstOrderSplines.resize(IonPos.size());
    // Read in k-points
    int numOrbs = orbitalSet->getOrbitalSetSize();
    int num = 0;
    if (root) {
      for (int iorb=0; iorb<N; iorb++) {
	int ti = SortBands[iorb].TwistIndex;
	PosType twist  = TwistAngles[ti];
	orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
	orbitalSet->MakeTwoCopies[iorb] = 
	  (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
	num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
      }
      PosType twist0 = TwistAngles[SortBands[0].TwistIndex];
      for (int i=0; i<OHMMS_DIM; i++)
	if (std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
	  orbitalSet->HalfG[i] = 1;
	else
	  orbitalSet->HalfG[i] = 0;
    }
    myComm->bcast(orbitalSet->kPoints);
    myComm->bcast(orbitalSet->MakeTwoCopies);
    myComm->bcast(orbitalSet->HalfG);

    // First, check to see if we have already read this in
    H5OrbSet set(H5FileName, spin, N);
    
    int nx, ny, nz, bi, ti;
    Array<complex<double>,3> rawData;
    Array<double,3>         splineData;
    TinyVector<int,3> mesh;

    // Find the orbital mesh size
    if (root) {
      HDFAttribIO<TinyVector<int,3> > h_mesh(mesh);
      h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
      h_mesh.read (H5FileID, "/electrons/mesh");
    }
    myComm->bcast(mesh);
    nx=mesh[0]; ny=mesh[1]; nz=mesh[2];
    splineData.resize(nx,ny,nz);

    Ugrid x_grid, y_grid, z_grid;
    BCtype_d xBC, yBC, zBC;

    if (orbitalSet->HalfG[0]) 
      { xBC.lCode = ANTIPERIODIC;    xBC.rCode = ANTIPERIODIC; }
    else
      { xBC.lCode = PERIODIC;        xBC.rCode = PERIODIC; }

    if (orbitalSet->HalfG[1]) 
      { yBC.lCode = ANTIPERIODIC;    yBC.rCode = ANTIPERIODIC; }
    else
      { yBC.lCode = PERIODIC;        yBC.rCode = PERIODIC; }

    if (orbitalSet->HalfG[2]) 
      { zBC.lCode = ANTIPERIODIC;    zBC.rCode = ANTIPERIODIC; }
    else
      { zBC.lCode = PERIODIC;        zBC.rCode = PERIODIC; }

    x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = nx;
    y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = ny;
    z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = nz;

    // Create the multiUBspline object
    orbitalSet->MultiSpline = 
      create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);

    if (HaveOrbDerivs) {
      orbitalSet->FirstOrderSplines.resize(IonPos.size());
      for (int ion=0; ion<IonPos.size(); ion++)
	for (int dir=0; dir<OHMMS_DIM; dir++)
	  orbitalSet->FirstOrderSplines[ion][dir] = 
	    create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
    }
    //////////////////////////////////////
    // Create the MuffinTin APW splines //
    //////////////////////////////////////
    orbitalSet->MuffinTins.resize(NumMuffinTins);
    for (int tin=0; tin<NumMuffinTins; tin++) {
      orbitalSet->MuffinTins[tin].Atom = tin;
      orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
      orbitalSet->MuffinTins[tin].set_lattice(Lattice);
      orbitalSet->MuffinTins[tin].init_APW 
	(MT_APW_rgrids[tin], MT_APW_lmax[tin], 
	 NumValenceOrbs);
    }

    for (int iat=0; iat<realOrbs.size(); iat++) {
      realOrbs[iat].set_num_bands(NumValenceOrbs);
      realOrbs[iat].allocate();
    }
           
    int iorb  = 0;
    int icore = 0;
    int ival = 0;

    int isComplex;
    if (root) {
      HDFAttribIO<int> h_isComplex(isComplex);
      h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");
    }
    myComm->bcast(isComplex);

    while (iorb < N) {
      bool isCore;
      if (root)  isCore = SortBands[iorb].IsCoreState;
      myComm->bcast (isCore);
      if (isCore) {
	app_error() << "Core states not supported by ES-HDF yet.\n";
	abort();
      }
      else {
	PosType twist;
	if (root) {
	  int ti   = SortBands[iorb].TwistIndex;
	  int bi   = SortBands[iorb].BandIndex;
	  double e = SortBands[iorb].Energy;
	  PosType k;
	  twist = TwistAngles[ti];
	  k = orbitalSet->PrimLattice.k_cart(twist);
	  fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
		   ti, bi, e, k[0], k[1], k[2], myComm->rank());

	  ostringstream path;
	  path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
	  string psiName = path.str() + "psi_r";

	  if (isComplex) {
	    HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
	    h_rawData.read(H5FileID, psiName.c_str());
	    if ((rawData.size(0) != nx) ||
		(rawData.size(1) != ny) ||
		(rawData.size(2) != nz)) {
	      fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
	      fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
	      abort();
	    }
	    PosType ru;
	    for (int ix=0; ix<nx; ix++) {
	      ru[0] = (RealType)ix / (RealType)nx;
	      for (int iy=0; iy<ny; iy++) {
		ru[1] = (RealType)iy / (RealType)ny;
		for (int iz=0; iz<nz; iz++) {
		  ru[2] = (RealType)iz / (RealType)nz;
		  double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
		  double s, c;
		  sincos(phi, &s, &c);
		  complex<double> phase(c,s);
		  complex<double> z = phase*rawData(ix,iy,iz);
		  splineData(ix,iy,iz) = z.real();
		}
	      }
	    }
	  }
	  else {
	    HDFAttribIO<Array<double,3> >  h_splineData(splineData);
	    h_splineData.read(H5FileID, psiName.c_str());
	    if ((splineData.size(0) != nx) ||
		(splineData.size(1) != ny) ||
		(splineData.size(2) != nz)) {
	      fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
	      fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
	      abort();
	    }
	  }
	}
	myComm->bcast(twist);
	myComm->bcast(splineData);
	set_multi_UBspline_3d_d 
	  (orbitalSet->MultiSpline, ival, splineData.data());

	// Read atomic orbital information
	for (int iat=0; iat<realOrbs.size(); iat++) {
	  app_log() << "Reading orbital " << iat << " for band " << ival << endl;
	  AtomicOrbital<double> &orb = realOrbs[iat];
	  //AtomicOrbital<complex<double> > &orb = realOrbs[iat];
	  Array<complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm), 
	    poly_coefs(orb.PolyOrder+1,orb.Numlm);
	  if (root) { 
	    int ti   = SortBands[iorb].TwistIndex;
	    int bi   = SortBands[iorb].BandIndex;
	    ostringstream path;
	    path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
	    ostringstream spline_path, poly_path;
	    spline_path << path.str() << "radial_spline_" << iat;
	    poly_path   << path.str() << "poly_coefs_"    << iat;
	    HDFAttribIO<Array<complex<double>,2> > h_radial_spline(radial_spline);
	    HDFAttribIO<Array<complex<double>,2> > h_poly_coefs(poly_coefs);
	    h_radial_spline.read(H5FileID, spline_path.str().c_str());
	    h_poly_coefs.read   (H5FileID, poly_path.str().c_str());
	  }
	  myComm->bcast(radial_spline);
	  myComm->bcast(poly_coefs);
	  realOrbs[iat].set_band (ival, radial_spline, poly_coefs, twist);
	}

      
	// Now read orbital derivatives if we have them
	if (HaveOrbDerivs) {
	  for (int ion=0; ion<IonPos.size(); ion++) 
	    for (int dim=0; dim<OHMMS_DIM; dim++) {
	      if (root) {
		int ti   = SortBands[iorb].TwistIndex;
		int bi   = SortBands[iorb].BandIndex;
		
		app_log() << "Reading orbital derivative for ion " << ion 
			  << " dim " << dim << " spin " << spin << " band "
			  << bi << " kpoint " << ti << endl;
		ostringstream path;
		path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/"
		     << "dpsi_" << ion << "_" << dim << "_r";
		string psiName = path.str();
		if (isComplex) {
		  HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
		  h_rawData.read(H5FileID, psiName.c_str());
		  if ((rawData.size(0) != nx) ||
		      (rawData.size(1) != ny) ||
		      (rawData.size(2) != nz)) {
		    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
		    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
		    abort();
		  }
		  PosType ru;
		  for (int ix=0; ix<nx; ix++) {
		    ru[0] = (RealType)ix / (RealType)nx;
		    for (int iy=0; iy<ny; iy++) {
		      ru[1] = (RealType)iy / (RealType)ny;
		      for (int iz=0; iz<nz; iz++) {
			ru[2] = (RealType)iz / (RealType)nz;
			double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
			double s, c;
			sincos(phi, &s, &c);
			complex<double> phase(c,s);
			complex<double> z = phase*rawData(ix,iy,iz);
			splineData(ix,iy,iz) = z.real();
		      }
		    }
		  }
		}
		else {
		  HDFAttribIO<Array<double,3> >  h_splineData(splineData);
		  h_splineData.read(H5FileID, psiName.c_str());
		  if ((splineData.size(0) != nx) ||
		      (splineData.size(1) != ny) ||
		      (splineData.size(2) != nz)) {
		    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
		    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
		    abort();
		  }
		}
	      }
	      myComm->bcast(splineData);
	      set_multi_UBspline_3d_d 
		(orbitalSet->FirstOrderSplines[ion][dim], ival, splineData.data());
	    }
	}
	

	
	// Now read muffin tin data
	for (int tin=0; tin<NumMuffinTins; tin++) {
	  // app_log() << "Reading data for muffin tin " << tin << endl;
	  PosType twist, k;
	  int lmax = MT_APW_lmax[tin];
	  int numYlm = (lmax+1)*(lmax+1);
	  Array<complex<double>,2> 
	    u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
	  Array<complex<double>,1> du_lm_dr (numYlm);
	  if (root) {
	    int ti   = SortBands[iorb].TwistIndex;
	    int bi   = SortBands[iorb].BandIndex;
	    twist = TwistAngles[ti];
	    k = orbitalSet->PrimLattice.k_cart(twist);
	    string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
	    string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
	    HDFAttribIO<Array<complex<double>,2> > h_u_lm_r(u_lm_r);
	    HDFAttribIO<Array<complex<double>,1> > h_du_lm_dr(du_lm_dr);
	    h_u_lm_r.read(H5FileID, uName.c_str());
	    h_du_lm_dr.read(H5FileID, duName.c_str());
	  }
	  myComm->bcast(u_lm_r);
	  myComm->bcast(du_lm_dr);
	  myComm->bcast(k);
	  double Z = (double)IonTypes(tin);
	  OrbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
	}
	ival++;
      } // valence state
      iorb++;
    }
    orbitalSet->AtomicOrbitals = realOrbs;
    for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
      orbitalSet->AtomicOrbitals[i].registerTimers();

    ExtendedMap_d[set] = orbitalSet->MultiSpline;
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

  string
  EinsplineSetBuilder::CoreStatePath(int ti, int cs)
  {
    string eigenstatesGroup;
    if (Version[0]==0 && Version[1]== 11) 
      eigenstatesGroup = "/eigenstates_3";
    else if (Version[0]==0 && Version[1]==20) 
      eigenstatesGroup = "/eigenstates";

    ostringstream groupPath;
    
    if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
      groupPath << eigenstatesGroup << "/twist_" 
		<< ti << "/core_state_" << cs << "/";
    else if (NumBands > 1)
      groupPath << eigenstatesGroup << "/twist/core_state_" << cs << "/";
    else 
      groupPath << eigenstatesGroup << "/twist/core_state/";

    return groupPath.str();
  }

  string
  EinsplineSetBuilder::MuffinTinPath(int ti, int bi, int tin)
  {
    ostringstream groupPath;
    if (NumMuffinTins > 0)
      groupPath << OrbitalPath(ti,bi) << "muffin_tin_" << tin << "/";
    else
      groupPath << OrbitalPath(ti,bi) << "muffin_tin/";
    return groupPath.str();
  }

  void
  EinsplineSetBuilder::ReadBands 
  (int spin, EinsplineSetExtended<complex<double> >* orbitalSet)
  {
    bool root = myComm->rank()==0;
    //bcastwith other stuff
    myComm->bcast(NumDistinctOrbitals);
    myComm->bcast (NumValenceOrbs);
    myComm->bcast (NumCoreOrbs);
    int N = NumDistinctOrbitals;

    orbitalSet->kPoints.resize(N);
    orbitalSet->MakeTwoCopies.resize(N);
    orbitalSet->StorageValueVector.resize(N);  orbitalSet->BlendValueVector.resize(N);
    orbitalSet->StorageLaplVector.resize(N);   orbitalSet->BlendLaplVector.resize(N);
    orbitalSet->StorageGradVector.resize(N);   orbitalSet->BlendGradVector.resize(N);
    orbitalSet->StorageHessVector.resize(N);
    orbitalSet->phase.resize(N);
    orbitalSet->eikr.resize(N);
    orbitalSet->NumValenceOrbs = NumValenceOrbs;
    orbitalSet->NumCoreOrbs    = NumCoreOrbs;

    if (root) {
      int numOrbs = orbitalSet->getOrbitalSetSize();
      int num = 0;
      for (int iorb=0; iorb<N; iorb++) {
	int ti = SortBands[iorb].TwistIndex;
	PosType twist  = TwistAngles[ti];
	orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
	orbitalSet->MakeTwoCopies[iorb] = 
	  (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
	num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
      }
    }
    myComm->bcast(orbitalSet->kPoints);
    myComm->bcast(orbitalSet->MakeTwoCopies);
    
    // First, check to see if we have already read this in
    H5OrbSet set(H5FileName, spin, N);
    // std::map<H5OrbSet,multi_UBspline_3d_z*>::iterator iter;
    // iter = ExtendedMap_z.find (set);
    // if (iter != ExtendedMap_z.end()) {
    //   cerr << "Using existing copy of multi_UBspline_3d_z for "
    // 	   << "thread number " << omp_get_thread_num() << ".\n";
    //   orbitalSet->MultiSpline = iter->second;
    //   return;
    // }

    int nx, ny, nz, bi, ti;
    Array<complex<double>,3> splineData, rawData;
    if (root) {
      // Find the orbital mesh size
      int i=0;
      while (SortBands[i].IsCoreState) i++;
      ti = SortBands[i].TwistIndex;
      bi = SortBands[i].BandIndex;
      string vectorName = OrbitalPath (ti, bi) + "eigenvector";
      HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
      h_rawData.read(H5FileID, vectorName.c_str());
      nx = rawData.size(0); 
      ny = rawData.size(1);
      nz = rawData.size(2);
      splineData.resize(nx-1, ny-1, nz-1);
      for (int ix=0; ix<(nx-1); ix++)
	for (int iy=0; iy<(ny-1); iy++)
	  for (int iz=0; iz<(nz-1); iz++)
	    splineData(ix,iy,iz) = rawData(ix,iy,iz);

      PosType twist, k;
      twist = TwistAngles[ti];
      k = orbitalSet->PrimLattice.k_cart(twist);
      double e = SortBands[i].Energy;
      // fprintf (stderr, "  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
      // 	       ti, bi, e, k[0], k[1], k[2], myComm->rank());   
    }

    TinyVector<int,3> nxyz(nx,ny,nz);
    myComm->bcast(nxyz);
    nx=nxyz[0]; ny=nxyz[1]; nz=nxyz[2];

    if (!root) splineData.resize(nx-1,ny-1,nz-1);

    myComm->bcast(splineData);
    Ugrid x_grid, y_grid, z_grid;
    BCtype_z xBC, yBC, zBC;

    xBC.lCode = PERIODIC;    xBC.rCode = PERIODIC;
    yBC.lCode = PERIODIC;    yBC.rCode = PERIODIC;
    zBC.lCode = PERIODIC;    zBC.rCode = PERIODIC;
    x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = nx-1;
    y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = ny-1;
    z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = nz-1;

    // Create the multiUBspline object
    orbitalSet->MultiSpline = 
      create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);

    set_multi_UBspline_3d_z (orbitalSet->MultiSpline, 0, splineData.data());

    //////////////////////////////////////
    // Create the MuffinTin APW splines //
    //////////////////////////////////////
    orbitalSet->MuffinTins.resize(NumMuffinTins);
    for (int tin=0; tin<NumMuffinTins; tin++) {
      orbitalSet->MuffinTins[tin].Atom = tin;
      orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
      orbitalSet->MuffinTins[tin].set_lattice(Lattice);
      orbitalSet->MuffinTins[tin].init_APW 
	(MT_APW_rgrids[tin], MT_APW_lmax[tin], 
	 NumValenceOrbs);
    }
           
    int iorb  = 0;
    int icore = 0;
    int ival = 0;
    while (iorb < N) {
      bool isCore;
      if (root)  isCore = SortBands[iorb].IsCoreState;
      myComm->bcast (isCore);
      if (isCore) {
	int atom, l, m=0;
	double rmax;
	Vector<double> g, r;
	PosType twist, k;
	if (root) {
	  ti   = SortBands[iorb].TwistIndex;
	  bi   = SortBands[iorb].BandIndex;
	  double e = SortBands[iorb].Energy;
	  twist = TwistAngles[ti];
	  k = orbitalSet->PrimLattice.k_cart(twist);
	  string atomName = CoreStatePath (ti, bi) + "atom";
	  string gName    = CoreStatePath (ti, bi) + "g";
	  string rMaxName = CoreStatePath (ti, bi) + "rmax";
	  string lName    = CoreStatePath (ti, bi) + "l";
	  string kName    = CoreStatePath (ti, bi) + "k";
	  string rName    = CoreStatePath (ti, bi) + "r";
	  HDFAttribIO<int> h_atom(atom), h_l(l);
	  HDFAttribIO<double> h_rmax(rmax);
	  HDFAttribIO<Vector<double> > h_g(g);
	  HDFAttribIO<Vector<double> > h_r(r);
	  h_atom.read (H5FileID, atomName.c_str());
	  h_l.read    (H5FileID,    lName.c_str());
	  h_rmax.read (H5FileID, rMaxName.c_str());
	  h_g.read    (H5FileID,   gName.c_str());
	  h_r.read    (H5FileID,   rName.c_str());
	  
	  fprintf (stderr, "  Core state:     ti=%3d  bi=%3d energy=%8.5f "
		   "k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
		   ti, bi, e, k[0], k[1], k[2], myComm->rank());
	}
	myComm->bcast (atom);  myComm->bcast(rmax);	
	myComm->bcast (l);     myComm->bcast (k);
	int ng = g.size();     myComm->bcast(ng);
	if (g.size() != ng) {
	  g.resize(ng);
	  r.resize(ng);
	}
	myComm->bcast (g);
	myComm->bcast (r);

	double Z = (double)IonTypes(atom);
	OrbitalSet->MuffinTins[atom].addCore (l, m, r, g, k, Z);
	icore++;
      }
      else {
	if (root) {
	  int ti   = SortBands[iorb].TwistIndex;
	  int bi   = SortBands[iorb].BandIndex;
	  double e = SortBands[iorb].Energy;
	  PosType twist, k;
	  twist = TwistAngles[ti];
	  k = orbitalSet->PrimLattice.k_cart(twist);
	  fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
		   ti, bi, e, k[0], k[1], k[2], myComm->rank());
	  
	  string vectorName = OrbitalPath (ti, bi) + "eigenvector";
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
	}
	myComm->bcast(splineData);
	set_multi_UBspline_3d_z 
	  (orbitalSet->MultiSpline, ival, splineData.data());
	
	// Now read muffin tin data
	for (int tin=0; tin<NumMuffinTins; tin++) {
	  // app_log() << "Reading data for muffin tin " << tin << endl;
	  PosType twist, k;
	  int lmax = MT_APW_lmax[tin];
	  int numYlm = (lmax+1)*(lmax+1);
	  Array<complex<double>,2> 
	    u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
	  Array<complex<double>,1> du_lm_dr (numYlm);
	  if (root) {
	    int ti   = SortBands[iorb].TwistIndex;
	    int bi   = SortBands[iorb].BandIndex;
	    twist = TwistAngles[ti];
	    k = orbitalSet->PrimLattice.k_cart(twist);
	    string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
	    string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
	    HDFAttribIO<Array<complex<double>,2> > h_u_lm_r(u_lm_r);
	    HDFAttribIO<Array<complex<double>,1> > h_du_lm_dr(du_lm_dr);
	    h_u_lm_r.read(H5FileID, uName.c_str());
	    h_du_lm_dr.read(H5FileID, duName.c_str());
	  }
	  myComm->bcast(u_lm_r);
	  myComm->bcast(du_lm_dr);
	  myComm->bcast(k);
	  double Z = (double)IonTypes(tin);
	  OrbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
	}
	ival++;
      } // valence state
      iorb++;
    }
    ExtendedMap_z[set] = orbitalSet->MultiSpline;
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
