//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
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

  std::map<H5OrbSet,SPOSetBase*,H5OrbSet>  EinsplineSetBuilder::SPOSetMap;
  std::map<TinyVector<int,4>,EinsplineSetBuilder::OrbType*,Int4less> EinsplineSetBuilder::OrbitalMap;
  std::map<H5OrbSet,multi_UBspline_3d_z*,H5OrbSet> EinsplineSetBuilder::ExtendedMap_z;
  std::map<H5OrbSet,multi_UBspline_3d_d*,H5OrbSet> EinsplineSetBuilder::ExtendedMap_d;

  EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, 
      PtclPoolType& psets, xmlNodePtr cur) 
    : XMLRoot(cur), TileFactor(1,1,1), TwistNum(0), LastSpinSet(-1), 
      NumOrbitalsRead(-1), NumMuffinTins(0), NumCoreStates(0),
      NumBands(0), NumElectrons(0), NumSpins(0), NumTwists(0),
      ParticleSets(psets), TargetPtcl(p), H5FileID(-1),
      Format(QMCPACK), makeRotations(false), MeshFactor(1.0),
      MeshSize(0,0,0)
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
  EinsplineSetBuilder::CheckLattice()
  {
    bool match=true;
    for (int i=0; i<OHMMS_DIM; i++)
      for (int j=0; j<OHMMS_DIM; j++) {
        RealType diff = SuperLattice(i,j) - TargetPtcl.Lattice.R(i,j);
        match = match && (std::fabs(diff) < 1.0e-6);
      }

    if(!match) 
    {
      ostringstream o;
      o.setf(std::ios::scientific, std::ios::floatfield);
      o.precision(6);
      o << "EinsplineSetBuilder::ReadOrbitalInfo_ESHDF \n"
        << "Mismatched supercell lattices.\n";
      o << " Lattice in ESHDF5 " << endl;
      o << SuperLattice << endl;
      o << " Lattice in xml" << endl;
      o << TargetPtcl.Lattice.R << endl;
      o << " Difference " << endl;
      o << SuperLattice-TargetPtcl.Lattice.R << endl;
      APP_ABORT(o.str());
    }
    return match;
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
    bool useGPU = false;
    attribs.add (H5FileName, "href");
    attribs.add (TileFactor, "tile");
    attribs.add (sortBands,  "sort");
    attribs.add (TileMatrix, "tilematrix");
    attribs.add (TwistNum,   "twistnum");
    attribs.add (sourceName, "source");
    attribs.add (MeshFactor, "meshfactor");
    attribs.add (useGPU,     "gpu");    
    attribs.put (XMLRoot);
    attribs.add (numOrbs,    "size");
    attribs.add (numOrbs,    "norbs");
    attribs.put (cur);

    ///////////////////////////////////////////////
    // Read occupation information from XML file //
    ///////////////////////////////////////////////
    cur = cur->children;
    int spinSet = -1;
    vector<int> Occ_Old(0,0);
    Occ.resize(0,0);
    bool NewOcc(false);
    while (cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "occupation") {
	string occ_mode("ground");
  occ_format="energy";
   particle_hole_pairs=0;
        OhmmsAttributeSet oAttrib;
        oAttrib.add(occ_mode,"mode");
        oAttrib.add(spinSet,"spindataset");
        oAttrib.add(occ_format,"format");
        oAttrib.add(particle_hole_pairs,"pairs");
        oAttrib.put(cur);
	if(occ_mode == "excited"){
          putContent(Occ,cur);
        } else if(occ_mode != "ground") {
	  app_error() << "Only ground state occupation currently supported "
		      << "in EinsplineSetBuilder.\n";
          APP_ABORT("EinsplineSetBuilder::createSPOSet");
	}
      }
      cur = cur->next;
    }
    if (Occ != Occ_Old){
      NewOcc=true;
      Occ_Old = Occ;
    }
    else NewOcc=false;


    H5OrbSet aset(H5FileName, spinSet, numOrbs);
    std::map<H5OrbSet,SPOSetBase*,H5OrbSet>::iterator iter;
    iter = SPOSetMap.find (aset);
    if ((iter != SPOSetMap.end() ) && (!NewOcc)) {
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
#ifdef QMC_CUDA
    else if (AtomicOrbitals.size() > 0) {
      if (UseRealOrbitals) 
	OrbitalSet = new EinsplineSetHybrid<double>;
      else
	OrbitalSet = new EinsplineSetHybrid<complex<double> >;
    }
#endif
    else {
      if (UseRealOrbitals) 
	OrbitalSet = new EinsplineSetExtended<double>;
      else
	OrbitalSet = new EinsplineSetExtended<complex<double> >;
    }

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
    OrbitalSet->GGt=dot(transpose(OrbitalSet->PrimLattice.G),
			OrbitalSet->PrimLattice.G);
    app_log() << "GGt = \n"
	      << OrbitalSet->GGt << endl;
    OrbitalSet->setOrbitalSetSize (numOrbs);
    OrbitalSet->BasisSetSize   = numOrbs;
    TileIons();
    
    if (HaveLocalizedOrbs) {
      EinsplineSetLocal *restrict orbitalSet = 
	dynamic_cast<EinsplineSetLocal*>(OrbitalSet);

#pragma omp critical(read_einspline_orbs)
      {
	if ((spinSet == LastSpinSet) && (numOrbs <= NumOrbitalsRead) && (!NewOcc) )
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
#ifdef QMC_CUDA
//	  if (true || useGPU) {
// 	    app_log() << "Copying einspline orbitals to GPU.\n";
// 	    create_multi_UBspline_3d_cuda 
// 	      (orbitalSet->MultiSpline, orbitalSet->CudaMultiSpline);
// 	    app_log() << "Successful copy.\n";
// 	    // Destroy original CPU spline
// 	    // HACK HACK HACK
// 	    //destroy_Bspline (orbitalSet->MultiSpline);
// 	    gpu::host_vector<CudaRealType> L_host(9), Linv_host(9);
// 	    orbitalSet->Linv_cuda.resize(9);
// 	    orbitalSet->L_cuda.resize(9);
// 	    for (int i=0; i<3; i++)
// 	      for (int j=0; j<3; j++) {
// 		L_host[i*3+j]    = (float)orbitalSet->PrimLattice.R(i,j);
// 		Linv_host[i*3+j] = (float)orbitalSet->PrimLattice.G(i,j);
// 	      }
// 	    orbitalSet->L_cuda    = L_host;
// 	    orbitalSet->Linv_cuda = Linv_host;
//	  }
#endif
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
#ifdef QMC_CUDA
// 	  if (useGPU) {
// 	    app_log() << "Copying einspline orbitals to GPU.\n";
// 	    create_multi_UBspline_3d_cuda (orbitalSet->MultiSpline,
// 	    				   orbitalSet->CudaMultiSpline);
// 	    app_log() << "Successful copy.\n";
// 	    // Destroy original CPU spline
// 	    // HACK HACK HACK
// 	    //destroy_Bspline (orbitalSet->MultiSpline);

// 	    gpu::host_vector<CudaRealType> L_host(9), Linv_host(9);
// 	    orbitalSet->Linv_cuda.resize(9);
// 	    orbitalSet->L_cuda.resize(9);
// 	    for (int i=0; i<3; i++)
// 	      for (int j=0; j<3; j++) {
// 		L_host[i*3+j]    = (float)orbitalSet->PrimLattice.R(i,j);
// 		Linv_host[i*3+j] = (float)orbitalSet->PrimLattice.G(i,j);
// 	      }
// 	    orbitalSet->L_cuda    = L_host;
// 	    orbitalSet->Linv_cuda = Linv_host;
// 	  }
#endif
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

    SPOSetMap[aset] = OrbitalSet;
    
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
#ifdef QMC_CUDA
    if (useGPU) {
      app_log() << "Initializing GPU data structures.\n";
      OrbitalSet->initGPU();
    }
#endif
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
    int numCopies = abs(det(TileMatrix));
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
	    if ((uSuper[0] >= -1.0e-4) && (uSuper[0] < 0.9999) &&
		(uSuper[1] >= -1.0e-4) && (uSuper[1] < 0.9999) &&
		(uSuper[2] >= -1.0e-4) && (uSuper[2] < 0.9999)) {
	      IonPos[index]= r;
	      IonTypes[index]= primTypes[iat];
	      index++;
	    }
	  }
    if (index != primPos.size()*numCopies) {
      
      app_error() << "The number of tiled ions, " << IonPos.size() 
		  << ", does not match the expected number of "
		  << primPos.size()*numCopies << " or the index "<< index <<".  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::TileIons()");
    }
    if (myComm->rank() == 0) {
      fprintf (stderr, "Supercell reduced ion positions = \n");
      for (int i=0; i<IonPos.size(); i++) {
	PosType u = SuperCell.toUnit(IonPos[i]);
	fprintf (stderr, "   %14.10f %14.10f %14.10f\n",
		 u[0], u[1], u[2]);
	//		 IonPos[i][0], IonPos[i][1], IonPos[i][2]);
      }
    }
  }

  
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
    TargetPtcl.setTwist(superFracs[TwistNum]);
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
    int numTwistsNeeded = abs(det(TileMatrix));
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
    UseRealOrbitals = false;
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


