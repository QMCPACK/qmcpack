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
/** @file common.cpp
 *
 * Instantiates the static data
 * Implements member functions of EinsplineSetBuilder
 * - EinsplineSetBuilder
 * -
*/
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
#include "Numerics/HDFSTLAttrib.h"
#include "ParticleIO/ESHDFParticleParser.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus {

  SPOSetBase*
  EinsplineSetBuilder::createSPOSet(xmlNodePtr cur) {

    OhmmsAttributeSet attribs;
    int numOrbs = 0;
    qafm=0;
    int sortBands(1);
    string sourceName;
#if defined(QMC_CUDA)
    string useGPU="yes";
#else
    string useGPU="no";
#endif
    attribs.add (H5FileName, "href");
    attribs.add (TileFactor, "tile");
    attribs.add (sortBands,  "sort");
    attribs.add (qafm,  "afmshift");
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
    if ((iter != SPOSetMap.end() ) && (!NewOcc) && (qafm==0)) {
      qafm=0;
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

    Timer mytimer;
    mytimer.restart();

    /////////////////////////////////////////////////////////////////
    // Read the basic orbital information, without reading all the //
    // orbitals themselves.                                        //
    /////////////////////////////////////////////////////////////////
    if (myComm->rank() == 0) 
      if (!ReadOrbitalInfo()) {
	app_error() << "Error reading orbital info from HDF5 file.  Aborting.\n";
        APP_ABORT("EinsplineSetBuilder::createSPOSet");
      }
    app_log() <<  "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << mytimer.elapsed() << endl;
    myComm->barrier();

    mytimer.restart();
    BroadcastOrbitalInfo();
    app_log() <<  "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << mytimer.elapsed() << endl;

    app_log().flush();

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
    else 
    {
      mytimer.restart();
      if (UseRealOrbitals) 
      { 
	EinsplineSetExtended<double> *restrict orbitalSet =
	  dynamic_cast<EinsplineSetExtended<double>* > (OrbitalSet);    

        OccupyBands(spinSet, sortBands);

	{ 
	  if (Format == ESHDF)
	    ReadBands_ESHDF(spinSet,orbitalSet);
	  else
	    ReadBands(spinSet, orbitalSet); 

          app_log() <<  "TIMER  EinsplineSetBuilder::ReadBands " << mytimer.elapsed() << endl;

//#ifdef QMC_CUDA
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
//#endif
	}
      }
      else 
      {
	EinsplineSetExtended<complex<double> > *restrict orbitalSet = 
	  dynamic_cast<EinsplineSetExtended<complex<double> >*>(OrbitalSet);
	OccupyBands(spinSet, sortBands);
	{ 
	  if (Format == ESHDF)
	    ReadBands_ESHDF(spinSet,orbitalSet);
	  else
	    ReadBands(spinSet, orbitalSet); 
//#ifdef QMC_CUDA
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
//#endif
	}
      }
      app_log() <<  "TIMER  EinsplineSetBuilder::ReadBands " << mytimer.elapsed() << endl;
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
    if (useGPU == "yes" || useGPU == "1")
    {
      app_log() << "Initializing GPU data structures.\n";
      OrbitalSet->initGPU();
    }
#endif
    return OrbitalSet;
  }
}


