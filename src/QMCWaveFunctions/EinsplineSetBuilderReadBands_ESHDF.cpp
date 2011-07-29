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

#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "Message/CommOperators.h"

#include <QMCWaveFunctions/einspline_helper.hpp>

namespace qmcplusplus {

  void EinsplineSetBuilder::ReadBands_ESHDF(int spin, EinsplineSetExtended<complex<double > >* orbitalSet)
  {

    ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF(EinsplineSetExtended<complex<double > >*");

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
    
    ///check mesh or ready for FFT grid
    bool havePsig=ReadGvectors_ESHDF();
    app_log() << "MeshSize = (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";

    int nx, ny, nz, bi, ti;
    nx=MeshSize[0]; ny=MeshSize[1]; nz=MeshSize[2];
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

    int isComplex=1;
    if (root) 
    {
      HDFAttribIO<int> h_isComplex(isComplex);
      h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");
    }
    myComm->bcast(isComplex);
    if (!isComplex) 
    {
      APP_ABORT("Expected complex orbitals in ES-HDF file, but found real ones.");
    }

    EinsplineSetBuilder::RotateBands_ESHDF(spin, orbitalSet);

    bool isCore = bcastSortBands(N,root);
    if(isCore)
    {
      APP_ABORT("Core states not supported by ES-HDF yet.");
    }

    /** For valence orbitals,
     * - extended orbitals either in G or in R
     * - localized orbitals
     */

    //this can potentially break
    Array<complex<double>,3> splineData(nx,ny,nz);

    if(havePsig)//perform FFT using FFTW
    {
      FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
      FFTplan = fftw_plan_dft_3d 
        (MeshSize[0], MeshSize[1], MeshSize[2],
         reinterpret_cast<fftw_complex*>(FFTbox.data()),
         reinterpret_cast<fftw_complex*>(FFTbox.data()),
         +1, FFTW_ESTIMATE);

      //this will be parallelized with OpenMP
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        Vector<complex<double> > cG;
        int ncg=0;
        int ti=SortBands[iorb].TwistIndex;
        if(root)
        {
          ostringstream path;
          path << "/electrons/kpoint_" << ti    //SortBands[iorb].TwistIndex 
            << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
          HDFAttribIO<Vector<complex<double> > >  h_cG(cG);
          h_cG.read (H5FileID, path.str().c_str());
          ncg=cG.size();
        }
        myComm->bcast(ncg);
        if(ncg != Gvecs[ti].size())
        {
          APP_ABORT("Failed : ncg != Gvecs[ti].size()");
        }

        if(!root) cG.resize(ncg);
        myComm->bcast(cG);

        unpack4fftw(cG,Gvecs[ti],MeshSize,FFTbox);
        fftw_execute (FFTplan);
        fix_phase_rotate(FFTbox,splineData,TwistAngles[ti]);
        set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
      }

      fftw_destroy_plan(FFTplan);
    }
    else
    {
      //this will be parallelized with OpenMP
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        //check dimension
        if(root)
        {
          ostringstream path;
          path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex 
            << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
          HDFAttribIO<Array<complex<double>,3> >  h_splineData(splineData);
          h_splineData.read(H5FileID, path.str().c_str());
        }
        myComm->bcast(splineData);
        set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
      }

      //return true;
    }

    //now localized orbitals
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      PosType twist=TwistAngles[SortBands[iorb].TwistIndex];
      // Read atomic orbital information
      for (int iat=0; iat<AtomicOrbitals.size(); iat++) 
      {
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
    }

    orbitalSet->AtomicOrbitals = AtomicOrbitals;
    for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
      orbitalSet->AtomicOrbitals[i].registerTimers();

    ExtendedMap_z[set] = orbitalSet->MultiSpline;
  }



  void EinsplineSetBuilder::ReadBands_ESHDF(int spin, EinsplineSetExtended<double>* orbitalSet)
  {
    ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF(EinsplineSetExtended<double>*");

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
      EinsplineSetBuilder::RotateBands_ESHDF(spin, orbitalSet);
    }
    myComm->bcast(orbitalSet->kPoints);
    myComm->bcast(orbitalSet->MakeTwoCopies);
    myComm->bcast(orbitalSet->HalfG);

    // First, check to see if we have already read this in
    H5OrbSet set(H5FileName, spin, N);
    
    bool havePsir=!ReadGvectors_ESHDF();

    app_log() << "MeshSize = (" << MeshSize[0] << ", " 
	      << MeshSize[1] << ", " << MeshSize[2] << ")\n";

    int nx, ny, nz, bi, ti;
    nx=MeshSize[0]; ny=MeshSize[1]; nz=MeshSize[2];
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

    for (int iat=0; iat<realOrbs.size(); iat++) 
    {
      realOrbs[iat].set_num_bands(NumValenceOrbs);
      realOrbs[iat].allocate();
    }
           
    int isComplex;
    if (root) {
      HDFAttribIO<int> h_isComplex(isComplex);
      h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");
    }
    myComm->bcast(isComplex);

    bool isCore = bcastSortBands(N,root);
    if(isCore)
    {
      APP_ABORT("Core states not supported by ES-HDF yet.");
    }

    //this is common
    Array<double,3> splineData(nx,ny,nz);
    if(havePsir)
    {
      if(isComplex)
      {
        app_log() << "   Reading complex psi_r and convert to real" << endl;
        Array<complex<double>,3> rawData;
        for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
        {
          if(root)
          {
            ostringstream path;
            path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex 
              << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
            HDFAttribIO<Array<complex<double>,3> >  h_splineData(rawData);
            h_splineData.read(H5FileID, path.str().c_str());
          }
          myComm->bcast(rawData);
          //multiply twist factor and project on the real
          fix_phase_c2r(rawData,splineData,TwistAngles[ti]);
          set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
        }
      }
      else
      {
        app_log() << "   Reading real psi_r" << endl;
        for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
        {
          if(root)
          {
            ostringstream path;
            path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex 
              << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
            HDFAttribIO<Array<double,3> >  h_splineData(splineData);
            h_splineData.read(H5FileID, path.str().c_str());
          }
          myComm->bcast(splineData);
          set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
        }
      }
    }
    else
    {
      FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
      FFTplan = fftw_plan_dft_3d 
        (MeshSize[0], MeshSize[1], MeshSize[2],
         reinterpret_cast<fftw_complex*>(FFTbox.data()),
         reinterpret_cast<fftw_complex*>(FFTbox.data()),
         +1, FFTW_ESTIMATE);

      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        Vector<complex<double> > cG;
        int ncg=0;
        int ti=SortBands[iorb].TwistIndex;
        if(root)
        {
          ostringstream path;
          path << "/electrons/kpoint_" << ti    //SortBands[iorb].TwistIndex 
            << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
          cout << "FFT on " <<  path.str() << endl;
          HDFAttribIO<Vector<complex<double> > >  h_cG(cG);
          h_cG.read (H5FileID, path.str().c_str());
          ncg=cG.size();
        }
        myComm->bcast(ncg);
        if(ncg != Gvecs[ti].size())
        {
          APP_ABORT("Failed : ncg != Gvecs[ti].size()");
        }

        if(!root) cG.resize(ncg);
        myComm->bcast(cG);
        unpack4fftw(cG,Gvecs[ti],MeshSize,FFTbox);
        fftw_execute (FFTplan);
        fix_phase_rotate(FFTbox,splineData,TwistAngles[ti]);
	set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
      }

      fftw_destroy_plan(FFTplan);
    }

    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      // Read atomic orbital information
      for (int iat=0; iat<realOrbs.size(); iat++) {
        app_log() << "Reading orbital " << iat << " for band " << ival << endl;
        AtomicOrbital<double> &orb = realOrbs[iat];
        //AtomicOrbital<complex<double> > &orb = realOrbs[iat];
        Array<complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm), 
          poly_coefs(orb.PolyOrder+1,orb.Numlm);

        int ti   = SortBands[iorb].TwistIndex;
        if (root) { 
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
        realOrbs[iat].set_band (ival, radial_spline, poly_coefs, TwistAngles[ti]);
      }
    }

    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      // Now read muffin tin data
      for (int tin=0; tin<NumMuffinTins; tin++) {
        // app_log() << "Reading data for muffin tin " << tin << endl;
        PosType twist, k;
        int lmax = MT_APW_lmax[tin];
        int numYlm = (lmax+1)*(lmax+1);
        Array<complex<double>,2> 
          u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
        Array<complex<double>,1> du_lm_dr (numYlm);
        int ti   = SortBands[iorb].TwistIndex;
        if (root) {
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
    }

    //FIX HaveOrbDerivs after debugging
//	// Now read orbital derivatives if we have them
//	if (HaveOrbDerivs) {
//	  for (int ion=0; ion<IonPos.size(); ion++) 
//	    for (int dim=0; dim<OHMMS_DIM; dim++) {
//	      if (root) {
//		int ti   = SortBands[iorb].TwistIndex;
//		int bi   = SortBands[iorb].BandIndex;
//		
//		app_log() << "Reading orbital derivative for ion " << ion 
//			  << " dim " << dim << " spin " << spin << " band "
//			  << bi << " kpoint " << ti << endl;
//		ostringstream path;
//		path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/"
//		     << "dpsi_" << ion << "_" << dim << "_r";
//		string psirName = path.str();
//		if (isComplex) {
//		  HDFAttribIO<Array<complex<double>,3> > h_rawData(rawData);
//		  h_rawData.read(H5FileID, psirName.c_str());
//		  if ((rawData.size(0) != nx) ||
//		      (rawData.size(1) != ny) ||
//		      (rawData.size(2) != nz)) {
//		    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
//		    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
//		    abort();
//		  }
//#pragma omp parallel for
//		  for (int ix=0; ix<nx; ix++) {
//		  PosType ru;
//		    ru[0] = (RealType)ix / (RealType)nx;
//		    for (int iy=0; iy<ny; iy++) {
//		      ru[1] = (RealType)iy / (RealType)ny;
//		      for (int iz=0; iz<nz; iz++) {
//			ru[2] = (RealType)iz / (RealType)nz;
//			double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
//			double s, c;
//			sincos(phi, &s, &c);
//			complex<double> phase(c,s);
//			complex<double> z = phase*rawData(ix,iy,iz);
//			splineData(ix,iy,iz) = z.real();
//		      }
//		    }
//		  }
//		}
//		else {
//		  HDFAttribIO<Array<double,3> >  h_splineData(splineData);
//		  h_splineData.read(H5FileID, psirName.c_str());
//		  if ((splineData.size(0) != nx) ||
//		      (splineData.size(1) != ny) ||
//		      (splineData.size(2) != nz)) {
//		    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
//		    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
//		    abort();
//		  }
//		}
//	      }
//	      myComm->bcast(splineData);
//	      set_multi_UBspline_3d_d 
//		(orbitalSet->FirstOrderSplines[ion][dim], ival, splineData.data());
//	    }
//	}
//	
//
//	
    orbitalSet->AtomicOrbitals = realOrbs;
    for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
      orbitalSet->AtomicOrbitals[i].registerTimers();
    
    ExtendedMap_d[set] = orbitalSet->MultiSpline;
  }

}

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: EinsplineSetBuilderESHDF.cpp 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
