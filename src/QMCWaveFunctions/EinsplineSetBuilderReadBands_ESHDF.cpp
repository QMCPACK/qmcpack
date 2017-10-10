//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "Message/CommOperators.h"
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>

namespace qmcplusplus
{

#if 0
void EinsplineSetBuilder::ReadBands_ESHDF(int spin, EinsplineSetExtended<std::complex<double > >* orbitalSet)
{
  update_token(__FILE__,__LINE__,"ReadBands_ESHDF:complex");
  ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF(EinsplineSetExtended<std::complex<double > >*");
  Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
  double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;
  c_prep.restart();
  bool root = myComm->rank()==0;
  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  // bcast other stuff
  myComm->bcast (NumDistinctOrbitals);
  myComm->bcast (NumValenceOrbs);
  myComm->bcast (NumCoreOrbs);
  int N = NumDistinctOrbitals;
  orbitalSet->kPoints.resize(N);
  orbitalSet->MakeTwoCopies.resize(N);
  orbitalSet->StorageValueVector.resize(N);
  orbitalSet->BlendValueVector.resize(N);
  orbitalSet->StorageLaplVector.resize(N);
  orbitalSet->BlendLaplVector.resize(N);
  orbitalSet->StorageGradVector.resize(N);
  orbitalSet->BlendGradVector.resize(N);
  orbitalSet->StorageHessVector.resize(N);
  orbitalSet->StorageGradHessVector.resize(N);
  orbitalSet->phase.resize(N);
  orbitalSet->eikr.resize(N);
  orbitalSet->NumValenceOrbs = NumValenceOrbs;
  orbitalSet->NumCoreOrbs    = NumCoreOrbs;
  // Read in k-points
  int numOrbs = orbitalSet->getOrbitalSetSize();
  int num = 0;
  if (root)
  {
    char s[1024];
    std::cerr << "#  Band    State   TwistIndex BandIndex Energy      Kx      Ky      Kz      K1      K2      K3    KmK " << std::endl;
    for (int iorb=0; iorb<N; iorb++)
    {
      int ti = SortBands[iorb].TwistIndex;
      int bi   = SortBands[iorb].BandIndex;
      double e = SortBands[iorb].Energy;
      int nd = (SortBands[iorb].MakeTwoCopies)?2:1;
      PosType twist  = TwistAngles[ti];
      orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
      orbitalSet->MakeTwoCopies[iorb] =
        (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
      sprintf (s, "%8d %8d %8d %8d %12.6f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6d\n",
          iorb, num, ti, bi, e, orbitalSet->kPoints[iorb][0], orbitalSet->kPoints[iorb][1], orbitalSet->kPoints[iorb][2],
          TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2],nd);
      std::cerr << s;
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
  nx=MeshSize[0];
  ny=MeshSize[1];
  nz=MeshSize[2];
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = PERIODIC;
  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;
  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;
  zBC.rCode = PERIODIC;
  x_grid.start = 0.0;
  x_grid.end = 1.0;
  x_grid.num = nx;
  y_grid.start = 0.0;
  y_grid.end = 1.0;
  y_grid.num = ny;
  z_grid.start = 0.0;
  z_grid.end = 1.0;
  z_grid.num = nz;
  // Create the multiUBspline object
  orbitalSet->MultiSpline =
    create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  //////////////////////////////////////
  // Create the MuffinTin APW splines //
  //////////////////////////////////////
  orbitalSet->MuffinTins.resize(NumMuffinTins);
  for (int tin=0; tin<NumMuffinTins; tin++)
  {
    orbitalSet->MuffinTins[tin].Atom = tin;
    orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
    orbitalSet->MuffinTins[tin].set_lattice(Lattice);
    orbitalSet->MuffinTins[tin].init_APW
    (MT_APW_rgrids[tin], MT_APW_lmax[tin],
     NumValenceOrbs);
  }
  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
  {
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
  bool isCore = bcastSortBands(spin,N,root);
  if(isCore)
  {
    APP_ABORT("Core states not supported by ES-HDF yet.");
  }
  t_prep += c_prep.elapsed();
  /** For valence orbitals,
   * - extended orbitals either in G or in R
   * - localized orbitals
   */
  //this can potentially break
  Array<ComplexType,3> splineData(nx,ny,nz);
  if(havePsig)//perform FFT using FFTW
  {
    c_init.restart();
    Array<ComplexType,3> FFTbox;
    FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
    fftw_plan FFTplan = fftw_plan_dft_3d
                        (MeshSize[0], MeshSize[1], MeshSize[2],
                         reinterpret_cast<fftw_complex*>(FFTbox.data()),
                         reinterpret_cast<fftw_complex*>(FFTbox.data()),
                         +1, FFTW_ESTIMATE);
    Vector<std::complex<double> > cG(MaxNumGvecs);
    //this will be parallelized with OpenMP
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      //Vector<std::complex<double> > cG;
      int ncg=0;
      int ti=SortBands[iorb].TwistIndex;
      c_h5.restart();
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
        HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
        h_cG.read (H5FileID, path.str().c_str());
        ncg=cG.size();
      }
      myComm->bcast(ncg);
      if(ncg != Gvecs[0].size())
      {
        APP_ABORT("Failed : ncg != Gvecs[0].size()");
      }
      if(!root)
        cG.resize(ncg);
      myComm->bcast(cG);
      t_h5 += c_h5.elapsed();
      c_unpack.restart();
      unpack4fftw(cG,Gvecs[0],MeshSize,FFTbox);
      t_unpack+= c_unpack.elapsed();
      c_fft.restart();
      fftw_execute (FFTplan);
      t_fft+= c_fft.elapsed();
      c_phase.restart();
      fix_phase_rotate_c2c(FFTbox,splineData,TwistAngles[ti]);
      t_phase+= c_phase.elapsed();
      c_spline.restart();
      set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
      t_spline+= c_spline.elapsed();
    }
    fftw_destroy_plan(FFTplan);
    t_init+=c_init.elapsed();
  }
  else
  {
    //this will be parallelized with OpenMP
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      //check dimension
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(splineData);
        h_splineData.read(H5FileID, path.str().c_str());
      }
      myComm->bcast(splineData);
      set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
    }
    //return true;
  }
  app_log() << "    READBANDS::PREP   = " << t_prep << std::endl;
  app_log() << "    READBANDS::H5     = " << t_h5 << std::endl;
  app_log() << "    READBANDS::UNPACK = " << t_unpack << std::endl;
  app_log() << "    READBANDS::FFT    = " << t_fft << std::endl;
  app_log() << "    READBANDS::PHASE  = " << t_phase << std::endl;
  app_log() << "    READBANDS::SPLINE = " << t_spline << std::endl;
  app_log() << "    READBANDS::SUM    = " << t_init << std::endl;
  //now localized orbitals
  for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
  {
    PosType twist=TwistAngles[SortBands[iorb].TwistIndex];
    // Read atomic orbital information
    for (int iat=0; iat<AtomicOrbitals.size(); iat++)
    {
      app_log() << "Reading orbital " << iat << " for band " << ival << std::endl;
      AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
      Array<std::complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm),
            poly_coefs(orb.PolyOrder+1,orb.Numlm);
      if (root)
      {
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
        std::ostringstream spline_path, poly_path;
        spline_path << path.str() << "radial_spline_" << iat;
        poly_path   << path.str() << "poly_coefs_"    << iat;
        HDFAttribIO<Array<std::complex<double>,2> > h_radial_spline(radial_spline);
        HDFAttribIO<Array<std::complex<double>,2> > h_poly_coefs(poly_coefs);
        h_radial_spline.read(H5FileID, spline_path.str().c_str());
        h_poly_coefs.read   (H5FileID, poly_path.str().c_str());
        // std::cerr << "radial_spline.size = (" << radial_spline.size(0)
        // 	 << ", " << radial_spline.size(1) << ")\n";
        // std::cerr << "poly_coefs.size = (" << poly_coefs.size(0)
        // 	 << ", " << poly_coefs.size(1) << ")\n";
      }
      myComm->bcast(radial_spline);
      myComm->bcast(poly_coefs);
      AtomicOrbitals[iat].set_band (ival, radial_spline, poly_coefs, twist);
    }
    // Now read muffin tin data
    for (int tin=0; tin<NumMuffinTins; tin++)
    {
      // app_log() << "Reading data for muffin tin " << tin << std::endl;
      PosType twist, k;
      int lmax = MT_APW_lmax[tin];
      int numYlm = (lmax+1)*(lmax+1);
      Array<std::complex<double>,2>
      u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
      Array<std::complex<double>,1> du_lm_dr (numYlm);
      if (root)
      {
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        std::string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
        std::string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
        HDFAttribIO<Array<std::complex<double>,2> > h_u_lm_r(u_lm_r);
        HDFAttribIO<Array<std::complex<double>,1> > h_du_lm_dr(du_lm_dr);
        h_u_lm_r.read(H5FileID, uName.c_str());
        h_du_lm_dr.read(H5FileID, duName.c_str());
      }
      myComm->bcast(u_lm_r);
      myComm->bcast(du_lm_dr);
      myComm->bcast(k);
      double Z = (double)IonTypes(tin);
      orbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
    }
  }
  orbitalSet->AtomicOrbitals = AtomicOrbitals;
  for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
    orbitalSet->AtomicOrbitals[i].registerTimers();
  //ExtendedMap_z[set] = orbitalSet->MultiSpline;
}



void EinsplineSetBuilder::ReadBands_ESHDF(int spin, EinsplineSetExtended<double>* orbitalSet)
{
  update_token(__FILE__,__LINE__,"ReadBands_ESHDF:double");
  ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF(EinsplineSetExtended<double>*");
  std::vector<AtomicOrbital<double> > realOrbs(AtomicOrbitals.size());
  for (int iat=0; iat<realOrbs.size(); iat++)
  {
    AtomicOrbital<std::complex<double> > &corb (AtomicOrbitals[iat]);
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
  orbitalSet->StorageValueVector.resize(N);
  orbitalSet->BlendValueVector.resize(N);
  orbitalSet->StorageLaplVector.resize(N);
  orbitalSet->BlendLaplVector.resize(N);
  orbitalSet->StorageGradVector.resize(N);
  orbitalSet->BlendGradVector.resize(N);
  orbitalSet->StorageHessVector.resize(N);
  orbitalSet->StorageGradHessVector.resize(N);
  orbitalSet->phase.resize(N);
  orbitalSet->eikr.resize(N);
  orbitalSet->NumValenceOrbs = NumValenceOrbs;
  orbitalSet->NumCoreOrbs    = NumCoreOrbs;
  orbitalSet->FirstOrderSplines.resize(IonPos.size());
  // Read in k-points
  int numOrbs = orbitalSet->getOrbitalSetSize();
  int num = 0;
  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  if (root)
  {
    for (int iorb=0; iorb<N; iorb++)
    {
      int ti = SortBands[iorb].TwistIndex;
      PosType twist  = TwistAngles[ti];
      orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
      orbitalSet->MakeTwoCopies[iorb] =
        (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
      num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
    }
    PosType twist0 = TwistAngles[SortBands[0].TwistIndex];
    for (int i=0; i<OHMMS_DIM; i++)
      if (std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)
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
  //int nx, ny, nz, bi, ti;
  int nx, ny, nz;
  nx=MeshSize[0];
  ny=MeshSize[1];
  nz=MeshSize[2];
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  if (orbitalSet->HalfG[0])
  {
    xBC.lCode = ANTIPERIODIC;
    xBC.rCode = ANTIPERIODIC;
  }
  else
  {
    xBC.lCode = PERIODIC;
    xBC.rCode = PERIODIC;
  }
  if (orbitalSet->HalfG[1])
  {
    yBC.lCode = ANTIPERIODIC;
    yBC.rCode = ANTIPERIODIC;
  }
  else
  {
    yBC.lCode = PERIODIC;
    yBC.rCode = PERIODIC;
  }
  if (orbitalSet->HalfG[2])
  {
    zBC.lCode = ANTIPERIODIC;
    zBC.rCode = ANTIPERIODIC;
  }
  else
  {
    zBC.lCode = PERIODIC;
    zBC.rCode = PERIODIC;
  }
  x_grid.start = 0.0;
  x_grid.end = 1.0;
  x_grid.num = nx;
  y_grid.start = 0.0;
  y_grid.end = 1.0;
  y_grid.num = ny;
  z_grid.start = 0.0;
  z_grid.end = 1.0;
  z_grid.num = nz;
  // Create the multiUBspline object
  orbitalSet->MultiSpline =
    create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  if (HaveOrbDerivs)
  {
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
  for (int tin=0; tin<NumMuffinTins; tin++)
  {
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
  if (root)
  {
    HDFAttribIO<int> h_isComplex(isComplex);
    h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");
  }
  myComm->bcast(isComplex);
  bool isCore = bcastSortBands(spin,N,root);
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
      app_log() << "   Reading complex psi_r and convert to real" << std::endl;
      Array<std::complex<double>,3> rawData;
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        int ti=SortBands[iorb].TwistIndex;
        if(root)
        {
          std::ostringstream path;
          path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
               << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
          HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(rawData);
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
      app_log() << "   Reading real psi_r" << std::endl;
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        if(root)
        {
          std::ostringstream path;
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
    Array<ComplexType,3> FFTbox;
    FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
    fftw_plan FFTplan = fftw_plan_dft_3d
                        (MeshSize[0], MeshSize[1], MeshSize[2],
                         reinterpret_cast<fftw_complex*>(FFTbox.data()),
                         reinterpret_cast<fftw_complex*>(FFTbox.data()),
                         +1, FFTW_ESTIMATE);
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      Vector<std::complex<double> > cG;
      int ncg=0;
      int ti=SortBands[iorb].TwistIndex;
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
        HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
        h_cG.read (H5FileID, path.str().c_str());
        ncg=cG.size();
      }
      myComm->bcast(ncg);
      if(ncg != Gvecs[0].size())
      {
        APP_ABORT("Failed : ncg != Gvecs[0].size()");
      }
      if(!root)
        cG.resize(ncg);
      myComm->bcast(cG);
      unpack4fftw(cG,Gvecs[0],MeshSize,FFTbox);
      fftw_execute (FFTplan);
      double phase_r,phase_i;
      compute_phase(FFTbox,TwistAngles[ti],phase_r,phase_i);
      fix_phase_rotate_c2r(FFTbox,splineData,TwistAngles[ti],phase_r,phase_i);
      set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
    }
    fftw_destroy_plan(FFTplan);
  }
  for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
  {
    // Read atomic orbital information
    for (int iat=0; iat<realOrbs.size(); iat++)
    {
      app_log() << "Reading orbital " << iat << " for band " << ival << std::endl;
      AtomicOrbital<double> &orb = realOrbs[iat];
      //AtomicOrbital<std::complex<double> > &orb = realOrbs[iat];
      Array<std::complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm),
            poly_coefs(orb.PolyOrder+1,orb.Numlm);
      int ti   = SortBands[iorb].TwistIndex;
      if (root)
      {
        int bi   = SortBands[iorb].BandIndex;
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        std::ostringstream spline_path, poly_path;
        spline_path << path.str() << "radial_spline_" << iat;
        poly_path   << path.str() << "poly_coefs_"    << iat;
        HDFAttribIO<Array<std::complex<double>,2> > h_radial_spline(radial_spline);
        HDFAttribIO<Array<std::complex<double>,2> > h_poly_coefs(poly_coefs);
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
    for (int tin=0; tin<NumMuffinTins; tin++)
    {
      // app_log() << "Reading data for muffin tin " << tin << std::endl;
      PosType twist, k;
      int lmax = MT_APW_lmax[tin];
      int numYlm = (lmax+1)*(lmax+1);
      Array<std::complex<double>,2>
      u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
      Array<std::complex<double>,1> du_lm_dr (numYlm);
      int ti   = SortBands[iorb].TwistIndex;
      if (root)
      {
        int bi   = SortBands[iorb].BandIndex;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        std::string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
        std::string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
        HDFAttribIO<Array<std::complex<double>,2> > h_u_lm_r(u_lm_r);
        HDFAttribIO<Array<std::complex<double>,1> > h_du_lm_dr(du_lm_dr);
        h_u_lm_r.read(H5FileID, uName.c_str());
        h_du_lm_dr.read(H5FileID, duName.c_str());
      }
      myComm->bcast(u_lm_r);
      myComm->bcast(du_lm_dr);
      myComm->bcast(k);
      double Z = (double)IonTypes(tin);
      orbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
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
//			  << bi << " kpoint " << ti << std::endl;
//		ostringstream path;
//		path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/"
//		     << "dpsi_" << ion << "_" << dim << "_r";
//		string psirName = path.str();
//		if (isComplex) {
//		  HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
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
  //ExtendedMap_d[set] = orbitalSet->MultiSpline;
}
#endif

bool EinsplineSetBuilder::ReadGvectors_ESHDF()
{
  update_token(__FILE__,__LINE__,"ReadGvectors_ESHDF");
  bool root=myComm->rank() ==0;
  //this is always ugly
  MeshSize = 0;
  int hasPsig=1;
#if defined(__bgp__)||(__bgq__)
  if(root)
  {
    hid_t gid=H5Dopen(H5FileID,"/electrons/kpoint_0/spin_0/state_0/psi_g");
    if(gid<0)
      hasPsig=0;
    H5Dclose(gid);
  }
  myComm->bcast(hasPsig);
#else
  if(root)
  {
    HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize);
    h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
    h_mesh.read (H5FileID, "/electrons/mesh");
  }
  myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
#endif
  if(hasPsig)
  {
    int nallowed=257;
    int allowed[] =
    {
      72,75,80,81,90,96,100,108,120,125,
      128,135,144,150,160,162,180,192,200,216,
      225,240,243,250,256,270,288,300,320,324,
      360,375,384,400,405,432,450,480,486,500,
      512,540,576,600,625,640,648,675,720,729,
      750,768,800,810,864,900,960,972,1000,1024,
      1080,1125,1152,1200,1215,1250,1280,1296,1350,1440,
      1458,1500,1536,1600,1620,1728,1800,1875,1920,1944,
      2000,2025,2048,2160,2187,2250,2304,2400,2430,2500,
      2560,2592,2700,2880,2916,3000,3072,3125,3200,3240,
      3375,3456,3600,3645,3750,3840,3888,4000,4050,4096,
      4320,4374,4500,4608,4800,4860,5000,5120,5184,5400,
      5625,5760,5832,6000,6075,6144,6250,6400,6480,6561,
      6750,6912,7200,7290,7500,7680,7776,8000,8100,8192,
      8640,8748,9000,9216,9375,9600,9720,10000,10125,10240,
      10368,10800,10935,11250,11520,11664,12000,12150,12288,12500,
      12800,12960,13122,13500,13824,14400,14580,15000,15360,15552,
      15625,16000,16200,16384,16875,17280,17496,18000,18225,18432,
      18750,19200,19440,19683,20000,20250,20480,20736,21600,21870,
      22500,23040,23328,24000,24300,24576,25000,25600,25920,26244,
      27000,27648,28125,28800,29160,30000,30375,30720,31104,31250,
      32000,32400,32768,32805,33750,34560,34992,36000,36450,36864,
      37500,38400,38880,39366,40000,40500,40960,41472,43200,43740,
      45000,46080,46656,46875,48000,48600,49152,50000,50625,51200,
      51840,52488,54000,54675,55296,56250,57600,58320,59049,60000,
      60750,61440,62208,62500,64000,64800,65536
    };
    int numk=0;
    MaxNumGvecs=0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int,3> maxIndex(0,0,0);
    Gvecs.resize(NumTwists);
    {
      int numg=0;
      if(root)
      {
        std::ostringstream Gpath;
        Gpath    << "/electrons/kpoint_0/gvectors";
        HDFAttribIO<std::vector<TinyVector<int,3> > > h_Gvecs(Gvecs[0]);
        h_Gvecs.read (H5FileID, Gpath.str().c_str());
        numg=Gvecs[0].size();
      }
      myComm->bcast(numg);
      if(!root)
        Gvecs[0].resize(numg);
      myComm->bcast(Gvecs[0]);
      MaxNumGvecs=Gvecs[0].size();
      for (int ig=0; ig<Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } //done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0*MeshFactor*maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0*MeshFactor*maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0*MeshFactor*maxIndex[2]);
    //only use 2^a 3^b 5^c where a>=2  up to 65536
    int *ix=std::lower_bound(allowed,allowed+nallowed,MeshSize[0]);
    int *iy=std::lower_bound(allowed,allowed+nallowed,MeshSize[1]);
    int *iz=std::lower_bound(allowed,allowed+nallowed,MeshSize[2]);
    MeshSize[0]=(MeshSize[0]>128)? *ix:(MeshSize[0]+MeshSize[0]%2);
    MeshSize[1]=(MeshSize[1]>128)? *iy:(MeshSize[1]+MeshSize[1]%2);
    MeshSize[2]=(MeshSize[2]>128)? *iz:(MeshSize[2]+MeshSize[2]%2);
    if(Version[0]<2)
    {
      //get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  ESHDF::Version " << Version << std::endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. Regenerate orbital files using updated QE." << std::endl;
      for(int k=0; k<DistinctTwists.size(); ++k)
      {
        int ik=DistinctTwists[k];
        if(ik==0)
          continue; //already done
        int numg=0;
        if(root)
        {
          std::ostringstream Gpath;
          Gpath    << "/electrons/kpoint_" << ik << "/gvectors";
          HDFAttribIO<std::vector<TinyVector<int,3> > > h_Gvecs(Gvecs[ik]);
          h_Gvecs.read (H5FileID, Gpath.str().c_str());
          numg=Gvecs[ik].size();
        }
        myComm->bcast(numg);
        if(numg==0)
        {
          //copy kpoint_0, default
          Gvecs[ik]=Gvecs[0];
        }
        else
        {
          if(numg !=  MaxNumGvecs)
          {
            std::ostringstream o;
            o<< "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
             << " This is not supported anymore. Rerun pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if(!root)
            Gvecs[ik].resize(numg);
          myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << std::endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << std::endl;
  app_log().flush();
  return hasPsig;
}


}

