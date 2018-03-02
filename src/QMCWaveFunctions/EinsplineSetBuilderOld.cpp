//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

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

namespace qmcplusplus
{

bool
EinsplineSetBuilder::ReadOrbitalInfo()
{
  update_token(__FILE__,__LINE__,"ReadOrbitalInfo");
  // Handle failed file open gracefully by temporarily replacing error handler
  H5E_auto_t old_efunc;
  void *old_efunc_data;
  H5Eget_auto(&old_efunc, &old_efunc_data);
  H5Eset_auto(NULL, NULL);
  H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  //  H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  if (H5FileID < 0)
  {
    app_error() << "Could not open HDF5 file \"" << H5FileName
                << "\" in EinsplineSetBuilder::ReadOrbitalInfo.  Aborting.\n";
    APP_ABORT("EinsplineSetBuilder::ReadOrbitalInfo");
  }
  H5Eset_auto(old_efunc,old_efunc_data);
  
  // Read format
  std::string format;
  HDFAttribIO<std::string> h_format(format);
  h_format.read(H5FileID, "/format");
  HDFAttribIO<TinyVector<int,3> > h_Version(Version);
  h_Version.read (H5FileID, "/version");
  app_log() << "  HDF5 orbital file version "
            << Version[0] << "." << Version[1] << "." << Version[2] << "\n";
  if (format.find("ES")<format.size())
  {
    Format = ESHDF;
    return ReadOrbitalInfo_ESHDF();
  }
  //////////////////////////////////////////////////
  // Read basic parameters from the orbital file. //
  //////////////////////////////////////////////////
  // Check the version
  if (Version[0]==0 && Version[1]== 11)
  {
    parameterGroup  = "/parameters_0";
    ionsGroup       = "/ions_2";
    eigenstatesGroup = "/eigenstates_3";
  }
  else
    if (Version[0]==0 && Version[1]==20)
    {
      parameterGroup  = "/parameters";
      ionsGroup       = "/ions";
      eigenstatesGroup = "/eigenstates";
    }
    else
    {
      std::ostringstream o;
      o << "Unknown HDF5 orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << "\n";
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
  CheckLattice();
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
            << ", muffin tins=" << NumMuffinTins << std::endl;
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
  for (int tin=0; tin<NumMuffinTins; tin++)
  {
    std::ostringstream MTstream;
    if (NumMuffinTins > 1)
      MTstream << parameterGroup << "/muffin_tins/muffin_tin_" << tin;
    else
      MTstream << parameterGroup << "/muffin_tins/muffin_tin";
    std::string MTgroup = MTstream.str();
    HDFAttribIO<int> h_lmax(MT_APW_lmax[tin]),
                h_num_radial_points(MT_APW_num_radial_points[tin]);
    HDFAttribIO<double> h_radius (MT_APW_radii[tin]);
    HDFAttribIO<TinyVector<double, OHMMS_DIM> > h_center (MT_centers[tin]);
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
  for (int ti=0; ti<NumTwists; ti++)
  {
    std::ostringstream path;
    if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
      path << eigenstatesGroup << "/twist_" << ti << "/twist_angle";
    else
      path << eigenstatesGroup << "/twist/twist_angle";
    TinyVector<double, OHMMS_DIM> TwistAngles_DP;
    HDFAttribIO<TinyVector<double, OHMMS_DIM> > h_Twist(TwistAngles_DP);
    h_Twist.read (H5FileID, path.str().c_str());
    TwistAngles[ti] = TwistAngles_DP;
    snprintf (buff, 1000, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n",
              TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    app_log() << buff;
  }
  //////////////////////////////////////////////////////////
  // If the density has not been set in TargetPtcl, and   //
  // the density is available, read it in and save it     //
  // in TargetPtcl.                                       //
  //////////////////////////////////////////////////////////
  if (!TargetPtcl.Density_G.size())
  {
    HDFAttribIO<std::vector<TinyVector<int,OHMMS_DIM> > >
    h_reduced_gvecs(TargetPtcl.DensityReducedGvecs);
    Array<double, OHMMS_DIM> Density_r_DP;
    HDFAttribIO<Array<double, OHMMS_DIM> >  h_density_r (Density_r_DP);
    h_reduced_gvecs.read (H5FileID, "/density/reduced_gvecs");
    h_density_r.read (H5FileID,     "/density/rho_r");
    TargetPtcl.Density_r = Density_r_DP;
    int numG = TargetPtcl.DensityReducedGvecs.size();
    // Convert primitive G-vectors to supercell G-vectors
    for (int iG=0; iG < numG; iG++)
      TargetPtcl.DensityReducedGvecs[iG] =
        dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
    app_log() << "  Read " << numG << " density G-vectors.\n";
    if (TargetPtcl.DensityReducedGvecs.size())
    {
      app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
      std::vector<std::complex<double> > Density_G_DP;
      HDFAttribIO<std::vector<std::complex<double> > > h_density_G (Density_G_DP);
      h_density_G.read (H5FileID, "/density/rho_G");
      TargetPtcl.Density_G.assign(Density_G_DP.begin(),Density_G_DP.end());
      if (!TargetPtcl.Density_G.size())
      {
        app_error() << "  Density reduced G-vectors defined, but not the"
                    << " density.\n";
        abort();
      }
    }
  }
  return true;
}


std::string
EinsplineSetBuilder::OrbitalPath(int ti, int bi)
{
  std::string eigenstatesGroup;
  if (Version[0]==0 && Version[1]== 11)
    eigenstatesGroup = "/eigenstates_3";
  else
    if (Version[0]==0 && Version[1]==20)
      eigenstatesGroup = "/eigenstates";
  std::ostringstream groupPath;
  if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
    groupPath << eigenstatesGroup << "/twist_"
              << ti << "/band_" << bi << "/";
  else
    if (NumBands > 1)
      groupPath << eigenstatesGroup << "/twist/band_" << bi << "/";
    else
      groupPath << eigenstatesGroup << "/twist/band/";
  return groupPath.str();
}

std::string
EinsplineSetBuilder::CoreStatePath(int ti, int cs)
{
  std::string eigenstatesGroup;
  if (Version[0]==0 && Version[1]== 11)
    eigenstatesGroup = "/eigenstates_3";
  else
    if (Version[0]==0 && Version[1]==20)
      eigenstatesGroup = "/eigenstates";
  std::ostringstream groupPath;
  if ((Version[0]==0 && Version[1]==11) || NumTwists > 1)
    groupPath << eigenstatesGroup << "/twist_"
              << ti << "/core_state_" << cs << "/";
  else
    if (NumBands > 1)
      groupPath << eigenstatesGroup << "/twist/core_state_" << cs << "/";
    else
      groupPath << eigenstatesGroup << "/twist/core_state/";
  return groupPath.str();
}

std::string
EinsplineSetBuilder::MuffinTinPath(int ti, int bi, int tin)
{
  std::ostringstream groupPath;
  if (NumMuffinTins > 0)
    groupPath << OrbitalPath(ti,bi) << "muffin_tin_" << tin << "/";
  else
    groupPath << OrbitalPath(ti,bi) << "muffin_tin/";
  return groupPath.str();
}

#if 0
void
EinsplineSetBuilder::ReadBands
(int spin, EinsplineSetExtended<std::complex<double> >* orbitalSet)
{
  update_token(__FILE__,__LINE__,"ReadBands:complex");
  bool root = myComm->rank()==0;
  //bcastwith other stuff
  myComm->bcast(NumDistinctOrbitals);
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
  orbitalSet->phase.resize(N);
  orbitalSet->eikr.resize(N);
  orbitalSet->NumValenceOrbs = NumValenceOrbs;
  orbitalSet->NumCoreOrbs    = NumCoreOrbs;

  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  if (root)
  {
    int numOrbs = orbitalSet->getOrbitalSetSize();
    int num = 0;
    for (int iorb=0; iorb<N; iorb++)
    {
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
  //   std::cerr << "Using existing copy of multi_UBspline_3d_z for "
  // 	   << "thread number " << omp_get_thread_num() << ".\n";
  //   orbitalSet->MultiSpline = iter->second;
  //   return;
  // }
  int nx, ny, nz, bi, ti;
  Array<std::complex<double>,3> splineData, rawData;
  if (root)
  {
    // Find the orbital mesh size
    int i=0;
    while (SortBands[i].IsCoreState)
      i++;
    ti = SortBands[i].TwistIndex;
    bi = SortBands[i].BandIndex;
    std::string vectorName = OrbitalPath (ti, bi) + "eigenvector";
    HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
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
  nx=nxyz[0];
  ny=nxyz[1];
  nz=nxyz[2];
  if (!root)
    splineData.resize(nx-1,ny-1,nz-1);
  myComm->bcast(splineData);
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
  x_grid.num = nx-1;
  y_grid.start = 0.0;
  y_grid.end = 1.0;
  y_grid.num = ny-1;
  z_grid.start = 0.0;
  z_grid.end = 1.0;
  z_grid.num = nz-1;
  // Create the multiUBspline object
  orbitalSet->MultiSpline =
    create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  set_multi_UBspline_3d_z (orbitalSet->MultiSpline, 0, splineData.data());
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
  int iorb  = 0;
  int icore = 0;
  int ival = 0;
  while (iorb < N)
  {
    bool isCore;
    if (root)
      isCore = SortBands[iorb].IsCoreState;
    myComm->bcast (isCore);
    if (isCore)
    {
      int atom, l, m=0;
      double rmax;
      Vector<double> g, r;
      PosType twist, k;
      if (root)
      {
        ti   = SortBands[iorb].TwistIndex;
        bi   = SortBands[iorb].BandIndex;
        double e = SortBands[iorb].Energy;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        std::string atomName = CoreStatePath (ti, bi) + "atom";
        std::string gName    = CoreStatePath (ti, bi) + "g";
        std::string rMaxName = CoreStatePath (ti, bi) + "rmax";
        std::string lName    = CoreStatePath (ti, bi) + "l";
        std::string kName    = CoreStatePath (ti, bi) + "k";
        std::string rName    = CoreStatePath (ti, bi) + "r";
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
      myComm->bcast (atom);
      myComm->bcast(rmax);
      myComm->bcast (l);
      myComm->bcast (k);
      int ng = g.size();
      myComm->bcast(ng);
      if (g.size() != ng)
      {
        g.resize(ng);
        r.resize(ng);
      }
      myComm->bcast (g);
      myComm->bcast (r);
      double Z = (double)IonTypes[atom];
      orbitalSet->MuffinTins[atom].addCore (l, m, r, g, k, Z);
      icore++;
    }
    else
    {
      if (root)
      {
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        double e = SortBands[iorb].Energy;
        PosType twist, k;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        char vs[256];
        sprintf (vs, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n",
                 ti, bi, e, k[0], k[1], k[2], myComm->rank());
        app_log() << vs << std::endl;

        std::string vectorName = OrbitalPath (ti, bi) + "eigenvector";
        HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
        h_rawData.read(H5FileID, vectorName.c_str());
        if ((rawData.size(0) != nx) ||
            (rawData.size(1) != ny) ||
            (rawData.size(2) != nz))
        {
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
        double Z = (double)IonTypes[tin];
        orbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
      }
      ival++;
    } // valence state
    iorb++;
  }
  //ExtendedMap_z[set] = orbitalSet->MultiSpline;
}

void
EinsplineSetBuilder::ReadBands
(int spin, EinsplineSetExtended<double>* orbitalSet)
{
  update_token(__FILE__,__LINE__,"ReadBands:double");
  std::vector<BandInfo>& SortBands(*FullBands[spin]);
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
  orbitalSet->phase.resize(N);
  orbitalSet->eikr.resize(N);
  orbitalSet->NumValenceOrbs = NumValenceOrbs;
  orbitalSet->NumCoreOrbs    = NumCoreOrbs;
  // Read in k-points
  int numOrbs = orbitalSet->getOrbitalSetSize();
  int num = 0;
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
  }
  myComm->bcast(orbitalSet->kPoints);
  myComm->bcast(orbitalSet->MakeTwoCopies);
  myComm->bcast(orbitalSet->HalfG);
  // First, check to see if we have already read this in
  H5OrbSet set(H5FileName, spin, N);
  int nx, ny, nz, bi, ti;
  Array<std::complex<double>,3> rawData;
  Array<double,3>         splineData;
  if (root)
  {
    // Find the orbital mesh size
    int i=0;
    while (SortBands[i].IsCoreState)
      i++;
    ti = SortBands[i].TwistIndex;
    bi = SortBands[i].BandIndex;
    std::string vectorName = OrbitalPath (ti, bi) + "eigenvector";
    HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
    h_rawData.read(H5FileID, vectorName.c_str());
    nx = rawData.size(0);
    ny = rawData.size(1);
    nz = rawData.size(2);
    splineData.resize(nx-1, ny-1, nz-1);
    PosType ru;
    for (int ix=0; ix<(nx-1); ix++)
    {
      ru[0] = (RealType)ix / (RealType)(nx-1);
      for (int iy=0; iy<(ny-1); iy++)
      {
        ru[1] = (RealType)iy / (RealType)(ny-1);
        for (int iz=0; iz<(nz-1); iz++)
        {
          ru[2] = (RealType)iz / (RealType)(nz-1);
          double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
          double s, c;
          sincos(phi, &s, &c);
          std::complex<double> phase(c,s);
          std::complex<double> z = phase*rawData(ix,iy,iz);
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
  nx=nxyz[0];
  ny=nxyz[1];
  nz=nxyz[2];
  if (!root)
    splineData.resize(nx-1,ny-1,nz-1);
  myComm->bcast(splineData);
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
  x_grid.num = nx-1;
  y_grid.start = 0.0;
  y_grid.end = 1.0;
  y_grid.num = ny-1;
  z_grid.start = 0.0;
  z_grid.end = 1.0;
  z_grid.num = nz-1;
  // Create the multiUBspline object
  orbitalSet->MultiSpline =
    create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  set_multi_UBspline_3d_d (orbitalSet->MultiSpline, 0, splineData.data());
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
  int iorb  = 0;
  int icore = 0;
  int ival = 0;
  while (iorb < N)
  {
    bool isCore;
    if (root)
      isCore = SortBands[iorb].IsCoreState;
    myComm->bcast (isCore);
    if (isCore)
    {
      int atom, l, m=0;
      double rmax;
      Vector<double> g, r;
      PosType twist, k;
      if (root)
      {
        ti   = SortBands[iorb].TwistIndex;
        bi   = SortBands[iorb].BandIndex;
        double e = SortBands[iorb].Energy;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        std::string atomName = CoreStatePath (ti, bi) + "atom";
        std::string gName    = CoreStatePath (ti, bi) + "g";
        std::string rMaxName = CoreStatePath (ti, bi) + "rmax";
        std::string lName    = CoreStatePath (ti, bi) + "l";
        std::string kName    = CoreStatePath (ti, bi) + "k";
        std::string rName    = CoreStatePath (ti, bi) + "r";
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
      myComm->bcast (atom);
      myComm->bcast(rmax);
      myComm->bcast (l);
      myComm->bcast (k);
      int ng = g.size();
      myComm->bcast(ng);
      if (g.size() != ng)
      {
        g.resize(ng);
        r.resize(ng);
      }
      myComm->bcast (g);
      myComm->bcast (r);
      double Z = (double)IonTypes[atom];
      orbitalSet->MuffinTins[atom].addCore (l, m, r, g, k, Z);
      icore++;
    }
    else
    {
      if (root)
      {
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        double e = SortBands[iorb].Energy;
        PosType twist, k;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n",
                 ti, bi, e, k[0], k[1], k[2], myComm->rank());
        std::string vectorName = OrbitalPath (ti, bi) + "eigenvector";
        HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
        h_rawData.read(H5FileID, vectorName.c_str());
        if ((rawData.size(0) != nx) ||
            (rawData.size(1) != ny) ||
            (rawData.size(2) != nz))
        {
          fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
          fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
          abort();
        }
        PosType ru;
        for (int ix=0; ix<(nx-1); ix++)
        {
          ru[0] = (RealType)ix / (RealType)(nx-1);
          for (int iy=0; iy<(ny-1); iy++)
          {
            ru[1] = (RealType)iy / (RealType)(ny-1);
            for (int iz=0; iz<(nz-1); iz++)
            {
              ru[2] = (RealType)iz / (RealType)(nz-1);
              double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
              double s, c;
              sincos(phi, &s, &c);
              std::complex<double> phase(c,s);
              std::complex<double> z = phase*rawData(ix,iy,iz);
              splineData(ix,iy,iz) = z.real();
            }
          }
        }
      }
      myComm->bcast(splineData);
      set_multi_UBspline_3d_d
      (orbitalSet->MultiSpline, ival, splineData.data());
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
        double Z = (double)IonTypes[tin];
        orbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
      }
      ival++;
    } // valence state
    iorb++;
  }
  //ExtendedMap_d[set] = orbitalSet->MultiSpline;
}
#endif

}

