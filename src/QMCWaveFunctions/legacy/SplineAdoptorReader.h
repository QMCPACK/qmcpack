//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SplineAdoptorReader.h
 *
 * The most general reader class for SplineAdoptor using the full single grid for the supercell
 * - SplineR2RAdoptor
 * - SplineC2CPackedAdoptor
 * - SplineC2RPackedAdoptor
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READER_H
#define QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READER_H
namespace qmcplusplus
{

/** base class to read data and manage spline tables
 *
 * Each SplineAdoptor needs a reader derived from BsplineReaderBase.
 * This base class handles common chores
 * - check_twists : read gvectors, set twists for folded bands if needed, and set the phase for the special K
 * - set_grid : create the basic grid and boundary conditions for einspline
 * Note that template is abused but it works.
 */
struct BsplineReaderBase
{
  EinsplineSetBuilder* mybuilder;
  Communicate* myComm;

  BsplineReaderBase(EinsplineSetBuilder* e):mybuilder(e)
  {
    myComm=mybuilder->getCommunicator();
  }

  /** copy minimal informatino from EinsplineSet to manage SplineAdoptor
   */
  template<typename SPE>
  void init(EinsplineSet* in, SPE* out)
  {
    out->PrimLattice=in->PrimLattice;
    out->SuperLattice=in->SuperLattice;
    out->GGt=in->GGt;
    out->setOrbitalSetSize(in->getOrbitalSetSize());
  }

  /** read gvectors and set the mesh, and prepare for einspline
   */
  template<typename GT, typename BCT>
  inline bool set_grid(const TinyVector<int,3>& halfg, GT* xyz_grid, BCT* xyz_bc)
  {
    bool havePsig=mybuilder->ReadGvectors_ESHDF();
    xyz_grid[0].start = 0.0;
    xyz_grid[0].end = 1.0;
    xyz_grid[0].num = mybuilder->MeshSize[0];
    xyz_grid[1].start = 0.0;
    xyz_grid[1].end = 1.0;
    xyz_grid[1].num = mybuilder->MeshSize[1];
    xyz_grid[2].start = 0.0;
    xyz_grid[2].end = 1.0;
    xyz_grid[2].num = mybuilder->MeshSize[2];
    for(int j=0; j<3; ++j)
    {
      if(halfg[j])
      {
        xyz_bc[j].lCode=ANTIPERIODIC;
        xyz_bc[j].rCode=ANTIPERIODIC;
      }
      else
      {
        xyz_bc[j].lCode=PERIODIC;
        xyz_bc[j].rCode=PERIODIC;
      }
    }
    return havePsig;
  }

  /** initialize twist-related data
   */
  template<typename SPE>
  inline void check_twists(EinsplineSet* orbitalSet, SPE* bspline)
  {
    //init(orbitalSet,bspline);
    bspline->PrimLattice=orbitalSet->PrimLattice;
    bspline->SuperLattice=orbitalSet->SuperLattice;
    bspline->GGt=orbitalSet->GGt;
    bspline->setOrbitalSetSize(orbitalSet->getOrbitalSetSize());
    int N = mybuilder->NumDistinctOrbitals;
    int NumValenceOrbs = mybuilder->NumValenceOrbs;
    bspline->resizeStorage(N,NumValenceOrbs);
    // Read in k-points
    int numOrbs=bspline->getOrbitalSetSize();
    int num = 0;
    const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
    for (int iorb=0; iorb<N; iorb++)
    {
      int ti = SortBands[iorb].TwistIndex;
      bspline->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(mybuilder->TwistAngles[ti]); //twist);
      bspline->MakeTwoCopies[iorb] =
        (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
      num += bspline->MakeTwoCopies[iorb] ? 2 : 1;
    }
    app_log() << "NumDistinctOrbitals " << mybuilder->NumDistinctOrbitals
              << " NumValenceOrbs " << mybuilder->NumValenceOrbs
              << " numOrbs = " << numOrbs << std::endl;
    bspline->HalfG=0;
    TinyVector<int,3> bconds=mybuilder->TargetPtcl.Lattice.BoxBConds;
    if(!bspline->is_complex)
    {
      //no k-point folding, single special k point (G, L ...)
      int ti=mybuilder->SortBands[0].TwistIndex;
      TinyVector<double,3> twist0 = mybuilder->TwistAngles[mybuilder->SortBands[0].TwistIndex];
      for (int i=0; i<3; i++)
        if (bconds[i] && ((std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)))
          bspline->HalfG[i] = 1;
        else
          bspline->HalfG[i] = 0;
      app_log() << "  TwistIndex = " << mybuilder->SortBands[0].TwistIndex << " TwistAngle " << twist0 << std::endl;
      app_log() <<"   HalfG = " << bspline->HalfG << std::endl;
    }
    app_log().flush();
  }

  /** return the path name in hdf5
   */
  inline std::string psi_g_path(int ti, int spin, int ib)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti
         << "/spin_" << spin << "/state_" << ib << "/psi_g";
    return path.str();
  }

  /** return the path name in hdf5
   */
  inline std::string psi_r_path(int ti, int spin, int ib)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti
         << "/spin_" << spin << "/state_" << ib << "/psi_r";
    return path.str();
  }

  /** read/bcast psi_g
   * @param ti twist index
   * @param spin spin index
   * @param ib band index
   * @param cG psi_g as stored in hdf5
   */
  void get_psi_g(int ti, int spin, int ib, Vector<std::complex<double> >& cG)
  {
    int ncg=0;
    if(myComm->rank()==0)
    {
      std::string path=psi_g_path(ti,spin,ib);
      HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
      h_cG.read (mybuilder->H5FileID, path.c_str());
      ncg=cG.size();
    }
    myComm->bcast(ncg);
    if(ncg != mybuilder->MaxNumGvecs)
    {
      APP_ABORT("Failed : ncg != MaxNumGvecs");
    }
    myComm->bcast(cG);
  }

  virtual ~BsplineReaderBase() {}

  /** create the actual spline sets
   */
  virtual SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalset)=0;
};

/** General SplineAdoptorReader to handle any unitcell
 */
template<typename SA>
struct SplineAdoptorReader: public BsplineReaderBase
{
  typedef SA adoptor_type;
  typedef typename adoptor_type::DataType    DataType;
  typedef typename adoptor_type::SplineType SplineType;

  SplineAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
  {}

  SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
  {
    ReportEngine PRE("SplineC2XAdoptorReader","create_spline_set(spin,SPE*)");
    Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
    double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;
    BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
    app_log() << "  AdoptorName = " << bspline->AdoptorName << std::endl;
    if(bspline->is_complex)
      app_log() << "  Using complex einspline table" << std::endl;
    else
      app_log() << "  Using real einspline table" << std::endl;
    //baseclass handles twists
    check_twists(orbitalSet,bspline);
    Ugrid xyz_grid[3];
    typename adoptor_type::BCType xyz_bc[3];
    bool havePsig=set_grid(bspline->HalfG,xyz_grid, xyz_bc);
    if(!havePsig)
    {
      APP_ABORT("EinsplineAdoptorReader needs psi_g. Set precision=\"double\".");
    }
    bspline->create_spline(xyz_grid,xyz_bc);
    int TwistNum = mybuilder->TwistNum;
    std::string splinefile
    =make_spline_filename(mybuilder->H5FileName,mybuilder->TileMatrix
                          ,spin,TwistNum,mybuilder->MeshSize);
    bool root=(myComm->rank() == 0);
    int foundspline=0;
    Timer now;
    if(root)
    {
      hdf_archive h5f;
      foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
      if(foundspline)
      {
        std::string aname("none");
        foundspline = h5f.read(aname,"adoptor_name");
        foundspline = (aname.find(bspline->KeyWord) != std::string::npos);
      }
      if(foundspline)
      {
        int sizeD=0;
        foundspline=h5f.read(sizeD,"sizeof");
        foundspline = (sizeD == sizeof(typename adoptor_type::DataType));
      }
      if(foundspline)
        bspline->read_splines(h5f);
    }
    myComm->bcast(foundspline);
    t_h5 = now.elapsed();
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << std::endl;
      chunked_bcast(myComm, bspline->MultiSpline);
      t_init+=now.elapsed();
    }
    else
    {
      int N = mybuilder->NumDistinctOrbitals;
      int NumValenceOrbs = mybuilder->NumValenceOrbs;
      const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
      int nx=mybuilder->MeshSize[0];
      int ny=mybuilder->MeshSize[1];
      int nz=mybuilder->MeshSize[2];
      /** For valence orbitals,
       * - extended orbitals either in G or in R
       * - localized orbitals
       */
      Array<DataType,3> splineData_r(nx,ny,nz),splineData_i(ny,ny,nz);
      if(havePsig)//perform FFT using FFTW
      {
        c_init.restart();
        Array<std::complex<double>,3> FFTbox;
        FFTbox.resize(nx, ny, nz);
        fftw_plan FFTplan = fftw_plan_dft_3d(nx, ny, nz,
                                             reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                             reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                             +1, FFTW_ESTIMATE);
        Vector<std::complex<double> > cG(mybuilder->MaxNumGvecs);
        //this will be parallelized with OpenMP
        for(int iorb=0; iorb<N; ++iorb)
        {
          int ti=SortBands[iorb].TwistIndex;
          get_psi_g(ti,spin,SortBands[iorb].BandIndex,cG);
          c_unpack.restart();
          unpack4fftw(cG,mybuilder->Gvecs[0],mybuilder->MeshSize,FFTbox);
          t_unpack+= c_unpack.elapsed();
          c_fft.restart();
          fftw_execute (FFTplan);
          t_fft+= c_fft.elapsed();
          c_phase.restart();
          if(bspline->is_complex)
            fix_phase_rotate_c2c(FFTbox,splineData_r, splineData_i,mybuilder->TwistAngles[ti]);
          else
            fix_phase_rotate_c2r(FFTbox,splineData_r, mybuilder->TwistAngles[ti]);
          t_phase+= c_phase.elapsed();
          c_spline.restart();
          bspline->set_spline(splineData_r.data(),splineData_i.data(),ti,iorb,0);
          t_spline+= c_spline.elapsed();
        }
        fftw_destroy_plan(FFTplan);
        t_init+=c_init.elapsed();
      }
      //else
      //{
      //  Array<std::complex<double>,3> rawData(nx,ny,nz);
      //  //this will be parallelized with OpenMP
      //  for(int iorb=0; iorb<N; ++iorb)
      //  {
      //    //check dimension
      //    if(root)
      //    {
      //      std::string path=psi_r_path(SortBands[iorb].TwistIndex,spin,SortBands[iorb].BandIndex);
      //      HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(rawData);
      //      h_splineData.read(mybuilder->H5FileID, path.c_str());
      //      simd::copy(splineData_r.data(),splineData_i.data(),rawData.data(),rawData.size());
      //    }
      //    myComm->bcast(splineData_r);
      //    myComm->bcast(splineData_i);
      //    bspline->set_spline(splineData_r.data(),splineData_i.data(),iorb);
      //  }
      //}
      if(qmc_common.save_wfs && root)
      {
        hdf_archive h5f;
        h5f.create(splinefile);
        h5f.write(bspline->AdoptorName,"adoptor_name");
        int sizeD=sizeof(typename adoptor_type::DataType);
        h5f.write(sizeD,"sizeof");
        bspline->write_splines(h5f);
      }
    }
    app_log() << "    READBANDS::PREP   = " << t_prep << std::endl;
    app_log() << "    READBANDS::H5     = " << t_h5 << std::endl;
    app_log() << "    READBANDS::UNPACK = " << t_unpack << std::endl;
    app_log() << "    READBANDS::FFT    = " << t_fft << std::endl;
    app_log() << "    READBANDS::PHASE  = " << t_phase << std::endl;
    app_log() << "    READBANDS::SPLINE = " << t_spline << std::endl;
    app_log() << "    READBANDS::SUM    = " << t_init << std::endl;
    return bspline;
  }
};
}
#endif
