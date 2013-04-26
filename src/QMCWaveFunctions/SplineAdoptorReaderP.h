//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file SplineAdoptorReader.h
 *
 * The most general reader class for SplineAdoptor using the full single grid for the supercell
 * - SplineR2RAdoptor
 * - SplineC2CPackedAdoptor
 * - SplineC2RPackedAdoptor
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READER_H
#define QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READER_H
#include <mpi/collectives.h>
#include <mpi/point2point.h>
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
              << " numOrbs = " << numOrbs << endl;
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
      app_log() << "  TwistIndex = " << mybuilder->SortBands[0].TwistIndex << " TwistAngle " << twist0 << endl;
      app_log() <<"   HalfG = " << bspline->HalfG << endl;
    }
    app_log().flush();
  }

  /** return the path name in hdf5
   */
  inline string psi_g_path(int ti, int spin, int ib)
  {
    ostringstream path;
    path << "/electrons/kpoint_" << ti
         << "/spin_" << spin << "/state_" << ib << "/psi_g";
    return path.str();
  }

  /** return the path name in hdf5
   */
  inline string psi_r_path(int ti, int spin, int ib)
  {
    ostringstream path;
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
  void get_psi_g(int ti, int spin, int ib, Vector<complex<double> >& cG)
  {
    int ncg=0;
    if(myComm->rank()==0)
    {
      string path=psi_g_path(ti,spin,ib);
      HDFAttribIO<Vector<complex<double> > >  h_cG(cG);
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

  TinyVector<int,3> MeshSize;
  Array<complex<double>,3> FFTbox;
  Array<double,3> splineData_r, splineData_i;
  vector<UBspline_3d_d*> spline_r;
  vector<UBspline_3d_d*> spline_i;
  BsplineSet<adoptor_type>* bspline;
  vector<int> OrbGroups;
  fftw_plan FFTplan;

  SplineAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e), FFTplan(NULL)
  {}

  ~SplineAdoptorReader()
  {
    for(int i=0; i<spline_r.size(); ++i)
      free(spline_r[i]);
    for(int i=0; i<spline_i.size(); ++i)
      free(spline_i[i]);
    fftw_destroy_plan(FFTplan);
  }

  SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
  {
    ReportEngine PRE("SplineC2XAdoptorReader","create_spline_set(spin,SPE*)");
    //Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
    //double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;
    bspline=new BsplineSet<adoptor_type>;
    app_log() << "  AdoptorName = " << bspline->AdoptorName << endl;
    if(bspline->is_complex)
      app_log() << "  Using complex einspline table" << endl;
    else
      app_log() << "  Using real einspline table" << endl;
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
    string splinefile
    =make_spline_filename(mybuilder->H5FileName,mybuilder->TileMatrix
                          ,spin,TwistNum,mybuilder->MeshSize);
    bool root=(myComm->rank() == 0);
    int foundspline=0;
    Timer now;
    if(root)
    {
      now.restart();
      hdf_archive h5f(myComm);
      foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
      if(foundspline)
      {
        string aname("none");
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
      app_log() << "  Time to read the table in " << splinefile << " = " << now.elapsed() << endl;;
    }
    myComm->bcast(foundspline);
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << endl;
      now.restart();
      chunked_bcast(myComm, bspline->MultiSpline);
      app_log() << "  SplineAdoptorReader bcast the full table " << now.elapsed() << " sec" << endl;
      app_log().flush();
    }
    else
    {
      int N = mybuilder->NumDistinctOrbitals;
      int NumValenceOrbs = mybuilder->NumValenceOrbs;
      int nx=mybuilder->MeshSize[0];
      int ny=mybuilder->MeshSize[1];
      int nz=mybuilder->MeshSize[2];
      if(havePsig)//perform FFT using FFTW
      {
        int np=std::min(N,myComm->size());
        OrbGroups.resize(np+1,0);
        FairDivideLow(N,np,OrbGroups);
        int ip=myComm->rank();
        int norbs_n=(ip<np)?OrbGroups[ip+1]-OrbGroups[ip]:0;
        TinyVector<double,3> start(0.0);
        TinyVector<double,3> end(1.0);
        splineData_r.resize(nx,ny,nz);
        if(bspline->is_complex)
          splineData_i.resize(nx,ny,nz);
        MeshSize=mybuilder->MeshSize;
        UBspline_3d_d* dummy=0;
        spline_r.resize(norbs_n+1);
        for(int i=0; i<spline_r.size(); ++i)
          spline_r[i]=einspline::create(dummy,start,end,MeshSize,bspline->HalfG);
        if(bspline->is_complex)
        {
          spline_i.resize(norbs_n+1);
          for(int i=0; i<spline_i.size(); ++i)
            spline_i[i]=einspline::create(dummy,start,end,MeshSize,bspline->HalfG);
        }
        FFTbox.resize(nx, ny, nz);
        FFTplan = fftw_plan_dft_3d(nx, ny, nz,
                                   reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                   reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                   +1, FFTW_ESTIMATE);
        //now.restart();
        //initialize_spline_pio_bcast(spin);
        //app_log() << "  SplineAdoptorReader initialize_spline_pio_bcast " << now.elapsed() << " sec" << endl;
        size_t ntot=size_t(nx*ny*nz)*size_t(N);
        if(ntot>>22) //Using 4M as the cutoff, candidate for autotuning
        {
          now.restart();
          initialize_spline_pio(spin);
          app_log() << "  SplineAdoptorReader initialize_spline_pio " << now.elapsed() << " sec" << endl;
        }
        else
        {
          now.restart();
          initialize_spline_slow(spin);
          app_log() << "  SplineAdoptorReader initialize_spline_slow " << now.elapsed() << " sec" << endl;
        }
      }
      else//why, don't know
        initialize_spline_psi_r(spin);
      if(qmc_common.save_wfs && root)
      {
        now.restart();
        hdf_archive h5f;
        h5f.create(splinefile);
        h5f.write(bspline->AdoptorName,"adoptor_name");
        int sizeD=sizeof(typename adoptor_type::DataType);
        h5f.write(sizeD,"sizeof");
        bspline->write_splines(h5f);
        app_log() << "  SplineAdoptorReader dump " << now.elapsed() << " sec" << endl;
      }
    }
    return bspline;
  }

  /** fft and spline cG
   * @param cG psi_g to be processed
   * @param ti twist index
   * @param iorb orbital index
   *
   * Perform FFT and spline to spline_r[iorb] and spline_i[iorb]
   */
  inline void fft_spline(Vector<complex<double> >& cG, int ti, int iorb)
  {
    unpack4fftw(cG,mybuilder->Gvecs[0],mybuilder->MeshSize,FFTbox);
    fftw_execute (FFTplan);
    if(bspline->is_complex)
    {
      fix_phase_rotate_c2c(FFTbox,splineData_r, splineData_i,mybuilder->TwistAngles[ti]);
      einspline::set(spline_r[iorb],splineData_r.data());
      einspline::set(spline_i[iorb],splineData_i.data());
    }
    else
    {
      fix_phase_rotate_c2r(FFTbox,splineData_r, mybuilder->TwistAngles[ti]);
      einspline::set(spline_r[iorb],splineData_r.data());
    }
  }


  void initialize_spline_slow(int spin)
  {
    int N = mybuilder->NumDistinctOrbitals;
    int NumValenceOrbs = mybuilder->NumValenceOrbs;
    int nx=mybuilder->MeshSize[0];
    int ny=mybuilder->MeshSize[1];
    int nz=mybuilder->MeshSize[2];
    Array<DataType,3> data_r(nx,ny,nz), data_i;
    if(bspline->is_complex)
      data_i.resize(ny,ny,nz);
    Vector<complex<double> > cG(mybuilder->MaxNumGvecs);
    const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
    //this will be parallelized with OpenMP
    for(int iorb=0; iorb<N; ++iorb)
    {
      int ti=SortBands[iorb].TwistIndex;
      get_psi_g(ti,spin,SortBands[iorb].BandIndex,cG);//bcast cG
      //fft_spline(cG,ti,0);
      //bspline->set_spline(spline_r[0],spline_i[0],iorb);
      unpack4fftw(cG,mybuilder->Gvecs[0],mybuilder->MeshSize,FFTbox);
      fftw_execute (FFTplan);
      if(bspline->is_complex)
        fix_phase_rotate_c2c(FFTbox,data_r, data_i,mybuilder->TwistAngles[ti]);
      else
        fix_phase_rotate_c2r(FFTbox,data_r, mybuilder->TwistAngles[ti]);
      bspline->set_spline(data_r.data(),data_i.data(),ti,iorb,0);
    }
  }

  void initialize_spline_pio(int spin)
  {
    int N = mybuilder->NumDistinctOrbitals;
    bool root=(myComm->rank()==0);
    bool foundit=true;
    int np=OrbGroups.size()-1;
    if(myComm->rank()<np)
    {
      int iorb_first=OrbGroups[myComm->rank()];
      int iorb_last =OrbGroups[myComm->rank()+1];
      hdf_archive h5f(myComm,false);
      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);
      Vector<complex<double> > cG(mybuilder->Gvecs[0].size());;
      const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=SortBands[iorb].TwistIndex;
        string s=psi_g_path(ti,spin,SortBands[iorb].BandIndex);
        foundit &= h5f.read(cG,s);
        fft_spline(cG,ti,ib);
      }
      if(root)
      {
        for(int iorb=OrbGroups[0],ib=0; iorb<OrbGroups[1]; ++iorb,++ib)
          bspline->set_spline(spline_r[ib],spline_i[ib],SortBands[iorb].TwistIndex,iorb,0);
      }
    }//root put splines to the big table
    myComm->barrier();
    myComm->bcast(foundit);
    if(!foundit)
      APP_ABORT("SplineAdoptorReader Failed to read band(s)");
    //mpi needs int
    int ng_big=static_cast<int>(spline_r[0]->coefs_size);
    //send back to zero  using synchronous send/recv
    for(int ip=1; ip<np; ++ip)
    {
      int remote=0;
      if(ip==myComm->rank())
      {
        for(int iorb=OrbGroups[ip],ib=0; iorb<OrbGroups[ip+1]; ++iorb,++ib)
        {
          mpi::send(*myComm,spline_r[ib]->coefs,ng_big,remote,iorb);
          if(bspline->is_complex)
            mpi::send(*myComm,spline_i[ib]->coefs,ng_big,remote,iorb+N);
        }
      }
      else
        if(root)
        {
          for(int iorb=OrbGroups[ip],ib=0; iorb<OrbGroups[ip+1]; ++iorb,++ib)
          {
            mpi::recv(*myComm,spline_r[ib]->coefs,ng_big,ip,iorb);
            if(bspline->is_complex)
              mpi::recv(*myComm,spline_i[ib]->coefs,ng_big,ip,iorb+N);
          }
          const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
          for(int iorb=OrbGroups[ip],ib=0; iorb<OrbGroups[ip+1]; ++iorb,++ib)
          {
            bspline->set_spline(spline_r[ib],spline_i[ib],SortBands[iorb].TwistIndex,iorb,0);
          }
        }
    }
    myComm->barrier();
    chunked_bcast(myComm, bspline->MultiSpline);
  }

  void initialize_spline_pio_bcast(int spin)
  {
    bool root=(myComm->rank()==0);
    int np=OrbGroups.size()-1;
    bool foundit=true;
    if(myComm->rank()<np)
    {
      int iorb_first=OrbGroups[myComm->rank()];
      int iorb_last =OrbGroups[myComm->rank()+1];
      hdf_archive h5f(myComm,false);
      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);
      Vector<complex<double> > cG(mybuilder->Gvecs[0].size());;
      const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=SortBands[iorb].TwistIndex;
        string s=psi_g_path(ti,spin,SortBands[iorb].BandIndex);
        foundit &= h5f.read(cG,s);
        fft_spline(cG,ti,ib);
      }
    }
    myComm->barrier();
    myComm->bcast(foundit);
    if(!foundit)
      APP_ABORT("SplineAdoptorReader Failed to read band(s)");
    int ng_big=static_cast<int>(spline_r[0]->coefs_size);
    //pointers to UBspline_3d_d for bcast without increaing mem
    UBspline_3d_d *dense_r, *dense_i;
    int iorb_target=(myComm->rank()<np)?(OrbGroups[myComm->rank()+1]-OrbGroups[myComm->rank()]):0;
    for(int ip=0; ip<np; ++ip)
    {
      const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
      for(int iorb=OrbGroups[ip], ib=0; iorb<OrbGroups[ip+1]; ++iorb, ++ib)
      {
        if(ip==myComm->rank())
        {
          dense_r=spline_r[ib];
          if(bspline->is_complex)
            dense_i=spline_i[ib];
        }
        else
        {
          //everyone else
          dense_r=spline_r[iorb_target];
          if(bspline->is_complex)
            dense_i=spline_i[iorb_target];
        }
        mpi::bcast(*myComm,dense_r->coefs,ng_big,ip);
        if(bspline->is_complex)
          mpi::bcast(*myComm,dense_i->coefs,ng_big,ip);
        bspline->set_spline(dense_r,dense_i,SortBands[iorb].TwistIndex,iorb,0);
      }
    }
  }

  void initialize_spline_psi_r(int spin)
  {
    //not used by may be enabled later
    int nx=mybuilder->MeshSize[0];
    int ny=mybuilder->MeshSize[1];
    int nz=mybuilder->MeshSize[2];
    Array<DataType,3> splineData_r(nx,ny,nz),splineData_i(ny,ny,nz);
    Array<complex<double>,3> rawData(nx,ny,nz);
    //this will be parallelized with OpenMP
    int N = mybuilder->NumDistinctOrbitals;
    const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
    for(int iorb=0; iorb<N; ++iorb)
    {
      //check dimension
      if(myComm->rank() ==0)
      {
        string path=psi_r_path(SortBands[iorb].TwistIndex,spin,SortBands[iorb].BandIndex);
        HDFAttribIO<Array<complex<double>,3> >  h_splineData(rawData);
        h_splineData.read(mybuilder->H5FileID, path.c_str());
        simd::copy(splineData_r.data(),splineData_i.data(),rawData.data(),rawData.size());
      }
      mpi::bcast(*myComm,splineData_r);
      if(bspline->is_complex)
        mpi::bcast(*myComm,splineData_i);
      bspline->set_spline(splineData_r.data(),splineData_i.data(),SortBands[iorb].TwistIndex,iorb,0);
    }
  }
};
}
#endif
