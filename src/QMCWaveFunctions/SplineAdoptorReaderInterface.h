//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file SplineAdoptorReaderInterface.h
 *
 * The most general reader class for SplineAdoptor using the full single grid for the supercell
 * - SplineR2RAdoptor
 * - SplineC2CPackedAdoptor
 * - SplineC2RPackedAdoptor
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READERINTERFACE_H
#define QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READERINTERFACE_H
#include <mpi/collectives.h>
#include <mpi/point2point.h>
namespace qmcplusplus
{
/** General SplineAdoptorReaderInterface to handle any unitcell
 */
template<typename SA>
struct SplineAdoptorReaderInterface: public BsplineReaderBase
{
  typedef SA adoptor_type;
  typedef typename adoptor_type::DataType    DataType;
  typedef typename adoptor_type::SplineType SplineType;

  Array<complex<double>,3> FFTbox;
  Array<double,3> splineData_r, splineData_i;
  vector<UBspline_3d_d*> spline_r;
  vector<UBspline_3d_d*> spline_i;
  BsplineSet<adoptor_type>* bspline;
  vector<int> OrbGroups;
  fftw_plan FFTplan;
  ESInterfaceBase* interface; 
 
  SplineAdoptorReaderInterface(EinsplineSetBuilder* e, ESInterfaceBase* intfc)
    : BsplineReaderBase(e), bspline(0),FFTplan(NULL), interface(intfc)
  {}

  ~SplineAdoptorReaderInterface()
  {
    clear();
  }

  void clear()
  {
    for(int i=0; i<spline_r.size(); ++i)
    {
      free(spline_r[i]->coefs);
      free(spline_r[i]);
    }
    for(int i=0; i<spline_i.size(); ++i)
    {
      if(spline_i[i]!=0)
      {
        free(spline_i[i]->coefs);
        free(spline_i[i]);
      }
    }
    spline_r.clear();
    spline_i.clear();
    if(FFTplan!=NULL) fftw_destroy_plan(FFTplan);
    FFTplan=NULL;
  }
  SPOSetBase* create_spline_set(int spin, const BandInfoGroup& bandgroup)
  {
    ReportEngine PRE("SplineC2XAdoptorReaderInterface","create_spline_set(spin,SPE*)");
    //Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
    //double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;
    bspline=new BsplineSet<adoptor_type>;
    app_log() << "  AdoptorName = " << bspline->AdoptorName << endl;
    if(bspline->is_complex)
      app_log() << "  Using complex einspline table" << endl;
    else
      app_log() << "  Using real einspline table" << endl;

    //baseclass handles twists
    check_twists(bspline,bandgroup);

    Ugrid xyz_grid[3];

    typename adoptor_type::BCType xyz_bc[3];
    bool havePsig=set_grid_interface(bspline->HalfG,xyz_grid, xyz_bc);
    if(!havePsig)
    {
      APP_ABORT("EinsplineAdoptorReader needs psi_g. Set precision=\"double\".");
    }
    bspline->create_spline(xyz_grid,xyz_bc);
    int TwistNum = mybuilder->TwistNum;
    ostringstream oo;
    oo<<bandgroup.myName << ".g"<<MeshSize[0]<<"x"<<MeshSize[1]<<"x"<<MeshSize[2]<<".h5";
    string splinefile= oo.str(); //bandgroup.myName+".h5";
    //=make_spline_filename(mybuilder->H5FileName,mybuilder->TileMatrix,spin,TwistNum,bandgroup.GroupID,MeshSize);
    bool root=(myComm->rank() == 0);
    int foundspline=0;
    Timer now;
    if(root)
    {
      now.restart();
      hdf_archive h5f(myComm); //Ray:  HWI
      foundspline=h5f.open(splinefile,H5F_ACC_RDONLY); //Ray:  HWI
      if(foundspline)
      {
        string aname("none");
        foundspline = h5f.read(aname,"adoptor_name"); //Ray:  HWI
        foundspline = (aname.find(bspline->KeyWord) != std::string::npos); //Ray: HWI
      }
      if(foundspline)
      {
        int sizeD=0;
        foundspline=h5f.read(sizeD,"sizeof"); //Ray:  HWI
        foundspline = (sizeD == sizeof(typename adoptor_type::DataType)); //Ray: HWI
      }
      if(foundspline)
        bspline->read_splines(h5f); //Ray: HWI
      app_log() << "  Time to read the table in " << splinefile << " = " << now.elapsed() << endl;;
    }
    myComm->bcast(foundspline);
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << endl;
      now.restart();
      chunked_bcast(myComm, bspline->MultiSpline);
      app_log() << "  SplineAdoptorReaderInterface bcast the full table " << now.elapsed() << " sec" << endl;
      app_log().flush();
    }
    else
    {
      app_log()<<"create_spline_set:  Didn't find the splines...\n";
      int N=bandgroup.getNumDistinctOrbitals();
      int nx=MeshSize[0];
      int ny=MeshSize[1];
      int nz=MeshSize[2];
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
        UBspline_3d_d* dummy=0;
        spline_r.resize(norbs_n+1);
        for(int i=0; i<spline_r.size(); ++i)
          spline_r[i]=einspline::create(dummy,start,end,MeshSize,bspline->HalfG);

        spline_i.resize(norbs_n+1,0);
        if(bspline->is_complex)
        {
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
        //app_log() << "  SplineAdoptorReaderInterface initialize_spline_pio_bcast " << now.elapsed() << " sec" << endl;
        size_t ntot=size_t(nx*ny*nz)*size_t(N);
        if(ntot>>22) //Using 4M as the cutoff, candidate for autotuning
        {
          now.restart();
          initialize_spline_pio(spin,bandgroup);
          app_log() << "  SplineAdoptorReaderInterface initialize_spline_pio " << now.elapsed() << " sec" << endl;
        }
        else
        {
          now.restart();
          initialize_spline_slow(spin,bandgroup);
          app_log() << "  SplineAdoptorReaderInterface initialize_spline_slow " << now.elapsed() << " sec" << endl;
        }
        fftw_destroy_plan(FFTplan);
        FFTplan=NULL;
      }
      else//why, don't know
        initialize_spline_psi_r(spin,bandgroup);
      if(qmc_common.save_wfs && root)
      {
        now.restart();
        hdf_archive h5f; //Ray:  HWI
        h5f.create(splinefile); //Ray:  HWI
        h5f.write(bspline->AdoptorName,"adoptor_name"); //Ray: HWI
        int sizeD=sizeof(typename adoptor_type::DataType); //Ray:  HWI
        h5f.write(sizeD,"sizeof"); //Ray: HWI
        bspline->write_splines(h5f); //Ray:  HWI
        app_log() << "  SplineAdoptorReaderInterface dump " << now.elapsed() << " sec" << endl;
      }
    }

    clear();
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
    unpack4fftw(cG,mybuilder->Gvecs[0],MeshSize,FFTbox);
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


  void initialize_spline_slow(int spin, const BandInfoGroup& bandgroup)
  {
    int N=bandgroup.getNumDistinctOrbitals();
    int nx=MeshSize[0];
    int ny=MeshSize[1];
    int nz=MeshSize[2];
    Array<DataType,3> data_r(nx,ny,nz), data_i;
    if(bspline->is_complex)
      data_i.resize(ny,ny,nz);
    Vector<complex<double> > cG(mybuilder->MaxNumGvecs);
    const vector<BandInfo>& cur_bands=bandgroup.myBands;
    //this will be parallelized with OpenMP
    for(int iorb=0; iorb<N; ++iorb)
    {
      int ti=cur_bands[iorb].TwistIndex;
      get_psi_g(ti,spin,cur_bands[iorb].BandIndex,cG);//bcast cG  //Ray:  HWI
      //fft_spline(cG,ti,0);
      //bspline->set_spline(spline_r[0],spline_i[0],iorb);
      unpack4fftw(cG,mybuilder->Gvecs[0],MeshSize,FFTbox);
      fftw_execute (FFTplan);
      if(bspline->is_complex)
        fix_phase_rotate_c2c(FFTbox,data_r, data_i,mybuilder->TwistAngles[ti]);
      else
        fix_phase_rotate_c2r(FFTbox,data_r, mybuilder->TwistAngles[ti]);

      bspline->set_spline(data_r.data(),data_i.data(),ti,iorb,0);
    }
  }

  void initialize_spline_pio(int spin, const BandInfoGroup& bandgroup)
  {
    int N=bandgroup.getNumDistinctOrbitals();
    bool root=(myComm->rank()==0);
    bool foundit=true;
    int np=OrbGroups.size()-1;
    if(myComm->rank()<np)
    {
      int iorb_first=OrbGroups[myComm->rank()];
      int iorb_last =OrbGroups[myComm->rank()+1];
//      hdf_archive h5f(myComm,false);  //Ray:  HWI 
//      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY); //Ray: HWI
      Vector<complex<double> > cG(mybuilder->Gvecs[0].size());;
      const vector<BandInfo>& cur_bands=bandgroup.myBands;
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=cur_bands[iorb].TwistIndex;
       // string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex); //Ray:  HWI (Handle with interface)
        //foundit &= h5f.read(cG,s); //RAY:  HWI  (Handle with interface).
        interface->getPsi_kspace(cG, spin, cur_bands[iorb].BandIndex, ti);
        fft_spline(cG,ti,ib);
      }
      if(root)
      {
        for(int iorb=OrbGroups[0],ib=0; iorb<OrbGroups[1]; ++iorb,++ib)
        {
          bspline->set_spline(spline_r[ib],spline_i[ib],cur_bands[iorb].TwistIndex,iorb,0);
        }
      }
    }//root put splines to the big table
    myComm->barrier();
    myComm->bcast(foundit);
    if(!foundit)
      APP_ABORT("SplineAdoptorReaderInterface Failed to read band(s)");
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
      else if(root)
      {
        for(int iorb=OrbGroups[ip],ib=0; iorb<OrbGroups[ip+1]; ++iorb,++ib)
        {
          mpi::recv(*myComm,spline_r[ib]->coefs,ng_big,ip,iorb);
          if(bspline->is_complex)
            mpi::recv(*myComm,spline_i[ib]->coefs,ng_big,ip,iorb+N);
        }
        const vector<BandInfo>& cur_bands=bandgroup.myBands;
        for(int iorb=OrbGroups[ip],ib=0; iorb<OrbGroups[ip+1]; ++iorb,++ib)
        {
          bspline->set_spline(spline_r[ib],spline_i[ib],cur_bands[iorb].TwistIndex,iorb,0);
        }
      }
    }
    myComm->barrier();
    chunked_bcast(myComm, bspline->MultiSpline);
  }

  void initialize_spline_pio_bcast(int spin, const BandInfoGroup& bandgroup)
  {
    bool root=(myComm->rank()==0);
    int np=OrbGroups.size()-1;
    bool foundit=true;
    const vector<BandInfo>& cur_bands=bandgroup.myBands;

    if(myComm->rank()<np)
    {
      int iorb_first=OrbGroups[myComm->rank()];
      int iorb_last =OrbGroups[myComm->rank()+1];
    //  hdf_archive h5f(myComm,false); //Ray:  HWI
   //   h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY); //Ray: HWI
      Vector<complex<double> > cG(mybuilder->Gvecs[0].size());;
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=cur_bands[iorb].TwistIndex;
    //    string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex); //Ray:  HWI
    //    foundit &= h5f.read(cG,s); //Ray:  HWI
        interface->getPsi_kspace(cG,spin,cur_bands[iorb].BandIndex,ti); 
        fft_spline(cG,ti,ib);
      }
    }
    myComm->barrier();
    myComm->bcast(foundit);
    if(!foundit)
      APP_ABORT("SplineAdoptorReaderInterface Failed to read band(s)");
    int ng_big=static_cast<int>(spline_r[0]->coefs_size);
    //pointers to UBspline_3d_d for bcast without increaing mem
    UBspline_3d_d *dense_r, *dense_i;
    int iorb_target=(myComm->rank()<np)?(OrbGroups[myComm->rank()+1]-OrbGroups[myComm->rank()]):0;
    for(int ip=0; ip<np; ++ip)
    {
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
        bspline->set_spline(dense_r,dense_i,cur_bands[iorb].TwistIndex,iorb,0);
      }
    }
  }

  void initialize_spline_psi_r(int spin, const BandInfoGroup& bandgroup)
  {
    //not used by may be enabled later
    int nx=MeshSize[0];
    int ny=MeshSize[1];
    int nz=MeshSize[2];
    Array<DataType,3> splineData_r(nx,ny,nz),splineData_i(ny,ny,nz);
    Array<complex<double>,3> rawData(nx,ny,nz);
    const vector<BandInfo>& cur_bands=bandgroup.myBands;
    //this will be parallelized with OpenMP
    int N=bandgroup.getNumDistinctOrbitals();
    for(int iorb=0; iorb<N; ++iorb)
    {
      //check dimension
      if(myComm->rank() ==0)
      {
        string path=psi_r_path(cur_bands[iorb].TwistIndex,spin,cur_bands[iorb].BandIndex);
        HDFAttribIO<Array<complex<double>,3> >  h_splineData(rawData); //Ray:  HWI
        h_splineData.read(mybuilder->H5FileID, path.c_str()); //Ray:  HWI
        simd::copy(splineData_r.data(),splineData_i.data(),rawData.data(),rawData.size());
      }
      mpi::bcast(*myComm,splineData_r);
      if(bspline->is_complex)
        mpi::bcast(*myComm,splineData_i);
      bspline->set_spline(splineData_r.data(),splineData_i.data(),cur_bands[iorb].TwistIndex,iorb,0);
    }
  }
  template<typename GT, typename BCT>
  inline bool set_grid_interface(const TinyVector<int,3>& halfg, GT* xyz_grid, BCT* xyz_bc)
  {
    //This sets MeshSize from the input file
    bool havePsig=mybuilder->ReadGvectors_Interface(interface);

    //If this MeshSize is not initialized, use the meshsize set by the input based on FFT grid and meshfactor
    if(MeshSize[0]==0) MeshSize=mybuilder->MeshSize;

    app_log() << "  Using meshsize=" << MeshSize 
            << "\n  vs input meshsize=" << mybuilder->MeshSize << endl;

    xyz_grid[0].start = 0.0;
    xyz_grid[0].end = 1.0;
    xyz_grid[0].num = MeshSize[0];
    xyz_grid[1].start = 0.0;
    xyz_grid[1].end = 1.0;
    xyz_grid[1].num = MeshSize[1];
    xyz_grid[2].start = 0.0;
    xyz_grid[2].end = 1.0;
    xyz_grid[2].num = MeshSize[2];
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
};
}
#endif
