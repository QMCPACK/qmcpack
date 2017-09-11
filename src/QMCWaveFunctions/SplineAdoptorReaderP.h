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
//                    Jeongnim Kim, jeongnim.kim@inte.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SplineAdoptorReader.h
 *
 * The most general reader class for SplineAdoptor using the full single grid for the supercell
 * - SplineR2RAdoptor
 * - SplineC2CPackedAdoptor
 * - SplineC2RPackedAdoptor
 * Each band is initialized with UBspline_3d_d and both real and imaginary parts are passed to the adoptor object
 * which will convert the data type to their internal precision. 
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READERP_H
#define QMCPLUSPLUS_EINSPLINE_BASE_ADOPTOR_READERP_H
#include <mpi/collectives.h>
#include <mpi/point2point.h>
namespace qmcplusplus
{
/** General SplineAdoptorReader to handle any unitcell
 */
template<typename SA>
struct SplineAdoptorReader: public BsplineReaderBase
{
  typedef SA adoptor_type;
  typedef typename adoptor_type::DataType    DataType;
  typedef typename adoptor_type::SplineType SplineType;

  Array<std::complex<double>,3> FFTbox;
  Array<double,3> splineData_r, splineData_i;
  double rotate_phase_r, rotate_phase_i;
  std::vector<UBspline_3d_d*> spline_r;
  std::vector<UBspline_3d_d*> spline_i;
  BsplineSet<adoptor_type>* bspline;
  std::vector<int> OrbGroups;
  fftw_plan FFTplan;

  SplineAdoptorReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e), bspline(0),FFTplan(NULL)
  {}

  ~SplineAdoptorReader()
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

  void export_MultiSpline(multi_UBspline_3d_z** target)
  {
    *target = new multi_UBspline_3d_z;
    multi_UBspline_3d_d* source_MultiSpline = (multi_UBspline_3d_d*) bspline->MultiSpline;

    (*target)->spcode = MULTI_U3D;
    (*target)->tcode  = DOUBLE_COMPLEX;
    (*target)->coefs = (complex_double*) source_MultiSpline->coefs;
    (*target)->x_stride = source_MultiSpline->x_stride/2;
    (*target)->y_stride = source_MultiSpline->y_stride/2;
    (*target)->z_stride = source_MultiSpline->z_stride/2;

    (*target)->x_grid.start     = source_MultiSpline->x_grid.start;
    (*target)->x_grid.end       = source_MultiSpline->x_grid.end;
    (*target)->x_grid.num       = source_MultiSpline->x_grid.num;
    (*target)->x_grid.delta     = source_MultiSpline->x_grid.delta;
    (*target)->x_grid.delta_inv = source_MultiSpline->x_grid.delta_inv;
    (*target)->y_grid.start     = source_MultiSpline->y_grid.start;
    (*target)->y_grid.end       = source_MultiSpline->y_grid.end;
    (*target)->y_grid.num       = source_MultiSpline->y_grid.num;
    (*target)->y_grid.delta     = source_MultiSpline->y_grid.delta;
    (*target)->y_grid.delta_inv = source_MultiSpline->y_grid.delta_inv;
    (*target)->z_grid.start     = source_MultiSpline->z_grid.start;
    (*target)->z_grid.end       = source_MultiSpline->z_grid.end;
    (*target)->z_grid.num       = source_MultiSpline->z_grid.num;
    (*target)->z_grid.delta     = source_MultiSpline->z_grid.delta;
    (*target)->z_grid.delta_inv = source_MultiSpline->z_grid.delta_inv;

    (*target)->xBC.lCode  = source_MultiSpline->xBC.lCode;
    (*target)->xBC.rCode  = source_MultiSpline->xBC.rCode;
    (*target)->xBC.lVal_r = source_MultiSpline->xBC.lVal;
    (*target)->xBC.lVal_i = source_MultiSpline->xBC.lVal;
    (*target)->xBC.rVal_r = source_MultiSpline->xBC.rVal;
    (*target)->xBC.rVal_i = source_MultiSpline->xBC.rVal;
    (*target)->yBC.lCode  = source_MultiSpline->yBC.lCode;
    (*target)->yBC.rCode  = source_MultiSpline->yBC.rCode;
    (*target)->yBC.lVal_r = source_MultiSpline->yBC.lVal;
    (*target)->yBC.lVal_i = source_MultiSpline->yBC.lVal;
    (*target)->yBC.rVal_r = source_MultiSpline->yBC.rVal;
    (*target)->yBC.rVal_i = source_MultiSpline->yBC.rVal;
    (*target)->zBC.lCode  = source_MultiSpline->zBC.lCode;
    (*target)->zBC.rCode  = source_MultiSpline->zBC.rCode;
    (*target)->zBC.lVal_r = source_MultiSpline->zBC.lVal;
    (*target)->zBC.lVal_i = source_MultiSpline->zBC.lVal;
    (*target)->zBC.rVal_r = source_MultiSpline->zBC.rVal;
    (*target)->zBC.rVal_i = source_MultiSpline->zBC.rVal;

    (*target)->num_splines =  source_MultiSpline->num_splines/2;
    (*target)->coefs_size = source_MultiSpline->coefs_size/2;
    // (*target)->lapl3 = (complex_double*) malloc (6*sizeof(double)*(*target)->z_stride);
  }

  void export_MultiSpline(multi_UBspline_3d_d** target)
  {
    *target = (multi_UBspline_3d_d*) bspline->MultiSpline;
  }

  SPOSetBase* create_spline_set(int spin, const BandInfoGroup& bandgroup)
  {
    ReportEngine PRE("SplineC2XAdoptorReader","create_spline_set(spin,SPE*)");
    //Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
    //double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;
    bspline=new BsplineSet<adoptor_type>;
    app_log() << "  AdoptorName = " << bspline->AdoptorName << std::endl;
    if(bspline->is_complex)
      app_log() << "  Using complex einspline table" << std::endl;
    else
      app_log() << "  Using real einspline table" << std::endl;
    if(bspline->is_soa_ready)
      app_log() << "  Can use SoA implementation for mGL" << std::endl;

    //baseclass handles twists
    check_twists(bspline,bandgroup);

    Ugrid xyz_grid[3];

    typename adoptor_type::BCType xyz_bc[3];
    bool havePsig=set_grid(bspline->HalfG,xyz_grid, xyz_bc);
    if(!havePsig)
    {
      APP_ABORT("EinsplineAdoptorReader needs psi_g. Set precision=\"double\".");
    }
    bspline->create_spline(xyz_grid,xyz_bc);
    int TwistNum = mybuilder->TwistNum;
    std::ostringstream oo;
    oo<<bandgroup.myName << ".g"<<MeshSize[0]<<"x"<<MeshSize[1]<<"x"<<MeshSize[2]<<".h5";
    std::string splinefile= oo.str(); //bandgroup.myName+".h5";
    //=make_spline_filename(mybuilder->H5FileName,mybuilder->TileMatrix,spin,TwistNum,bandgroup.GroupID,MeshSize);
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
      {
        foundspline=bspline->read_splines(h5f);
        if(foundspline) app_log() << "  Time to read the table in " << splinefile << " = " << now.elapsed() << std::endl;
      }
    }
    myComm->bcast(foundspline);
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << std::endl;
      now.restart();
      chunked_bcast(myComm, bspline->MultiSpline);
      app_log() << "  SplineAdoptorReader bcast the full table " << now.elapsed() << " sec" << std::endl;
      app_log().flush();
    }
    else
    {
      int N=bandgroup.getNumDistinctOrbitals();
      int nx=MeshSize[0];
      int ny=MeshSize[1];
      int nz=MeshSize[2];
      if(havePsig)//perform FFT using FFTW
      {
        FFTbox.resize(nx, ny, nz);
        FFTplan = fftw_plan_dft_3d(nx, ny, nz,
                                   reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                   reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                   +1, FFTW_ESTIMATE);
        splineData_r.resize(nx,ny,nz);
        if(bspline->is_complex) splineData_i.resize(nx,ny,nz);

        bool usingSerialIO=(myComm->size()==1);

        int np=std::min(N,myComm->size());
        OrbGroups.resize(np+1,0);
        FairDivideLow(N,np,OrbGroups);
        int ip=myComm->rank();
        int norbs_n=(ip<np)?OrbGroups[ip+1]-OrbGroups[ip]:0;
        if(usingSerialIO)  norbs_n=0;

        TinyVector<double,3> start(0.0);
        TinyVector<double,3> end(1.0);
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

        if(usingSerialIO)
        {
          now.restart();
          initialize_spline_serial(spin, bandgroup);
          app_log() << "  SplineAdoptorReader initialize_spline_serial " << now.elapsed() << " sec" << std::endl;
        }
        else
        {
          //now.restart();
          //initialize_spline_pio_bcast(spin);
          //app_log() << "  SplineAdoptorReader initialize_spline_pio_bcast " << now.elapsed() << " sec" << std::endl;
          size_t ntot=size_t(nx*ny*nz)*size_t(N);
          //if(ntot>>22) //Using 4M as the cutoff, candidate for autotuning
          {
            now.restart();
            initialize_spline_pio(spin,bandgroup);
            //initialize_spline_pio_reduce(spin,bandgroup);
            app_log() << "  SplineAdoptorReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;
          }
        }

        fftw_destroy_plan(FFTplan);
        FFTplan=NULL;
      }
      else//why, don't know
        initialize_spline_psi_r(spin,bandgroup);
      if(qmc_common.save_wfs && root)
      {
        now.restart();
        hdf_archive h5f;
        h5f.create(splinefile);
        h5f.write(bspline->AdoptorName,"adoptor_name");
        int sizeD=sizeof(typename adoptor_type::DataType);
        h5f.write(sizeD,"sizeof");
        bspline->write_splines(h5f);
        app_log() << "  SplineAdoptorReader dump " << now.elapsed() << " sec" << std::endl;
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
  inline void fft_spline(Vector<std::complex<double> >& cG, int ti, int iorb)
  {
    unpack4fftw(cG,mybuilder->Gvecs[0],MeshSize,FFTbox);
    fftw_execute (FFTplan);
    if(bspline->is_complex)
    {
      fix_phase_rotate_c2c(FFTbox,splineData_r, splineData_i,mybuilder->TwistAngles[ti], rotate_phase_r, rotate_phase_i);
      einspline::set(spline_r[iorb],splineData_r.data());
      einspline::set(spline_i[iorb],splineData_i.data());
    }
    else
    {
      fix_phase_rotate_c2r(FFTbox,splineData_r, mybuilder->TwistAngles[ti], rotate_phase_r, rotate_phase_i);
      einspline::set(spline_r[iorb],splineData_r.data());
    }
  }


#if 0
  void initialize_spline_slow(int spin, const BandInfoGroup& bandgroup)
  {
    int N=bandgroup.getNumDistinctOrbitals();
    int nx=MeshSize[0];
    int ny=MeshSize[1];
    int nz=MeshSize[2];
    Array<DataType,3> data_r(nx,ny,nz), data_i;
    if(bspline->is_complex)
      data_i.resize(ny,ny,nz);
    Vector<std::complex<double> > cG(mybuilder->MaxNumGvecs);
    const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
    //this will be parallelized with OpenMP
    for(int iorb=0; iorb<N; ++iorb)
    {
      int ti=cur_bands[iorb].TwistIndex;
      get_psi_g(ti,spin,cur_bands[iorb].BandIndex,cG);//bcast cG
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

#endif

  /** initialize the set to minimize the memroy use
   */
  bool initialize_spline_serial(int spin, const BandInfoGroup& bandgroup)
  {
    hdf_archive h5f(myComm,false);
    h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);

    const size_t N=bandgroup.getNumDistinctOrbitals();
    Vector<std::complex<double> > cG(mybuilder->MaxNumGvecs);
    const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
    bool foundit=true;
    for(size_t iorb=0; iorb<N; ++iorb)
    {
      int ti=cur_bands[iorb].TwistIndex;
      std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex);
      foundit &= h5f.read(cG,s);
      get_psi_g(ti,spin,cur_bands[iorb].BandIndex,cG);//bcast cG
      fft_spline(cG,ti,0);
      bspline->set_spline(spline_r[0],spline_i[0],cur_bands[iorb].TwistIndex,iorb,0);
    }
    return foundit;
  }

  void initialize_spline_pio_reduce(int spin, const BandInfoGroup& bandgroup)
  {
    int Nbands=bandgroup.getNumDistinctOrbitals();
    const int Nprocs=myComm->size();
    const int Nbandgroups=std::min(Nbands,Nprocs);
    Communicate band_group_comm(*myComm, Nbandgroups);
    std::vector<int> band_groups(Nbandgroups+1,0);
    FairDivideLow(Nbands,Nbandgroups,band_groups);
    int iorb_first=band_groups[band_group_comm.getGroupID()];
    int iorb_last =band_groups[band_group_comm.getGroupID()+1];
    if(band_group_comm.isGroupLeader())
    {
      hdf_archive h5f(&band_group_comm,false);
      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);
      Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());
      const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=cur_bands[iorb].TwistIndex;
        std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex);
        if(!h5f.read(cG,s)) APP_ABORT("SplineAdoptorReader Failed to read band(s) from h5!\n");
        fft_spline(cG,ti,ib);
        bspline->set_spline(spline_r[ib],spline_i[ib],cur_bands[iorb].TwistIndex,iorb,0);
      }
      chunked_reduce(band_group_comm.GroupLeaderComm, bspline->MultiSpline);
    }
    myComm->barrier();
    chunked_bcast(myComm, bspline->MultiSpline);
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
      hdf_archive h5f(myComm,false);
      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);
      Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());
      const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=cur_bands[iorb].TwistIndex;
        std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex);
        foundit &= h5f.read(cG,s);
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
      else if(root)
      {
        for(int iorb=OrbGroups[ip],ib=0; iorb<OrbGroups[ip+1]; ++iorb,++ib)
        {
          mpi::recv(*myComm,spline_r[ib]->coefs,ng_big,ip,iorb);
          if(bspline->is_complex)
            mpi::recv(*myComm,spline_i[ib]->coefs,ng_big,ip,iorb+N);
        }
        const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
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
    const std::vector<BandInfo>& cur_bands=bandgroup.myBands;

    if(myComm->rank()<np)
    {
      int iorb_first=OrbGroups[myComm->rank()];
      int iorb_last =OrbGroups[myComm->rank()+1];
      hdf_archive h5f(myComm,false);
      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);
      Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());
      for(int ib=0, iorb=iorb_first; iorb<iorb_last; ib++, iorb++)
      {
        int ti=cur_bands[iorb].TwistIndex;
        std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex);
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
    Array<std::complex<double>,3> rawData(nx,ny,nz);
    const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
    //this will be parallelized with OpenMP
    int N=bandgroup.getNumDistinctOrbitals();
    for(int iorb=0; iorb<N; ++iorb)
    {
      //check dimension
      if(myComm->rank() ==0)
      {
        std::string path=psi_r_path(cur_bands[iorb].TwistIndex,spin,cur_bands[iorb].BandIndex);
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(rawData);
        h_splineData.read(mybuilder->H5FileID, path.c_str());
        simd::copy(splineData_r.data(),splineData_i.data(),rawData.data(),rawData.size());
      }
      mpi::bcast(*myComm,splineData_r);
      if(bspline->is_complex)
        mpi::bcast(*myComm,splineData_i);
      bspline->set_spline(splineData_r.data(),splineData_i.data(),cur_bands[iorb].TwistIndex,iorb,0);
    }
  }
};
}
#endif
