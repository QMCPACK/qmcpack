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
    
    
/** @file SplineHybridAdoptorReader.h
 *
 * The most general reader class for SplineAdoptor using the full single grid for the supercell
 * - SplineR2RAdoptor
 * - SplineC2CPackedAdoptor
 * - SplineC2RPackedAdoptor
 * Each band is initialized with UBspline_3d_d and both real and imaginary parts are passed to the adoptor object
 * which will convert the data type to their internal precision. 
 */
#ifndef QMCPLUSPLUS_EINSPLINE_HYBRID_ADOPTOR_READERP_H
#define QMCPLUSPLUS_EINSPLINE_HYBRID_ADOPTOR_READERP_H
#include <mpi/collectives.h>
#include <mpi/point2point.h>
#include <Numerics/Quadrature.h>
#include <Numerics/Bessel.h>
#include <QMCWaveFunctions/BsplineFactory/HybridAdoptorBase.h>

//#include <QMCHamiltonians/Ylm.h>
//#define PRINT_RADIAL
//#define REPORT_MISMATCH

namespace qmcplusplus
{

template<typename ST, typename LT>
struct Gvectors
{
  typedef TinyVector<ST,3> PosType;
  typedef std::complex<ST> ValueType;

  const LT& Lattice;
  const std::vector<TinyVector<int,3> >& gvecs;
  ST                                     gmag_max;
  std::vector<PosType>                   gvecs_cart; //Cartesian.
  std::vector<ST>                        gmag;
  std::vector<aligned_vector<ST> >       YlmG;
  const size_t NumGvecs;
  int mylmax;

  Gvectors(const std::vector<TinyVector<int,3> >& gvecs_in, const LT& Lattice_in, const TinyVector<int,3>& HalfG):
  gvecs(gvecs_in), Lattice(Lattice_in), NumGvecs(gvecs.size()), gmag(0.0), mylmax(-1)
  {
    gvecs_cart.resize(NumGvecs);
    gmag.resize(NumGvecs);
    #pragma omp parallel for
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      TinyVector<ST,3> gvec_shift;
      gvec_shift=gvecs[ig]+HalfG*0.5;
      gvecs_cart[ig]=Lattice.k_cart(gvec_shift);
      gmag[ig]=std::sqrt(dot(gvecs_cart[ig],gvecs_cart[ig]));
      if (gmag[ig]>gmag_max) gmag_max=gmag[ig];
    }
  }

  void calc_YlmG(const int lmax)
  {
    mylmax=lmax;
    SoaSphericalTensor<ST> Ylm(lmax);
    const int lm_tot=(lmax+1)*(lmax+1);
    YlmG.resize(lm_tot);
    for(size_t lm=0; lm<lm_tot; lm++)
      YlmG[lm].resize(NumGvecs);
    #pragma omp parallel for
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      PosType Ghat;
      aligned_vector<ST> Ylm_vals(lm_tot);
      if (gmag[ig]==0)
        Ghat=PosType(0.0,0.0,1.0);
      else
        Ghat=gvecs_cart[ig]/gmag[ig];
      Ylm.evaluateV(Ghat[0], Ghat[1], Ghat[2], Ylm_vals.data());
      for(size_t lm=0; lm<lm_tot; lm++)
        YlmG[lm][ig]=Ylm_vals[lm];
    }
    //std::cout << "Calculated " << NumGvecs << " YlmG!" << std::endl;
  }

  void calc_jlm_G(const int lmax, ST& r, std::vector<aligned_vector<ST> >& j_lm_G)
  {
    if(lmax>mylmax)
    {
      app_error() << "Current Gvectors only has Ylm_G up to lmax = " << mylmax << " but asking for " << lmax << std::endl;
      abort();
    }

    // allocate space for j_lm_G[lm][ig]
    const int lm_tot=(lmax+1)*(lmax+1);
    j_lm_G.resize(lm_tot);
    for(size_t lm=0; lm<lm_tot; lm++)
      j_lm_G[lm].resize(NumGvecs);

    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      std::vector<double> jl_vals(lmax+1,0.0);
      bessel_steed_array_cpu(lmax, gmag[ig]*r, jl_vals.data());
      for(size_t l=0; l<=lmax; l++)
        for(size_t lm=l*l; lm<(l+1)*(l+1); lm++)
          j_lm_G[lm][ig]=jl_vals[l];
    }
  }

  template<typename PT>
  void calc_phase_shift(const PT& pos, aligned_vector<ValueType>& phase_shift)
  {
    phase_shift.resize(NumGvecs);
    #pragma omp parallel for
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      ST s,c;
      sincos(dot(gvecs_cart[ig],pos),&s,&c);
      phase_shift[ig]=ValueType(c,s);
    }
  }

  template<typename PT>
  ValueType evaluate_psi_r(const Vector<std::complex<double> >& cG, const PT& pos)
  {
    assert(cG.size()==NumGvecs);
    std::complex<ST> val(0.0,0.0);
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      ST s,c;
      sincos(dot(gvecs_cart[ig],pos),&s,&c);
      ValueType pw0(c,s);
      val+=cG[ig]*pw0;
    }
    return val;
  }

  template<typename PT>
  void evaluate_psi_r(const Vector<std::complex<double> >& cG, const PT& pos, ValueType& phi, ValueType& d2phi)
  {
    assert(cG.size()==NumGvecs);
    d2phi=phi=0.0;
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      ST s,c;
      sincos(dot(gvecs_cart[ig],pos),&s,&c);
      ValueType pw0(c,s);
      phi+=cG[ig]*pw0;
      d2phi+=cG[ig]*pw0*(-dot(gvecs_cart[ig],gvecs_cart[ig]));
    }
  }

  double evaluate_KE(const Vector<std::complex<double> >& cG)
  {
    assert(cG.size()==NumGvecs);
    double KE=0;
    for(size_t ig=0; ig<NumGvecs; ig++)
      KE+=dot(gvecs_cart[ig],gvecs_cart[ig])*(cG[ig].real()*cG[ig].real()+cG[ig].imag()*cG[ig].imag());
    return KE/2.0;
  }

};


/** General SplineHybridAdoptorReader to handle any unitcell
 */
template<typename SA>
struct SplineHybridAdoptorReader: public BsplineReaderBase
{
  typedef SA adoptor_type;
  typedef typename adoptor_type::DataType    DataType;
  typedef typename adoptor_type::SplineType SplineType;

  Array<std::complex<double>,3> FFTbox;
  Array<double,3> splineData_r, splineData_i;
  double rotate_phase_r, rotate_phase_i;
  UBspline_3d_d* spline_r;
  UBspline_3d_d* spline_i;
  BsplineSet<adoptor_type>* bspline;
  fftw_plan FFTplan;
#ifdef REPORT_MISMATCH
  std::vector<std::vector<double> > mismatch_energy_AO_to_PW;
#endif

  SplineHybridAdoptorReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e), spline_r(NULL), spline_i(NULL), bspline(0), FFTplan(NULL)
  {}

  ~SplineHybridAdoptorReader()
  {
    clear();
  }

  void clear()
  {
    einspline::destroy(spline_r);
    einspline::destroy(spline_i);
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

    // set info for Hybrid
    bspline->set_info(*(mybuilder->SourcePtcl), mybuilder->TargetPtcl, mybuilder->Super2Prim);
    initialize_atomic_centers(bspline->AtomicCenters);

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
      //chunked_bcast(myComm, bspline->MultiSpline);
      bspline->bcast_tables(myComm);
      app_log() << "  SplineHybridAdoptorReader bcast the full table " << now.elapsed() << " sec" << std::endl;
      app_log().flush();
    }
    else
    {
      bspline->flush_zero();

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

        TinyVector<double,3> start(0.0);
        TinyVector<double,3> end(1.0);
        spline_r=einspline::create(spline_r,start,end,MeshSize,bspline->HalfG);
        if(bspline->is_complex)
          spline_i=einspline::create(spline_i,start,end,MeshSize,bspline->HalfG);

        now.restart();
        //initialize_spline_pio(spin,bandgroup);
        initialize_spline_pio_reduce(spin,bandgroup);
        app_log() << "  SplineHybridAdoptorReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;

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
        app_log() << "  SplineHybridAdoptorReader dump " << now.elapsed() << " sec" << std::endl;
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
   * Perform FFT and spline to spline_r and spline_i
   */
  inline void fft_spline(Vector<std::complex<double> >& cG, int ti)
  {
    unpack4fftw(cG,mybuilder->Gvecs[0],MeshSize,FFTbox);
    fftw_execute (FFTplan);
    if(bspline->is_complex)
    {
      fix_phase_rotate_c2c(FFTbox,splineData_r, splineData_i,mybuilder->TwistAngles[ti], rotate_phase_r, rotate_phase_i);
      einspline::set(spline_r,splineData_r.data());
      einspline::set(spline_i,splineData_i.data());
    }
    else
    {
      fix_phase_rotate_c2r(FFTbox,splineData_r, mybuilder->TwistAngles[ti], rotate_phase_r, rotate_phase_i);
      einspline::set(spline_r,splineData_r.data());
    }
  }

  /** initialize basic parameters of atomic orbitals */
  void initialize_atomic_centers(std::vector<AtomicOrbitalSoA<DataType> >& centers)
  {
    const auto& ACInfo=mybuilder->AtomicCentersInfo;
    // load atomic center info only when it is not initialized
    if(centers.size()==0)
    {
      bool success=true;
      app_log() << "Reading atomic center info for hybrid representation" << std::endl;
      for(int center_idx=0; center_idx<ACInfo.Ncenters; center_idx++)
      {
        const int my_GroupID = ACInfo.GroupID[center_idx];
        if(ACInfo.lmax[center_idx]<0)
        {
          app_error() << "Hybrid representation needs parameter 'lmax' for atom " << center_idx << std::endl;
          success=false;
        }
        if(ACInfo.cutoff[center_idx]<0)
        {
          app_error() << "Hybrid representation needs parameter 'cutoff_radius' for atom " << center_idx << std::endl;
          success=false;
        }
        if(ACInfo.spline_radius[center_idx]<0)
        {
          app_error() << "Hybrid representation needs parameter 'spline_radius' for atom " << center_idx << std::endl;
          success=false;
        }
        if(ACInfo.spline_npoints[center_idx]<0)
        {
          app_error() << "Hybrid representation needs parameter 'spline_npoints' for atom " << center_idx << std::endl;
          success=false;
        }
        double max_allowed_cutoff=ACInfo.spline_radius[center_idx]-2.0*ACInfo.spline_radius[center_idx]/(ACInfo.spline_npoints[center_idx]-1);
        if(success && ACInfo.cutoff[center_idx]>max_allowed_cutoff)
        {
          app_error() << "Hybrid representation requires cutoff_radius<=" << max_allowed_cutoff
                      << " calculated by spline_radius-2*spline_radius/(spline_npoints-1) for atom " << center_idx << std::endl;
          success=false;
        }
      }
      if(!success) abort();

#ifdef REPORT_MISMATCH
      mismatch_energy_AO_to_PW.resize(ACInfo.Ncenters);
#endif
      for(int center_idx=0; center_idx<ACInfo.Ncenters; center_idx++)
      {
        AtomicOrbitalSoA<DataType> oneCenter(ACInfo.lmax[center_idx]);
        oneCenter.set_info(ACInfo.ion_pos[center_idx], ACInfo.cutoff[center_idx],
                           ACInfo.spline_radius[center_idx], ACInfo.spline_npoints[center_idx]);
        centers.push_back(oneCenter);
      }
    }

#ifdef REPORT_MISMATCH
    for(int center_idx=0; center_idx<ACInfo.Ncenters; center_idx++)
      mismatch_energy_AO_to_PW[center_idx].resize(std::max(ACInfo.spline_npoints[center_idx]-3,0),0.0);
#endif
  }

  /** initialize construct atomic orbital radial functions from plane waves */
  inline void create_atomic_centers_Gspace(Vector<std::complex<double> >& cG, Communicate& band_group_comm, int iorb)
  {
    std::vector<AtomicOrbitalSoA<DataType> >& centers=bspline->AtomicCenters;
    //distribute atomic centers over processor groups
    int Ncenters=centers.size();
    const int Nprocs=band_group_comm.size();
    const int Ncentergroups=std::min(Ncenters,Nprocs);
    Communicate center_group_comm(band_group_comm, Ncentergroups);
    std::vector<int> center_groups(Ncentergroups+1,0);
    FairDivideLow(Ncenters,Ncentergroups,center_groups);
    int center_first=center_groups[center_group_comm.getGroupID()];
    int center_last =center_groups[center_group_comm.getGroupID()+1];

    // no workload case
    if(!center_group_comm.isGroupLeader()) return;

    typedef typename EinsplineSetBuilder::UnitCellType UnitCellType;

    // prepare Gvecs Ylm(G)
    Gvectors<double, UnitCellType> Gvecs(mybuilder->Gvecs[0], mybuilder->PrimCell, bspline->HalfG);
    // if(band_group_comm.isGroupLeader()) std::cout << "print band=" << iorb << " KE=" << Gvecs.evaluate_KE(cG) << std::endl;
    const int lmax_limit=7;
    Gvecs.calc_YlmG(lmax_limit);
    std::vector<std::complex<double> > i_power;
    // rotate phase is introduced here.
    std::complex<double> i_temp(rotate_phase_r, rotate_phase_i);
    for(size_t l=0; l<=lmax_limit; l++)
    {
      for(size_t lm=l*l; lm<(l+1)*(l+1); lm++)
        i_power.push_back(i_temp);
      i_temp*=std::complex<double>(0.0,1.0);
    }

    app_log() << "Transforming band " << iorb << " centers from " << center_first << " to " << center_last-1 << " on Rank 0" << std::endl;
    for(int center_idx=center_first; center_idx<center_last; center_idx++)
    {
      AtomicOrbitalSoA<DataType>& mycenter=centers[center_idx];
      const double cutoff_radius = mycenter.cutoff;
      const double spline_radius = mycenter.spline_radius;
      const int   spline_npoints = mycenter.spline_npoints;
      const int             lmax = mycenter.lmax;
      const double delta = spline_radius/static_cast<double>(spline_npoints-1);
      const int lm_tot=(lmax+1)*(lmax+1);

#ifdef PRINT_RADIAL
      char fname[64];
      sprintf(fname, "band_%d_center_%d_pw.dat", iorb, center_idx);
      FILE *fout_pw  = fopen (fname, "w");
      sprintf(fname, "band_%d_center_%d_spline_v.dat", iorb, center_idx);
      FILE *fout_spline_v  = fopen (fname, "w");
      sprintf(fname, "band_%d_center_%d_spline_g.dat", iorb, center_idx);
      FILE *fout_spline_g  = fopen (fname, "w");
      sprintf(fname, "band_%d_center_%d_spline_l.dat", iorb, center_idx);
      FILE *fout_spline_l  = fopen (fname, "w");
      fprintf(fout_pw, "# r vals(lm)\n");
      fprintf(fout_spline_v, "# r vals(lm)\n");
      fprintf(fout_spline_g, "# r grads(lm)\n");
      fprintf(fout_spline_l, "# r lapls(lm)\n");
#endif

      std::vector<std::vector<std::complex<double> > > all_vals;
      all_vals.resize(spline_npoints);
      aligned_vector<std::complex<double> > phase_shift;
      Gvecs.calc_phase_shift(mycenter.pos, phase_shift);

      #pragma omp parallel for
      for(int ip=0; ip<spline_npoints; ip++)
      {
        double r=delta*static_cast<double>(ip);
        std::vector<std::complex<double> > vals(lm_tot, std::complex<double>(0.0,0.0));
        std::vector<aligned_vector<double> > j_lm_G;
        Gvecs.calc_jlm_G(lmax, r, j_lm_G);

        for(int lm=0; lm<lm_tot; lm++)
          for(size_t ig=0; ig<Gvecs.NumGvecs; ig++)
          {
            vals[lm]+=cG[ig]*phase_shift[ig]*j_lm_G[lm][ig]*Gvecs.YlmG[lm][ig];
          }

        all_vals[ip]=vals;
        for(int lm=0; lm<lm_tot; lm++)
          all_vals[ip][lm]*=4.0*M_PI*i_power[lm];
      }
      //app_log() << "Building band " << iorb << " at center " << center_idx << std::endl;

      #pragma omp parallel for
      for(int lm=0; lm<lm_tot; lm++)
      {
        aligned_vector<double> splineData_r(spline_npoints);
        UBspline_1d_d* atomic_spline_r;
        for(size_t ip=0; ip<spline_npoints; ip++)
          splineData_r[ip]=real(all_vals[ip][lm]);
        atomic_spline_r=einspline::create(atomic_spline_r, 0.0, spline_radius, spline_npoints, splineData_r.data(), ((lm==0)||(lm>3)));
        if(!bspline->is_complex)
        {
          mycenter.set_spline(atomic_spline_r,lm,iorb);
          einspline::destroy(atomic_spline_r);
        }
        else
        {
          aligned_vector<double> splineData_i(spline_npoints);
          UBspline_1d_d* atomic_spline_i;
          for(size_t ip=0; ip<spline_npoints; ip++)
            splineData_i[ip]=imag(all_vals[ip][lm]);
          atomic_spline_i=einspline::create(atomic_spline_i, 0.0, spline_radius, spline_npoints, splineData_i.data(), ((lm==0)||(lm>3)));
          int iband=bspline->BandIndexMap.size()>0?bspline->BandIndexMap[iorb]:iorb;
          mycenter.set_spline(atomic_spline_r,lm,iband*2);
          mycenter.set_spline(atomic_spline_i,lm,iband*2+1);
          einspline::destroy(atomic_spline_r);
          einspline::destroy(atomic_spline_i);
        }
      }

#ifdef PRINT_RADIAL
      // write to file for plotting
      for(int ip=0; ip<spline_npoints-1; ip++)
      {
        double r=delta*static_cast<double>(ip);
        mycenter.SplineInst->evaluate_vgl(r,mycenter.localV,mycenter.localG,mycenter.localL);
        fprintf(fout_pw, "%15.10lf  ", r);
        fprintf(fout_spline_v, "%15.10lf  ", r);
        fprintf(fout_spline_g, "%15.10lf  ", r);
        fprintf(fout_spline_l, "%15.10lf  ", r);
        for(int lm=0; lm<lm_tot; lm++)
        {
          fprintf(fout_pw, "%15.10lf  %15.10lf  ", all_vals[ip][lm].real(), all_vals[ip][lm].imag());
          fprintf(fout_spline_v, "%15.10lf  %15.10lf  ", mycenter.localV[lm*mycenter.Npad+iorb*2], mycenter.localV[lm*mycenter.Npad+iorb*2+1]);
          fprintf(fout_spline_g, "%15.10lf  %15.10lf  ", mycenter.localG[lm*mycenter.Npad+iorb*2], mycenter.localG[lm*mycenter.Npad+iorb*2+1]);
          fprintf(fout_spline_l, "%15.10lf  %15.10lf  ", mycenter.localL[lm*mycenter.Npad+iorb*2], mycenter.localL[lm*mycenter.Npad+iorb*2+1]);
        }
        fprintf(fout_pw, "\n");
        fprintf(fout_spline_v, "\n");
        fprintf(fout_spline_g, "\n");
        fprintf(fout_spline_l, "\n");
      }
#endif

#ifdef PRINT_RADIAL
      fclose(fout_pw);
      fclose(fout_spline_v);
      fclose(fout_spline_g);
      fclose(fout_spline_l);
#endif

#ifdef REPORT_MISMATCH
      const int total_grid_size=mismatch_energy_AO_to_PW[center_idx].size();
      Quadrature3D<double> quad(7);
      #pragma omp parallel for
      for(int ip=0; ip<total_grid_size; ip++)
      {
        double r=delta*static_cast<double>(ip+1);
        double AO_to_PW(0.0);
        int quad_counter=0;
        // sum up quadratures
        for(int j=0; j<quad.nk; j++)
        {
          typedef TinyVector<DataType,3> PosType;
          PosType pos_quad=quad.xyz_m[j]*r+mycenter.pos;
          double PW_ene, AO_ene, PW_norm2, AO_norm2;

          // Evaluate PW.
          std::complex<double> PW_psi, PW_d2psi;
          Gvecs.evaluate_psi_r(cG, pos_quad, PW_psi, PW_d2psi);
          PW_psi*=std::complex<double>(rotate_phase_r, rotate_phase_i);
          PW_d2psi*=std::complex<double>(rotate_phase_r, rotate_phase_i);
          if(bspline->is_complex)
          {
            PW_ene  =PW_psi.real()*PW_d2psi.real()+PW_psi.imag()*PW_d2psi.imag();
            PW_norm2=PW_psi.real()*  PW_psi.real()+PW_psi.imag()*  PW_psi.imag();
          }
          else
          {
            PW_ene  =PW_psi.real()*PW_d2psi.real();
            PW_norm2=PW_psi.real()*  PW_psi.real();
          }

          // Evaluate AO. This is a lazy implementation.
          // The evaluation goes over all the bands althought the results of only one band is used.
          auto& myV=bspline->myV;
          auto& myG=bspline->myG;
          auto& myL=bspline->myL;
          #pragma omp critical
          mycenter.evaluate_vgl(r, quad.xyz_m[j]*r, myV, myG, myL);
          if(bspline->is_complex)
          {
            AO_ene  =myV[iorb*2]*myL[iorb*2]+myV[iorb*2+1]*myL[iorb*2+1];
            AO_norm2=myV[iorb*2]*myV[iorb*2]+myV[iorb*2+1]*myV[iorb*2+1];
          }
          else
          {
            AO_ene  =myV[iorb]*myL[iorb];
            AO_norm2=myV[iorb]*myV[iorb];
          }

          //collect mismatch energy
          if(AO_norm2>std::numeric_limits<float>::epsilon() && PW_norm2>std::numeric_limits<float>::epsilon())
          {
            quad_counter++;
            AO_to_PW+=std::fabs(AO_ene-PW_ene)/2.0;
          }

          //std::cout << "debug PW KE " << PW_ene << std::endl;
          //std::cout << "debug AO KE " << AO_ene << std::endl;
        }
        //std::cout << "debug iorb=" << iorb << " center=" << center_idx << " ip=" << ip << " r=" << r << " mismatch_energy_AO_to_PW=" << AO_to_PW/(quad_counter>0?quad_counter:1) << std::endl;
        mismatch_energy_AO_to_PW[center_idx][ip]+=AO_to_PW/(quad_counter>0?quad_counter:1);
      }
#endif
    }
    // collect atomic orbitals from all the centers
    bspline->reduce_atomic_tables(center_group_comm.GroupLeaderComm);
#ifdef REPORT_MISMATCH
    // collect mismatch_energy
    for(int center_idx=0; center_idx<Ncenters; center_idx++)
      mpi::reduce(*center_group_comm.GroupLeaderComm,mismatch_energy_AO_to_PW[center_idx]);
#endif
  }


  /** initialize the splines
   */
  void initialize_spline_pio_reduce(int spin, const BandInfoGroup& bandgroup)
  {
    //distribute bands over processor groups
    int Nbands=bandgroup.getNumDistinctOrbitals();
    const int Nprocs=myComm->size();
    const int Nbandgroups=std::min(Nbands,Nprocs);
    Communicate band_group_comm(*myComm, Nbandgroups);
    std::vector<int> band_groups(Nbandgroups+1,0);
    FairDivideLow(Nbands,Nbandgroups,band_groups);
    int iorb_first=band_groups[band_group_comm.getGroupID()];
    int iorb_last =band_groups[band_group_comm.getGroupID()+1];

    app_log() << "Start transforming 3D B-Splines and atomic orbitals for hybrid representation." << std::endl;
    hdf_archive h5f(&band_group_comm,false);
    Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());
    const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
    if(band_group_comm.isGroupLeader())
      h5f.open(mybuilder->H5FileName,H5F_ACC_RDONLY);
    for(int iorb=iorb_first; iorb<iorb_last; iorb++)
    {
      if(band_group_comm.isGroupLeader())
      {
        int ti=cur_bands[iorb].TwistIndex;
        std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex);
        if(!h5f.read(cG,s)) APP_ABORT("SplineHybridAdoptorReader Failed to read band(s) from h5!\n");
        fft_spline(cG,ti);
        bspline->set_spline(spline_r,spline_i,cur_bands[iorb].TwistIndex,iorb,0);
      }
      band_group_comm.bcast(rotate_phase_r);
      band_group_comm.bcast(rotate_phase_i);
      band_group_comm.bcast(cG);
      create_atomic_centers_Gspace(cG, band_group_comm, iorb);
    }

    myComm->barrier();
    Timer now;
    now.restart();
    if(band_group_comm.isGroupLeader())
    {
      bspline->reduce_tables(band_group_comm.GroupLeaderComm);
#ifdef REPORT_MISMATCH
      // collect mismatch_energy
      for(int center_idx=0; center_idx<bspline->AtomicCenters.size(); center_idx++)
        mpi::reduce(*band_group_comm.GroupLeaderComm,mismatch_energy_AO_to_PW[center_idx]);
#endif
    }
    bspline->bcast_tables(myComm);
    app_log() << "  Time to reduce and bcast the table = " << now.elapsed() << std::endl;

#ifdef REPORT_MISMATCH
    if(!myComm->rank()) print_mismatch_scan(spin);
#endif
  }


#ifdef REPORT_MISMATCH
  /** print the mismatch energy between plane waves and AO */
  void print_mismatch_scan(const int spin)
  {
    const auto& mySpecies=mybuilder->SourcePtcl->mySpecies;
    int Nspecies=mySpecies.size();
    std::vector<int> atom_counter(Nspecies,0),lmax_per_species(Nspecies,0.0);
    std::vector<double> delta_per_species(Nspecies,0.0);
    std::vector<std::vector<double> > mismatch_ene_species_AO_to_PW(Nspecies);
    const auto& AtomicCenters=bspline->AtomicCenters;
    // set proper size
    for(int center_idx=0; center_idx<AtomicCenters.size(); center_idx++)
    {
      const int my_GroupID = mybuilder->AtomicCentersInfo.GroupID[center_idx];
      if(mismatch_energy_AO_to_PW[center_idx].size()>mismatch_ene_species_AO_to_PW[my_GroupID].size())
      {
        mismatch_ene_species_AO_to_PW[my_GroupID].resize(mismatch_energy_AO_to_PW[center_idx].size(),0.0);
        delta_per_species[my_GroupID]=AtomicCenters[center_idx].spline_radius/static_cast<double>(AtomicCenters[center_idx].spline_npoints-1);
        lmax_per_species[my_GroupID]=AtomicCenters[center_idx].lmax;
      }
    }
    // average per species
    for(int center_idx=0; center_idx<AtomicCenters.size(); center_idx++)
    {
      //app_log() << " Yedebug center=" << center_idx << std::endl;
      //for(int j=0; j<mismatch_energy_AO_to_PW[center_idx].size(); j++)
      //  app_log() << " ip=" << j << " dis=" << (j+1)*AtomicCenters[center_idx].spline_radius/static_cast<double>(AtomicCenters[center_idx].spline_npoints-1)
      //            << " mismatch=" << mismatch_energy_AO_to_PW[center_idx][j] << std::endl;
      const int my_GroupID = mybuilder->AtomicCentersInfo.GroupID[center_idx];
      for(int ip=0; ip<mismatch_energy_AO_to_PW[center_idx].size(); ip++)
        mismatch_ene_species_AO_to_PW[my_GroupID][ip]+=mismatch_energy_AO_to_PW[center_idx][ip];
      atom_counter[my_GroupID]++;
    }
    // print per species
    for(int ispecies=0; ispecies<Nspecies; ispecies++)
    {
      std::ostringstream filename;
      filename << "AO_to_PW.spin_" << spin << ".tw" << mybuilder->TwistNum << "." << mySpecies.speciesName[ispecies]
               << ".lmax" << lmax_per_species[ispecies] << ".dat";
      std::ofstream outfile(filename.str().c_str());
      outfile << "# ip r abs(Ediff)" << std::endl;
      for(int ip=0; ip<mismatch_ene_species_AO_to_PW[ispecies].size(); ip++)
        outfile << ip << "  " << std::fixed << std::setprecision(4) << (ip+1)*delta_per_species[ispecies] << "  "
                << std::setprecision(9)
                << mismatch_ene_species_AO_to_PW[ispecies][ip]/(atom_counter[ispecies]>0?atom_counter[ispecies]:1)
                << std::endl;
      outfile.close();
    }
  }
#endif


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
