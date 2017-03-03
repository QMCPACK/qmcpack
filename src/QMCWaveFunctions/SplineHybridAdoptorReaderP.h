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
#define PRINT_RADIAL

namespace qmcplusplus
{

template<typename ST, typename LT>
struct Gvectors
{
  typedef TinyVector<ST,3> PosType;

  const LT& Lattice;
  const std::vector<TinyVector<int,3> >& gvecs;
  ST                                     gmag_max;
  std::vector<PosType>                   gvecs_cart; //Cartesian.
  std::vector<ST>                        gmag;
  std::vector<aligned_vector<ST> >       YlmG;
  const size_t NumGvecs;
  int mylmax;

  Gvectors(const std::vector<TinyVector<int,3> >& gvecs_in, const LT& Lattice_in):
  gvecs(gvecs_in), Lattice(Lattice_in), NumGvecs(gvecs.size()), gmag(0.0), mylmax(-1)
  {
    gvecs_cart.resize(NumGvecs);
    gmag.resize(NumGvecs);
    #pragma omp parallel for
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      gvecs_cart[ig]=Lattice.k_cart(gvecs[ig]);
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
  void calc_phase_shift(const PT& pos, aligned_vector<std::complex<ST> >& phase_shift)
  {
    phase_shift.resize(NumGvecs);
    #pragma omp parallel for
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      ST s,c;
      sincos(dot(gvecs_cart[ig],pos),&s,&c);
      phase_shift[ig]=std::complex<ST>(c,s);
    }
  }

  template<typename PT>
  std::complex<ST> evaluate_psi_r(const Vector<std::complex<double> >& cG, const PT& pos)
  {
    std::complex<ST> val(0.0,0.0);
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      ST s,c;
      sincos(dot(gvecs_cart[ig],pos),&s,&c);
      std::complex<ST> pw0(c,s);
      val+=cG[ig]*pw0;
    }
    return val;
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
  std::vector<UBspline_3d_d*> spline_r;
  std::vector<UBspline_3d_d*> spline_i;
  BsplineSet<adoptor_type>* bspline;
  std::vector<int> OrbGroups;
  fftw_plan FFTplan;

  SplineHybridAdoptorReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e), bspline(0),FFTplan(NULL)
  {}

  ~SplineHybridAdoptorReader()
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

    // set info for Hybrid
    bspline->set_info(*(mybuilder->SourcePtcl), mybuilder->TargetPtcl, mybuilder->Super2Prim);

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
        if(foundspline) app_log() << "  Time to read the table in " << splinefile << " = " << now.elapsed() << std::endl;;
      }
    }
    myComm->bcast(foundspline);
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << std::endl;
      now.restart();
      chunked_bcast(myComm, bspline->MultiSpline);
      app_log() << "  SplineHybridAdoptorReader bcast the full table " << now.elapsed() << " sec" << std::endl;
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

        initialize_atomic_centers(bspline->AtomicCenters);

        if(usingSerialIO)
        {
          now.restart();
          initialize_spline_serial(spin, bandgroup);
          app_log() << "  SplineHybridAdoptorReader initialize_spline_serial " << now.elapsed() << " sec" << std::endl;
        }
        else
        {
          //now.restart();
          //initialize_spline_pio_bcast(spin);
          //app_log() << "  SplineHybridAdoptorReader initialize_spline_pio_bcast " << now.elapsed() << " sec" << std::endl;
          size_t ntot=size_t(nx*ny*nz)*size_t(N);
          //if(ntot>>22) //Using 4M as the cutoff, candidate for autotuning
          {
            now.restart();
            initialize_spline_pio(spin,bandgroup);
            app_log() << "  SplineHybridAdoptorReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;
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

  inline int checkout_parameter_index(SpeciesSet& mySpecies, const std::string& name)
  {
    const int index=mySpecies.findAttribute(name);
    if(index<0)
    {
      app_error() << "Parameter " << name << " needed for each atom species to use hybrid representation." << std::endl;
      abort();
    }
    return index;
  }

  void initialize_atomic_centers(std::vector<AtomicOrbitalSoA<DataType> >& centers)
  {
    // load atomic center info only when it is not initialized
    if(centers.size()==0)
    {
      ParticleSet& PrimSourcePtcl=mybuilder->PrimSourcePtcl;
      SpeciesSet& mySpecies=PrimSourcePtcl.mySpecies;

      int cutoff_radius_ind=checkout_parameter_index(mySpecies,"cutoff_radius");
      int spline_radius_ind=checkout_parameter_index(mySpecies,"spline_radius");
      int spline_npoints_ind=checkout_parameter_index(mySpecies,"spline_npoints");
      int lmax_ind=checkout_parameter_index(mySpecies,"lmax");

      app_log() << "Reading atomic center info for hybrid representation" << std::endl;
      for(int center_idx=0; center_idx<PrimSourcePtcl.R.size(); center_idx++)
      {
        const int my_GroupID = PrimSourcePtcl.GroupID[center_idx];
        double cutoff_radius = mySpecies(cutoff_radius_ind, my_GroupID);
        double spline_radius = mySpecies(spline_radius_ind, my_GroupID);
        int   spline_npoints = mySpecies(spline_npoints_ind, my_GroupID);
        int             lmax = mySpecies(lmax_ind, my_GroupID);
        const int npad=bspline->myV.size();
        AtomicOrbitalSoA<DataType> oneCenter(lmax,npad);
        oneCenter.set_info(PrimSourcePtcl.R[center_idx], cutoff_radius, spline_radius, spline_npoints);
        centers.push_back(oneCenter);
      }
    }

    // allocate the memory hosting the spline coefficients.
    for(int center_idx=0; center_idx<centers.size(); center_idx++)
      centers[center_idx].create_spline();
  }

  inline void create_atomic_centers_Gspace(Vector<std::complex<double> >& cG, int ti, int iorb)
  {
    typedef typename EinsplineSetBuilder::UnitCellType UnitCellType;
    //Quadrature3D<double> quad(5);

    // prepare Gvecs Ylm(G)
    Gvectors<double, UnitCellType> Gvecs(mybuilder->Gvecs[0], mybuilder->PrimSourcePtcl.Lattice);
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

#if 0
    //checking Ylm
    const int lm_limit_tot=(lmax_limit+1)*(lmax_limit+1);
    SoaSphericalTensor<double> myYlm(lmax_limit,true);
    for(int j=0; j<quad.nk; j++)
    {
      std::vector<double> Ylm_vals(lm_limit_tot);
      myYlm.evaluateV(quad.xyz_m[j][0], quad.xyz_m[j][1], quad.xyz_m[j][2], Ylm_vals.data());
      TinyVector<double,3> rotate(quad.xyz_m[j][2], quad.xyz_m[j][0], quad.xyz_m[j][1]);
      app_log() << "*****************************" << std::endl;
      app_log() << "checking quad " << quad.xyz_m[j] << std::endl;
      for(int l=0; l<=lmax_limit; l++)
      {
        app_log() << "l,m " << l << ",0 " << Ylm_vals[l*(l+1)] << " " << real(Ylm(l,0,rotate)) << std::endl;
        for(int m=1; m<=l; m++)
        {
          if(std::abs(Ylm_vals[l*(l+1)-m]-std::sqrt(2)*(m%2?-1:1)*imag(Ylm(l,m,rotate)))>1e-10) 
            app_log() << "l,m " << l << "," << -m << " " << Ylm_vals[l*(l+1)-m] << " " << std::sqrt(2)*(m%2?-1:1)*imag(Ylm(l,m,rotate)) << std::endl;
          if(std::abs(Ylm_vals[l*(l+1)+m]-std::sqrt(2)*(m%2?-1:1)*real(Ylm(l,m,rotate)))>1e-10) 
            app_log() << "l,m " << l << "," << m << " " << Ylm_vals[l*(l+1)+m] << " " << std::sqrt(2)*(m%2?-1:1)*real(Ylm(l,m,rotate)) << std::endl;
        }
      }
    }
    abort();
#endif

    std::vector<AtomicOrbitalSoA<DataType> >& centers=bspline->AtomicCenters;
    for(int center_idx=0; center_idx<centers.size(); center_idx++)
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
      app_log() << "Building band " << iorb << " at center " << center_idx << std::endl;
#if 0
      app_log() << "checking band " << iorb << " center " << center_idx << " at " << mycenter.pos << std::endl;
      SoaSphericalTensor<double> Ylm(lmax);
      for(int j=1; j<quad.nk; j++)
      {
        std::vector<double> Ylm_vals(lm_tot);
        Ylm.evaluateV(quad.xyz_m[j][0], quad.xyz_m[j][1], quad.xyz_m[j][2], Ylm_vals.data());
        app_log() << "checking quad " << quad.xyz_m[j] << std::endl;
        //print out error in each direction
        for(int ip=0; ip<spline_npoints; ip++)
        {
          double r=delta*static_cast<double>(ip);
          std::complex<double> psi_ref=Gvecs.evaluate_psi_r(cG,quad.xyz_m[j]*r+mycenter.pos);
          std::complex<double> psi_sum(0.0,0.0);
          app_log() << " quad " << j << " r " << r << "  " << real(psi_ref) << "  " << imag(psi_ref);
          for(int l=0; l<=lmax; l++)
          {
            for(int lm=l*l; lm<(l+1)*(l+1); lm++)
              psi_sum+=all_vals[ip][lm]*Ylm_vals[lm];
            app_log() << "  " << real(psi_sum-psi_ref) << "  " << imag(psi_sum-psi_ref);
          }
          app_log() << std::endl;
        }
      }
      if( iorb==1 && center_idx==0 ) abort();
#endif

      #pragma omp parallel for
      for(int lm=0; lm<lm_tot; lm++)
      {
        aligned_vector<double> splineData_r(spline_npoints);
        UBspline_1d_d* spline_r;
        for(size_t ip=0; ip<spline_npoints; ip++)
          splineData_r[ip]=real(all_vals[ip][lm]);
        spline_r=einspline::create(spline_r, 0.0, spline_radius, spline_npoints, ((lm==0)||(lm>3)));
        einspline::set(spline_r,splineData_r.data());
        if(!bspline->is_complex)
        {
          mycenter.set_spline(spline_r,lm,iorb);
          einspline::destroy(spline_r);
        }
        else
        {
          aligned_vector<double> splineData_i(spline_npoints);
          UBspline_1d_d* spline_i;
          for(size_t ip=0; ip<spline_npoints; ip++)
            splineData_i[ip]=imag(all_vals[ip][lm]);
          spline_i=einspline::create(spline_i, 0.0, spline_radius, spline_npoints, ((lm==0)||(lm>3)));
          einspline::set(spline_i,splineData_i.data());
          int iband=bspline->BandIndexMap.size()>0?bspline->BandIndexMap[iorb]:iorb;
          mycenter.set_spline(spline_r,lm,iband*2);
          mycenter.set_spline(spline_i,lm,iband*2+1);
          einspline::destroy(spline_r);
          einspline::destroy(spline_i);
        }
      }

#ifdef PRINT_RADIAL
      // write to file for plotting
      for(int ip=0; ip<spline_npoints; ip++)
      {
        double r=delta*static_cast<double>(ip);
        fprintf(fout_pw, "%15.10lf  ", r);
        fprintf(fout_spline_v, "%15.10lf  ", r);
        fprintf(fout_spline_g, "%15.10lf  ", r);
        fprintf(fout_spline_l, "%15.10lf  ", r);
        //einspline::evaluate(mycenter.MultiSpline,r,mycenter.localV);
        einspline::evaluate(mycenter.MultiSpline,r,mycenter.localV,mycenter.localG,mycenter.localL);
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

#if 0
      // check u(r)
      SoaSphericalTensor<double> Ylm(lmax);
      std::vector<double> Ylm_vals(lm_tot);
      Ylm.evaluateV(0.0,1.0,0.0, Ylm_vals.data());
      std::cout << "# iorb " << iorb << std::endl;
      for(int ip=0; ip<spline_npoints; ip++)
      {
        double r=delta*static_cast<double>(ip);
        std::complex<double> psi_ref(0.0,0.0);
        for(int lm=0; lm<lm_tot; lm++)
          psi_ref+=all_vals[ip][lm]*Ylm_vals[lm];
        std::cout << "x_y_z " << 10.0 << " " << 10.0+r << " " << 10.0 << " : " << psi_ref << std::endl;
      }
#endif
#if 0
      SoaSphericalTensor<double> Ylm(lmax);
      std::vector<double> Ylm_vals(lm_tot);
      Ylm.evaluateV(0.0,1.0,0.0, Ylm_vals.data());
      double mydelta=2.0/100;
      std::cout << "# iorb " << iorb << std::endl;
      for(int ip=0; ip<121; ip++)
      {
        double r=mydelta*ip;
        TinyVector<double,3> mypos(10.0,10.0+r,10.0);
        std::complex<double> psi_ref_0(0.0,0.0), psi_ref_1(0.0,0.0);
        //psi_ref=Gvecs.evaluate_psi_r(cG,mypos)*i_power[0];
        einspline::evaluate(mycenter.MultiSpline,r,mycenter.localV);
        for(int lm=0; lm<lm_tot; lm++)
        {
          //psi_ref+=std::complex<double>(mycenter.localV[lm*mycenter.Npad+iorb*2]*Ylm_vals[lm],mycenter.localV[lm*mycenter.Npad+iorb*2+1]*Ylm_vals[lm]);
          psi_ref_0+=std::complex<double>(mycenter.localV[lm*mycenter.Npad+0]*Ylm_vals[lm],mycenter.localV[lm*mycenter.Npad+1]*Ylm_vals[lm]);
          psi_ref_1+=std::complex<double>(mycenter.localV[lm*mycenter.Npad+2]*Ylm_vals[lm],mycenter.localV[lm*mycenter.Npad+3]*Ylm_vals[lm]);
          std::cout << "debug lm " << Ylm_vals[lm] << " " << mycenter.localV[lm*mycenter.Npad+0] << " " << mycenter.localV[lm*mycenter.Npad+1]
                                                   << " " << mycenter.localV[lm*mycenter.Npad+2] << " " << mycenter.localV[lm*mycenter.Npad+3] << std::endl;
        }
        std::cout << "x_y_z " << 10.0 << " " << 10.0+r << " " << 10.0 << " : " << psi_ref_0 << psi_ref_1 << std::endl;
      }
#endif

#ifdef PRINT_RADIAL
      fclose(fout_pw);
      fclose(fout_spline_v);
      fclose(fout_spline_g);
      fclose(fout_spline_l);
#endif
    }
  }

  inline void create_atomic_centers_Rspace(Vector<std::complex<double> >& cG, int ti, int iorb)
  {
    typedef typename EinsplineSetBuilder::UnitCellType UnitCellType;

    Quadrature3D<double> quad(6);
    Gvectors<double, UnitCellType> Gvecs(mybuilder->Gvecs[0], mybuilder->PrimSourcePtcl.Lattice);

    std::vector<AtomicOrbitalSoA<DataType> >& centers=bspline->AtomicCenters;
    //#pragma omp parallel for
    for(int center_idx=0; center_idx<centers.size(); center_idx++)
    {
      AtomicOrbitalSoA<DataType>& mycenter=centers[center_idx];
      const double cutoff_radius = mycenter.cutoff;
      const double spline_radius = mycenter.spline_radius;
      const int   spline_npoints = mycenter.spline_npoints;
      const int             lmax = mycenter.lmax;
      double delta = spline_radius/static_cast<double>(spline_npoints-1);

      char fname[20];
      sprintf(fname, "band_%d_center_%d", iorb, center_idx);
      FILE *fout_pw  = fopen (fname, "w");
      fprintf(fout_pw, "# r vals(lm)\n");
      std::vector<std::vector<std::complex<double> > > all_vals;
      all_vals.resize(spline_npoints);

      SoaSphericalTensor<double> Ylm(lmax);
      Array<double,2> Ylm_vals(quad.nk,Ylm.size());
      for(int j=0; j<quad.nk; j++)
        Ylm.evaluateV(quad.xyz_m[j][0], quad.xyz_m[j][1], quad.xyz_m[j][2], &Ylm_vals(j,0));
      // vector splineData_r/i for each grid point
      #pragma omp parallel for
      for(int ip=0; ip<spline_npoints; ip++)
      {
        double r=delta*static_cast<double>(ip);
        std::vector<std::complex<double> > vals(Ylm.size(), std::complex<double>(0.0,0.0));
        // sum up quadratures
        for(int j=0; j<quad.nk; j++)
        {
          std::complex<double> psi=Gvecs.evaluate_psi_r(cG, quad.xyz_m[j]*r+mycenter.pos);
          for(int lm=0; lm<Ylm.size(); lm++)
          {
            vals[lm]+=Ylm_vals(j,lm)*psi*quad.weight_m[j];
          }
        }
        all_vals[ip]=vals;
        for(int lm=0; lm<Ylm.size(); lm++)
          all_vals[ip][lm]*=4.0*M_PI;
      }
      app_log() << "debug band " << iorb << " center " << center_idx << std::endl;

#if 0
      app_log() << "checking band " << iorb << " center " << center_idx << " at " << mycenter.pos << std::endl;
      for(int j=1; j<quad.nk; j++)
      {
        //print out error in each direction
        for(int ip=0; ip<spline_npoints; ip++)
        {
          double r=delta*static_cast<double>(ip);
          std::complex<double> psi_ref=Gvecs.evaluate_psi_r(cG,quad.xyz_m[j]*r+mycenter.pos);
          std::complex<double> psi_sum(0.0,0.0);
          app_log() << " quad " << j << " r " << r << "  " << real(psi_ref) << "  " << imag(psi_ref);
          for(int l=0; l<=lmax; l++)
          {
            for(int lm=l*l; lm<(l+1)*(l+1); lm++)
              psi_sum+=all_vals[ip][lm]*Ylm_vals(j,lm);
            app_log() << "  " << real(psi_sum-psi_ref) << "  " << imag(psi_sum-psi_ref);
          }
          app_log() << std::endl;
        }
        if( iorb==0 && center_idx==1 ) abort();
      }
#endif

      for(int ip=0; ip<spline_npoints; ip++)
      {
        fprintf(fout_pw, "%15.10lf  ", delta*static_cast<double>(ip));
        for(int lm=0; lm<Ylm.size(); lm++)
          fprintf(fout_pw, "%15.10lf  %15.10lf  ", all_vals[ip][lm].real(), all_vals[ip][lm].imag());
        fprintf(fout_pw, "\n");
      }
      // fill it in the big table N bands
      // push into class.
      fclose(fout_pw);
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
      create_atomic_centers_Gspace(cG,ti,iorb);
      bspline->set_spline(spline_r[0],spline_i[0],cur_bands[iorb].TwistIndex,iorb,0);
    }
    return foundit;
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
      Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());;
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
      APP_ABORT("SplineHybridAdoptorReader Failed to read band(s)");
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
      Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());;
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
      APP_ABORT("SplineHybridAdoptorReader Failed to read band(s)");
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
