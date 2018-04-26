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

namespace qmcplusplus
{

template<typename ST, typename LT>
struct Gvectors
{
  typedef TinyVector<ST,3> PosType;
  typedef std::complex<ST> ValueType;

  const LT& Lattice;
  std::vector<PosType>                   gvecs_cart; //Cartesian.
  std::vector<ST>                        gmag;
  const size_t NumGvecs;

  Gvectors(const std::vector<TinyVector<int,3> >& gvecs_in,
           const LT& Lattice_in,
           const TinyVector<int,3>& HalfG,
           size_t first, size_t last):
  Lattice(Lattice_in), NumGvecs(last-first)
  {
    gvecs_cart.resize(NumGvecs);
    gmag.resize(NumGvecs);
    #pragma omp parallel for
    for(size_t ig=0; ig<NumGvecs; ig++)
    {
      TinyVector<ST,3> gvec_shift;
      gvec_shift=gvecs_in[ig+first]+HalfG*0.5;
      gvecs_cart[ig]=Lattice.k_cart(gvec_shift);
      gmag[ig]=std::sqrt(dot(gvecs_cart[ig],gvecs_cart[ig]));
    }
  }

  template<typename YLM_ENGINE, typename VVT>
  void calc_Ylm_G(const size_t ig, YLM_ENGINE &Ylm, VVT &YlmG) const
  {
    PosType Ghat(0.0,0.0,1.0);
    if(gmag[ig]>0)
      Ghat=gvecs_cart[ig]/gmag[ig];
    Ylm.evaluateV(Ghat[0], Ghat[1], Ghat[2], YlmG.data());
  }

  template<typename VVT>
  inline void calc_jlm_G(const int lmax, ST &r, const size_t ig, VVT &j_lm_G) const
  {
    bessel_steed_array_cpu(lmax, gmag[ig]*r, j_lm_G.data());
    for(size_t l=lmax; l>0; l--)
      for(size_t lm=l*l; lm<(l+1)*(l+1); lm++)
        j_lm_G[lm]=j_lm_G[l];
  }

  template<typename PT, typename VT>
  inline void calc_phase_shift(const PT& RSoA, const size_t ig, VT &phase_shift_real, VT &phase_shift_imag) const
  {
    const ST* restrict px=RSoA.data(0);
    const ST* restrict py=RSoA.data(1);
    const ST* restrict pz=RSoA.data(2);
    ST* restrict v_r=phase_shift_real.data();
    ST* restrict v_i=phase_shift_imag.data();
    const ST &gv_x=gvecs_cart[ig][0];
    const ST &gv_y=gvecs_cart[ig][1];
    const ST &gv_z=gvecs_cart[ig][2];

    #pragma omp simd aligned(px,py,pz,v_r,v_i)
    for(size_t iat=0; iat<RSoA.size(); iat++)
      sincos(px[iat]*gv_x+py[iat]*gv_y+pz[iat]*gv_z,v_i+iat,v_r+iat);
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

  SplineHybridAdoptorReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e), spline_r(NULL), spline_i(NULL),
      bspline(0), FFTplan(NULL)
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
      h5f.close();
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
        h5f.close();
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
    auto& ACInfo=mybuilder->AtomicCentersInfo;
    // load atomic center info only when it is not initialized
    if(centers.size()==0)
    {
      bool success=true;
      app_log() << "Reading atomic center info for hybrid representation" << std::endl;
      for(int center_idx=0; center_idx<ACInfo.Ncenters; center_idx++)
      {
        const int my_GroupID = ACInfo.GroupID[center_idx];
        if(ACInfo.cutoff[center_idx]<0)
        {
          app_error() << "Hybrid orbital representation needs parameter 'cutoff_radius' for atom " << center_idx << std::endl;
          success=false;
        }

        if(ACInfo.inner_cutoff[center_idx]<0)
        {
          const double inner_cutoff = std::max(ACInfo.cutoff[center_idx]-0.3, 0.0);
          app_log() << "Hybrid orbital representation setting 'inner_cutoff' to "
                    << inner_cutoff << " for group " << my_GroupID << " as atom " << center_idx << std::endl;
          // overwrite the inner_cutoff of all the atoms of the same species
          for(int id=0; id<ACInfo.Ncenters; id++)
            if(my_GroupID==ACInfo.GroupID[id])
              ACInfo.inner_cutoff[id] = inner_cutoff;
        }
        else if(ACInfo.inner_cutoff[center_idx]>ACInfo.cutoff[center_idx])
        {
          app_error() << "Hybrid orbital representation 'inner_cutoff' must be smaller than 'spline_radius' for atom " << center_idx << std::endl;
          success=false;
        }

        if(ACInfo.cutoff[center_idx]>0)
        {
          if(ACInfo.lmax[center_idx]<0)
          {
            app_error() << "Hybrid orbital representation needs parameter 'lmax' for atom " << center_idx << std::endl;
            success=false;
          }

          if(ACInfo.spline_radius[center_idx]<0 && ACInfo.spline_npoints[center_idx]<0)
          {
            app_log() << "Parameters 'spline_radius' and 'spline_npoints' for group " << my_GroupID
                      << " as atom " << center_idx << " are not specified." << std::endl;
            const double delta = std::min(0.02, ACInfo.cutoff[center_idx]/4.0);
            const int n_grid_point = std::ceil((ACInfo.cutoff[center_idx] + 1e-4) / delta) + 3;
            for(int id=0; id<ACInfo.Ncenters; id++)
              if(my_GroupID==ACInfo.GroupID[id])
              {
                ACInfo.spline_npoints[id] = n_grid_point;
                ACInfo.spline_radius[id] = (n_grid_point-1) * delta;
              }
            app_log() << "  Based on default grid point distance " << delta << std::endl;
            app_log() << "  Setting 'spline_npoints' to " << ACInfo.spline_npoints[center_idx] << std::endl;
            app_log() << "  Setting 'spline_radius' to " << ACInfo.spline_radius[center_idx] << std::endl;
          }
          else
          {
            if(ACInfo.spline_radius[center_idx]<0)
            {
              app_error() << "Hybrid orbital representation needs parameter 'spline_radius' for atom " << center_idx << std::endl;
              success=false;
            }

            if(ACInfo.spline_npoints[center_idx]<0)
            {
              app_error() << "Hybrid orbital representation needs parameter 'spline_npoints' for atom " << center_idx << std::endl;
              success=false;
            }
          }

          // check maximally allowed cutoff_radius
          double max_allowed_cutoff=ACInfo.spline_radius[center_idx]-2.0*ACInfo.spline_radius[center_idx]/(ACInfo.spline_npoints[center_idx]-1);
          if(success && ACInfo.cutoff[center_idx]>max_allowed_cutoff)
          {
            app_error() << "Hybrid orbital representation requires cutoff_radius<=" << max_allowed_cutoff
                        << " calculated by spline_radius-2*spline_radius/(spline_npoints-1) for atom " << center_idx << std::endl;
            success=false;
          }
        }
        else
        {
          // no atomic regions for this atom type
          ACInfo.spline_radius[center_idx] = 0.0;
          ACInfo.spline_npoints[center_idx] = 0;
          ACInfo.lmax[center_idx] = 0;
        }
      }
      if(!success) abort();

      for(int center_idx=0; center_idx<ACInfo.Ncenters; center_idx++)
      {
        AtomicOrbitalSoA<DataType> oneCenter(ACInfo.lmax[center_idx]);
        oneCenter.set_info(ACInfo.ion_pos[center_idx], ACInfo.cutoff[center_idx],
                           ACInfo.inner_cutoff[center_idx], ACInfo.spline_radius[center_idx],
                           ACInfo.non_overlapping_radius[center_idx], ACInfo.spline_npoints[center_idx]);
        centers.push_back(oneCenter);
      }
    }
  }

  /** initialize construct atomic orbital radial functions from plane waves */
  inline void create_atomic_centers_Gspace(Vector<std::complex<double> >& cG, Communicate& band_group_comm, int iorb)
  {
    //distribute G-vectors over processor groups
    const int Ngvecs=mybuilder->Gvecs[0].size();
    const int Nprocs=band_group_comm.size();
    const int Ngvecgroups=std::min(Ngvecs,Nprocs);
    Communicate gvec_group_comm(band_group_comm, Ngvecgroups);
    std::vector<int> gvec_groups(Ngvecgroups+1,0);
    FairDivideLow(Ngvecs,Ngvecgroups,gvec_groups);
    const int gvec_first=gvec_groups[gvec_group_comm.getGroupID()];
    const int gvec_last=gvec_groups[gvec_group_comm.getGroupID()+1];

    // prepare Gvecs Ylm(G)
    typedef typename EinsplineSetBuilder::UnitCellType UnitCellType;
    Gvectors<double, UnitCellType> Gvecs(mybuilder->Gvecs[0], mybuilder->PrimCell, bspline->HalfG,
                                         gvec_first, gvec_last);
    // if(band_group_comm.isGroupLeader()) std::cout << "print band=" << iorb << " KE=" << Gvecs.evaluate_KE(cG) << std::endl;

    std::vector<AtomicOrbitalSoA<DataType> >& centers=bspline->AtomicCenters;
    app_log() << "Transforming band " << iorb << " on Rank 0" << std::endl;
    // collect atomic centers by group
    std::vector<int> uniq_species;
    for(int center_idx=0; center_idx<centers.size(); center_idx++)
    {
      auto& ACInfo = mybuilder->AtomicCentersInfo;
      const int my_GroupID = ACInfo.GroupID[center_idx];
      int found_idx = -1;
      for(size_t idx=0; idx<uniq_species.size(); idx++)
        if(my_GroupID==uniq_species[idx])
        {
          found_idx = idx; break;
        }
      if(found_idx<0) uniq_species.push_back(my_GroupID);
    }
    // construct group list
    std::vector<std::vector<int> > group_list(uniq_species.size());
    for(int center_idx=0; center_idx<centers.size(); center_idx++)
    {
      auto& ACInfo = mybuilder->AtomicCentersInfo;
      const int my_GroupID = ACInfo.GroupID[center_idx];
      for(size_t idx=0; idx<uniq_species.size(); idx++)
        if(my_GroupID==uniq_species[idx])
        {
          group_list[idx].push_back(center_idx);
          break;
        }
    }

    for(int group_idx=0; group_idx<group_list.size(); group_idx++)
    {
      const auto &mygroup = group_list[group_idx];
      const double spline_radius = centers[mygroup[0]].spline_radius;
      const int   spline_npoints = centers[mygroup[0]].spline_npoints;
      const int             lmax = centers[mygroup[0]].lmax;
      const double delta = spline_radius/static_cast<double>(spline_npoints-1);
      const int lm_tot=(lmax+1)*(lmax+1);
      const size_t natoms = mygroup.size();
      const int policy = lm_tot>natoms?0:1;

      std::vector<std::complex<double> > i_power(lm_tot);
      // rotate phase is introduced here.
      std::complex<double> i_temp(rotate_phase_r, rotate_phase_i);
      for(size_t l=0; l<=lmax; l++)
      {
        for(size_t lm=l*l; lm<(l+1)*(l+1); lm++)
          i_power[lm]=i_temp;
        i_temp*=std::complex<double>(0.0,1.0);
      }

      std::vector<Matrix<double> > all_vals(natoms);
      std::vector<std::vector<aligned_vector<double> > > vals_local(spline_npoints*omp_get_max_threads());
      VectorSoaContainer<double,3> myRSoA(natoms);
      for(size_t idx=0; idx<natoms; idx++)
      {
        all_vals[idx].resize(spline_npoints, lm_tot*2);
        all_vals[idx]=0.0;
        myRSoA(idx)=centers[mygroup[idx]].pos;
      }

      #pragma omp parallel
      {
        const size_t tid = omp_get_thread_num();
        const size_t nt = omp_get_num_threads();

        for(int ip=0; ip<spline_npoints; ip++)
        {
          const size_t ip_idx=tid*spline_npoints+ip;
          if(policy==1)
          {
            vals_local[ip_idx].resize(lm_tot*2);
            for(size_t lm=0; lm<lm_tot*2; lm++)
            {
              auto &vals = vals_local[ip_idx][lm];
              vals.resize(natoms);
              std::fill(vals.begin(),vals.end(),0.0);
            }
          }
          else
          {
            vals_local[ip_idx].resize(natoms*2);
            for(size_t iat=0; iat<natoms*2; iat++)
            {
              auto &vals = vals_local[ip_idx][iat];
              vals.resize(lm_tot);
              std::fill(vals.begin(),vals.end(),0.0);
            }
          }
        }

        const size_t size_pw_tile = 32;
        const size_t num_pw_tiles = (Gvecs.NumGvecs+size_pw_tile-1)/size_pw_tile;
        aligned_vector<double> j_lm_G(lm_tot,0.0);
        std::vector<aligned_vector<double> > phase_shift_r(size_pw_tile);
        std::vector<aligned_vector<double> > phase_shift_i(size_pw_tile);
        std::vector<aligned_vector<double> > YlmG(size_pw_tile);
        for(size_t ig=0; ig<size_pw_tile; ig++)
        {
          phase_shift_r[ig].resize(natoms);
          phase_shift_i[ig].resize(natoms);
          YlmG[ig].resize(lm_tot);
        }
        SoaSphericalTensor<double> Ylm(lmax);

        #pragma omp for
        for(size_t tile_id=0; tile_id<num_pw_tiles; tile_id++)
        {
          const size_t ig_first = tile_id*size_pw_tile;
          const size_t ig_last  = std::min((tile_id+1)*size_pw_tile,Gvecs.NumGvecs);
          for(size_t ig=ig_first; ig<ig_last; ig++)
          {
            const size_t ig_local=ig-ig_first;
            // calculate phase shift for all the centers of this group
            Gvecs.calc_phase_shift(myRSoA, ig, phase_shift_r[ig_local], phase_shift_i[ig_local]);
            Gvecs.calc_Ylm_G(ig, Ylm, YlmG[ig_local]);
          }

          for(int ip=0; ip<spline_npoints; ip++)
          {
            double r=delta*static_cast<double>(ip);
            const size_t ip_idx=tid*spline_npoints+ip;

            for(size_t ig=ig_first; ig<ig_last; ig++)
            {
              const size_t ig_local=ig-ig_first;
              // calculate spherical bessel function
              Gvecs.calc_jlm_G(lmax, r, ig, j_lm_G);
              for(size_t lm=0; lm<lm_tot; lm++)
                j_lm_G[lm]*=YlmG[ig_local][lm];

              const double cG_r=cG[ig+gvec_first].real();
              const double cG_i=cG[ig+gvec_first].imag();
              if(policy==1)
              {
                for(size_t lm=0; lm<lm_tot; lm++)
                {
                  double* restrict vals_r = vals_local[ip_idx][lm*2].data();
                  double* restrict vals_i = vals_local[ip_idx][lm*2+1].data();
                  const double* restrict ps_r_ptr = phase_shift_r[ig_local].data();
                  const double* restrict ps_i_ptr = phase_shift_i[ig_local].data();
                  double cG_j_r=cG_r*j_lm_G[lm];
                  double cG_j_i=cG_i*j_lm_G[lm];
                  #pragma omp simd aligned(vals_r,vals_i,ps_r_ptr,ps_i_ptr)
                  for(size_t idx=0; idx<natoms; idx++)
                  {
                    const double ps_r=ps_r_ptr[idx];
                    const double ps_i=ps_i_ptr[idx];
                    vals_r[idx]+=cG_j_r*ps_r-cG_j_i*ps_i;
                    vals_i[idx]+=cG_j_i*ps_r+cG_j_r*ps_i;
                  }
                }
              }
              else
              {
                for(size_t idx=0; idx<natoms; idx++)
                {
                  double* restrict vals_r = vals_local[ip_idx][idx*2].data();
                  double* restrict vals_i = vals_local[ip_idx][idx*2+1].data();
                  const double* restrict j_lm_G_ptr = j_lm_G.data();
                  double cG_ps_r=cG_r*phase_shift_r[ig_local][idx]-cG_i*phase_shift_i[ig_local][idx];
                  double cG_ps_i=cG_i*phase_shift_r[ig_local][idx]+cG_r*phase_shift_i[ig_local][idx];
                  #pragma omp simd aligned(vals_r,vals_i,j_lm_G_ptr)
                  for(size_t lm=0; lm<lm_tot; lm++)
                  {
                    const double jlm=j_lm_G_ptr[lm];
                    vals_r[lm]+=cG_ps_r*jlm;
                    vals_i[lm]+=cG_ps_i*jlm;
                  }
                }
              }
            }
          }
        }

        #pragma omp for collapse(2)
        for(int ip=0; ip<spline_npoints; ip++)
          for(size_t idx=0; idx<natoms; idx++)
          {
            double *vals = all_vals[idx][ip];
            for(size_t tid=0; tid<nt; tid++)
              for(size_t lm=0; lm<lm_tot; lm++)
              {
                double vals_th_r, vals_th_i;
                const size_t ip_idx=tid*spline_npoints+ip;
                if(policy==1)
                {
                  vals_th_r = vals_local[ip_idx][lm*2][idx];
                  vals_th_i = vals_local[ip_idx][lm*2+1][idx];
                }
                else
                {
                  vals_th_r = vals_local[ip_idx][idx*2][lm];
                  vals_th_i = vals_local[ip_idx][idx*2+1][lm];
                }
                const double real_tmp = 4.0*M_PI*i_power[lm].real();
                const double imag_tmp = 4.0*M_PI*i_power[lm].imag();
                vals[lm]        += vals_th_r*real_tmp - vals_th_i*imag_tmp;
                vals[lm+lm_tot] += vals_th_i*real_tmp + vals_th_r*imag_tmp;
              }
          }
      }
      //app_log() << "Building band " << iorb << " at center " << center_idx << std::endl;

      for(size_t idx=0; idx<natoms; idx++)
      {
        // reduce all_vals
        band_group_comm.reduce_in_place(all_vals[idx].data(),all_vals[idx].size());
        if(!band_group_comm.isGroupLeader()) continue;
        #pragma omp parallel for
        for(int lm=0; lm<lm_tot; lm++)
        {
          auto &mycenter = centers[mygroup[idx]];
          aligned_vector<double> splineData_r(spline_npoints);
          UBspline_1d_d* atomic_spline_r;
          for(size_t ip=0; ip<spline_npoints; ip++)
            splineData_r[ip]=all_vals[idx][ip][lm];
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
              splineData_i[ip]=all_vals[idx][ip][lm+lm_tot];
            atomic_spline_i=einspline::create(atomic_spline_i, 0.0, spline_radius, spline_npoints, splineData_i.data(), ((lm==0)||(lm>3)));
            int iband=bspline->BandIndexMap.size()>0?bspline->BandIndexMap[iorb]:iorb;
            mycenter.set_spline(atomic_spline_r,lm,iband*2);
            mycenter.set_spline(atomic_spline_i,lm,iband*2+1);
            einspline::destroy(atomic_spline_r);
            einspline::destroy(atomic_spline_i);
          }
        }
      }

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
          fprintf(fout_pw, "%15.10lf  %15.10lf  ", all_vals[center_idx][ip][lm].real(), all_vals[center_idx][ip][lm].imag());
          fprintf(fout_spline_v, "%15.10lf  %15.10lf  ", mycenter.localV[lm*mycenter.Npad+iorb*2], mycenter.localV[lm*mycenter.Npad+iorb*2+1]);
          fprintf(fout_spline_g, "%15.10lf  %15.10lf  ", mycenter.localG[lm*mycenter.Npad+iorb*2], mycenter.localG[lm*mycenter.Npad+iorb*2+1]);
          fprintf(fout_spline_l, "%15.10lf  %15.10lf  ", mycenter.localL[lm*mycenter.Npad+iorb*2], mycenter.localL[lm*mycenter.Npad+iorb*2+1]);
        }
        fprintf(fout_pw, "\n");
        fprintf(fout_spline_v, "\n");
        fprintf(fout_spline_g, "\n");
        fprintf(fout_spline_l, "\n");
      }
      fclose(fout_pw);
      fclose(fout_spline_v);
      fclose(fout_spline_g);
      fclose(fout_spline_l);
#endif
    }
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
        double total_norm = compute_norm(cG);
        if(std::abs(total_norm-1.0)>PW_COEFF_NORM_TOLERANCE)
        {
          std::cerr << "The orbital " << iorb << " has a wrong norm " << total_norm
                    << ", computed from plane wave coefficients!" << std::endl
                    << "This may indicate a problem with the HDF5 library versions used "
                    << "during wavefunction conversion or read." << std::endl;
          APP_ABORT("SplineHybridAdoptorReader Wrong orbital norm!");
        }
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
    if(band_group_comm.isGroupLeader())
    {
      now.restart();
      bspline->gather_tables(band_group_comm.GroupLeaderComm);
      app_log() << "  Time to gather the table = " << now.elapsed() << std::endl;
    }
    now.restart();
    bspline->bcast_tables(myComm);
    app_log() << "  Time to bcast the table = " << now.elapsed() << std::endl;
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
