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
    
    
/** @file SplineMixedAdoptorReader.h
 */

#ifndef QMCPLUSPLUS_EINSPLINE_SPLINE_MIXED_ADOPTOR_PARALLEL_READER_H
#define QMCPLUSPLUS_EINSPLINE_SPLINE_MIXED_ADOPTOR_PARALLEL_READER_H

#include <mpi/point2point.h>
#include <mpi/collectives.h>

namespace qmcplusplus
{

/** reader for a EinsplineAdoptor with truncation
 * @tparam SA spline adoptor, SplineMixedAdoptor, SplineOpenAdoptor
 */
template<typename SA>
struct SplineMixedAdoptorReader: public BsplineReaderBase
{
  typedef SA adoptor_type;
  typedef typename adoptor_type::SplineType SplineType;

  TinyVector<int,3> coarse_mesh;
  TinyVector<int,3> coarse_stride;
  Array<std::complex<double>,3> FFTbox;
  Array<double,3> bigD_r, bigD_i;
  Array<double,3> smallD_r, smallD_i;
  std::vector<UBspline_3d_d*> spline_r;
  std::vector<UBspline_3d_d*> spline_i;
  BsplineSet<adoptor_type>* bspline;
  std::vector<int> OrbGroups;
  fftw_plan FFTplan;
  bool use_imaginary;

  SplineMixedAdoptorReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e),bspline(0),FFTplan(NULL),use_imaginary(false)
  {}

  ~SplineMixedAdoptorReader()
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
    APP_ABORT("SplineMixedAdoptorReader::export_MultiSpline(multi_UBspline_3d_z* target) should not be called");
  }

  void export_MultiSpline(multi_UBspline_3d_d** target)
  {
    *target = (multi_UBspline_3d_d*) bspline->MultiSpline;
  }

  /** fft and spline a band
   * @param ti twist index
   * @param iorb orbital index
   */
  void fft_spline(Vector<std::complex<double> >& cG, int ti, int iorb)
  {
    unpack4fftw(cG,mybuilder->Gvecs[0],mybuilder->MeshSize,FFTbox);
    fftw_execute (FFTplan);
    //fix_phase_rotate_c2r(FFTbox,bigD,TwistAngle);
    if(bspline->is_complex)
      fix_phase_rotate_c2c(FFTbox,bigD_r, bigD_i,mybuilder->TwistAngles[ti]);
    else
      fix_phase_rotate_c2r(FFTbox,bigD_r,mybuilder->TwistAngles[ti]);
    for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=coarse_stride[0])
      for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=coarse_stride[1])
        for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=coarse_stride[2])
          smallD_r(i,j,k)=bigD_r(i2,j2,k2);
    einspline::set(spline_r[2*iorb],bigD_r.data());
    einspline::set(spline_r[2*iorb+1],smallD_r.data());
    if(bspline->is_complex)
    {
      for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=coarse_stride[0])
        for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=coarse_stride[1])
          for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=coarse_stride[2])
            smallD_i(i,j,k)=bigD_i(i2,j2,k2);
      einspline::set(spline_i[2*iorb],bigD_i.data());
      einspline::set(spline_i[2*iorb+1],smallD_i.data());
    }
  }

  SPOSetBase* create_spline_set(int spin, const BandInfoGroup& bandgroup)
  {
    ReportEngine PRE("SplineOpenAdoptorReader","create_spline_set(int, EinsplineSet*)");
    //typedef typename adoptor_type::real_type spline_data_type;
    bspline=new BsplineSet<adoptor_type>;
    app_log() << "  AdoptorName = " << bspline->AdoptorName << std::endl;
    if(bspline->is_complex)
      app_log() << "  Using complex einspline table" << std::endl;
    else
      app_log() << "  Using real einspline table" << std::endl;

    check_twists(bspline,bandgroup);
    Ugrid xyz_grid[3];
    typename adoptor_type::BCType xyz_bc[3];

    bool havePsig=set_grid(bspline->HalfG,xyz_grid, xyz_bc);
    if(!havePsig)
    {
      APP_ABORT("Need psi_g with truncate=\"yes\"");
    }
    TinyVector<int,3> bconds=mybuilder->TargetPtcl.Lattice.BoxBConds;
    bool use_cartesian= (mybuilder->TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN);

    //coarse mesh only along the open direction
    for(int j=0; j<3; ++j)
    {
      coarse_mesh[j]=(bconds[j])?MeshSize[j]:MeshSize[j]/2;
      coarse_stride[j]=(bconds[j])?1:2;
    }
    TinyVector<double,3> start(0.0);
    TinyVector<double,3> end(1.0);
    TinyVector<double,3> maxL(bspline->PrimLattice.R(0,0)
                              ,bspline->PrimLattice.R(1,1),bspline->PrimLattice.R(2,2));
    if(use_cartesian)
      end=maxL;
    //determine the bonding box
    double buffer=mybuilder->BufferLayer; //this has to be option
    TinyVector<double,3> lower=end;
    TinyVector<double,3> upper=start;
    const Vector<TinyVector<double,3> >& IonPos(mybuilder->IonPos);;
    if(use_cartesian)
    {
      for(int i=0; i<IonPos.size(); ++i)
      {
        for(int j=0; j<3; ++j)
        {
          lower[j]=std::min(IonPos[i][j],lower[j]);
          upper[j]=std::max(IonPos[i][j],upper[j]);
        }
      }
      //add buffer
      for(int j=0; j<3; ++j)
      {
        lower[j]=std::max(lower[j]-buffer,start[j]);
        upper[j]=std::min(upper[j]+buffer,end[j]);
      }
    }
    else
    {
      for(int i=0; i<IonPos.size(); ++i)
      {
        TinyVector<double,3> pos=bspline->PrimLattice.toUnit(IonPos[i]);
        for(int j=0; j<3; ++j)
        {
          if(!bconds[j])
          {
            lower[j]=std::min(pos[j],lower[j]);
            upper[j]=std::max(pos[j],upper[j]);
          }
        }
      }
      for(int j=0; j<3; ++j)
      {
        if(bconds[j])
        {
          lower[j]=0.0;
          upper[j]=1.0;
        }
        else
        {
          lower[j]=std::max(lower[j]-buffer/maxL[j],start[j]);
          upper[j]=std::min(upper[j]+buffer/maxL[j],end[j]);
        }
      }
    }
    int N = mybuilder->NumDistinctOrbitals;
    //bspline->create_spline(MeshSize,N,fullgrid);
    bspline->create_spline(MeshSize,coarse_mesh,N,N);
    int np=std::min(N,myComm->size());
    OrbGroups.resize(np+1,0);
    FairDivideLow(N,np,OrbGroups);
    int ip=myComm->rank();
    int norbs_n=(ip<np)?OrbGroups[ip+1]-OrbGroups[ip]:0;
    UBspline_3d_d* dummy=0;
    spline_r.resize(norbs_n*2+2,0);
    for(int i=0; i<spline_r.size()/2; ++i)
    {
      spline_r[2*i]=einspline::create(dummy,start,end,MeshSize,bspline->HalfG);
      spline_r[2*i+1]=einspline::create(dummy,start,end,coarse_mesh,bspline->HalfG);
    }
    spline_i.resize(norbs_n*2+2,0);
    if(bspline->is_complex)
    {
      use_imaginary=true;
      for(int i=0; i<spline_i.size()/2; ++i)
      {
        spline_i[2*i]=einspline::create(dummy,start,end,MeshSize,bspline->HalfG);
        spline_i[2*i+1]=einspline::create(dummy,start,end,coarse_mesh,bspline->HalfG);
      }
    }
    app_log() << "  Original Mesh " << MeshSize << std::endl;
    app_log() << "  Coarse Mesh " << coarse_mesh << std::endl;
    if(use_cartesian)
      app_log() << "  Using Cartesian grids for open systems. " << std::endl;
    else
      app_log() << "  Using Primitive-cell grids for mixed systems. " << std::endl;
    app_log() << "  Using buffer layer for the small box= " << buffer << " bohr " << std::endl;
    app_log() << "  Adding a small box" << "\n  LowerBound " << lower << "\n  UpperBound " << upper << std::endl;
    app_log().flush();
    bspline->add_box(spline_r[0],lower,upper);
    int foundspline=0;
    std::string splinefile=make_spline_filename(mybuilder->H5FileName,spin,mybuilder->TwistNum,MeshSize);
    Timer now;
    if(myComm->rank()==0)
    {
      now.restart();
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
        foundspline=bspline->read_splines(h5f);
    }
    myComm->bcast(foundspline);
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << std::endl;
      now.restart();
      chunked_bcast(myComm, bspline->MultiSpline);
      chunked_bcast(myComm, bspline->smallBox);
      app_log() << "   Bcast Time for the big table = " << now.elapsed() << std::endl;
    }
    else
    {
      app_log() << "Perform FFT+spline and dump bspline tables to " << splinefile << std::endl;
      FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
      FFTplan = fftw_plan_dft_3d
                (MeshSize[0], MeshSize[1], MeshSize[2],
                 reinterpret_cast<fftw_complex*>(FFTbox.data()),
                 reinterpret_cast<fftw_complex*>(FFTbox.data()),
                 +1, FFTW_ESTIMATE);
      bigD_r.resize(MeshSize[0],MeshSize[1],MeshSize[2]);
      smallD_r.resize(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      if(bspline->is_complex)
      {
        bigD_i.resize(MeshSize[0],MeshSize[1],MeshSize[2]);
        smallD_i.resize(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      }
      Timer clock;
      initialize_spline_pio(spin,bandgroup);
      app_log() << "Time with PIO " << clock.elapsed() << std::endl;
      //clock.restart();
      //initialize_spline_pio_bcast(spin);
      //app_log() << "Time with PIO/Bcast " << clock.elapsed() << std::endl;
      //clock.restart();
      //initialize_spline_slow(spin);
      //app_log() << "Time slow method " << clock.elapsed() << std::endl;
      clock.restart();
      if(qmc_common.save_wfs && myComm->rank()==0)
      {
        clock.restart();
        hdf_archive h5f;
        h5f.create(splinefile);
        h5f.write(bspline->AdoptorName,"adoptor_name");
        int sizeD=sizeof(typename adoptor_type::DataType);
        h5f.write(sizeD,"sizeof");
        bspline->write_splines(h5f);
        app_log() << "Time write TABLE " << clock.elapsed() << std::endl;
      }
    }
    app_log() << "TIME READBANDS " << now.elapsed() << std::endl;

    clear();
    return bspline;
  }

  void initialize_spline_slow(int spin, const BandInfoGroup& bandgroup)
  {
    Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());
    int N = mybuilder->NumDistinctOrbitals;
    const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      int ti=cur_bands[iorb].TwistIndex;
      get_psi_g(ti,spin,cur_bands[iorb].BandIndex,cG);
      fft_spline(cG,ti,0);
      bspline->set_spline(spline_r[1],spline_i[1],ti,ival,0);//level0
      bspline->set_spline(spline_r[0],spline_i[0],ti,ival,1);//level1
    }
  }

  /** initialize big spline using parallel read
   *
   * - initialize a number of bands per MPI task
   * - send bspline coefficients to the root
   * - root complete the table including data transpose
   * - bcast the big table
   * Minimize temporary memory use and buffering.
   */
  void initialize_spline_pio(int spin, const BandInfoGroup& bandgroup)
  {
    int N = mybuilder->NumDistinctOrbitals;
    int np=OrbGroups.size()-1;
    bool root=(myComm->rank()==0);
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
      if(root)
      {
        for(int iorb=OrbGroups[0],ib2=0; iorb<OrbGroups[1]; ++iorb,ib2+=2)
        {
          bspline->set_spline(spline_r[ib2],  spline_i[ib2],  cur_bands[iorb].TwistIndex, iorb,1);
          bspline->set_spline(spline_r[ib2+1],spline_i[ib2+1],cur_bands[iorb].TwistIndex, iorb,0);
        }
      }
    }//root put splines to the big table
    myComm->barrier();
    myComm->bcast(foundit);
    if(!foundit)
      APP_ABORT("SplineMixedAdoptorReader Failed to read band(s)");
    //mpi needs int
    int ng_big=static_cast<int>(spline_r[0]->coefs_size);
    int ng_small=static_cast<int>(spline_r[1]->coefs_size);
    //send back to zero  using synchronous send/recv
    for(int ip=1; ip<np; ++ip)
    {
      int remote=0;
      if(ip==myComm->rank())
      {
        for(int iorb=OrbGroups[ip],ib2=0; iorb<OrbGroups[ip+1]; ++iorb,ib2+=2)
        {
          mpi::send(*myComm,spline_r[ib2]->coefs,ng_big,remote,iorb);
          mpi::send(*myComm,spline_r[ib2+1]->coefs,ng_small,remote,iorb+N);
          if(bspline->is_complex)
          {
            mpi::send(*myComm,spline_i[ib2]->coefs,ng_big,remote,iorb+2*N);
            mpi::send(*myComm,spline_i[ib2+1]->coefs,ng_small,remote,iorb+3*N);
          }
        }
      }
      else
        if(root)
        {
          for(int iorb=OrbGroups[ip],ib2=0; iorb<OrbGroups[ip+1]; ++iorb,ib2+=2)
          {
            mpi::recv(*myComm,spline_r[ib2]->coefs,ng_big,ip,iorb);
            mpi::recv(*myComm,spline_r[ib2+1]->coefs,ng_small,ip,iorb+N);
            if(bspline->is_complex)
            {
              mpi::recv(*myComm,spline_i[ib2]->coefs,ng_big,ip,iorb+2*N);
              mpi::recv(*myComm,spline_i[ib2+1]->coefs,ng_small,ip,iorb+3*N);
            }
          }
          for(int iorb=OrbGroups[ip],ib2=0; iorb<OrbGroups[ip+1]; ++iorb,ib2+=2)
          {
            bspline->set_spline(spline_r[ib2],  spline_i[ib2],  cur_bands[iorb].TwistIndex, iorb,1);
            bspline->set_spline(spline_r[ib2+1],spline_i[ib2+1],cur_bands[iorb].TwistIndex, iorb,0);
          }
        }
    }
    myComm->barrier();
    chunked_bcast(myComm, bspline->MultiSpline);
    chunked_bcast(myComm, bspline->smallBox);
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
      APP_ABORT("SplineMixedAdoptorReader Failed to read band(s)");
    int ng_big=static_cast<int>(spline_r[0]->coefs_size);
    int ng_small=static_cast<int>(spline_r[1]->coefs_size);
    //pointers to UBspline_3d_d for bcast without increaing mem
    UBspline_3d_d *dense_r=0, *dense_i=0;
    UBspline_3d_d *coarse_r=0, *coarse_i=0;
    int iorb_target=(myComm->rank()<np)?(OrbGroups[myComm->rank()+1]-OrbGroups[myComm->rank()])*2:0;
    for(int ip=0; ip<np; ++ip)
    {
      for(int iorb=OrbGroups[ip], ib2=0; iorb<OrbGroups[ip+1]; ++iorb, ib2+=2)
      {
        if(ip==myComm->rank())
        {
          dense_r=spline_r[ib2];
          coarse_r=spline_r[ib2+1];
          if(bspline->is_complex)
          {
            dense_i=spline_i[ib2];
            coarse_i=spline_i[ib2+1];
          }
        }
        else
        {
          //everyone else
          dense_r=spline_r[iorb_target];
          coarse_r=spline_r[iorb_target+1];
          if(bspline->is_complex)
          {
            dense_i=spline_i[iorb_target];
            coarse_i=spline_i[iorb_target+1];
          }
        }
        myComm->barrier();
        mpi::bcast(*myComm,dense_r->coefs,ng_big,ip);
        mpi::bcast(*myComm,coarse_r->coefs,ng_small,ip);
        if(bspline->is_complex)
        {
          mpi::bcast(*myComm,dense_i->coefs,ng_big,ip);
          mpi::bcast(*myComm,coarse_i->coefs,ng_small,ip);
        }
        myComm->barrier();
        bspline->set_spline(coarse_r,coarse_i,cur_bands[iorb].TwistIndex, iorb,0);
        bspline->set_spline(dense_r,dense_i,  cur_bands[iorb].TwistIndex, iorb,1);
      }
    }
  }
};

//  template<typename T>
//    struct SplineOpenAdoptorReader: BsplineReaderBase
//    {
//      typedef SplineOpenAdoptor<T,double,3> adoptor_type;
//
//      SplineOpenAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
//      {}
//
//      SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
//      {
//        ReportEngine PRE("SplineOpenAdoptorReader","create_spline_set(int, EinsplineSet*)");
//
//        typedef T spline_data_ype;
//
//        BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
//        init(orbitalSet,bspline);
//        int N=bspline->getOrbitalSetSize();
//
//        bspline->resizeStorage(N,N);
//
//        std::string H5FileName(mybuilder->H5FileName);
//        H5OrbSet set(H5FileName, spin, N);
//
//        bool havePsir=!(mybuilder->ReadGvectors_ESHDF());
//
//        bool fullgrid=false;
//
//        TinyVector<int,3> MeshSize=mybuilder->MeshSize;
//        TinyVector<int,3> coarse_mesh(MeshSize[0]/2,MeshSize[1]/2,MeshSize[2]/2);
//
//        TinyVector<double,3> start(0.0);
//        TinyVector<double,3> end(1.0);
//        for(int i=0; i<3; ++i) end(i)=mybuilder->SuperLattice(i,i);
//
//        //create dense and coarse UBspline_3d_d
//        UBspline_3d_d* dense=0;
//        dense=einspline::create(dense,start,end,MeshSize,PERIODIC);
//        UBspline_3d_d* coarse=0;
//        coarse=einspline::create(coarse,start,end,coarse_mesh,PERIODIC);
//
//        //determine the bonding box
//        TinyVector<double,3> lower=end;
//        TinyVector<double,3> upper=start;
//        const Vector<TinyVector<double,OHMMS_DIM> >& IonPos(mybuilder->IonPos);;
//        for(int i=0; i<IonPos.size(); ++i)
//        {
//          for(int j=0; j<3; ++j)
//          {
//            lower[j]=std::min(IonPos[i][j],lower[j]);
//            upper[j]=std::max(IonPos[i][j],upper[j]);
//          }
//        }
//        //add 2 bohr
//        for(int j=0; j<3; ++j)
//        {
//          lower[j]=std::max(lower[j]-2.0,start[j]);
//          upper[j]=std::min(upper[j]+2.0,end[j]);
//        }
//
//        bspline->create_spline(MeshSize,N,fullgrid);
//
//        app_log() << "Original Mesh " << MeshSize << std::endl;
//        app_log() << "Coarse Mesh " << coarse_mesh << std::endl;
//
//        app_log() << "  Adding a small box" << "\n  LowerBound " << lower << "\n  UpperBound " << upper << std::endl;
//
//        app_log().flush();
//
//        bspline->add_box(dense,lower,upper);
//
//        int foundspline=0;
//
//        std::string splinefile=make_spline_filename(H5FileName,spin,0,MeshSize);
//
//        Timer now;
//
//        if(myComm->rank()==0)
//        {
//          hdf_archive h5f;
//          foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
//          if(foundspline)
//          {
//            std::string aname("none");
//            foundspline = h5f.read(aname,"adoptor_name");
//            foundspline = (aname.find(bspline->KeyWord) != std::string::npos);
//          }
//          if(foundspline)
//          {
//            TinyVector<double,3> lower_in(end);
//            TinyVector<double,3> upper_in(start);
//            h5f.read(lower_in,"lower_bound");
//            h5f.read(upper_in,"upper_bound");
//
//            lower_in-=lower; upper_in-=upper;
//            if(dot(lower_in,lower_in)<1e-12 &&dot(upper_in,upper_in)<1e-12)
//            {
//              einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
//              einspline_engine<typename adoptor_type::SplineType> smalltable(bspline->smallBox);
//              foundspline=h5f.read(bigtable,"spline_0");
//              foundspline=h5f.read(smalltable,"spline_1");
//            }
//            else
//            {
//              app_log() << "  The upper/lower bound of the input is different from the current value."<< std::endl;
//              foundspline=0;
//            }
//          }
//        }
//
//        myComm->bcast(foundspline);
//
//        if(foundspline)
//        {
//          app_log() << "Use existing bspline tables in " << splinefile << std::endl;
//          chunked_bcast(myComm, bspline->MultiSpline);
//          chunked_bcast(myComm, bspline->smallBox);
//        }
//        else
//        {
//          app_log() << "Perform FFT+spline and dump bspline tables to " << splinefile << std::endl;
//          Array<std::complex<double>,3> FFTbox;
//          FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
//          fftw_plan FFTplan = fftw_plan_dft_3d
//            (MeshSize[0], MeshSize[1], MeshSize[2],
//             reinterpret_cast<fftw_complex*>(FFTbox.data()),
//             reinterpret_cast<fftw_complex*>(FFTbox.data()),
//             +1, FFTW_ESTIMATE);
//
//          Array<double,3> bigD(MeshSize[0],MeshSize[1],MeshSize[2]);
//          Array<double,3> smallD(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
//
//          TinyVector<double,3> TwistAngle(0.0,0.0,0.0);
//          int ti=0;
//          int ncg=mybuilder->Gvecs[ti].size();
//          Vector<std::complex<double> > cG(ncg);
//
//          for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
//          {
//            get_psi_g(ti,spin,mybuilder->SortBands[iorb].BandIndex,cG);
//            unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
//            fftw_execute (FFTplan);
//            fix_phase_rotate_c2r(FFTbox,bigD,TwistAngle);
//
//            for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=2)
//              for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=2)
//                for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=2)
//                  smallD(i,j,k)=bigD(i2,j2,k2);
//
//            einspline::set(dense,bigD.data());
//            einspline::set(coarse,smallD.data());
//
//            bspline->init_spline(dense,coarse,ival);
//          }
//
//          free(coarse);
//          free(dense);
//          fftw_destroy_plan(FFTplan);
//
//          if(myComm->rank()==0)
//          {
//            hdf_archive h5f;
//            h5f.create(splinefile);
//            einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
//            einspline_engine<typename adoptor_type::SplineType> smalltable(bspline->smallBox);
//            std::string aname("EinsplineOpenAdoptor");
//            h5f.write(aname,"adoptor_name");
//            h5f.write(lower,"lower_bound");
//            h5f.write(upper,"upper_bound");
//            h5f.write(bigtable,"spline_0");
//            h5f.write(smalltable,"spline_1");
//          }
//        }
//
//        app_log() << "TIME READBANDS " << now.elapsed() << std::endl;
//
//        return bspline;
//      }
//    };
}
#endif

