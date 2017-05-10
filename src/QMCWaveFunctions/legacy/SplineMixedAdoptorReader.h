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
    
    
/** @file SplineMixedAdoptorReader.h
 */

#ifndef QMCPLUSPLUS_EINSPLINE_SPLINE_MIXED_ADOPTOR_READER_H
#define QMCPLUSPLUS_EINSPLINE_SPLINE_MIXED_ADOPTOR_READER_H

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

  SplineMixedAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
  {}

  SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
  {
    ReportEngine PRE("SplineOpenAdoptorReader","create_spline_set(int, EinsplineSet*)");
    //typedef typename adoptor_type::real_type spline_data_type;
    BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
    app_log() << "  AdoptorName = " << bspline->AdoptorName << std::endl;
    if(bspline->is_complex)
      app_log() << "  Using complex einspline table" << std::endl;
    else
      app_log() << "  Using real einspline table" << std::endl;
    check_twists(orbitalSet,bspline);
    Ugrid xyz_grid[3];
    typename adoptor_type::BCType xyz_bc[3];
    bool havePsig=set_grid(bspline->HalfG,xyz_grid, xyz_bc);
    if(!havePsig)
    {
      APP_ABORT("Need psi_g with truncate=\"yes\"");
    }
    TinyVector<int,3> bconds=mybuilder->TargetPtcl.Lattice.BoxBConds;
    bool use_cartesian= (mybuilder->TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN);
    TinyVector<int,3> MeshSize=mybuilder->MeshSize;
    TinyVector<int,3> coarse_mesh;
    TinyVector<int,3> coarse_stride;
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
    //create dense and coarse UBspline_3d_d
    UBspline_3d_d* dense_r=0;
    UBspline_3d_d* dense_i=0;
    dense_r=einspline::create(dense_r,start,end,MeshSize,bspline->HalfG);//PERIODIC);
    UBspline_3d_d* coarse_r=0;
    UBspline_3d_d* coarse_i=0;
    coarse_r=einspline::create(coarse_r,start,end,coarse_mesh,bspline->HalfG);//PERIODIC);
    if(bspline->is_complex)
    {
      dense_i=einspline::create(dense_i,start,end,MeshSize,bspline->HalfG);//PERIODIC);
      coarse_i=einspline::create(coarse_i,start,end,coarse_mesh,bspline->HalfG);//PERIODIC);
    }
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
    //bspline->create_spline(MeshSize,N,fullgrid);
    int N = mybuilder->NumDistinctOrbitals;
    bspline->create_spline(MeshSize,coarse_mesh,N,N);
    app_log() << "  Original Mesh " << MeshSize << std::endl;
    app_log() << "  Coarse Mesh " << coarse_mesh << std::endl;
    if(use_cartesian)
      app_log() << "  Using Cartesian grids for open systems. " << std::endl;
    else
      app_log() << "  Using Primitive-cell grids for mixed systems. " << std::endl;
    app_log() << "  Using buffer layer for the small box= " << buffer << " bohr " << std::endl;
    app_log() << "  Adding a small box" << "\n  LowerBound " << lower << "\n  UpperBound " << upper << std::endl;
    app_log().flush();
    bspline->add_box(dense_r,lower,upper);
    int foundspline=0;
    std::string splinefile=make_spline_filename(mybuilder->H5FileName,spin,mybuilder->TwistNum,MeshSize);
    Timer now;
    if(myComm->rank()==0)
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
        foundspline=bspline->read_splines(h5f);
    }
    myComm->bcast(foundspline);
    if(foundspline)
    {
      app_log() << "Use existing bspline tables in " << splinefile << std::endl;
      chunked_bcast(myComm, bspline->MultiSpline);
      chunked_bcast(myComm, bspline->smallBox);
    }
    else
    {
      app_log() << "Perform FFT+spline and dump bspline tables to " << splinefile << std::endl;
      Array<std::complex<double>,3> FFTbox;
      FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
      fftw_plan FFTplan = fftw_plan_dft_3d
                          (MeshSize[0], MeshSize[1], MeshSize[2],
                           reinterpret_cast<fftw_complex*>(FFTbox.data()),
                           reinterpret_cast<fftw_complex*>(FFTbox.data()),
                           +1, FFTW_ESTIMATE);
      Array<double,3> bigD_r(MeshSize[0],MeshSize[1],MeshSize[2]);
      Array<double,3> smallD_r(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      Array<double,3> bigD_i, smallD_i;
      if(bspline->is_complex)
      {
        bigD_i.resize(MeshSize[0],MeshSize[1],MeshSize[2]);
        smallD_i.resize(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      }
      Vector<std::complex<double> > cG(mybuilder->MaxNumGvecs);
      const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        int ti=SortBands[iorb].TwistIndex;
        get_psi_g(ti,spin,mybuilder->SortBands[iorb].BandIndex,cG);
        unpack4fftw(cG,mybuilder->Gvecs[0],MeshSize,FFTbox);
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
        einspline::set(dense_r,bigD_r.data());
        einspline::set(coarse_r,smallD_r.data());
        if(bspline->is_complex)
        {
          for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=coarse_stride[0])
            for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=coarse_stride[1])
              for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=coarse_stride[2])
                smallD_i(i,j,k)=bigD_i(i2,j2,k2);
          einspline::set(dense_i,bigD_i.data());
          einspline::set(coarse_i,smallD_i.data());
        }
        bspline->set_spline(coarse_r,coarse_i,ti,ival,0);//level=0, full grid
        bspline->set_spline(dense_r, dense_i, ti,ival,1);//level=1, box grid
      }
      free(coarse_r);
      free(dense_r);
      if(bspline->is_complex)
      {
        free(coarse_i);
        free(dense_i);
      }
      fftw_destroy_plan(FFTplan);
      if(qmc_common.save_wfs && myComm->rank()==0)
      {
        hdf_archive h5f;
        h5f.create(splinefile);
        h5f.write(bspline->AdoptorName,"adoptor_name");
        int sizeD=sizeof(typename adoptor_type::DataType);
        h5f.write(sizeD,"sizeof");
        bspline->write_splines(h5f);
      }
    }
    app_log() << "TIME READBANDS " << now.elapsed() << std::endl;
    return bspline;
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

