/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EINSPLINE_SPLINE_MIXED_ADOPTOR_READER_H
#define QMCPLUSPLUS_EINSPLINE_SPLINE_MIXED_ADOPTOR_READER_H

namespace qmcplusplus {

  template<typename SA>
    struct SplineMixedAdoptorReader: BsplineReaderBase
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
        init(orbitalSet,bspline);
        int N=bspline->getOrbitalSetSize();

        bspline->resizeStorage(N,N);

        TinyVector<int,3> bconds=mybuilder->TargetPtcl.Lattice.BoxBConds;

        bspline->HalfG=0;
        if(!bspline->is_complex)
        {//no k-point folding, single special k point (G, L ...)
          int ti=mybuilder->SortBands[0].TwistIndex;
          TinyVector<double,3> twist0 = mybuilder->TwistAngles[mybuilder->SortBands[0].TwistIndex];
          for (int i=0; i<3; i++)
            if (bconds[i] && std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
              bspline->HalfG[i] = 1;
            else
              bspline->HalfG[i] = 0;
          app_log() << "  TwistIndex = " << mybuilder->SortBands[0].TwistIndex << " TwistAngle " << twist0 << endl;
          app_log() <<"   HalfG = " << bspline->HalfG << endl;
        }

        string H5FileName(mybuilder->H5FileName);
        H5OrbSet set(H5FileName, spin, N);
    
        bool havePsir=!(mybuilder->ReadGvectors_ESHDF());

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
        if(use_cartesian) end=maxL;

        //create dense and coarse UBspline_3d_d 
        UBspline_3d_d* dense=0;
        dense=einspline::create(dense,start,end,MeshSize,bspline->HalfG);//PERIODIC);
        UBspline_3d_d* coarse=0;
        coarse=einspline::create(coarse,start,end,coarse_mesh,bspline->HalfG);//PERIODIC);

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
          for(int j=0;j<3; ++j)
          {
            if(bconds[j])
            {
              lower[j]=0.0; upper[j]=1.0;
            }
            else
            {
              lower[j]=std::max(lower[j]-buffer/maxL[j],start[j]);
              upper[j]=std::min(upper[j]+buffer/maxL[j],end[j]);
            }
          }
        }

        //bspline->create_spline(MeshSize,N,fullgrid);
        bspline->create_spline(MeshSize,coarse_mesh,N);

        app_log() << "  Original Mesh " << MeshSize << endl;
        app_log() << "  Coarse Mesh " << coarse_mesh << endl;

        if(use_cartesian)
          app_log() << "  Using Cartesian grids for open systems. " << endl;
        else
          app_log() << "  Using Primitive-cell grids for open systems. " << endl;

        app_log() << "  Using buffer layer for the small box= " << buffer << " bohr " << endl;
        app_log() << "  Adding a small box" << "\n  LowerBound " << lower << "\n  UpperBound " << upper << endl;

        app_log().flush();

        bspline->add_box(dense,lower,upper);

        int foundspline=0;

        string splinefile=make_spline_filename(H5FileName,spin,mybuilder->TwistNum,MeshSize);

        Timer now;

        if(myComm->rank()==0)
        {
          hdf_archive h5f;
          foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
          if(foundspline)
          {
            string aname("none");
            foundspline = h5f.read(aname,"adoptor_name");
            foundspline = (aname.find(bspline->KeyWord) != std::string::npos);
          }
          if(foundspline)
          {
            TinyVector<double,3> lower_in(end);
            TinyVector<double,3> upper_in(start);
            h5f.read(lower_in,"lower_bound");
            h5f.read(upper_in,"upper_bound");

            lower_in-=lower; upper_in-=upper;
            if(dot(lower_in,lower_in)<1e-12 &&dot(upper_in,upper_in)<1e-12)
            {
              einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
              einspline_engine<typename adoptor_type::SplineType> smalltable(bspline->smallBox);
              foundspline=h5f.read(bigtable,"spline_0");
              foundspline=h5f.read(smalltable,"spline_1");
            }
            else
            {
              app_log() << "  The upper/lower bound of the input is different from the current value."<< endl;
              foundspline=0;
            }
          }
        }

        myComm->bcast(foundspline);

        if(foundspline)
        {
          app_log() << "Use existing bspline tables in " << splinefile << endl;
          chunked_bcast(myComm, bspline->MultiSpline); 
          chunked_bcast(myComm, bspline->smallBox); 
        }
        else
        {
          app_log() << "Perform FFT+spline and dump bspline tables to " << splinefile << endl;
          Array<complex<double>,3> FFTbox;
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

          //TinyVector<double,3> TwistAngle(0.0,0.0,0.0);
          int ncg=mybuilder->Gvecs[0].size();
          Vector<complex<double> > cG(ncg);

          const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
          for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
          {
            int ti=SortBands[iorb].TwistIndex;
            get_psi_g(ti,spin,mybuilder->SortBands[iorb].BandIndex,cG);
            unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
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

            if(bspline->is_complex)
            {
              einspline::set(dense,bigD_r.data());
              einspline::set(coarse,smallD_r.data());
              bspline->set_spline(dense,coarse,2*ival);

              for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=coarse_stride[0])
                for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=coarse_stride[1])
                  for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=coarse_stride[2])
                    smallD_i(i,j,k)=bigD_i(i2,j2,k2);

              einspline::set(dense,bigD_i.data());
              einspline::set(coarse,smallD_i.data());
              bspline->set_spline(dense,coarse,2*ival+1);
            }
            else
            {
              einspline::set(dense,bigD_r.data());
              einspline::set(coarse,smallD_r.data());
              bspline->set_spline(dense,coarse,ival);
            }
          }

          free(coarse);
          free(dense);
          fftw_destroy_plan(FFTplan);

          if(myComm->rank()==0)
          {
            hdf_archive h5f;
            h5f.create(splinefile);
            einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
            einspline_engine<typename adoptor_type::SplineType> smalltable(bspline->smallBox);
            h5f.write(bspline->AdoptorName,"adoptor_name");
            h5f.write(lower,"lower_bound");
            h5f.write(upper,"upper_bound");
            h5f.write(bigtable,"spline_0");
            h5f.write(smalltable,"spline_1");
          }
        }

        app_log() << "TIME READBANDS " << now.elapsed() << endl;

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
//        string H5FileName(mybuilder->H5FileName);
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
//        app_log() << "Original Mesh " << MeshSize << endl;
//        app_log() << "Coarse Mesh " << coarse_mesh << endl;
//
//        app_log() << "  Adding a small box" << "\n  LowerBound " << lower << "\n  UpperBound " << upper << endl;
//
//        app_log().flush();
//
//        bspline->add_box(dense,lower,upper);
//
//        int foundspline=0;
//
//        string splinefile=make_spline_filename(H5FileName,spin,0,MeshSize);
//
//        Timer now;
//
//        if(myComm->rank()==0)
//        {
//          hdf_archive h5f;
//          foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
//          if(foundspline)
//          {
//            string aname("none");
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
//              app_log() << "  The upper/lower bound of the input is different from the current value."<< endl;
//              foundspline=0;
//            }
//          }
//        }
//
//        myComm->bcast(foundspline);
//
//        if(foundspline)
//        {
//          app_log() << "Use existing bspline tables in " << splinefile << endl;
//          chunked_bcast(myComm, bspline->MultiSpline); 
//          chunked_bcast(myComm, bspline->smallBox); 
//        }
//        else
//        {
//          app_log() << "Perform FFT+spline and dump bspline tables to " << splinefile << endl;
//          Array<complex<double>,3> FFTbox;
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
//          Vector<complex<double> > cG(ncg);
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
//            string aname("EinsplineOpenAdoptor");
//            h5f.write(aname,"adoptor_name");
//            h5f.write(lower,"lower_bound");
//            h5f.write(upper,"upper_bound");
//            h5f.write(bigtable,"spline_0");
//            h5f.write(smalltable,"spline_1");
//          }
//        }
//
//        app_log() << "TIME READBANDS " << now.elapsed() << endl;
//
//        return bspline;
//      }
//    };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: EinsplineSetBuilderESHDF.cpp 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
