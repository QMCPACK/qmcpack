/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EINSPLINE_R2RADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_R2RADOPTOR_H

namespace qmcplusplus {

  /** adoptor class to match ST real spline with TT real SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineR2RAdoptor: public SplineAdoptorBase<ST,D>
    {
      typedef typename einspline_traits<ST,D>::SplineType SplineType;
      typedef typename einspline_traits<ST,D>::BCType     BCType;
      typedef typename SplineAdoptorBase<ST,D>::PointType PointType;

      using SplineAdoptorBase<ST,D>::HalfG;
      using SplineAdoptorBase<ST,D>::GGt;
      using SplineAdoptorBase<ST,D>::PrimLattice;

      using SplineAdoptorBase<ST,D>::myV;
      using SplineAdoptorBase<ST,D>::myL;
      using SplineAdoptorBase<ST,D>::myG;
      using SplineAdoptorBase<ST,D>::myH;
      using SplineAdoptorBase<ST,D>::myGH;
      SplineType *MultiSpline;

      SplineR2RAdoptor(): MultiSpline(0)
      {
        this->is_complex=false;
        this->AdoptorName="SplineR2RAdoptor";
        this->KeyWord="R2R";
      }

      void resizeStorage(int n, int nv)
      {
        SplineAdoptorBase<ST,D>::init_base(n);
        myV.resize(n);
        myL.resize(n);
        myG.resize(n);
        myH.resize(n);
        myGH.resize(n);
      }

    template<typename GT, typename BCT>
      void create_spline(GT& xyz_g, BCT& xyz_bc)
      {
        cout << this->AdoptorName << ":create_spline(xyz_g,xyz_bc)" << endl;
        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
        MultiSpline=einspline::create(MultiSpline,xyz_g,xyz_bc,myV.size());
      }

      void create_spline(TinyVector<int,D>& mesh, int n)
      {
        Ugrid xyz_grid[D];
        BCType xyz_bc[D];
        for(int i=0; i<D; ++i)
        {
          xyz_grid[i].start = 0.0;  xyz_grid[i].end = 1.0;  
          xyz_grid[i].num = mesh[i];
          xyz_bc[i].lCode=xyz_bc[i].rCode=(HalfG[i])? ANTIPERIODIC:PERIODIC;
        }

        SplineType* dummy=0;
        MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);
      }

      inline bool isready()
      {
        return true;
      }

      inline void set_spline(int ival, ST* psi_r, ST* psi_i)
      {
        einspline::set(MultiSpline, ival,psi_r);
      }
      /** convert postion in PrimLattice unit and return sign */
      inline int convertPos(const PointType& r, PointType& ru)
      {
        ru=PrimLattice.toUnit(r);
        int bc_sign=0;
        for (int i=0; i<D; i++) {
          ST img = std::floor(ru[i]);
          ru[i] -= img;
          bc_sign += HalfG[i] * (int)img;
        }
        return bc_sign;
      }

      /** assign myV to psi
       */
      template<typename VV>
        inline void assign_v(const PointType& r, int bc_sign, VV& psi) 
        {
          if (bc_sign & 1) 
            for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(-myV[j]);
          else
            for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(myV[j]);
        }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          PointType ru;
          int bc_sign=convertPos(r,ru);
          einspline::evaluate(MultiSpline,ru,myV);
          assign_v(r,bc_sign,psi);
        }

      /** assign internal data to psi's
       */
      template<typename VV, typename GV>
        inline void assign_vgl(const PointType& r, int bc_sign, VV& psi, GV& dpsi, VV& d2psi) 
        {
          const int N=psi.size();
          const Tensor<ST,D> gConv(PrimLattice.G);
          if (bc_sign & 1) 
          {
            const ST minus_one=-1.0;
            for(int j=0; j<N; ++j) psi[j]=-myV[j];
            for(int j=0; j<N; ++j) dpsi[j]=minus_one*dot(gConv,myG[j]);
            for(int j=0; j<N; ++j) d2psi[j]=-trace(myH[j],GGt);
          }
          else
          {
            for(int j=0; j<N; ++j) psi[j]=myV[j];
            for(int j=0; j<N; ++j) dpsi[j]=dot(gConv,myG[j]);
            for(int j=0; j<N; ++j) d2psi[j]=trace(myH[j],GGt);
          }
        }

      template<typename VV, typename GV>
        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          PointType ru;
          int bc_sign=convertPos(r,ru);

          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
          assign_vgl(r,bc_sign,psi,dpsi,d2psi);
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {}
    };

//  template<typename T>
//    struct SplineR2RAdoptorReader: BsplineReaderBase
//  {
//    typedef SplineR2RAdoptor<T,double,3> adoptor_type;
//    SplineR2RAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
//      {}
//
//    SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
//    {
//      ReportEngine PRE("SplineR2RAdoptorReader","create_spline_set(int,EinsplineSet*)");
//
//      typedef T spline_data_type;
//
//      BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
//
//      init(orbitalSet,bspline);
//
//      int N=bspline->getOrbitalSetSize();
//      bspline->resizeStorage(N,N);
//
//      int num = 0;
//      bool root=myComm->rank()==0;
//
//      int ti=mybuilder->SortBands[0].TwistIndex;
//      TinyVector<double,3> twist0 = mybuilder->TwistAngles[mybuilder->SortBands[0].TwistIndex];
//      for (int i=0; i<3; i++)
//        if (std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
//          bspline->HalfG[i] = 1;
//        else
//          bspline->HalfG[i] = 0;
//
//      app_log() << "  TwistIndex = " << mybuilder->SortBands[0].TwistIndex << " TwistAngle " << twist0 << endl;
//      app_log() <<"   HalfG = " << bspline->HalfG << endl;
//
//      string H5FileName(mybuilder->H5FileName);
//
//      H5OrbSet set(H5FileName, spin, N);
//
//      bool havePsir=!(mybuilder->ReadGvectors_ESHDF());
//
//      TinyVector<int,3> MeshSize=mybuilder->MeshSize;
//
//      app_log() << "MeshSize = (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
//
//      int nx, ny, nz;
//      nx=MeshSize[0]; ny=MeshSize[1]; nz=MeshSize[2];
//
//      bspline->allocate(MeshSize,N);
//
//      string splinefile=make_spline_filename(H5FileName,spin,mybuilder->TwistNum,MeshSize);
//      int foundspline=0;
//      Timer now;
//      if(root)
//      {
//        hdf_archive h5f;
//        foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
//        if(foundspline)
//        {
//          string aname("none");
//          foundspline = h5f.read(aname,"adoptor_name");
//          foundspline = (aname.find(bspline->KeyWord) != std::string::npos);
//        }
//        if(foundspline)
//        {
//          einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
//          foundspline=h5f.read(bigtable,"spline_0");
//        }
//      }
//
//      myComm->bcast(foundspline);
//      if(foundspline)
//      {
//        app_log() << "Use existing bspline tables in " << splinefile << endl;
//        chunked_bcast(myComm, bspline->MultiSpline); 
//      }
//      else
//      {
//        Array<complex<double>,3> FFTbox;
//        FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
//        fftw_plan FFTplan = fftw_plan_dft_3d 
//          (MeshSize[0], MeshSize[1], MeshSize[2],
//           reinterpret_cast<fftw_complex*>(FFTbox.data()),
//           reinterpret_cast<fftw_complex*>(FFTbox.data()),
//           +1, FFTW_ESTIMATE);
//        //Array<spline_data_type,3> splineData(nx,ny,nz);
//        Array<spline_data_type,3> splineData(nx,ny,nz);
//
//        int ncg=mybuilder->Gvecs[ti].size();
//        Vector<complex<double> > cG(ncg);
//
//        for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
//        {
//          int ti=mybuilder->SortBands[iorb].TwistIndex;
//
//          get_psi_g(ti,spin,mybuilder->SortBands[iorb].BandIndex,cG);
//          unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
//          fftw_execute (FFTplan);
//          fix_phase_rotate_c2r(FFTbox,splineData,twist0);
//          bspline->set_spline(ival,splineData.data(),0);
//          //einspline::set(bspline->MultiSpline, ival,splineData.data());
//        }
//
//        fftw_destroy_plan(FFTplan);
//
//        if(root)
//        {
//          hdf_archive h5f;
//          h5f.create(splinefile);
//          einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
//          string aname("EinsplineR2RAdoptor");
//          h5f.write(aname,"adoptor_name");
//          h5f.write(bigtable,"spline_0");
//        }
//      }
//
//      app_log() << "TIME READBANDS " << now.elapsed() << endl;
//
//      return bspline;
//    }
//
//    //ExtendedMap_d[set] = orbitalSet->MultiSpline;
//  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: EinsplineSetBuilderESHDF.cpp 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
