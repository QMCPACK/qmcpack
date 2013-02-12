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
namespace qmcplusplus {
  template<typename SA>
    struct SplineAdoptorReader: BsplineReaderBase
  {
    typedef SA adoptor_type;
    typedef typename adoptor_type::real_type  spline_real_type;
    typedef typename adoptor_type::value_type spline_value_type;
    typedef typename adoptor_type::SplineType SplineType;

    SplineAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
    {}

    SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
    {

      ReportEngine PRE("SplineC2XAdoptorReader","create_spline_set(spin,SPE*)");

      Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
      double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;

      BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
      init(orbitalSet,bspline);


      int N = mybuilder->NumDistinctOrbitals;
      int NumValenceOrbs = mybuilder->NumValenceOrbs;

      bspline->resizeStorage(N,NumValenceOrbs);
      int numOrbs=bspline->getOrbitalSetSize();

      // Read in k-points
      //int numOrbs = bspline->getOrbitalSetSize();
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

      app_log() << "  AdoptorName = " << bspline->AdoptorName << endl;

      if(!bspline->is_complex)
      {//no k-point folding, single special k point (G, L ...)
        int ti=mybuilder->SortBands[0].TwistIndex;
        TinyVector<double,3> twist0 = mybuilder->TwistAngles[mybuilder->SortBands[0].TwistIndex];
        for (int i=0; i<3; i++)
          if (std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
            bspline->HalfG[i] = 1;
          else
            bspline->HalfG[i] = 0;
        app_log() << "  TwistIndex = " << mybuilder->SortBands[0].TwistIndex << " TwistAngle " << twist0 << endl;
        app_log() <<"   HalfG = " << bspline->HalfG << endl;
      }

      // First, check to see if we have already read this in
      string H5FileName(mybuilder->H5FileName);
      H5OrbSet set(H5FileName, spin, N);

      bool havePsig=mybuilder->ReadGvectors_ESHDF();

      TinyVector<int,3> MeshSize=mybuilder->MeshSize;
      app_log() << "MeshSize = (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";

      int nx, ny, nz, bi, ti;
      nx=MeshSize[0]; ny=MeshSize[1]; nz=MeshSize[2];

      Ugrid xyz_grid[3];
      typename adoptor_type::BCType xyz_bc[3];
      xyz_grid[0].start = 0.0;  xyz_grid[0].end = 1.0;  xyz_grid[0].num = nx;
      xyz_grid[1].start = 0.0;  xyz_grid[1].end = 1.0;  xyz_grid[1].num = ny;
      xyz_grid[2].start = 0.0;  xyz_grid[2].end = 1.0;  xyz_grid[2].num = nz;

      for(int j=0; j<3; ++j)
      {
        if(bspline->HalfG[j])
        {
          xyz_bc[j].lCode=ANTIPERIODIC; xyz_bc[j].rCode=ANTIPERIODIC;
        }
        else
        {
          xyz_bc[j].lCode=PERIODIC; xyz_bc[j].rCode=PERIODIC;
        }
      }

      bspline->create_spline(xyz_grid,xyz_bc);

      int TwistNum = mybuilder->TwistNum;

      //app_log() << "  TEST TWIST " << TargetPtcl.getTwist() << endl;
      //string splinefile=make_spline_filename(H5FileName,TwistNum,MeshSize);
      string splinefile=make_spline_filename(H5FileName,mybuilder->TileMatrix,spin,TwistNum,MeshSize);

      bool root=(myComm->rank() == 0);

      int foundspline=0;
      Timer now;
      if(root)
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
          einspline_engine<SplineType> bigtable(bspline->MultiSpline);
          foundspline=h5f.read(bigtable,"spline_0");
        }
      }
      myComm->bcast(foundspline);
      t_h5 = now.elapsed();

      if(foundspline)
      {
        app_log() << "Use existing bspline tables in " << splinefile << endl;
        chunked_bcast(myComm, bspline->MultiSpline); 
        t_init+=now.elapsed();
      }
      else
      {

        /** For valence orbitals,
         * - extended orbitals either in G or in R
         * - localized orbitals
         */
        Array<spline_real_type,3> splineData_r(nx,ny,nz),splineData_i(ny,ny,nz);

        if(havePsig)//perform FFT using FFTW
        {
          c_init.restart();
          Array<complex<double>,3> FFTbox;
          FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
          fftw_plan FFTplan = fftw_plan_dft_3d 
            (MeshSize[0], MeshSize[1], MeshSize[2],
             reinterpret_cast<fftw_complex*>(FFTbox.data()),
             reinterpret_cast<fftw_complex*>(FFTbox.data()),
             +1, FFTW_ESTIMATE);

          Vector<complex<double> > cG(mybuilder->MaxNumGvecs);

          //this will be parallelized with OpenMP
          for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
          {
            int ti=SortBands[iorb].TwistIndex;
            get_psi_g(ti,spin,SortBands[iorb].BandIndex,cG);

            c_unpack.restart();
            unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
            t_unpack+= c_unpack.elapsed();

            c_fft.restart();
            fftw_execute (FFTplan);
            t_fft+= c_fft.elapsed();

            c_phase.restart();
            if(bspline->is_complex)
              fix_phase_rotate_c2c(FFTbox,splineData_r, splineData_i,mybuilder->TwistAngles[ti]);
            else
              fix_phase_rotate_c2r(FFTbox,splineData_r, mybuilder->TwistAngles[ti]);
            t_phase+= c_phase.elapsed();

            c_spline.restart();
            bspline->set_spline(ival,splineData_r.data(),splineData_i.data());
            t_spline+= c_spline.elapsed();
          }

          fftw_destroy_plan(FFTplan);
          t_init+=c_init.elapsed();
        }
        else
        {
          Array<complex<double>,3> rawData(nx,ny,nz);
          //this will be parallelized with OpenMP
          for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
          {
            //check dimension
            if(root)
            {
              string path=psi_r_path(SortBands[iorb].TwistIndex,spin,SortBands[iorb].BandIndex);
              HDFAttribIO<Array<complex<double>,3> >  h_splineData(rawData);
              h_splineData.read(mybuilder->H5FileID, path.c_str());
              simd::copy(splineData_r.data(),splineData_i.data(),rawData.data(),rawData.size());
            }
            myComm->bcast(splineData_r);
            myComm->bcast(splineData_i);
            bspline->set_spline(ival,splineData_r.data(),splineData_i.data());
          }
        }

        if(root)
        {
          hdf_archive h5f;
          h5f.create(splinefile);
          einspline_engine<SplineType> bigtable(bspline->MultiSpline);
          h5f.write(bspline->AdoptorName,"adoptor_name");
          h5f.write(bigtable,"spline_0");
        }
      }

      app_log() << "    READBANDS::PREP   = " << t_prep << endl;
      app_log() << "    READBANDS::H5     = " << t_h5 << endl;
      app_log() << "    READBANDS::UNPACK = " << t_unpack << endl;
      app_log() << "    READBANDS::FFT    = " << t_fft << endl;
      app_log() << "    READBANDS::PHASE  = " << t_phase << endl;
      app_log() << "    READBANDS::SPLINE = " << t_spline << endl;
      app_log() << "    READBANDS::SUM    = " << t_init << endl;

      return bspline;
    }
  };
}
#endif
