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
#include "Utilities/ProgressReportEngine.h"
#include "Message/CommOperators.h"
#include <QMCWaveFunctions/einspline_helper.hpp>
#include <spline/einspline_util.hpp>
#include <fftw3.h>

namespace qmcplusplus {

  /** adoptor class to match ST real spline with TT real SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineR2RAdoptor
    {
      typedef typename einspline_traits<ST,D>::SplineType SplineType;
      typedef typename einspline_traits<ST,D>::BCType     BCType;
      typedef typename einspline_traits<ST,D>::DataType   DataType;
      typedef CrystalLattice<ST,D>                        UnitCellType;
      typedef TinyVector<ST,D>                            PointType;
      typedef typename OrbitalSetTraits<ST>::ValueVector_t      StorageValueVector_t;
      typedef typename OrbitalSetTraits<ST>::GradVector_t       StorageGradVector_t;
      typedef typename OrbitalSetTraits<ST>::HessVector_t       StorageHessVector_t;
      typedef typename OrbitalSetTraits<ST>::GradHessVector_t   StorageGradHessVector_t;

      SplineType *MultiSpline;
      Tensor<ST,D> GGt;
      TinyVector<int,D> HalfG;
      UnitCellType PrimLattice;

      //these are not needed for R2R adoptor and will be removed
      UnitCellType SuperLattice;
      vector<TinyVector<ST,D> > kPoints;
      vector<bool>      MakeTwoCopies;

      // Temporary storage for Eispline calls
      StorageValueVector_t myV, myL;
      StorageGradVector_t      myG;
      StorageHessVector_t      myH;
      StorageGradHessVector_t  myGH;

      void resizeStorage(int n, int nv)
      {
        kPoints.resize(n);
        MakeTwoCopies.resize(n);
        myV.resize(n);
        myL.resize(n);
        myG.resize(n);
        myH.resize(n);
        myGH.resize(n);
      }

      void create_spline(TinyVector<int,D>& mesh, int n)
      {
        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
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

      ///** return sign */
      inline void convertPos(const PointType& r, PointType& ru, int& sign)
      {
        ru=PrimLattice.toUnit(r);
        sign=0;
        for (int i=0; i<D; i++) {
          ST img = std::floor(ru[i]);
          ru[i] -= img;
          sign += HalfG[i] * (int)img;
        }
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          int phase;
          TinyVector<ST,D> ru;
          convertPos(r,ru,phase);

          einspline::evaluate(MultiSpline,ru,myV);
          if (phase & 1) 
            for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(-myV[j]);
          else
            for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(myV[j]);
        }

      template<typename VV, typename GV>
        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          int phase;
          TinyVector<ST,D> ru;
          convertPos(r,ru,phase);

          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);

          const int N=psi.size();
          const Tensor<ST,D> gConv(PrimLattice.G);
          if (phase & 1) 
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

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {}
    };

  template<typename T>
    struct SplineR2RAdoptorReader: BsplineReaderBase
  {
    typedef SplineR2RAdoptor<T,double,3> adoptor_type;
    SplineR2RAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
      {}

    SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
    {
      ReportEngine PRE("SplineR2RAdoptorReader","create_spline_set(int,EinsplineSet*)");

      typedef T spline_data_type;

      BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;

      init(orbitalSet,bspline);

      int N=bspline->getOrbitalSetSize();
      bspline->resizeStorage(N,N);

      int num = 0;
      bool root=myComm->rank()==0;

      int ti=mybuilder->SortBands[0].TwistIndex;
      TinyVector<double,3> twist0 = mybuilder->TwistAngles[mybuilder->SortBands[0].TwistIndex];
      for (int i=0; i<3; i++)
        if (std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
          bspline->HalfG[i] = 1;
        else
          bspline->HalfG[i] = 0;

      app_log() << "  TwistIndex = " << mybuilder->SortBands[0].TwistIndex << " TwistAngle " << twist0 << endl;
      app_log() <<"   HalfG = " << bspline->HalfG << endl;

      string H5FileName(mybuilder->H5FileName);

      H5OrbSet set(H5FileName, spin, N);

      bool havePsir=!(mybuilder->ReadGvectors_ESHDF());

      TinyVector<int,3> MeshSize=mybuilder->MeshSize;

      app_log() << "MeshSize = (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";

      int nx, ny, nz;
      nx=MeshSize[0]; ny=MeshSize[1]; nz=MeshSize[2];

      bspline->allocate(MeshSize,N);

      string splinefile=make_spline_filename(H5FileName,spin,mybuilder->TwistNum,MeshSize);
      int foundspline=0;
      Timer now;
      if(root)
      {
        hdf_archive h5f;
        foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
        if(foundspline)
        {
          einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
          foundspline=h5f.read(bigtable,"spline_0");
        }
      }

      myComm->bcast(foundspline);
      if(foundspline)
      {
        app_log() << "Use existing bspline tables in " << splinefile << endl;
        chunked_bcast(myComm, bspline->MultiSpline); 
      }
      else
      {
        Array<complex<double>,3> FFTbox;
        FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
        fftw_plan FFTplan = fftw_plan_dft_3d 
          (MeshSize[0], MeshSize[1], MeshSize[2],
           reinterpret_cast<fftw_complex*>(FFTbox.data()),
           reinterpret_cast<fftw_complex*>(FFTbox.data()),
           +1, FFTW_ESTIMATE);
        //Array<spline_data_type,3> splineData(nx,ny,nz);
        Array<double,3> splineData(nx,ny,nz);

        int ncg=mybuilder->Gvecs[ti].size();
        Vector<complex<double> > cG(ncg);

        for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
        {
          int ti=mybuilder->SortBands[iorb].TwistIndex;

          get_psi_g(ti,spin,mybuilder->SortBands[iorb].BandIndex,cG);
          unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
          fftw_execute (FFTplan);
          fix_phase_rotate_c2r(FFTbox,splineData,twist0);
          einspline::set(bspline->MultiSpline, ival,splineData.data());
        }

        fftw_destroy_plan(FFTplan);

        if(root)
        {
          hdf_archive h5f;
          h5f.create(splinefile);
          einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
          string aname("EinsplineR2RAdoptor");
          h5f.write(aname,"bspline_type");
          h5f.write(bigtable,"spline_0");
        }
      }

      app_log() << "TIME READBANDS " << now.elapsed() << endl;

      return bspline;
    }

    //ExtendedMap_d[set] = orbitalSet->MultiSpline;
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: EinsplineSetBuilderESHDF.cpp 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
