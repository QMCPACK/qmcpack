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

#ifndef QMCPLUSPLUS_EINSPLINE_BUILDER_ESHDF_BIG_H
#define QMCPLUSPLUS_EINSPLINE_BUILDER_ESHDF_BIG_H
#include "Utilities/ProgressReportEngine.h"
#include "Message/CommOperators.h"
#include <QMCWaveFunctions/einspline_helper.hpp>
#include <fftw3.h>

namespace qmcplusplus {

  template<typename SPE>
  void EinsplineSetBuilder::ReadBands_ESHDF_Big(int spin, SPE* orbitalSet)
  {
    ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF_Big(spin,SPE*,need2convert)");

    bool root = myComm->rank()==0;
    // bcast other stuff
    myComm->bcast (NumDistinctOrbitals);
    myComm->bcast (NumValenceOrbs);
    myComm->bcast (NumCoreOrbs);
    int N = NumDistinctOrbitals;

    typedef typename SPE::DataType spline_data_type;
    //typedef double spline_data_type;
    orbitalSet->PrimLattice=Lattice;
    orbitalSet->resizeStorage(N,NumValenceOrbs);

    int numOrbs = orbitalSet->getOrbitalSetSize();

    //int num = 0;
    //vector<int> twist_id(N);
    ////twist-related things can be completely ignored but will keep them for now
    //if (root) {
    //  for (int iorb=0; iorb<N; iorb++) {
    //    int ti = SortBands[iorb].TwistIndex;
    //    twist_id[iorb]=ti;
    //    PosType twist  = TwistAngles[ti];
    //    orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
    //    orbitalSet->MakeTwoCopies[iorb] = 
    //      (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
    //    num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
    //  }
    //  PosType twist0 = TwistAngles[SortBands[0].TwistIndex];
    //  for (int i=0; i<OHMMS_DIM; i++)
    //    if (std::fabs(std::fabs(twist0[i]) - 0.5) < 1.0e-8)
    //      orbitalSet->HalfG[i] = 1;
    //    else
    //      orbitalSet->HalfG[i] = 0;
    //  //EinsplineSetBuilder::RotateBands_ESHDF(spin, orbitalSet);
    //}
    //myComm->bcast(orbitalSet->kPoints);
    //myComm->bcast(orbitalSet->MakeTwoCopies);
    //myComm->bcast(orbitalSet->HalfG);
    //myComm->bcast(twist_id);

    // First, check to see if we have already read this in
    H5OrbSet set(H5FileName, spin, N);
    
    bool havePsir=!ReadGvectors_ESHDF();

    bool fullgrid=false;

    TinyVector<int,3> coarse_mesh(MeshSize[0]/2,MeshSize[1]/2,MeshSize[2]/2);
    TinyVector<double,3> start(0.0);
    TinyVector<double,3> end(1.0);
    for(int i=0; i<3; ++i) end(i)=SuperLattice(i,i);

    //create dense and coarse UBspline_3d_d 
    UBspline_3d_d* dense=0;
    dense=einspline::create(dense,start,end,MeshSize,PERIODIC);
    UBspline_3d_d* coarse=0;
    coarse=einspline::create(coarse,start,end,coarse_mesh,PERIODIC);

    //determine the bonding box
    TinyVector<double,3> lower=end;
    TinyVector<double,3> upper=start;
    for(int i=0; i<IonPos.size(); ++i)
    {
      for(int j=0; j<3; ++j)
      {
        lower[j]=std::min(IonPos[i][j],lower[j]);
        upper[j]=std::max(IonPos[i][j],upper[j]);
      }
    }
    //add 2 bohr
    for(int j=0; j<3; ++j)
    {
      lower[j]=std::max(lower[j]-2.0,start[j]);
      upper[j]=std::min(upper[j]+2.0,end[j]);
    }

    orbitalSet->create_spline(MeshSize,NumValenceOrbs,fullgrid);

    app_log() << "  Adding a small box" 
      << "\n  LowerBound " << lower 
      << "\n  UpperBound " << upper << endl;
    orbitalSet->add_box(dense,lower,upper);

    Timer now;
    {
      Array<ComplexType,3> FFTbox;
      FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
      fftw_plan FFTplan = fftw_plan_dft_3d 
        (MeshSize[0], MeshSize[1], MeshSize[2],
         reinterpret_cast<fftw_complex*>(FFTbox.data()),
         reinterpret_cast<fftw_complex*>(FFTbox.data()),
         +1, FFTW_ESTIMATE);

      Array<double,3> bigD(MeshSize[0],MeshSize[1],MeshSize[2]);
      Array<double,3> smallD(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);

      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        Vector<complex<double> > cG;
        int ncg=0;
        int ti=0;
        //int ti=SortBands[iorb].TwistIndex;
        if(root)
        {
          ostringstream path;
          path << "/electrons/kpoint_" << ti    //SortBands[iorb].TwistIndex 
            << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
          HDFAttribIO<Vector<complex<double> > >  h_cG(cG);
          h_cG.read (H5FileID, path.str().c_str());
          ncg=cG.size();
        }
        myComm->bcast(ncg);
        if(ncg != Gvecs[ti].size())
        {
          APP_ABORT("Failed : ncg != Gvecs[ti].size()");
        }

        if(!root) cG.resize(ncg);
        myComm->bcast(cG);

        unpack4fftw(cG,Gvecs[ti],MeshSize,FFTbox);
        fftw_execute (FFTplan);
        fix_phase_rotate_c2r(FFTbox,bigD,TwistAngles[ti]);

        for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=2)
          for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=2)
            for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=2)
              smallD(i,j,k)=bigD(i2,j2,k2);

        einspline::set(dense,bigD.data());
        einspline::set(coarse,smallD.data());

        orbitalSet->init_spline(dense,coarse,ival);
      }

      free(coarse);
      free(dense);
      fftw_destroy_plan(FFTplan);
    }

    app_log() << "TIME READBANDS " << now.elapsed() << endl;
  }
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: EinsplineSetBuilderESHDF.cpp 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
