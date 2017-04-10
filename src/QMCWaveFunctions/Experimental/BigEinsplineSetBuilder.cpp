//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_EINSPLINE_BUILDER_ESHDF_BIG_H
#define QMCPLUSPLUS_EINSPLINE_BUILDER_ESHDF_BIG_H
#include "Utilities/ProgressReportEngine.h"
#include <QMCWaveFunctions/einspline_helper.hpp>
#include <spline/einspline_util.hpp>
#include <fftw3.h>

namespace qmcplusplus
{

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
  //    if (std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)
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
  for(int i=0; i<3; ++i)
    end(i)=SuperLattice(i,i);
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
            << "\n  UpperBound " << upper << std::endl;
  orbitalSet->add_box(dense,lower,upper);
  int foundspline=0;
  std::string splinefile=make_spline_filename(H5FileName,spin,0,MeshSize);
  Timer now;
  if(root)
  {
    hdf_archive h5f;
    foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
    if(foundspline)
    {
      TinyVector<double,3> lower_in(end);
      TinyVector<double,3> upper_in(start);
      h5f.read(lower_in,"lower_bound");
      h5f.read(upper_in,"upper_bound");
      lower_in-=lower;
      upper_in-=upper;
      if(dot(lower_in,lower_in)<1e-12 &&dot(upper_in,upper_in)<1e-12)
      {
        einspline_engine<typename SPE::SplineType> bigtable(orbitalSet->MultiSpline);
        einspline_engine<typename SPE::SplineType> smalltable(orbitalSet->smallBox);
        foundspline=h5f.read(bigtable,"spline_0");
        foundspline=h5f.read(smalltable,"spline_1");
      }
      else
      {
        app_log() << "  The upper/lower bound of the input is different from the current value."<< std::endl;
        foundspline=0;
      }
    }
  }
  myComm->bcast(foundspline);
  if(foundspline)
  {
    app_log() << "Use existing bspline tables in " << splinefile << std::endl;
    chunked_bcast(myComm, orbitalSet->MultiSpline);
    chunked_bcast(myComm, orbitalSet->smallBox);
  }
  else
  {
    app_log() << "Perform FFT+spline and dump bspline tables to " << splinefile << std::endl;
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
      Vector<std::complex<double> > cG;
      int ncg=0;
      int ti=0;
      //int ti=SortBands[iorb].TwistIndex;
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti    //SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
        HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
        h_cG.read (H5FileID, path.str().c_str());
        ncg=cG.size();
      }
      myComm->bcast(ncg);
      if(ncg != Gvecs[ti].size())
      {
        APP_ABORT("Failed : ncg != Gvecs[ti].size()");
      }
      if(!root)
        cG.resize(ncg);
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
    if(root)
    {
      hdf_archive h5f;
      h5f.create(splinefile);
      einspline_engine<typename SPE::SplineType> bigtable(orbitalSet->MultiSpline);
      einspline_engine<typename SPE::SplineType> smalltable(orbitalSet->smallBox);
      std::string aname("EinsplineOpenAdoptor");
      h5f.write(aname,"bspline_type");
      h5f.write(lower,"lower_bound");
      h5f.write(upper,"upper_bound");
      h5f.write(bigtable,"spline_0");
      h5f.write(smalltable,"spline_1");
    }
  }
  app_log() << "TIME READBANDS " << now.elapsed() << std::endl;
}
}
#endif

