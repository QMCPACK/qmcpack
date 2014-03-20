/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////
/** @file MultiGridBsplineSetReader.h
 */

#ifndef QMCPLUSPLUS_EINSPLINE_SPLINE_MIXEDGRID_ADOPTOR_PARALLEL_READER_H
#define QMCPLUSPLUS_EINSPLINE_SPLINE_MIXEDGRID_ADOPTOR_PARALLEL_READER_H

#include <mpi/point2point.h>
#include <mpi/collectives.h>

namespace qmcplusplus
{

/** reader for MultiGridBsplineSEt
 */
template<typename SA>
struct MultiGridBsplineSetReader: public BsplineReaderBase
{
  typedef MultiGridBsplineSet<SA> ThisSPOSetType;

  static const int EXTENDED=1;
  static const int LOCALIZED=0;

  /** if true, need to spline both real and imaginary
   */
  bool use_complex;
  /** dilation factor dense/corase only use one for all the directions
   */
  int GridFactor;
  /** actual SPOSet to be created and returned */
  ThisSPOSetType* thisSPOSet;

  TinyVector<int,3> coarse_mesh;
  TinyVector<int,3> coarse_stride;
  /** raw data on the original grid */
  Array<double,3> dense_r, dense_i;
  Array<double,3> coarse_r, coarse_i;
  vector<UBspline_3d_d*> spline_r;
  vector<UBspline_3d_d*> spline_i;
  Array<complex<double>,3> FFTbox;
  fftw_plan FFTplan;
  vector<int> OrbGroups;

  MultiGridBsplineSetReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e),use_complex(false),GridFactor(2),thisSPOSet(0),FFTplan(NULL)
  {}

  ~MultiGridBsplineSetReader()
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
      free(spline_i[i]->coefs);
      free(spline_i[i]);
    }
    spline_r.clear();
    spline_i.clear();
    if(FFTplan!=NULL) fftw_destroy_plan(FFTplan);
    FFTplan=NULL;
  }

  /** fft and spline a band
   * @param ti twist index
   * @param iorb orbital index
   */
  void fft_spline(Vector<complex<double> >& cG, int ti, int iorb)
  {
    unpack4fftw(cG,mybuilder->Gvecs[0],mybuilder->MeshSize,FFTbox);
    fftw_execute (FFTplan);
    if(use_complex)
      fix_phase_rotate_c2c(FFTbox,dense_r, dense_i,mybuilder->TwistAngles[ti]);
    else
      fix_phase_rotate_c2r(FFTbox,dense_r,mybuilder->TwistAngles[ti]);

    for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=coarse_stride[0])
      for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=coarse_stride[1])
        for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=coarse_stride[2])
          coarse_r(i,j,k)=dense_r(i2,j2,k2);

    einspline::set(spline_r[2*iorb+LOCALIZED],dense_r.data());
    einspline::set(spline_r[2*iorb+EXTENDED],coarse_r.data());

    if(use_complex)
    {
      for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=coarse_stride[0])
        for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=coarse_stride[1])
          for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=coarse_stride[2])
            coarse_i(i,j,k)=dense_i(i2,j2,k2);
      einspline::set(spline_i[2*iorb+LOCALIZED],dense_i.data());
      einspline::set(spline_i[2*iorb+EXTENDED],coarse_i.data());
    }
  }

  SPOSetBase* create_spline_set(int spin, const BandInfoGroup& bandgroup)
  {
    ReportEngine PRE("MultiGridBsplineSetReader","create_spline_set(int, EinsplineSet*)");

    thisSPOSet=new ThisSPOSetType;
    app_log() << "  AdoptorName = MixedGridBsplineSet "<< endl;

    use_complex=thisSPOSet->is_complex;
    if(use_complex)
      app_log() << "  Using complex einspline table" << endl;
    else
      app_log() << "  Using real einspline table" << endl;

    check_twists(thisSPOSet->Extended,bandgroup);
    check_twists(thisSPOSet->Localized,bandgroup);

    typename ThisSPOSetType::bspline_type* extended=thisSPOSet->Extended;
    typename ThisSPOSetType::bspline_type* localized=thisSPOSet->Localized;

    Ugrid xyz_grid[3];
    typename SA::BCType xyz_bc[3];

    bool havePsig=set_grid(localized->HalfG,xyz_grid, xyz_bc);

    if(!havePsig)
    {
      APP_ABORT("Need psi_g. Must be old HDF5. Regenerate it with pwscf/pw2qmcpack\n");
    }
    for(int j=0; j<3; ++j)
    {
      coarse_mesh[j]=MeshSize[j]/GridFactor;
      coarse_stride[j]=GridFactor;
    }

    localized->create_spline(xyz_grid,xyz_bc);

    for(int j=0; j<3; ++j) xyz_grid[j].num=coarse_mesh[j];
    extended->create_spline(xyz_grid,xyz_bc);

    {//create temporary single-value spline objects
      TinyVector<double,3> start(0.0);
      TinyVector<double,3> end(1.0);
      int norbs_n=1; //using serial version for now
      UBspline_3d_d* dummy=0; //dummy argument to handle c++-c template calls
      spline_r.resize(norbs_n*2,0);
      for(int i=0; i<spline_r.size()/2; ++i)
      {
        spline_r[2*i+LOCALIZED]=einspline::create(dummy,start,end,MeshSize,extended->HalfG);
        spline_r[2*i+EXTENDED]=einspline::create(dummy,start,end,coarse_mesh,extended->HalfG);
      }
      if(use_complex)
      {
        spline_i.resize(norbs_n*2,0);
        for(int i=0; i<spline_i.size()/2; ++i)
        {
          spline_i[2*i+LOCALIZED]=einspline::create(dummy,start,end,MeshSize,extended->HalfG);
          spline_i[2*i+EXTENDED]=einspline::create(dummy,start,end,coarse_mesh,extended->HalfG);
        }
      }
    }

    app_log() << "  Original Mesh " << MeshSize << endl;
    app_log() << "  Coarse Mesh " << coarse_mesh << endl;
    app_log().flush();


    /** every MPI node does this: to be parallelized */
    {
      FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
      FFTplan = fftw_plan_dft_3d
                (MeshSize[0], MeshSize[1], MeshSize[2],
                 reinterpret_cast<fftw_complex*>(FFTbox.data()),
                 reinterpret_cast<fftw_complex*>(FFTbox.data()),
                 +1, FFTW_ESTIMATE);
      dense_r.resize(MeshSize[0],MeshSize[1],MeshSize[2]);
      coarse_r.resize(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      if(use_complex)
      {
        dense_i.resize(MeshSize[0],MeshSize[1],MeshSize[2]);
        coarse_i.resize(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      }

      Vector<complex<double> > cG(mybuilder->Gvecs[0].size());
      int N = mybuilder->NumDistinctOrbitals;
      const vector<BandInfo>& cur_bands=bandgroup.myBands;
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        int ti=cur_bands[iorb].TwistIndex;
        cout << "Getting PW coefficients for " << iorb << " at twist=" << ti << endl;
        get_psi_g(ti,spin,cur_bands[iorb].BandIndex,cG);
        fft_spline(cG,ti,0); 
        thisSPOSet->set_spline(spline_r[EXTENDED],spline_i[EXTENDED],ti,ival,EXTENDED);//localized/dense
        thisSPOSet->set_spline(spline_r[LOCALIZED],spline_i[LOCALIZED],ti,ival,LOCALIZED);//extended/coarse
      }

    }

    clear();

    return thisSPOSet;
  }

};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: EinsplineSetBuilderESHDF.cpp 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
