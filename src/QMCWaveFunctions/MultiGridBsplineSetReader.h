//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  typedef typename SA::SingleSplineType SingleSplineType;

  static const int EXTENDED=1;
  static const int LOCALIZED=0;

  /** if true, need to spline both real and imaginary
   */
  bool use_complex;
  /** actual SPOSet to be created and returned */
  ThisSPOSetType* thisSPOSet;

  TinyVector<int,3> coarse_mesh;
  TinyVector<int,3> coarse_stride;
  /** raw data on the original grid */
  Array<double,3> dense_r, dense_i;
  Array<double,3> coarse_r, coarse_i;
  std::vector<SingleSplineType*> spline_r;
  std::vector<SingleSplineType*> spline_i;
  Array<std::complex<double>,3> FFTbox;
  fftw_plan FFTplan;
  std::vector<int> OrbGroups;

  MultiGridBsplineSetReader(EinsplineSetBuilder* e)
    : BsplineReaderBase(e),use_complex(false),thisSPOSet(0),FFTplan(NULL)
  { }

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
    *target = new multi_UBspline_3d_z;
    APP_ABORT("MultiGridBsplineSetReader::export_MultiSpline(multi_UBspline_3d_z* target) not ready");
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
    app_log() << "  AdoptorName = MixedGridBsplineSet "<< std::endl;

    use_complex=thisSPOSet->is_complex;
    if(use_complex)
      app_log() << "  Using complex einspline table" << std::endl;
    else
      app_log() << "  Using real einspline table" << std::endl;

    check_twists(thisSPOSet->Extended,bandgroup);

    //create the coarse grid
    Ugrid xyz_grid[3];
    typename SA::BCType xyz_bc[3];

    typename ThisSPOSetType::bspline_type* extended=thisSPOSet->Extended;

    //setting grid based on the FFT grid & meshfactor
    bool havePsig=set_grid(extended->HalfG, xyz_grid, xyz_bc); 

    if(!havePsig)
    {
      APP_ABORT("Need psi_g. Must be old HDF5. Regenerate it with pwscf/pw2qmcpack\n");
    }

    int vol_factor=1;
    for(int j=0; j<3; ++j)
    {
      vol_factor *=GridFactor;
      coarse_mesh[j]=MeshSize[j]/GridFactor;
      coarse_stride[j]=GridFactor;
    }

    for(int j=0; j<3; ++j) xyz_grid[j].num=coarse_mesh[j];
    extended->create_spline(xyz_grid,xyz_bc);

    //create temporary single-value spline objects: dense/coarse pairs
    {
      int norbs_n=1; //using serial version for now to do the block allocation

      TinyVector<double,3> start(0.0);
      TinyVector<double,3> end(1.0);

      SingleSplineType* dummy=0; //dummy argument to handle c++-c template calls
      spline_r.resize(norbs_n*2,0);
      spline_i.resize(norbs_n*2,0);

      for(int i=0; i<norbs_n; ++i)
      {
        spline_r[2*i+LOCALIZED]=einspline::create(dummy,start,end,MeshSize,extended->HalfG);
        spline_r[2*i+EXTENDED]=einspline::create(dummy,start,end,coarse_mesh,extended->HalfG);
      }
      if(use_complex)
      {
        for(int i=0; i<norbs_n; ++i)
        {
          spline_i[2*i+LOCALIZED]=einspline::create(dummy,start,end,MeshSize,extended->HalfG);
          spline_i[2*i+EXTENDED]=einspline::create(dummy,start,end,coarse_mesh,extended->HalfG);
        }
      }
    }

    app_log().flush();

    //time to do the dense part
    for(int j=0; j<3; ++j) xyz_grid[j].num=MeshSize[j];
    int tableindex=thisSPOSet->myTableIndex=mybuilder->myTableIndex;
    const ParticleSet& ions=mybuilder->TargetPtcl.DistTables[tableindex]->origin();
    int ncenters=ions.getTotalNum()/(static_cast<int>(round(ions.Lattice.Volume/ions.PrimitiveLattice.Volume)));
    size_t loc_data_size=0;
    {
      const double onethird=1.0/3.0;
      double rc=ions.PrimitiveLattice.SimulationCellRadius;
      double pixels=static_cast<double>(MeshSize[0])*static_cast<double>(MeshSize[1])*static_cast<double>(MeshSize[2]);
      double rc_pixel=rc/std::pow(pixels,onethird);
      double n2=Rcut*0.5/rc_pixel;

      typedef QMCTraits::PosType pos_type;
      pos_type delta(n2/static_cast<double>(MeshSize[0]),n2/static_cast<double>(MeshSize[1]),n2/static_cast<double>(MeshSize[2]));

      app_log() << "  Estimated SubDomain " << 2.0*delta << std::endl;
      TinyVector<double,3> lower(0.0);
      TinyVector<double,3> upper(1.0);

      //prepare to add sub domains
      thisSPOSet->resizeSubDomains(ions.getTotalNum(),ncenters);
      copy(ions.PCID.begin(),ions.PCID.end(),thisSPOSet->PCID.begin());

      std::vector<int> boxes(ncenters,-1); 
      char s[1024];
      for(int i=0; i<ions.getTotalNum(); ++i)
      {
        int pcid=ions.PCID[i];
        if(boxes[pcid]<0)
        {
          pos_type u=ions.PrimitiveLattice.toUnit(ions.R[i]);
          lower=u-delta;
          upper=u+delta;
          sprintf(s,"SubDomain %4d at %10.5e %10.5e %10.5e\n  Bounds = [%10.5e,%10.5e) x [%10.5e,%10.5e) x [%10.5e,%10.5e) \n",
              pcid,u[0],u[1],u[2],lower[0],upper[0],lower[1],upper[1],lower[2],upper[2]);
          app_log()<< s;
          if(einspline::outOfBound(lower) || einspline::outOfBound(upper) || !einspline::validRange(lower,upper))
          {
            APP_ABORT("  Choose right-handed cell and place the atoms so that the subdomain is within [0,1)^3 of the supercell ");
          }
          loc_data_size+=thisSPOSet->setSubDomain(pcid,spline_r[LOCALIZED],lower,upper);
          check_twists(thisSPOSet->Localized[pcid],bandgroup);
          boxes[pcid]=pcid;
        }
      }
      app_log().flush();
    }

    {
      size_t org_data_size=thisSPOSet->sizeOfExtended()*vol_factor;
      app_log() << "Bspline Memory Use in MB: Original " << (org_data_size>>20)
        << " Global " << (thisSPOSet->sizeOfExtended()>>20)
        << " Local  " << (loc_data_size>>20) 
        << "\n  Saving factor= " << static_cast<double>(org_data_size)/static_cast<double>(thisSPOSet->sizeOfExtended()+loc_data_size)
        << std::endl;
    }

    //APP_ABORT("DONE");

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

      Vector<std::complex<double> > cG(mybuilder->Gvecs[0].size());
      int N = mybuilder->NumDistinctOrbitals;
      const std::vector<BandInfo>& cur_bands=bandgroup.myBands;
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        int ti=cur_bands[iorb].TwistIndex;

        get_psi_g(ti,spin,cur_bands[iorb].BandIndex,cG);
        fft_spline(cG,ti,0); 

        //assign the extended orbital
        thisSPOSet->set_spline(spline_r[EXTENDED],spline_i[EXTENDED],ti,ival,-1);

        for(int ic=0; ic<ncenters; ++ic)
          thisSPOSet->set_spline(spline_r[LOCALIZED],spline_i[LOCALIZED],ti,ival,ic);
      }
    }

    clear();

    return thisSPOSet;
  }

};
}
#endif

