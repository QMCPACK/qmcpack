//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_EINSPLINE_SPLINEOPENADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_SPLINEOPENADOPTOR_H
#include "Utilities/ProgressReportEngine.h"
#include <QMCWaveFunctions/einspline_helper.hpp>
#include <spline/einspline_util.hpp>
#include <fftw3.h>

namespace qmcplusplus
{

/** adoptor class for the non-periodic systems
 *
 * Support two-level grid system for big systems
 * - MultiSpline : full cell with the coarse grid
 * - smallBox    : small cell with the original grid
 * No coordinate transform is needed for the open BC with a orthorhombic cell
 */
template<typename ST, typename TT, unsigned D>
struct SplineOpenAdoptor
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
  SplineType *smallBox;

  TinyVector<ST,D> Lower;
  TinyVector<ST,D> Upper;
  TinyVector<ST,D> L;
  TinyVector<ST,D> InvL;

  Tensor<ST,D> GGt;
  UnitCellType PrimLattice;
  UnitCellType SuperLattice;

  ///used for testing only
  bool UseFullGrid;
  ///grid transform
  GridConvert<ST> gTransform;

  // Temporary storage for Eispline calls
  StorageValueVector_t myV, myL;
  StorageGradVector_t      myG;
  StorageGradHessVector_t  myGH;

  SplineOpenAdoptor(): MultiSpline(0), smallBox(0)
  {
  }

  void resizeStorage(int n, int nv)
  {
    for(int i=0; i<D; ++i)
      L[i]=SuperLattice.R(i,i);
    for(int i=0; i<D; ++i)
      InvL[i]=1.0/L[i];
    myV.resize(n);
    myL.resize(n);
    myG.resize(n);
  }

  /** create MultiSpline for the full cell with a coarse grid
   * @param mesh original mesh (e.g. FFT)
   * @param n number of states
   * @param samegrid if true, use full grid
   */
  void create_spline(TinyVector<int,D>& mesh, int n, bool samegrid)
  {
    UseFullGrid=samegrid;
    Ugrid xyz_grid[D];
    BCType xyz_bc[D];
    for(int i=0; i<D; ++i)
    {
      L[i]=PrimLattice.R(i,i);
      InvL[i]=1.0/L[i];
      xyz_grid[i].start = 0.0;
      xyz_grid[i].end = L[i]; //1.0;
      xyz_grid[i].num = (samegrid)?mesh[i]:mesh[i]/2;
      xyz_bc[i].lCode=xyz_bc[i].rCode=PERIODIC;
      gTransform.BaseOffset[i]=0;
      gTransform.BaseN[i]=xyz_grid[i].num+3;
    }
    SplineType* dummy=0;
    MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);
  }

  /** create smallBox
   * @param dense single-spline for the full grid
   * @param lower Cartesian lower bound
   * @param upper Cartesian upper bound
   */
  template <typename UBspline, typename PT>
  void add_box(UBspline* dense, PT& lower, PT& upper)
  {
    if(smallBox==0)
      gTransform.create(smallBox,dense,lower,upper,MultiSpline->num_splines);
    Lower=lower;//+2.0*gTransform.Delta;
    Upper=upper;//-2.0*gTransform.Delta;
  }

  /** set ival-th state
   * @param dense a single-spline on the full grid
   * @param coarse a single-spline on the half grid
   * @param ival the state to be initialized
   */
  template <typename UBspline>
  void init_spline(UBspline* dense, UBspline* coarse, int ival)
  {
    if(UseFullGrid)
      einspline::set(MultiSpline,ival,dense,  gTransform.BaseOffset,gTransform.BaseN);
    else
      einspline::set(MultiSpline,ival,coarse, gTransform.BaseOffset,gTransform.BaseN);
    einspline::set(smallBox,ival,dense,gTransform.Offset,gTransform.N);
  }

  inline bool isready()
  {
    return true;
  }

  ///not needed for molecular system
  inline void convertPos(const PointType& r, PointType& ru)
  {
    for(int i=0; i<D; ++i)
      ru[i] = r[i]-L[i]*std::floor(r[i]*InvL[i]);
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    TinyVector<ST,D> ru;
    convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate(smallBox,ru,myV);
    else
      einspline::evaluate(MultiSpline,ru,myV);
    for (int j=0; j<psi.size(); j++)
      psi[j]=static_cast<TT>(myV[j]);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    TinyVector<ST,D> ru;
    convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate_vgl(smallBox,ru,myV,myG,myL);
    else
      einspline::evaluate_vgl(MultiSpline,ru,myV,myG,myL);
    const int N=psi.size();
    for(int j=0; j<N; ++j)
      psi[j]=myV[j];
    for(int j=0; j<N; ++j)
      dpsi[j]=myG[j];
    for(int j=0; j<N; ++j)
      d2psi[j]=myL[j];
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {}
};

template<typename T>
struct SplineOpenAdoptorReader: BsplineReaderBase
{
  typedef SplineOpenAdoptor<T,double,3> adoptor_type;

  SplineOpenAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
  {}

  SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
  {
    ReportEngine PRE("SplineOpenAdoptorReader","create_spline_set(int, EinsplineSet*)");
    typedef T spline_data_ype;
    BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
    init(orbitalSet,bspline);
    int N=bspline->getOrbitalSetSize();
    bspline->resizeStorage(N,N);
    std::string H5FileName(mybuilder->H5FileName);
    H5OrbSet set(H5FileName, spin, N);
    bool havePsir=!(mybuilder->ReadGvectors_ESHDF());
    bool fullgrid=false;
    TinyVector<int,3> MeshSize=mybuilder->MeshSize;
    TinyVector<int,3> coarse_mesh(MeshSize[0]/2,MeshSize[1]/2,MeshSize[2]/2);
    TinyVector<double,3> start(0.0);
    TinyVector<double,3> end(1.0);
    for(int i=0; i<3; ++i)
      end(i)=mybuilder->SuperLattice(i,i);
    //create dense and coarse UBspline_3d_d
    UBspline_3d_d* dense=0;
    dense=einspline::create(dense,start,end,MeshSize,PERIODIC);
    UBspline_3d_d* coarse=0;
    coarse=einspline::create(coarse,start,end,coarse_mesh,PERIODIC);
    //determine the bonding box
    TinyVector<double,3> lower=end;
    TinyVector<double,3> upper=start;
    const Vector<TinyVector<double,OHMMS_DIM> >& IonPos(mybuilder->IonPos);;
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
    bspline->create_spline(MeshSize,N,fullgrid);
    app_log() << "Original Mesh " << MeshSize << std::endl;
    app_log() << "Coarse Mesh " << coarse_mesh << std::endl;
    app_log() << "  Adding a small box" << "\n  LowerBound " << lower << "\n  UpperBound " << upper << std::endl;
    app_log().flush();
    bspline->add_box(dense,lower,upper);
    int foundspline=0;
    std::string splinefile=make_spline_filename(H5FileName,spin,0,MeshSize);
    Timer now;
    if(myComm->rank()==0)
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
          einspline_engine<typename adoptor_type::SplineType> bigtable(bspline->MultiSpline);
          einspline_engine<typename adoptor_type::SplineType> smalltable(bspline->smallBox);
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
      Array<double,3> bigD(MeshSize[0],MeshSize[1],MeshSize[2]);
      Array<double,3> smallD(coarse_mesh[0],coarse_mesh[1],coarse_mesh[2]);
      TinyVector<double,3> TwistAngle(0.0,0.0,0.0);
      int ti=0;
      int ncg=mybuilder->Gvecs[ti].size();
      Vector<std::complex<double> > cG(ncg);
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        get_psi_g(ti,spin,mybuilder->SortBands[iorb].BandIndex,cG);
        unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
        fftw_execute (FFTplan);
        fix_phase_rotate_c2r(FFTbox,bigD,TwistAngle);
        for(int i=0,i2=0; i<coarse_mesh[0]; ++i,i2+=2)
          for(int j=0,j2=0; j<coarse_mesh[1]; ++j,j2+=2)
            for(int k=0,k2=0; k<coarse_mesh[2]; ++k,k2+=2)
              smallD(i,j,k)=bigD(i2,j2,k2);
        einspline::set(dense,bigD.data());
        einspline::set(coarse,smallD.data());
        bspline->init_spline(dense,coarse,ival);
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
        std::string aname("EinsplineOpenAdoptor");
        h5f.write(aname,"bspline_type");
        h5f.write(lower,"lower_bound");
        h5f.write(upper,"upper_bound");
        h5f.write(bigtable,"spline_0");
        h5f.write(smalltable,"spline_1");
      }
    }
    app_log() << "TIME READBANDS " << now.elapsed() << std::endl;
    return bspline;
  }
};
}
#endif

