//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@intel.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file SplineC2RSoA.h
 *
 * Adoptor classes to handle complex-to-(real,complex) with arbitrary precision
 */
#ifndef QMCPLUSPLUS_EINSPLINE_C2R_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_C2R_ADOPTOR_H

#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorBase.h"
namespace qmcplusplus
{

/** adoptor class to match std::complex<ST> spline with TT real SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 *
 * Requires temporage storage and multiplication of phase vectors
 */
template<typename ST, typename TT>
struct SplineC2RSoA: public SplineAdoptorBase<ST,3>
{
  static const int D=3;
  using BaseType=SplineAdoptorBase<ST,3>;
  using SplineType=typename bspline_traits<ST,3>::SplineType;
  using BCType=typename bspline_traits<ST,3>::BCType;
  using DataType=ST;
  using PointType=typename BaseType::PointType;
  using SingleSplineType=typename BaseType::SingleSplineType;

  using vContainer_type=Vector<ST,aligned_allocator<ST> >;
  using gContainer_type=VectorSoaContainer<ST,3>;
  using hContainer_type=VectorSoaContainer<ST,6>;

  using BaseType::first_spo;
  using BaseType::last_spo;
  using BaseType::GGt;
  using BaseType::PrimLattice;
  using BaseType::kPoints;
  using BaseType::MakeTwoCopies;
  using BaseType::offset_cplx;
  using BaseType::offset_real;

  ///number of complex bands
  int nComplexBands;
  ///number of points of the original grid
  int BaseN[3];
  ///offset of the original grid, always 0
  int BaseOffset[3];
  ///multi bspline set
  MultiBspline<ST>* SplineInst;
  ///expose the pointer to reuse the reader and only assigned with create_spline
  SplineType* MultiSpline;

  vContainer_type  mKK;
  VectorSoaContainer<ST,3>  myKcart;

  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;

  SplineC2RSoA(): BaseType(),SplineInst(nullptr), MultiSpline(nullptr)
  {
    this->is_complex=true;
    this->is_soa_ready=true;
    this->AdoptorName="SplineC2RSoAAdoptor";
    this->KeyWord="SplineC2RSoA";
  }

  ///** copy the base property */
  //SplineC2RSoA(BaseType& rhs): BaseType(rhs)
  //{
  //  this->is_complex=true;
  //  this->AdoptorName="SplineC2RSoA";
  //  this->KeyWord="C2RSoA";
  //}

  SplineC2RSoA(const SplineC2RSoA& a):
    SplineAdoptorBase<ST,3>(a),SplineInst(a.SplineInst),MultiSpline(nullptr),
    nComplexBands(a.nComplexBands),mKK(a.mKK), myKcart(a.myKcart)
  {
    const size_t n=a.myL.size();
    myV.resize(n); myG.resize(n); myL.resize(n); myH.resize(n);
  }

  ~SplineC2RSoA()
  {
    if(MultiSpline != nullptr) delete SplineInst;
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseType::init_base(n);
    size_t npad=getAlignedSize<ST>(2*n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);
  }

  void bcast_tables(Communicate* comm)
  {
    chunked_bcast(comm, MultiSpline);
  }

  void gather_tables(Communicate* comm)
  {
    if(comm->size()==1) return;
    const int Nbands = kPoints.size();
    const int Nbandgroups = comm->size();
    std::vector<int> offset(Nbandgroups+1,0);
    FairDivideLow(Nbands,Nbandgroups,offset);

#ifdef QMC_CUDA
    for(size_t ib=0; ib<offset.size(); ib++)
      offset[ib] = offset[ib]*2;
    gatherv(comm, MultiSpline, MultiSpline->z_stride, offset);
#else
    // complex bands
    int gid=1;
    offset_cplx.resize(Nbandgroups+1,0);
    for(int ib=0; ib<Nbands; ++ib)
    {
      if(ib==offset[gid]) gid++;
      if(MakeTwoCopies[ib])
        offset_cplx[gid]++;
    }
    for(int bg=0; bg<Nbandgroups; ++bg)
      offset_cplx[bg+1] = offset_cplx[bg+1]*2+offset_cplx[bg];
    gatherv(comm, MultiSpline, MultiSpline->z_stride, offset_cplx);

    // real bands
    gid=1;
    offset_real.resize(Nbandgroups+1,0);
    for(int ib=0; ib<Nbands; ++ib)
    {
      if(ib==offset[gid]) gid++;
      if(!MakeTwoCopies[ib])
        offset_real[gid]++;
    }
    offset_real[0]=nComplexBands*2;
    for(int bg=0; bg<Nbandgroups; ++bg)
      offset_real[bg+1] = offset_real[bg+1]*2+offset_real[bg];
    gatherv(comm, MultiSpline, MultiSpline->z_stride, offset_real);
#endif
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    resize_kpoints();
    SplineInst=new MultiBspline<ST>();
    SplineInst->create(xyz_g,xyz_bc,myV.size());
    MultiSpline=SplineInst->spline_m;

    for(size_t i=0; i<D; ++i)
    {
      BaseOffset[i]=0;
      BaseN[i]=xyz_g[i].num+3;
    }
    qmc_common.memory_allocated += SplineInst->sizeInByte();
  }

  inline void flush_zero()
  {
    SplineInst->flush_zero();
  }

  /** remap kPoints to pack the double copy */
  inline void resize_kpoints()
  {
    nComplexBands=this->remap_kpoints();
    int nk=kPoints.size();
    mKK.resize(nk);
    myKcart.resize(nk);
    for(size_t i=0; i<nk; ++i)
    {
      mKK[i]=-dot(kPoints[i],kPoints[i]);
      myKcart(i)=kPoints[i];
    }
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
#ifdef QMC_CUDA
    // GPU code needs the old ordering.
    int iband=ispline;
#else
    int iband=this->BandIndexMap[ispline];
#endif
    SplineInst->copy_spline(spline_r,2*iband  ,BaseOffset, BaseN);
    SplineInst->copy_spline(spline_i,2*iband+1,BaseOffset, BaseN);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    Vector<ST> v_r(psi_r,0), v_i(psi_i,0);
#ifdef QMC_CUDA
    // GPU code needs the old ordering.
    int iband=ispline;
#else
    int iband=this->BandIndexMap[ispline];
#endif
    SplineInst->set(2*iband  ,v_r);
    SplineInst->set(2*iband+1,v_i);
  }


  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i,
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
  }

  bool read_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o<<"spline_" << SplineAdoptorBase<ST,D>::MyIndex;
    einspline_engine<SplineType> bigtable(SplineInst->spline_m);
    return h5f.read(bigtable,o.str().c_str());//"spline_0");
  }

  bool write_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o<<"spline_" << SplineAdoptorBase<ST,D>::MyIndex;
    einspline_engine<SplineType> bigtable(SplineInst->spline_m);
    return h5f.write(bigtable,o.str().c_str());//"spline_0");
  }

  TT evaluate_dot(const ParticleSet& P, int iat, const TT* restrict arow, ST* scratch, bool compute_spline=true)
  {
    Vector<ST> vtmp(scratch,myV.size());
    const PointType& r=P.activeR(iat);
    if(compute_spline)
    {
      PointType ru(PrimLattice.toUnit_floor(r));
      SplineInst->evaluate(ru,vtmp);
    }

    TT res=TT();
    const size_t N=kPoints.size();
    const ST x=r[0], y=r[1], z=r[2];
    const ST* restrict kx=myKcart.data(0);
    const ST* restrict ky=myKcart.data(1);
    const ST* restrict kz=myKcart.data(2);

    const TT* restrict arow_s=arow+first_spo;
    #pragma omp simd reduction(+:res)
    for (size_t j=0; j<nComplexBands; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;
      const ST val_r=vtmp[jr];
      const ST val_i=vtmp[ji];
      ST s, c;
      sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
      res+=arow_s[jr] * (val_r*c-val_i*s);
      res+=arow_s[ji] * (val_i*c+val_r*s);
    }
    const TT* restrict arow_c=arow+first_spo+nComplexBands;
    #pragma omp simd reduction(+:res)
    for (size_t j=nComplexBands; j<N; j++)
    {
      const ST val_r=vtmp[2*j  ];
      const ST val_i=vtmp[2*j+1];
      ST s, c;
      sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
      res+=arow_c[j]*(val_r*c-val_i*s);
    }
    return res;
  }

  template<typename VV>
  inline void assign_v(const PointType& r, const vContainer_type& myV, VV& psi)
  {
    const size_t N=kPoints.size();
    const ST x=r[0], y=r[1], z=r[2];
    const ST* restrict kx=myKcart.data(0);
    const ST* restrict ky=myKcart.data(1);
    const ST* restrict kz=myKcart.data(2);
#if defined(USE_VECTOR_ML)
    {//reuse myH
      ST* restrict KdotR=myH.data(0);
      ST* restrict CosV=myH.data(1);
      ST* restrict SinV=myH.data(2);
      #pragma omp simd
      for(size_t j=0; j<N; ++j)
        KdotR[j]=-(x*kx[j]+y*ky[j]+z*kz[j]);

      eval_e2iphi(N,KdotR,CosV,SinV);

      #pragma omp simd
      for (size_t j=0,psiIndex=first_spo; j<nComplexBands; j++, psiIndex+=2)
      {
        const ST val_r=myV[2*j  ];
        const ST val_i=myV[2*j+1];
        psi[psiIndex  ] = val_r*CosV[j]-val_i*SinV[j];
        psi[psiIndex+1] = val_i*CosV[j]+val_r*SinV[j];
      }
      #pragma omp simd
      for (size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands; j<N; j++,psiIndex++)
      {
        const ST val_r=myV[2*j  ];
        const ST val_i=myV[2*j+1];
        psi[psiIndex  ] = val_r*CosV[j]-val_i*SinV[j];
      }
    }
#else
    {
      TT* restrict psi_s=psi.data()+first_spo;
      #pragma omp simd
      for (size_t j=0; j<nComplexBands; j++)
      {
        ST s, c;
        const size_t jr=j<<1;
        const size_t ji=jr+1;
        const ST val_r=myV[jr];
        const ST val_i=myV[ji];
        sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
        psi_s[jr] = val_r*c-val_i*s;
        psi_s[ji] = val_i*c+val_r*s;
      }
    }

    {
      TT* restrict psi_s=psi.data()+first_spo+nComplexBands;
      #pragma omp simd
      for (size_t j=nComplexBands; j<N; j++)
      {
        ST s, c;
        const ST val_r=myV[2*j  ];
        const ST val_i=myV[2*j+1];
        sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
        psi_s[j] = val_r*c-val_i*s;
      }
    }
#endif
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate(ru,myV);
    assign_v(r,myV,psi);
  }

  template<typename VM>
  inline void evaluateValues(const VirtualParticleSet& VP, VM& psiM)
  {
    const size_t m=psiM.cols();
    for(int iat=0; iat<VP.getTotalNum(); ++iat)
    {
      Vector<TT> psi(psiM[iat],m);
      evaluate_v(VP,iat,psi);
    }
  }

  inline size_t estimateMemory(const int nP) { return 0; }

  /** assign_vgl
   */
  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    CONSTEXPR ST zero(0);
    CONSTEXPR ST two(2);
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);
    const ST x=r[0], y=r[1], z=r[2];
    const ST symGG[6]={GGt[0],GGt[1]+GGt[3],GGt[2]+GGt[6],GGt[4],GGt[5]+GGt[7],GGt[8]};

    const ST* restrict k0=myKcart.data(0); ASSUME_ALIGNED(k0);
    const ST* restrict k1=myKcart.data(1); ASSUME_ALIGNED(k1);
    const ST* restrict k2=myKcart.data(2); ASSUME_ALIGNED(k2);

    const ST* restrict g0=myG.data(0); ASSUME_ALIGNED(g0);
    const ST* restrict g1=myG.data(1); ASSUME_ALIGNED(g1);
    const ST* restrict g2=myG.data(2); ASSUME_ALIGNED(g2);
    const ST* restrict h00=myH.data(0); ASSUME_ALIGNED(h00);
    const ST* restrict h01=myH.data(1); ASSUME_ALIGNED(h01);
    const ST* restrict h02=myH.data(2); ASSUME_ALIGNED(h02);
    const ST* restrict h11=myH.data(3); ASSUME_ALIGNED(h11);
    const ST* restrict h12=myH.data(4); ASSUME_ALIGNED(h12);
    const ST* restrict h22=myH.data(5); ASSUME_ALIGNED(h22);

    const size_t N=kPoints.size();
    const size_t nsplines=myL.size();
#if defined(PRECOMPUTE_L)
    for(size_t j=0; j<nsplines; ++j)
    {
      myL[j]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
    }
#endif

    #pragma omp simd
    for (size_t j=0; j<nComplexBands; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[jr]+g01*g1[jr]+g02*g2[jr];
      const ST dY_r = g10*g0[jr]+g11*g1[jr]+g12*g2[jr];
      const ST dZ_r = g20*g0[jr]+g21*g1[jr]+g22*g2[jr];

      const ST dX_i = g00*g0[ji]+g01*g1[ji]+g02*g2[ji];
      const ST dY_i = g10*g0[ji]+g11*g1[ji]+g12*g2[ji];
      const ST dZ_i = g20*g0[ji]+g21*g1[ji]+g22*g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;

#if defined(PRECOMPUTE_L)
      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#else
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],symGG);
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],symGG);
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif


      const size_t psiIndex=first_spo+jr;
      //this will be fixed later
      psi[psiIndex  ]=c*val_r-s*val_i;
      psi[psiIndex+1]=c*val_i+s*val_r;
      d2psi[psiIndex  ]=c*lap_r-s*lap_i;
      d2psi[psiIndex+1]=c*lap_i+s*lap_r;
      //this will go way with Determinant
      dpsi[psiIndex  ][0]=c*gX_r-s*gX_i;
      dpsi[psiIndex  ][1]=c*gY_r-s*gY_i;
      dpsi[psiIndex  ][2]=c*gZ_r-s*gZ_i;
      dpsi[psiIndex+1][0]=c*gX_i+s*gX_r;
      dpsi[psiIndex+1][1]=c*gY_i+s*gY_r;
      dpsi[psiIndex+1][2]=c*gZ_i+s*gZ_r;
    }

    #pragma omp simd
    for (size_t j=nComplexBands; j<N; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[jr]+g01*g1[jr]+g02*g2[jr];
      const ST dY_r = g10*g0[jr]+g11*g1[jr]+g12*g2[jr];
      const ST dZ_r = g20*g0[jr]+g21*g1[jr]+g22*g2[jr];

      const ST dX_i = g00*g0[ji]+g01*g1[ji]+g02*g2[ji];
      const ST dY_i = g10*g0[ji]+g11*g1[ji]+g12*g2[ji];
      const ST dZ_i = g20*g0[ji]+g21*g1[ji]+g22*g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;

      const size_t psiIndex=first_spo+nComplexBands+j;
      psi[psiIndex  ]=c*val_r-s*val_i;
      //this will be fixed later
      dpsi[psiIndex  ][0]=c*gX_r-s*gX_i;
      dpsi[psiIndex  ][1]=c*gY_r-s*gY_i;
      dpsi[psiIndex  ][2]=c*gZ_r-s*gZ_i;

#if defined(PRECOMPUTE_L)
      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#else
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],symGG);
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],symGG);
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif
      d2psi[psiIndex  ]=c*lap_r-s*lap_i;
    }
  }

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  template<typename VV, typename GV>
  inline void assign_vgl_from_l(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    CONSTEXPR ST two(2);
    const ST x=r[0], y=r[1], z=r[2];

    const ST* restrict k0=myKcart.data(0); ASSUME_ALIGNED(k0);
    const ST* restrict k1=myKcart.data(1); ASSUME_ALIGNED(k1);
    const ST* restrict k2=myKcart.data(2); ASSUME_ALIGNED(k2);

    const ST* restrict g0=myG.data(0); ASSUME_ALIGNED(g0);
    const ST* restrict g1=myG.data(1); ASSUME_ALIGNED(g1);
    const ST* restrict g2=myG.data(2); ASSUME_ALIGNED(g2);

    const size_t N=kPoints.size();
    const size_t nsplines=myL.size();

    #pragma omp simd
    for (size_t j=0; j<nComplexBands; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g0[jr];
      const ST dY_r = g1[jr];
      const ST dZ_r = g2[jr];

      const ST dX_i = g0[ji];
      const ST dY_i = g1[ji];
      const ST dZ_i = g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;

      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);

      //this will be fixed later
      const size_t psiIndex=first_spo+jr;
      psi[psiIndex  ]=c*val_r-s*val_i;
      psi[psiIndex+1]=c*val_i+s*val_r;
      d2psi[psiIndex  ]=c*lap_r-s*lap_i;
      d2psi[psiIndex+1]=c*lap_i+s*lap_r;
      //this will go way with Determinant
      dpsi[psiIndex  ][0]=c*gX_r-s*gX_i;
      dpsi[psiIndex  ][1]=c*gY_r-s*gY_i;
      dpsi[psiIndex  ][2]=c*gZ_r-s*gZ_i;
      dpsi[psiIndex+1][0]=c*gX_i+s*gX_r;
      dpsi[psiIndex+1][1]=c*gY_i+s*gY_r;
      dpsi[psiIndex+1][2]=c*gZ_i+s*gZ_r;
    }

    const size_t nComputed=2*nComplexBands;
    #pragma omp simd
    for (size_t j=nComplexBands; j<N; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g0[jr];
      const ST dY_r = g1[jr];
      const ST dZ_r = g2[jr];

      const ST dX_i = g0[ji];
      const ST dY_i = g1[ji];
      const ST dZ_i = g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;
      const size_t psiIndex=first_spo+nComplexBands+j;
      psi[psiIndex  ]=c*val_r-s*val_i;
      //this will be fixed later
      dpsi[psiIndex  ][0]=c*gX_r-s*gX_i;
      dpsi[psiIndex  ][1]=c*gY_r-s*gY_i;
      dpsi[psiIndex  ][2]=c*gZ_r-s*gZ_i;

      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
      d2psi[psiIndex  ]=c*lap_r-s*lap_i;
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl(r,psi,dpsi,d2psi);
  }

  /** identical to assign_vgl but the output container is SoA container
   */
  template<typename VGL>
  inline void assign_vgl_soa(const PointType& r, VGL& vgl)
  {
    const int N=kPoints.size();
    CONSTEXPR ST zero(0);
    CONSTEXPR ST two(2);
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);
    const ST x=r[0], y=r[1], z=r[2];
    const ST symGG[6]={GGt[0],GGt[1]+GGt[3],GGt[2]+GGt[6],GGt[4],GGt[5]+GGt[7],GGt[8]};

    const ST* restrict k0=myKcart.data(0); ASSUME_ALIGNED(k0);
    const ST* restrict k1=myKcart.data(1); ASSUME_ALIGNED(k1);
    const ST* restrict k2=myKcart.data(2); ASSUME_ALIGNED(k2);

    const ST* restrict g0=myG.data(0); ASSUME_ALIGNED(g0);
    const ST* restrict g1=myG.data(1); ASSUME_ALIGNED(g1);
    const ST* restrict g2=myG.data(2); ASSUME_ALIGNED(g2);
    const ST* restrict h00=myH.data(0); ASSUME_ALIGNED(h00);
    const ST* restrict h01=myH.data(1); ASSUME_ALIGNED(h01);
    const ST* restrict h02=myH.data(2); ASSUME_ALIGNED(h02);
    const ST* restrict h11=myH.data(3); ASSUME_ALIGNED(h11);
    const ST* restrict h12=myH.data(4); ASSUME_ALIGNED(h12);
    const ST* restrict h22=myH.data(5); ASSUME_ALIGNED(h22);
    const size_t nsplines=myL.size();
#if defined(PRECOMPUTE_L)
    #pragma omp simd
    for(size_t j=0; j<nsplines; ++j)
    {
      myL[j]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
    }
#endif

    TT* restrict psi =vgl.data(0)+first_spo; ASSUME_ALIGNED(psi);
    TT* restrict vl_x=vgl.data(1)+first_spo; ASSUME_ALIGNED(vl_x);
    TT* restrict vl_y=vgl.data(2)+first_spo; ASSUME_ALIGNED(vl_y);
    TT* restrict vl_z=vgl.data(3)+first_spo; ASSUME_ALIGNED(vl_z);
    TT* restrict vl_l=vgl.data(4)+first_spo; ASSUME_ALIGNED(vl_l);

    #pragma omp simd
    for (size_t j=0; j<nComplexBands; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[jr]+g01*g1[jr]+g02*g2[jr];
      const ST dY_r = g10*g0[jr]+g11*g1[jr]+g12*g2[jr];
      const ST dZ_r = g20*g0[jr]+g21*g1[jr]+g22*g2[jr];

      const ST dX_i = g00*g0[ji]+g01*g1[ji]+g02*g2[ji];
      const ST dY_i = g10*g0[ji]+g11*g1[ji]+g12*g2[ji];
      const ST dZ_i = g20*g0[ji]+g21*g1[ji]+g22*g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;

#if defined(PRECOMPUTE_L)
      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#else
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],symGG);
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],symGG);
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif

      //this will be fixed later
      psi[jr]=c*val_r-s*val_i;
      psi[ji]=c*val_i+s*val_r;
      vl_x[jr]=c*gX_r-s*gX_i;
      vl_x[ji]=c*gX_i+s*gX_r;
      vl_y[jr]=c*gY_r-s*gY_i;
      vl_y[ji]=c*gY_i+s*gY_r;
      vl_z[jr]=c*gZ_r-s*gZ_i;
      vl_z[ji]=c*gZ_i+s*gZ_r;
      vl_l[jr]=c*lap_r-s*lap_i;
      vl_l[ji]=c*lap_i+s*lap_r;
    }

    #pragma omp simd
    for (size_t j=nComplexBands; j<N; j++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[jr]+g01*g1[jr]+g02*g2[jr];
      const ST dY_r = g10*g0[jr]+g11*g1[jr]+g12*g2[jr];
      const ST dZ_r = g20*g0[jr]+g21*g1[jr]+g22*g2[jr];

      const ST dX_i = g00*g0[ji]+g01*g1[ji]+g02*g2[ji];
      const ST dY_i = g10*g0[ji]+g11*g1[ji]+g12*g2[ji];
      const ST dZ_i = g20*g0[ji]+g21*g1[ji]+g22*g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;
#if defined(PRECOMPUTE_L)
      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#else
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],symGG);
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],symGG);
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif
      const size_t psiIndex=first_spo+nComplexBands+j;
      //const size_t psiIndex=j+nComplexBands;
      psi [psiIndex]=c*val_r-s*val_i;
      vl_x[psiIndex]=c*gX_r-s*gX_i;
      vl_y[psiIndex]=c*gY_r-s*gY_i;
      vl_z[psiIndex]=c*gZ_r-s*gZ_i;
      vl_l[psiIndex]=c*lap_r-s*lap_i;
    }
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& vgl)
  {
    const PointType& r=P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl_soa(r,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);
    const ST x=r[0], y=r[1], z=r[2];

    const ST* restrict k0=myKcart.data(0);
    const ST* restrict k1=myKcart.data(1);
    const ST* restrict k2=myKcart.data(2);

    const ST* restrict g0=myG.data(0);
    const ST* restrict g1=myG.data(1);
    const ST* restrict g2=myG.data(2);
    const ST* restrict h00=myH.data(0);
    const ST* restrict h01=myH.data(1);
    const ST* restrict h02=myH.data(2);
    const ST* restrict h11=myH.data(3);
    const ST* restrict h12=myH.data(4);
    const ST* restrict h22=myH.data(5);

    const size_t N=kPoints.size();

    #pragma omp simd
    for (size_t j=0; j<nComplexBands; j++)
    {
      int jr=j<<1;
      int ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[jr]+g01*g1[jr]+g02*g2[jr];
      const ST dY_r = g10*g0[jr]+g11*g1[jr]+g12*g2[jr];
      const ST dZ_r = g20*g0[jr]+g21*g1[jr]+g22*g2[jr];

      const ST dX_i = g00*g0[ji]+g01*g1[ji]+g02*g2[ji];
      const ST dY_i = g10*g0[ji]+g11*g1[ji]+g12*g2[ji];
      const ST dZ_i = g20*g0[ji]+g21*g1[ji]+g22*g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;

      const size_t psiIndex=first_spo+jr;

      psi[psiIndex]  =c*val_r-s*val_i;
      dpsi[psiIndex][0]  =c*gX_r-s*gX_i;
      dpsi[psiIndex][1]  =c*gY_r-s*gY_i;
      dpsi[psiIndex][2]  =c*gZ_r-s*gZ_i;

      psi[psiIndex+1]=c*val_i+s*val_r;
      dpsi[psiIndex+1][0]=c*gX_i+s*gX_r;
      dpsi[psiIndex+1][1]=c*gY_i+s*gY_r;
      dpsi[psiIndex+1][2]=c*gZ_i+s*gZ_r;

      const ST h_xx_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g00,g01,g02,g00,g01,g02)+kX*(gX_i+dX_i);
      const ST h_xy_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g00,g01,g02,g10,g11,g12)+kX*(gY_i+dY_i);
      const ST h_xz_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g00,g01,g02,g20,g21,g22)+kX*(gZ_i+dZ_i);
      const ST h_yx_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g10,g11,g12,g00,g01,g02)+kY*(gX_i+dX_i);
      const ST h_yy_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g10,g11,g12,g10,g11,g12)+kY*(gY_i+dY_i);
      const ST h_yz_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g10,g11,g12,g20,g21,g22)+kY*(gZ_i+dZ_i);
      const ST h_zx_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g20,g21,g22,g00,g01,g02)+kZ*(gX_i+dX_i);
      const ST h_zy_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g20,g21,g22,g10,g11,g12)+kZ*(gY_i+dY_i);
      const ST h_zz_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g20,g21,g22,g20,g21,g22)+kZ*(gZ_i+dZ_i);

      const ST h_xx_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g00,g01,g02,g00,g01,g02)-kX*(gX_r+dX_r);
      const ST h_xy_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g00,g01,g02,g10,g11,g12)-kX*(gY_r+dY_r);
      const ST h_xz_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g00,g01,g02,g20,g21,g22)-kX*(gZ_r+dZ_r);
      const ST h_yx_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g10,g11,g12,g00,g01,g02)-kY*(gX_r+dX_r);
      const ST h_yy_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g10,g11,g12,g10,g11,g12)-kY*(gY_r+dY_r);
      const ST h_yz_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g10,g11,g12,g20,g21,g22)-kY*(gZ_r+dZ_r);
      const ST h_zx_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g20,g21,g22,g00,g01,g02)-kZ*(gX_r+dX_r);
      const ST h_zy_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g20,g21,g22,g10,g11,g12)-kZ*(gY_r+dY_r);
      const ST h_zz_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g20,g21,g22,g20,g21,g22)-kZ*(gZ_r+dZ_r);

      grad_grad_psi[psiIndex][0]=c*h_xx_r-s*h_xx_i;
      grad_grad_psi[psiIndex][1]=c*h_xy_r-s*h_xy_i;
      grad_grad_psi[psiIndex][2]=c*h_xz_r-s*h_xz_i;
      grad_grad_psi[psiIndex][3]=c*h_yx_r-s*h_yx_i;
      grad_grad_psi[psiIndex][4]=c*h_yy_r-s*h_yy_i;
      grad_grad_psi[psiIndex][5]=c*h_yz_r-s*h_yz_i;
      grad_grad_psi[psiIndex][6]=c*h_zx_r-s*h_zx_i;
      grad_grad_psi[psiIndex][7]=c*h_zy_r-s*h_zy_i;
      grad_grad_psi[psiIndex][8]=c*h_zz_r-s*h_zz_i;

      grad_grad_psi[psiIndex+1][0]=c*h_xx_i+s*h_xx_r;
      grad_grad_psi[psiIndex+1][1]=c*h_xy_i+s*h_xy_r;
      grad_grad_psi[psiIndex+1][2]=c*h_xz_i+s*h_xz_r;
      grad_grad_psi[psiIndex+1][3]=c*h_yx_i+s*h_yx_r;
      grad_grad_psi[psiIndex+1][4]=c*h_yy_i+s*h_yy_r;
      grad_grad_psi[psiIndex+1][5]=c*h_yz_i+s*h_yz_r;
      grad_grad_psi[psiIndex+1][6]=c*h_zx_i+s*h_zx_r;
      grad_grad_psi[psiIndex+1][7]=c*h_zy_i+s*h_zy_r;
      grad_grad_psi[psiIndex+1][8]=c*h_zz_i+s*h_zz_r;
    }

    #pragma omp simd
    for (size_t j=nComplexBands; j<N; j++)
    {
      int jr=j<<1;
      int ji=jr+1;

      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val_r=myV[jr];
      const ST val_i=myV[ji];

      //phase
      ST s, c;
      sincos(-(x*kX+y*kY+z*kZ),&s,&c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[jr]+g01*g1[jr]+g02*g2[jr];
      const ST dY_r = g10*g0[jr]+g11*g1[jr]+g12*g2[jr];
      const ST dZ_r = g20*g0[jr]+g21*g1[jr]+g22*g2[jr];

      const ST dX_i = g00*g0[ji]+g01*g1[ji]+g02*g2[ji];
      const ST dY_i = g10*g0[ji]+g11*g1[ji]+g12*g2[ji];
      const ST dZ_i = g20*g0[ji]+g21*g1[ji]+g22*g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r=dX_r+val_i*kX;
      const ST gY_r=dY_r+val_i*kY;
      const ST gZ_r=dZ_r+val_i*kZ;
      const ST gX_i=dX_i-val_r*kX;
      const ST gY_i=dY_i-val_r*kY;
      const ST gZ_i=dZ_i-val_r*kZ;

      const size_t psiIndex=first_spo+nComplexBands+j;

      psi[psiIndex]  =c*val_r-s*val_i;
      dpsi[psiIndex][0]  =c*gX_r-s*gX_i;
      dpsi[psiIndex][1]  =c*gY_r-s*gY_i;
      dpsi[psiIndex][2]  =c*gZ_r-s*gZ_i;

      const ST h_xx_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g00,g01,g02,g00,g01,g02)+kX*(gX_i+dX_i);
      const ST h_xy_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g00,g01,g02,g10,g11,g12)+kX*(gY_i+dY_i);
      const ST h_xz_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g00,g01,g02,g20,g21,g22)+kX*(gZ_i+dZ_i);
      const ST h_yx_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g10,g11,g12,g00,g01,g02)+kY*(gX_i+dX_i);
      const ST h_yy_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g10,g11,g12,g10,g11,g12)+kY*(gY_i+dY_i);
      const ST h_yz_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g10,g11,g12,g20,g21,g22)+kY*(gZ_i+dZ_i);
      const ST h_zx_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g20,g21,g22,g00,g01,g02)+kZ*(gX_i+dX_i);
      const ST h_zy_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g20,g21,g22,g10,g11,g12)+kZ*(gY_i+dY_i);
      const ST h_zz_r=v_m_v(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],g20,g21,g22,g20,g21,g22)+kZ*(gZ_i+dZ_i);

      const ST h_xx_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g00,g01,g02,g00,g01,g02)-kX*(gX_r+dX_r);
      const ST h_xy_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g00,g01,g02,g10,g11,g12)-kX*(gY_r+dY_r);
      const ST h_xz_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g00,g01,g02,g20,g21,g22)-kX*(gZ_r+dZ_r);
      const ST h_yx_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g10,g11,g12,g00,g01,g02)-kY*(gX_r+dX_r);
      const ST h_yy_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g10,g11,g12,g10,g11,g12)-kY*(gY_r+dY_r);
      const ST h_yz_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g10,g11,g12,g20,g21,g22)-kY*(gZ_r+dZ_r);
      const ST h_zx_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g20,g21,g22,g00,g01,g02)-kZ*(gX_r+dX_r);
      const ST h_zy_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g20,g21,g22,g10,g11,g12)-kZ*(gY_r+dY_r);
      const ST h_zz_i=v_m_v(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],g20,g21,g22,g20,g21,g22)-kZ*(gZ_r+dZ_r);

      grad_grad_psi[psiIndex][0]=c*h_xx_r-s*h_xx_i;
      grad_grad_psi[psiIndex][1]=c*h_xy_r-s*h_xy_i;
      grad_grad_psi[psiIndex][2]=c*h_xz_r-s*h_xz_i;
      grad_grad_psi[psiIndex][3]=c*h_yx_r-s*h_yx_i;
      grad_grad_psi[psiIndex][4]=c*h_yy_r-s*h_yy_i;
      grad_grad_psi[psiIndex][5]=c*h_yz_r-s*h_yz_i;
      grad_grad_psi[psiIndex][6]=c*h_zx_r-s*h_zx_i;
      grad_grad_psi[psiIndex][7]=c*h_zy_r-s*h_zy_i;
      grad_grad_psi[psiIndex][8]=c*h_zz_r-s*h_zz_i;
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgh(r,psi,dpsi,grad_grad_psi);
  }
};

}
#endif
