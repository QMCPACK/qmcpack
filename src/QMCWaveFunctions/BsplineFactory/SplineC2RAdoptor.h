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

#include <Numerics/VectorViewer.h>
#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>

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
  using PointType=typename BaseType::PointType; 
  using SingleSplineType=typename BaseType::SingleSplineType;

  using vContainer_type=aligned_vector<ST>;
  using gContainer_type=VectorSoaContainer<ST,3>;
  using hContainer_type=VectorSoaContainer<ST,6>;

  using BaseType::first_spo;
  using BaseType::last_spo;
  using BaseType::GGt;
  using BaseType::PrimLattice;
  using BaseType::kPoints;
  using BaseType::MakeTwoCopies;

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

  //vContainer_type  KdotR;
  //vContainer_type  CosV;
  //vContainer_type  SinV;
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
    this->AdoptorName="SplineC2RSoA";
    this->KeyWord="C2RSoA";
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
    int iband=this->BandIndexMap[ispline];
    SplineInst->copy_spline(spline_r,2*iband  ,BaseOffset, BaseN);
    SplineInst->copy_spline(spline_i,2*iband+1,BaseOffset, BaseN);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    VectorViewer<ST> v_r(psi_r,0), v_i(psi_i,0);
    int iband=this->BandIndexMap[ispline];
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

  template<typename VV>
  inline void assign_v(const PointType& r, int bc_sign, VV& psi)
  {
    const size_t N=kPoints.size();
    const ST x=r[0], y=r[1], z=r[2];
    const ST* restrict kx=myKcart.data(0);
    const ST* restrict ky=myKcart.data(1);
    const ST* restrict kz=myKcart.data(2);
#if defined(USE_VECTOR_ML)
    {
      for(size_t j=0; j<N; ++j)
        KdotR[j]=-(x*kx[j]+y*ky[j]+z*kz[j]);
      eval_e2iphi(nk,KdotR.data(),CosV.data(),SinV.data());
    }
    for (size_t j=0,psiIndex=first_spo; j<nComplexBands; j++, psiIndex+=2)
    {
      const ST val_r=myV[2*j  ];
      const ST val_i=myV[2*j+1];
      psi[psiIndex  ] = val_r*CosV[j]-val_i*SinV[j];
      psi[psiIndex+1] = val_i*CosV[j]+val_r*SinV[j];
    }

    for (size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands; j<N; j++,psiIndex++)
    {
      const ST val_r=myV[2*j  ];
      const ST val_i=myV[2*j+1];
      psi[psiIndex  ] = val_r*CosV[j]-val_i*SinV[j];
    }
#else
    ST s, c;
#pragma omp simd private(s,c)
    for (size_t j=0,psiIndex=first_spo; j<nComplexBands; j++, psiIndex+=2)
    {
      const ST val_r=myV[2*j  ];
      const ST val_i=myV[2*j+1];
      sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
      psi[psiIndex  ] = val_r*c-val_i*s;
      psi[psiIndex+1] = val_i*c+val_r*s;
    }

#pragma omp simd private(s,c)
    for (size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands; j<N; j++,psiIndex++)
    {
      const ST val_r=myV[2*j  ];
      const ST val_i=myV[2*j+1];
      sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
      psi[psiIndex  ] = val_r*c-val_i*s;
    }
#endif
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate(ru,myV);
    assign_v(r,0,psi);
  }

  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    CONSTEXPR ST zero(0);
    CONSTEXPR ST two(2);
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
    const size_t nsplines=myL.size();
#if defined(PRECOMPUTE_L)
    for(size_t j=0; j<nsplines; ++j)
    {
      myL[j]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],GGt.data());
    }
#endif

#pragma omp simd
    for (size_t j=0, psiIndex=first_spo; j<nComplexBands; j++,psiIndex+=2)
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
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],GGt.data());
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],GGt.data());
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif  

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

    const size_t nComputed=2*nComplexBands;
#pragma omp simd
    for (size_t j=nComplexBands,psiIndex=first_spo+nComputed; j<N; j++,psiIndex++)
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
      psi[psiIndex  ]=c*val_r-s*val_i;
      //this will be fixed later
      dpsi[psiIndex  ][0]=c*gX_r-s*gX_i;
      dpsi[psiIndex  ][1]=c*gY_r-s*gY_i;
      dpsi[psiIndex  ][2]=c*gZ_r-s*gZ_i;

#if defined(PRECOMPUTE_L)
      const ST lap_r=myL[jr]+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=myL[ji]+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#else
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],GGt.data());
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],GGt.data());
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif
      d2psi[psiIndex  ]=c*lap_r-s*lap_i;
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl(r,0,psi,dpsi,d2psi);
  }

  /** identical to assign_vgl but the output container is SoA container
   */
  template<typename VV, typename GL>
  inline void assign_vgl_soa(const PointType& r, int bc_sign, VV& psi, GL& dpsi)
  {
    const int N=kPoints.size();
    CONSTEXPR ST zero(0);
    CONSTEXPR ST two(2);
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
    const size_t nsplines=myL.size();
#if defined(PRECOMPUTE_L)
    for(size_t j=0; j<nsplines; ++j)
    {
      myL[j]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],GGt.data());
    }
#endif

    TT* restrict vl_x=dpsi.data(0);
    TT* restrict vl_y=dpsi.data(1);
    TT* restrict vl_z=dpsi.data(2);
    TT* restrict vl_l=dpsi.data(3);
#pragma omp simd
    for (size_t j=0, psiIndex=first_spo; j<nComplexBands; j++,psiIndex+=2)
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
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],GGt.data());
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],GGt.data());
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif
      
      //this will be fixed later
      psi[psiIndex   ]=c*val_r-s*val_i;
      psi[psiIndex+1 ]=c*val_i+s*val_r;
      vl_x[psiIndex  ]=c*gX_r-s*gX_i;
      vl_x[psiIndex+1]=c*gX_i+s*gX_r;
      vl_y[psiIndex  ]=c*gY_r-s*gY_i;
      vl_y[psiIndex+1]=c*gY_i+s*gY_r;
      vl_z[psiIndex  ]=c*gZ_r-s*gZ_i;
      vl_z[psiIndex+1]=c*gZ_i+s*gZ_r;
      vl_l[psiIndex  ]=c*lap_r-s*lap_i;
      vl_l[psiIndex+1]=c*lap_i+s*lap_r;
    }

#pragma omp simd
    for (size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands; j<N; j++,psiIndex++)
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
      const ST lcart_r=SymTrace(h00[jr],h01[jr],h02[jr],h11[jr],h12[jr],h22[jr],GGt.data());
      const ST lcart_i=SymTrace(h00[ji],h01[ji],h02[ji],h11[ji],h12[ji],h22[ji],GGt.data());
      const ST lap_r=lcart_r+mKK[j]*val_r+two*(kX*dX_i+kY*dY_i+kZ*dZ_i);
      const ST lap_i=lcart_i+mKK[j]*val_i-two*(kX*dX_r+kY*dY_r+kZ*dZ_r);
#endif

      psi[psiIndex   ]=c*val_r-s*val_i;
      vl_x[psiIndex  ]=c*gX_r-s*gX_i;
      vl_y[psiIndex  ]=c*gY_r-s*gY_i;
      vl_z[psiIndex  ]=c*gZ_r-s*gZ_i;
      vl_l[psiIndex  ]=c*lap_r-s*lap_i;
    }
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VV, typename GL>
  inline void evaluate_vgl_combo(const PointType& r, VV& psi, GL& dpsi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl_soa(r,0,psi,dpsi);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    //missing
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    //missing
  }
};

#if 0
/** adoptor class to match std::complex<ST> spline with TT real SPOs
 *
 * Reimplement SplineC2RPackedAdoptor with internal grouping
 * This is to debug SplineC2RSoA and not is intended for use in the production run
 */
template<typename ST, typename TT, unsigned D>
struct SplineC2RAdoptor: public SplineAdoptorBase<ST,D>
{
  using SplineType      = typename einspline_traits<ST,D>::SplineType;
  using BCType          = typename einspline_traits<ST,D>::BCType;
  using PointType       = typename SplineAdoptorBase<ST,D>::PointType;
  using SingleSplineType= typename SplineAdoptorBase<ST,D>::SingleSplineType;

  using SplineAdoptorBase<ST,D>::first_spo;
  using SplineAdoptorBase<ST,D>::last_spo;
  using SplineAdoptorBase<ST,D>::GGt;
  using SplineAdoptorBase<ST,D>::PrimLattice;
  using SplineAdoptorBase<ST,D>::kPoints;
  using SplineAdoptorBase<ST,D>::MakeTwoCopies;

  typename OrbitalSetTraits<ST>::ValueVector_t     myV;
  typename OrbitalSetTraits<ST>::ValueVector_t     myL;
  typename OrbitalSetTraits<ST>::GradVector_t      myG;
  typename OrbitalSetTraits<ST>::HessVector_t      myH;
  typename OrbitalSetTraits<ST>::GradHessVector_t  myGH;

  ///Actual spline table, multi_bspline_3d_(d,s)
  SplineType *MultiSpline;

  ///number of points of the original grid
  int BaseN[3];
  ///offset of the original grid, always 0
  int BaseOffset[3];
  ///number of complex bands
  int nComplexBands;

  std::vector<ST>   KdotR;
  std::vector<ST>   CosV;
  std::vector<ST>   SinV;
  std::vector<ST>   mKK;
  std::vector<Tensor<ST,D> >  KK; //k^k

  SplineC2RAdoptor():MultiSpline(nullptr)
  {
    this->is_complex=true;
    this->AdoptorName="SplineC2RAdoptor";
    this->KeyWord="C2RPacked";
  }

  virtual ~SplineC2RAdoptor() {}

  inline void resizeStorage(int n, int nvals)
  {
    SplineAdoptorBase<ST,D>::init_base(n);
    myV.resize(2*n);
    myL.resize(2*n);
    myG.resize(2*n);
    myH.resize(2*n);
    CosV.resize(n);
    SinV.resize(n);
    KdotR.resize(n);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    //Shadow=new SplineC2RSoA<ST,TT>(*this);
    //Shadow->resizeStorage(this->nunique_orbitals,this->nunique_orbitals);
    //Shadow->create_spline(xyz_g,xyz_bc);
   
    resize_kk();
    MultiSpline=einspline::create(MultiSpline,xyz_g,xyz_bc,myV.size());
    for(int i=0; i<D; ++i)
    {
      BaseOffset[i]=0;
      BaseN[i]=xyz_g[i].num+3;
    }
    qmc_common.memory_allocated += MultiSpline->coefs_size*sizeof(ST);
  }

  inline void resize_kk()
  {
    nComplexBands=this->remap_kpoints();
    mKK.resize(kPoints.size());
    for(size_t i=0; i<kPoints.size(); ++i)
      mKK[i]=-dot(kPoints[i],kPoints[i]);
    KK.resize(kPoints.size());
    for(size_t i=0; i<kPoints.size(); ++i)
      KK[i]=outerProduct(kPoints[i],kPoints[i]);
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    int iband=this->BandIndexMap[ispline];
    einspline::set(MultiSpline, 2*iband  ,spline_r, BaseOffset, BaseN);
    einspline::set(MultiSpline, 2*iband+1,spline_i, BaseOffset, BaseN);
    //Shadow->set_spline(spline_r,spline_i,twist,ispline,level);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    int iband=this->BandIndexMap[ispline];
    einspline::set(MultiSpline, 2*iband,   psi_r);
    einspline::set(MultiSpline, 2*iband+1, psi_i);
  }

  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i, 
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
    int iband=this->BandIndexMap[ispline];
    if(mKK.empty()) resize_kk();
    einspline::set(MultiSpline, 2*iband,  spline_r, offset_l, mesh_l);
    einspline::set(MultiSpline, 2*iband+1,spline_i, offset_l, mesh_l);
  }

  bool read_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o<<"spline_" << SplineAdoptorBase<ST,D>::MyIndex;
    einspline_engine<SplineType> bigtable(MultiSpline);
    return h5f.read(bigtable,o.str().c_str());//"spline_0");
  }

  bool write_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o<<"spline_" << SplineAdoptorBase<ST,D>::MyIndex;
    einspline_engine<SplineType> bigtable(MultiSpline);
    return h5f.write(bigtable,o.str().c_str());//"spline_0");
  }

  template<typename VV>
  inline void assign_v(const PointType& r, int bc_sign, VV& psi)
  {
    const int N=kPoints.size();
    for(int j=0; j<N; ++j)
      KdotR[j]=-dot(r,kPoints[j]);
    eval_e2iphi(N,KdotR.data(),CosV.data(),SinV.data());
    //for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
    for(size_t j=0,psiIndex=first_spo; j<nComplexBands; ++j, psiIndex+=2)
    {
      const size_t jr=j<<1;
      psi[psiIndex  ] = static_cast<TT>(myV[jr  ]*CosV[j]-myV[jr+1]*SinV[j]);
      psi[psiIndex+1] = static_cast<TT>(myV[jr+1]*CosV[j]+myV[jr  ]*SinV[j]);
    }
    for (size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands;
        j<N; ++j, psiIndex++)
    {
      const size_t jr=j<<1;
      psi[psiIndex  ] = static_cast<TT>(myV[jr]*CosV[j]-myV[jr+1]*SinV[j]);
    }
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate(MultiSpline,ru,myV);
    assign_v(r,0,psi);
    //Shadow->evaluate_v(r,psi);
  }

  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    const size_t N=kPoints.size();
    CONSTEXPR ST zero(0);
    CONSTEXPR ST two(2);
    for(size_t j=0; j<N; ++j)
      KdotR[j]=-dot(r,kPoints[j]);
    eval_e2iphi(N,KdotR.data(),CosV.data(),SinV.data());
    //for(size_t j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
    for(size_t j=0,psiIndex=first_spo; j<nComplexBands; ++j, psiIndex+=2)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;
      const auto myG_r= dot(PrimLattice.G, myG[jr]);
      const auto myG_i= dot(PrimLattice.G, myG[ji]);
      const auto g_r=myG_r+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const auto g_i=myG_i-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      psi[psiIndex   ] =CosV[j]*myV[jr]-SinV[j]*myV[ji];
      psi[psiIndex+1 ] =CosV[j]*myV[ji]+SinV[j]*myV[jr];
      dpsi[psiIndex  ] =CosV[j]*g_r-SinV[j]*g_i; // multiply phase
      dpsi[psiIndex+1] =CosV[j]*g_i+SinV[j]*g_r;
      const auto myL_r= trace(myH[jr],GGt)+mKK[j]*myV[jr]+two*dot(kPoints[j],myG_i);
      const auto myL_i= trace(myH[ji],GGt)+mKK[j]*myV[ji]-two*dot(kPoints[j],myG_r);
      d2psi[psiIndex  ]=CosV[j]*myL_r-SinV[j]*myL_i;
      d2psi[psiIndex+1]=CosV[j]*myL_i+SinV[j]*myL_r;
    }

    for (size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands;
        j<N; ++j, psiIndex++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;
      const auto myG_r= dot(PrimLattice.G, myG[jr]);
      const auto myG_i= dot(PrimLattice.G, myG[ji]);
      const auto g_r=myG_r+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const auto g_i=myG_i-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      psi[psiIndex   ] =CosV[j]*myV[jr]-SinV[j]*myV[ji];
      dpsi[psiIndex  ] =CosV[j]*g_r-SinV[j]*g_i; // multiply phase
      const auto myL_r= trace(myH[jr],GGt)+mKK[j]*myV[jr]+two*dot(kPoints[j],myG_i);
      const auto myL_i= trace(myH[ji],GGt)+mKK[j]*myV[ji]-two*dot(kPoints[j],myG_r);
      d2psi[psiIndex ] =CosV[j]*myL_r-SinV[j]*myL_i;
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {

    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgl(r,0,psi,dpsi,d2psi);
    //Shadow->evaluate_vgl(r,psi,dpsi,d2psi);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    //convert to Cartesian G and Hess
    const size_t N=kPoints.size();
    for (size_t j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (size_t j=0; j<2*N; j++)
      myH[j] = dot(myH[j],GGt);
    ST s,c;
    //can easily make three independent loops
    for(size_t j=0,psiIndex=first_spo; j<nComplexBands; ++j, psiIndex+=2)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;
      const PointType g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const PointType g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      //kk=outerProduct(kPoints[j],kPoints[j]); // \f$kk=k^k \f$
      const Tensor<ST,D> h_r=myH[jr]-myV[jr]*KK[j]+outerProductSymm(kPoints[j],myG[ji]);
      const Tensor<ST,D> h_i=myH[ji]-myV[ji]*KK[j]-outerProductSymm(kPoints[j],myG[jr]);
      sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
      psi[psiIndex]=c*myV[jr]-s*myV[ji];
      psi[psiIndex+1]=s*myV[jr]+c*myV[ji];

      for(size_t idim=0; idim<D; ++idim)
        dpsi[psiIndex][idim]=c*g_r[idim]-s*g_i[idim];
      for(size_t idim=0; idim<D; ++idim)
        dpsi[psiIndex+1][idim]=c*g_i[idim]+s*g_r[idim];

      for(size_t t=0; t<D*D; ++t)
        grad_grad_psi[psiIndex](t)=c*h_r(t)-s*h_i(t);
      for(size_t t=0; t<D*D; ++t)
        grad_grad_psi[psiIndex+1](t)=s*h_r(t)+c*h_i(t);
    }
    for(size_t j=nComplexBands,psiIndex=first_spo+2*nComplexBands; j<N; ++j, psiIndex++)
    {
      const size_t jr=j<<1;
      const size_t ji=jr+1;
      auto g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      auto g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      //kk=outerProduct(kPoints[j],kPoints[j]); // \f$kk=k^k \f$
      auto h_r=myH[jr]-myV[jr]*KK[j]+outerProductSymm(kPoints[j],myG[ji]);
      auto h_i=myH[ji]-myV[ji]*KK[j]-outerProductSymm(kPoints[j],myG[jr]);
      sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
      psi[psiIndex]=c*myV[jr]-s*myV[ji];
      for(size_t idim=0; idim<D; ++idim)
        dpsi[psiIndex][idim]=c*g_r[idim]-s*g_i[idim];
      for(size_t t=0; t<D*D; ++t)
        grad_grad_psi[psiIndex](t)=c*h_r(t)-s*h_i(t);
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgh(r,0,psi,dpsi,grad_grad_psi);
  }
};
#endif

}
#endif
