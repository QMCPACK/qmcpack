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
    
    
/** @file SplineC2CSoA.h
 *
 * Adoptor classes to handle complex-to-(real,complex) with arbitrary precision
 */
#ifndef QMCPLUSPLUS_EINSPLINE_C2C_SOA_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_C2C_SOA_ADOPTOR_H

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
struct SplineC2CSoA: public SplineAdoptorBase<ST,3>
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

  SplineC2CSoA(): BaseType(),SplineInst(nullptr), MultiSpline(nullptr)
  {
    this->is_complex=true;
    this->AdoptorName="SplineC2CSoA";
    this->KeyWord="C2RSoA";
  }

  ///** copy the base property */
  //SplineC2CSoA(BaseType& rhs): BaseType(rhs)
  //{
  //  this->is_complex=true;
  //  this->AdoptorName="SplineC2CSoA";
  //  this->KeyWord="C2RSoA";
  //}

  SplineC2CSoA(const SplineC2CSoA& a):
    SplineAdoptorBase<ST,3>(a),SplineInst(a.SplineInst),MultiSpline(nullptr), 
    nComplexBands(a.nComplexBands),mKK(a.mKK), myKcart(a.myKcart)
  {
    const size_t n=a.myL.size();
    myV.resize(n); myG.resize(n); myL.resize(n); myH.resize(n);
  }

  ~SplineC2CSoA() 
  { 
    if(MultiSpline != nullptr) delete SplineInst;
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseType::init_base(n);
    myV.resize(2*n);
    myG.resize(2*n);
    myL.resize(2*n);
    myH.resize(2*n);
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
    const size_t nk=kPoints.size();
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
    SplineInst->copy_spline(spline_r,2*ispline  ,BaseOffset, BaseN);
    SplineInst->copy_spline(spline_i,2*ispline+1,BaseOffset, BaseN);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    VectorViewer<ST> v_r(psi_r,0), v_i(psi_i,0);
    SplineInst->set(2*ispline  ,v_r);
    SplineInst->set(2*ispline+1,v_i);
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
    typedef typename psi::value_type ComplexT;
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
    for (size_t j=0,psiIndex=first_spo; j<N; j++, psiIndex++)
    {
      const ST val_r=myV[2*j  ];
      const ST val_i=myV[2*j+1];
      psi[psiIndex]=ComplexT( val_r*CosV[j]-val_i*SinV[j], val_i*CosV[j]+val_r*SinV[j]);
    }
#else
    ST s, c;
#pragma omp simd private(s,c)
    for (size_t j=0,psiIndex=first_spo; j<N; j++, psiIndex++)
    {
      const ST val_r=myV[2*j  ];
      const ST val_i=myV[2*j+1];
      sincos(-(x*kx[j]+y*ky[j]+z*kz[j]),&s,&c);
      psi[psiIndex  ] = ComplexT(val_r*c-val_i*s,val_i*c+val_r*s);
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
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    //DO NOT USE THESE PLEASE
  }

  /** identical to assign_vgl but the output container is SoA container
   */
  template<typename VV, typename GL>
  inline void assign_vgl_soa(const PointType& r, int bc_sign, VV& psi, GL& dpsi)
  {
    typedef typename VV::value_type ComplexT;

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

    ComplexT* restrict vg_x=dpsi.data(0); ASSUME_ALIGNED(vl_x);
    ComplexT* restrict vg_y=dpsi.data(1); ASSUME_ALIGNED(vl_y);
    ComplexT* restrict vg_z=dpsi.data(2); ASSUME_ALIGNED(vl_z);
    ComplexT* restrict vl_l=dpsi.data(3); ASSUME_ALIGNED(vl_l);

    for (size_t j=0, psiIndex=first_spo; j<N; j++,psiIndex++)
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
      psi[psiIndex ]=ComplexT(c*val_r-s*val_i,c*val_i+s*val_r);
      vg_x[psiIndex]=ComplexT(c*gX_r -s*gX_i, c*gX_i +s*gX_r);
      vg_y[psiIndex]=ComplexT(c*gY_r -s*gY_i, c*gY_i +s*gY_r);
      vg_z[psiIndex]=ComplexT(c*gZ_r -s*gZ_i, c*gZ_i +s*gZ_r);
      vl_l[psiIndex]=ComplexT(c*lap_r-s*lap_i,c*lap_i+s*lap_r);
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
    assign_vgl_soa(r,psi,dpsi);
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

}
#endif
