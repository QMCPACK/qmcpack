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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EINSPLINE_R2RSOA_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_R2RSOA_ADOPTOR_H

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
struct SplineR2RSoA: public SplineAdoptorBase<ST,3>
{
  static const int D=3;
  bool IsGamma;
  using BaseType=SplineAdoptorBase<ST,3>;
  using SplineType=typename bspline_traits<ST,3>::SplineType;
  using BCType=typename bspline_traits<ST,3>::BCType;
  using PointType=typename BaseType::PointType;
  using SingleSplineType=typename BaseType::SingleSplineType;

  using vContainer_type=Vector<ST,aligned_allocator<ST> >;
  using gContainer_type=VectorSoaContainer<ST,3>;
  using hContainer_type=VectorSoaContainer<ST,6>;

  using BaseType::first_spo;
  using BaseType::last_spo;
  using SplineAdoptorBase<ST,D>::HalfG;
  using BaseType::GGt;
  using BaseType::PrimLattice;
  using BaseType::kPoints;
  using BaseType::offset_cplx;
  using BaseType::offset_real;

  ///number of points of the original grid
  int BaseN[3];
  ///offset of the original grid, always 0
  int BaseOffset[3];
  ///multi bspline set
  MultiBspline<ST>* SplineInst;
  ///expose the pointer to reuse the reader and only assigned with create_spline
  SplineType* MultiSpline;

  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;

  SplineR2RSoA(): BaseType(),SplineInst(nullptr), MultiSpline(nullptr)
  {
    this->is_complex=false;
    this->is_soa_ready=true;
    this->AdoptorName="SplineR2RSoAAdoptor";
    this->KeyWord="SplineR2RSoA";
  }

  ///** copy the base property */
  //SplineC2CSoA(BaseType& rhs): BaseType(rhs)
  //{
  //  this->is_complex=true;
  //  this->AdoptorName="SplineC2CSoA";
  //  this->KeyWord="C2RSoA";
  //}

  SplineR2RSoA(const SplineR2RSoA& a):
    SplineAdoptorBase<ST,3>(a),SplineInst(a.SplineInst),MultiSpline(nullptr)
  {
    const size_t n=a.myV.size();
    myV.resize(n); myG.resize(n); myL.resize(n); myH.resize(n);
  }

  ~SplineR2RSoA()
  {
    if(MultiSpline != nullptr) delete SplineInst;
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseType::init_base(n);
    const size_t npad=getAlignedSize<ST>(n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);

    IsGamma=( (HalfG[0]==0) && (HalfG[1]==0) && (HalfG[2]==0));
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
    offset_real.resize(Nbandgroups+1,0);
    FairDivideLow(Nbands,Nbandgroups,offset_real);
    gatherv(comm, MultiSpline, MultiSpline->z_stride, offset_real);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
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

  void create_spline(TinyVector<int,D>& mesh, int n)
  {
    Ugrid xyz_grid[D];
    BCType xyz_bc[D];
    for(int i=0; i<D; ++i)
    {
      xyz_grid[i].start = 0.0;
      xyz_grid[i].end = 1.0;
      xyz_grid[i].num = mesh[i];
      xyz_bc[i].lCode=xyz_bc[i].rCode=(HalfG[i])? ANTIPERIODIC:PERIODIC;
      BaseOffset[i]=0;
      BaseN[i]=xyz_grid[i].num+3;
    }
    SplineInst=new MultiBspline<ST>();
    SplineInst->create(xyz_grid,xyz_bc,n);
    MultiSpline=SplineInst->spline_m;
    qmc_common.memory_allocated += MultiSpline->coefs_size*sizeof(ST);
  }

  inline void flush_zero()
  {
    SplineInst->flush_zero();
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    SplineInst->copy_spline(spline_r, ispline, BaseOffset, BaseN);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    Vector<ST> v_r(psi_r,0);
    SplineInst->set(ispline, v_r);
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

  /** convert postion in PrimLattice unit and return sign */
  inline int convertPos(const PointType& r, PointType& ru)
  {
    ru=PrimLattice.toUnit(r);
    int bc_sign=0;
    for(int i=0; i<D; i++)
      if( -std::numeric_limits<ST>::epsilon() < ru[i] && ru[i] < 0 )
        ru[i] = ST(0.0);
      else
      {
        ST img = std::floor(ru[i]);
        ru[i] -= img;
        bc_sign += HalfG[i] * (int)img;
      }
    return bc_sign;
  }

  template<typename VV>
  inline void assign_v(int bc_sign, const vContainer_type& myV, VV& psi)
  {
    if (bc_sign & 1)
      for(size_t psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=-myV[j];
    else
      for(size_t psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=myV[j];
  }

  inline TT evaluate_dot(const ParticleSet& P, const int iat, const TT* restrict arow, ST* scratch, bool compute_spline=true)
  {
    Vector<ST> vtmp(scratch,myV.size());
    PointType ru;
    int bc_sign=convertPos(P.activeR(iat),ru);
    if(compute_spline) SplineInst->evaluate(ru,vtmp);

    TT res=TT();
    if (bc_sign & 1)
      #pragma omp simd reduction(+:res)
      for(size_t psiIndex=first_spo; psiIndex<last_spo; ++psiIndex)
        res -= vtmp[psiIndex-first_spo]*arow[psiIndex];
    else
      #pragma omp simd reduction(+:res)
      for(size_t psiIndex=first_spo; psiIndex<last_spo; ++psiIndex)
        res += vtmp[psiIndex-first_spo]*arow[psiIndex];
    return res;
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.activeR(iat);
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate(ru,myV);
    assign_v(bc_sign,myV,psi);
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

  template<typename VV, typename GV>
  inline void assign_vgl(int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);
    const ST symGG[6]={GGt[0],GGt[1]+GGt[3],GGt[2]+GGt[6],GGt[4],GGt[5]+GGt[7],GGt[8]};

    const ST* restrict g0=myG.data(0);
    const ST* restrict g1=myG.data(1);
    const ST* restrict g2=myG.data(2);
    const ST* restrict h00=myH.data(0);
    const ST* restrict h01=myH.data(1);
    const ST* restrict h02=myH.data(2);
    const ST* restrict h11=myH.data(3);
    const ST* restrict h12=myH.data(4);
    const ST* restrict h22=myH.data(5);

    if (bc_sign & 1)
    {
      #pragma omp simd
      for(size_t psiIndex=first_spo; psiIndex<last_spo; ++psiIndex)
      {
        const size_t j=psiIndex-first_spo;
        psi[psiIndex]=-myV[j];
        dpsi[psiIndex][0]=-(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        dpsi[psiIndex][1]=-(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        dpsi[psiIndex][2]=-(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        d2psi[psiIndex]=-SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
      }
    }
    else
    {
      #pragma omp simd
      for(size_t psiIndex=first_spo; psiIndex<last_spo; ++psiIndex)
      {
        const size_t j=psiIndex-first_spo;
        psi[psiIndex]=myV[j];
        dpsi[psiIndex][0]=(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        dpsi[psiIndex][1]=(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        dpsi[psiIndex][2]=(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        d2psi[psiIndex]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
      }
    }
  }

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  template<typename VV, typename GV>
  inline void assign_vgl_from_l(int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    const ST* restrict g0=myG.data(0);
    const ST* restrict g1=myG.data(1);
    const ST* restrict g2=myG.data(2);

    if (bc_sign & 1)
    {
      #pragma omp simd
      for(int psiIndex=first_spo; psiIndex<last_spo; ++psiIndex)
      {
        const size_t j=psiIndex-first_spo;
        psi[psiIndex]=-myV[j];
        dpsi[psiIndex][0]=-g0[j];
        dpsi[psiIndex][1]=-g1[j];
        dpsi[psiIndex][2]=-g2[j];
        d2psi[psiIndex]=-myL[j];
      }
    }
    else
    {
      #pragma omp simd
      for(int psiIndex=first_spo; psiIndex<last_spo; ++psiIndex)
      {
        const size_t j=psiIndex-first_spo;
        psi[psiIndex]=myV[j];
        dpsi[psiIndex][0]=g0[j];
        dpsi[psiIndex][1]=g1[j];
        dpsi[psiIndex][2]=g2[j];
        d2psi[psiIndex]=myL[j];
      }
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.activeR(iat);
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl(bc_sign,psi,dpsi,d2psi);
  }

  /** identical to assign_vgl but the output container is SoA container
   */
  template<typename VGL>
  inline void assign_vgl_soa(int bc_sign, VGL& vgl)
  {
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);
    const ST symGG[6]={GGt[0],GGt[1]+GGt[3],GGt[2]+GGt[6],GGt[4],GGt[5]+GGt[7],GGt[8]};

    const ST* restrict g0=myG.data(0);
    const ST* restrict g1=myG.data(1);
    const ST* restrict g2=myG.data(2);
    const ST* restrict h00=myH.data(0);
    const ST* restrict h01=myH.data(1);
    const ST* restrict h02=myH.data(2);
    const ST* restrict h11=myH.data(3);
    const ST* restrict h12=myH.data(4);
    const ST* restrict h22=myH.data(5);

    TT* restrict psi =vgl.data(0)+first_spo; ASSUME_ALIGNED(psi);
    TT* restrict vg_x=vgl.data(1)+first_spo; ASSUME_ALIGNED(vg_x);
    TT* restrict vg_y=vgl.data(2)+first_spo; ASSUME_ALIGNED(vg_y);
    TT* restrict vg_z=vgl.data(3)+first_spo; ASSUME_ALIGNED(vg_z);
    TT* restrict vl_l=vgl.data(4)+first_spo; ASSUME_ALIGNED(vl_l);

    const size_t N=last_spo-first_spo;
    if (bc_sign & 1)
    {
      #pragma omp simd
      for(int j=0; j<N; ++j)
      {
        psi [j]=-myV[j];
        vg_x[j]=-(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        vg_y[j]=-(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        vg_z[j]=-(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        vl_l[j]=-SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
      }
    }
    else
    {
      #pragma omp simd
      for(int j=0; j<N; ++j)
      {
        psi [j]=myV[j];
        vg_x[j]=(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        vg_y[j]=(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        vg_z[j]=(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        vl_l[j]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
      }
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
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl_soa(bc_sign,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const ST cone = (bc_sign &1)? -1:1;
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);

    const ST* restrict g0=myG.data(0);
    const ST* restrict g1=myG.data(1);
    const ST* restrict g2=myG.data(2);
    const ST* restrict h00=myH.data(0);
    const ST* restrict h01=myH.data(1);
    const ST* restrict h02=myH.data(2);
    const ST* restrict h11=myH.data(3);
    const ST* restrict h12=myH.data(4);
    const ST* restrict h22=myH.data(5);

    const size_t N=last_spo-first_spo;

    #pragma omp simd
    for (size_t j=0; j<N; ++j)
    {
      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00*g0[j]+g01*g1[j]+g02*g2[j];
      const ST dY_r = g10*g0[j]+g11*g1[j]+g12*g2[j];
      const ST dZ_r = g20*g0[j]+g21*g1[j]+g22*g2[j];

      const size_t psiIndex=j+first_spo;
      psi[psiIndex]    =cone*myV[j];
      dpsi[psiIndex][0]=cone*dX_r;
      dpsi[psiIndex][1]=cone*dY_r;
      dpsi[psiIndex][2]=cone*dZ_r;

      const ST h_xx_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g00,g01,g02,g00,g01,g02);
      const ST h_xy_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g00,g01,g02,g10,g11,g12);
      const ST h_xz_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g00,g01,g02,g20,g21,g22);
      const ST h_yx_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g10,g11,g12,g00,g01,g02);
      const ST h_yy_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g10,g11,g12,g10,g11,g12);
      const ST h_yz_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g10,g11,g12,g20,g21,g22);
      const ST h_zx_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g20,g21,g22,g00,g01,g02);
      const ST h_zy_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g20,g21,g22,g10,g11,g12);
      const ST h_zz_r=v_m_v(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],g20,g21,g22,g20,g21,g22);

      grad_grad_psi[psiIndex][0]=cone*h_xx_r;
      grad_grad_psi[psiIndex][1]=cone*h_xy_r;
      grad_grad_psi[psiIndex][2]=cone*h_xz_r;
      grad_grad_psi[psiIndex][3]=cone*h_yx_r;
      grad_grad_psi[psiIndex][4]=cone*h_yy_r;
      grad_grad_psi[psiIndex][5]=cone*h_yz_r;
      grad_grad_psi[psiIndex][6]=cone*h_zx_r;
      grad_grad_psi[psiIndex][7]=cone*h_zy_r;
      grad_grad_psi[psiIndex][8]=cone*h_zz_r;
    }
  }


  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.activeR(iat);
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgh(bc_sign,psi,dpsi,grad_grad_psi);
  }
};


}
#endif

