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
struct SplineR2RSoA: public SplineAdoptorBase<ST,3>
{
  static const int D=3;
  bool IsGamma;
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
  using SplineAdoptorBase<ST,D>::HalfG;
  using BaseType::GGt;
  using BaseType::PrimLattice;

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
    this->KeyWord="R2RSoA";
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


  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    SplineInst->copy_spline(spline_r,ispline  ,BaseOffset, BaseN);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    VectorViewer<ST> v_r(psi_r,0);
    SplineInst->set(ispline  ,v_r);
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
    {
      ST img = std::floor(ru[i]);
      ru[i] -= img;
      bc_sign += HalfG[i] * (int)img;
    }
    return bc_sign;
  }

  template<typename VV>
  inline void assign_v(int bc_sign, VV& psi)
  {
    if (bc_sign & 1)
      for(size_t psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=-myV[j];
    else
      for(size_t psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=myV[j];
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate(ru,myV);
    assign_v(bc_sign,psi);
  }

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
#pragma simd
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      {
        psi[psiIndex]=-myV[j];
        dpsi[psiIndex][0]=-(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        dpsi[psiIndex][1]=-(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        dpsi[psiIndex][2]=-(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        d2psi[psiIndex]=-SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
      }
    }
    else
    {
#pragma simd
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      {
        psi[psiIndex]=myV[j];
        dpsi[psiIndex][0]=(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        dpsi[psiIndex][1]=(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        dpsi[psiIndex][2]=(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        d2psi[psiIndex]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],symGG);
      }
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
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
  inline void evaluate_vgl_combo(const PointType& r, VGL& vgl)
  {
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl_soa(bc_sign,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
#if 0
    const ST g00=PrimLattice.G(0), g01=PrimLattice.G(1), g02=PrimLattice.G(2),
             g10=PrimLattice.G(3), g11=PrimLattice.G(4), g12=PrimLattice.G(5),
             g20=PrimLattice.G(6), g21=PrimLattice.G(7), g22=PrimLattice.G(8);
    const ST symGG[6]={GGt[0],GGt[1]+GGt[3],GGt[2]+GGt[6],GGt[4],GGt[5]+GGt[7],GGt[8]};

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

    ST* restrict psi =psi.data(0)+first_spo; ASSUME_ALIGNED(psi);
    ST* restrict vg_x=dpsi.data(0)+first_spo; ASSUME_ALIGNED(vg_x);
    ST* restrict vg_y=dpsi.data(1)+first_spo; ASSUME_ALIGNED(vg_y);
    ST* restrict vg_z=dpsi.data(2)+first_spo; ASSUME_ALIGNED(vg_z);
    ST* restrict gg_xx=grad_grad_psi.data(0)+first_spo; ASSUME_ALIGNED(gg_xx);
    ST* restrict gg_xy=grad_grad_psi.data(1)+first_spo; ASSUME_ALIGNED(gg_xy);
    ST* restrict gg_xz=grad_grad_psi.data(2)+first_spo; ASSUME_ALIGNED(gg_xz);
    ST* restrict gg_yx=grad_grad_psi.data(3)+first_spo; ASSUME_ALIGNED(gg_yx);
    ST* restrict gg_yy=grad_grad_psi.data(4)+first_spo; ASSUME_ALIGNED(gg_yy);
    ST* restrict gg_yz=grad_grad_psi.data(5)+first_spo; ASSUME_ALIGNED(gg_yz);
    ST* restrict gg_zx=grad_grad_psi.data(6)+first_spo; ASSUME_ALIGNED(gg_zx);
    ST* restrict gg_zy=grad_grad_psi.data(7)+first_spo; ASSUME_ALIGNED(gg_zy);
    ST* restrict gg_zz=grad_grad_psi.data(8)+first_spo; ASSUME_ALIGNED(gg_zz);

    const ST cone = (bc_sign &1)? -1:1;
#pragma simd 
    for (size_t j=0; j<N; ++j)
    {
      const ST kX=k0[j];
      const ST kY=k1[j];
      const ST kZ=k2[j];
      const ST val=myV[j];
      const ST kkV=mKK[j]*val;

      const ST gX = g00*g0[j]+g01*g1[j]+g02*g2[j];
      const ST gY = g10*g0[j]+g11*g1[j]+g12*g2[j];
      const ST gZ = g20*g0[j]+g21*g1[j]+g22*g2[j];

      psi[j]  =cone*val;
      vg_x[j] =cone*gX;
      vg_y[j] =cone*gY;
      vg_z[j] =cone*gZ;
      gg_xx[j]=cone*(h00[j] + kkV + kX*gX); 
      gg_xy[j]=cone*(h01[j] + kkV + kX*gY); 
      gg_xz[j]=cone*(h02[j] + kkV + kX*gZ);
      gg_yx[j]=cone*(h01[j] + kkV + kY*gX); 
      gg_yy[j]=cone*(h11[j] + kkV + kY*gY); 
      gg_yz[j]=cone*(h12[j] + kkV + kY*gZ);
      gg_zx[j]=cone*(h02[j] + kkV + kz*gX); 
      gg_zy[j]=cone*(h12[j] + kkV + kz*gY); 
      gg_zz[j]=cone*(h22[j] + kkV + kz*gZ);
    }
#endif
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgh(bc_sign,psi,dpsi,grad_grad_psi);
  }
};


}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: SplineR2RAdoptor.h 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
