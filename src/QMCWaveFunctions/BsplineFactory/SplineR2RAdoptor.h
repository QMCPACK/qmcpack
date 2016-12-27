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
    this->AdoptorName="SplineR2RSoA";
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
    SplineInst->create(xyz_g,xyz_bc,n);
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
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      {
        psi[psiIndex]=-myV[j];
        dpsi[psiIndex][0]=-(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        dpsi[psiIndex][1]=-(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        dpsi[psiIndex][2]=-(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        d2psi[psiIndex]=-SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],GGt.data());
      }
    }
    else
    {
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      {
        psi[psiIndex]=myV[j];
        dpsi[psiIndex][0]=(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        dpsi[psiIndex][1]=(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        dpsi[psiIndex][2]=(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        d2psi[psiIndex]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],GGt.data());
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
  template<typename VV, typename GL>
  inline void assign_vgl_soa(int bc_sign, VV& psi, GL& dpsi)
  {
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

    TT* restrict vg_x=dpsi.data(0); ASSUME_ALIGNED(vg_x);
    TT* restrict vg_y=dpsi.data(1); ASSUME_ALIGNED(vg_y);
    TT* restrict vg_z=dpsi.data(2); ASSUME_ALIGNED(vg_z);
    TT* restrict vl_l=dpsi.data(3); ASSUME_ALIGNED(vl_l);

    if (bc_sign & 1)
    {
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      {
        psi[psiIndex]=-myV[j];
        vg_x[psiIndex]=-(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        vg_y[psiIndex]=-(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        vg_z[psiIndex]=-(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        vl_l[psiIndex]=-SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],GGt.data());
      }
    }
    else
    {
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      {
        psi[psiIndex]=myV[j];
        vg_x[psiIndex]=(g00*g0[j]+g01*g1[j]+g02*g2[j]);
        vg_y[psiIndex]=(g10*g0[j]+g11*g1[j]+g12*g2[j]);
        vg_z[psiIndex]=(g20*g0[j]+g21*g1[j]+g22*g2[j]);
        vl_l[psiIndex]=SymTrace(h00[j],h01[j],h02[j],h11[j],h12[j],h22[j],GGt.data());
      }
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
    PointType ru;
    int bc_sign=convertPos(r,ru);
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    assign_vgl_soa(bc_sign,psi,dpsi);
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

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5260 $   $Date: 2011-06-18 07:45:58 -0400 (Sat, 18 Jun 2011) $
 * $Id: SplineR2RAdoptor.h 5260 2011-06-18 11:45:58Z jeongnim.kim $
 ***************************************************************************/
