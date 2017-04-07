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
    
    
#ifndef QMCPLUSPLUS_EINSPLINE_R2RADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_R2RADOPTOR_H

#include <Numerics/VectorViewer.h>
#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>

namespace qmcplusplus
{

/** adoptor class to match ST real spline with TT real SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 */
template<typename ST, typename TT, unsigned D>
struct SplineR2RAdoptor: public SplineAdoptorBase<ST,D>
{
  typedef typename einspline_traits<ST,D>::SplineType SplineType;
  typedef typename einspline_traits<ST,D>::BCType     BCType;
  typedef typename SplineAdoptorBase<ST,D>::PointType PointType;
  typedef typename SplineAdoptorBase<ST,D>::SingleSplineType SingleSplineType;

  using SplineAdoptorBase<ST,D>::first_spo;
  using SplineAdoptorBase<ST,D>::last_spo;
  using SplineAdoptorBase<ST,D>::HalfG;
  using SplineAdoptorBase<ST,D>::GGt;
  using SplineAdoptorBase<ST,D>::PrimLattice;

  typename OrbitalSetTraits<ST>::ValueVector_t     myV;
  typename OrbitalSetTraits<ST>::ValueVector_t     myL;
  typename OrbitalSetTraits<ST>::GradVector_t      myG;
  typename OrbitalSetTraits<ST>::HessVector_t      myH;
  typename OrbitalSetTraits<ST>::GradHessVector_t  myGH;

  SplineType *MultiSpline;

  ///number of points of the original grid
  int BaseN[3];
  ///offset of the original grid, always 0
  int BaseOffset[3];

  SplineR2RAdoptor(): MultiSpline(0)
  {
    this->is_complex=false;
    this->AdoptorName="SplineR2RAdoptor";
    this->KeyWord="R2R";
  }

  void resizeStorage(int n, int nv)
  {
    SplineAdoptorBase<ST,D>::init_base(n);
    myV.resize(n);
    myL.resize(n);
    myG.resize(n);
    myH.resize(n);
    myGH.resize(n);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    MultiSpline=einspline::create(MultiSpline,xyz_g,xyz_bc,myV.size());

    qmc_common.memory_allocated += MultiSpline->coefs_size*sizeof(ST);

    for(int i=0; i<D; ++i)
    {
      BaseOffset[i]=0;
      BaseN[i]=xyz_g[i].num+3;
    }
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
    SplineType* dummy=0;
    MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);

    qmc_common.memory_allocated += MultiSpline->coefs_size*sizeof(ST);
  }

  inline bool isready()
  {
    return true;
  }

  inline void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    einspline::set(MultiSpline, ispline,psi_r);
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    einspline::set(MultiSpline, ispline,spline_r, BaseOffset,BaseN);
  }

  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i, 
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
    einspline::set(MultiSpline, ispline,  spline_r, offset_l, mesh_l);
  }

  bool read_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o<<"spline_" << SplineAdoptorBase<ST,D>::MyIndex;
    einspline_engine<SplineType> bigtable(MultiSpline);
    return h5f.read(bigtable,o.str().c_str()); //"spline_0");
  }

  bool write_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o<<"spline_" << SplineAdoptorBase<ST,D>::MyIndex;
    einspline_engine<SplineType> bigtable(MultiSpline);
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

  /** assign myV to psi
   */
  template<typename VV>
  inline void assign_v(const PointType& r, int bc_sign, VV& psi)
  {
    if (bc_sign & 1)
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=static_cast<TT>(-myV[j]);
    else
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=static_cast<TT>(myV[j]);
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.R[iat];
    PointType ru;
    int bc_sign=convertPos(r,ru);
    einspline::evaluate(MultiSpline,ru,myV);
    assign_v(r,bc_sign,psi);
  }

  /** assign internal data to psi's
   */
  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    if (bc_sign & 1)
    {
      const ST minus_one=-1.0;
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=-myV[j];
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        dpsi[psiIndex]=minus_one*dot(PrimLattice.G,myG[j]);
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        d2psi[psiIndex]=-trace(myH[j],GGt);
    }
    else
    {
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=myV[j];
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        dpsi[psiIndex]=dot(PrimLattice.G,myG[j]);
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        d2psi[psiIndex]=trace(myH[j],GGt);
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.R[iat];
    PointType ru;
    int bc_sign=convertPos(r,ru);
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgl(r,bc_sign,psi,dpsi,d2psi);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    if (bc_sign & 1)
    {
      const ST minus_one=-1.0;
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=-myV[j];
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        dpsi[psiIndex]=minus_one*dot(PrimLattice.G,myG[j]);
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        grad_grad_psi[psiIndex]=minus_one*dot(myH[j],GGt);
    }
    else
    {
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        psi[psiIndex]=myV[j];
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        dpsi[psiIndex]=dot(PrimLattice.G,myG[j]);
      for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
        grad_grad_psi[psiIndex]=dot(myH[j],GGt);
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.R[iat];
    PointType ru;
    int bc_sign=convertPos(r,ru);
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgh(r,bc_sign,psi,dpsi,grad_grad_psi);
  }

  template<typename VGL>
  void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& dpsi)
  {
    const PointType& r=P.R[iat];
  }
};


}
#endif

