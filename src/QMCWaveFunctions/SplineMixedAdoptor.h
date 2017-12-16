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
    
    
/** @file EsinplineMixedAdoptor.h
 *
 * SplineAdoptor to save memory for the vacuum for 3D systems
 *
 * SplineAdoptor for mixed grid. Supports full-grid orbitals and truncated orbitals.
 * - EinsplineMixedAdoptor: for slabs or wires
 * - EinsplineOpenAdoptor : for molecular systems with an orthorhombic cell
 */
#ifndef QMCPLUSPLUS_EINSPLINE_MIXED_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_MIXED_ADOPTOR_H

namespace qmcplusplus
{

template<typename SMA>
struct hdf_dual_grid
{
  static bool read(SMA* bspline, hdf_archive& h5f)
  {
    TinyVector<double,3> lower_in;
    TinyVector<double,3> upper_in;
    h5f.read(lower_in,"lower_bound");
    h5f.read(upper_in,"upper_bound");
    bool foundspline=true;
    lower_in-=bspline->Lower;
    upper_in-=bspline->Upper;
    if(dot(lower_in,lower_in)<1e-12 &&dot(upper_in,upper_in)<1e-12)
    {
      einspline_engine<typename SMA::SplineType> bigtable(bspline->MultiSpline);
      einspline_engine<typename SMA::SplineType> smalltable(bspline->smallBox);
      foundspline=h5f.read(bigtable,"spline_0");
      foundspline=h5f.read(smalltable,"spline_1");
    }
    else
    {
      app_log() << "  The upper/lower bound of the input is different from the current value."<< std::endl;
      foundspline=false;
    }
    return foundspline;
  }

  static bool write(SMA* bspline, hdf_archive& h5f)
  {
    einspline_engine<typename SMA::SplineType> bigtable(bspline->MultiSpline);
    einspline_engine<typename SMA::SplineType> smalltable(bspline->smallBox);
    TinyVector<double,3> lower(bspline->Lower);
    TinyVector<double,3> upper(bspline->Upper);
    int doneit=1;
    doneit=h5f.write(lower,"lower_bound");
    doneit=h5f.write(upper,"upper_bound");
    doneit=h5f.write(bigtable,"spline_0");
    doneit=h5f.write(smalltable,"spline_1");
    return doneit;
  }

};
/** adoptor class to match ST real spline with TT real SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 */
template<typename ST, typename TT, unsigned D>
struct SplineMixedAdoptor: public SplineR2RAdoptor<ST,TT,D>
{
  typedef typename einspline_traits<ST,D>::SplineType SplineType;
  typedef typename einspline_traits<ST,D>::BCType     BCType;
  typedef typename einspline_traits<ST,D>::DataType   DataType;
  typedef typename SplineAdoptorBase<ST,D>::PointType PointType;
  typedef typename SplineAdoptorBase<ST,D>::SingleSplineType SingleSplineType;
  typedef SplineMixedAdoptor<ST,TT,D> ThisType;

  using SplineAdoptorBase<ST,D>::HalfG;
  using SplineAdoptorBase<ST,D>::GGt;
  using SplineAdoptorBase<ST,D>::PrimLattice;

  typename OrbitalSetTraits<ST>::ValueVector_t     myV;
  typename OrbitalSetTraits<ST>::ValueVector_t     myL;
  typename OrbitalSetTraits<ST>::GradVector_t      myG;
  typename OrbitalSetTraits<ST>::HessVector_t      myH;
  typename OrbitalSetTraits<ST>::GradHessVector_t  myGH;

  SplineType *MultiSpline;
  SplineType *smallBox;
  TinyVector<ST,D> Lower;
  TinyVector<ST,D> Upper;
  GridConvert<ST> gTransform;

  ///used for testing only, to be removed
  bool UseFullGrid;

  // Temporary storage for Eispline calls
  SplineMixedAdoptor(): MultiSpline(0), smallBox(0)
  {
    this->AdoptorName="SplineR2RMixedAdoptor";
    this->AdoptorName="R2RMixed";
  }

  void create_spline(TinyVector<int,D>& mesh_0, TinyVector<int,D>& mesh_1
                     , int n, int nval)
  {
    UseFullGrid=(mesh_0[0]==mesh_1[0]);
    for(int i=1; i<D; ++i)
      UseFullGrid &= (mesh_0[i]==mesh_1[i]);
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    Ugrid xyz_grid[D];
    BCType xyz_bc[D];
    for(int i=0; i<D; ++i)
    {
      xyz_grid[i].start = 0.0;
      xyz_grid[i].end = 1.0;
      xyz_grid[i].num = mesh_1[i]; //(samegrid)?mesh[i]:mesh[i]/2;
      xyz_bc[i].lCode=xyz_bc[i].rCode=(HalfG[i])?ANTIPERIODIC:PERIODIC;
      gTransform.BaseOffset[i]=0;
      gTransform.BaseN[i]=xyz_grid[i].num+3;
    }
    SplineType* dummy=0;
    if(this->is_complex)
      MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,2*n);
    else
      MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);
  }

  bool read_splines(hdf_archive& h5f)
  {
    return hdf_dual_grid<ThisType>::read(this,h5f);
  }

  bool write_splines(hdf_archive& h5f)
  {
    return hdf_dual_grid<ThisType>::write(this,h5f);
  }

  template <typename PT>
  void add_box(SingleSplineType* dense, PT& lower, PT& upper)
  {
    if(smallBox==0)
    {
      gTransform.create(smallBox,dense,lower,upper,MultiSpline->num_splines);
    }
    Lower=lower;
    Upper=upper;
  }

  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    if(level==0)
      einspline::set(MultiSpline,ispline,spline_r, gTransform.BaseOffset,gTransform.BaseN);
    else
      einspline::set(smallBox,   ispline,spline_r, gTransform.Offset, gTransform.N);
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.R[iat];
    PointType ru;
    int bc_sign=this->convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate(smallBox,ru,myV);
    else
      einspline::evaluate(MultiSpline,ru,myV);
    this->assign_v(r,bc_sign,psi);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.R[iat];
    PointType ru;
    int bc_sign=this->convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate_vgh(smallBox,ru,myV,myG,myH);
    else
      einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    this->assign_vgl(r,bc_sign,psi,dpsi,d2psi);
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.R[iat];
    PointType ru;
    int bc_sign=this->convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate_vgh(smallBox,ru,myV,myG,myH);
    else
      einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    this->assign_vgh(r,bc_sign,psi,dpsi,grad_grad_psi);
  }

  template<typename VV, typename GL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VV& psi, GL& dpsi)
  {
    const PointType& r=P.R[iat];
  }
};

/** adoptor class for the non-periodic systems
 *
 * Support two-level grid system for big systems
 * - MultiSpline : full cell with the coarse grid
 * - smallBox    : small cell with the original grid
 * No coordinate transform is needed for the open BC with a orthorhombic cell
 */
template<typename ST, typename TT, unsigned D>
struct SplineOpenAdoptor: public SplineAdoptorBase<ST,D>
{
  typedef typename einspline_traits<ST,D>::SplineType SplineType;
  typedef typename einspline_traits<ST,D>::BCType     BCType;
  typedef typename einspline_traits<ST,D>::DataType   DataType;
  typedef typename SplineAdoptorBase<ST,D>::PointType PointType;
  typedef typename SplineAdoptorBase<ST,D>::SingleSplineType SingleSplineType;
  typedef SplineOpenAdoptor<ST,TT,D> ThisType;

  using SplineAdoptorBase<ST,D>::first_spo;
  using SplineAdoptorBase<ST,D>::last_spo;
  using SplineAdoptorBase<ST,D>::SuperLattice;

  typename OrbitalSetTraits<ST>::ValueVector_t     myV;
  typename OrbitalSetTraits<ST>::ValueVector_t     myL;
  typename OrbitalSetTraits<ST>::GradVector_t      myG;
  typename OrbitalSetTraits<ST>::HessVector_t      myH;
  typename OrbitalSetTraits<ST>::GradHessVector_t  myGH;

  SplineType *MultiSpline;
  SplineType *smallBox;

  TinyVector<ST,D> Lower;
  TinyVector<ST,D> Upper;
  TinyVector<ST,D> L;
  TinyVector<ST,D> InvL;

  ///used for testing only
  bool UseFullGrid;
  ///grid transform
  GridConvert<ST> gTransform;

  SplineOpenAdoptor(): MultiSpline(0), smallBox(0)
  {
    this->is_complex=false;
    this->AdoptorName="SplineOpenAdoptor";
    this->KeyWord="Open";
  }

  void resizeStorage(int n, int nv)
  {
    SplineAdoptorBase<ST,D>::init_base(n);
    myV.resize(n);
    myL.resize(n);
    myG.resize(n);
    myH.resize(n);
  }

  /** create MultiSpline for the full cell with a coarse grid
   * @param mesh_0 original mesh (e.g. FFT)
   * @param mesh_1 coarse mesh (e.g. FFT)
   * @param n number of distict einsplines
   * @param nval number of valence els
   */
  void create_spline(TinyVector<int,D>& mesh_0, TinyVector<int,D>& mesh_1
                     , int n, int nval)
  {
    UseFullGrid=(mesh_0[0]==mesh_1[0]);
    for(int i=1; i<D; ++i)
      UseFullGrid &= (mesh_0[i]==mesh_1[i]);
    Ugrid xyz_grid[D];
    BCType xyz_bc[D];
    for(int i=0; i<D; ++i)
    {
      L[i]=SuperLattice.R(i,i);
      InvL[i]=1.0/L[i];
      xyz_grid[i].start = 0.0;
      xyz_grid[i].end = L[i]; //1.0;
      xyz_grid[i].num = mesh_1[i];
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
  template <typename PT>
  void add_box(SingleSplineType* dense, PT& lower, PT& upper)
  {
    if(smallBox==0)
      gTransform.create(smallBox,dense,lower,upper,MultiSpline->num_splines);
    Lower=lower;
    Upper=upper;
  }

  /** set ival-th state
   * @param dense a single-spline on the full grid
   * @param coarse a single-spline on the half grid
   * @param ival the state to be initialized
   */
  //template <typename UBspline>
  //  void set_spline(UBspline* dense, UBspline* coarse, int ival)
  //  {
  //    if(UseFullGrid)
  //      einspline::set(MultiSpline,ival,dense,  gTransform.BaseOffset,gTransform.BaseN);
  //    else
  //      einspline::set(MultiSpline,ival,coarse, gTransform.BaseOffset,gTransform.BaseN);

  //    einspline::set(smallBox,ival,dense,gTransform.Offset,gTransform.N);
  //  }
  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    if(level==0)
      einspline::set(MultiSpline,ispline,spline_r, gTransform.BaseOffset,gTransform.BaseN);
    else
      einspline::set(smallBox,ispline,spline_r,gTransform.Offset,gTransform.N);
  }

  bool read_splines(hdf_archive& h5f)
  {
    return hdf_dual_grid<ThisType>::read(this,h5f);
  }

  bool write_splines(hdf_archive& h5f)
  {
    return hdf_dual_grid<ThisType>::write(this,h5f);
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
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.R[iat];
    TinyVector<ST,D> ru;
    convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate(smallBox,ru,myV);
    else
      einspline::evaluate(MultiSpline,ru,myV);
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      psi[psiIndex]=static_cast<TT>(myV[j]);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.R[iat];
    TinyVector<ST,D> ru;
    convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate_vgl(smallBox,ru,myV,myG,myL);
    else
      einspline::evaluate_vgl(MultiSpline,ru,myV,myG,myL);
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      psi[psiIndex]=myV[j];
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      dpsi[psiIndex]=myG[j];
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
      d2psi[psiIndex]=myL[j];
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.R[iat];
    TinyVector<ST,D> ru;
    convertPos(r,ru);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate_vgh(smallBox,ru,myV,myG,myH);
    else
      einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    const int N=psi.size();
    for(int j=0; j<N; ++j)
      psi[j]=myV[j];
    for(int j=0; j<N; ++j)
      dpsi[j]=myG[j];
    for(int j=0; j<N; ++j)
      grad_grad_psi[j]=myH[j];
  }

  template<typename VV, typename GL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VV& psi, GL& dpsi)
  {
    const PointType& r=P.R[iat];
  }
};

}
#endif
