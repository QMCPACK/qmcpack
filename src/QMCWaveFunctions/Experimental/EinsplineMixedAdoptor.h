//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file EsinplineMixedAdoptor.h
 *
 * SplineAdoptor to save memory for the vacuum for 3D systems
 *
 * A template parameter for BsplineSet<SplineAdoptor>
 * - EinsplineMixedAdoptor: for slabs or wires
 * - EinsplineOpenAdoptor : for molecular systems with an orthorhombic cell
 */
#ifndef QMCPLUSPLUS_EINSPLINE_MIXED_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_MIXED_ADOPTOR_H

#include <spline/einspline_util.hpp>
namespace qmcplusplus
{

/** adoptor class to match ST real spline with TT real SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 */
template<typename ST, typename TT, unsigned D>
struct SplineMixedAdoptor
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
  Tensor<ST,D> GGt;
  UnitCellType PrimLattice;
  GridConvert<ST> gTransform;
  //these are not needed for R2R adoptor and will be removed
  UnitCellType SuperLattice;

  ///used for testing only, to be removed
  bool UseFullGrid;

  // Temporary storage for Eispline calls
  StorageValueVector_t myV, myL;
  StorageGradVector_t      myG;
  StorageHessVector_t      myH;
  StorageGradHessVector_t  myGH;

  SplineMixedAdoptor(): MultiSpline(0), smallBox(0)
  {
  }

  void resizeStorage(int n, int nv)
  {
    myV.resize(n);
    myL.resize(n);
    myG.resize(n);
    myH.resize(n);
    myGH.resize(n);
  }

  void create_spline(TinyVector<int,D>& mesh, int n, bool samegrid)
  {
    UseFullGrid=samegrid;
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    Ugrid xyz_grid[D];
    BCType xyz_bc[D];
    for(int i=0; i<D; ++i)
    {
      xyz_grid[i].start = 0.0;
      xyz_grid[i].end = 1.0;
      xyz_grid[i].num = (samegrid)?mesh[i]:mesh[i]/2;
      xyz_bc[i].lCode=xyz_bc[i].rCode=PERIODIC;
      gTransform.BaseOffset[i]=0;
      gTransform.BaseN[i]=xyz_grid[i].num+3;
    }
    SplineType* dummy=0;
    MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);
  }

  template <typename UBspline, typename PT>
  void add_box(UBspline* dense, PT& lower, PT& upper)
  {
    if(smallBox==0)
      gTransform.create(smallBox,dense,lower,upper,MultiSpline->num_splines);
    Lower=lower+2.0*gTransform.Delta;
    Upper=upper-2.0*gTransform.Delta;
  }

  template <typename UBspline>
  void init_spline(UBspline* dense, UBspline* coarse, int ival)
  {
    if(UseFullGrid)
      einspline::set(MultiSpline,ival,dense, gTransform.BaseOffset,  gTransform.BaseN);
    else
      einspline::set(MultiSpline,ival,coarse, gTransform.BaseOffset,  gTransform.BaseN);
    einspline::set(smallBox,   ival,dense,  gTransform.Offset,gTransform.N);
  }

  inline bool isready()
  {
    return true;
  }

  ///** return sign */
  inline void convertPos(const PointType& r, PointType& ru, int& sign)
  {
    ru=PrimLattice.toUnit(r);
    for (int i=0; i<D; i++)
      ru[i]-=std::floor(ru[i]);
    sign=0; //need to add ANTIPERIODIC
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    int phase;
    TinyVector<ST,D> ru;
    convertPos(r,ru,phase);
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
    int phase;
    TinyVector<ST,D> ru;
    convertPos(r,ru,phase);
    if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
      einspline::evaluate_vgh(smallBox,ru,myV,myG,myH);
    else
      einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    const int N=psi.size();
    const Tensor<ST,D> gConv(PrimLattice.G);
    for(int j=0; j<N; ++j)
      psi[j]=myV[j];
    for(int j=0; j<N; ++j)
      dpsi[j]=dot(gConv,myG[j]);
    for(int j=0; j<N; ++j)
      d2psi[j]=trace(myH[j],GGt);
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {}
};

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

}
#endif
