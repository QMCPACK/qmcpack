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
    
    
/** @file EsinplineAdoptor2.h
 *
 * Adoptor classes for BsplineSet<SplineAdoptor>  using packed storage for einspline
 */
#ifndef QMCPLUSPLUS_EINSPLINE_ADOPTOR_PACKED_H
#define QMCPLUSPLUS_EINSPLINE_ADOPTOR_PACKED_H
//#include <spline/einspline_engine.hpp>
namespace qmcplusplus
{

/** adoptor class to match std::complex<ST> spline stored in a packed array with std::complex<TT> SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 */
template<typename ST, typename TT, unsigned D>
struct SplineC2CAdoptorPacked
{
  typedef ST                                                      real_type;
  typedef std::complex<ST>                                             value_type;
  typedef typename einspline_traits<real_type,D>::SplineType      SplineType;
  typedef typename einspline_traits<real_type,D>::BCType          BCType;
  typedef typename OrbitalSetTraits<real_type>::ValueVector_t     StorageValueVector_t;
  typedef typename OrbitalSetTraits<real_type>::GradVector_t      StorageGradVector_t;
  typedef typename OrbitalSetTraits<real_type>::HessVector_t      StorageHessVector_t;
  typedef typename OrbitalSetTraits<real_type>::GradHessVector_t  StorageGradHessVector_t;

  typedef CrystalLattice<ST,D> UnitCellType;
  typedef TinyVector<ST,D>     PointType;

  SplineType          *MultiSpline;
  UnitCellType        SuperLattice;
  UnitCellType        PrimLattice;
  TinyVector<int,D>   HalfG;
  std::vector<bool>        MakeTwoCopies;
  Tensor<real_type,D> GGt;
  std::vector<PointType>   kPoints;

  std::vector<real_type> phase;
  std::vector<value_type> eikr;

  StorageValueVector_t     myV;
  StorageValueVector_t     myL;
  StorageGradVector_t      myG;
  StorageHessVector_t      myH;
  StorageGradHessVector_t  myGH;

  inline void resizeStorage(int n, int nvals)
  {
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
    myV.resize(2*n);
    myL.resize(2*n);
    myG.resize(2*n);
    myH.resize(2*n);
  }

  inline bool isready()
  {
    return true;
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    PointType ru(PrimLattice.toUnit(r));
    for (int i=0; i<D; i++)
      ru[i] -= std::floor (ru[i]);
    einspline::evaluate(MultiSpline,ru,myV);
    //computePhases(r);
    //simd::multiadd(psi.size(),eikr.data(),psi.data());
    register ST s,c;
    TT* restrict t_ptr=reinterpret_cast<TT*>(psi.data());
    for(int j=0,jr=0; j<psi.size(); ++j, jr+=2)
    {
      sincos(-dot(r,kPoints[j]),&s,&c);
      t_ptr[jr]  =c*myV[jr]-s*myV[jr+1];
      t_ptr[jr+1]=s*myV[jr]+c*myV[jr+1];
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    PointType ru(PrimLattice.toUnit(r));
    for (int i=0; i<D; i++)
      ru[i] -= std::floor (ru[i]);
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    const int N=kPoints.size();
    for (int j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (int j=0; j<2*N; j++)
      myL[j] = trace(myH[j],GGt);
    const ST two=2.0;
    ST s,c;
    TinyVector<ST,D> g_r, g_i;
    //can easily make three independent loops
    for (int j=0,jr=0,ji=1; j<N; j++,jr+=2,ji+=2)
    {
      g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      ST kk=-dot(kPoints[j],kPoints[j]);
      myL[jr]+=kk*myV[jr]+two*dot(kPoints[j],myG[ji]);
      myL[ji]+=kk*myV[ji]-two*dot(kPoints[j],myG[jr]);
      sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
      psi[j]= std::complex<TT>(c*myV[jr]-s*myV[ji],c*myV[ji]+s*myV[jr]);
      d2psi[j]= std::complex<TT>(c*myL[jr]-s*myL[ji],c*myL[ji]+s*myL[jr]);
      for(int idim=0; idim<D; ++idim)
        dpsi[j][idim]= std::complex<TT>(c*g_r[idim]-s*g_i[idim], c*g_i[idim]+s*g_r[idim]);
      //complex<ST> e_mikr(c,s);
      //convert(e_mikr * myV[j], psi[j]);
      //convert(e_mikr*(-myV[j]*ck + myG[j]), dpsi[j]);
      //convert(e_mikr*(-myV[j]*kk - two*dot(ck,myG[j]) + myL[j]), d2psi[j]);
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
  }
};

/** adoptor class to match std::complex<ST> spline with TT real SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 *
 * Requires temporage storage and multiplication of phase vectors
 */
template<typename ST, typename TT, unsigned D>
struct SplineC2RAdoptorPacked
{
  typedef ST                                                    real_type;
  typedef std::complex<ST>                                           value_type;

  typedef typename einspline_traits<real_type,D>::SplineType           SplineType;
  typedef typename einspline_traits<real_type,D>::BCType               BCType;
  typedef typename OrbitalSetTraits<real_type>::ValueVector_t          StorageValueVector_t;
  typedef typename OrbitalSetTraits<real_type>::GradVector_t           StorageGradVector_t;
  typedef typename OrbitalSetTraits<real_type>::HessVector_t           StorageHessVector_t;
  typedef typename OrbitalSetTraits<real_type>::GradHessVector_t       StorageGradHessVector_t;

  typedef CrystalLattice<ST,D> UnitCellType;
  typedef TinyVector<ST,D> PointType;

  SplineType          *MultiSpline;
  UnitCellType        SuperLattice;
  UnitCellType        PrimLattice;
  TinyVector<int,D>   HalfG;
  Tensor<real_type,D> GGt;
  std::vector<PointType>   kPoints;
  std::vector<bool>        MakeTwoCopies;
  std::vector<real_type>   CosV;
  std::vector<real_type>   SinV;
  std::vector<real_type>   mKK;

  // Temporary storage for Eispline calls
  StorageValueVector_t     myV;
  StorageValueVector_t     myL;
  StorageGradVector_t      myG;
  StorageHessVector_t      myH;
  StorageGradHessVector_t  myGH;

  SplineC2RAdoptorPacked():MultiSpline(0) { }

  virtual ~SplineC2RAdoptorPacked() {}

  inline void resizeStorage(int n, int nvals)
  {
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
    myV.resize(2*n);
    myL.resize(2*n);
    myG.resize(2*n);
    myH.resize(2*n);
    CosV.resize(n);
    SinV.resize(n);
  }

  inline bool isready()
  {
    mKK.resize(kPoints.size());
    for(int i=0; i<kPoints.size(); ++i)
      mKK[i]=-dot(kPoints[i],kPoints[i]);
    return true;
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    PointType ru(PrimLattice.toUnit(r));
    for (int i=0; i<D; i++)
      ru[i] -= std::floor (ru[i]);
    einspline::evaluate(MultiSpline,ru,myV);
    int N=kPoints.size();
    ST phase[N];
    for(int j=0; j<N; ++j)
      phase[j]=-dot(r,kPoints[j]);
    eval_e2iphi(N,phase,CosV.data(),SinV.data());
    //for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
    int psiIndex = 0;
    for (int j=0,jr=0; j<N; j++,jr+=2)
    {
      psi[psiIndex] = static_cast<TT>(myV[jr]*CosV[j]-myV[jr+1]*SinV[j]);
      psiIndex++;
      if (MakeTwoCopies[j])
      {
        psi[psiIndex] = static_cast<TT>(myV[jr+1]*CosV[j]+myV[jr]*SinV[j]);
        psiIndex++;
      }
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    PointType ru(PrimLattice.toUnit(r));
    for (int i=0; i<D; i++)
      ru[i] -= std::floor (ru[i]);
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    const int N=kPoints.size();
    for (int j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (int j=0; j<2*N; j++)
      myL[j] = trace(myH[j],GGt);
    const ST zero=0.0;
    const ST two=2.0;
    for(int j=0; j<N; ++j)
      sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
    int psiIndex=0;
    TinyVector<ST,D> g_r, g_i;
    for (int j=0,jr=0,ji=1; j<N; j++,jr+=2,ji+=2)
    {
      g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      myL[jr]+=mKK[j]*myV[jr]+two*dot(kPoints[j],myG[ji]);
      myL[ji]+=mKK[j]*myV[ji]-two*dot(kPoints[j],myG[jr]);
      psi[psiIndex]=CosV[j]*myV[jr]-SinV[j]*myV[ji];
      dpsi[psiIndex]=CosV[j]*g_r-SinV[j]*g_i; // multiply phase
      d2psi[psiIndex]=CosV[j]*myL[jr]-SinV[j]*myL[ji];
      ++psiIndex;
      if(MakeTwoCopies[j])
      {
        psi[psiIndex]=CosV[j]*myV[ji]+SinV[j]*myV[jr];
        dpsi[psiIndex]=CosV[j]*g_i+SinV[j]*g_r;
        d2psi[psiIndex]=CosV[j]*myL[ji]+SinV[j]*myL[jr];
        ++psiIndex;
      }
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
  }
};

}
#endif
