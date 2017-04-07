//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SplineC2XAdoptor.h
 *
 * Adoptor classes to handle complex-to-(real,complex) with arbitrary precision
 */
#ifndef QMCPLUSPLUS_EINSPLINE_C2X_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_C2X_ADOPTOR_H

#include <spline/einspline_engine.hpp>

namespace qmcplusplus
{
//template<typename T, unsigned D>
//  inline Tensor<T,D> dot(const Tensor<T,D>& a, const Tensor<T,D>& b)
//  {
//    Tensor<T,D> res;
//    for(int i=0; i<D; ++i)
//      for(int j=0; j<D; ++j)
//      {
//        T s=0;
//        for(int k=0; k<D; ++k)
//            s+=a(i,k)*b(k,j);
//        res(i,j)=s;
//      }
//    return res;
//  }


/** adoptor class to match std::complex<ST> spline stored in a packed array with std::complex<TT> SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 */
template<typename ST, typename TT, unsigned D>
struct SplineC2CPackedAdoptor: public SplineAdoptorBase<ST,D>
{
  typedef typename einspline_traits<ST,D>::SplineType  SplineType;
  typedef typename einspline_traits<ST,D>::BCType      BCType;
  typedef typename SplineAdoptorBase<ST,D>::PointType  PointType;
  typedef typename SplineAdoptorBase<ST,D>::SingleSplineType SingleSplineType;

  using SplineAdoptorBase<ST,D>::first_spo;
  using SplineAdoptorBase<ST,D>::last_spo;
  using SplineAdoptorBase<ST,D>::GGt;
  using SplineAdoptorBase<ST,D>::PrimLattice;
  using SplineAdoptorBase<ST,D>::kPoints;

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

  //vector<ST> phase;

  SplineC2CPackedAdoptor():MultiSpline(0)
  {
    this->is_complex=true;
    this->AdoptorName="SplineC2CPackedAdoptor";
    this->KeyWord="C2CPacked";
  }

  inline void resizeStorage(int n, int nvals)
  {
    SplineAdoptorBase<ST,D>::init_base(n);
    myV.resize(2*n);
    myL.resize(2*n);
    myG.resize(2*n);
    myH.resize(2*n);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    MultiSpline=einspline::create(MultiSpline,xyz_g,xyz_bc,myV.size());
    for(int i=0; i<D; ++i)
    {
      BaseOffset[i]=0;
      BaseN[i]=xyz_g[i].num+3;
    }
    qmc_common.memory_allocated += MultiSpline->coefs_size*sizeof(ST);
  }

  inline void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    einspline::set(MultiSpline, 2*ispline, psi_r);
    einspline::set(MultiSpline, 2*ispline+1, psi_i);
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    einspline::set(MultiSpline, 2*ispline,  spline_r, BaseOffset, BaseN);
    einspline::set(MultiSpline, 2*ispline+1,spline_i, BaseOffset, BaseN);
  }

  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i, 
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
    einspline::set(MultiSpline, 2*ispline,  spline_r, offset_l, mesh_l);
    einspline::set(MultiSpline, 2*ispline+1,spline_i, offset_l, mesh_l);
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
    return h5f.write(bigtable,o.str().c_str()); //"spline_0");
  }

  /** assign myV to psi
   *
   * Taken out for the derived classes
   */
  template<typename VV>
  inline void assign_v(const PointType& r, int bc_sign, VV& psi)
  {
    register ST s,c;
    //TT* restrict t_ptr=reinterpret_cast<TT*>(psi.data())+(first_spo<<1);
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
    {
      int jr=j<<1;
      sincos(-dot(r,kPoints[j]),&s,&c);
      psi[psiIndex]= std::complex<TT>(c*myV[jr]-s*myV[jr+1],s*myV[jr]+c*myV[jr+1]);
      //t_ptr[jr  ]=c*myV[jr]-s*myV[jr+1];
      //t_ptr[jr+1]=s*myV[jr]+c*myV[jr+1];
    }
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.R[iat];
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate(MultiSpline,ru,myV);
    assign_v(r,0,psi);
    ////computePhases(r);
    ////simd::multiadd(psi.size(),eikr.data(),psi.data());
  }

  /** assign internal data to psi,dpsi,d2psi
   */
  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    const int N=kPoints.size();
    for (int j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (int j=0; j<2*N; j++)
      myL[j] = trace(myH[j],GGt);
    const ST two=2.0;
    ST s,c;
    PointType g_r, g_i;
    //can easily make three independent loops
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
    {
      int jr=j<<1;
      int ji=jr+1;
      g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      ST kk=-dot(kPoints[j],kPoints[j]);
      myL[jr]+=kk*myV[jr]+two*dot(kPoints[j],myG[ji]);
      myL[ji]+=kk*myV[ji]-two*dot(kPoints[j],myG[jr]);
      sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
      psi[psiIndex]= std::complex<TT>(c*myV[jr]-s*myV[ji],c*myV[ji]+s*myV[jr]);
      d2psi[psiIndex]= std::complex<TT>(c*myL[jr]-s*myL[ji],c*myL[ji]+s*myL[jr]);
      for(int idim=0; idim<D; ++idim)
        dpsi[psiIndex][idim]= std::complex<TT>(c*g_r[idim]-s*g_i[idim], c*g_i[idim]+s*g_r[idim]);
    }
    //complex<ST> e_mikr(c,s);
    //convert(e_mikr * myV[j], psi[j]);
    //convert(e_mikr*(-myV[j]*ck + myG[j]), dpsi[j]);
    //convert(e_mikr*(-myV[j]*kk - two*dot(ck,myG[j]) + myL[j]), d2psi[j]);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.R[iat];
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgl(r,0,psi,dpsi,d2psi);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    //convert to Cartesian G and Hess
    const int N=kPoints.size();
    for (int j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (int j=0; j<2*N; j++)
      myH[j] = dot(myH[j],GGt);
    ST s,c;
    PointType g_r, g_i;
    Tensor<ST,D> kk,h_r,h_i;
    //can easily make three independent loops
    for(int psiIndex=first_spo,j=0; psiIndex<last_spo; ++psiIndex,++j)
    {
      int jr=j<<1;
      int ji=jr+1;
      g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
      psi[psiIndex]= std::complex<TT>(c*myV[jr]-s*myV[ji],c*myV[ji]+s*myV[jr]);
      for(int idim=0; idim<D; ++idim)
        dpsi[psiIndex][idim]= std::complex<TT>(c*g_r[idim]-s*g_i[idim], c*g_i[idim]+s*g_r[idim]);
      kk=outerProduct(kPoints[j],kPoints[j]); // \f$kk=k^k \f$
      h_r=myH[jr]-myV[jr]*kk+outerProductSymm(kPoints[j],myG[ji]); //kdotg_i;
      h_i=myH[ji]-myV[ji]*kk-outerProductSymm(kPoints[j],myG[jr]); //kdotg_r;
      for(int t=0; t<D*D; ++t)
        grad_grad_psi[psiIndex](t)= std::complex<TT>(c*h_r(t)-s*h_i(t), c*h_i(t)+s*h_r(t));
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.R[iat];
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgh(r,0,psi,dpsi,grad_grad_psi);
  }

  template<typename VGL>
  void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& dpsi)
  {
    const PointType& r=P.R[iat];
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
struct SplineC2RPackedAdoptor: public SplineAdoptorBase<ST,D>
{
  typedef typename einspline_traits<ST,D>::SplineType  SplineType;
  typedef typename einspline_traits<ST,D>::BCType      BCType;
  typedef typename SplineAdoptorBase<ST,D>::PointType         PointType;
  typedef typename SplineAdoptorBase<ST,D>::SingleSplineType SingleSplineType;

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

  std::vector<ST>   KdotR;
  std::vector<ST>   CosV;
  std::vector<ST>   SinV;
  std::vector<ST>   mKK;
  std::vector<Tensor<ST,D> >  KK; //k^k

  SplineC2RPackedAdoptor():MultiSpline(0)
  {
    this->is_complex=true;
    this->AdoptorName="SplineC2RPackedAdoptor";
    this->KeyWord="C2RPacked";
  }

  virtual ~SplineC2RPackedAdoptor() {}

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
    mKK.resize(kPoints.size());
    for(int i=0; i<kPoints.size(); ++i)
      mKK[i]=-dot(kPoints[i],kPoints[i]);
    KK.resize(kPoints.size());
    for(int i=0; i<kPoints.size(); ++i)
      KK[i]=outerProduct(kPoints[i],kPoints[i]);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    //if(mKK.empty()) resize_kk();
    einspline::set(MultiSpline, 2*ispline, psi_r);
    einspline::set(MultiSpline, 2*ispline+1, psi_i);
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    //if(mKK.empty()) resize_kk();
    einspline::set(MultiSpline, 2*ispline,spline_r, BaseOffset, BaseN);
    einspline::set(MultiSpline, 2*ispline+1,spline_i, BaseOffset, BaseN);
  }

  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i, 
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
    if(mKK.empty()) resize_kk();
    einspline::set(MultiSpline, 2*ispline,  spline_r, offset_l, mesh_l);
    einspline::set(MultiSpline, 2*ispline+1,spline_i, offset_l, mesh_l);
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
    int psiIndex=first_spo;
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

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.R[iat];
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate(MultiSpline,ru,myV);
    assign_v(r,0,psi);
  }

  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, int bc_sign, VV& psi, GV& dpsi, VV& d2psi)
  {
    const int N=kPoints.size();
    for (int j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (int j=0; j<2*N; j++)
      myL[j] = trace(myH[j],GGt);
    const ST zero(0);
    const ST two(2);
    for(int j=0; j<N; ++j)
      KdotR[j]=-dot(r,kPoints[j]);
    eval_e2iphi(N,KdotR.data(),CosV.data(),SinV.data());
    //for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
    int psiIndex=first_spo;
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

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.R[iat];
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgl(r,0,psi,dpsi,d2psi);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, int bc_sign, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    //convert to Cartesian G and Hess
    const int N=kPoints.size();
    for (int j=0; j<2*N; j++)
      myG[j] = dot(PrimLattice.G, myG[j]);
    for (int j=0; j<2*N; j++)
      myH[j] = dot(myH[j],GGt);
    ST s,c;
    PointType g_r, g_i;
    Tensor<ST,D> h_r,h_i;
    //can easily make three independent loops
    int psiIndex=first_spo;
    for (int j=0,jr=0,ji=1; j<N; j++,jr+=2,ji+=2)
    {
      g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
      //kk=outerProduct(kPoints[j],kPoints[j]); // \f$kk=k^k \f$
      h_r=myH[jr]-myV[jr]*KK[j]+outerProductSymm(kPoints[j],myG[ji]);
      h_i=myH[ji]-myV[ji]*KK[j]-outerProductSymm(kPoints[j],myG[jr]);
      sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
      psi[psiIndex]=c*myV[jr]-s*myV[ji];
      for(int idim=0; idim<D; ++idim)
        dpsi[psiIndex][idim]=c*g_r[idim]-s*g_i[idim];
      for(int t=0; t<D*D; ++t)
        grad_grad_psi[psiIndex](t)=c*h_r(t)-s*h_i(t);
      ++psiIndex;
      if(MakeTwoCopies[j])
      {
        psi[psiIndex]=s*myV[jr]+c*myV[ji];
        for(int idim=0; idim<D; ++idim)
          dpsi[psiIndex][idim]=c*g_i[idim]+s*g_r[idim];
        for(int t=0; t<D*D; ++t)
          grad_grad_psi[psiIndex](t)=s*h_r(t)+c*h_i(t);
        ++psiIndex;
      }
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.R[iat];
    PointType ru(PrimLattice.toUnit_floor(r));
    einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
    assign_vgh(r,0,psi,dpsi,grad_grad_psi);
  }

  template<typename VGL>
  void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& dpsi)
  {
    const PointType& r=P.R[iat];
  }

};

//  /** adoptor class to match std::complex<ST> spline with std::complex<TT> SPOs, just references
//   * @tparam ST precision of spline
//   * @tparam TT precision of SPOs
//   * @tparam D dimension
//   *
//   * This is a reference implementation but NOT USED. Use Packed version
//   */
//  template<typename ST, typename TT, unsigned D>
//    struct SplineC2CAdoptor
//    {
//      typedef ST                                                       real_type;
//      typedef std::complex<ST>                                              value_type;
//      typedef typename einspline_traits<value_type,D>::SplineType      SplineType;
//      typedef typename einspline_traits<value_type,D>::BCType          BCType;
//      typedef typename OrbitalSetTraits<value_type>::ValueVector_t     StorageValueVector_t;
//      typedef typename OrbitalSetTraits<value_type>::GradVector_t      StorageGradVector_t;
//      typedef typename OrbitalSetTraits<value_type>::HessVector_t      StorageHessVector_t;
//      typedef typename OrbitalSetTraits<value_type>::GradHessVector_t  StorageGradHessVector_t;
//
//      typedef CrystalLattice<ST,D> UnitCellType;
//      typedef TinyVector<ST,D>     PointType;
//
//      SplineType          *MultiSpline;
//      UnitCellType        SuperLattice;
//      UnitCellType        PrimLattice;
//      std::vector<PointType>   kPoints;
//      TinyVector<int,D>   HalfG;
//      std::vector<bool>        MakeTwoCopies;
//      Tensor<real_type,D> GGt;
//
//      std::vector<real_type> phase;
//      std::vector<value_type> eikr;
//
//      StorageValueVector_t     myV;
//      StorageValueVector_t     myL;
//      StorageGradVector_t      myG;
//      StorageHessVector_t      myH;
//      StorageGradHessVector_t  myGH;
//
//      inline void resizeStorage(int n, int nvals)
//      {
//        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
//
//        kPoints.resize(n);
//        MakeTwoCopies.resize(n);
//        myV.resize(n);
//        myL.resize(n);
//        myG.resize(n);
//        myH.resize(n);
//        myH.resize(n);
//      }
//
//      inline bool isready()
//      {
//        return true;
//      }
//
//      template<typename VV>
//        inline void evaluate_v(const PointType& r, VV& psi)
//        {
//          PointType ru(PrimLattice.toUnit(r));
//          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
//          einspline::evaluate(MultiSpline,ru,myV);
//
//          const int N=myV.size();
//          //computePhases(r);
//          //simd::multiadd(psi.size(),eikr.data(),psi.data());
//          register ST s,c;
//          for(int j=0; j<psi.size(); ++j)
//          {
//            sincos(-dot(r,kPoints[j]),&s,&c);
//            psi[j]= std::complex<TT>(static_cast<TT>(c*myV[j].real()-s*myV[j].imag()),
//                static_cast<TT>(c*myV[j].imag()+s*myV[j].real()));
//          }
//
//        }
//
//      template<typename VV, typename GV>
//        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
//        {
//          PointType ru(PrimLattice.toUnit(r));
//          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
//          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
//
//          const int N=kPoints.size();
//          for (int j=0; j<N; j++) myG[j] = dot(PrimLattice.G, myG[j]);
//          for (int j=0; j<N; j++) myL[j] = trace(myH[j],GGt);
//
//          const ST two=2.0;
//          TinyVector<std::complex<ST>,D> ck;
//          register ST s,c;
//          for (int j=0; j<psi.size(); j++)
//          {
//            ST kk=dot(kPoints[j],kPoints[j]);
//            for(int i=0; i<D; ++i) ck[i]= std::complex<ST>(0.0,kPoints[j][i]);
//            sincos(-dot(r,kPoints[j]),&s,&c);
//            std::complex<ST> e_mikr(c,s);
//            convert(e_mikr * myV[j], psi[j]);
//            convert(e_mikr*(-myV[j]*ck + myG[j]), dpsi[j]);
//            convert(e_mikr*(-myV[j]*kk - two*dot(ck,myG[j]) + myL[j]), d2psi[j]);
//          }
//        }
//
//      template<typename VV, typename GV, typename GGV>
//        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
//        {
//        }
//    };
//
//  /** adoptor class to match std::complex<ST> spline with TT real SPOs
//   * @tparam ST precision of spline
//   * @tparam TT precision of SPOs
//   * @tparam D dimension
//   *
//   * Requires temporage storage and multiplication of phase vectors
//   */
//  template<typename ST, typename TT, unsigned D>
//    struct SplineC2RAdoptor
//    {
//      typedef ST                                                    real_type;
//      typedef std::complex<ST>                                           value_type;
//
//      typedef typename einspline_traits<value_type,D>::SplineType  SplineType;
//      typedef typename einspline_traits<value_type,D>::BCType      BCType;
//      typedef typename OrbitalSetTraits<value_type>::ValueVector_t          StorageValueVector_t;
//      typedef typename OrbitalSetTraits<value_type>::GradVector_t           StorageGradVector_t;
//      typedef typename OrbitalSetTraits<value_type>::HessVector_t           StorageHessVector_t;
//      typedef typename OrbitalSetTraits<value_type>::GradHessVector_t       StorageGradHessVector_t;
//
//      typedef CrystalLattice<ST,D> UnitCellType;
//      typedef TinyVector<ST,D> PointType;
//
//      SplineType*         MultiSpline;
//      UnitCellType        SuperLattice;
//      UnitCellType        PrimLattice;
//      TinyVector<int,D>   HalfG;
//      Tensor<real_type,D> GGt;
//      std::vector<PointType>   kPoints;
//      std::vector<bool>        MakeTwoCopies;
//      std::vector<real_type>   CosV;
//      std::vector<real_type>   SinV;
//      std::vector<value_type>  mKK;
//
//      // Temporary storage for Eispline calls
//      StorageValueVector_t     myV;
//      StorageValueVector_t     myL;
//      StorageGradVector_t      myG;
//      StorageHessVector_t      myH;
//      StorageGradHessVector_t  myGH;
//
//      SplineC2RAdoptor():MultiSpline(0) { }
//
//      virtual ~SplineC2RAdoptor(){}
//
//      inline void resizeStorage(int n, int nvals)
//      {
//
//        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
//
//        HalfG=0;
//        kPoints.resize(n);
//        MakeTwoCopies.resize(n);
//
//        myV.resize(n);
//        myL.resize(n);
//        myG.resize(n);
//        myH.resize(n);
//        CosV.resize(n);
//        SinV.resize(n);
//      }
//
//      inline bool isready()
//      {
//        mKK.resize(kPoints.size());
//        for(int i=0; i<kPoints.size(); ++i) mKK[i]=-dot(kPoints[i],kPoints[i]);
//        return true;
//      }
//
//      template<typename VV>
//        inline void evaluate_v(const PointType& r, VV& psi)
//        {
//          PointType ru(PrimLattice.toUnit(r));
//          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
//          einspline::evaluate(MultiSpline,ru,myV);
//          const int N=kPoints.size();
//          for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
//
//          int psiIndex = 0;
//          for (int j=0; j<N; j++)
//          {
//            psi[psiIndex] = static_cast<TT>(myV[j].real()*CosV[j]-myV[j].imag()*SinV[j]);
//            //psi[psiIndex] = static_cast<TT>(myV[j].real());
//            psiIndex++;
//            if (MakeTwoCopies[j])
//            {
//              //psi[psiIndex] = static_cast<TT>(myV[j].imag());
//              psi[psiIndex] = static_cast<TT>(myV[j].imag()*CosV[j]+myV[j].real()*SinV[j]);
//              psiIndex++;
//            }
//          }
//        }
//
//      template<typename VV, typename GV>
//      inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
//        {
//          PointType ru(PrimLattice.toUnit(r));
//          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
//
//          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
//          const int N=kPoints.size();
//          for (int j=0; j<N; j++) myG[j] = dot(PrimLattice.G, myG[j]);
//          for (int j=0; j<N; j++) myL[j] = trace(myH[j],GGt);
//
//          const ST zero=0.0;
//          const ST two=2.0;
//          TinyVector<std::complex<ST>,D> ck;
//          ST s,c;
//          for(int j=0; j<N; ++j)
//          {
//            for (int n=0; n<D; n++)ck[n] = std::complex<ST>(zero,-kPoints[j][n]);
//            sincos (-dot(r,kPoints[j]), &s, &c);
//            std::complex<ST> e_mikr(c,s);
//            myV[j]  = e_mikr*myV[j];
//            myL[j]  = -dot(kPoints[j],kPoints[j])*myV[j]+e_mikr*(two*dot(ck,myG[j]) + myL[j]);
//            myG[j]  = myV[j]*ck+e_mikr*myG[j];
//          }
//
//          //const std::complex<ST> eye (0.0, 1.0);
//          //ST s,c;
//          //for(int j=0; j<N; ++j)
//          //{
//          //  std::complex<ST> u = myV[j];
//          //  TinyVector<std::complex<ST>,D> gradu = myG[j];
//          //  std::complex<ST> laplu = myL[j];
//          //  PointType k = kPoints[j];
//          //  TinyVector<std::complex<ST>,D> ck;
//          //  for (int n=0; n<D; n++)	  ck[n] = k[n];
//          //  sincos (-dot(r,k), &s, &c);
//          //  std::complex<ST> e_mikr (c,s);
//          //  myV[j]  = e_mikr*u;
//          //  myG[j]  = e_mikr*(-eye*u*ck + gradu);
//          //  myL[j]  = e_mikr*(-dot(k,k)*u - two*eye*dot(ck,gradu) + laplu);
//          //}
//
//          int psiIndex=0;
//          for(int j=0; j<N; ++j)
//          {
//            psi[psiIndex]=static_cast<TT>(myV[j].real());
//            for(int n=0; n<D; ++n) dpsi[psiIndex][n]=static_cast<TT>(myG[j][n].real());
//            d2psi[psiIndex]=static_cast<TT>(myL[j].real());
//            ++psiIndex;
//            if(MakeTwoCopies[j])
//            {
//              psi[psiIndex]=static_cast<TT>(myV[j].imag());
//              for(int n=0; n<D; ++n) dpsi[psiIndex][n]=static_cast<TT>(myG[j][n].imag());
//              d2psi[psiIndex]=static_cast<TT>(myL[j].imag());
//              ++psiIndex;
//            }
//          }
//        }
//
//      template<typename VV, typename GV, typename GGV>
//        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
//        {
//        }
//    };
//
}
#endif
