//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file EsinplineAdoptor.h
 *
 * Adoptor classes and BsplineSet<SplineAdoptor>
 * BsplineSet<SplineAdoptor> is a SPOSetBase class to work with determinant classes
 * SplineAdoptor provides these functions
 * - evaluate_v    value only
 * - evaluate_vgl  vgl
 * - evaluate_vgh  vgh
 * Specializations are implemented  in Spline*Adoptor.h and include
 * - SplineC2RAdoptor<ST,TT,D> : real wavefunction using complex einspline, tiling
 * - SplineC2CAdoptor<ST,TT,D> : complex wavefunction using complex einspline, tiling
 * - SplineR2RAdoptor<ST,TT,D> : real wavefunction using real einspline, a single twist
 * where ST (TT) is the precision of the einspline (SPOSetBase).
 *
 * typedefs and data members are duplicated for each adoptor class.
 * @todo Specalization and optimization for orthorhombic cells to use vgl not vgh
 */
#ifndef QMCPLUSPLUS_EINSPLINE_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_ADOPTOR_H

#include <Lattice/CrystalLattice.h>
#include <spline/einspline_engine.hpp>
#include <spline/einspline_util.hpp>
#include <simd/allocator.hpp>
#include <Numerics/VectorViewer.h>

namespace qmcplusplus
{

/** einspline trait class equivalent to  MultiOrbitalTraits
 * @tparam ST spline datatype
 * @tparam D dimension
 * @tparam TT target datatype
 */
template<typename ST, unsigned D>
struct einspline_traits { };

/** specialization for 3D double */
template<>
struct einspline_traits<double,3>
{
  typedef multi_UBspline_3d_d SplineType;
  typedef UBspline_3d_d       SingleSplineType;
  typedef BCtype_d            BCType;
  typedef double              DataType;
};

/** specialization for 3D float */
template<>
struct einspline_traits<float,3>
{
  typedef multi_UBspline_3d_s SplineType;
  typedef UBspline_3d_s       SingleSplineType;
  typedef BCtype_s            BCType;
  typedef float               DataType;
};

#if 0
/** specialization for 3D std::complex<double> */
template<>
struct einspline_traits<std::complex<double>,3>
{
  typedef multi_UBspline_3d_z SplineType;
  typedef UBspline_3d_z       SingleSplineType;
  typedef BCtype_z            BCType;
  typedef std::complex<double>     DataType;
};

/** specialization for 3D std::complex<float> */
template<>
struct einspline_traits<std::complex<float>,3>
{
  typedef multi_UBspline_3d_c SplineType;
  typedef UBspline_3d_c       SingleSplineType;
  typedef BCtype_c            BCType;
  typedef std::complex<float>      DataType;
};
#endif

/** symmetric outer product
 * @param v a vector
 * @param w a vector
 * @return \f$v^w+w^v\f$
 *
 * Used for the gradient coming from the phase and gradient
 */
template<typename T, unsigned D>
inline Tensor<T,D> outerProductSymm(const TinyVector<T,D>& v, const TinyVector<T,D>& w)
{
  Tensor<T,D> res;
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
      res(i,j)=v(i)*w(j)+v(j)*w(i);
  return res;
}

//inline void computePhases(const PointType& r)
//{
//  for (int i=0; i<kPoints.size(); i++) phase[i] = -dot(r, kPoints[i]);
//  eval_e2iphi(kPoints.size(),phase.data(),eikr.data());
//}

/** base class any SplineAdoptor
 *
 * This handles SC and twist and declare storage for einspline
 */
template<typename ST, unsigned D>
struct SplineAdoptorBase
{
#if (__cplusplus >= 201103L)
  using PointType=TinyVector<ST,D>;
  using SingleSplineType=UBspline_3d_d;
  using DataType=ST; 
#else
  typedef TinyVector<ST,D> PointType;
  typedef UBspline_3d_d SingleSplineType;
  typedef ST DataType;
#endif
  ///true if the computed values are complex
  bool is_complex;
  ///true, if it has only one k point and Gamma
  bool is_gamma_only;
  ///true, if it has only one k point and Gamma
  bool is_soa_ready;
  ///Index of this adoptor, when multiple adoptors are used for NUMA or distributed cases
  size_t MyIndex;
  ///number of unique orbitals
  size_t nunique_orbitals;
  ///first index of the SPOs this Spline handles
  size_t first_spo;
  ///last index of the SPOs this Spline handles
  size_t last_spo;
  ///sign bits at the G/2 boundaries
  TinyVector<int,D>          HalfG;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST,D>               GGt;
  CrystalLattice<ST,D>       SuperLattice;
  CrystalLattice<ST,D>       PrimLattice;
  /// flags to unpack sin/cos
  std::vector<bool>               MakeTwoCopies;
  ///kpoints for each unique orbitals
  std::vector<TinyVector<ST,D> >  kPoints;
  ///remap band
  aligned_vector<int> BandIndexMap;

  ///name of the adoptor
  std::string AdoptorName;
  ///keyword used to match hdf5
  std::string KeyWord;

  SplineAdoptorBase()
    :is_complex(false),is_gamma_only(false), is_soa_ready(false),
    MyIndex(0),nunique_orbitals(0),first_spo(0),last_spo(0)
  { }

#if (__cplusplus >= 201103L)
  SplineAdoptorBase(const SplineAdoptorBase& rhs)=default;
#endif

  inline void init_base(int n)
  {
    nunique_orbitals=n;
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
  }

  ///remap kpoints to group general kpoints & special kpoints
  int remap_kpoints()
  {
    std::vector<TinyVector<ST,D> >  k_copy(kPoints);
    const int nk=kPoints.size();
    BandIndexMap.resize(nk);
    int nCB=0;
    //two pass
    for(int i=0; i<nk; ++i)
    {
      if(MakeTwoCopies[i]) 
      {
        kPoints[nCB]=k_copy[i];
        BandIndexMap[i]=nCB++;
      }
    }
    int nRealBands=nCB;
    for(int i=0; i<nk; ++i)
    {
      if(!MakeTwoCopies[i]) 
      {
        kPoints[nRealBands]=k_copy[i];
        BandIndexMap[i]=nRealBands++;
      }
    }
    return nCB; //return the number of complex bands
  }
};


/** BsplineSet<SplineAdoptor>, a SPOSetBase
 * @tparam SplineAdoptor implements evaluation functions that matched the storage requirements.
 *
 * Equivalent to EinsplineSetExtended<Storage>
 * Storage is now handled by SplineAdoptor class that is specialized for precision, storage etc.
 * @todo Make SplineAdoptor be a member not the base class. This is needed
 * to make MultiBsplineSet (TBD) which has multiple SplineAdoptors for distributed
 * cases.
 */
template<typename SplineAdoptor>
struct BsplineSet: public SPOSetBase, public SplineAdoptor
{
  typedef typename SplineAdoptor::SplineType SplineType;
  typedef typename SplineAdoptor::PointType  PointType;
  typedef typename SplineAdoptor::DataType  DataType;

  ///** default constructor */
  //BsplineSet() { }

  SPOSetBase* makeClone() const
  {
    return new BsplineSet<SplineAdoptor>(*this);
  }

  /** set_spline to the big table
   * @param psi_r starting address of real part of psi(ispline)
   * @param psi_i starting address of imaginary part of psi(ispline)
   * @param twist twist id, reserved to sorted adoptor, ignored
   * @param ispline index of this spline function
   * @param level refinement level
   *
   * Each adoptor handles the map with respect to the twist, state index and refinement level
   */
  template<typename CT>
  void set_spline(CT* spline_r, CT* spline_i, int twist, int ispline, int level)
  {
    SplineAdoptor::set_spline(spline_r,spline_i,twist,ispline,level);
  }

  inline ValueType RATIO(const ParticleSet& P, int iat, const ValueType* restrict arow)
  {
    //this is just an example how to resuse t_logpsi
    int ip=omp_get_thread_num()*2;
    // YYYY: need to fix
    //return SplineAdoptor::evaluate_dot(P,iat,arow,reinterpret_cast<DataType*>(t_logpsi[ip]));
    return ValueType();
  }

  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    SplineAdoptor::evaluate_v(P,iat,psi);
  }

  inline void evaluateValues(const ParticleSet& P, ValueMatrix_t& psiM)
  {
    const size_t m=psiM.cols();
    for(int iat=0; iat<P.getTotalNum(); ++iat)
    {
      ValueVector_t psi(psiM[iat],m);
      SplineAdoptor::evaluate_v(P,iat,psi);
    }
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    SplineAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);

#if 0
    //debug GL combo
    CONSTEXPR double eps=std::numeric_limits<float>::epsilon();
    ValueVector_t psi_copy(psi);
    GLVector_t gl(psi.size());
    SplineAdoptor::evaluate_vgl_combo(P,iat,psi_copy,gl);
    auto gradX=gl.data(0);
    auto gradY=gl.data(1);
    auto gradZ=gl.data(2);
    auto lap=gl.data(3);
    double v_err=0, g_err=0, l_err=0;
    for(size_t i=0; i<psi.size(); ++i)
    {
      v_err+=std::abs(psi[i]-psi_copy[i]);
      double dx=std::abs(dpsi[i][0]-gradX[i]);
      double dy=std::abs(dpsi[i][1]-gradY[i]);
      double dz=std::abs(dpsi[i][2]-gradZ[i]);
      g_err+=std::sqrt(dx*dx+dy*dy+dz*dz);
      l_err+=std::abs(d2psi[i]-lap[i]);
    }
    if(v_err>eps || g_err > eps || l_err>eps)
      std::cout << "ERROR " << v_err << " " << g_err << " " << l_err << std::endl;
#endif
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
  {
    SplineAdoptor::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
  }

  void resetParameters(const opt_variables_type& active)
  { }

  void resetTargetParticleSet(ParticleSet& e)
  { }

  void setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
    BasisSetSize=norbs;
    //SplineAdoptor::first_spo=0;
    //SplineAdoptor::last_spo=norbs;
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      VectorViewer<value_type> v(logdet[i],OrbitalSetSize);
      VectorViewer<grad_type> g(dlogdet[i],OrbitalSetSize);
      VectorViewer<value_type> l(d2logdet[i],OrbitalSetSize);
      SplineAdoptor::evaluate_vgl(P,iat,v,g,l);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    typedef HessMatrix_t::value_type hess_type;
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      VectorViewer<value_type> v(logdet[i],OrbitalSetSize);
      VectorViewer<grad_type> g(dlogdet[i],OrbitalSetSize);
      VectorViewer<hess_type> h(grad_grad_logdet[i],OrbitalSetSize);
      SplineAdoptor::evaluate_vgh(P,iat,v,g,h);
    }
  }

  /** einspline does not need any other state data */
  void evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl, bool newp)
  {
    SplineAdoptor::evaluate_vgl_combo(P,iat,vgl);
  }

};

}
#endif
