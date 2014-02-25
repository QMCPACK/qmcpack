/////////////////////////////////////////////////////////////////nspli
//jjjkkkkkk/
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
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

/** specialization for 3D complex<double> */
template<>
struct einspline_traits<complex<double>,3>
{
  typedef multi_UBspline_3d_z SplineType;
  typedef UBspline_3d_z       SingleSplineType;
  typedef BCtype_z            BCType;
  typedef complex<double>     DataType;
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

/** specialization for 3D complex<float> */
template<>
struct einspline_traits<complex<float>,3>
{
  typedef multi_UBspline_3d_c SplineType;
  typedef UBspline_3d_c       SingleSplineType;
  typedef BCtype_c            BCType;
  typedef complex<float>      DataType;
};

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
  typedef TinyVector<ST,D>   PointType;
  typedef UBspline_3d_d      SingleSplineType;
  ///true, if this for complex, each derived class has to set this
  bool is_complex;
  ///Index of this adoptor, when multiple adoptors are used for NUMA or distributed cases
  int MyIndex;
  ///number of unique orbitals
  int nunique_orbitals;
  ///first index of the SPOs this Spline handles
  int first_spo;
  ///last index of the SPOs this Spline handles
  int last_spo;
  ///sign bits at the G/2 boundaries
  TinyVector<int,D>          HalfG;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST,D>               GGt;
  CrystalLattice<ST,D>       SuperLattice;
  CrystalLattice<ST,D>       PrimLattice;
  /// flags to unpack sin/cos
  vector<bool>               MakeTwoCopies;
  /// kpoints for each unique orbitals
  vector<TinyVector<ST,D> >  kPoints;

  ///name of the adoptor
  string AdoptorName;
  ///keyword used to match hdf5
  string KeyWord;

  typedef typename einspline_traits<ST,D>::DataType   DataType;
  typename OrbitalSetTraits<ST>::ValueVector_t     myV;
  typename OrbitalSetTraits<ST>::ValueVector_t     myL;
  typename OrbitalSetTraits<ST>::GradVector_t      myG;
  typename OrbitalSetTraits<ST>::HessVector_t      myH;
  typename OrbitalSetTraits<ST>::GradHessVector_t  myGH;

  SplineAdoptorBase()
    :is_complex(false),MyIndex(0),nunique_orbitals(0),first_spo(0),last_spo(0)
  {
  }

  inline void init_base(int n)
  {
    nunique_orbitals=n;
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
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

  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    SplineAdoptor::evaluate_v(P.R[iat],psi);
  }

  inline void evaluateValues(const ParticleSet& P, ValueMatrix_t& psiM)
  {
    ValueVector_t psi(psiM.cols());
    for(int iat=0; iat<P.getTotalNum(); ++iat)
    {
      SplineAdoptor::evaluate_v(P.R[iat],psi);
      std::copy(psi.begin(),psi.end(),psiM[iat]);
    }
  }

  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    SplineAdoptor::evaluate_vgl(P.R[iat],psi,dpsi,d2psi);
  }
  inline void evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
  {
    SplineAdoptor::evaluate_vgh(P.R[iat],psi,dpsi,grad_grad_psi);
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
      SplineAdoptor::evaluate_vgl(P.R[iat],v,g,l);
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
      SplineAdoptor::evaluate_vgh(P.R[iat],v,g,h);
    }
  }

};

}
#endif
