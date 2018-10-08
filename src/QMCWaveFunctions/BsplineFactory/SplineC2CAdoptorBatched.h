//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

// This seems like a waste of effort
///\file SplineC2CAdoptorCUDA.h
#ifndef QMCPLUSPLUS_BATCHED_C2C_SOA_ADOPTOR_H
#define QMCPLUSPLUS_BATCHED_C2C_SOA_ADOPTOR_H

#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>
#include "QMCWaveFunctions/BsplineFactory/SplineC2CAdoptor.h"
#include "spline/einspline_util.hpp"

//#define USE_VECTOR_ML 1

namespace qmcplusplus
{

/** adoptor class to match std::complex<ST> spline with TT Complex SPOs
 * when using a cuda device
 *
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 *
 * Requires temporage storage and multiplication of phase vectors
 */
template<typename ST, typename TT>
struct SplineC2CAdoptorBatched: public SplineC2CAdoptor<ST, TT>,
				public SplineAdoptorBatched<BsplineDeviceCUDA, ST, OHMMS_DIM>
{
  static const int D= OHMMS_DIM;
  using BaseType=SplineC2CAdoptor<ST,TT>;
  using SplineType=typename bspline_traits<ST,3>::SplineType;
  using BCType=typename bspline_traits<ST,3>::BCType;
  using DataType=ST;
  using PointType=typename BaseType::PointType;
  using SingleSplineType=typename BaseType::SingleSplineType;

  using vContainer_type=Vector<ST,aligned_allocator<ST> >;
  using gContainer_type=VectorSoaContainer<ST,3>;
  using hContainer_type=VectorSoaContainer<ST,6>;

  using BaseType::first_spo;
  using BaseType::last_spo;
  using BaseType::GGt;
  using BaseType::PrimLattice;
  using BaseType::kPoints;
  using BaseType::MakeTwoCopies;
  using BaseType::offset_cplx;
  using BaseType::offset_real;


  SplineC2CAdoptorBatched(): BaseType()
  {
    this->is_complex=true;
    this->is_soa_ready=true;
    this->AdoptorName="SplineC2CAdoptorBatched";
    this->KeyWord="SplineC2CAdoptorBatched";
  }

  ///** copy the base property */
  //SplineC2CSoA(BaseType& rhs): BaseType(rhs)
  //{
  //  this->is_complex=true;
  //  this->AdoptorName="SplineC2CSoA";
  //  this->KeyWord="C2RSoA";
  //}

  // SplineC2CAdoptorBatched(const SplineC2CAdoptorBatched& a):
  //   SplineAdoptorBatched<ST,3>(a),SplineInst(a.SplineInst),MultiSpline(nullptr),
  //   mKK(a.mKK), myKcart(a.myKcart)
  // {
  //   const size_t n=a.myL.size();
  //   myV.resize(n); myG.resize(n); myL.resize(n); myH.resize(n);
  // }


  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseType::resizeStorage(n, nvals);
  }

  void bcast_tables(Communicate* comm)
  {
    chunked_bcast(comm, BaseType::MultiSpline);
  }

  void gather_tables(Communicate* comm)
  {
    BaseType::gather_tables(comm);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    BaseType::create_spline(xyz_g, xyz_bc);
  }

  inline void flush_zero()
  {
    BaseType::flush_zero();
  }

  /** remap kPoints to pack the double copy */
  inline void resize_kpoints()
  {
    BaseType::resize_kpoints();
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    BaseType::set_spline(spline_r, spline_i, twist, ispline, level);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    Vector<ST> v_r(psi_r,0), v_i(psi_i,0);
    BaseType::SplineInst->set(2*ispline  ,v_r);
    BaseType::SplineInst->set(2*ispline+1,v_i);
  }


  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i,
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
  }

  // bool read_splines(hdf_archive& h5f)
  // {
  //   std::ostringstream o;
  //   o<<"spline_" << SplineAdoptor<ST,D>::MyIndex;
  //   einspline_engine<SplineType> bigtable(SplineInst->spline_m);
  //   return h5f.read(bigtable,o.str().c_str());//"spline_0");
  // }

  // bool write_splines(hdf_archive& h5f)
  // {
  //   std::ostringstream o;
  //   o<<"spline_" << SplineAdoptor<ST,D>::MyIndex;
  //   einspline_engine<SplineType> bigtable(BaseType::SplineInst->spline_m);
  //   return h5f.write(bigtable,o.str().c_str());//"spline_0");
  // }

  inline std::complex<TT> 
    evaluate_dot(const ParticleSet& P, const int iat, const std::complex<TT>* restrict arow, ST* scratch,
        bool compute_spline=true)
  {
    return BaseType::evaluate_dot(P, iat, arow, scratch,
			compute_spline);
  }

  template<typename VV>
  inline void assign_v(const PointType& r, const vContainer_type& myV, VV& psi)
  {
    BaseType::assign_v(r, myV, psi);
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    BaseType::evaluate_v(P, iat, psi);
  }

  template<typename VM>
  inline void evaluateValues(const VirtualParticleSet& VP, VM& psiM)
  {
    BaseType::evaluateValues(VP, psiM);
  }

  inline size_t estimateMemory(const int nP) { return 0; }

  /** assign_vgl
   */
  template<typename VV, typename GV>
  inline void assign_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    BaseType::assign_vgl(r, psi, dpsi, d2psi);
  }

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  template<typename VV, typename GV>
  inline void assign_vgl_from_l(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    BaseType::assign_vgl_from_l(r, psi, dpsi, d2psi);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    BaseType::evaluate_vgl(P, iat, psi, dpsi, d2psi);
  }

  /** identical to assign_vgl but the output container is SoA container
   */
  template<typename VGL>
  inline void assign_vgl_soa(const PointType& r, VGL& vgl)
  {
    BaseType::assign_vgl_soa(r, vgl);
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& vgl)
  {
    BaseType::evaluate_vgl_combo(P, iat, vgl);
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    BaseType::assign_vgh(r, psi, dpsi, grad_grad_psi);
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    BaseType::evaluate_vgh(P, iat, psi, dpsi, grad_grad_psi);
  }
};

}
#endif
