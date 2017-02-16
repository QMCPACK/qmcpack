//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory 
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory

//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file HybridCplxAdoptor.h
 *
 * Adoptor classes to handle complex hybrid orbitals with arbitrary precision
 */
#ifndef QMCPLUSPLUS_HYBRID_CPLX_SOA_ADOPTOR_H
#define QMCPLUSPLUS_HYBRID_CPLX_SOA_ADOPTOR_H

#include <QMCWaveFunctions/BsplineFactory/HybridAdoptorBase.h>
namespace qmcplusplus
{

/** adoptor class to match 
 *
 */
template<typename BaseAdoptor>
struct HybridCplxSoA: public BaseAdoptor, public HybridAdoptorBase<typename BaseAdoptor::DataType>
{
  using ST               = typename BaseAdoptor::DataType;
  using PointType        = typename BaseAdoptor::PointType;
  using SingleSplineType = typename BaseAdoptor::SingleSplineType;

  using BaseAdoptor::PrimLattice;
  using BaseAdoptor::myV;
  using BaseAdoptor::myG;
  using BaseAdoptor::myL;
  using BaseAdoptor::myH;
  using BaseAdoptor::SplineInst;

  HybridCplxSoA(): BaseAdoptor()
  {
    this->is_complex=true;
    this->is_soa_ready=true;
    this->AdoptorName="Hybrid"+this->AdoptorName;
    this->KeyWord="Hybrid"+this->KeyWord;
  }

  HybridCplxSoA(const HybridCplxSoA& a):
    BaseAdoptor(a)
  {
  }

  ~HybridCplxSoA() 
  { 
    //if(MultiSpline != nullptr) delete SplineInst;
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseAdoptor::resizeStorage(n, nvals); 
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    BaseAdoptor::create_spline(xyz_g, xyz_bc);
  }

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    BaseAdoptor::set_spline(spline_r, spline_i, twist, ispline, level);
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
    BaseAdoptor::set_spline(psi_r, psi_i, twist, ispline, level);
  }


  inline void set_spline_domain(SingleSplineType* spline_r, SingleSplineType* spline_i, 
      int twist, int ispline, const int* offset_l, const int* mesh_l)
  {
    BaseAdoptor::set_spline_domain(spline_r, spline_i, twist, ispline, offset_l, mesh_l);
  }

  bool read_splines(hdf_archive& h5f)
  {
    BaseAdoptor::read_splines(h5f);
  }

  bool write_splines(hdf_archive& h5f)
  {
    BaseAdoptor::write_splines(h5f);
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate(ru,myV);
    BaseAdoptor::assign_v(r,psi);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    BaseAdoptor::assign_vgl(r,psi,dpsi,d2psi);
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const PointType& r, VGL& vgl)
  {
    PointType ru(PrimLattice.toUnit_floor(r));
    SplineInst->evaluate_vgh(ru,myV,myG,myH);
    BaseAdoptor::assign_vgl_soa(r,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    //missing
  }
};

}
#endif
