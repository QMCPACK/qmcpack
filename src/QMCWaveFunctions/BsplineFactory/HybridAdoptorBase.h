//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

/** @file HybridAdoptorBase.h
 *
 * Hybrid adoptor base class
 */
#ifndef QMCPLUSPLUS_HYBRID_ADOPTOR_BASE_H
#define QMCPLUSPLUS_HYBRID_ADOPTOR_BASE_H

#include <Particle/DistanceTableData.h>
#include <QMCWaveFunctions/lcao/SoaSphericalTensor.h>
#include <spline2/MultiBspline.hpp>

namespace qmcplusplus
{

template<typename ST>
struct AtomicOrbitalSoA
{
  static const int D=3;
  using AtomicSplineType=typename bspline_traits<ST,1>::SplineType;
  using AtomicBCType=typename bspline_traits<ST,1>::BCType;
  using AtomicSingleSplineType=UBspline_1d_d;
  using PointType=TinyVector<ST,D>;
  using value_type=ST;

  using vContainer_type=aligned_vector<ST>;

  // near core cutoff
  ST rmin;
  // far from core cutoff, rmin_sqrt>=rmin
  ST rmin_sqrt;
  ST cutoff, spline_radius;
  int spline_npoints, BaseN;
  int NumBands, Npad;
  PointType pos;
  const int lmax, lm_tot;
  SoaSphericalTensor<ST> Ylm;
  vContainer_type l_vals;
  vContainer_type r_power_minus_l;
  AtomicSplineType* MultiSpline;
  MultiBspline1D<ST>* SplineInst;

  vContainer_type localV, localG, localL;

  AtomicOrbitalSoA(int Lmax):
  Ylm(Lmax), MultiSpline(nullptr), SplineInst(nullptr), lmax(Lmax),
  lm_tot((Lmax+1)*(Lmax+1))
  {
    r_power_minus_l.resize(lm_tot);
    l_vals.resize(lm_tot);
    for(int l=0; l<=lmax; l++)
      for(int m=-l; m<=l; m++)
        l_vals[l*(l+1)+m] = l;
    rmin = std::exp(std::log(std::numeric_limits<ST>::min())/std::max(Lmax,1));
    rmin = std::max(rmin,std::numeric_limits<ST>::epsilon());
    rmin_sqrt=std::max(rmin,std::sqrt(std::numeric_limits<ST>::epsilon()));
  }

  inline void resizeStorage(size_t Nb)
  {
    NumBands=Nb;
    Npad=getAlignedSize<ST>(Nb);
    localV.resize(Npad*lm_tot);
    localG.resize(Npad*lm_tot);
    localL.resize(Npad*lm_tot);
    create_spline();
    qmc_common.memory_allocated += SplineInst->sizeInByte();
  }

  void bcast_tables(Communicate* comm)
  {
    chunked_bcast(comm, MultiSpline);
  }

  void reduce_tables(Communicate* comm)
  {
    chunked_reduce(comm, MultiSpline);
  }

  template<typename PT, typename VT>
  inline void set_info(const PT& R, const VT& cutoff_in, const VT& spline_radius_in, const int& spline_npoints_in)
  {
    pos[0]=R[0];
    pos[1]=R[1];
    pos[2]=R[2];
    cutoff=cutoff_in;
    spline_radius=spline_radius_in;
    spline_npoints=spline_npoints_in;
    BaseN=spline_npoints+2;
  }

  inline void create_spline()
  {
    AtomicBCType bc;
    bc.lCode = FLAT;
    bc.rCode = NATURAL;
    Ugrid grid;
    grid.start = 0.0;
    grid.end   = spline_radius;
    grid.num   = spline_npoints;
    SplineInst = new MultiBspline1D<ST>();
    SplineInst->create(grid, bc, lm_tot*Npad);
    MultiSpline=&(SplineInst->spline_m);
  }

  inline void flush_zero()
  {
    SplineInst->flush_zero();
  }

  inline void set_spline(AtomicSingleSplineType* spline, int lm, int ispline)
  {
    SplineInst->copy_spline(spline, lm*Npad+ispline, 0, BaseN);
  }

  bool read_splines(hdf_archive& h5f)
  {
    einspline_engine<AtomicSplineType> bigtable(MultiSpline);
    int lmax_in, spline_npoints_in;
    ST spline_radius_in;
    bool success=true;
    success = success && h5f.read(lmax_in, "l_max");
    success = success && h5f.read(spline_radius_in, "spline_radius");
    success = success && h5f.read(spline_npoints_in, "spline_npoints");
    if(lmax_in!=lmax) return false;
    if(spline_radius_in!=spline_radius) return false;
    if(spline_npoints_in!=spline_npoints) return false;
    return success && h5f.read(bigtable,"radial_spline");
  }

  bool write_splines(hdf_archive& h5f)
  {
    bool success=true;
    success = success && h5f.write(spline_radius, "spline_radius");
    success = success && h5f.write(spline_npoints, "spline_npoints");
    success = success && h5f.write(lmax, "l_max");
    success = success && h5f.write(pos, "position");
    einspline_engine<AtomicSplineType> bigtable(MultiSpline);
    success = success && h5f.write(bigtable,"radial_spline");
    return success;
  }

  //evaluate only V
  template<typename VV>
  inline void evaluate_v(const ST& r, const PointType& dr, VV& myV)
  {
    if (r>std::numeric_limits<ST>::epsilon())
      Ylm.evaluateV(dr[0]/r, dr[1]/r, dr[2]/r);
    else
      Ylm.evaluateV(0,0,1);
    const ST* restrict Ylm_v=Ylm[0];

    CONSTEXPR ST czero(0);
    ST* restrict val=myV.data();
    ST* restrict local_val=localV.data();
    std::fill(myV.begin(),myV.end(),czero);

    SplineInst->evaluate(r,localV);

    for(size_t lm=0; lm<lm_tot; lm++)
    {
      #pragma omp simd aligned(val,local_val)
      for(size_t ib=0; ib<myV.size(); ib++)
        val[ib]+=Ylm_v[lm]*local_val[ib];
      local_val+=Npad;
    }
  }

  //evaluate VGL
  template<typename VV, typename GV>
  inline void evaluate_vgl(const ST& r, const PointType& dr, VV& myV, GV& myG, VV& myL)
  {
    ST drx, dry, drz, rhatx, rhaty, rhatz, rinv;
    if (r>rmin)
    {
      rinv=1.0/r;
      drx=dr[0];
      dry=dr[1];
      drz=dr[2];
      rhatx=drx*rinv;
      rhaty=dry*rinv;
      rhatz=drz*rinv;
    }
    else
    {
      rinv=0;
      drx=dr[0];
      dry=dr[1];
      drz=dr[2];
    }

    Ylm.evaluateVGL(drx, dry, drz);
    const ST* restrict Ylm_v=Ylm[0];
    const ST* restrict Ylm_gx=Ylm[1];
    const ST* restrict Ylm_gy=Ylm[2];
    const ST* restrict Ylm_gz=Ylm[3];

    ST* restrict g0=myG.data(0);
    ST* restrict g1=myG.data(1);
    ST* restrict g2=myG.data(2);
    CONSTEXPR ST czero(0), cone(1), chalf(0.5);
    std::fill(myV.begin(),myV.end(),czero);
    std::fill(g0,g0+Npad,czero);
    std::fill(g1,g1+Npad,czero);
    std::fill(g2,g2+Npad,czero);
    std::fill(myL.begin(),myL.end(),czero);
    ST* restrict val=myV.data();
    ST* restrict lapl=myL.data();
    ST* restrict local_val=localV.data();
    ST* restrict local_grad=localG.data();
    ST* restrict local_lapl=localL.data();

    SplineInst->evaluate_vgl(r,localV,localG,localL);

    if(r>rmin_sqrt)
    {
      // far from core
      r_power_minus_l[0]=cone;
      ST r_power_temp=cone;
      for(int l=1; l<=lmax; l++)
      {
        r_power_temp*=rinv;
        for(int m=-l, lm=l*l; m<=l; m++,lm++)
          r_power_minus_l[lm]=r_power_temp;
      }

      for(size_t lm=0; lm<lm_tot; lm++)
      {
        const ST& l_val=l_vals[lm];
        const ST& r_power=r_power_minus_l[lm];
        const ST Ylm_rescale=Ylm_v[lm]*r_power;
        #pragma omp simd aligned(val,g0,g1,g2,lapl,local_val,local_grad,local_lapl)
        for(size_t ib=0; ib<myV.size(); ib++)
        {
          const ST local_v=local_val[ib];
          const ST local_g=local_grad[ib];
          const ST local_l=local_lapl[ib];
          // value
          const ST Vpart = l_val*rinv*local_v;
          val[ib] += Ylm_rescale*local_v;

          // grad
          g0[ib] += local_g * rhatx * Ylm_rescale + local_v * Ylm_gx[lm] * r_power - Vpart * rhatx * Ylm_rescale;
          g1[ib] += local_g * rhaty * Ylm_rescale + local_v * Ylm_gy[lm] * r_power - Vpart * rhaty * Ylm_rescale;
          g2[ib] += local_g * rhatz * Ylm_rescale + local_v * Ylm_gz[lm] * r_power - Vpart * rhatz * Ylm_rescale;

          // laplacian
          ST rhat_dot_G = ( rhatx*Ylm_gx[lm] + rhaty*Ylm_gy[lm] + rhatz*Ylm_gz[lm] ) * r_power;
          lapl[ib] += (local_l + ( local_g * (static_cast<ST>(2)-l_val) - Vpart ) * rinv) * Ylm_rescale
                    + (local_g - Vpart ) * rhat_dot_G;
        }
        local_val+=Npad;
        local_grad+=Npad;
        local_lapl+=Npad;
      }
    }
    else if(r>rmin)
    {
      // the possibility of reaching here is very very low
      std::cout << "Warning: an electron is very close to an ion, distance=" << r << " be careful!" << std::endl;
      // near core, kill divergence in the laplacian
      r_power_minus_l[0]=cone;
      ST r_power_temp=cone;
      for(int l=1; l<=lmax; l++)
      {
        r_power_temp*=rinv;
        for(int m=-l, lm=l*l; m<=l; m++,lm++)
          r_power_minus_l[lm]=r_power_temp;
      }

      for(size_t lm=0; lm<lm_tot; lm++)
      {
        const ST& l_val=l_vals[lm];
        const ST& r_power=r_power_minus_l[lm];
        const ST Ylm_rescale=Ylm_v[lm]*r_power;
        #pragma omp simd aligned(val,g0,g1,g2,lapl,local_val,local_grad,local_lapl)
        for(size_t ib=0; ib<myV.size(); ib++)
        {
          const ST local_v=local_val[ib];
          const ST local_g=local_grad[ib];
          const ST local_l=local_lapl[ib];
          // value
          const ST Vpart = Ylm_rescale*local_v;
          val[ib] += Vpart;

          // grad
          g0[ib] += local_g * rhatx * Ylm_rescale + local_v * Ylm_gx[lm] * r_power - l_val * rhatx * Vpart * rinv;
          g1[ib] += local_g * rhaty * Ylm_rescale + local_v * Ylm_gy[lm] * r_power - l_val * rhaty * Vpart * rinv;
          g2[ib] += local_g * rhatz * Ylm_rescale + local_v * Ylm_gz[lm] * r_power - l_val * rhatz * Vpart * rinv;

          // laplacian
          ST rhat_dot_G = (Ylm_gx[lm] * rhatx + Ylm_gy[lm] * rhaty + Ylm_gz[lm] * rhatz ) * r_power * r;
          lapl[ib] += local_l * (cone - chalf *l_val) * ( static_cast<ST>(3) * Ylm_rescale + rhat_dot_G );
        }
        local_val+=Npad;
        local_grad+=Npad;
        local_lapl+=Npad;
      }
    }
    else
    {
      std::cout << "Warning: an electron is on top of an ion!" << std::endl;
      // strictly zero

      #pragma omp simd aligned(val,lapl,local_val,local_lapl)
      for(size_t ib=0; ib<myV.size(); ib++)
      {
        // value
        val[ib] = Ylm_v[0]*local_val[ib];

        // laplacian
        lapl[ib] = local_lapl[ib] * static_cast<ST>(3) * Ylm_v[0];
      }
      local_val+=Npad;
      local_grad+=Npad;
      local_lapl+=Npad;
      if(lm_tot>0)
      {
        //std::cout << std::endl;
        for(size_t lm=1; lm<4; lm++)
        {
          #pragma omp simd aligned(g0,g1,g2,local_grad)
          for(size_t ib=0; ib<myV.size(); ib++)
          {
            const ST local_g=local_grad[ib];
            // grad
            g0[ib] += local_g * Ylm_gx[lm];
            g1[ib] += local_g * Ylm_gy[lm];
            g2[ib] += local_g * Ylm_gz[lm];
          }
          local_grad+=Npad;
        }
      }
    }
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const ST& r, const PointType& dr, VV& myV, GV& myG, HT& myH)
  {
    //Needed to do tensor product here
    APP_ABORT("AtomicOrbitalSoA::evaluate_vgh");
  }
};

/** adoptor class to match
 *
 */
template<typename ST>
struct HybridAdoptorBase
{
  static const int D=3;
  using PointType=typename AtomicOrbitalSoA<ST>::PointType;

  // atomic centers
  std::vector<AtomicOrbitalSoA<ST> > AtomicCenters;
  ///table index
  int myTableID;
  //mapping supercell to primitive cell
  std::vector<int> Super2Prim;
  // r, dr for distance table
  DistanceTableData::RealType dist_r;
  DistanceTableData::PosType dist_dr;
  // local r, dr
  ST r;
  PointType dr;
  PointType r_image;

  HybridAdoptorBase() { }

  void set_info(const ParticleSet& ions, ParticleSet& els, const std::vector<int>& mapping)
  {
    myTableID=els.addTable(ions,DT_SOA);
    Super2Prim=mapping;
  }

  inline void resizeStorage(size_t Nb)
  {
    for(int ic=0; ic<AtomicCenters.size(); ic++)
      AtomicCenters[ic].resizeStorage(Nb);
  }

  void bcast_tables(Communicate* comm)
  {
    for(int ic=0; ic<AtomicCenters.size(); ic++)
      AtomicCenters[ic].bcast_tables(comm);
  }

  void reduce_atomic_tables(Communicate* comm)
  {
    for(int ic=0; ic<AtomicCenters.size(); ic++)
      AtomicCenters[ic].reduce_tables(comm);
  }

  inline void flush_zero()
  {
    for(int ic=0; ic<AtomicCenters.size(); ic++)
      AtomicCenters[ic].flush_zero();
  }

  bool read_splines(hdf_archive& h5f)
  {
    bool success=true;
    size_t ncenter;

    success = success && h5f.push("atomic_centers",false);
    success = success && h5f.read(ncenter,"number_of_centers");
    if(!success) return success;
    if(ncenter!=AtomicCenters.size()) success=false;
    // read splines of each center
    for(int ic=0; ic<AtomicCenters.size(); ic++)
    {
      std::ostringstream gname;
      gname << "center_" << ic;
      success = success && h5f.push(gname.str().c_str(),false);
      success = success && AtomicCenters[ic].read_splines(h5f);
      h5f.pop();
    }
    h5f.pop();
    return success;
  }

  bool write_splines(hdf_archive& h5f)
  {
    bool success=true;
    int ncenter=AtomicCenters.size();
    success = success && h5f.push("atomic_centers",true);
    success = success && h5f.write(ncenter,"number_of_centers");
    // write splines of each center
    for(int ic=0; ic<AtomicCenters.size(); ic++)
    {
      std::ostringstream gname;
      gname << "center_" << ic;
      success = success && h5f.push(gname.str().c_str(),true);
      success = success && AtomicCenters[ic].write_splines(h5f);
      h5f.pop();
    }
    h5f.pop();
    return success;
  }

  inline int get_bc_sign(const PointType& r, TinyVector<int,D>& HalfG)
  {
    int bc_sign=0;
    for(int i=0; i<D; i++)
    {
      ST img = round(r[i]-r_image[i]);
      bc_sign += HalfG[i] * (int)img;
    }
    return bc_sign;
  }

  //evaluate only V
  template<typename VV>
  inline bool evaluate_v(const ParticleSet& P, const int iat, VV& myV)
  {
    bool inAtom=false;
    const auto* ei_dist=P.DistTables[myTableID];
    const int center_idx=ei_dist->get_first_neighbor(iat, dist_r, dist_dr, P.activePtcl==iat);
    if(center_idx<0) abort();
    r=dist_r;
    dr[0]=-dist_dr[0];
    dr[1]=-dist_dr[1];
    dr[2]=-dist_dr[2];
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    r_image=myCenter.pos+dr;
    if ( r < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_v(r, dr, myV);
    }
    return inAtom;
  }

  //evaluate only VGL
  template<typename VV, typename GV>
  inline bool evaluate_vgl(const ParticleSet& P, const int iat, VV& myV, GV& myG, VV& myL)
  {
    bool inAtom=false;
    const auto* ei_dist=P.DistTables[myTableID];
    const int center_idx=ei_dist->get_first_neighbor(iat, dist_r, dist_dr, P.activePtcl==iat);
    if(center_idx<0) abort();
    r=dist_r;
    dr[0]=-dist_dr[0];
    dr[1]=-dist_dr[1];
    dr[2]=-dist_dr[2];
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    r_image=myCenter.pos+dr;
    if ( r < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_vgl(r, dr, myV, myG, myL);
    }
    return inAtom;
  }

  //evaluate only VGH
  template<typename VV, typename GV, typename HT>
  inline bool evaluate_vgh(const ParticleSet& P, const int iat, VV& myV, GV& myG, HT& myH)
  {
    bool inAtom=false;
    const auto* ei_dist=P.DistTables[myTableID];
    const int center_idx=ei_dist->get_first_neighbor(iat, dist_r, dist_dr, P.activePtcl==iat);
    if(center_idx<0) abort();
    r=dist_r;
    dr[0]=-dist_dr[0];
    dr[1]=-dist_dr[1];
    dr[2]=-dist_dr[2];
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    r_image=myCenter.pos+dr;
    if ( r < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_vgh(r, dr, myV, myG, myH);
    }
    return inAtom;
  }
};

}
#endif
