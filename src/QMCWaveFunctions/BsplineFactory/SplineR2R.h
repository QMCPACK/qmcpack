//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPLINE_R2R_H
#define QMCPLUSPLUS_SPLINE_R2R_H

#include <memory>
#include "BsplineFactory/BsplineSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "spline2/MultiBsplineBase.hpp"
#include "Utilities/FairDivide.h"
#include <ResourceHandle.h>
#include "SplineOMPTargetMultiWalkerMem.h"

namespace qmcplusplus
{
/** class to match ST real spline with BsplineSet::ValueType (real) SPOs
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of the sign of the real part of the phase
 * Internal storage ST type arrays are aligned and padded.
 */
template<typename ST>
class SplineR2R : public BsplineSet
{
public:
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;
  // types for evaluation results
  using TT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector;
  using BsplineSet::GradVector;
  using BsplineSet::HessVector;
  using BsplineSet::ValueVector;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

  template<typename DT>
  using OffloadVector = Vector<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using OffloadPosVector = VectorSoaContainer<DT, 3, OffloadAllocator<DT>>;

private:
  /// if true, use OpenMP offload computation
  const bool use_offload_;
  /// timer for offload portion
  NewTimer& offload_timer_;
  /// if true, gamma point calculation
  bool IsGamma;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///multi bspline set
  std::shared_ptr<MultiBsplineBase<ST>> SplineInst;
  /// const offload copy of GGt
  std::shared_ptr<OffloadVector<ST>> GGt_offload;
  /// const offload copy of GPrimLattice_G
  std::shared_ptr<OffloadVector<ST>> PrimLattice_G_offload;
  /// crowd resource
  ResourceHandle<SplineOMPTargetMultiWalkerMem<ST, TT>> mw_mem_handle_;

  ///Copy of original splines for orbital rotation
  std::shared_ptr<std::vector<ST>> coef_copy_;

  ///thread private ratios for reduction when using nested threading, numVP x numThread
  Matrix<TT> ratios_private;

protected:
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineR2R(const std::string& my_name, bool use_offload = false);
  SplineR2R(const SplineR2R& in);
  virtual std::string getClassName() const override { return "SplineR2R"; }
  virtual std::string getKeyword() const override { return "SplineR2R"; }
  bool isComplex() const override { return false; };
  bool isRotationSupported() const override { return true; }
  virtual bool isOMPoffload() const override { return use_offload_; }

  void createResource(ResourceCollection& collection) const override
  {
    auto resource_index = collection.addResource(std::make_unique<SplineOMPTargetMultiWalkerMem<ST, TT>>());
  }

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const override
  {
    assert(this == &spo_list.getLeader());
    auto& phi_leader          = spo_list.getCastedLeader<SplineR2R>();
    phi_leader.mw_mem_handle_ = collection.lendResource<SplineOMPTargetMultiWalkerMem<ST, TT>>();
  }

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const override
  {
    assert(this == &spo_list.getLeader());
    auto& phi_leader = spo_list.getCastedLeader<SplineR2R>();
    collection.takebackResource(phi_leader.mw_mem_handle_);
  }

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<SplineR2R>(*this); }

  /// Store an original copy of the spline coefficients for orbital rotation
  void storeParamsBeforeRotation() override;

  /*
     Implements orbital rotations via [1,2].
     Should be called by RotatedSPOs::apply_rotation()

     This implementation requires that NSPOs > Nelec. In other words,
     if you want to run a orbopt wfn, you must include some virtual orbitals!

     Some results (using older Berkeley branch) were published in [3].

     [1] Filippi & Fahy, JCP 112, (2000)
     [2] Toulouse & Umrigar, JCP 126, (2007)
     [3] Townsend et al., PRB 102, (2020)
  */
  void applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy) override;

  inline void resizeStorage(size_t n) override
  {
    init_base(n);
    const size_t npad = getAlignedSize<ST>(n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);
    mygH.resize(npad);

    IsGamma = ((HalfG[0] == 0) && (HalfG[1] == 0) && (HalfG[2] == 0));
  }

  void bcast_tables(Communicate* comm) { chunked_bcast(comm, SplineInst->getSplinePtr()); }

  void gather_tables(Communicate* comm)
  {
    if (comm->size() == 1)
      return;
    const int Nbands      = kPoints.size();
    const int Nbandgroups = comm->size();
    offset.resize(Nbandgroups + 1, 0);
    FairDivideLow(Nbands, Nbandgroups, offset);
    gatherv(comm, SplineInst->getSplinePtr(), SplineInst->getSplinePtr()->z_stride, offset);
  }

  template<typename BCT>
  void create_spline(const Ugrid xyz_g[3], const BCT& xyz_bc)
  {
    GGt = dot(transpose(PrimLattice.G), PrimLattice.G);
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  /// this routine can not be called from threaded region
  void finalizeConstruction() override;

  inline void flush_zero() { SplineInst->flush_zero(); }

  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level);

  bool read_splines(hdf_archive& h5f);

  bool write_splines(hdf_archive& h5f);

  /** convert position in PrimLattice unit and return sign */
  inline int convertPos(const PointType& r, PointType& ru) const
  {
    ru          = PrimLattice.toUnit(r);
    int bc_sign = 0;
    for (int i = 0; i < D; i++)
      if (-std::numeric_limits<ST>::epsilon() < ru[i] && ru[i] < 0)
        ru[i] = ST(0.0);
      else
      {
        ST img = std::floor(ru[i]);
        ru[i] -= img;
        bc_sign += HalfG[i] * (int)img;
      }
    return bc_sign;
  }

  void assign_v(int bc_sign, const vContainer_type& myV, ValueVector& psi, int first, int last) const;

  void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override;

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<TT>& ratios) override;

  void mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                            const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                            const RefVector<ValueVector>& psi_list,
                            const std::vector<const ValueType*>& invRow_ptr_list,
                            std::vector<std::vector<ValueType>>& ratios_list) const override;

  void assign_vgl(int bc_sign, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi, int first, int last) const;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(int bc_sign, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  void evaluateVGL(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   ValueVector& d2psi) override;

  void mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                      const RefVectorWithLeader<ParticleSet>& P_list,
                                      int iat,
                                      const std::vector<const ValueType*>& invRow_ptr_list,
                                      OffloadMWVGLArray& phi_vgl_v,
                                      std::vector<ValueType>& ratios,
                                      std::vector<GradType>& grads) const override;

  void assign_vgh(int bc_sign, ValueVector& psi, GradVector& dpsi, HessVector& grad_grad_psi, int first, int last)
      const;

  void evaluateVGH(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override;

  void assign_vghgh(int bc_sign,
                    ValueVector& psi,
                    GradVector& dpsi,
                    HessVector& grad_grad_psi,
                    GGGVector& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const;

  void evaluateVGHGH(const ParticleSet& P,
                     const int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override;

  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend class BsplineReader;
};

extern template class SplineR2R<float>;
extern template class SplineR2R<double>;

} // namespace qmcplusplus
#endif
