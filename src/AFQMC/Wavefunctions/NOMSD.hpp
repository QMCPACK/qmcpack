//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_NOMSD_H
#define QMCPLUSPLUS_AFQMC_NOMSD_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/config.h"
#include "mpi3/shm/mutex.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Walkers/WalkerConfig.hpp"
#include "AFQMC/Memory/buffer_managers.h"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"


namespace qmcplusplus
{
namespace afqmc
{
/*
 * Class that implements a multi-Slater determinant trial wave-function.
 * Single determinant wfns are also allowed. 
 * No relation between different determinants in the expansion is assumed.
 * Designed for non-orthogonal MSD expansions. 
 * For particle-hole orthogonal MSD wfns, use FastMSD.
 */
template<class devPsiT>
class NOMSD : public AFQMCInfo
{
  // Note:
  // if number_of_devices > 0, nextra should always be 0,
  // so code doesn't need to be portable in places guarded by if(nextra>0)

  // allocators
  using Allocator        = device_allocator<ComplexType>;
  using Allocator_shared = localTG_allocator<ComplexType>;

  // type defs
  using pointer              = typename Allocator::pointer;
  using const_pointer        = typename Allocator::const_pointer;
  using pointer_shared       = typename Allocator_shared::pointer;
  using const_pointer_shared = typename Allocator_shared::const_pointer;

  using buffer_alloc_type     = DeviceBufferManager::template allocator_t<ComplexType>;
  using shm_buffer_alloc_type = LocalTGBufferManager::template allocator_t<ComplexType>;

  using CVector      = boost::multi::array<ComplexType, 1, Allocator>;
  using CMatrix      = boost::multi::array<ComplexType, 2, Allocator>;
  using CTensor      = boost::multi::array<ComplexType, 3, Allocator>;
  using CVector_ref  = boost::multi::array_ref<ComplexType, 1, pointer>;
  using CMatrix_ref  = boost::multi::array_ref<ComplexType, 2, pointer>;
  using CMatrix_ptr  = boost::multi::array_ptr<ComplexType, 2, pointer>;
  using CMatrix_cref = boost::multi::array_ref<const ComplexType, 2, const_pointer>;
  using CTensor_ref  = boost::multi::array_ref<ComplexType, 3, pointer>;
  using CTensor_cref = boost::multi::array_ref<const ComplexType, 3, const_pointer>;
  using shmCVector   = boost::multi::array<ComplexType, 1, Allocator_shared>;
  using shmCMatrix   = boost::multi::array<ComplexType, 2, Allocator_shared>;
  using shared_mutex = boost::mpi3::shm::mutex;

  using stdCVector  = boost::multi::array<ComplexType, 1>;
  using stdCMatrix  = boost::multi::array<ComplexType, 2>;
  using stdCTensor  = boost::multi::array<ComplexType, 3>;
  using mpi3CVector = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;

  using stdCMatrix_ref = boost::multi::array_ref<ComplexType, 2>;

  using StaticVector  = boost::multi::static_array<ComplexType, 1, buffer_alloc_type>;
  using StaticMatrix  = boost::multi::static_array<ComplexType, 2, buffer_alloc_type>;
  using Static3Tensor = boost::multi::static_array<ComplexType, 3, buffer_alloc_type>;

  using StaticSHMVector = boost::multi::static_array<ComplexType, 1, shm_buffer_alloc_type>;
  using StaticSHMMatrix = boost::multi::static_array<ComplexType, 2, shm_buffer_alloc_type>;

public:
  template<class MType>
  NOMSD(AFQMCInfo& info,
        xmlNodePtr cur,
        afqmc::TaskGroup_& tg_,
        SlaterDetOperations&& sdet_,
        HamiltonianOperations&& hop_,
        std::vector<ComplexType>&& ci_,
        std::vector<MType>&& orbs_,
        WALKER_TYPES wlk,
        ValueType nce,
        int targetNW = 1)
      : AFQMCInfo(info),
        TG(tg_),
        alloc_(), // right now device_allocator is default constructible
        alloc_shared_(make_localTG_allocator<ComplexType>(TG)),
        buffer_manager(),
        shm_buffer_manager(),
        SDetOp(std::move(sdet_)),
        HamOp(std::move(hop_)),
        ci(std::move(ci_)),
        OrbMats(move_vector<devPsiT>(std::move(orbs_))),
        RefOrbMats({0, 0}, shared_allocator<ComplexType>{TG.Node()}),
        mutex(std::make_unique<shared_mutex>(TG.TG_local())),
        walker_type(wlk),
        nspins((walker_type == COLLINEAR) ? (2) : (1)),
        number_of_references(-1),
        NuclearCoulombEnergy(nce),
        last_number_extra_tasks(-1),
        last_task_index(-1),
        local_group_comm(),
        shmbuff_for_G(nullptr)
  {
    compact_G_for_vbias     = (ci.size() == 1); // this should be input, since it is determined by HOps
    transposed_G_for_vbias_ = HamOp.transposed_G_for_vbias();
    transposed_G_for_E_     = HamOp.transposed_G_for_E();
    transposed_vHS_         = HamOp.transposed_vHS();

    excitedState = false;
    std::string rediag("");
    std::string excited_file("");
    std::string svd_Gf("no");
    std::string svd_O("no");
    std::string svd_Gm("no");
    int i_ = -1, a_ = -1;
    nbatch    = ((number_of_devices() > 0) ? -1 : 0);
    nbatch_qr = ((number_of_devices() > 0) ? -1 : 0);
    if (NMO > 1024 || NAEA > 512)
      nbatch_qr = 0;

    ParameterSet m_param;
    m_param.add(number_of_references, "number_of_references");
    m_param.add(number_of_references, "nrefs");
    m_param.add(excited_file, "excited");
    // generalize this to multi-particle excitations, how do I read a list of integers???
    m_param.add(i_, "i");
    m_param.add(a_, "a");
    m_param.add(rediag, "rediag");
    m_param.add(svd_Gf, "svd_with_Gfull");
    m_param.add(svd_Gm, "svd_with_Gmix");
    m_param.add(svd_O, "svd_with_Ovlp");
    if (TG.TG_local().size() == 1)
      m_param.add(nbatch, "nbatch");
    if (TG.TG_local().size() == 1)
      m_param.add(nbatch_qr, "nbatch_qr");
    m_param.put(cur);

    if (omp_get_num_threads() > 1 && (nbatch == 0))
    {
      app_log() << " WARNING!!!: Found OMP_NUM_THREADS > 1 with nbatch=0.\n"
                << "             This will lead to low performance. Set nbatch. \n";
    }
    std::transform(svd_Gf.begin(), svd_Gf.end(), svd_Gf.begin(), (int (*)(int))tolower);
    if (svd_Gf == "yes" || svd_Gf == "true")
      useSVD_in_Gfull = true;
    std::transform(svd_Gm.begin(), svd_Gm.end(), svd_Gm.begin(), (int (*)(int))tolower);
    if (svd_Gm == "yes" || svd_Gm == "true")
      useSVD_in_Gmix = true;
    std::transform(svd_O.begin(), svd_O.end(), svd_O.begin(), (int (*)(int))tolower);
    if (svd_O == "yes" || svd_O == "true")
      useSVD_in_Ovlp = true;


    if (excited_file != "" && i_ >= 0 && a_ >= 0)
    {
      if (i_ < NMO && a_ < NMO)
      {
        if (i_ >= NAEA || a_ < NAEA)
          APP_ABORT(" Errors: Inconsistent excited orbitals for alpha electrons. \n");
        excitedState        = true;
        maxOccupExtendedMat = {a_, NAEB};
        numExcitations      = {1, 0};
        excitations.push_back({i_, a_});
      }
      else if (i_ >= NMO && a_ >= NMO)
      {
        if (i_ >= NMO + NAEB || a_ < NMO + NAEB)
          APP_ABORT(" Errors: Inconsistent excited orbitals for beta electrons. \n");
        excitedState        = true;
        maxOccupExtendedMat = {NAEA, a_ - NMO};
        numExcitations      = {0, 1};
        excitations.push_back({i_ - NMO, a_ - NMO});
      }
      else
      {
        APP_ABORT(" Errors: Inconsistent excited orbitals. \n");
      }
      stdCTensor excitedOrbMat_;
      readWfn(excited_file, excitedOrbMat_, NMO, maxOccupExtendedMat.first, maxOccupExtendedMat.second);
      excitedOrbMat = excitedOrbMat_;
    }

    std::transform(rediag.begin(), rediag.end(), rediag.begin(), (int (*)(int))tolower);
    if (rediag == "yes" || rediag == "true")
      recompute_ci();
  }

  ~NOMSD() {}

  NOMSD(NOMSD const& other) = delete;
  NOMSD& operator=(NOMSD const& other) = delete;
  NOMSD(NOMSD&& other)                 = default;
  NOMSD& operator=(NOMSD&& other) = delete;

  int local_number_of_cholesky_vectors() const { return HamOp.local_number_of_cholesky_vectors(); }
  int global_number_of_cholesky_vectors() const { return HamOp.global_number_of_cholesky_vectors(); }
  int global_origin_cholesky_vector() const { return HamOp.global_origin_cholesky_vector(); }
  bool distribution_over_cholesky_vectors() const { return HamOp.distribution_over_cholesky_vectors(); }

  int size_of_G_for_vbias() const { return dm_size(!compact_G_for_vbias); }

  bool transposed_G_for_vbias() const { return transposed_G_for_vbias_; }
  bool transposed_G_for_E() const { return transposed_G_for_E_; }
  bool transposed_vHS() const { return transposed_vHS_; }
  WALKER_TYPES getWalkerType() const { return walker_type; }

  template<class Vec>
  void vMF(Vec&& v);

  SlaterDetOperations* getSlaterDetOperations() { return std::addressof(SDetOp); }
  HamiltonianOperations* getHamiltonianOperations() { return std::addressof(HamOp); }

  /*
     * local contribution to vbias for the Green functions in G 
     * G: [size_of_G_for_vbias()][nW]
     * v: [local # Chol. Vectors][nW]
     */
  template<class MatG, class MatA>
  void vbias(const MatG& G, MatA&& v, double a = 1.0)
  {
    if (transposed_G_for_vbias_)
    {
      assert(std::get<0>(G.sizes()) == std::get<1>(v.sizes()));
      assert(std::get<1>(G.sizes()) == size_of_G_for_vbias());
    }
    else
    {
      assert(std::get<0>(G.sizes()) == size_of_G_for_vbias());
      assert(std::get<1>(G.sizes()) == std::get<1>(v.sizes()));
    }
    assert(std::get<0>(v.sizes()) == HamOp.local_number_of_cholesky_vectors());
    if (ci.size() == 1)
    {
      // HamOp expects a compact Gc with alpha/beta components
      HamOp.vbias(G, std::forward<MatA>(v), a);
    }
    else
    {
      if (walker_type == CLOSED)
        HamOp.vbias(G, std::forward<MatA>(v), a); // factor of 2 now in HamOps
      else if (walker_type == NONCOLLINEAR)
        HamOp.vbias(G, std::forward<MatA>(v), a);
      else
      {
        if(transposed_G_for_vbias_)
          APP_ABORT(" Error in NOMSD::vbias: transposed_G_for_vbias_ should be false. \n");
        // HamOp expects either alpha or beta, so must be called twice
        HamOp.vbias(G.sliced(0, NMO * NMO), std::forward<MatA>(v), a, 0.0);
        HamOp.vbias(G.sliced(NMO * NMO, 2 * NMO * NMO), std::forward<MatA>(v), a, 1.0);
      }
    }
    TG.local_barrier();
  }

  /*
     * local contribution to vHS for the Green functions in G 
     * X: [# chol vecs][nW]
     * v: [NMO^2][nW]
     */
  template<class MatX, class MatA>
  void vHS(MatX&& X, MatA&& v, double a = 1.0)
  {
    assert(std::get<0>(X.sizes()) == HamOp.local_number_of_cholesky_vectors());
    if (transposed_vHS_)
      assert(std::get<1>(X.sizes()) == std::get<0>(v.sizes()));
    else
      assert(std::get<1>(X.sizes()) == std::get<1>(v.sizes()));
    HamOp.vHS(std::forward<MatX>(X), std::forward<MatA>(v), a);
    TG.local_barrier();
  }

  /*
     * Calculates the local energy and overlaps of all the walkers in the set and stores
     * them in the wset data
     */
  template<class WlkSet>
  void Energy(WlkSet& wset)
  {
    int nw = wset.size();
    StaticVector ovlp(iextensions<1u>{nw}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    StaticMatrix eloc({nw, 3}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    Energy(wset, eloc, ovlp);
    TG.local_barrier();
    if (TG.getLocalTGRank() == 0)
    {
      wset.setProperty(OVLP, ovlp);
      wset.setProperty(E1_, eloc(eloc.extension(), 0));
      wset.setProperty(EXX_, eloc(eloc.extension(), 1));
      wset.setProperty(EJ_, eloc(eloc.extension(), 2));
    }
    TG.local_barrier();
  }

  /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
  template<class WlkSet, class Mat, class TVec>
  void Energy(const WlkSet& wset, Mat&& E, TVec&& Ov)
  {
    if (TG.getNGroupsPerTG() > 1)
      Energy_distributed(wset, std::forward<Mat>(E), std::forward<TVec>(Ov));
    else
      Energy_shared(wset, std::forward<Mat>(E), std::forward<TVec>(Ov));
  }

  /*
     * Calculates the mixed density matrix for all walkers in the walker set. 
     * Options:
     *  - compact:   If true (default), returns compact form with Dim: [NEL*NMO], 
     *                 otherwise returns full form with Dim: [NMO*NMO]. 
     *  - transpose: If false (default), returns standard form with Dim: [XXX][nW]
     *                 otherwise returns the transpose with Dim: [nW][XXX}
     */
  template<class WlkSet, class MatG>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, bool compact = true, bool transpose = false)
  {
    int nw = wset.size();
    StaticVector ovlp(iextensions<1u>{nw}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    MixedDensityMatrix(wset, std::forward<MatG>(G), ovlp, compact, transpose);
  }

  template<class WlkSet, class MatG, class TVec>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact = true, bool transpose = false)
  {
    if (nbatch != 0)
      MixedDensityMatrix_batched(wset, std::forward<MatG>(G), std::forward<TVec>(Ov), compact, transpose);
    else
      MixedDensityMatrix_shared(wset, std::forward<MatG>(G), std::forward<TVec>(Ov), compact, transpose);
  }

  /*
     * Calculates the density matrix with respect to a given Reference
     * for all walkers in the walker set. 
     */
  template<class WlkSet, class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix(const WlkSet& wset,
                     MatA&& RefA,
                     MatB&& RefB,
                     MatG&& G,
                     TVec&& Ov,
                     bool herm,
                     bool compact,
                     bool transposed)
  {
    if (nbatch != 0)
      DensityMatrix_batched(wset, std::forward<MatA>(RefA), std::forward<MatB>(RefB), std::forward<MatG>(G),
                            std::forward<TVec>(Ov), herm, compact, transposed);
    else
      DensityMatrix_shared(wset, std::forward<MatA>(RefA), std::forward<MatB>(RefB), std::forward<MatG>(G),
                           std::forward<TVec>(Ov), herm, compact, transposed);
  }

  template<class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix(std::vector<MatA>& Left,
                     std::vector<MatB>& Right,
                     std::vector<MatG>& G,
                     TVec&& Ov,
                     double LogOverlapFactor,
                     bool herm,
                     bool compact)
  {
    if (nbatch != 0)
      DensityMatrix_batched(Left, Right, G, std::forward<TVec>(Ov), LogOverlapFactor, herm, compact);
    else
      DensityMatrix_shared(Left, Right, G, std::forward<TVec>(Ov), LogOverlapFactor, herm, compact);
  }

  /*
     * Calculates the walker averaged density matrix.
     * Options:
     *  - free_projection: If false (default), assumes using phaseless approximation
     *                       otherwise assumes using free projection.
     */
  template<class WlkSet, class MatG, class CVec1, class CVec2, class Mat1, class Mat2>
  void WalkerAveragedDensityMatrix(const WlkSet& wset,
                                   CVec1& wgt,
                                   MatG& G,
                                   CVec2& denom,
                                   Mat1&& Ovlp,
                                   Mat2&& DMsum,
                                   bool free_projection                          = false,
                                   boost::multi::array_ref<ComplexType, 3>* Refs = nullptr,
                                   boost::multi::array<ComplexType, 2>* detR     = nullptr)
  {
    //      if(nbatch != 0)
    //        WalkerAveragedDensityMatrix_batched(wset,wgt,G,denom,free_projection,Refs,detR);
    //      else
    // having problems with underflow with large (back) projection times and multidets,
    // mainly from the normalization coming out of the orthonormalization (detR)
    // writing specialized version for single det which doesn;t have this issues
    if (ci.size() == 1)
      WalkerAveragedDensityMatrix_shared_single_det(wset, wgt, G, denom, Ovlp, DMsum, free_projection, Refs, detR);
    else
      WalkerAveragedDensityMatrix_shared(wset, wgt, G, denom, Ovlp, DMsum, free_projection, Refs, detR);
  }

  /*
     * Calculates the mixed density matrix for all walkers in the walker set
     *   with a format consistent with (and expected by) the vbias routine.
     * This is implementation dependent, so this density matrix should ONLY be used
     * in conjunction with vbias. 
     */
  template<class WlkSet, class MatG>
  void MixedDensityMatrix_for_vbias(const WlkSet& wset, MatG&& G)
  {
    int nw = wset.size();
    StaticVector ovlp(iextensions<1u>{nw}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    MixedDensityMatrix(wset, std::forward<MatG>(G), ovlp, compact_G_for_vbias, transposed_G_for_vbias_);
  }

  /*
     * Calculates the overlaps of all walkers in the set. Returns values in arrays. 
     */
  template<class WlkSet, class TVec>
  void Overlap(const WlkSet& wset, TVec&& Ov)
  {
    if (nbatch != 0)
      Overlap_batched(wset, std::forward<TVec>(Ov));
    else
      Overlap_shared(wset, std::forward<TVec>(Ov));
  }

  /*
     * Calculates the overlaps of all walkers in the set. Updates values in wset. 
     */
  template<class WlkSet>
  void Overlap(WlkSet& wset)
  {
    int nw = wset.size();
    StaticVector ovlp(iextensions<1u>{nw}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    Overlap(wset, ovlp);
    TG.local_barrier();
    if (TG.getLocalTGRank() == 0)
    {
      wset.setProperty(OVLP, ovlp);
    }
    TG.local_barrier();
  }

  /*
     * Orthogonalizes the Slater matrices of all walkers in the set.  
     * Options:
     *  - bool importanceSamplingt(default=true): use algorithm appropriate for importance sampling. 
     *         This means that the determinant of the R matrix in the QR decomposition is ignored.
     *         If false, add the determinant of R to the weight of the walker. 
     */
  template<class WlkSet>
  void Orthogonalize(WlkSet& wset, bool impSamp)
  {
    if (not excitedState && (nbatch_qr != 0))
      Orthogonalize_batched(wset, impSamp);
    else
      Orthogonalize_shared(wset, impSamp);
  }

  /*
     * Orthogonalizes the Slater matrix of a walker in an excited state calculation.
     */
  template<class Mat>
  void OrthogonalizeExcited(Mat&& A, SpinTypes spin, double LogOverlapFactor);


  /*
     * Returns the number of reference Slater Matrices needed for back propagation.  
     */
  int number_of_references_for_back_propagation() const
  {
    if (number_of_references > 0)
      return number_of_references;
    else
      return ((walker_type == COLLINEAR) ? OrbMats.size() / 2 : OrbMats.size());
  }

  ComplexType getReferenceWeight(int i) const { return ci[i]; }

  /*
     * Returns the reference Slater Matrices needed for back propagation.  
     */
  template<class Mat, class Ptr = ComplexType*>
  void getReferencesForBackPropagation(Mat&& A)
  {
    static_assert(std::decay<Mat>::type::dimensionality == 2, "Wrong dimensionality");
    int ndet = number_of_references_for_back_propagation();
    assert(A.size() == ndet);
    if (RefOrbMats.size() == 0)
    {
      TG.Node().barrier(); // for safety
      int nrow(NMO * ((walker_type == NONCOLLINEAR) ? 2 : 1));
      int ncol(NAEA + ((walker_type == CLOSED) ? 0 : NAEB)); //careful here, spins are stored contiguously
      RefOrbMats.reextent({ndet, nrow * ncol});
      TG.Node().barrier(); // for safety
      if (TG.Node().root())
      {
        if (walker_type != COLLINEAR)
        {
          for (int i = 0; i < ndet; ++i)
          {
            boost::multi::array_ref<ComplexType, 2> A_(to_address(RefOrbMats[i].origin()), {nrow, ncol});
            ma::Matrix2MAREF('H', OrbMats[i], A_);
          }
        }
        else
        {
          for (int i = 0; i < ndet; ++i)
          {
            boost::multi::array_ref<ComplexType, 2> A_(to_address(RefOrbMats[i].origin()), {NMO, NAEA});
            ma::Matrix2MAREF('H', OrbMats[2 * i], A_);
            boost::multi::array_ref<ComplexType, 2> B_(A_.origin() + A_.num_elements(), {NMO, NAEB});
            ma::Matrix2MAREF('H', OrbMats[2 * i + 1], B_);
          }
        }
      }                    // TG.Node().root()
      TG.Node().barrier(); // for safety
    }
    assert(std::get<0>(RefOrbMats.sizes()) == ndet);
    assert(std::get<1>(RefOrbMats.sizes()) == std::get<1>(A.sizes()));
    auto&& RefOrbMats_(boost::multi::static_array_cast<ComplexType, ComplexType*>(RefOrbMats));
    auto&& A_(boost::multi::static_array_cast<ComplexType, Ptr>(A));
    using std::copy_n;
    int n0, n1;
    std::tie(n0, n1) = FairDivideBoundary(TG.getLocalTGRank(), int(std::get<1>(A.sizes())), TG.getNCoresPerTG());
    for (int i = 0; i < ndet; i++)
      copy_n(RefOrbMats_[i].origin() + n0, n1 - n0, A_[i].origin() + n0);
    TG.TG_local().barrier();
  }

protected:
  TaskGroup_& TG;

  Allocator alloc_;
  Allocator_shared alloc_shared_;

  DeviceBufferManager buffer_manager;
  LocalTGBufferManager shm_buffer_manager;

  //SlaterDetOperations_shared<ComplexType> SDetOp;
  SlaterDetOperations SDetOp;

  HamiltonianOperations HamOp;

  std::vector<ComplexType> ci;

  // eventually switched from CMatrix to SMHSparseMatrix(node)
  std::vector<devPsiT> OrbMats;
  mpi3CMatrix RefOrbMats;

  std::unique_ptr<shared_mutex> mutex;

  // in both cases below: closed_shell=0, UHF/ROHF=1, GHF=2
  WALKER_TYPES walker_type;
  int nspins;

  int number_of_references;

  int nbatch;
  int nbatch_qr;
  bool useSVD_in_Ovlp  = false;
  bool useSVD_in_Gmix  = false;
  bool useSVD_in_Gfull = false;

  bool compact_G_for_vbias;

  // in the 3 cases, true means [nwalk][...], false means [...][nwalk]
  bool transposed_G_for_vbias_;
  bool transposed_G_for_E_;
  bool transposed_vHS_;

  ValueType NuclearCoulombEnergy;

  // not elegant, but reasonable for now
  int last_number_extra_tasks;
  int last_task_index;

  // shared_communicator for parallel work within TG_local()
  //std::unique_ptr<shared_communicator> local_group_comm;
  shared_communicator local_group_comm;
  std::unique_ptr<mpi3CVector> shmbuff_for_G;

  // excited states
  bool excitedState;
  std::vector<std::pair<int, int>> excitations;
  CTensor excitedOrbMat;
  CMatrix extendedMatAlpha;
  CMatrix extendedMatBeta;
  std::pair<int, int> maxOccupExtendedMat;
  std::pair<int, int> numExcitations;

  void recompute_ci();

  /* 
     * Computes the density matrix with respect to a given reference. 
     * Intended to be used in combination with the energy evaluation routine.
     * G and Ov are expected to be in shared memory.
     */
  template<class WlkSet, class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix_shared(const WlkSet& wset,
                            MatA&& RefsA,
                            MatB&& RefsB,
                            MatG&& G,
                            TVec&& Ov,
                            bool herm,
                            bool compact,
                            bool transposed);

  template<class WlkSet, class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix_batched(const WlkSet& wset,
                             MatA&& RefsA,
                             MatB&& RefsB,
                             MatG&& G,
                             TVec&& Ov,
                             bool herm,
                             bool compact,
                             bool transposed);

  template<class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix_shared(std::vector<MatA>& Left,
                            std::vector<MatB>& Right,
                            std::vector<MatG>& G,
                            TVec&& Ov,
                            double LogOverlapFactor,
                            bool herm,
                            bool compact);

  template<class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix_batched(std::vector<MatA>& Left,
                             std::vector<MatB>& Right,
                             std::vector<MatG>& G,
                             TVec&& Ov,
                             double LogOverlapFactor,
                             bool herm,
                             bool compact);

  template<class MatSM, class MatG, class TVec>
  void MixedDensityMatrix_for_E_from_SM(const MatSM& SM, MatG&& G, TVec&& Ov, int nd, double LogOverlapFactor);

  /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
  template<class WlkSet, class Mat, class TVec>
  void Energy_shared(const WlkSet& wset, Mat&& E, TVec&& Ov);

  template<class WlkSet, class MatG, class TVec>
  void MixedDensityMatrix_shared(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact = true, bool transpose = false);

  template<class WlkSet, class MatG, class TVec>
  void MixedDensityMatrix_batched(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact = true, bool transpose = false);

  template<class WlkSet, class TVec>
  void Overlap_shared(const WlkSet& wset, TVec&& Ov);

  template<class WlkSet, class TVec>
  void Overlap_batched(const WlkSet& wset, TVec&& Ov);

  template<class WlkSet>
  void Orthogonalize_batched(WlkSet& wset, bool impSamp);

  template<class WlkSet>
  void Orthogonalize_shared(WlkSet& wset, bool impSamp);

  /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
  template<class WlkSet, class Mat, class TVec>
  void Energy_distributed(const WlkSet& wset, Mat&& E, TVec&& Ov)
  {
    if (ci.size() == 1)
      Energy_distributed_singleDet(wset, std::forward<Mat>(E), std::forward<TVec>(Ov));
    else
      Energy_distributed_multiDet(wset, std::forward<Mat>(E), std::forward<TVec>(Ov));
  }

  template<class WlkSet, class MatG, class CVec1, class CVec2, class Mat1, class Mat2>
  void WalkerAveragedDensityMatrix_batched(const WlkSet& wset,
                                           CVec1& wgt,
                                           MatG& G,
                                           CVec2& denom,
                                           Mat1&& Ovlp,
                                           Mat2&& DMsum,
                                           bool free_projection,
                                           boost::multi::array_ref<ComplexType, 3>* Refs,
                                           boost::multi::array<ComplexType, 2>* detR);

  template<class WlkSet, class MatG, class CVec1, class CVec2, class Mat1, class Mat2>
  void WalkerAveragedDensityMatrix_shared(const WlkSet& wset,
                                          CVec1& wgt,
                                          MatG& G,
                                          CVec2& denom,
                                          Mat1&& Ovlp,
                                          Mat2&& DMsum,
                                          bool free_projection,
                                          boost::multi::array_ref<ComplexType, 3>* Refs,
                                          boost::multi::array<ComplexType, 2>* detR);

  template<class WlkSet, class MatG, class CVec1, class CVec2, class Mat1, class Mat2>
  void WalkerAveragedDensityMatrix_shared_single_det(const WlkSet& wset,
                                                     CVec1& wgt,
                                                     MatG& G,
                                                     CVec2& denom,
                                                     Mat1&& Ovlp,
                                                     Mat2&& DMsum,
                                                     bool free_projection,
                                                     boost::multi::array_ref<ComplexType, 3>* Refs,
                                                     boost::multi::array<ComplexType, 2>* detR);

  template<class WlkSet, class Mat, class TVec>
  void Energy_distributed_singleDet(const WlkSet& wset, Mat&& E, TVec&& Ov);

  template<class WlkSet, class Mat, class TVec>
  void Energy_distributed_multiDet(const WlkSet& wset, Mat&& E, TVec&& Ov);

  int dm_size(bool full) const
  {
    switch (walker_type)
    {
    case CLOSED: // closed-shell RHF
      return (full) ? (NMO * NMO) : (NAEA * NMO);
      break;
    case COLLINEAR:
      return (full) ? (2 * NMO * NMO) : ((NAEA + NAEB) * NMO);
      break;
    case NONCOLLINEAR:
      return (full) ? (4 * NMO * NMO) : (NAEA * 2 * NMO);
      break;
    default:
      APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
      return -1;
    }
  }
  // dimensions for each component of the DM.
  std::pair<int, int> dm_dims(bool full, SpinTypes sp = Alpha) const
  {
    using arr = std::pair<int, int>;
    switch (walker_type)
    {
    case CLOSED: // closed-shell RHF
      return (full) ? (arr{NMO, NMO}) : (arr{NAEA, NMO});
      break;
    case COLLINEAR:
      return (full) ? (arr{NMO, NMO}) : ((sp == Alpha) ? (arr{NAEA, NMO}) : (arr{NAEB, NMO}));
      break;
    case NONCOLLINEAR:
      return (full) ? (arr{2 * NMO, 2 * NMO}) : (arr{NAEA, 2 * NMO});
      break;
    default:
      APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
      return arr{-1, -1};
    }
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#include "AFQMC/Wavefunctions/NOMSD.icc"

#endif
