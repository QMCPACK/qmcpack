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

#ifndef QMCPLUSPLUS_AFQMC_PHMSD_HPP
#define QMCPLUSPLUS_AFQMC_PHMSD_HPP

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/config.h"
#include "mpi3/shm/mutex.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/array_of_sequences.hpp"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"

#include "AFQMC/Wavefunctions/phmsd_helpers.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"


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
 * NOTE: Optimization note: CLOSED and NONCOLLINEAR calculations with a unique reference
 *                          only need a single set of unique overlaps/determinants/energies.
 *                          Fix this to improve performance!!!! 
 * THERE IS A PROBLEM WITH CLOSED SHELL CALCULATIONS!!!!
 * INCONSISTENCY WHEN SPIN DEPENDENT QUANTITIES ARE REUESTED, e.g. G. 
 * SOLUTION: For CLOSED and NONCOLLINEAR, force a single reference PHMSD.
 *           Then, only calculate alpha component of things and assume beta==alpha when needed,
 *           even if this is not true.
 */
class PHMSD : public AFQMCInfo
{
  // allocators
  using Allocator          = std::allocator<ComplexType>;     //device_allocator<ComplexType>;
  using SPAllocator        = std::allocator<SPComplexType>;   //device_allocator<ComplexType>;
  using Allocator_shared   = shared_allocator<ComplexType>;   //localTG_allocator<ComplexType>;
  using SPAllocator_shared = shared_allocator<SPComplexType>; //localTG_allocator<ComplexType>;

  // type defs
  using pointer              = typename Allocator::pointer;
  using const_pointer        = typename Allocator::const_pointer;
  using pointer_shared       = typename Allocator_shared::pointer;
  using const_pointer_shared = typename Allocator_shared::const_pointer;

  using CVector       = boost::multi::array<ComplexType, 1, Allocator>;
  using CMatrix       = boost::multi::array<ComplexType, 2, Allocator>;
  using SPCMatrix     = boost::multi::array<SPComplexType, 2, SPAllocator>;
  using CTensor       = boost::multi::array<ComplexType, 3, Allocator>;
  using CVector_ref   = boost::multi::array_ref<ComplexType, 1, pointer>;
  using CMatrix_ref   = boost::multi::array_ref<ComplexType, 2, pointer>;
  using CMatrix_cref  = boost::multi::array_ref<const ComplexType, 2, const_pointer>;
  using CTensor_ref   = boost::multi::array_ref<ComplexType, 3, pointer>;
  using CTensor_cref  = boost::multi::array_ref<const ComplexType, 3, const_pointer>;
  using shmCVector    = boost::multi::array<ComplexType, 1, Allocator_shared>;
  using shmCMatrix    = boost::multi::array<ComplexType, 2, Allocator_shared>;
  using shmC3Tensor   = boost::multi::array<ComplexType, 3, Allocator_shared>;
  using shmSPCMatrix  = boost::multi::array<SPComplexType, 2, SPAllocator_shared>;
  using shmSPC3Tensor = boost::multi::array<SPComplexType, 3, SPAllocator_shared>;
  using shmC4Tensor   = boost::multi::array<ComplexType, 4, Allocator_shared>;
  using shared_mutex  = boost::mpi3::shm::mutex;
  using index_aos     = ma::sparse::array_of_sequences<int, int, shared_allocator<int>, ma::sparse::is_root>;

  using stdCVector  = boost::multi::array<ComplexType, 1>;
  using stdCMatrix  = boost::multi::array<ComplexType, 2>;
  using stdCTensor  = boost::multi::array<ComplexType, 3>;
  using mpi3CVector = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;

public:
  PHMSD(AFQMCInfo& info,
        xmlNodePtr cur,
        afqmc::TaskGroup_& tg_,
        HamiltonianOperations&& hop_,
        std::map<int, int>&& acta2mo_,
        std::map<int, int>&& actb2mo_,
        ph_excitations<int, ComplexType>&& abij_,
        index_aos&& beta_coupled_to_unique_alpha__,
        index_aos&& alpha_coupled_to_unique_beta__,
        std::vector<PsiT_Matrix>&& orbs_,
        WALKER_TYPES wlk,
        ValueType nce,
        int targetNW = 1)
      : AFQMCInfo(info),
        TG(tg_),
        SDetOp(SlaterDetOperations_shared<ComplexType>(((wlk != NONCOLLINEAR) ? (NMO) : (2 * NMO)),
                                                       ((wlk != NONCOLLINEAR) ? (NAEA) : (NAEA + NAEB)))),
        HamOp(std::move(hop_)),
        acta2mo(std::move(acta2mo_)),
        actb2mo(std::move(actb2mo_)),
        abij(std::move(abij_)),
        OrbMats(std::move(orbs_)),
        RefOrbMats({0, 0}, shared_allocator<ComplexType>{TG.Node()}),
        number_of_references(-1),
        shmbuff_for_E(nullptr),
        mutex(std::make_unique<shared_mutex>(TG.TG_local())),
        walker_type(wlk),
        NuclearCoulombEnergy(nce),
        last_number_extra_tasks(-1),
        last_task_index(-1),
        local_group_comm(),
        shmbuff_for_G(nullptr),
        maxn_unique_confg(std::max(abij.number_of_unique_excitations()[0], abij.number_of_unique_excitations()[1])),
        maxnactive(std::max(OrbMats[0].size(), OrbMats.back().size())),
        max_exct_n(std::max(abij.maximum_excitation_number()[0], abij.maximum_excitation_number()[1])),
        unique_overlaps({2, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        unique_Etot({2, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        QQ0inv0({1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        QQ0inv1({1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        GA2D0_shm({1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        GB2D0_shm({1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        local_ov  ({2, static_cast<boost::multi::size_t>(maxn_unique_confg)}),
        local_etot({2, static_cast<boost::multi::size_t>(maxn_unique_confg)}),
        local_QQ0inv0({static_cast<boost::multi::size_t>(OrbMats[0].size()), NAEA}),
        local_QQ0inv1({static_cast<boost::multi::size_t>(OrbMats.back().size()), NAEB}),
        Qwork({2 * static_cast<boost::multi::size_t>(max_exct_n), static_cast<boost::multi::size_t>(max_exct_n)}),
        Gwork({NAEA, static_cast<boost::multi::size_t>(maxnactive)}),
        Ovmsd({1, 1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        Emsd({1, 1, 1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        QQ0A({1, 1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        QQ0B({1, 1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        GrefA({1, 1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        GrefB({1, 1, 1}, shared_allocator<ComplexType>{TG.TG_local()}),
        KEright({1, 1, 1}, shared_allocator<SPComplexType>{TG.TG_local()}),
        KEleft({1, 1}, shared_allocator<SPComplexType>{TG.TG_local()}),
        det_couplings{std::move(beta_coupled_to_unique_alpha__), std::move(alpha_coupled_to_unique_beta__)},
        req_Gsend(MPI_REQUEST_NULL),
        req_Grecv(MPI_REQUEST_NULL),
        req_SMsend(MPI_REQUEST_NULL),
        req_SMrecv(MPI_REQUEST_NULL)
  {
    /* To me, PHMSD is not compatible with walker_type=CLOSED unless
       * the MSD expansion is symmetric with respect to spin. For this, 
       * it is better to write a specialized class that assumes either spin symmetry
       * or e.g. Perfect Pairing.
       */
    if (walker_type == CLOSED)
      APP_ABORT("Error: PHMSD requires walker_type != CLOSED.\n");

    compact_G_for_vbias     = true;
    transposed_G_for_vbias_ = HamOp.transposed_G_for_vbias();
    transposed_G_for_E_     = HamOp.transposed_G_for_E();
    transposed_vHS_         = HamOp.transposed_vHS();
    fast_ph_energy          = HamOp.fast_ph_energy();

    excitedState = false;
    std::string excited_file("");
    int i_ = -1, a_ = -1;
    ParameterSet m_param;
    m_param.add(number_of_references, "number_of_references");
    m_param.add(number_of_references, "nrefs");
    m_param.add(excited_file, "excited");
    // generalize this to multi-particle excitations, how do I read a list of integers???
    m_param.add(i_, "i");
    m_param.add(a_, "a");
    m_param.put(cur);

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
      readWfn(excited_file, excitedOrbMat, NMO, maxOccupExtendedMat.first, maxOccupExtendedMat.second);
    }
  }

  ~PHMSD()
  {
    if (req_SMrecv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_SMrecv);
    if (req_SMsend != MPI_REQUEST_NULL)
      MPI_Request_free(&req_SMsend);
    if (req_Grecv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_Grecv);
    if (req_Gsend != MPI_REQUEST_NULL)
      MPI_Request_free(&req_Gsend);
  }

  PHMSD(PHMSD const& other) = delete;
  PHMSD& operator=(PHMSD const& other) = delete;
  PHMSD(PHMSD&& other)                 = default;
  PHMSD& operator=(PHMSD&& other) = delete;

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
    assert(std::get<0>(v.sizes()) == HamOp.local_number_of_cholesky_vectors());
    double scl = (walker_type == COLLINEAR) ? 0.5 : 1.0;
    if (transposed_G_for_vbias_)
    {
      assert(std::get<0>(G.sizes()) == std::get<1>(v.sizes()));
      assert(std::get<1>(G.sizes()) == size_of_G_for_vbias());
      HamOp.vbias(G(G.extension(), {0, long(OrbMats[0].size() * NMO)}), std::forward<MatA>(v), scl * a, 0.0);
      if (walker_type == COLLINEAR) {
        APP_ABORT(" Error in PHMSD::vbias: transposed_G_for_vbias_ should be false. \n");
        HamOp.vbias(G(G.extension(), {long(OrbMats[0].size() * NMO), std::get<1>(G.sizes())}),                                       std::forward<MatA>(v), scl * a, 1.0);
      }
    }
    else
    {
      assert(G.size() == size_of_G_for_vbias());
      assert(std::get<1>(G.sizes()) == std::get<1>(v.sizes()));
      HamOp.vbias(G.sliced(0, OrbMats[0].size() * NMO), std::forward<MatA>(v), scl * a, 0.0);
      if (walker_type == COLLINEAR)
        HamOp.vbias(G.sliced(OrbMats[0].size() * NMO, G.size()), std::forward<MatA>(v), scl * a, 1.0);
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
    if (ovlp.num_elements() != nw)
      ovlp.reextent(iextensions<1u>{nw});
    if (std::get<0>(eloc.sizes()) != nw || std::get<1>(eloc.sizes()) != 3)
      eloc.reextent({nw, 3});
    Energy(wset, eloc, ovlp);
    TG.local_barrier();
    if (TG.getLocalTGRank() == 0)
    {
      int p = 0;
      for (typename WlkSet::iterator it = wset.begin(); it != wset.end(); ++it, ++p)
      {
        *it->overlap() = ovlp[p];
        *it->E1()      = eloc[p][0];
        *it->EXX()     = eloc[p][1];
        *it->EJ()      = eloc[p][2];
      }
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
    if (ovlp.num_elements() != nw)
      ovlp.reextent(iextensions<1u>{nw});
    MixedDensityMatrix(wset, std::forward<MatG>(G), ovlp, compact, transpose);
  }

  template<class WlkSet, class MatG, class TVec>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact = true, bool transpose = false);

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
    /*
      if(nbatch != 0)
        DensityMatrix_batched(wset,std::forward<MatA>(RefA),std::forward<MatB>(RefB),
                                std::forward<MatG>(G),std::forward<TVec>(Ov),
                                herm,compact,transposed);
      else
*/
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
    /*
      if(nbatch != 0)
        DensityMatrix_batched(Left,Right,G,std::forward<TVec>(Ov),LogOverlapFactor,
                                herm,compact);
      else
*/
    DensityMatrix_shared(Left, Right, G, std::forward<TVec>(Ov), LogOverlapFactor, herm, compact);
  }

  template<class WlkSet, class MatG, class CVec1, class CVec2, class Mat1, class Mat2>
  void WalkerAveragedDensityMatrix(const WlkSet& wset,
                                   CVec1& wgt,
                                   MatG& G,
                                   CVec2& denom,
                                   Mat1&& Ovlp,
                                   Mat2&& DMsum,
                                   bool free_projection                          = false,
                                   boost::multi::array_ref<ComplexType, 3>* Refs = nullptr,
                                   boost::multi::array<ComplexType, 2>* detR     = nullptr);

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
    if (ovlp.num_elements() != nw)
      ovlp.reextent(iextensions<1u>{nw});
    MixedDensityMatrix(wset, std::forward<MatG>(G), ovlp, compact_G_for_vbias, transposed_G_for_vbias_);
  }

  /*
     * Calculates the overlaps of all walkers in the set. Returns values in arrays. 
     */
  template<class WlkSet, class TVec>
  void Overlap(const WlkSet& wset, TVec&& Ov);

  /*
     * Calculates the overlaps of all walkers in the set. Updates values in wset. 
     */
  template<class WlkSet>
  void Overlap(WlkSet& wset)
  {
    int nw = wset.size();
    if (ovlp.num_elements() != nw)
      ovlp.reextent(iextensions<1u>{nw});
    Overlap(wset, ovlp);
    TG.local_barrier();
    if (TG.getLocalTGRank() == 0)
    {
      int p = 0;
      for (typename WlkSet::iterator it = wset.begin(); it != wset.end(); ++it, ++p)
        *it->overlap() = ovlp[p];
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
  void Orthogonalize(WlkSet& wset, bool impSamp);

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
      return abij.number_of_configurations();
  }

  ComplexType getReferenceWeight(int i) const { return std::get<2>(*abij.configuration(i)); }

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
      int ncol(NAEA + NAEB); //careful here, spins are stored contiguously
      RefOrbMats.reextent({ndet, nrow * ncol});
      TG.Node().barrier(); // for safety
      if (TG.Node().root())
      {
        boost::multi::array<ComplexType, 2> OA_({
			static_cast<boost::multi::size_t>(std::get<1>(OrbMats[0].sizes())),
			static_cast<boost::multi::size_t>(std::get<0>(OrbMats[0].sizes()))
		});
        boost::multi::array<ComplexType, 2> OB_({0, 0});
        if (OrbMats.size() > 1)
          OB_.reextent({
            static_cast<boost::multi::size_t>(std::get<1>(OrbMats[1].sizes())),
            static_cast<boost::multi::size_t>(std::get<0>(OrbMats[1].sizes()))
          });
        ma::Matrix2MAREF('H', OrbMats[0], OA_);
        if (OrbMats.size() > 1)
          ma::Matrix2MAREF('H', OrbMats[1], OB_);
        std::vector<int> Ac(NAEA);
        std::vector<int> Bc(NAEB);
        for (int i = 0; i < ndet; ++i)
        {
          auto c(abij.configuration(i));
          abij.get_configuration(0, std::get<0>(*c), Ac);
          abij.get_configuration(1, std::get<1>(*c), Bc);
          boost::multi::array_ref<ComplexType, 2> A_(to_address(RefOrbMats[i].origin()), {NMO, NAEA});
          boost::multi::array_ref<ComplexType, 2> B_(A_.origin() + A_.num_elements(), {NMO, NAEB});
          for (int i = 0, ia = 0; i < NMO; ++i)
            for (int a = 0; a < NAEA; ++a, ia++)
              A_[i][a] = OA_[i][Ac[a]];
          if (OrbMats.size() > 1)
          {
            for (int i = 0, ia = 0; i < NMO; ++i)
              for (int a = 0; a < NAEB; ++a, ia++)
                B_[i][a] = OB_[i][Bc[a]];
          }
          else
          {
            for (int i = 0, ia = 0; i < NMO; ++i)
              for (int a = 0; a < NAEB; ++a, ia++)
                B_[i][a] = OA_[i][Bc[a]];
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

  //SlaterDetOperations_shared<ComplexType> SDetOp;
  SlaterDetOperations SDetOp;

  HamiltonianOperations HamOp;

  std::map<int, int> acta2mo;
  std::map<int, int> actb2mo;

  ph_excitations<int, ComplexType> abij;

  // eventually switched from CMatrix to SMHSparseMatrix(node)
  std::vector<PsiT_Matrix> OrbMats;
  mpi3CMatrix RefOrbMats;
  int number_of_references;

  std::unique_ptr<shmCVector> shmbuff_for_E;

  std::unique_ptr<shared_mutex> mutex;

  // in both cases below: closed_shell=0, UHF/ROHF=1, GHF=2
  WALKER_TYPES walker_type;

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
  std::unique_ptr<shmCVector> shmbuff_for_G;

  // shared memory arrays for temporary calculations
  bool fast_ph_energy;
  size_t maxn_unique_confg; // maximum number of unque configurations
  size_t maxnactive;        // maximum number of states in active space
  size_t max_exct_n;        // maximum excitation number (number of electrons excited simultaneously)
  // used by OVerlap and MixedDensityMatrix
  shmCMatrix unique_overlaps;
  shmCMatrix unique_Etot;
  shmCMatrix QQ0inv0; // Q * inv(Q0)
  shmCMatrix QQ0inv1; // Q * inv(Q0)
  shmCMatrix GA2D0_shm;
  shmCMatrix GB2D0_shm;
  boost::multi::array<ComplexType, 2> local_ov;
  boost::multi::array<ComplexType, 2> local_etot;
  boost::multi::array<ComplexType, 2> local_QQ0inv0;
  boost::multi::array<ComplexType, 2> local_QQ0inv1;
  boost::multi::array<ComplexType, 2> Qwork;
  boost::multi::array<ComplexType, 2> Gwork;
  // used by Energy_shared
  boost::multi::array<ComplexType, 1> wgt;
  boost::multi::array<ComplexType, 1> opSpinEJ;
  shmC3Tensor Ovmsd; // [nspins][maxn_unique_confg][nwalk]
  shmC4Tensor Emsd;  // [nspins][maxn_unique_confg][nwalk][3]
  shmC3Tensor QQ0A;  // [nwalk][NAOA][NAEA]
  shmC3Tensor QQ0B;  // [nwalk][NAOB][NAEB]
  shmC3Tensor GrefA; // [nwalk][NAOA][NMO]
  shmC3Tensor GrefB; // [nwalk][NAOB][NMO]
  shmSPC3Tensor KEright;
  shmSPCMatrix KEleft;


  // array of sequence structure storing the list of connected alpha/beta configurations
  std::array<index_aos, 2> det_couplings;

  // excited states
  bool excitedState;
  std::vector<std::pair<int, int>> excitations;
  boost::multi::array<ComplexType, 3> excitedOrbMat;
  CMatrix extendedMatAlpha;
  CMatrix extendedMatBeta;
  std::pair<int, int> maxOccupExtendedMat;
  std::pair<int, int> numExcitations;

  // buffers and work arrays, careful here!!!
  CVector ovlp, localGbuff, ovlp2;
  CMatrix eloc, eloc2, eloc3;

  MPI_Request req_Gsend, req_Grecv;
  MPI_Request req_SMsend, req_SMrecv;

  /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
  template<class WlkSet, class Mat, class TVec>
  void Energy_shared(const WlkSet& wset, Mat&& E, TVec&& Ov);

  /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
  template<class WlkSet, class Mat, class TVec>
  void Energy_distributed(const WlkSet& wset, Mat&& E, TVec&& Ov);

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

  /*
    template<class WlkSet, class MatA, class MatB, class MatG, class TVec>
    void DensityMatrix_batched(const WlkSet& wset, MatA&& RefsA, MatB&& RefsB, MatG&& G,
                                TVec&& Ov, bool herm, bool compact, bool transposed);
*/

  template<class MatA, class MatB, class MatG, class TVec>
  void DensityMatrix_shared(std::vector<MatA>& Left,
                            std::vector<MatB>& Right,
                            std::vector<MatG>& G,
                            TVec&& Ov,
                            double LogOverlapFactor,
                            bool herm,
                            bool compact);

  /*
    template<class MatA, class MatB, class MatG, class TVec>
    void DensityMatrix_batched(std::vector<MatA>& Left, std::vector<MatB>& Right, std::vector<MatG>& G,
                        TVec&& Ov, double LogOverlapFactor, bool herm, bool compact);
*/

  int dm_size(bool full) const
  {
    switch (walker_type)
    {
    case CLOSED: // closed-shell RHF
      return (full) ? (NMO * NMO) : (OrbMats[0].size() * NMO);
      break;
    case COLLINEAR:
      return (full) ? (2 * NMO * NMO) : ((OrbMats[0].size() + OrbMats.back().size()) * NMO);
      break;
    case NONCOLLINEAR:
      return (full) ? (4 * NMO * NMO) : ((OrbMats[0].size()) * 2 * NMO);
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
      return (full) ? (arr{NMO, NMO}) : (arr{OrbMats[0].size(), NMO});
      break;
    case COLLINEAR:
      return (full) ? (arr{NMO, NMO})
                    : ((sp == Alpha) ? (arr{OrbMats[0].size(), NMO}) : (arr{OrbMats.back().size(), NMO}));
      break;
    case NONCOLLINEAR:
      return (full) ? (arr{2 * NMO, 2 * NMO}) : (arr{OrbMats[0].size(), 2 * NMO});
      break;
    default:
      APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
      return arr{-1, -1};
    }
  }
  std::pair<int, int> dm_dims_ref(bool full, SpinTypes sp = Alpha) const
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

#include "AFQMC/Wavefunctions/PHMSD.icc"

#endif
