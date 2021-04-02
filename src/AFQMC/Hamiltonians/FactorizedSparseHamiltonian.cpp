#include <cstdlib>
#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <numeric>
#include <functional>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include "Configuration.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/HSPotential_Helpers.h"
#include "AFQMC/Hamiltonians/generateHijkl.hpp"
#include "AFQMC/Matrix/array_partition.hpp"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"
#include "AFQMC/HamiltonianOperations/SparseTensorIO.hpp"

namespace qmcplusplus
{
namespace afqmc
{
boost::multi::array<ComplexType, 1> FactorizedSparseHamiltonian::halfRotatedHij(WALKER_TYPES type,
                                                                                PsiT_Matrix* Alpha,
                                                                                PsiT_Matrix* Beta)
{
  check_wavefunction_consistency(type, Alpha, Beta, NMO, NAEA, NAEB);
#if defined(QMC_COMPLEX)
  return rotateHij(type, Alpha, Beta, OneBodyHamiltonian::H1);
#else
  boost::multi::array<ComplexType, 2> H1_(OneBodyHamiltonian::H1);
  return rotateHij(type, Alpha, Beta, H1_);
#endif
}

SpVType_shm_csr_matrix FactorizedSparseHamiltonian::calculateHSPotentials(double cut,
                                                                          TaskGroup_& TGprop,
                                                                          boost::multi::array<ComplexType, 2>& vn0)
{
  using Alloc = shared_allocator<SPValueType>;
  if (TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: HSPotential not implemented with distributed Hamiltonian. \n");

  vn0.reextent({NMO, NMO});
  std::fill_n(vn0.data_elements(), NMO * NMO, ComplexType(0));
  for (int i = 0, cnt = 0; i < NMO; i++)
    for (int l = i; l < NMO; l++, cnt++)
    {
      if (cnt % TG.Global().size() != TG.Global().rank())
        continue;
      ValueType vl = ValueType(0);
      for (int j = 0; j < NMO; j++)
        vl += H(i, j, j, l);
      vn0[i][l] -= 0.5 * vl;
      if (i != l)
        vn0[l][i] -= 0.5 * ma::conj(vl);
    }
  TG.Global().all_reduce_in_place_n(vn0.origin(), vn0.num_elements(), std::plus<>());

  if (TG.getNumberOfTGs() > 1)
  {
    APP_ABORT(" Finish HS. \n");
    return SpVType_shm_csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
  }
  else
  {
    // you always need to count since some vectors might be empty
    auto nnz_per_cv = HamHelper::count_nnz_per_cholvec(cut, TG, V2_fact, NMO);

    // count vectors and make mapping from 2*n/2*n+1 to actual value
    std::vector<int> map_;
    map_.reserve(nnz_per_cv.size());
    int cnt = 0;
    for (auto& n : nnz_per_cv)
      map_.emplace_back((n > 0) ? (cnt++) : (-1));

    // partition and
    std::size_t cv0, cvN;
    if (TGprop.getNGroupsPerTG() == 1)
    { // Spvn is not distributed
      cv0 = 0;
      cvN = nnz_per_cv.size();
      if (TG.Global().root())
      {
        app_log() << std::endl << "Partition of cholesky vector over nodes in TG: ";
        app_log() << std::count_if(nnz_per_cv.begin(), nnz_per_cv.begin() + cvN, [](std::size_t i) { return i > 0; });
        app_log() << std::endl;
        app_log() << "Number of terms in Cholesky Matrix per node in TG: ";
        app_log() << std::accumulate(nnz_per_cv.begin(), nnz_per_cv.begin() + cvN, std::size_t(0));
        app_log() << std::endl;
      }
    }
    else
    {
      std::vector<std::size_t> cv_boundaries(TGprop.getNGroupsPerTG() + 1);
      simple_matrix_partition<TaskGroup_, std::size_t, double> split(V2_fact.size(0), nnz_per_cv.size(), cut);
      // no need for all cores to do this
      if (TG.Global().root())
        split.partition(TGprop, false, nnz_per_cv, cv_boundaries);
      TG.Global().broadcast_n(cv_boundaries.begin(), cv_boundaries.size());
      cv0 = cv_boundaries[TGprop.getLocalGroupNumber()];
      cvN = cv_boundaries[TGprop.getLocalGroupNumber() + 1];
      // no need for all cores to do this
      if (TG.Global().root())
      {
        app_log() << std::endl << "Partition of cholesky vector over nodes in TG: ";
        for (int i = 0; i < TGprop.getNGroupsPerTG(); i++)
          app_log() << std::count_if(nnz_per_cv.begin() + cv_boundaries[i], nnz_per_cv.begin() + cv_boundaries[i + 1],
                                     [](std::size_t i) { return i > 0; })
                    << " ";
        app_log() << std::endl;
        app_log() << "Number of terms in Cholesky Matrix per node in TG: ";
        for (int i = 0; i < TGprop.getNGroupsPerTG(); i++)
          app_log() << std::accumulate(nnz_per_cv.begin() + cv_boundaries[i], nnz_per_cv.begin() + cv_boundaries[i + 1],
                                       std::size_t(0))
                    << " ";
        app_log() << std::endl << std::endl;
      }
    }

    auto nnz_per_ik = HamHelper::count_nnz_per_ik(cut, TG, V2_fact, NMO, cv0, cvN);

    std::size_t nvec =
        std::count_if(nnz_per_cv.begin() + cv0, nnz_per_cv.begin() + cvN, [](std::size_t const& i) { return i > 0; });

    std::size_t cv_origin =
        std::count_if(nnz_per_cv.begin(), nnz_per_cv.begin() + cv0, [](std::size_t const& i) { return i > 0; });

    // can build csr directly since cores work on non-overlapping rows
    // and can use emplace_back
    SpVType_shm_csr_matrix csr(tp_ul_ul{NMO * NMO, nvec}, tp_ul_ul{0, cv_origin}, nnz_per_ik, Alloc(TG.Node()));

    HamHelper::generateHSPotential(csr, map_, cut, TG, V2_fact, NMO, cv0, cvN);
    TG.node_barrier();

    return csr;
  }
}

SpCType_shm_csr_matrix FactorizedSparseHamiltonian::halfRotatedHijkl(WALKER_TYPES type,
                                                                     bool addCoulomb,
                                                                     TaskGroup_& TGHam,
                                                                     PsiT_Matrix_t<SPComplexType>* Alpha,
                                                                     PsiT_Matrix_t<SPComplexType>* Beta,
                                                                     RealType const cut)
{
  check_wavefunction_consistency(type, Alpha, Beta, NMO, NAEA, NAEB);
  using Alloc = shared_allocator<SPComplexType>;
  assert(Alpha != nullptr);
  if (type == COLLINEAR)
    assert(Beta != nullptr);
  std::size_t nr = Alpha->size(0) * Alpha->size(1);
  if (type == COLLINEAR)
    nr = (Alpha->size(0) + Beta->size(0)) * Alpha->size(1);
  if (TGHam.getNGroupsPerTG() > 1)
  {
    using tvec = std::vector<std::tuple<int, int, SPComplexType>>;
    tvec tmat;
    tmat.reserve(100000); // reserve some space
    rotateHijkl<tvec>(factorizedHalfRotationType, type, addCoulomb, TG, tmat, Alpha, Beta, V2_fact, cut,
                      maximum_buffer_size, false, false);
    TG.node_barrier();
    return csr::shm::construct_distributed_csr_matrix_from_distributed_containers<SpCType_shm_csr_matrix>(tmat, nr, nr,
                                                                                                          TGHam);
  }
  else
  {
    using ucsr_mat = SpCType_shm_csr_matrix::base;
    ucsr_mat ucsr(tp_ul_ul{nr, nr}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
    if (TG.getTotalNodes() > 1)
      rotateHijkl<ucsr_mat>(factorizedHalfRotationType, type, addCoulomb, TG, ucsr, Alpha, Beta, V2_fact, cut,
                            maximum_buffer_size, true, true);
    else
      rotateHijkl_single_node<ucsr_mat>(factorizedHalfRotationType, type, addCoulomb, TG, ucsr, Alpha, Beta, V2_fact,
                                        cut, maximum_buffer_size, true);
    TG.node_barrier();
    return csr::shm::construct_csr_matrix_from_distributed_ucsr<SpCType_shm_csr_matrix, TaskGroup_>(std::move(ucsr),
                                                                                                    TG);
  }
}

SpVType_shm_csr_matrix FactorizedSparseHamiltonian::generateHijkl(
    WALKER_TYPES type,
    bool addCoulomb,
    TaskGroup_& TGwfn,
    std::map<IndexType, std::pair<bool, IndexType>>& occ_a,
    std::map<IndexType, std::pair<bool, IndexType>>& occ_b,
    RealType const cut)
{
  using Alloc = shared_allocator<SPValueType>;
  int nel     = 0;
  for (auto& i : occ_a)
    if (i.second.first)
      nel++;
  if (type != CLOSED)
    for (auto& i : occ_b)
      if (i.second.first)
        nel++;
  std::size_t nr = nel * NMO;
  if (type == NONCOLLINEAR)
    nr *= 2;
  if (TGwfn.getNGroupsPerTG() > 1)
  {
    APP_ABORT("Finish \n");
    return SpVType_shm_csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
  }
  else
  {
    using ucsr_mat = SpVType_shm_csr_matrix::base;
    using wrapper  = csr::matrix_emplace_wrapper<ucsr_mat>;
    using namespace std::placeholders;
    ucsr_mat ucsr(tp_ul_ul{nr, nr}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
    wrapper ucsr_wrapper(ucsr, TG.Node());
    HamHelper::generateHijkl(type, addCoulomb, TG, ucsr_wrapper, NMO, occ_a, occ_b, *this, true, true, cut);
    ucsr_wrapper.push_buffer();
    TG.node_barrier();
    return csr::shm::construct_csr_matrix_from_distributed_ucsr<SpVType_shm_csr_matrix, TaskGroup_>(std::move(ucsr),
                                                                                                    TG);
  }
}

HamiltonianOperations FactorizedSparseHamiltonian::getHamiltonianOperations(bool pureSD,
                                                                            bool addCoulomb,
                                                                            WALKER_TYPES type,
                                                                            std::vector<PsiT_Matrix>& PsiT,
                                                                            double cutvn,
                                                                            double cutv2,
                                                                            TaskGroup_& TGprop,
                                                                            TaskGroup_& TGwfn,
                                                                            hdf_archive& dump)
{
  if (type == COLLINEAR)
    assert(PsiT.size() % 2 == 0);

#if defined(MIXED_PRECISION)
  auto PsiTsp(csr::shm::CSRvector_to_single_precision<PsiT_Matrix_t<SPComplexType>>(PsiT));
#else
  auto& PsiTsp(PsiT);
#endif

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if (TGwfn.Global().root())
    write_hdf = !dump.closed();
  //    if(TGwfn.Global().root()) write_hdf = (dump.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  boost::multi::array<ComplexType, 2> vn0({NMO, NMO});
  auto Spvn(calculateHSPotentials(cutvn, TGprop, vn0));
  auto Spvnview(csr::shm::local_balanced_partition(Spvn, TGprop));

  // trick: the last node always knows what the total # of chol vecs is
  int global_ncvecs = 0;
  if (TG.getNodeID() == TG.getTotalNodes() - 1 && TG.getCoreID() == 0)
    global_ncvecs = Spvn.global_origin()[1] + Spvn.size(1);
  global_ncvecs = (TG.Global() += global_ncvecs);

  ValueType E0 = OneBodyHamiltonian::NuclearCoulombEnergy + OneBodyHamiltonian::FrozenCoreEnergy;

  // dense one body hamiltonian
  auto H1 = getH1();

  // several posibilities
  int ndet = ((type != COLLINEAR) ? (PsiT.size()) : (PsiT.size() / 2));
  // SparseTensor<Integrals_Type, SpvnT_Type>
  if (ndet == 1)
  {
    std::vector<boost::multi::array<ComplexType, 1>> hij;
    hij.reserve(ndet);
    hij.emplace_back(halfRotatedHij(type, &PsiT[0], ((type == COLLINEAR) ? (&PsiT[1]) : (&PsiT[0]))));
    std::vector<SpCType_shm_csr_matrix> SpvnT;
    using matrix_view = typename SpCType_shm_csr_matrix::template matrix_view<int>;
    std::vector<matrix_view> SpvnTview;
    SpvnT.emplace_back(
        sparse_rotate::halfRotateCholeskyMatrixForBias<SPComplexType>(type, TGprop, &PsiTsp[0],
                                                                      ((type == COLLINEAR) ? (&PsiTsp[1])
                                                                                           : (&PsiTsp[0])),
                                                                      Spvn, cutv2));
    SpvnTview.emplace_back(csr::shm::local_balanced_partition(SpvnT[0], TGprop));

    // in single determinant, SpvnT is the half rotated cholesky matrix
    if (pureSD)
    {
      // in pureSD: V2 is ValueType
      using sparse_ham = SparseTensor<ValueType, ComplexType>;

      if (type != CLOSED)
        APP_ABORT(" Error: pureSD only implemented for closed shell system right now. \n");

      std::map<IndexType, std::pair<bool, IndexType>> occ_a;
      for (int i = 0; i < NMO; i++)
        occ_a[i] = {false, 0};
      auto nel = PsiT[0].size(0);
      for (std::size_t i = 0; i < nel; i++)
      {
        if (PsiT[0].num_non_zero_elements(i) != 1)
          APP_ABORT(" Error: PsiT is not of pureSD form: Too many nnz. \n");
        auto val = PsiT[0].non_zero_values_data(i)[0];
        if (std::abs(val - ValueType(1.0)) > 1e-8)
          APP_ABORT(" Error: PsiT is not of pureSD form: Coeff != 1.0. \n");
        auto orb   = PsiT[0].non_zero_indices2_data(i)[0];
        occ_a[orb] = {true, IndexType(i)};
      }

      std::vector<SpVType_shm_csr_matrix> V2;
      V2.reserve(ndet);
      V2.emplace_back(generateHijkl(type, addCoulomb, TGwfn, occ_a, occ_a, cutv2));
      std::vector<SpVType_shm_csr_matrix::template matrix_view<int>> V2view;
      V2view.reserve(ndet);
      V2view.emplace_back(csr::shm::local_balanced_partition(V2[0], TGwfn));

      if (write_hdf)
        writeSparseTensor(dump, type, NMO, NAEA, NAEB, TGprop, TGwfn, H1, V2, Spvn, vn0, E0, global_ncvecs, 12);

      return HamiltonianOperations(sparse_ham(TGwfn.TG_local(), type, std::move(H1), std::move(hij), std::move(V2),
                                              std::move(V2view), std::move(Spvn), std::move(Spvnview), std::move(vn0),
                                              std::move(SpvnT), std::move(SpvnTview), E0, global_ncvecs));
    }
    else
    {
      // in this case: V2 is ComplexType since it is half rotated
      using sparse_ham = SparseTensor<ComplexType, ComplexType>;

      std::vector<SpCType_shm_csr_matrix> V2;
      V2.reserve(ndet);
      V2.emplace_back(halfRotatedHijkl(type, addCoulomb, TGwfn, &PsiTsp[0],
                                       ((type == COLLINEAR) ? (&PsiTsp[1]) : (&PsiTsp[0])), cutv2));
      std::vector<SpCType_shm_csr_matrix::template matrix_view<int>> V2view;
      V2view.reserve(ndet);
      V2view.emplace_back(csr::shm::local_balanced_partition(V2[0], TGwfn));

      if (write_hdf)
        writeSparseTensor(dump, type, NMO, NAEA, NAEB, TGprop, TGwfn, H1, V2, Spvn, vn0, E0, global_ncvecs, 22);

      return HamiltonianOperations(sparse_ham(TGwfn.TG_local(), type, std::move(H1), std::move(hij), std::move(V2),
                                              std::move(V2view), std::move(Spvn), std::move(Spvnview), std::move(vn0),
                                              std::move(SpvnT), std::move(SpvnTview), E0, global_ncvecs));
    }
  }
  else if (addCoulomb)
  {
    // in multi determinant with addCoulomb, SpvnT is transposed(Spvn)

    std::vector<boost::multi::array<ComplexType, 1>> hij;
    hij.reserve(ndet);
    int skp = ((type == COLLINEAR) ? 1 : 0);
    for (int n = 0, nd = 0; n < ndet; ++n, nd += (skp + 1))
      hij.emplace_back(halfRotatedHij(type, &PsiT[nd], &PsiT[nd + skp]));
    std::vector<SpVType_shm_csr_matrix> SpvnT;
    using matrix_view = typename SpVType_shm_csr_matrix::template matrix_view<int>;
    std::vector<matrix_view> SpvnTview;
    SpvnT.emplace_back(csr::shm::transpose<SpVType_shm_csr_matrix>(Spvn));
    SpvnTview.emplace_back(csr::shm::local_balanced_partition(SpvnT[0], TGprop));

    if (pureSD)
    {
      // in pureSD: V2 is ValueType
      using sparse_ham = SparseTensor<ValueType, ValueType>;

      return HamiltonianOperations{};
      //return HamiltonianOperations(sparse_ham(std::move(H1),std::move(hij),std::move(V2),
      //    std::move(V2view),std::move(Spvn),std::move(Spvnview),
      //    std::move(vn0),std::move(SpvnT),std::move(SpvnTview),E0,global_ncvecs));
    }
    else
    {
      // in this case: V2 is ComplexType since it is half rotated
      using sparse_ham = SparseTensor<ComplexType, ValueType>;

      std::vector<SpCType_shm_csr_matrix> V2;
      V2.reserve(ndet);
      for (int n = 0, nd = 0; n < ndet; ++n, nd += (skp + 1))
        V2.emplace_back(halfRotatedHijkl(type, addCoulomb, TGwfn, &PsiTsp[nd], &PsiTsp[nd + skp], cutv2));
      std::vector<SpCType_shm_csr_matrix::template matrix_view<int>> V2view;
      V2view.reserve(ndet);
      for (auto& v : V2)
        V2view.emplace_back(csr::shm::local_balanced_partition(v, TGwfn));

      if (write_hdf)
        writeSparseTensor(dump, type, NMO, NAEA, NAEB, TGprop, TGwfn, H1, V2, Spvn, vn0, E0, global_ncvecs, 21);

      return HamiltonianOperations(sparse_ham(TGwfn.TG_local(), type, std::move(H1), std::move(hij), std::move(V2),
                                              std::move(V2view), std::move(Spvn), std::move(Spvnview), std::move(vn0),
                                              std::move(SpvnT), std::move(SpvnTview), E0, global_ncvecs));
    }
  }
  else
  { // multideterminant for PHMSD

    assert(type == CLOSED);
    std::vector<boost::multi::array<ComplexType, 1>> hij;
    hij.reserve(ndet);
    for (int n = 0, nd = 0; n < ndet; ++n, nd++)
      hij.emplace_back(halfRotatedHij(type, &PsiT[nd], &PsiT[nd]));
    std::vector<SpCType_shm_csr_matrix> SpvnT;
    using matrix_view = typename SpCType_shm_csr_matrix::template matrix_view<int>;
    std::vector<matrix_view> SpvnTview;
    SpvnT.reserve(ndet);
    for (int n = 0; n < ndet; ++n)
    {
      SpvnT.emplace_back(sparse_rotate::halfRotateCholeskyMatrixForBias<SPComplexType>(type, TGprop, &PsiTsp[n],
                                                                                       &PsiTsp[n], Spvn, cutv2));
      SpvnTview.emplace_back(csr::shm::local_balanced_partition(SpvnT[n], TGprop));
    }

    if (pureSD)
    {
      // in pureSD: V2 is ValueType
      using sparse_ham = SparseTensor<ValueType, ValueType>;

      return HamiltonianOperations{};
      //return HamiltonianOperations(sparse_ham(std::move(H1),std::move(hij),std::move(V2),
      //    std::move(V2view),std::move(Spvn),std::move(Spvnview),
      //    std::move(vn0),std::move(SpvnT),std::move(SpvnTview),E0,global_ncvecs));
    }
    else
    {
      // in this case: V2 is ComplexType since it is half rotated
      using sparse_ham = SparseTensor<ComplexType, ComplexType>;

      std::vector<SpCType_shm_csr_matrix> V2;
      V2.reserve(ndet);
      for (int n = 0, nd = 0; n < ndet; ++n, nd++)
        V2.emplace_back(halfRotatedHijkl(type, addCoulomb, TGwfn, &PsiTsp[nd], &PsiTsp[nd], cutv2));
      std::vector<SpCType_shm_csr_matrix::template matrix_view<int>> V2view;
      V2view.reserve(ndet);
      for (auto& v : V2)
        V2view.emplace_back(csr::shm::local_balanced_partition(v, TGwfn));

      if (write_hdf)
        writeSparseTensor(dump, type, NMO, NAEA, NAEB, TGprop, TGwfn, H1, V2, Spvn, vn0, E0, global_ncvecs, 21);

      return HamiltonianOperations(sparse_ham(TGwfn.TG_local(), type, std::move(H1), std::move(hij), std::move(V2),
                                              std::move(V2view), std::move(Spvn), std::move(Spvnview), std::move(vn0),
                                              std::move(SpvnT), std::move(SpvnTview), E0, global_ncvecs));
    }
  }
}

} // namespace afqmc
} // namespace qmcplusplus
