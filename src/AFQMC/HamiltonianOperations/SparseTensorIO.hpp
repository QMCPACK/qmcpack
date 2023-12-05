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

#ifndef QMCPLUSPLUS_AFQMC_SPARSETENSORIO_HPP
#define QMCPLUSPLUS_AFQMC_SPARSETENSORIO_HPP

#include <fstream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"

#include "AFQMC/HamiltonianOperations/SparseTensor.hpp"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"

namespace qmcplusplus
{
namespace afqmc
{
template<typename T1, typename T2>
SparseTensor<T1, T2> loadSparseTensor(hdf_archive& dump,
                                      WALKER_TYPES type,
                                      int NMO,
                                      int NAEA,
                                      int NAEB,
                                      std::vector<PsiT_Matrix>& PsiT,
                                      TaskGroup_& TGprop,
                                      TaskGroup_& TGwfn,
                                      RealType cutvn,
                                      RealType cutv2)
{
#if defined(MIXED_PRECISION)
  using SpT1 = typename to_single_precision<T1>::value_type;
  using SpT2 = typename to_single_precision<T2>::value_type;
#else
  using SpT1 = T1;
  using SpT2 = T2;
#endif

  // NEEDS TO BE FIXED FOR SP CASE
  using T1_shm_csr_matrix = ma::sparse::csr_matrix<SpT1, int, std::size_t, shared_allocator<SpT1>, ma::sparse::is_root>;
  using T2_shm_csr_matrix = ma::sparse::csr_matrix<SpT2, int, std::size_t, shared_allocator<SpT2>, ma::sparse::is_root>;

  std::vector<int> dims(10);
  int ndet = (type == COLLINEAR ? PsiT.size() / 2 : PsiT.size());
  ValueType E0;
  int global_ncvecs = 0;
  int V2_nrows, V2_ncols, Spvn_nrows, Spvn_ncols;

  // read from HDF
  dump.push("HamiltonianOperations", false);
  dump.push("SparseTensor", false);

  if (TGwfn.Global().root())
  {
    if (!dump.readEntry(dims, "dims"))
    {
      app_error() << " Error in loadSparseTensor: Problems reading dataset. \n";
      APP_ABORT("");
    }
    if (dims[0] != NMO)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: NMO. \n";
      APP_ABORT("");
    }
    if (dims[1] != NAEA)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: NAEA. \n";
      APP_ABORT("");
    }
    if (dims[2] != NAEB)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: NAEB. \n";
      APP_ABORT("");
    }
    if (dims[3] != ndet)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: ndet. \n";
      APP_ABORT("");
    }
    if (type == CLOSED && dims[4] != 1)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    if (type == COLLINEAR && dims[4] != 2)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    if (type == NONCOLLINEAR && dims[4] != 3)
    {
      app_error() << " Error in loadSparseTensor: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    std::vector<ValueType> et;
    if (!dump.readEntry(et, "E0"))
    {
      app_error() << " Error in loadSparseTensor: Problems reading dataset. \n";
      APP_ABORT("");
    }
    E0 = et[0];
  }
  TGwfn.Global().broadcast_n(dims.data(), 9);
  TGwfn.Global().broadcast_value(E0);
  V2_nrows   = dims[5];
  V2_ncols   = dims[6];
  Spvn_nrows = dims[7];
  Spvn_ncols = dims[8];

  // read 1-body hamiltonian and exchange potential (v0)
  boost::multi::array<ComplexType, 2> H1({NMO, NMO});
  boost::multi::array<ComplexType, 2> v0({NMO, NMO});
  if (TGwfn.Global().root())
  {
    boost::multi::array<ValueType, 2> H1_({NMO, NMO});
    if (!dump.readEntry(H1_, "H1"))
    {
      app_error() << " Error in loadSparseTensor: Problems reading dataset. \n";
      APP_ABORT("");
    }
    copy_n_cast(H1_.origin(), NMO * NMO, to_address(H1.origin()));
    if (!dump.readEntry(v0, "v0"))
    {
      app_error() << " Error in loadSparseTensor: Problems reading dataset. \n";
      APP_ABORT("");
    }
  }
  TGwfn.Global().broadcast_n(H1.origin(), H1.num_elements());
  TGwfn.Global().broadcast_n(v0.origin(), v0.num_elements());

  // read half-rotated exchange matrix
  std::vector<T1_shm_csr_matrix> V2;
  V2.reserve(ndet);
  for (int i = 0; i < ndet; i++)
  {
    dump.push(std::string("HalfRotatedHijkl_") + std::to_string(i), false);
    V2.emplace_back(
        csr_hdf5::unstructured_distributed_CSR_from_HDF<T1_shm_csr_matrix>(dump, TGwfn, V2_nrows, V2_ncols, cutv2));
    dump.pop();
  }

  // read Cholesky matrix
  dump.push(std::string("SparseCholeskyMatrix"), false);
  SpVType_shm_csr_matrix Spvn(
      csr_hdf5::column_distributed_CSR_from_HDF<SpVType_shm_csr_matrix>(dump, TGprop, Spvn_nrows, Spvn_ncols, cutv2));
  dump.pop();
  // leave hdf_archive in group it started from
  dump.pop();
  dump.pop();

  // rotated 1 body hamiltonians
  std::vector<boost::multi::array<ComplexType, 1>> hij;
  hij.reserve(ndet);
  int skp = ((type == COLLINEAR) ? 1 : 0);
  for (int n = 0, nd = 0; n < ndet; ++n, nd += (skp + 1))
  {
    check_wavefunction_consistency(type, &PsiT[nd], &PsiT[nd + skp], NMO, NAEA, NAEB);
    hij.emplace_back(rotateHij(type, &PsiT[nd], &PsiT[nd + skp], H1));
  }

  // setup views
  using matrix_view = typename T1_shm_csr_matrix::template matrix_view<int>;
  std::vector<matrix_view> V2view;
  V2view.reserve(ndet);
  for (auto& v : V2)
    V2view.emplace_back(csr::shm::local_balanced_partition(v, TGwfn));

  auto Spvnview(csr::shm::local_balanced_partition(Spvn, TGprop));

  if (ndet == 1)
  {
    std::vector<T2_shm_csr_matrix> SpvnT;
    // MAM: chech that T2 is SpComplexType
    //std::vector<SpCType_shm_csr_matrix> SpvnT;
    using matrix_view = typename T2_shm_csr_matrix::template matrix_view<int>;
    std::vector<matrix_view> SpvnTview;
#if defined(MIXED_PRECISION)
    auto PsiTsp(csr::shm::CSRvector_to_single_precision<PsiT_Matrix_t<SPComplexType>>(PsiT));
#else
    auto& PsiTsp(PsiT);
#endif
    SpvnT.emplace_back(
        sparse_rotate::halfRotateCholeskyMatrixForBias<T2>(type, TGprop, &PsiTsp[0],
                                                           ((type == COLLINEAR) ? (&PsiTsp[1]) : (&PsiTsp[0])), Spvn,
                                                           cutv2));
    SpvnTview.emplace_back(csr::shm::local_balanced_partition(SpvnT[0], TGprop));

    return SparseTensor<T1, T2>(TGwfn.TG_local(), type, std::move(H1), std::move(hij), std::move(V2), std::move(V2view),
                                std::move(Spvn), std::move(Spvnview), std::move(v0), std::move(SpvnT),
                                std::move(SpvnTview), E0, Spvn_ncols);
  }
  else
  {
    // problem here!!!!!
    // don't know how to tell if this is NOMSD or PHMSD!!!
    // whether to rotate or not! That's the question!
    std::vector<T2_shm_csr_matrix> SpvnT;
    //std::vector<SpVType_shm_csr_matrix> SpvnT;
    using matrix_view = typename T2_shm_csr_matrix::template matrix_view<int>;
    std::vector<matrix_view> SpvnTview;
    SpvnT.emplace_back(csr::shm::transpose<T2_shm_csr_matrix>(Spvn));
    SpvnTview.emplace_back(csr::shm::local_balanced_partition(SpvnT[0], TGprop));

    return SparseTensor<T1, T2>(TGwfn.TG_local(), type, std::move(H1), std::move(hij), std::move(V2), std::move(V2view),
                                std::move(Spvn), std::move(Spvnview), std::move(v0), std::move(SpvnT),
                                std::move(SpvnTview), E0, global_ncvecs);
  }
}

// single writer right now
template<class shm_mat1, class shm_mat2>
inline void writeSparseTensor(hdf_archive& dump,
                              WALKER_TYPES type,
                              int NMO,
                              int NAEA,
                              int NAEB,
                              TaskGroup_& TGprop,
                              TaskGroup_& TGwfn,
                              boost::multi::array<ValueType, 2>& H1,
                              std::vector<shm_mat1> const& v2,
                              shm_mat2 const& Spvn,
                              boost::multi::array<ComplexType, 2>& v0,
                              ValueType E0,
                              int gncv,
                              int code)
{
  assert(v2.size() > 0);
  if (TGwfn.Global().root())
  {
    dump.push("HamiltonianOperations");
    dump.push("SparseTensor");
    std::vector<int>
        dims{NMO, NAEA, NAEB, int(v2.size()), type, int(v2[0].size(0)), int(v2[0].size(1)), int(Spvn.size(0)), gncv};
    dump.write(dims, "dims");
    std::vector<int> dm{code};
    dump.write(dm, "type");
    std::vector<ValueType> et{E0};
    dump.write(et, "E0");
    dump.write(H1, "H1");
    dump.write(v0, "v0");
  }

  for (int i = 0; i < v2.size(); i++)
  {
    if (TGwfn.Global().root())
      dump.push(std::string("HalfRotatedHijkl_") + std::to_string(i));
    csr_hdf5::write_distributed_CSR_to_HDF(v2[i], dump, TGwfn);
    if (TGwfn.Global().root())
      dump.pop();
  }

  if (TGprop.Global().root())
    dump.push("SparseCholeskyMatrix");
  csr_hdf5::write_distributed_CSR_to_HDF(Spvn, dump, TGprop);

  if (TGwfn.Global().root())
  {
    dump.pop();
    dump.pop();
    dump.pop();
  }
  TGwfn.Global().barrier();
}

} // namespace afqmc
} // namespace qmcplusplus

#endif
