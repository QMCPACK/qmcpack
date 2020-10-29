#include <cstdlib>
#include <memory>
#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <numeric>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"

#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Utilities/readHeader.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Matrix/hdf5_readers_new.hpp"
#include "AFQMC/Matrix/array_partition.hpp"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/coo_matrix.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// it is possible to subdivide the rows within equivalent nodes and communicate at the end
FactorizedSparseHamiltonian::shm_csr_matrix read_V2fact(hdf_archive& dump,
                                                        TaskGroup_& TG,
                                                        int nread,
                                                        int NMO,
                                                        int nvecs,
                                                        double cutoff1bar,
                                                        int int_blocks)
{
  using counter = qmcplusplus::afqmc::sparse_matrix_element_counter;
  using Alloc   = shared_allocator<SPValueType>;
  using ucsr_matrix =
      ma::sparse::ucsr_matrix<SPValueType, int, std::size_t, shared_allocator<SPValueType>, ma::sparse::is_root>;

  int min_i = 0;
  int max_i = nvecs;

  int nrows           = NMO * NMO;
  bool distribute_Ham = (TG.getNGroupsPerTG() < TG.getTotalNodes());
  std::vector<IndexType> row_counts(nrows);

  app_log() << " Reading sparse factorized two-electron integrals." << std::endl;
  // calculate column range that belong to this node
  if (distribute_Ham)
  {
    // count number of non-zero elements along each column (Chol Vec)
    std::vector<IndexType> col_counts(nvecs);
    csr_hdf5::multiple_reader_global_count(dump, counter(false, nrows, nvecs, 0, nrows, 0, nvecs, cutoff1bar),
                                           col_counts, TG, nread);

    std::vector<IndexType> sets(TG.getNumberOfTGs() + 1);
    simple_matrix_partition<TaskGroup_, IndexType, RealType> split(nrows, nvecs, cutoff1bar);
    if (TG.getCoreID() < nread)
      split.partition_over_TGs(TG, false, col_counts, sets);

    if (TG.getGlobalRank() == 0)
    {
      app_log() << " Partitioning of (factorized) Hamiltonian Vectors: ";
      for (int i = 0; i <= TG.getNumberOfTGs(); i++)
        app_log() << sets[i] << " ";
      app_log() << std::endl;
      app_log() << " Number of terms in each partitioning: ";
      for (int i = 0; i < TG.getNumberOfTGs(); i++)
        app_log() << accumulate(col_counts.begin() + sets[i], col_counts.begin() + sets[i + 1], 0) << " ";
      app_log() << std::endl;
    }

    TG.Node().broadcast(sets.begin(), sets.end());
    min_i = sets[TG.getTGNumber()];
    max_i = sets[TG.getTGNumber() + 1];

    csr_hdf5::multiple_reader_local_count(dump, counter(true, nrows, nvecs, 0, nrows, min_i, max_i, cutoff1bar),
                                          row_counts, TG, nread);
  }
  else
  {
    // should be faster if ham is not distributed
    csr_hdf5::multiple_reader_global_count(dump, counter(true, nrows, nvecs, 0, nrows, 0, nvecs, cutoff1bar),
                                           row_counts, TG, nread);
  }

  ucsr_matrix ucsr(tp_ul_ul{nrows, max_i - min_i}, tp_ul_ul{0, min_i}, row_counts, Alloc(TG.Node()));
  csr::matrix_emplace_wrapper<ucsr_matrix> csr_wrapper(ucsr, TG.Node());

  using mat_map = qmcplusplus::afqmc::matrix_map;
  csr_hdf5::multiple_reader_hdf5_csr<SPValueType, int>(csr_wrapper,
                                                       mat_map(false, true, nrows, nvecs, 0, nrows, min_i, max_i,
                                                               cutoff1bar),
                                                       dump, TG, nread);
  csr_wrapper.push_buffer();
  TG.node_barrier();

  return FactorizedSparseHamiltonian::shm_csr_matrix(std::move(ucsr));
}

} // namespace afqmc

} // namespace qmcplusplus
