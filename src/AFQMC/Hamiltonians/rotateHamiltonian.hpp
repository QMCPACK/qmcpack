//////////////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HPP
#define QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HPP

#include <vector>
#include <tuple>
#include <mpi.h>
#include <algorithm>
#include <numeric>

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"
#include "AFQMC/Utilities/afqmc_TTI.hpp"
#include "AFQMC/Utilities/myTimer.h"

#include "AFQMC/Hamiltonians/rotateHamiltonian_Helper2.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// due to generalized Slater matrices, the conditions change
template<class PsiT_Type>
inline void check_wavefunction_consistency(WALKER_TYPES type, PsiT_Type* A, PsiT_Type* B, int NMO, int NAEA, int NAEB)
{
  if (type == CLOSED)
  {
    if (A->size(1) != NMO || A->size(0) < NAEA)
    {
      app_error() << " Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency(): wfn_type=0, NMO, "
                     "NAEA, A.rows, A.cols: "
                  << NMO << " " << NAEA << " " << A->size(0) << " " << A->size(1) << std::endl;
      APP_ABORT(" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency().\n");
    }
  }
  else if (type == COLLINEAR)
  {
    if (A->size(1) != NMO || A->size(0) < NAEA || B->size(1) != NMO || B->size(0) < NAEB)
    {
      app_error() << " Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency(): wfn_type=1, NMO, "
                     "NAEA, NAEB, A.rows, A.cols, B.rows, B.cols: "
                  << NMO << " " << NAEA << " " << NAEB << " " << A->size(0) << " " << A->size(1) << " " << B->size(0)
                  << " " << B->size(1) << std::endl;
      APP_ABORT(" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency().\n");
    }
  }
  else if (type == NONCOLLINEAR)
  {
    if (A->size(1) != 2 * NMO || A->size(0) < (NAEB + NAEA))
    {
      app_error() << " Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency(): wfn_type=1, NMO, "
                     "NAEA, NAEB, A.rows, A.cols: "
                  << NMO << " " << NAEA << " " << NAEB << " " << A->size(0) << " " << A->size(1) << std::endl;
      APP_ABORT(" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency().\n");
    }
  }
  else
  {
    app_error() << " Error: Unacceptable walker_type in check_wavefunction_consistency(): " << type << std::endl;
    APP_ABORT(" Error: Unacceptable walker_type in check_wavefunction_consistency(). \n");
  }
}

template<class Mat>
inline boost::multi::array<ComplexType, 1> rotateHij(WALKER_TYPES walker_type,
                                                     PsiT_Matrix* Alpha,
                                                     PsiT_Matrix* Beta,
                                                     Mat const& H1)
{
  assert(Alpha != nullptr);
  int NAEA = Alpha->size(0);
  int NMO  = Alpha->size(1);

  boost::multi::array<ComplexType, 1> N(iextensions<1u>{1});
  const ComplexType one  = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0);

  // 1-body part
  if (walker_type == CLOSED || walker_type == NONCOLLINEAR)
  {
    N.reextent(iextensions<1u>{NAEA * NMO});
    boost::multi::array_ref<ComplexType, 2> N_(N.origin(), {NAEA, NMO});

    ma::product(*Alpha, H1, N_);
    ma::scal(ComplexType(2.0), N);
  }
  else if (walker_type == COLLINEAR)
  {
    assert(Beta != nullptr);
    int NAEB = Beta->size(0);

    N.reextent(iextensions<1u>{(NAEA + NAEB) * NMO});
    boost::multi::array_ref<ComplexType, 2> NA_(N.origin(), {NAEA, NMO});
    boost::multi::array_ref<ComplexType, 2> NB_(N.origin() + NAEA * NMO, {NAEB, NMO});

    ma::product(*Alpha, H1, NA_);
    ma::product(*Beta, H1, NB_);
  }

  return N;
}

template<class Container = std::vector<std::tuple<int, int, SPComplexType>>,
         class PsiT_Type = PsiT_Matrix_t<SPComplexType>>
inline void rotateHijkl(std::string& type,
                        WALKER_TYPES walker_type,
                        bool addCoulomb,
                        TaskGroup_& TG,
                        Container& Vijkl,
                        PsiT_Type* Alpha,
                        PsiT_Type* Beta,
                        SpVType_shm_csr_matrix const& V2_fact,
                        const RealType cut,
                        int maximum_buffer_size,
                        bool reserve_to_fit_ = true,
                        bool global_reserve  = true)
{
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool distribute_Ham = TG.getNumberOfTGs() > 1;
  if (distribute_Ham)
    APP_ABORT(" Distributed V2_fact not yet implemented. \n");

  int NAEA = Alpha->size(0);
  int NMO  = Alpha->size(1);
  int NAEB = NAEA;
  if (walker_type == COLLINEAR)
    NAEB = Beta->size(0);

  // <ab||kl> = sum_n Qk(a,n) * Rl(b,n) - Rl(a,n)*Qk(b,n),
  // where:
  //   Qk(a,n) = sum_i ma::conj(Amat(i,a)) * V2_fact(ik,n)
  //   Rl(a,n) = sum_i ma::conj(Amat(i,a)) * ma::conj(V2_fact(li,n))
  // For real build, Qk=Rk
  //
  // For parallelization, distribute (k,l) pairs over nodes.
  // Build ahead of time Qk/Rl matrices in shared memory to reduce memory/setup time.
  // Assemble integrals in parallel and fully distributed.
  // Collect on all nodes.
  //    - For distributed hamiltonians, you do not need to worry about keeping contiguous
  //    segments of the hamiltonian. Only that the distribution over nodes is roughly equal.
  //    Write a simple algorithm that balances the number of terms in a TG.
  //

  bool sparseQk = (type == "SD" || type == "SS");
  bool sparseRl = (type == "SS");

  app_log() << " Calculating half-rotated Hamiltonian using ";
  if (type == "SD")
    app_log() << "Sparse x Dense";
  else if (type == "SS")
    app_log() << "Sparse x Sparse";
  else if (type == "DD")
    app_log() << "Dense x Dense";
  app_log() << " matrix multiplication. \n";

  using shmSpVector = SPComplexVector<shared_allocator<SPComplexType>>;
  using shmSpMatrix = SPComplexMatrix<shared_allocator<SPComplexType>>;

  int NMO2 = (walker_type == CLOSED) ? NMO : 2 * NMO;
  int ngrp = std::min(NMO2, nnodes);
  std::vector<int> M_split(ngrp + 1);

  // split the orbitals among processors:  ngrp/2 + ngrp%2 for alpha, ngrp/2 for beta
  M_split[0] = 0;
  if (walker_type == CLOSED)
  {
    FairDivide(NMO2, ngrp, M_split);
  }
  else if (walker_type == COLLINEAR)
  {
    if (ngrp == 1)
      APP_ABORT(" Error: Current implementation of rotateHijkl requires at least 2 nodes.\n");
    // This should be 2/3-1/3 partitioning between alpha/beta
    // to balance the resulting matrix elements better
    int nalpha = ngrp / 2;
    for (int i = 0; i < nalpha; i++)
    {
      int m0, m1;
      std::tie(m0, m1) = FairDivideBoundary(i, NMO, nalpha);
      M_split[i + 1]   = m1;
    }
    assert(M_split[nalpha] == NMO);
    for (int i = 0, nbeta = ngrp / 2 + ngrp % 2; i < nbeta; i++)
    {
      int m0, m1;
      std::tie(m0, m1)        = FairDivideBoundary(i, NMO, nbeta);
      M_split[i + nalpha + 1] = NMO + m1;
    }
    assert(M_split[ngrp] == 2 * NMO);
  }
  else if (walker_type == NONCOLLINEAR)
  {
    APP_ABORT(" Finish. \n");
  }

  // Construct your set of Q(k,a,m), R(l,b,m)
  int l0        = (nodeid < ngrp) ? (M_split[nodeid]) : (-1);
  int lN        = (nodeid < ngrp) ? (M_split[nodeid + 1]) : (-1);
  bool amIAlpha = true; // for simplicity, subset of bands must be either elpha or beta, not mixed
  if (l0 < NMO && (lN - 1) < NMO)
    amIAlpha = true;
  else if (l0 >= NMO && lN >= NMO)
    amIAlpha = false;
  else
  {
    std::cerr << "l0, lN, nodeid, ngrp, NMO: " << l0 << " " << lN << " " << nodeid << " " << ngrp << " " << NMO
              << std::endl;
    std::cerr << " Error: Current algorithm requires an even number of processors. \n\n\n";
    APP_ABORT(" Error: Current algorithm requires an even number of processors. \n\n\n");
  }

  // use new communicator in case addCoulomb is false.
  int key;
  if (l0 < 0)
    key = 0;
  if (addCoulomb)
    key = 1; // single communicator with all working nodes
  else
    key = (amIAlpha ? 1 : 2); // 2 communicators, one for each spin (if there are 2)
  boost::mpi3::communicator comm(TG.Cores().split(key));
  //boost::mpi3::communicator comm(TG.Cores().split(0));

  int norb    = lN - l0;
  int maxnorb = 0;
  for (int i = 0; i < ngrp; i++)
    maxnorb = std::max(maxnorb, M_split[i + 1] - M_split[i]);
  int nvec = V2_fact.size(1);
  // must gather over heads of TG to get nchol per TG and total # chol vecs
  // Rl(k, a, m), k:[0:NMO2], a:[0:NAEA], m:[0:nvec]

  app_log() << " Approximate memory usage for half-rotated Hamiltonian construction: \n"
            << "   max. number of orbital in a node: " << maxnorb << "\n"
            << "   Qk/Rl matrices size (assuming dense) each = maxnorb * nup * ncholvecs complex numbers = "
            << sizeof(SPComplexType) * maxnorb * NAEA * nvec / 1024.0 / 1024.0 << " MB \n"
            << "   Maximum size of communication buffer: " << maximum_buffer_size << " MB" << std::endl;

  const int nrow = norb * ((amIAlpha) ? NAEA : NAEB);
  const int ncol = nvec;
  int dummy_nrow = nrow, dummy_ncol = ncol;
  if (nodeid >= ngrp || sparseQk)
    dummy_nrow = dummy_ncol = 0;

  using Alloc = shared_allocator<SPComplexType>;
  // global_origin is not set correctly, careful not to rely on it
  SpCType_shm_csr_matrix SpQk(tp_ul_ul{nrow, ncol}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
  SpCType_shm_csr_matrix SpRl(tp_ul_ul{ncol, nrow}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));

  shmSpMatrix Qk({dummy_nrow, dummy_ncol}, shared_allocator<SPComplexType>{TG.Node()});
  if (coreid == 0)
    std::fill_n(Qk.origin(), Qk.num_elements(), SPComplexType(0.0));
  dummy_nrow = nrow;
  dummy_ncol = ncol;
  if (sparseRl or nodeid >= ngrp)
    dummy_nrow = dummy_ncol = 0;
  shmSpMatrix Rl({dummy_ncol, dummy_nrow}, shared_allocator<SPComplexType>{TG.Node()});
  if (coreid == 0)
    std::fill_n(Rl.origin(), Rl.num_elements(), SPComplexType(0.0));

  if (distribute_Ham)
  {
    APP_ABORT(" Finish THIS (43)!!! \n\n\n");
  }
  else
  {
    //   Q(k,a,n) = sum_i ma::conj(Amat(i,a)) * V2_fact(ik,n)
    //   R(l,a,n) = sum_i ma::conj(Amat(i,a)) * ma::conj(V2_fact(li,n))

    // Construct SpQk[k,n,nvec]
    if (sparseQk)
    {
      sparse_rotate::halfRotateCholeskyMatrix(walker_type, TG, l0, lN, SpQk, Alpha, Beta, V2_fact, false, false, cut,
                                              true);
      SpQk.remove_empty_spaces(); // just in case
    }
    else
      ma_rotate::halfRotateCholeskyMatrix(walker_type, TG, l0, lN, Qk, Alpha, Beta, V2_fact, false, false, cut);

#if defined(QMC_COMPLEX)
    // Construct SpRl[nvec,k,a]
    if (sparseRl)
    {
      // since Rl is transposed, I can't emplace_back on the csr matrix.
      // In this case, I need to use a temporary ucsr with an emplace_wrapper
      using ucsr_matrix = ma::sparse::ucsr_matrix<SPComplexType, int, std::size_t, shared_allocator<SPComplexType>,
                                                  ma::sparse::is_root>;
      ucsr_matrix ucsr(tp_ul_ul{ncol, nrow}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
      csr::matrix_emplace_wrapper<ucsr_matrix> ucsr_wrapper(ucsr, TG.Node());
      sparse_rotate::halfRotateCholeskyMatrix(walker_type, TG, l0, lN, ucsr_wrapper, Alpha, Beta, V2_fact, true, true,
                                              cut, true);
      ucsr_wrapper.push_buffer(); // push any remaining elements in temporary buffer
      SpRl = std::move(ucsr);
      SpRl.remove_empty_spaces(); // just in case
    }
    else
      ma_rotate::halfRotateCholeskyMatrix(walker_type, TG, l0, lN, Rl, Alpha, Beta, V2_fact, true, true, cut);
#else
    if (sparseRl)
    {
      if (sparseQk)
      {
        SpRl = csr::shm::transpose<SpCType_shm_csr_matrix>(SpQk);
        SpRl.remove_empty_spaces(); // just in case
      }
      else
      {
        app_error() << " Error: Incorrect matrix setup in createHamiltonianForGeneralDeterminant. sparseRl=True, "
                       "sparseQk=False."
                    << std::endl;
        APP_ABORT(" Error: Incorrect matrix setup in createHamiltonianForGeneralDeterminant. sparseRl=True, "
                  "sparseQk=False. \n");
      }
    }
    else
    {
      if (sparseQk)
      {
        csr::shm::transpose(SpQk, Rl);
      }
      else
      {
        // Qk[norb*NAEA,nvec]
        // Rl[nvec,norb*NAEA]
        int n0_, n1_, sz_ = Qk.size();
        std::tie(n0_, n1_) = FairDivideBoundary(coreid, sz_, ncores);
        if (n1_ - n0_ > 0)
          ma::transpose(Qk.sliced(n0_, n1_), Rl(Rl.extension(0), {n0_, n1_}));
      }
    }
#endif
  }

  TG.node_barrier();
  // let the maximum message be maximum_buffer_size MB
  int maxnt = std::max(1, static_cast<int>(std::floor(maximum_buffer_size * 1024.0 * 1024.0 / sizeof(SPComplexType))));

  // control the size of the MPI messages.
  // communicate blocks of Qk up to maxnt terms
  std::vector<int> nkbounds;           // local bounds for communication
  std::vector<int> Qknum(comm.size()); // number of blocks per node
  std::vector<int> Qksizes;            // number of terms and number of k vectors in block for all nodes
  if (coreid == 0)
  {
    int n_ = 0, ntcnt = 0, n0 = 0;
    for (int i = 0; i < norb; i++)
    {
      int ntt;
      if (sparseQk)
        ntt = SpQk.num_non_zero_elements(i);
      else
        ntt = nvec * NAEA;
      assert(ntt < maxnt);
      if (ntcnt + ntt > maxnt)
      {
        nkbounds.push_back(ntcnt);
        nkbounds.push_back(i - n0);
        ntcnt = ntt;
        n0    = i;
        n_++;
      }
      else
      {
        ntcnt += ntt;
      }
    }
    if (ntcnt > 0)
    {
      // push last block
      n_++;
      nkbounds.push_back(ntcnt);
      nkbounds.push_back(norb - n0);
    }

    MPI_Allgather(&n_, 1, MPI_INT, Qknum.data(), 1, MPI_INT, comm.get());

    int ntt = std::accumulate(Qknum.begin(), Qknum.end(), 0);
    Qksizes.resize(2 * ntt);

    std::vector<int> cnts(comm.size());
    std::vector<int> disp(comm.size());
    int cnt = 0;
    for (int i = 0; i < comm.size(); i++)
    {
      cnts[i] = Qknum[i] * 2;
      disp[i] = cnt;
      cnt += cnts[i];
    }
    MPI_Allgatherv(nkbounds.data(), nkbounds.size(), MPI_INT, Qksizes.data(), cnts.data(), disp.data(), MPI_INT,
                   comm.get());
  }

  MPI_Bcast(Qknum.data(), comm.size(), MPI_INT, 0, TG.Node().get());
  int ntt = std::accumulate(Qknum.begin(), Qknum.end(), 0);
  if (coreid != 0)
    Qksizes.resize(2 * ntt);
  MPI_Bcast(Qksizes.data(), Qksizes.size(), MPI_INT, 0, TG.Node().get());

  // store {nterms,nk} for all nodes
  // use it to know communication pattern

  int maxnk      = 0; // maximum number of k vectors in communication block
  long maxqksize = 0; // maximum size of communication block
  for (int i = 0; i < ntt; i++)
  {
    if (Qksizes[2 * i] > maxqksize)
      maxqksize = Qksizes[2 * i];
    if (Qksizes[2 * i + 1] > maxnk)
      maxnk = Qksizes[2 * i + 1];
  }

  app_log() << "   Maximum number of (k,l) pairs communicated in a block: " << maxnk << "\n"
            << "   Temporary integral matrix Ta: "
            << norb * NAEA * maxnk * NAEA * sizeof(SPComplexType) / 1024.0 / 1024.0 << " MB " << std::endl;

  // temporary shared memory space for local "dense" result
  shmSpVector Ta_shmbuff(iextensions<1u>{norb * NAEA * maxnk * NAEA}, shared_allocator<SPComplexType>{TG.Node()});

  // setup working sparse matrix
  dummy_nrow = maxnk * NAEA;
  dummy_ncol = nvec;
  shmSpVector tQk_shmbuff(iextensions<1u>{1}, shared_allocator<SPComplexType>{TG.Node()});
  SpCType_shm_csr_matrix SptQk(tp_ul_ul{maxnk * NAEA, nvec}, tp_ul_ul{0, 0}, 0, Alloc(TG.Node()));
  if (sparseQk)
  {
    std::size_t sz_ = std::ceil(maxqksize / SptQk.size(0));
    SptQk.reserve(sz_);
  }
  else
    tQk_shmbuff.reextent(iextensions<1u>{maxnk * NAEA * nvec});
  if (sparseQk)
    dummy_nrow = dummy_ncol = 0;

  myTimer Timer_;

  SPComplexType EJX(0.0);
  if (addCoulomb)
    EJX = SPComplexType(1.0);

  TG.Node().barrier();
  if (reserve_to_fit_)
  {
    // count and resize container
    std::vector<std::size_t> sz_local;
    if (walker_type == CLOSED)
      sz_local.resize(NMO * NAEA);
    else if (walker_type == COLLINEAR)
      sz_local.resize(NMO * (NAEA + NAEB));
    else if (walker_type == NONCOLLINEAR)
      sz_local.resize(2 * NMO * (NAEA + NAEB));

    int K0_ = (key == 2 ? NMO : 0);
    int nn0 = (key == 2 ? ngrp / 2 : 0);
    for (int nn = 0, nb = 0, nkcum = 0; nn < comm.size(); nn++)
    {
      // just checking
      assert(nkcum + K0_ == M_split[nn + nn0]);
      if (M_split[nn + nn0 + 1] == M_split[nn + nn0])
        continue;
      int nblk = Qknum[nn];
#ifndef NDEBUG
      long ntermscum = 0;
#endif
      for (int bi = 0; bi < nblk; bi++, nb++)
      {
        int nterms = Qksizes[2 * nb];     // number of terms in block
        int nk     = Qksizes[2 * nb + 1]; // number of k-blocks in block
        int k0     = nkcum + K0_;         // first value of k in block
        nkcum += nk;
        int kN   = nkcum + K0_;              // last+1 value
        int NEL0 = (k0 < NMO) ? NAEA : NAEB; // number of electrons in this spin block
        assert(nk > 0 && nk <= maxnk);       // just checking

        boost::multi::array_ref<SPComplexType, 2> tQk(to_address(tQk_shmbuff.origin()), {nk * NEL0, nvec});

        Timer_.reset("T0");
        Timer_.start("T0");
        if (sparseQk)
        {
          if (coreid == 0)
          {
            if (nn == comm.rank())
            {
              auto ka0 = (k0 - M_split[nn + nn0]) * NEL0;
              auto kaN = (k0 - M_split[nn + nn0] + nk) * NEL0;
              auto n0  = *SpQk.pointers_begin(ka0);
              auto n1  = *SpQk.pointers_end(kaN);
              int nt_  = static_cast<int>(n1 - n0);
              assert(ntermscum == n0);
              assert(nt_ == nterms);
              std::copy(to_address(SpQk.non_zero_values_data(ka0)), to_address(SpQk.non_zero_values_data(kaN)),
                        to_address(SptQk.non_zero_values_data()));
              std::copy(to_address(SpQk.non_zero_indices2_data(ka0)), to_address(SpQk.non_zero_indices2_data(kaN)),
                        to_address(SptQk.non_zero_indices2_data()));
              for (int i = 0, j = ka0; i <= nk * NEL0; i++, j++)
              {
                SptQk.pointers_begin()[i] = SpQk.pointers_begin()[j] - n0;
                SptQk.pointers_end()[i]   = SpQk.pointers_end()[j] - n0;
              }
            }
            comm.broadcast_value(nterms, nn);
            comm.broadcast_n(to_address(SptQk.pointers_begin()), SptQk.size(0), nn);
            comm.broadcast_n(to_address(SptQk.pointers_end()), SptQk.size(0), nn);
            comm.broadcast_n(to_address(SptQk.non_zero_values_data()), nterms, nn);
            comm.broadcast_n(to_address(SptQk.non_zero_indices2_data()), nterms, nn);
          }
          TG.node_barrier();
#ifndef NDEBUG
          // for safety, keep track of sum
          ntermscum += static_cast<long>(nterms);
#endif
        }
        else
        {
          if (coreid == 0)
          {
            if (nn == comm.rank())
              std::copy(Qk.origin() + bi * maxnk * NEL0 * nvec, Qk.origin() + (bi * maxnk + nk) * NEL0 * nvec,
                        tQk.origin());
            comm.broadcast_n(tQk.origin(), nk * NEL0 * nvec, nn);
          }
          TG.node_barrier();
        }
        Timer_.stop("T0");
        app_log() << " Loop: " << nn << "/" << comm.size() << " " << bi << "/" << nblk
                  << " communication: " << Timer_.total("T0") << " ";

        boost::multi::array_ref<SPComplexType, 2> Ta(to_address(Ta_shmbuff.origin()), {nk * NEL0, nrow});

        Timer_.reset("T0");
        Timer_.start("T0");
        if (type == "SD")
          count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0, kN, l0, lN, NMO, NAEA, NAEB,
                        SptQk[{0, std::size_t(nk * NEL0)}], Rl, Ta, cut);
        else if (type == "DD")
          count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0, kN, l0, lN, NMO, NAEA, NAEB, tQk, Rl, Ta, cut);
        Timer_.stop("T0");
        app_log() << " QxR: " << Timer_.total("T0") << std::endl;
        //      else if(type == "SS")
      }
    }

    std::size_t tot_sz_local = std::accumulate(sz_local.begin(), sz_local.end(), std::size_t(0));

    std::vector<std::size_t> sz_global(sz_local.size());
    TG.Global().all_reduce_n(sz_local.begin(), sz_local.size(), sz_global.begin(), std::plus<>());
    std::size_t tot_sz_global = std::accumulate(sz_global.begin(), sz_global.end(), std::size_t(0));

    std::vector<std::size_t> sz_node(sz_local.size());
    TG.Node().all_reduce_n(sz_local.begin(), sz_local.size(), sz_node.begin(), std::plus<>());
    std::size_t tot_sz_node = std::accumulate(sz_node.begin(), sz_node.end(), std::size_t(0));

    std::size_t sz_node_min, sz_node_max, sz_local_min, sz_local_max;
    TG.Global().all_reduce_n(&tot_sz_node, 1, &sz_node_min, boost::mpi3::min<>());
    TG.Global().all_reduce_n(&tot_sz_node, 1, &sz_node_max, boost::mpi3::max<>());
    TG.Global().all_reduce_n(&tot_sz_local, 1, &sz_local_min, boost::mpi3::min<>());
    TG.Global().all_reduce_n(&tot_sz_local, 1, &sz_local_max, boost::mpi3::max<>());

    app_log() << "  Number of terms in Vijkl: \n"
              << "    Local: (min/max)" << sz_local_min << " "
              << sz_local_min * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB  -  " << sz_local_max
              << " " << sz_local_max * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB \n"
              << "    Node (min/max): " << sz_node_min << " "
              << sz_node_min * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB   -   " << sz_node_max
              << " " << sz_node_max * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB \n"
              << "    Global: " << tot_sz_global << " "
              << tot_sz_global * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB" << std::endl
              << std::endl;

    if (global_reserve)
      reserve_to_fit(Vijkl, sz_global);
    else
      reserve_to_fit(Vijkl, sz_local);
  }

  // now calculate fully distributed matrix elements
  int K0_ = (key == 2 ? NMO : 0);
  int nn0 = (key == 2 ? ngrp / 2 : 0);
  for (int nn = 0, nb = 0, nkcum = 0; nn < comm.size(); nn++)
  {
    // just checking
    assert(nkcum + K0_ == M_split[nn + nn0]);
    if (M_split[nn + nn0 + 1] == M_split[nn + nn0])
      continue;
    int nblk = Qknum[nn];
#ifndef NDEBUG
    long ntermscum = 0;
#endif
    for (int bi = 0; bi < nblk; bi++, nb++)
    {
      int nterms = Qksizes[2 * nb];     // number of terms in block
      int nk     = Qksizes[2 * nb + 1]; // number of k-blocks in block
      int k0     = nkcum + K0_;         // first value of k in block
      nkcum += nk;
      int kN   = nkcum + K0_;              // last+1 value
      int NEL0 = (k0 < NMO) ? NAEA : NAEB; // number of electrons in this spin block
      assert(nk > 0 && nk <= maxnk);       // just checking

      boost::multi::array_ref<SPComplexType, 2> tQk(to_address(tQk_shmbuff.origin()), {nk * NEL0, nvec});

      Timer_.reset("T0");
      Timer_.start("T0");
      if (sparseQk)
      {
        if (coreid == 0)
        {
          if (nn == comm.rank())
          {
            auto ka0 = (k0 - M_split[nn + nn0]) * NEL0;
            auto kaN = (k0 - M_split[nn + nn0] + nk) * NEL0;
            auto n0  = *SpQk.pointers_begin(ka0);
            auto n1  = *SpQk.pointers_end(kaN);
            int nt_  = static_cast<int>(n1 - n0);
            assert(ntermscum == n0);
            assert(nt_ == nterms);
            std::copy(to_address(SpQk.non_zero_values_data(ka0)), to_address(SpQk.non_zero_values_data(kaN)),
                      to_address(SptQk.non_zero_values_data()));
            std::copy(to_address(SpQk.non_zero_indices2_data(ka0)), to_address(SpQk.non_zero_indices2_data(kaN)),
                      to_address(SptQk.non_zero_indices2_data()));
            for (int i = 0, j = ka0; i <= nk * NEL0; i++, j++)
            {
              SptQk.pointers_begin()[i] = SpQk.pointers_begin()[j] - n0;
              SptQk.pointers_end()[i]   = SpQk.pointers_end()[j] - n0;
            }
          }
          comm.broadcast_value(nterms, nn);
          comm.broadcast_n(to_address(SptQk.pointers_begin()), SptQk.size(0), nn);
          comm.broadcast_n(to_address(SptQk.pointers_end()), SptQk.size(0), nn);
          comm.broadcast_n(to_address(SptQk.non_zero_values_data()), nterms, nn);
          comm.broadcast_n(to_address(SptQk.non_zero_indices2_data()), nterms, nn);
        }
        TG.node_barrier();
#ifndef NDEBUG
        // for safety, keep track of sum
        ntermscum += static_cast<long>(nterms);
#endif
      }
      else
      {
        if (coreid == 0)
        {
          if (nn == comm.rank())
            std::copy(Qk.origin() + bi * maxnk * NEL0 * nvec, Qk.origin() + (bi * maxnk + nk) * NEL0 * nvec,
                      tQk.origin());
          comm.broadcast_n(tQk.origin(), nk * NEL0 * nvec, nn);
        }
        TG.node_barrier();
      }
      app_log() << " Loop: " << nn << "/" << comm.size() << " " << bi << "/" << nblk
                << " communication: " << Timer_.total("T0") << " ";

      boost::multi::array_ref<SPComplexType, 2> Ta(to_address(Ta_shmbuff.origin()), {nk * NEL0, nrow});

      Timer_.reset("T0");
      Timer_.start("T0");
      if (type == "SD")
        Qk_x_Rl(walker_type, EJX, TG, k0, kN, l0, lN, NMO, NAEA, NAEB, SptQk[{0, std::size_t(nk * NEL0)}], Rl, Ta,
                Vijkl, cut);
      else if (type == "DD")
        Qk_x_Rl(walker_type, EJX, TG, k0, kN, l0, lN, NMO, NAEA, NAEB, tQk, Rl, Ta, Vijkl, cut);
      Timer_.stop("T0");
      app_log() << " QxR: " << Timer_.total("T0") << std::endl;
      //      else if(type == "SS")
    }
  }

  TG.node_barrier();
}

template<class Container = std::vector<std::tuple<int, int, SPComplexType>>,
         class PsiT_Type = PsiT_Matrix_t<SPComplexType>>
inline void rotateHijkl_single_node(std::string& type,
                                    WALKER_TYPES walker_type,
                                    bool addCoulomb,
                                    TaskGroup_& TG,
                                    Container& Vijkl,
                                    PsiT_Type* Alpha,
                                    PsiT_Type* Beta,
                                    SpVType_shm_csr_matrix const& V2_fact,
                                    const RealType cut,
                                    int maximum_buffer_size,
                                    bool reserve_to_fit_ = true)
{
  int nnodes = TG.getTotalNodes();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if (nnodes != 1)
    APP_ABORT(" Error: rotateHijkl_single_node requres nnodes==1. \n");

  bool distribute_Ham = TG.getNumberOfTGs() > 1;
  if (distribute_Ham)
    APP_ABORT(" Distributed V2_fact not yet implemented. \n");

  int NAEA = Alpha->size(0);
  int NMO  = Alpha->size(1);
  int NAEB = NAEA;
  if (walker_type == COLLINEAR)
    NAEB = Beta->size(0);

  // <ab||kl> = sum_n Qk(a,n) * Rl(b,n) - Rl(a,n)*Qk(b,n),
  // where:
  //   Qk(a,n) = sum_i ma::conj(Amat(i,a)) * V2_fact(ik,n)
  //   Rl(a,n) = sum_i ma::conj(Amat(i,a)) * ma::conj(V2_fact(li,n))
  // For real build, Qk=Rk
  //
  // For parallelization, distribute (k,l) pairs over nodes.
  // Build ahead of time Qk/Rl matrices in shared memory to reduce memory/setup time.
  // Assemble integrals in parallel and fully distributed.
  // Collect on all nodes.
  //    - For distributed hamiltonians, you do not need to worry about keeping contiguous
  //    segments of the hamiltonian. Only that the distribution over nodes is roughly equal.
  //    Write a simple algorithm that balances the number of terms in a TG.
  //

  bool sparseQk = (type == "SD" || type == "SS");

  app_log() << " Calculating half-rotated Hamiltonian using ";
  if (type == "SD")
    app_log() << "Sparse x Dense";
  else if (type == "SS")
    app_log() << "Sparse x Sparse";
  else if (type == "DD")
    app_log() << "Dense x Dense";
  app_log() << " matrix multiplication. \n";

  int NEL = NAEA + (walker_type == COLLINEAR ? NAEB : 0);

  int nvec = V2_fact.size(1);

  app_log() << " Approximate memory usage for half-rotated Hamiltonian construction: \n"
            << "   max. number of orbital in a node: " << NMO << "\n"
            << "   Qk/Rl matrices size (assuming dense) each = maxnorb * nup * ncholvecs complex numbers = "
            << sizeof(SPComplexType) * NMO * NAEA * nvec / 1024.0 / 1024.0 << " MB \n"
            << "   Maximum size of communication buffer: " << maximum_buffer_size << " MB" << std::endl;


  using Allocator   = shared_allocator<SPComplexType>;
  using shmSpVector = SPComplexVector<shared_allocator<SPComplexType>>;
  using shmSpMatrix = SPComplexMatrix<shared_allocator<SPComplexType>>;
  Allocator alloc{TG.Node()};

  // allocate Qk or SpQk
  // sadly no std::optional in c++14!!!
  int nx = (sparseQk ? 0 : 1);
  shmSpMatrix Qk({nx * NMO * NEL, nx * nvec}, shared_allocator<SPComplexType>{TG.Node()});
  if (coreid == 0)
    std::fill_n(Qk.origin(), Qk.num_elements(), SPComplexType(0.0));
  nx = (sparseQk ? 1 : 0);
  SpCType_shm_csr_matrix SpQk(tp_ul_ul{nx * NMO * NEL, nx * nvec}, tp_ul_ul{0, 0}, 0, alloc);

  shmSpMatrix Rl({nvec, NMO * NEL}, shared_allocator<SPComplexType>{TG.Node()});
  if (coreid == 0)
    std::fill_n(Rl.origin(), Rl.num_elements(), SPComplexType(0.0));

  {
    //   Q(k,a,n) = sum_i ma::conj(Amat(i,a)) * V2_fact(ik,n)
    //   R(l,a,n) = sum_i ma::conj(Amat(i,a)) * ma::conj(V2_fact(li,n))

    // Construct SpQk[k,n,nvec]
    int NMO2 = (walker_type != CLOSED) ? 2 * NMO : NMO;
    if (sparseQk)
    {
      sparse_rotate::halfRotateCholeskyMatrix(walker_type, TG, 0, NMO2, SpQk, Alpha, Beta, V2_fact, false, false, cut,
                                              true);
      SpQk.remove_empty_spaces(); // just in case
    }
    else
      ma_rotate::halfRotateCholeskyMatrix(walker_type, TG, 0, NMO2, Qk, Alpha, Beta, V2_fact, false, false, cut);

#if defined(QMC_COMPLEX)
    ma_rotate::halfRotateCholeskyMatrix(walker_type, TG, 0, NMO2, Rl, Alpha, Beta, V2_fact, true, true, cut);
#else
    if (sparseQk)
    {
      csr::shm::transpose(SpQk, Rl);
    }
    else
    {
      // Qk[norb*NAEA,nvec]
      // Rl[nvec,norb*NAEA]
      int n0_, n1_, sz_ = std::get<0>(Qk.sizes());
      std::tie(n0_, n1_) = FairDivideBoundary(coreid, sz_, ncores);
      if (n1_ - n0_ > 0)
        ma::transpose(Qk.sliced(n0_, n1_), Rl(Rl.extension(0), {n0_, n1_}));
    }
#endif
  }

  TG.node_barrier();

  // set maxnk based on size of Ta
  int maxnk = std::max(1,
                       static_cast<int>(std::floor(1.0 * maximum_buffer_size /
                                                   (NMO * NAEA * NAEA * sizeof(SPComplexType) / 1024.0 / 1024.0))));
  maxnk     = std::min(maxnk, NMO);

  app_log() << "   Temporary integral matrix Ta: "
            << NMO * NAEA * maxnk * NAEA * sizeof(SPComplexType) / 1024.0 / 1024.0 << " MB " << std::endl;

  shmSpVector Ta_shmbuff(iextensions<1u>{NMO * NAEA * maxnk * NAEA}, shared_allocator<SPComplexType>{TG.Node()});
  myTimer Timer_;

  SPComplexType EJX(0.0);
  if (addCoulomb)
    EJX = SPComplexType(1.0);
  int nblk = (NMO + maxnk - 1) / maxnk;

  TG.Node().barrier();
  if (reserve_to_fit_)
  {
    // count and resize container
    std::vector<std::size_t> sz_local;
    if (walker_type == CLOSED)
      sz_local.resize(NMO * NAEA);
    else if (walker_type == COLLINEAR)
      sz_local.resize(NMO * (NAEA + NAEB));
    else if (walker_type == NONCOLLINEAR)
      sz_local.resize(2 * NMO * (NAEA + NAEB));

    for (int bi = 0; bi < nblk; bi++)
    {
      int k0 = bi * maxnk;
      int kN = std::min(k0 + maxnk, NMO);
      int nk = kN - k0;
      { // alpha-alpha
        boost::multi::array_ref<SPComplexType, 2> Ta(to_address(Ta_shmbuff.origin()), {nk * NAEA, NAEA * NMO});
        if (type == "SD")
          count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0, kN, 0, NMO, NMO, NAEA, NAEB,
                        SpQk[{size_t(k0 * NAEA), std::size_t(kN * NAEA)}], Rl(Rl.extension(0), {0, NAEA * NMO}), Ta,
                        cut);
        else if (type == "DD")
        {
          count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0, kN, 0, NMO, NMO, NAEA, NAEB,
                        Qk.sliced(size_t(k0 * NAEA), size_t(kN * NAEA)), Rl(Rl.extension(0), {0, NAEA * NMO}), Ta, cut);
        }
      }
      TG.Node().barrier();
      if (walker_type == COLLINEAR)
      { // beta-beta
        boost::multi::array_ref<SPComplexType, 2> Ta(to_address(Ta_shmbuff.origin()), {nk * NAEB, NAEB * NMO});
        if (type == "SD")
          count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0 + NMO, kN + NMO, NMO, 2 * NMO, NMO, NAEA, NAEB,
                        SpQk[{size_t(NAEA * NMO + k0 * NAEB), std::size_t(NAEA * NMO + kN * NAEB)}],
                        Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta, cut);
        else if (type == "DD")
        {
          count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0 + NMO, kN + NMO, NMO, 2 * NMO, NMO, NAEA, NAEB,
                        Qk.sliced(NAEA * NMO + k0 * NAEB, NAEA * NMO + kN * NAEB),
                        Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta, cut);
        }
        TG.Node().barrier();
        if (addCoulomb)
        { // alpha-beta
          if (type == "SD")
            count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0, kN, NMO, 2 * NMO, NMO, NAEA, NAEB,
                          SpQk[{size_t(k0 * NAEA), std::size_t(kN * NAEA)}],
                          Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta, cut);
          else if (type == "DD")
          {
            boost::multi::array_ref<SPComplexType, 2> Ta_(to_address(Ta_shmbuff.origin()), {nk * NAEA, NAEB * NMO});
            count_Qk_x_Rl(walker_type, EJX, TG, sz_local, k0, kN, NMO, 2 * NMO, NMO, NAEA, NAEB,
                          Qk.sliced(k0 * NAEA, kN * NAEA), Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta_,
                          cut);
          }
        }
        TG.Node().barrier();
      }
    }

    std::size_t tot_sz_local = std::accumulate(sz_local.begin(), sz_local.end(), std::size_t(0));

    std::vector<std::size_t> sz_global(sz_local.size());
    TG.Global().all_reduce_n(sz_local.begin(), sz_local.size(), sz_global.begin(), std::plus<>());
    std::size_t tot_sz_global = std::accumulate(sz_global.begin(), sz_global.end(), std::size_t(0));

    std::vector<std::size_t> sz_node(sz_local.size());
    TG.Node().all_reduce_n(sz_local.begin(), sz_local.size(), sz_node.begin(), std::plus<>());
    std::size_t tot_sz_node = std::accumulate(sz_node.begin(), sz_node.end(), std::size_t(0));

    std::size_t sz_node_min, sz_node_max, sz_local_min, sz_local_max;
    TG.Global().all_reduce_n(&tot_sz_node, 1, &sz_node_min, boost::mpi3::min<>());
    TG.Global().all_reduce_n(&tot_sz_node, 1, &sz_node_max, boost::mpi3::max<>());
    TG.Global().all_reduce_n(&tot_sz_local, 1, &sz_local_min, boost::mpi3::min<>());
    TG.Global().all_reduce_n(&tot_sz_local, 1, &sz_local_max, boost::mpi3::max<>());

    app_log() << "  Number of terms in Vijkl: \n"
              << "    Local: (min/max)" << sz_local_min << " "
              << sz_local_min * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB  -  " << sz_local_max
              << " " << sz_local_max * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB \n"
              << "    Node (min/max): " << sz_node_min << " "
              << sz_node_min * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB   -   " << sz_node_max
              << " " << sz_node_max * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB \n"
              << "    Global: " << tot_sz_global << " "
              << tot_sz_global * (sizeof(SPComplexType) + sizeof(int)) / 1024.0 / 1024.0 << " MB" << std::endl
              << std::endl;

    reserve_to_fit(Vijkl, sz_global);
  }
  TG.Node().barrier();

  // now calculate fully distributed matrix elements
  for (int bi = 0; bi < nblk; bi++)
  {
    int k0 = bi * maxnk;
    int kN = std::min(k0 + maxnk, NMO);
    int nk = kN - k0;
    { // alpha-alpha
      boost::multi::array_ref<SPComplexType, 2> Ta(to_address(Ta_shmbuff.origin()), {nk * NAEA, NAEA * NMO});
      if (type == "SD")
        Qk_x_Rl(walker_type, EJX, TG, k0, kN, 0, NMO, NMO, NAEA, NAEB,
                SpQk[{size_t(k0 * NAEA), std::size_t(kN * NAEA)}], Rl(Rl.extension(0), {0, NAEA * NMO}), Ta, Vijkl,
                cut);
      else if (type == "DD")
        Qk_x_Rl(walker_type, EJX, TG, k0, kN, 0, NMO, NMO, NAEA, NAEB, Qk.sliced(size_t(k0 * NAEA), size_t(kN * NAEA)),
                Rl(Rl.extension(0), {0, NAEA * NMO}), Ta, Vijkl, cut);
    }
    TG.Node().barrier();
    if (walker_type == COLLINEAR)
    { // beta-beta
      boost::multi::array_ref<SPComplexType, 2> Ta(to_address(Ta_shmbuff.origin()), {nk * NAEB, NAEB * NMO});
      if (type == "SD")
        Qk_x_Rl(walker_type, EJX, TG, k0 + NMO, kN + NMO, NMO, 2 * NMO, NMO, NAEA, NAEB,
                SpQk[{size_t(NAEA * NMO + k0 * NAEB), std::size_t(NAEA * NMO + kN * NAEB)}],
                Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta, Vijkl, cut);
      else if (type == "DD")
        Qk_x_Rl(walker_type, EJX, TG, k0 + NMO, kN + NMO, NMO, 2 * NMO, NMO, NAEA, NAEB,
                Qk.sliced(NAEA * NMO + k0 * NAEB, NAEA * NMO + kN * NAEB),
                Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta, Vijkl, cut);
      TG.Node().barrier();
      if (addCoulomb)
      { // alpha-beta
        boost::multi::array_ref<SPComplexType, 2> Ta_(to_address(Ta_shmbuff.origin()), {nk * NAEA, NAEB * NMO});
        if (type == "SD")
          Qk_x_Rl(walker_type, EJX, TG, k0, kN, NMO, 2 * NMO, NMO, NAEA, NAEB,
                  SpQk[{size_t(k0 * NAEA), std::size_t(kN * NAEA)}],
                  Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta, Vijkl, cut);
        else if (type == "DD")
          Qk_x_Rl(walker_type, EJX, TG, k0, kN, NMO, 2 * NMO, NMO, NAEA, NAEB, Qk.sliced(k0 * NAEA, kN * NAEA),
                  Rl(Rl.extension(0), {NAEA * NMO, (NAEA + NAEB) * NMO}), Ta_, Vijkl, cut);
      }
      TG.Node().barrier();
    }
  }

  TG.node_barrier();
}

} // namespace afqmc

} // namespace qmcplusplus

#endif
