#include <cstdlib>
#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <utility>
#include <vector>
#include <numeric>
#include <functional>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include "Configuration.h"
#include "hdf/hdf_multi.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Utilities/kp_utilities.hpp"
#include "KPFactorizedHamiltonian.h"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"
//#include "AFQMC/HamiltonianOperations/KP3IndexFactorizationIO.hpp"

namespace qmcplusplus
{
namespace afqmc
{
#if defined(QMC_COMPLEX)
HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations(bool pureSD,
                                                                        bool addCoulomb,
                                                                        WALKER_TYPES type,
                                                                        std::vector<PsiT_Matrix>& PsiT,
                                                                        double cutvn,
                                                                        double cutv2,
                                                                        TaskGroup_& TGprop,
                                                                        TaskGroup_& TGwfn,
                                                                        hdf_archive& hdf_restart)
{
  if (TG.TG_local().size() > 1 || not(batched == "yes" || batched == "true"))
    return getHamiltonianOperations_shared(pureSD, addCoulomb, type, PsiT, cutvn, cutv2, TGprop, TGwfn, hdf_restart);
  else
    return getHamiltonianOperations_batched(pureSD, addCoulomb, type, PsiT, cutvn, cutv2, TGprop, TGwfn, hdf_restart);
}

HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations_shared(bool pureSD,
                                                                               bool addCoulomb,
                                                                               WALKER_TYPES type,
                                                                               std::vector<PsiT_Matrix>& PsiT,
                                                                               double cutvn,
                                                                               double cutv2,
                                                                               TaskGroup_& TGprop,
                                                                               TaskGroup_& TGwfn,
                                                                               hdf_archive& hdf_restart)
{
  using shmIMatrix    = boost::multi::array<int, 2, shared_allocator<int>>;
  using shmCVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using shmCMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using shmCTensor    = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using shmSpMatrix   = boost::multi::array<SPComplexType, 2, shared_allocator<SPComplexType>>;
  using shmSpTensor   = boost::multi::array<SPComplexType, 3, shared_allocator<SPComplexType>>;
  using IVector       = boost::multi::array<int, 1>;
  using CMatrix       = boost::multi::array<ComplexType, 2>;
  using SpMatrix      = boost::multi::array<SPComplexType, 2>;
  using SpMatrix_ref  = boost::multi::array_ref<SPComplexType, 2>;
  using Sp3Tensor_ref = boost::multi::array_ref<SPComplexType, 3>;

  if (TGprop.TG() != TGwfn.TG())
  {
    app_error() << " Error: KPFactorizedHamiltonian requires nnodes to be the same in Wavefunction \n"
                << "        and Propagator xml blocks." << std::endl;
    APP_ABORT("Error: Inconsistent nnodes in KPFactorizedHamiltonian \n");
  }

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if (TGwfn.Global().root())
    write_hdf = !hdf_restart.closed();
  //  if(TGwfn.Global().root()) write_hdf = (hdf_restart.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  if (type == COLLINEAR)
    assert(PsiT.size() % 2 == 0);
  int nspins = ((type != COLLINEAR) ? 1 : 2);
  int ndet   = PsiT.size() / nspins;
  int npol   = ((type == NONCOLLINEAR) ? 2 : 1);

  if (ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  long nkpts;
  hdf_archive dump(TGwfn.Global());
  // right now only Node.root() reads
  if (TG.Node().root())
  {
    if (!dump.open(fileName, H5F_ACC_RDONLY))
    {
      app_error() << " Error opening integral file in THCHamiltonian. \n";
      APP_ABORT("");
    }
    dump.push("Hamiltonian", false);
  }

  std::vector<int> Idata(8);
  if (TG.Global().root())
  {
    if (!dump.readEntry(Idata, "dims"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading dims. \n";
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(Idata.begin(), 8, 0);
  nkpts = Idata[2];
  app_log() << " nkpts: " << nkpts << std::endl;

  // partition Q over nodes if distributed Q

  IVector nmo_per_kp(iextensions<1u>{nkpts});
  IVector nchol_per_kp(iextensions<1u>{nkpts});
  IVector kminus(iextensions<1u>{nkpts});
  IVector Qmap(iextensions<1u>{nkpts});
  shmIMatrix QKtok2({nkpts, nkpts}, shared_allocator<int>{TG.Node()});
  ValueType E0;
  if (TG.Global().root())
  {
    if (!dump.readEntry(nmo_per_kp, "NMOPerKP"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading NMOPerKP. \n";
      APP_ABORT("");
    }
    if (!dump.readEntry(nchol_per_kp, "NCholPerKP"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading NCholPerKP. \n";
      APP_ABORT("");
    }
    if (!dump.readEntry(kminus, "MinusK"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading MinusK. \n";
      APP_ABORT("");
    }
    for (int k = 0; k < nkpts; k++)
    {
      if (kminus[k] < k)
        nchol_per_kp[k] = nchol_per_kp[kminus[k]];
    }
    if (!dump.readEntry(QKtok2, "QKTok2"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading QKTok2. \n";
      APP_ABORT("");
    }
    std::vector<RealType> E_(2);
    if (!dump.readEntry(E_, "Energies"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading Energies. \n";
      APP_ABORT("");
    }
    E0 = E_[0] + E_[1];
    if (nmo_per_kp.size() != nkpts || nchol_per_kp.size() != nkpts || kminus.size() != nkpts ||
        std::get<0>(QKtok2.sizes()) != nkpts || std::get<1>(QKtok2.sizes()) != nkpts)
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Inconsistent dimension (NMOPerKP,NCholPerKP,QKtTok2): " << nkpts << " " << nmo_per_kp.size()
                  << " " << nchol_per_kp.size() << " " << kminus.size() << " " << std::get<0>(QKtok2.sizes()) << " "
                  << std::get<1>(QKtok2.sizes()) << std::endl;
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(&E0, 1, 0);
  TG.Global().broadcast_n(nmo_per_kp.origin(), nmo_per_kp.size(), 0);
  TG.Global().broadcast_n(nchol_per_kp.origin(), nchol_per_kp.size(), 0);
  TG.Global().broadcast_n(kminus.origin(), kminus.size(), 0);
  if (TG.Node().root())
    TG.Cores().broadcast_n(to_address(QKtok2.origin()), QKtok2.num_elements(), 0);
  TG.Node().barrier();

  int number_of_symmetric_Q = 0;
  int global_origin(0);
  // Defines behavior over Q vector:
  //   <0: Ignore (handled by another TG)
  //    0: Calculate, without rho^+ contribution
  //   >0: Calculate, with rho^+ contribution. LQKbln data located at Qmap[Q]-1
  std::fill_n(Qmap.origin(), Qmap.num_elements(), -1);
  {
    int ngrp(TGwfn.getNGroupsPerTG());
    int ig(TGwfn.getLocalGroupNumber());
    int work(0);
    // assign Q/Qm pairs of vectors to groups round-robin
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (kminus[Q] == Q)
      {
        if (work % ngrp == ig)
          Qmap[Q] = 1 + (number_of_symmetric_Q++);
        if (work % ngrp < ig)
          global_origin += 2 * nchol_per_kp[Q];
        work++;
      }
      else if (Q < kminus[Q])
      {
        if (work % ngrp == ig)
        {
          Qmap[Q]         = 0;
          Qmap[kminus[Q]] = 0;
        }
        if (work % ngrp < ig)
          global_origin += 4 * nchol_per_kp[Q];
        work++;
      }
    }
    if (work < ngrp)
      APP_ABORT(" Error: Too many nodes in group (nnodes) for given number of kpoints. \n");
  }
  // new communicator over nodes that share the same set of Q
  auto Qcomm(TG.Global().split(TGwfn.getLocalGroupNumber(), TG.Global().rank()));
  auto Qcomm_roots(Qcomm.split(TGwfn.Node().rank(), Qcomm.rank()));

  int nmo_max   = *std::max_element(nmo_per_kp.begin(), nmo_per_kp.end());
  int nchol_max = *std::max_element(nchol_per_kp.begin(), nchol_per_kp.end());
  shmCTensor H1({nkpts, npol * nmo_max, npol * nmo_max}, shared_allocator<ComplexType>{TG.Node()});
  std::vector<shmSpMatrix> LQKikn;
  LQKikn.reserve(nkpts);
  for (int Q = 0; Q < nkpts; Q++)
    if (Qmap[Q] >= 0 && Q <= kminus[Q])
      LQKikn.emplace_back(
          shmSpMatrix({nkpts, nmo_max * nmo_max * nchol_per_kp[Q]}, shared_allocator<SPComplexType>{TG.Node()}));
    else
      LQKikn.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{TG.Node()}));

  if (TG.Node().root())
  {
    // now read H1_kpQ
    for (int Q = 0; Q < nkpts; Q++)
    {
      // until double_hyperslabs work!
      boost::multi::array<ComplexType, 2> h1({npol * nmo_per_kp[Q], npol * nmo_per_kp[Q]});
      if (!dump.readEntry(h1, std::string("H1_kp") + std::to_string(Q)))
      {
        app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                    << " Problems reading /Hamiltonian/H1_kp" << Q << ". \n";
        APP_ABORT("");
      }
      // H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}) = h1;
      // using add to get raw pointer dispatch, otherwise matrix copy is going to sync
      ma::add(ComplexType(1.0), h1, ComplexType(0.0), h1, H1[Q]({0, npol * nmo_per_kp[Q]}, {0, npol * nmo_per_kp[Q]}));
    }
    // read LQ
    dump.push("KPFactorized", false);

    for (int Q = 0; Q < nkpts; Q++)
    {
      using ma::conj;
      if (Qmap[Q] >= 0 && Q <= kminus[Q])
      {
        if (!dump.readEntry(LQKikn[Q], std::string("L") + std::to_string(Q)))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Problems reading /Hamiltonian/KPFactorized/L" << Q << ". \n";
          APP_ABORT("");
        }
        if (std::get<0>(LQKikn[Q].sizes()) != nkpts || std::get<1>(LQKikn[Q].sizes()) != nmo_max * nmo_max * nchol_per_kp[Q])
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Problems reading /Hamiltonian/KPFactorized/L" << Q << ". \n"
                      << " Unexpected dimensins: " << std::get<0>(LQKikn[Q].sizes()) << " " << std::get<1>(LQKikn[Q].sizes()) << std::endl;
          APP_ABORT("");
        }
      }
    }
    dump.pop();
  }
  TG.Node().barrier();

  // calculate vn0
  shmCTensor vn0({nkpts, nmo_max, nmo_max}, shared_allocator<ComplexType>{TG.Node()});

  // generate nocc_per_kp using PsiT and nmo_per_kp
  shmIMatrix nocc_per_kp({ndet, nspins * nkpts}, shared_allocator<int>{TG.Node()});
  if (TG.Node().root())
  {
    if (type == COLLINEAR)
    {
      for (int i = 0; i < ndet; i++)
      {
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[2 * i], nocc_per_kp[i]({0, nkpts})))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[2 * i + 1], nocc_per_kp[i]({nkpts, 2 * nkpts})))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
      }
    }
    else
    {
      for (int i = 0; i < ndet; i++)
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[i], nocc_per_kp[i], npol == 2))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
    }
  }
  TG.Node().barrier();
  int nocc_max = *std::max_element(to_address(nocc_per_kp.origin()),
                                   to_address(nocc_per_kp.origin()) + nocc_per_kp.num_elements());

  /* half-rotate LQ and H1:
   * Given that PsiT = H(SM),
   * h[K][a][k] = sum_i PsiT[K][a][i] * h[K][i][k]
   * L[Q][K][a][k][n] = sum_i PsiT[K][a][i] * L[Q][K][i][k][n]
   * L[Q][K][l][b][n] = sum_i PsiT[K][b][] * L[Q][K][l][k][n]*
   * LQKak has a special transposition to facilitate computations
   * of the energy, and they are stored with padding to max linear dimension
   * LQKank[Q][K][...] = LQKank[Q][K][a][n][k] = LQKakn[Q][K][a][k][n]
   */
  std::vector<shmSpMatrix> LQKank;
  LQKank.reserve(ndet * nspins * nkpts); // storing 2 components for Q=0, since it is not assumed symmetric
  shmCMatrix haj({ndet * nkpts, (type == COLLINEAR ? 2 : 1) * nocc_max * npol * nmo_max},
                 shared_allocator<ComplexType>{TG.Node()});
  if (TG.Node().root())
    std::fill_n(haj.origin(), haj.num_elements(), ComplexType(0.0));
  int ank_max = nocc_max * nchol_max * nmo_max;
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < nkpts; Q++)
      if (Qmap[Q] >= 0)
        LQKank.emplace_back(shmSpMatrix({nkpts, npol * ank_max}, shared_allocator<SPComplexType>{TG.Node()}));
      else
        LQKank.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{TG.Node()}));
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < nkpts; Q++)
        if (Qmap[Q] >= 0)
          LQKank.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<SPComplexType>{TG.Node()}));
        else
          LQKank.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{TG.Node()}));
    }
  }
  for (int nd = 0, nt = 0, nq0 = 0; nd < ndet; nd++, nq0 += nkpts * nspins)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] < 0)
        continue;
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % TG.Node().size() == TG.Node().rank())
        {
          std::fill_n(to_address(LQKank[nq0 + Q][K].origin()), LQKank[nq0 + Q][K].num_elements(), SPComplexType(0.0));
          if (type == COLLINEAR)
          {
            std::fill_n(to_address(LQKank[nq0 + nkpts + Q][K].origin()), LQKank[nq0 + nkpts + Q][K].num_elements(),
                        SPComplexType(0.0));
          }
        }
      }
    }
  }

  // NOTE: LQKbnl is indexed by the K index of 'b', L[Q][Kb]
  std::vector<shmSpMatrix> LQKbnl;
  LQKbnl.reserve(ndet * nspins *
                 number_of_symmetric_Q); // storing 2 components for Q=0, since it is not assumed symmetric
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < number_of_symmetric_Q; Q++)
      LQKbnl.emplace_back(shmSpMatrix({nkpts, npol * ank_max}, shared_allocator<SPComplexType>{TG.Node()}));
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < number_of_symmetric_Q; Q++)
        LQKbnl.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<SPComplexType>{TG.Node()}));
    }
  }
  for (int nd = 0, nt = 0, nq0 = 0; nd < ndet; nd++, nq0 += number_of_symmetric_Q * nspins)
  {
    for (int Q = 0; Q < number_of_symmetric_Q; Q++)
    {
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % TG.Node().size() == TG.Node().rank())
        {
          std::fill_n(to_address(LQKbnl[nq0 + Q][K].origin()), LQKbnl[nq0 + Q][K].num_elements(), SPComplexType(0.0));
          if (type == COLLINEAR)
            std::fill_n(to_address(LQKbnl[nq0 + number_of_symmetric_Q + Q][K].origin()),
                        LQKbnl[nq0 + number_of_symmetric_Q + Q][K].num_elements(), SPComplexType(0.0));
        }
      }
    }
  }
  TG.Node().barrier();

  int Q0 = -1; // if K=(0,0,0) exists, store index here
  for (int Q = 0; Q < nkpts; Q++)
  {
    if (kminus[Q] == Q)
    {
      bool found = true;
      for (int KI = 0; KI < nkpts; KI++)
        if (KI != QKtok2[Q][KI])
        {
          found = false;
          break;
        }
      if (found)
      {
        Q0 = Q;
        break;
      }
    }
  }
  if (Q0 < 0)
    APP_ABORT(" Error: Could not find Q=0. \n");

  boost::multi::array<SPComplexType, 2> buff({npol * nmo_max, nchol_max});
  int nt = 0;
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int K = 0; K < nkpts; K++, nt++)
    {
      if (nt % TG.Global().size() == TG.Global().rank())
      {
        // haj and add half-transformed right-handed rotation for Q=0
        int na = nocc_per_kp[nd][K];
        int nb = (nspins == 2 ? nocc_per_kp[nd][nkpts + K] : na);
        int ni = nmo_per_kp[K];
        if (type == COLLINEAR)
        {
          { // Alpha
            auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], K);
            assert(std::get<0>(Psi.sizes()) == na);
            boost::multi::array_ref<ComplexType, 2> haj_r(to_address(haj[nd * nkpts + K].origin()), {na, ni});
            if (na > 0)
              ma::product(Psi, H1[K]({0, ni}, {0, ni}), haj_r);
          }
          { // Beta
            auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], K);
            assert(std::get<0>(Psi.sizes()) == nb);
            boost::multi::array_ref<ComplexType, 2> haj_r(to_address(haj[nd * nkpts + K].origin()) + na * ni, {nb, ni});
            if (nb > 0)
              ma::product(Psi, H1[K]({0, ni}, {0, ni}), haj_r);
          }
        }
        else
        {
          RealType scl = (type == CLOSED ? 2.0 : 1.0);
          auto Psi     = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[nd], K, npol == 2);
          assert(std::get<0>(Psi.sizes()) == na);
          boost::multi::array_ref<ComplexType, 2> haj_r(to_address(haj[nd * nkpts + K].origin()), {na, npol * ni});
          if (na > 0)
            ma::product(ComplexType(scl), Psi, H1[K]({0, npol * ni}, {0, npol * ni}), ComplexType(0.0), haj_r);
        }
      }
    }
  }
  // Generate LQKank
  for (int nd = 0, nq0 = 0; nd < ndet; nd++, nq0 += nkpts * nspins)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] < 0)
        continue;
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % Qcomm.size() == Qcomm.rank())
        {
          // add half-transformed right-handed rotation for Q=0
          int Qm    = kminus[Q];
          int QK    = QKtok2[Q][K];
          int na    = nocc_per_kp[nd][K];
          int nb    = (nspins == 2 ? nocc_per_kp[nd][nkpts + K] : na);
          int ni    = nmo_per_kp[K];
          int nk    = nmo_per_kp[QK];
          int nchol = nchol_per_kp[Q];
          if (type == COLLINEAR)
          {
            { // Alpha
              auto Psi = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], K);
              assert(std::get<0>(Psi.sizes()) == nocc_per_kp[nd][K]);
              if (Q <= Qm)
              {
                Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {ni, nk, nchol});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {na, nchol, nk});
                ma_rotate::getLank(Psi, Likn, Lank, buff);
              }
              else
              {
                Sp3Tensor_ref Lkin(to_address(LQKikn[Qm][QK].origin()), {nk, ni, nchol});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {na, nchol, nk});
                ma_rotate::getLank_from_Lkin(Psi, Lkin, Lank, buff);
              }
            }
            { // Beta
              auto Psi = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], K);
              assert(std::get<0>(Psi.sizes()) == nb);
              if (Q <= Qm)
              {
                Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {ni, nk, nchol});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + nkpts + Q][K].origin()), {nb, nchol, nk});
                ma_rotate::getLank(Psi, Likn, Lank, buff);
              }
              else
              {
                Sp3Tensor_ref Lkin(to_address(LQKikn[Qm][QK].origin()), {nk, ni, nchol});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + nkpts + Q][K].origin()), {nb, nchol, nk});
                ma_rotate::getLank_from_Lkin(Psi, Lkin, Lank, buff);
              }
            }
          }
          else
          {
            auto Psi = get_PsiK<SpMatrix>(nmo_per_kp, PsiT[nd], K, npol == 2);
            assert(std::get<0>(Psi.sizes()) == na);
            if (Q <= Qm)
            {
              Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {ni, nk, nchol});
              Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {na, nchol, npol * nk});
              ma_rotate::getLank(Psi, Likn, Lank, buff, npol == 2);
            }
            else
            {
              Sp3Tensor_ref Lkin(to_address(LQKikn[Qm][QK].origin()), {nk, ni, nchol});
              Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {na, nchol, npol * nk});
              ma_rotate::getLank_from_Lkin(Psi, Lkin, Lank, buff, npol == 2);
            }
          }
        }
      }
    }
  }

  // now generate LQKbnl if Q==(-Q)
  for (int nd = 0, nq0 = 0; nd < ndet; nd++, nq0 += number_of_symmetric_Q * nspins)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] <= 0)
        continue;
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % Qcomm.size() == Qcomm.rank())
        {
          // careful with subtle redefinition of na,nb,... here
          int QK    = QKtok2[Q][K];
          int na    = nocc_per_kp[nd][QK];
          int nb    = (nspins == 2 ? nocc_per_kp[nd][nkpts + QK] : na);
          int ni    = nmo_per_kp[K];
          int nk    = nmo_per_kp[QK];
          int nchol = nchol_per_kp[Q];
          Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {ni, nk, nchol});
          // NOTE: LQKbnl is indexed by the K index of 'b', L[Q][Kb]
          if (type == COLLINEAR)
          {
            { // Alpha
              auto PsiQK = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], QK);
              Sp3Tensor_ref Lbnl(to_address(LQKbnl[nq0 + Qmap[Q] - 1][QK].origin()), {na, nchol, ni});
              ma_rotate::getLank_from_Lkin(PsiQK, Likn, Lbnl, buff);
            }
            { // Beta
              auto PsiQK = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], QK);
              assert(std::get<0>(PsiQK.sizes()) == nb);
              Sp3Tensor_ref Lbnl(to_address(LQKbnl[nq0 + number_of_symmetric_Q + Qmap[Q] - 1][QK].origin()),
                                 {nb, nchol, ni});
              ma_rotate::getLank_from_Lkin(PsiQK, Likn, Lbnl, buff);
            }
          }
          else
          {
            auto PsiQK = get_PsiK<SpMatrix>(nmo_per_kp, PsiT[nd], QK, npol == 2);
            assert(std::get<0>(PsiQK.sizes()) == na);
            Sp3Tensor_ref Lbnl(to_address(LQKbnl[nq0 + Qmap[Q] - 1][QK].origin()), {na, nchol, npol * ni});
            ma_rotate::getLank_from_Lkin(PsiQK, Likn, Lbnl, buff, npol == 2);
          }
        }
      }
    }
  }
  Qcomm.barrier();
  if (TG.Node().root())
  {
    TG.Cores().all_reduce_in_place_n(to_address(haj.origin()), haj.num_elements(), std::plus<>());
    for (int Q = 0; Q < LQKank.size(); Q++)
      Qcomm_roots.all_reduce_in_place_n(to_address(LQKank[Q].origin()), LQKank[Q].num_elements(), std::plus<>());
    for (int Q = 0; Q < LQKbnl.size(); Q++)
      Qcomm_roots.all_reduce_in_place_n(to_address(LQKbnl[Q].origin()), LQKbnl[Q].num_elements(), std::plus<>());
    std::fill_n(to_address(vn0.origin()), vn0.num_elements(), ComplexType(0.0));
  }
  // need to broadcast haj from root of Qcomm with Qsym[0]>=0, to all other ones
  // NOTE NOTE NOTE
  TG.Node().barrier();

  // calculate vn0(I,L) = -0.5 sum_K sum_j sum_n L[0][K][i][j][n] ma::conj(L[0][K][l][j][n])
  for (int Q = 0; Q < nkpts; Q++)
  {
    if (Qmap[Q] < 0)
      continue;
    for (int K = 0; K < nkpts; K++)
    {
      if (K % TG.Node().size() == TG.Node().rank())
      {
        int QK = QKtok2[Q][K];
        int Qm = kminus[Q];
        if (Q <= Qm)
        {
          boost::multi::array_ref<SPComplexType, 2> Likn(to_address(LQKikn[Q][K].origin()),
                                                         {nmo_per_kp[K], nmo_per_kp[QK] * nchol_per_kp[Q]});
          using ma::H;
#if defined(MIXED_PRECISION)
          boost::multi::array<SPComplexType, 2> v1_({nmo_per_kp[K], nmo_per_kp[K]});
          ma::product(SPComplexType(-0.5), Likn, H(Likn), SPComplexType(0.0), v1_);
          boost::multi::array<ComplexType, 2> v2_(v1_);
          ma::add(ComplexType(1.0), v2_, ComplexType(1.0), vn0[K]({0, nmo_per_kp[K]}, {0, nmo_per_kp[K]}),
                  vn0[K]({0, nmo_per_kp[K]}, {0, nmo_per_kp[K]}));
#else
          ma::product(-0.5, Likn, H(Likn), 1.0, vn0[K]({0, nmo_per_kp[K]}, {0, nmo_per_kp[K]}));
#endif
        }
        else
        {
          int QmK = QKtok2[Qm][K];
          boost::multi::array_ref<SPComplexType, 3> Lkin(to_address(LQKikn[Qm][QK].origin()),
                                                         {nmo_per_kp[QK], nmo_per_kp[K], nchol_per_kp[Qm]});
          boost::multi::array<SPComplexType, 3> buff3D({nmo_per_kp[K], nmo_per_kp[QK], nchol_per_kp[Qm]});
          using ma::conj;
          for (int i = 0; i < nmo_per_kp[K]; i++)
            for (int k = 0; k < nmo_per_kp[QK]; k++)
              for (int n = 0; n < nchol_per_kp[Qm]; n++)
                buff3D[i][k][n] = ma::conj(Lkin[k][i][n]);
          boost::multi::array_ref<SPComplexType, 2> L_(to_address(buff3D.origin()),
                                                       {nmo_per_kp[K], nmo_per_kp[QK] * nchol_per_kp[Qm]});
          using ma::H;
#if defined(MIXED_PRECISION)
          boost::multi::array<SPComplexType, 2> v1_({nmo_per_kp[K], nmo_per_kp[K]});
          ma::product(SPComplexType(-0.5), L_, H(L_), SPComplexType(0.0), v1_);
          boost::multi::array<ComplexType, 2> v2_(v1_);
          ma::add(ComplexType(1.0), v2_, ComplexType(1.0), vn0[K]({0, nmo_per_kp[K]}, {0, nmo_per_kp[K]}),
                  vn0[K]({0, nmo_per_kp[K]}, {0, nmo_per_kp[K]}));
#else
          ma::product(-0.5, L_, H(L_), 1.0, vn0[K]({0, nmo_per_kp[K]}, {0, nmo_per_kp[K]}));
#endif
        }
      }
    }
    // need sync here to avoid having multiple Q's overwritting each other
    // either this or you need local storage
    TG.Node().barrier();
  }
  TG.Node().barrier();

  if (TG.Node().root())
  {
    dump.pop();
    dump.close();
  }

  int global_ncvecs = 2 * std::accumulate(nchol_per_kp.begin(), nchol_per_kp.end(), 0);

  std::vector<RealType> gQ(nkpts);
  if (nsampleQ > 0)
  {
    app_log() << " Sampling EXX energy using distribution over Q vector obtained from "
              << " trial energy. \n";

    if (npol == 2)
      APP_ABORT("Error: nsampleQ>0 not yet implemented for noncollinear.\n\n\n");

    RealType scl = (type == CLOSED ? 2.0 : 1.0);
    size_t nqk   = 0;
    for (int Q = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Qmap[Q] < 0)
        continue;
      int Qm = kminus[Q];
      for (int Ka = 0; Ka < nkpts; ++Ka)
      {
        int Kk = QKtok2[Q][Ka];
        int Kb = Kk;
        int Kl = QKtok2[Qm][Kb];
        if ((Ka != Kl) || (Kb != Kk))
          APP_ABORT(" Error: Problems with EXX.\n");
        if ((nqk++) % Qcomm.size() == Qcomm.rank())
        {
          int nchol = nchol_per_kp[Q];
          int nl    = nmo_per_kp[Kl];
          int nb    = nocc_per_kp[0][Kb];
          int nk    = nmo_per_kp[Kk];
          int na    = nocc_per_kp[0][Ka];

          if (na == 0 || nb == 0)
            continue;

          SpMatrix_ref Lank(to_address(LQKank[Q][Ka].origin()), {na * nchol, nk});
          auto bnl_ptr(to_address(LQKank[Qm][Kb].origin()));
          if (Qmap[Q] > 0)
            bnl_ptr = to_address(LQKbnl[Qmap[Q] - 1][Kb].origin());
          SpMatrix_ref Lbnl(bnl_ptr, {nb * nchol, nl});

          SpMatrix Tban({nb, na * nchol});
          Sp3Tensor_ref T3ban(Tban.origin(), {nb, na, nchol});
          SpMatrix Tabn({na, nb * nchol});
          Sp3Tensor_ref T3abn(Tabn.origin(), {na, nb, nchol});

          auto Gal = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[0], Ka, npol == 2);
          auto Gbk = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[0], Kb, npol == 2);
          for (int a = 0; a < na; ++a)
            for (int l = 0; l < nl; ++l)
              Gal[a][l] = ma::conj(Gal[a][l]);
          for (int b = 0; b < nb; ++b)
            for (int k = 0; k < nk; ++k)
              Gbk[b][k] = ma::conj(Gbk[b][k]);

          ma::product(Gal, ma::T(Lbnl), Tabn);
          ma::product(Gbk, ma::T(Lank), Tban);

          ComplexType E_(0.0);
          for (int a = 0; a < na; ++a)
            for (int b = 0; b < nb; ++b)
              E_ += static_cast<SPComplexType>(ma::dot(T3abn[a][b], T3ban[b][a]));
          gQ[Q] -= scl * 0.5 * real(E_);
        }
        if (type == COLLINEAR)
        {
          APP_ABORT(" Finish UHF.\n ");
        }
      }
    }
    TG.Global().all_reduce_in_place_n(gQ.begin(), nkpts, std::plus<>());
    RealType E_ = std::accumulate(gQ.begin(), gQ.end(), RealType(0.0));
    for (auto& v : gQ)
      v /= E_;
    app_log() << " EXX: " << E_ << std::endl;
    for (auto v : gQ)
    {
      if (v < 0.0)
        APP_ABORT(" Error: g(Q) < 0.0, implement shift to g(Q). \n")
    }
  }

  //  if(write_hdf)
  //    writeKP3IndexFactorization(hdf_restart,type,NMO,NAEA,NAEB,TGprop,TGwfn,
  //                            nmo_per_kp,nchol_per_kp,kminus,QKtok2,H1,LQKikn,
  //                            vn0,nsampleQ,gQ,E0,global_ncvecs);

  return HamiltonianOperations(
      KP3IndexFactorization(TGwfn.TG_local(), type, std::move(nmo_per_kp), std::move(nchol_per_kp), std::move(kminus),
                            std::move(nocc_per_kp), std::move(QKtok2), std::move(H1), std::move(haj), std::move(LQKikn),
                            std::move(LQKank), std::move(LQKbnl), std::move(Qmap), std::move(vn0), std::move(gQ),
                            nsampleQ, E0, global_origin, global_ncvecs));
}

HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations_batched(bool pureSD,
                                                                                bool addCoulomb,
                                                                                WALKER_TYPES type,
                                                                                std::vector<PsiT_Matrix>& PsiT,
                                                                                double cutvn,
                                                                                double cutv2,
                                                                                TaskGroup_& TGprop,
                                                                                TaskGroup_& TGwfn,
                                                                                hdf_archive& hdf_restart)
{
  // For now doing setup in CPU and moving structures to GPU in HamOps constructor
  using shmIMatrix    = boost::multi::array<int, 2, shared_allocator<int>>;
  using shmCVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using shmCMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using shmCTensor    = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using stdCTensor    = boost::multi::array<ComplexType, 3>;
  using shmSpMatrix   = boost::multi::array<SPComplexType, 2, shared_allocator<SPComplexType>>;
  using shmSpTensor   = boost::multi::array<SPComplexType, 3, shared_allocator<SPComplexType>>;
  using IVector       = boost::multi::array<int, 1>;
  using CMatrix       = boost::multi::array<ComplexType, 2>;
  using SpMatrix      = boost::multi::array<SPComplexType, 2>;
  using SpMatrix_ref  = boost::multi::array_ref<SPComplexType, 2>;
  using Sp3Tensor_ref = boost::multi::array_ref<SPComplexType, 3>;
  using Sp4Tensor_ref = boost::multi::array_ref<SPComplexType, 4>;

  if (TGprop.TG() != TGwfn.TG())
  {
    app_error() << " Error: KPFactorizedHamiltonian requires nnodes to be the same in Wavefunction \n"
                << "        and Propagator xml blocks." << std::endl;
    APP_ABORT("Error: Inconsistent nnodes in KPFactorizedHamiltonian \n");
  }

  if (TG.TG_local().size() > 1)
    APP_ABORT(" Error: KPFactorizedHamiltonian::getHamiltonianOperations_batched expects ncores=1. \n");

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if (TGwfn.Global().root())
    write_hdf = (not hdf_restart.closed());
  TGwfn.Global().broadcast_value(write_hdf);

  if (type == COLLINEAR)
    assert(PsiT.size() % 2 == 0);
  int nspins = ((type != COLLINEAR) ? 1 : 2);
  int ndet   = PsiT.size() / nspins;
  int npol   = ((type == NONCOLLINEAR) ? 2 : 1);

  if (ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  auto Qcomm(TG.Global().split(TGwfn.getLocalGroupNumber(), TG.Global().rank()));
  auto distNode(TG.Node().split(TGwfn.getLocalGroupNumber(), TG.Node().rank()));
  auto Qcomm_roots(Qcomm.split(distNode.rank(), Qcomm.rank()));

  long nkpts;
  hdf_archive dump(TGwfn.Global());
  // right now only Node.root() reads
  if (distNode.root())
  {
    if (!dump.open(fileName, H5F_ACC_RDONLY))
    {
      app_error() << " Error opening integral file in THCHamiltonian. \n";
      APP_ABORT("");
    }
    dump.push("Hamiltonian", false);
  }

  std::vector<int> Idata(8);
  if (TG.Global().root())
  {
    if (!dump.readEntry(Idata, "dims"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading dims. \n";
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(Idata.begin(), 8, 0);
  nkpts = Idata[2];
  app_log() << " nkpts: " << nkpts << std::endl;

  // partition Q over nodes if distributed Q

  IVector nmo_per_kp(iextensions<1u>{nkpts});
  IVector nchol_per_kp(iextensions<1u>{nkpts});
  IVector kminus(iextensions<1u>{nkpts});
  IVector Qmap(iextensions<1u>{nkpts});
  shmIMatrix QKtok2({nkpts, nkpts}, shared_allocator<int>{TG.Node()});
  ValueType E0;
  if (TG.Global().root())
  {
    if (!dump.readEntry(nmo_per_kp, "NMOPerKP"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading NMOPerKP. \n";
      APP_ABORT("");
    }
    if (!dump.readEntry(nchol_per_kp, "NCholPerKP"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading NCholPerKP. \n";
      APP_ABORT("");
    }
    if (!dump.readEntry(kminus, "MinusK"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading MinusK. \n";
      APP_ABORT("");
    }
    for (int k = 0; k < nkpts; k++)
    {
      if (kminus[k] < k)
        nchol_per_kp[k] = nchol_per_kp[kminus[k]];
    }
    if (!dump.readEntry(QKtok2, "QKTok2"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading QKTok2. \n";
      APP_ABORT("");
    }
    std::vector<RealType> E_(2);
    if (!dump.readEntry(E_, "Energies"))
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Problems reading Energies. \n";
      APP_ABORT("");
    }
    E0 = E_[0] + E_[1];
    if (nmo_per_kp.size() != nkpts || nchol_per_kp.size() != nkpts || kminus.size() != nkpts ||
        std::get<0>(QKtok2.sizes()) != nkpts || std::get<1>(QKtok2.sizes()) != nkpts)
    {
      app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                  << " Inconsistent dimension (NMOPerKP,NCholPerKP,QKtTok2): " << nkpts << " " << nmo_per_kp.size()
                  << " " << nchol_per_kp.size() << " " << kminus.size() << " " << std::get<0>(QKtok2.sizes()) << " "
                  << std::get<1>(QKtok2.sizes()) << std::endl;
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(&E0, 1, 0);
  TG.Global().broadcast_n(nmo_per_kp.origin(), nmo_per_kp.size(), 0);
  TG.Global().broadcast_n(nchol_per_kp.origin(), nchol_per_kp.size(), 0);
  TG.Global().broadcast_n(kminus.origin(), kminus.size(), 0);
  if (TG.Node().root())
    TG.Cores().broadcast_n(to_address(QKtok2.origin()), QKtok2.num_elements(), 0);
  TG.Node().barrier();

  // Defines behavior over Q vector:
  //   <0: Ignore (handled by another TG)
  //    0: Calculate, without rho^+ contribution
  //   >0: Calculate, with rho^+ contribution. LQKbln data located at Qmap[Q]-1
  int number_of_symmetric_Q = 0;
  int global_origin(0);
  std::fill_n(Qmap.origin(), Qmap.num_elements(), -1);
  {
    int ngrp(TGwfn.getNGroupsPerTG());
    int ig(TGwfn.getLocalGroupNumber());
    int work(0);
    // assign Q/Qm pairs of vectors to groups round-robin
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (kminus[Q] == Q)
      {
        if (work % ngrp == ig)
          Qmap[Q] = 1 + (number_of_symmetric_Q++);
        if (work % ngrp < ig)
          global_origin += 2 * nchol_per_kp[Q];
        work++;
      }
      else if (Q < kminus[Q])
      {
        if (work % ngrp == ig)
        {
          Qmap[Q]         = 0;
          Qmap[kminus[Q]] = 0;
        }
        if (work % ngrp < ig)
          global_origin += 4 * nchol_per_kp[Q];
        work++;
      }
    }
    if (work < ngrp)
      APP_ABORT(" Error: Too many nodes in group (nnodes) for given number of kpoints. \n");
  }

  int nmo_max   = *std::max_element(nmo_per_kp.begin(), nmo_per_kp.end());
  int nchol_max = *std::max_element(nchol_per_kp.begin(), nchol_per_kp.end());
  shmCTensor H1({nkpts, npol * nmo_max, npol * nmo_max}, shared_allocator<ComplexType>{TG.Node()});
  std::vector<shmSpMatrix> LQKikn;
  LQKikn.reserve(nkpts);
  for (int Q = 0; Q < nkpts; Q++)
    if (Qmap[Q] >= 0 && Q <= kminus[Q])
      LQKikn.emplace_back(
          shmSpMatrix({nkpts, nmo_max * nmo_max * nchol_max}, shared_allocator<SPComplexType>{distNode}));
    else // Q > kminus[Q]
      LQKikn.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{distNode}));

  if (TG.Node().root())
  {
    // now read H1_kpQ
    for (int Q = 0; Q < nkpts; Q++)
    {
      // until double_hyperslabs work!
      boost::multi::array<ComplexType, 2> h1({npol * nmo_per_kp[Q], npol * nmo_per_kp[Q]});
      if (!dump.readEntry(h1, std::string("H1_kp") + std::to_string(Q)))
      {
        app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                    << " Problems reading /Hamiltonian/H1_kp" << Q << ". \n";
        APP_ABORT("");
      }
      //H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}) = h1;
      ma::add(ComplexType(1.0), h1, ComplexType(0.0), h1, H1[Q]({0, npol * nmo_per_kp[Q]}, {0, npol * nmo_per_kp[Q]}));
    }
  }
  if (distNode.root())
  {
    for (auto& v : LQKikn)
      std::fill_n(to_address(v.origin()), v.num_elements(), SPComplexType(0.0));
    // read LQ
    dump.push("KPFactorized", false);
    // read in compact form and transform to padded
    SpMatrix L_({1, 1});
    for (int Q = 0; Q < nkpts; Q++)
    {
      using ma::conj;
      int nchol = nchol_per_kp[Q];
      if (Qmap[Q] >= 0 && Q <= kminus[Q])
      {
        if (!dump.readEntry(L_, std::string("L") + std::to_string(Q)))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Problems reading /Hamiltonian/KPFactorized/L" << Q << ". \n";
          APP_ABORT("");
        }
        assert(L_.size() == nkpts);
        Sp4Tensor_ref L2(to_address(LQKikn[Q].origin()), {nkpts, nmo_max, nmo_max, nchol_max});
        for (int K = 0; K < nkpts; ++K)
        {
          int QK = QKtok2[Q][K];
          int ni = nmo_per_kp[K];
          int nk = nmo_per_kp[QK];
          Sp3Tensor_ref L1(to_address(L_[K].origin()), {ni, nk, nchol});
          for (int i = 0; i < ni; i++)
            for (int k = 0; k < nk; k++)
              copy_n(L1[i][k].origin(), nchol, L2[K][i][k].origin());
        }
      }
    }
    dump.pop();
  }
  TG.Node().barrier();

  // calculate vn0
  shmCTensor vn0({nkpts, nmo_max, nmo_max}, shared_allocator<ComplexType>{TG.Node()});

  // generate nocc_per_kp using PsiT and nmo_per_kp
  shmIMatrix nocc_per_kp({ndet, nspins * nkpts}, shared_allocator<int>{TG.Node()});
  TG.Node().barrier();
  if (TG.Node().root())
  {
    if (type == COLLINEAR)
    {
      for (int i = 0; i < ndet; i++)
      {
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[2 * i], nocc_per_kp[i]({0, nkpts})))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[2 * i + 1], nocc_per_kp[i]({nkpts, 2 * nkpts})))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
      }
    }
    else
    {
      for (int i = 0; i < ndet; i++)
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[i], nocc_per_kp[i], npol == 2))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
    }
  }
  TG.Node().barrier();
  int nocc_max = *std::max_element(to_address(nocc_per_kp.origin()),
                                   to_address(nocc_per_kp.origin()) + nocc_per_kp.num_elements());

  int nocc_tot = std::accumulate(to_address(nocc_per_kp.origin()),
                                 to_address(nocc_per_kp.origin()) + nocc_per_kp.num_elements(), 0);
  app_log() << " Total number of electrons: " << nocc_tot << std::endl;

  /* half-rotate LQ and H1:
   * Given that PsiT = H(SM),
   * h[K][a][k] = sum_i PsiT[K][a][i] * h[K][i][k]
   * L[Q][K][a][k][n] = sum_i PsiT[K][a][i] * L[Q][K][i][k][n]
   * Both permutations are stores, akn and ank, for performance reasons.
   */
  std::vector<shmSpMatrix> LQKank;
  LQKank.reserve(ndet * nspins * nkpts); // storing 2 components for Q=0, since it is not assumed symmetric
  shmCMatrix haj({ndet * nkpts, (type == COLLINEAR ? 2 : 1) * nocc_max * npol * nmo_max},
                 shared_allocator<ComplexType>{TG.Node()});
  if (TG.Node().root())
    std::fill_n(to_address(haj.origin()), haj.num_elements(), ComplexType(0.0));
  int ank_max = nocc_max * nchol_max * nmo_max;
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < nkpts; Q++)
      if (Qmap[Q] >= 0)
        LQKank.emplace_back(shmSpMatrix({nkpts, npol * ank_max}, shared_allocator<SPComplexType>{distNode}));
      else
        LQKank.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{distNode}));
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < nkpts; Q++)
        if (Qmap[Q] >= 0)
          LQKank.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<SPComplexType>{distNode}));
        else
          LQKank.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{distNode}));
    }
  }
  if (distNode.root())
    for (auto& v : LQKank)
      std::fill_n(to_address(v.origin()), v.num_elements(), SPComplexType(0.0));

  std::vector<shmSpMatrix> LQKakn;
  LQKakn.reserve(ndet * nspins * nkpts);
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] >= 0)
        LQKakn.emplace_back(shmSpMatrix({nkpts, npol * ank_max}, shared_allocator<SPComplexType>{distNode}));
      else
        LQKakn.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{distNode}));
    }
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < nkpts; Q++)
      {
        if (Qmap[Q] >= 0)
          LQKakn.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<SPComplexType>{distNode}));
        else
          LQKakn.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{distNode}));
      }
    }
  }
  if (distNode.root())
    for (auto& v : LQKakn)
      std::fill_n(to_address(v.origin()), v.num_elements(), SPComplexType(0.0));
  // NOTE: LQKbnl and LQKbln are indexed by the K index of 'b', L[Q][Kb]
  std::vector<shmSpMatrix> LQKbnl;
  LQKbnl.reserve(ndet * nspins *
                 number_of_symmetric_Q); // storing 2 components for Q=0, since it is not assumed symmetric
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < number_of_symmetric_Q; Q++)
      LQKbnl.emplace_back(shmSpMatrix({nkpts, npol * ank_max}, shared_allocator<SPComplexType>{distNode}));
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < number_of_symmetric_Q; Q++)
        LQKbnl.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<SPComplexType>{distNode}));
    }
  }
  if (distNode.root())
    for (auto& v : LQKbnl)
      std::fill_n(to_address(v.origin()), v.num_elements(), SPComplexType(0.0));

  std::vector<shmSpMatrix> LQKbln;
  LQKbln.reserve(ndet * nspins *
                 number_of_symmetric_Q); // storing 2 components for Q=0, since it is not assumed symmetric
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < number_of_symmetric_Q; Q++)
    {
      LQKbln.emplace_back(shmSpMatrix({nkpts, npol * ank_max}, shared_allocator<SPComplexType>{distNode}));
    }
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < number_of_symmetric_Q; Q++)
      {
        LQKbln.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<SPComplexType>{distNode}));
      }
    }
  }
  if (distNode.root())
    for (auto& v : LQKbln)
      std::fill_n(to_address(v.origin()), v.num_elements(), SPComplexType(0.0));

  int Q0 = -1; // if K=(0,0,0) exists, store index here
  for (int Q = 0; Q < nkpts; Q++)
  {
    if (kminus[Q] == Q)
    {
      bool found = true;
      for (int KI = 0; KI < nkpts; KI++)
        if (KI != QKtok2[Q][KI])
        {
          found = false;
          break;
        }
      if (found)
      {
        Q0 = Q;
        break;
      }
    }
  }
  if (Q0 < 0)
    APP_ABORT(" Error: Could not find Q=0. \n");

  TG.Node().barrier();
  boost::multi::array<SPComplexType, 2> buff({npol * nmo_max, nchol_max});
  int nt = 0;
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int K = 0; K < nkpts; K++, nt++)
    {
      if (nt % TG.Global().size() == TG.Global().rank())
      {
        // haj and add half-transformed right-handed rotation for Q=0
        int na = nocc_per_kp[nd][K];
        int nb = (nspins == 2 ? nocc_per_kp[nd][nkpts + K] : na);
        int ni = nmo_per_kp[K];
        int nk = ni;
        if (type == COLLINEAR)
        {
          { // Alpha
            auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], K);
            assert(Psi.size() == na);
            boost::multi::array_ref<ComplexType, 2> haj_r(to_address(haj[nd * nkpts + K].origin()),
                                                          {nocc_max, nmo_max});
            if (na > 0)
              ma::product(Psi, H1[K]({0, ni}, {0, ni}), haj_r({0, na}, {0, ni}));
          }
          { // Beta
            auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], K);
            assert(Psi.size() == nb);
            boost::multi::array_ref<ComplexType, 2> haj_r(to_address(haj[nd * nkpts + K].origin()) + nocc_max * nmo_max,
                                                          {nocc_max, nmo_max});
            if (nb > 0)
              ma::product(Psi, H1[K]({0, ni}, {0, ni}), haj_r({0, nb}, {0, ni}));
          }
        }
        else
        {
          RealType scl = (type == CLOSED ? 2.0 : 1.0);
          auto Psi     = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[nd], K, npol == 2);
          assert(Psi.size() == na);
          boost::multi::array_ref<ComplexType, 2> haj_r(to_address(haj[nd * nkpts + K].origin()),
                                                        {nocc_max, npol * nmo_max});
          if (na > 0)
            ma::product(ComplexType(scl), Psi, H1[K]({0, npol * ni}, {0, npol * ni}), ComplexType(0.0),
                        haj_r({0, na}, {0, npol * ni}));
        }
      }
    }
  }
  for (int nd = 0, nq0 = 0; nd < ndet; nd++, nq0 += nkpts * nspins)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] < 0)
        continue;
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % Qcomm.size() == Qcomm.rank())
        {
          // add half-transformed right-handed rotation for Q=0
          int Qm = kminus[Q];
          int QK = QKtok2[Q][K];
          int na = nocc_per_kp[nd][K];
          int nb = (nspins == 2 ? nocc_per_kp[nd][nkpts + K] : na);
          int ni = nmo_per_kp[K];
          int nk = nmo_per_kp[QK];
          if (type == COLLINEAR)
          {
            { // Alpha
              auto Psi = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], K);
              assert(Psi.size() == na);
              if (Q <= Qm)
              {
                Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {nmo_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lakn(to_address(LQKakn[nq0 + Q][K].origin()), {nocc_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {nocc_max, nchol_max, nmo_max});
                ma_rotate_padded::getLakn_Lank(Psi, Likn, Lakn, Lank);
              }
              else
              {
                Sp3Tensor_ref Lkin(to_address(LQKikn[Qm][QK].origin()), {nmo_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {nocc_max, nchol_max, nmo_max});
                Sp3Tensor_ref Lakn(to_address(LQKakn[nq0 + Q][K].origin()), {nocc_max, nmo_max, nchol_max});
                ma_rotate_padded::getLakn_Lank_from_Lkin(Psi, Lkin, Lakn, Lank, buff);
              }
            }
            { // Beta
              auto Psi = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], K);
              assert(Psi.size() == nb);
              if (Q <= Qm)
              {
                Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {nmo_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lakn(to_address(LQKakn[nq0 + nkpts + Q][K].origin()), {nocc_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + nkpts + Q][K].origin()), {nocc_max, nchol_max, nmo_max});
                ma_rotate_padded::getLakn_Lank(Psi, Likn, Lakn, Lank);
              }
              else
              {
                Sp3Tensor_ref Lkin(to_address(LQKikn[Qm][QK].origin()), {nmo_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lakn(to_address(LQKakn[nq0 + nkpts + Q][K].origin()), {nocc_max, nmo_max, nchol_max});
                Sp3Tensor_ref Lank(to_address(LQKank[nq0 + nkpts + Q][K].origin()), {nocc_max, nchol_max, nmo_max});
                ma_rotate_padded::getLakn_Lank_from_Lkin(Psi, Lkin, Lakn, Lank, buff);
              }
            }
          }
          else
          {
            auto Psi = get_PsiK<SpMatrix>(nmo_per_kp, PsiT[nd], K, npol == 2);
            assert(Psi.size() == na);
            if (Q <= Qm)
            {
              Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {nmo_max, nmo_max, nchol_max});
              Sp3Tensor_ref Lakn(to_address(LQKakn[nq0 + Q][K].origin()), {nocc_max, npol * nmo_max, nchol_max});
              Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {nocc_max, nchol_max, npol * nmo_max});
              ma_rotate_padded::getLakn_Lank(Psi, Likn, Lakn, Lank, npol == 2);
            }
            else
            {
              Sp3Tensor_ref Lkin(to_address(LQKikn[Qm][QK].origin()), {nmo_max, nmo_max, nchol_max});
              Sp3Tensor_ref Lakn(to_address(LQKakn[nq0 + Q][K].origin()), {nocc_max, npol * nmo_max, nchol_max});
              Sp3Tensor_ref Lank(to_address(LQKank[nq0 + Q][K].origin()), {nocc_max, nchol_max, npol * nmo_max});
              ma_rotate_padded::getLakn_Lank_from_Lkin(Psi, Lkin, Lakn, Lank, buff, npol == 2);
            }
          }
        }
      }
    }
  }

  // now generate LQKbnl if Q==(-Q)
  for (int nd = 0, nq0 = 0; nd < ndet; nd++, nq0 += number_of_symmetric_Q * nspins)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] <= 0)
        continue;
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % Qcomm.size() == Qcomm.rank())
        {
          // careful with subtle redefinition of na,nb,... here
          int QK = QKtok2[Q][K];
          int na = nocc_per_kp[nd][QK];
          int nb = (nspins == 2 ? nocc_per_kp[nd][nkpts + QK] : na);
          Sp3Tensor_ref Likn(to_address(LQKikn[Q][K].origin()), {nmo_max, nmo_max, nchol_max});
          if (type == COLLINEAR)
          {
            { // Alpha
              auto PsiQK = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], QK);
              assert(PsiQK.size() == na);
              Sp3Tensor_ref Lbnl(to_address(LQKbnl[nq0 + Qmap[Q] - 1][QK].origin()), {nocc_max, nchol_max, nmo_max});
              Sp3Tensor_ref Lbln(to_address(LQKbln[nq0 + Qmap[Q] - 1][QK].origin()), {nocc_max, nmo_max, nchol_max});
              ma_rotate_padded::getLakn_Lank_from_Lkin(PsiQK, Likn, Lbln, Lbnl, buff);
            }
            { // Beta
              auto PsiQK = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], QK);
              assert(PsiQK.size() == nb);
              Sp3Tensor_ref Lbnl(to_address(LQKbnl[nq0 + number_of_symmetric_Q + Qmap[Q] - 1][QK].origin()),
                                 {nocc_max, nchol_max, nmo_max});
              Sp3Tensor_ref Lbln(to_address(LQKbln[nq0 + number_of_symmetric_Q + Qmap[Q] - 1][QK].origin()),
                                 {nocc_max, nmo_max, nchol_max});
              ma_rotate_padded::getLakn_Lank_from_Lkin(PsiQK, Likn, Lbln, Lbnl, buff);
            }
          }
          else
          {
            auto PsiQK = get_PsiK<SpMatrix>(nmo_per_kp, PsiT[nd], QK, npol == 2);
            assert(PsiQK.size() == na);
            Sp3Tensor_ref Lbnl(to_address(LQKbnl[nq0 + Qmap[Q] - 1][QK].origin()),
                               {nocc_max, nchol_max, npol * nmo_max});
            Sp3Tensor_ref Lbln(to_address(LQKbln[nq0 + Qmap[Q] - 1][QK].origin()),
                               {nocc_max, npol * nmo_max, nchol_max});
            ma_rotate_padded::getLakn_Lank_from_Lkin(PsiQK, Likn, Lbln, Lbnl, buff, npol == 2);
          }
        }
      }
    }
  }
  Qcomm.barrier();

  if (TG.Node().root())
  {
    TG.Cores().all_reduce_in_place_n(to_address(haj.origin()), haj.num_elements(), std::plus<>());
    std::fill_n(to_address(vn0.origin()), vn0.num_elements(), ComplexType(0.0));
  }

  if (distNode.root())
  {
    for (int Q = 0; Q < LQKank.size(); Q++)
      Qcomm_roots.all_reduce_in_place_n(to_address(LQKank[Q].origin()), LQKank[Q].num_elements(), std::plus<>());

    for (int Q = 0; Q < LQKakn.size(); Q++)
      Qcomm_roots.all_reduce_in_place_n(to_address(LQKakn[Q].origin()), LQKakn[Q].num_elements(), std::plus<>());

    for (int Q = 0; Q < LQKbnl.size(); Q++)
      Qcomm_roots.all_reduce_in_place_n(to_address(LQKbnl[Q].origin()), LQKbnl[Q].num_elements(), std::plus<>());

    for (int Q = 0; Q < LQKbln.size(); Q++)
      Qcomm_roots.all_reduce_in_place_n(to_address(LQKbln[Q].origin()), LQKbln[Q].num_elements(), std::plus<>());
  }
  TG.Node().barrier();

  // local storage seems necessary
  stdCTensor vn0_({nkpts, nmo_max, nmo_max}, ComplexType(0.0));

  // calculate vn0(I,L) = -0.5 sum_K sum_j sum_n L[0][K][i][j][n] ma::conj(L[0][K][l][j][n])
  nt = 0;
  for (int Q = 0; Q < nkpts; Q++)
  {
    if (Qmap[Q] < 0)
      continue;
    for (int K = 0; K < nkpts; K++)
    {
      if (nt % Qcomm.size() == Qcomm.rank())
      {
        int QK = QKtok2[Q][K];
        int Qm = kminus[Q];
        if (Q <= Qm)
        {
          boost::multi::array_ref<SPComplexType, 2> Likn(to_address(LQKikn[Q][K].origin()),
                                                         {nmo_max, nmo_max * nchol_max});
          using ma::H;
#if defined(MIXED_PRECISION)
          boost::multi::array<SPComplexType, 2> v1_({nmo_max, nmo_max});
          ma::product(SPComplexType(-0.5), Likn, H(Likn), SPComplexType(0.0), v1_);
          using std::copy_n;
          boost::multi::array<ComplexType, 2> v2_(v1_);
          ma::add(ComplexType(1.0), v2_, ComplexType(1.0), vn0_[K], vn0_[K]);
#else
          ma::product(-0.5, Likn, H(Likn), 1.0, vn0_[K]);
#endif
        }
        else
        {
          int QmK = QKtok2[Qm][K];
          boost::multi::array_ref<SPComplexType, 3> Lkin(to_address(LQKikn[Qm][QK].origin()),
                                                         {nmo_max, nmo_max, nchol_max});
          boost::multi::array<SPComplexType, 3> buff3D({nmo_max, nmo_max, nchol_max});
          using ma::conj;
          for (int i = 0; i < nmo_max; i++)
            for (int k = 0; k < nmo_max; k++)
              for (int n = 0; n < nchol_max; n++)
                buff3D[i][k][n] = ma::conj(Lkin[k][i][n]);
          boost::multi::array_ref<SPComplexType, 2> L_(to_address(buff3D.origin()), {nmo_max, nmo_max * nchol_max});
          using ma::H;
#if defined(MIXED_PRECISION)
          boost::multi::array<SPComplexType, 2> v1_({nmo_max, nmo_max});
          ma::product(SPComplexType(-0.5), L_, H(L_), SPComplexType(0.0), v1_);
          boost::multi::array<ComplexType, 2> v2_(v1_);
          ma::add(ComplexType(1.0), v2_, ComplexType(1.0), vn0_[K], vn0_[K]);
#else
          ma::product(-0.5, L_, H(L_), 1.0, vn0_[K]);
#endif
        }
      }
    }
  }
  TG.Global().all_reduce_in_place_n(vn0_.origin(), vn0_.num_elements(), std::plus<>());
  copy_n(vn0_.origin(), vn0_.num_elements(), vn0.origin());

  //  TG.Node().barrier();
  //  if(TG.Node().root())
  //    TG.Cores().all_reduce_in_place_n(to_address(vn0.origin()),vn0.num_elements(),std::plus<>());

  if (TG.Node().root())
  {
    dump.pop();
    dump.close();
  }

  int global_ncvecs = 2 * std::accumulate(nchol_per_kp.begin(), nchol_per_kp.end(), 0);

  std::vector<RealType> gQ(nkpts);
  if (nsampleQ > 0)
  {
    app_log() << " Sampling EXX energy using distribution over Q vector obtained from "
              << " trial energy. \n";

    if (npol == 2)
      APP_ABORT("Error: nsampleQ>0 not yet implemented for noncollinear.\n\n\n");

    RealType scl = (type == CLOSED ? 2.0 : 1.0);
    size_t nqk   = 0;
    for (int Q = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Qmap[Q] < 0)
        continue;
      int Qm = kminus[Q];
      for (int Ka = 0; Ka < nkpts; ++Ka)
      {
        int Kk = QKtok2[Q][Ka];
        int Kb = Kk;
        int Kl = QKtok2[Qm][Kb];
        if ((Ka != Kl) || (Kb != Kk))
          APP_ABORT(" Error: Problems with EXX.\n");
        if ((nqk++) % Qcomm.size() == Qcomm.rank())
        {
          int nchol = nchol_per_kp[Q];
          int nl    = nmo_per_kp[Kl];
          int nb    = nocc_per_kp[0][Kb];
          int nk    = nmo_per_kp[Kk];
          int na    = nocc_per_kp[0][Ka];

          if (na == 0 || nb == 0)
            continue;

          SpMatrix_ref Lank(to_address(LQKank[Q][Ka].origin()), {na * nchol_max, nmo_max});
          auto bnl_ptr(to_address(LQKank[Qm][Kb].origin()));
          if (Qmap[Q] > 0)
            bnl_ptr = to_address(LQKbnl[Qmap[Q] - 1][Kb].origin());
          SpMatrix_ref Lbnl(bnl_ptr, {nb * nchol_max, nmo_max});

          SpMatrix Tban({nb, na * nchol_max});
          Sp3Tensor_ref T3ban(Tban.origin(), {nb, na, nchol_max});
          SpMatrix Tabn({na, nb * nchol_max});
          Sp3Tensor_ref T3abn(Tabn.origin(), {na, nb, nchol_max});

          auto Gal = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[0], Ka, npol == 2);
          auto Gbk = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[0], Kb, npol == 2);
          for (int a = 0; a < na; ++a)
            for (int l = 0; l < nl; ++l)
              Gal[a][l] = ma::conj(Gal[a][l]);
          for (int b = 0; b < nb; ++b)
            for (int k = 0; k < nk; ++k)
              Gbk[b][k] = ma::conj(Gbk[b][k]);

          ma::product(Gal, ma::T(Lbnl), Tabn);
          ma::product(Gbk, ma::T(Lank), Tban);

          ComplexType E_(0.0);
          for (int a = 0; a < na; ++a)
            for (int b = 0; b < nb; ++b)
              E_ += ma::dot(T3abn[a][b], T3ban[b][a]);
          gQ[Q] -= scl * 0.5 * real(E_);
        }
        if (type == COLLINEAR)
        {
          APP_ABORT(" Finish UHF.\n ");
        }
      }
    }
    TG.Global().all_reduce_in_place_n(gQ.begin(), nkpts, std::plus<>());
    RealType E_ = std::accumulate(gQ.begin(), gQ.end(), RealType(0.0));
    for (auto& v : gQ)
      v /= E_;
    app_log() << " EXX: " << E_ << std::endl;
    for (auto v : gQ)
    {
      if (v < 0.0)
        APP_ABORT(" Error: g(Q) < 0.0, implement shift to g(Q). \n")
    }
  }

  //  if(write_hdf)
  //    writeKP3IndexFactorization(hdf_restart,type,NMO,NAEA,NAEB,TGprop,TGwfn,
  //                            nmo_per_kp,nchol_per_kp,kminus,QKtok2,H1,LQKikn,
  //                            vn0,nsampleQ,gQ,E0,global_ncvecs);

  if (ooc == "yes" || ooc == "true")
  {
    return HamiltonianOperations(
        KP3IndexFactorization_batched<shmSpMatrix>(type, TG, std::move(nmo_per_kp), std::move(nchol_per_kp),
                                                   std::move(kminus), std::move(nocc_per_kp), std::move(QKtok2),
                                                   std::move(H1), std::move(haj), std::move(LQKikn), std::move(LQKank),
                                                   std::move(LQKakn), std::move(LQKbnl), std::move(LQKbln),
                                                   std::move(Qmap), std::move(vn0), std::move(gQ), nsampleQ, E0,
                                                   device_allocator<ComplexType>{}, global_origin, global_ncvecs,
                                                   memory));
  }
  else
  {
    using devSpMatrix = boost::multi::array<SPComplexType, 2, device_allocator<SPComplexType>>;
    return HamiltonianOperations(
        KP3IndexFactorization_batched<devSpMatrix>(type, TG, std::move(nmo_per_kp), std::move(nchol_per_kp),
                                                   std::move(kminus), std::move(nocc_per_kp), std::move(QKtok2),
                                                   std::move(H1), std::move(haj), std::move(LQKikn), std::move(LQKank),
                                                   std::move(LQKakn), std::move(LQKbnl), std::move(LQKbln),
                                                   std::move(Qmap), std::move(vn0), std::move(gQ), nsampleQ, E0,
                                                   device_allocator<ComplexType>{}, global_origin, global_ncvecs,
                                                   memory));
  }
}

#else

HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations(bool pureSD,
                                                                        bool addCoulomb,
                                                                        WALKER_TYPES type,
                                                                        std::vector<PsiT_Matrix>& PsiT,
                                                                        double cutvn,
                                                                        double cutv2,
                                                                        TaskGroup_& TGprop,
                                                                        TaskGroup_& TGwfn,
                                                                        hdf_archive& hdf_restart)
{
  APP_ABORT(" Error: KPFactorizedHamiltonian requires complex build (QMC_COMPLEX=1). \n");
  return HamiltonianOperations();
}

HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations_shared(bool pureSD,
                                                                               bool addCoulomb,
                                                                               WALKER_TYPES type,
                                                                               std::vector<PsiT_Matrix>& PsiT,
                                                                               double cutvn,
                                                                               double cutv2,
                                                                               TaskGroup_& TGprop,
                                                                               TaskGroup_& TGwfn,
                                                                               hdf_archive& hdf_restart)
{
  APP_ABORT(" Error: KPFactorizedHamiltonian requires complex build (QMC_COMPLEX=1). \n");
  return HamiltonianOperations();
}

HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations_batched(bool pureSD,
                                                                                bool addCoulomb,
                                                                                WALKER_TYPES type,
                                                                                std::vector<PsiT_Matrix>& PsiT,
                                                                                double cutvn,
                                                                                double cutv2,
                                                                                TaskGroup_& TGprop,
                                                                                TaskGroup_& TGwfn,
                                                                                hdf_archive& hdf_restart)
{
  APP_ABORT(" Error: KPFactorizedHamiltonian requires complex build (QMC_COMPLEX=1). \n");
  return HamiltonianOperations();
}


#endif

} // namespace afqmc
} // namespace qmcplusplus
