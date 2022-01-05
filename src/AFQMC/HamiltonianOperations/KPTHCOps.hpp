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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KPTHCOPS_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KPTHCOPS_HPP

#include <fstream>
#include <mutex>

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "mpi3/shm/mutex.hpp"
#include "mpi3/shared_communicator.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "type_traits/complex_help.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"
#include "AFQMC/Wavefunctions/phmsd_helpers.hpp"
#include "AFQMC/Utilities/myTimer.h"

namespace qmcplusplus
{
namespace afqmc
{
class KPTHCOps
{
#if defined(MIXED_PRECISION)
  using SpC = typename to_single_precision<ComplexType>::value_type;
#else
  using SpC = ComplexType;
#endif

  using IVector       = boost::multi::array<int, 1>;
  using CMatrix       = boost::multi::array<ComplexType, 2>;
  using CMatrix_cref  = boost::multi::array_cref<ComplexType, 2>;
  using CMatrix_ref   = boost::multi::array_ref<ComplexType, 2>;
  using SpMatrix_cref = boost::multi::array_cref<SPComplexType, 2>;
  using SpMatrix_ref  = boost::multi::array_ref<SPComplexType, 2>;
  using SpVector_ref  = boost::multi::array_ref<SPComplexType, 1>;
  using C3Tensor      = boost::multi::array<ComplexType, 3>;
  using SpMatrix      = boost::multi::array<SPComplexType, 2>;
  using Sp3Tensor     = boost::multi::array<SPComplexType, 3>;
  using Sp3Tensor_ref = boost::multi::array_ref<SPComplexType, 3>;
  using shmCVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using shmCMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using shmIMatrix    = boost::multi::array<int, 2, shared_allocator<int>>;
  using shmC3Tensor   = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using shmSpVector   = boost::multi::array<SPComplexType, 1, shared_allocator<SPComplexType>>;
  using shmSpMatrix   = boost::multi::array<SPComplexType, 2, shared_allocator<SPComplexType>>;
  using shmSp3Tensor  = boost::multi::array<SPComplexType, 3, shared_allocator<SPComplexType>>;
  using communicator  = boost::mpi3::shared_communicator;
  using shared_mutex  = boost::mpi3::shm::mutex;
  using this_t        = KPTHCOps;

public:
  static const HamiltonianTypes HamOpType = KPTHC;
  HamiltonianTypes getHamType() const { return HamOpType; }

  /*
     * NAOA/NAOB stands for number of active orbitals alpha/beta (instead of active electrons)
     */
  KPTHCOps(communicator& c_,
           WALKER_TYPES type,
           IVector&& nopk_,
           IVector&& ncholpQ_,
           IVector&& kminus_,
           shmIMatrix&& nelpk_,
           shmIMatrix&& QKToK2_,
           shmIMatrix&& QKToG_,
           shmC3Tensor&& h1_,
           shmCMatrix&& haj_,
           std::vector<shmSpMatrix>&& rotmuv_,
           shmCMatrix&& rotpiu_,
           std::vector<shmSpMatrix>&& rotpau_,
           std::vector<shmSpMatrix>&& lun_,
           shmSpMatrix&& piu_,
           std::vector<shmCMatrix>&& pau_,
           shmC3Tensor&& vn0_,
           ValueType e0_,
           int gncv,
           std::vector<shmSpMatrix>&& vik,
           std::vector<shmSpMatrix>&& vak,
           bool verbose = false)
      : comm(std::addressof(c_)),
        walker_type(type),
        global_nCV(gncv),
        H1(std::move(h1_)),
        haj(std::move(haj_)),
        nopk(std::move(nopk_)),
        ncholpQ(std::move(ncholpQ_)),
        kminus(std::move(kminus_)),
        nelpk(std::move(nelpk_)),
        QKToK2(std::move(QKToK2_)),
        QKToG(std::move(QKToG_)),
        rotMuv(std::move(rotmuv_)),
        rotPiu(std::move(rotpiu_)),
        rotcPua(std::move(rotpau_)),
        LQGun(std::move(lun_)),
        Piu(std::move(piu_)),
        cPua(std::move(pau_)),
        vn0(std::move(vn0_)),
        E0(e0_),
        KK2Q({nopk.size(), nopk.size()}, shared_allocator<SPComplexType>{c_}),
        SM_TMats({1, 1}, shared_allocator<SPComplexType>{c_}),
        TMats({1, 1}),
        mutex(0),
        LQKikn(std::move(vik)),
        LQKank(std::move(vak))
  {
    nGpk.resize(nopk.size());
    nG2pk.resize(nopk.size());
    int nGG(0);
    for (int Q = 0; Q < nopk.size(); Q++)
    {
      nGpk[Q]  = *std::max_element(QKToG[Q].begin(), QKToG[Q].end()) + 1;
      nG2pk[Q] = nGG;
      nGG += nGpk[Q] * nGpk[Q];
    }
    app_log() << "  KPTHCOps: Found " << nGG << " {Q,dq} pairs. \n";
    int nu    = Piu.size(1);
    int rotnu = rotPiu.size(1);
    int nkpts = nopk.size();
    local_nCV = std::accumulate(ncholpQ.begin(), ncholpQ.end(), 0);
    mutex.reserve(ncholpQ.size());
    for (int nQ = 0; nQ < ncholpQ.size(); nQ++)
      mutex.emplace_back(std::make_unique<shared_mutex>(*comm));
    // transposing until Fionn fixes ordering
    if (comm->root())
    {
      for (int Q = 0; Q < nkpts; Q++)
      {
        int nG = nGpk[Q];
        boost::multi::array<SPComplexType, 2> T({nG * rotnu, nG * rotnu});
        std::copy_n(to_address(rotMuv[Q].origin()), rotMuv[Q].num_elements(), T.origin());
        boost::multi::array_ref<SPComplexType, 4> Muv(to_address(rotMuv[Q].origin()), {nG, nG, rotnu, rotnu});
        boost::multi::array_ref<SPComplexType, 4> T_(T.origin(), {nG, rotnu, nG, rotnu});
        for (int G1 = 0; G1 < nG; G1++)
          for (int G2 = 0; G2 < nG; G2++)
            for (int u = 0; u < rotnu; u++)
              for (int v = 0; v < rotnu; v++)
                Muv[G1][G2][u][v] = T_[G1][u][G2][v];
      }
    }
    /*
      // storing full matrix for now, switch to upper triangular later
      Muv.reserve(nGG);
      for(int i=0; i<nGG; i++)
        Muv.emplace_back(shmSpMatrix({nu,nu},shared_allocator<SPComplexType>{c_}));
      comm->barrier();
      // split over  Global later
      for(int Q=0, nq=0; Q<nkpts; ++Q) {
        int nG = nGpk[Q];
        for(int G1=0; G1<nG; ++G1) {
          //for(int G2=G1; G2<nG; ++G2) {
          for(int G2=0; G2<nG; ++G2, ++nq) {
            if(nq%comm->size() == comm->rank())
              ma::product(LQGun[Q].sliced(nu*G1,nu*(G1+1)),ma::H(LQGun[Q].sliced(nu*G2,nu*(G2+1))),Muv[nq]);
          }
        }
      }
*/
    if (comm->root())
    {
      for (int KI = 0; KI < nkpts; KI++)
        for (int KK = 0; KK < nkpts; KK++)
        {
          KK2Q[KI][KK] = -1;
          for (int Q = 0; Q < nkpts; Q++)
            if (QKToK2[Q][KI] == KK)
            {
              KK2Q[KI][KK] = Q;
              break;
            }
          assert(KK2Q[KI][KK] >= 0);
        }
    }
    comm->barrier();
  }

  ~KPTHCOps() {}

  KPTHCOps(KPTHCOps const& other) = delete;
  KPTHCOps& operator=(KPTHCOps const& other) = delete;

  KPTHCOps(KPTHCOps&& other) = default;
  KPTHCOps& operator=(KPTHCOps&& other) = default;

  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(TaskGroup_& TG,
                                                                 boost::multi::array<ComplexType, 1> const& vMF)
  {
    int nkpts = nopk.size();
    int NMO   = std::accumulate(nopk.begin(), nopk.end(), 0);
    // in non-collinear case with SO, keep SO matrix here and add it
    // for now, stay collinear
    boost::multi::array<ComplexType, 2> P1({NMO, NMO});

    // making a copy of vMF since it will be modified
    set_shm_buffer(vMF.size(0));
    boost::multi::array_ref<ComplexType, 1> vMF_(to_address(SM_TMats.origin()), {vMF.size(0)});

    boost::multi::array_ref<ComplexType, 1> P1D(to_address(P1.origin()), {NMO * NMO});
    std::fill_n(P1D.origin(), P1D.num_elements(), ComplexType(0));
    vHS(vMF_, P1D);
    TG.TG().all_reduce_in_place_n(P1D.origin(), P1D.num_elements(), std::plus<>());

    // add H1 + vn0 and symmetrize
    using ma::conj;

    for (int K = 0, nk0 = 0; K < nkpts; ++K)
    {
      for (int i = 0, I = nk0; i < nopk[K]; i++, I++)
      {
        P1[I][I] += H1[K][i][i] + vn0[K][i][i];
        for (int j = i + 1, J = I + 1; j < nopk[K]; j++, J++)
        {
          P1[I][J] += H1[K][i][j] + vn0[K][i][j];
          P1[J][I] += H1[K][j][i] + vn0[K][j][i];
          // This is really cutoff dependent!!!
          if (std::abs(P1[I][J] - conj(P1[J][I])) > 1e-6)
          {
            app_error() << " WARNING in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
            app_error() << I << " " << J << " " << P1[I][J] << " " << P1[j][i] << " " << H1[K][i][j] << " "
                        << H1[K][j][i] << " " << vn0[K][i][j] << " " << vn0[K][j][i] << std::endl;
            //APP_ABORT("Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n");
          }
          P1[I][J] = 0.5 * (P1[I][J] + conj(P1[J][I]));
          P1[J][I] = conj(P1[I][J]);
        }
      }
      nk0 += nopk[K];
    }
    return P1;
  }

  template<class Mat, class MatB>
  void energy(Mat&& E, MatB const& G, int k, bool addH1 = true, bool addEJ = true, bool addEXX = true)
  {
    MatB* Kr(nullptr);
    MatB* Kl(nullptr);
    energy(E, G, k, Kl, Kr, addH1, addEJ, addEXX);
  }

  // Kl and Kr must be in shared memory for this to work correctly
  template<class Mat, class MatB, class MatC, class MatD>
  void energy(Mat&& E,
              MatB const& Gc,
              int nd,
              MatC* KEleft,
              MatD* KEright,
              bool addH1  = true,
              bool addEJ  = true,
              bool addEXX = true)
  {
    using ma::conj;
    using ma::H;
    using ma::T;
    if (nd > 0)
      APP_ABORT(" Error: KPTHC not yet implemented for multiple references.\n");
    static_assert(E.dimensionality == 2, "Wrong dimensionality");
    static_assert(Gc.dimensionality == 2, "Wrong dimensionality");
    assert(E.size(0) == Gc.size(0));
    assert(E.size(1) == 3);
    assert(nd >= 0 && nd < nelpk.size());

    SPComplexType one(1.0);
    int nkpts     = nopk.size();
    int nwalk     = Gc.size(0);
    int nspin     = (walker_type == COLLINEAR ? 2 : 1);
    int nmo_tot   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int nmo_max   = *std::max_element(nopk.begin(), nopk.end());
    int nocca_tot = std::accumulate(nelpk[nd].begin(), nelpk[nd].begin() + nkpts, 0);
    int nG_max    = *std::max_element(nGpk.begin(), nGpk.end());
    int nocca_max = *std::max_element(nelpk[nd].begin(), nelpk[nd].begin() + nkpts);
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());
    int noccb_tot = 0;
    if (walker_type == COLLINEAR)
      noccb_tot = std::accumulate(nelpk[nd].begin() + nkpts, nelpk[nd].begin() + 2 * nkpts, 0);

    for (int n = 0; n < nwalk; n++)
      std::fill_n(E[n].origin(), 3, ComplexType(0.));

    assert(Gc.num_elements() == nwalk * (nocca_tot + noccb_tot) * nmo_tot);
    boost::multi::array_cref<ComplexType, 3> G3Da(to_address(Gc.origin()), {nwalk, nocca_tot, nmo_tot});
    boost::multi::array_cref<ComplexType, 3> G3Db(to_address(Gc.origin()) + G3Da.num_elements() * (nspin - 1),
                                                  {nwalk, noccb_tot, nmo_tot});

    if (addH1)
    {
      int na = 0, nj = 0, nb = 0;
      for (int n = 0; n < nwalk; n++)
        E[n][0] = E0;
      for (int K = 0; K < nkpts; ++K)
      {
        boost::multi::array_ref<ComplexType, 2> haj_K(to_address(haj[nd * nkpts + K].origin()),
                                                      {nelpk[nd][K], nopk[K]});
        for (int n = 0; n < nwalk; ++n)
        {
          ComplexType E_(0.0);
          for (int a = 0; a < nelpk[nd][K]; ++a)
            //E_ += ma::dot(G3Da[n][na+a]({nj,nj+nopk[K]}),haj_K[a]);
            E_ += ma::dot(G3Da[n][na + a].sliced(nj, nj + nopk[K]), haj_K[a]);
          E[n][0] += E_;
        }
        na += nelpk[nd][K];
        if (walker_type == COLLINEAR)
        {
          boost::multi::array_ref<ComplexType, 2> haj_Kb(haj_K.origin() + haj_K.num_elements(),
                                                         {nelpk[nd][nkpts + K], nopk[K]});
          for (int n = 0; n < nwalk; ++n)
          {
            ComplexType E_(0.0);
            for (int a = 0; a < nelpk[nd][nkpts + K]; ++a)
              E_ += ma::dot(G3Db[n][nb + a]({nj, nj + nopk[K]}), haj_Kb[a]);
            E[n][0] += E_;
          }
          nb += nelpk[nd][nkpts + K];
        }
        nj += nopk[K];
      }
    }

    if (not addEXX and not addEJ)
      return;

    RealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
    SPComplexType *Krptr, *Klptr;
    int getKr           = KEright != nullptr;
    int getKl           = KEleft != nullptr;
    int rotnu           = rotPiu.size(1);
    int nGu             = std::accumulate(nGpk.begin(), nGpk.end(), 0) * rotnu;
    size_t memory_needs = 0;
    if (addEXX)
      memory_needs += nkpts * nkpts * rotnu * rotnu;
    if (addEJ)
    {
      if (not getKr)
        memory_needs += nwalk * nGu;
      if (not getKl)
        memory_needs += nwalk * nGu;
    }
    set_shm_buffer(memory_needs);
    size_t cnt = 0;

    // messy
    long Knr = 0, Knc = 0;
    if (addEJ)
    {
      Knr = nwalk;
      Knc = nGu;
      if (getKr)
      {
        assert(KEright->size(0) == nwalk && KEright->size(1) == nGu);
        assert(KEright->stride(0) == KEright->size(1));
        Krptr = to_address(KEright->origin());
      }
      else
      {
        Krptr = to_address(SM_TMats.origin()) + cnt;
        cnt += nwalk * nGu;
      }
      if (getKl)
      {
        assert(KEleft->size(0) == nwalk && KEleft->size(1) == nGu);
        assert(KEleft->stride(0) == KEleft->size(1));
        Klptr = to_address(KEleft->origin());
      }
      else
      {
        Klptr = to_address(SM_TMats.origin()) + cnt;
        cnt += nwalk * nGu;
      }
      if (comm->root())
        std::fill_n(Krptr, Knr * Knc, SPComplexType(0.0));
      if (comm->root())
        std::fill_n(Klptr, Knr * Knc, SPComplexType(0.0));
    }
    else if (getKr or getKl)
    {
      APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
    }
    SpMatrix_ref Kl(Klptr, {Knr, Knc});
    SpMatrix_ref Kr(Krptr, {Knr, Knc});
    comm->barrier();

    Timer.reset("T0");
    Timer.reset("T1");
    Timer.reset("T2");
    Timer.reset("T3");
    if (addEXX)
    {
      size_t local_memory_needs = nocca_tot * rotnu + rotnu;
      if (TMats.num_elements() < local_memory_needs)
        TMats.reextent({local_memory_needs, 1});

      // Fuv[k1][k2][u][v] = sum_a_l rotcPua[u][k1][a] * G[k1][a][k2][l] rotPiu[k2][l][v]
      boost::multi::array_ref<SPComplexType, 4> Fuv(to_address(SM_TMats.origin()) + cnt, {nkpts, nkpts, rotnu, rotnu});
      cnt += Fuv.num_elements();

      size_t cnt_local = 0;
      boost::multi::array_ref<SPComplexType, 2> TAv(TMats.origin() + cnt_local, {nocca_tot, rotnu});
      cnt_local += TAv.num_elements();

      // avoiding vectors for now
      SpVector_ref Zu(TMats.origin() + cnt_local, {rotnu});
      cnt_local += Zu.num_elements();
      std::fill_n(Zu.origin(), Zu.num_elements(), SPComplexType(0.0));

      for (int n = 0; n < nwalk; ++n)
      {
        size_t nqk = 1; // start count at 1 to "offset" the calcuation of E1 done at root
        Timer.start("T0");
        for (int Kl = 0, nl0 = 0; Kl < nkpts; ++Kl)
        {
          if ((nqk++) % comm->size() == comm->rank())
          {
            ma::product(G3Da[n]({0, nocca_tot}, {nl0, nl0 + nopk[Kl]}), rotPiu({nl0, nl0 + nopk[Kl]}, {0, rotnu}), TAv);
            // write with BatchedGEMM later
            for (int Ka = 0, na0 = 0; Ka < nkpts; ++Ka)
            {
              ma::product(rotcPua[nd]({0, rotnu}, {na0, na0 + nelpk[nd][Ka]}),
                          TAv({na0, na0 + nelpk[nd][Ka]}, {0, rotnu}), Fuv[Ka][Kl]);
              na0 += nelpk[nd][Ka];
            }
          }
          nl0 += nopk[Kl];
        }
        comm->barrier();
        Timer.stop("T0");

        Timer.start("T2");
        // NOT OPTIMAL: FIX!!!
        if (addEJ)
        {
          nqk = 0;
          for (int Q = 0, nqGa = 0; Q < nkpts; ++Q)
          { // momentum conservation index
            int nG = nGpk[Q];
            for (int Ga = 0; Ga < nG; ++Ga, nqGa += rotnu)
            {
              if ((nqk++) % comm->size() == comm->rank())
              {
                boost::multi::array_ref<SPComplexType, 4> Muv(to_address(rotMuv[Q].origin()), {nG, nG, rotnu, rotnu});
                auto&& KlQa(Kl[n].sliced(nqGa, nqGa + rotnu));
                auto&& KrQa(Kr[n].sliced(nqGa, nqGa + rotnu));
                std::fill_n(KlQa.origin(), KlQa.num_elements(), SPComplexType(0.0));
                std::fill_n(KrQa.origin(), KrQa.num_elements(), SPComplexType(0.0));
                // Kl(n,Q,Ga,u) = sum_K_in_Ga F(K,Q[K],u,u)
                for (int Ka = 0; Ka < nkpts; ++Ka)
                {
                  if (QKToG[Q][Ka] != Ga)
                    continue;
                  int Kk = QKToK2[Q][Ka];
                  auto f_(Fuv[Ka][Kk].origin());
                  auto ku_(KlQa.origin());
                  for (int u = 0; u < rotnu; u++, ku_++, f_ += (rotnu + 1))
                    (*ku_) += (*f_);
                }
                // Kr(n,Q,Ga,u) = sum_Gl sum_v M(Q,Ga,Gl)(u,v) sum_K_in_Gl F(Q[K],K,u,u)
                for (int Gl = 0; Gl < nG; ++Gl)
                {
                  std::fill_n(Zu.origin(), Zu.num_elements(), SPComplexType(0.0));
                  for (int Kl = 0; Kl < nkpts; ++Kl)
                  {
                    if (QKToG[Q][Kl] != Gl)
                      continue;
                    int Kb = QKToK2[Q][Kl];
                    auto f_(Fuv[Kb][Kl].origin());
                    auto zu_(Zu.origin());
                    for (int u = 0; u < rotnu; u++, zu_++, f_ += (rotnu + 1))
                      (*zu_) += (*f_);
                  }
                  ma::product(one, Muv[Ga][Gl], Zu, one, KrQa);
                }
                /*
                  for(int Gb=0; Gb<nG; ++Gb) {
                    std::fill_n(Zu.origin(),Zu.num_elements(),SPComplexType(0.0));
                    for(int Kb=0; Kb<nkpts; ++Kb) {
                      if(QKToG[Q][Kb] != Gb) continue;
                      int Kl = QKToK2[kminus[Q]][Kb];
                      auto f_(Fuv[Kb][Kl].origin());
                      auto zu_(Zu.origin());
                      for(int u=0; u<rotnu; u++, zu_++, f_ += (rotnu+1))
                        (*zu_) += (*f_);
                    }
                    ma::product(one,Muv[Ga][Gb],Zu,one,KrQa);
                  }
*/
              }
            }
          }
        }
        Timer.stop("T2");
        Timer.start("T1");
        // FIX parallelization!!!
        int bsz = 256;
        int nbu = (rotnu + bsz - 1) / bsz;
        nqk     = 0;
        for (int Q = 0; Q < nkpts; ++Q)
        { // momentum conservation index
          int nG = nGpk[Q];
          boost::multi::array_ref<SPComplexType, 4> Muv(to_address(rotMuv[Q].origin()), {nG, nG, rotnu, rotnu});
          for (int K1 = 0; K1 < nkpts; ++K1)
          {
            for (int K2 = 0; K2 < nkpts; ++K2)
            {
              if ((nqk++) % comm->size() == comm->rank())
              {
                int QK1 = QKToK2[Q][K1];
                int QK2 = QKToK2[Q][K2];
                // EXX += sum_u_v Muv[u][v] * Fuv[K1][K2][u][v] * Fuv[QK2][QK1][v][u]
                ComplexType E_(0.0);
                //                  int nq = nG2pk[Q] + QKToG[Q][K1]*nG + QKToG[Q][K2];
                auto F1_(to_address(Fuv[K1][K2].origin()));
                auto muv_(Muv[QKToG[Q][K1]][QKToG[Q][K1]].origin());
                auto F2_(to_address(Fuv[QK2][QK1].origin()));
                for (int u = 0; u < rotnu; ++u, ++F2_)
                {
                  auto Fv_(F2_);
                  for (int v = 0; v < rotnu; ++v, ++muv_, ++F1_, Fv_ += rotnu)
                    E_ += (*F1_) * (*muv_) * (*Fv_);
                }
                E[n][1] -= 0.5 * scl * E_;
              }
            }
          }
        }
        Timer.stop("T1");
        if (walker_type == COLLINEAR)
        {
          APP_ABORT("Error: Finish UHF in KPTHC.\n");
        }
      }
    }

    Timer.start("T2");
    if (addEJ)
    {
      if (not addEXX)
      {
        // calculate Kr
        APP_ABORT(" Error: Finish addEJ and not addEXX");
      }
      RealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
      size_t nqk   = 0; // start count at 1 to "offset" the calcuation of E1 done at root
      for (int n = 0; n < nwalk; ++n)
      {
        for (int Q = 0; Q < nkpts; ++Q)
        { // momentum conservation index
          if ((nqk++) % comm->size() == comm->rank())
          {
            int nc0 = std::accumulate(nGpk.begin(), nGpk.begin() + Q, 0) * rotnu;
            E[n][2] +=
                0.5 * scl * scl * ma::dot(Kl[n]({nc0, nc0 + rotnu * nGpk[Q]}), Kr[n]({nc0, nc0 + rotnu * nGpk[Q]}));
          }
        }
      }
    }
    Timer.stop("T2");
    app_log() << " E time: " << Timer.total("T0") << " " << Timer.total("T1") << " " << Timer.total("T2") << "\n";
  }

  template<class MatE, class MatO, class MatG, class MatQ, class MatB, class index_aos>
  void fast_energy(MatE&& E,
                   MatO&& Ov,
                   MatG const& GrefA,
                   MatG const& GrefB,
                   MatQ const& QQ0A,
                   MatQ const& QQ0B,
                   MatB&& Qwork,
                   ph_excitations<int, ComplexType> const& abij,
                   std::array<index_aos, 2> const& det_couplings)
  {}

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    boost::multi::array_ref<ComplexType, 2> X_(to_address(X.origin()), {X.size(0), 1});
    boost::multi::array_ref<ComplexType, 2> v_(to_address(v.origin()), {1, v.size(0)});
    vHS(X_, v_, a, c);
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    int nkpts = nopk.size();
    int nu    = Piu.size(1);
    int nwalk = X.size(1);
    assert(v.size(0) == nwalk);
    int nspin     = (walker_type == COLLINEAR ? 2 : 1);
    int nmo_tot   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int nmo_max   = *std::max_element(nopk.begin(), nopk.end());
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());
    assert(X.num_elements() == nwalk * 2 * local_nCV);
    assert(v.num_elements() == nwalk * nmo_tot * nmo_tot);
    SPComplexType one(1.0, 0.0);
    SPComplexType im(0.0, 1.0);

    size_t local_memory_needs = nmo_max * nmo_max * nwalk + nwalk * nu * (nmo_max + 1);
    if (TMats.num_elements() < local_memory_needs)
      TMats.reextent({local_memory_needs, 1});
    size_t cnt = 0;
    SpMatrix_ref Twu(TMats.origin(), {nwalk, nu});
    cnt += Twu.num_elements();
    auto vik_ptr(TMats.origin() + cnt);
    cnt += nwalk * nmo_max * nmo_max;
    auto Qniu_ptr(TMats.origin() + cnt);
    cnt += nwalk * nu * nmo_max;

    Sp3Tensor_ref v3D(to_address(v.origin()), {nwalk, nmo_tot, nmo_tot});

    // "rotate" X
    //  XIJ = 0.5*a*(Xn+ -i*Xn-), XJI = 0.5*a*(Xn+ +i*Xn-)
    for (int Q = 0, nq = 0; Q < nkpts; ++Q)
    {
      int nc0, ncN;
      std::tie(nc0, ncN) = FairDivideBoundary(comm->rank(), ncholpQ[Q], comm->size());
      auto Xnp           = to_address(X[nq + nc0].origin());
      auto Xnm           = to_address(X[nq + ncholpQ[Q] + nc0].origin());
      for (int n = nc0; n < ncN; ++n)
      {
        for (int nw = 0; nw < nwalk; ++nw, ++Xnp, ++Xnm)
        {
          ComplexType Xnp_ = 0.5 * a * ((*Xnp) - im * (*Xnm));
          *Xnm             = 0.5 * a * ((*Xnp) + im * (*Xnm));
          *Xnp             = Xnp_;
        }
      }
      nq += 2 * ncholpQ[Q];
    }
    // scale v by 'c': assuming contiguous data
    {
      size_t i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(size_t(comm->rank()), size_t(v.num_elements()), size_t(comm->size()));
      auto v_          = to_address(v.origin()) + i0;
      for (size_t i = i0; i < iN; ++i, ++v_)
        *v_ *= c;
    }
    comm->barrier();
    using ma::H;
    using ma::T;
    size_t nqk = 0;
    for (int Q = 0, nc0 = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Q % comm->size() == comm->rank())
      {
        int nchol = ncholpQ[Q];
        int nG    = nGpk[Q];
        for (int G = 0; G < nG; ++G)
        {
          // combine n+ and n- by maping X({n+,n-},w) -> X(n,{w+,w-})
          ma::product(ma::T(X.sliced(nc0, nc0 + nchol)), ma::T(LQGun[Q].sliced(nu * G, nu * (G + 1))), Twu);
          for (int K = 0; K < nkpts; ++K)
          {
            if (QKToG[Q][K] != G)
              continue;
            int ni  = nopk[K];
            int ni0 = std::accumulate(nopk.begin(), nopk.begin() + K, 0);
            int nk  = nopk[QKToK2[Q][K]];
            int nk0 = std::accumulate(nopk.begin(), nopk.begin() + QKToK2[Q][K], 0);
            SpMatrix_ref vik(vik_ptr, {nwalk * ni, nk});
            Sp3Tensor_ref vik3D(vik_ptr, {nwalk, ni, nk});

            SpMatrix_ref Qniu(Qniu_ptr, {nwalk * ni, nu});

            auto Qniu_(Qniu_ptr);
            for (int n = 0; n < nwalk; ++n)
            {
              auto p_(to_address(Piu[ni0].origin()));
              for (int i = 0; i < ni; ++i)
              {
                auto Tu(Twu[n].origin());
                for (int u = 0; u < nu; ++u, ++p_, ++Qniu_, ++Tu)
                  (*Qniu_) = (*Tu) * conj(*p_);
              }
            }
            // v[nw][i(in K)][k(in Q(K))] += sum_u_n conj(Piu(I,u)) Piu(J,u) LQGun(u,n) X[Q][n+][nw]
            ma::product(Qniu, ma::T(Piu.sliced(nk0, nk0 + nk)), vik);

            // it is possible to add the second half here by calculating the (Q*,K*) pair that maps
            // to the JI term corresponding to this (Q,K) pair. Not doing it for now
            auto vik_(vik.origin());
            for (int n = 0; n < nwalk; n++)
            {
              auto v_(v3D[n][ni0].origin() + nk0);
              for (int i = 0; i < ni; i++, vik_ += nk, v_ += nmo_tot)
                ma::axpy(nk, one, vik_, 1, v_, 1);
            }
          }
        }
      }
      nc0 += 2 * ncholpQ[Q];
    }
    comm->barrier();
    // adding second half. sync here to avoid need for locks
    for (int Q = 0, nc0 = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Q % comm->size() == comm->rank())
      {
        int nchol = ncholpQ[Q];
        int nG    = nGpk[Q];
        for (int G = 0; G < nG; ++G)
        {
          // combine n+ and n- by maping X({n+,n-},w) -> X(n,{w+,w-})
          ma::product(ma::T(X.sliced(nc0 + nchol, nc0 + 2 * nchol)), ma::H(LQGun[Q].sliced(nu * G, nu * (G + 1))), Twu);
          for (int K = 0; K < nkpts; ++K)
          {
            if (QKToG[Q][K] != G)
              continue;
            int ni  = nopk[K];
            int ni0 = std::accumulate(nopk.begin(), nopk.begin() + K, 0);
            int nk  = nopk[QKToK2[Q][K]];
            int nk0 = std::accumulate(nopk.begin(), nopk.begin() + QKToK2[Q][K], 0);
            SpMatrix_ref vki(vik_ptr, {nwalk * nk, ni});
            Sp3Tensor_ref vki3D(vik_ptr, {nwalk, nk, ni});

            SpMatrix_ref Qnku(Qniu_ptr, {nwalk * nk, nu});

            auto Qnku_(Qniu_ptr);
            for (int n = 0; n < nwalk; ++n)
            {
              auto p_(to_address(Piu[nk0].origin()));
              for (int k = 0; k < nk; ++k)
              {
                auto Tu(Twu[n].origin());
                for (int u = 0; u < nu; ++u, ++p_, ++Qnku_, ++Tu)
                  (*Qnku_) = (*Tu) * conj(*p_);
              }
            }
            // v[nw][i(in K)][k(in Q(K))] += sum_u_n conj(Piu(I,u)) Piu(J,u) LQGun(u,n) X[Q][n+][nw]
            ma::product(Qnku, ma::T(Piu.sliced(ni0, ni0 + ni)), vki);

            // it is possible to add the second half here by calculating the (Q*,K*) pair that maps
            // to the JI term corresponding to this (Q,K) pair. Not doing it for now
            auto vki_(vki.origin());
            for (int n = 0; n < nwalk; n++)
            {
              auto v_(v3D[n][nk0].origin() + ni0);
              for (int k = 0; k < nk; k++, vki_ += ni, v_ += nmo_tot)
                ma::axpy(ni, one, vki_, 1, v_, 1);
            }
          }
        }
      }
      nc0 += 2 * ncholpQ[Q];
    }
    comm->barrier();
    // do I need to "rotate" back, can be done if necessary
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vbias(MatA const& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    boost::multi::array_cref<ComplexType, 2> G_(to_address(G.origin()), {1, G.size(0)});
    boost::multi::array_ref<ComplexType, 2> v_(to_address(v.origin()), {v.size(0), 1});
    vbias(G_, v_, a, c, k);
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vbias(MatA const& G, MatB&& v, double a = 1., double c = 0., int nd = 0)
  {
    if (comm->size() > 1)
      APP_ABORT(" Error: only ncore==1 for now. \n");

    int nkpts = nopk.size();
    int nu    = Piu.size(1);
    assert(nd >= 0 && nd < nelpk.size());
    int nwalk = G.size(0);
    assert(v.size(0) == 2 * local_nCV);
    assert(v.size(1) == nwalk);
    int nspin     = (walker_type == COLLINEAR ? 2 : 1);
    int nmo_tot   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int nmo_max   = *std::max_element(nopk.begin(), nopk.end());
    int nocca_tot = std::accumulate(nelpk[nd].begin(), nelpk[nd].begin() + nkpts, 0);
    int nocca_max = *std::max_element(nelpk[nd].begin(), nelpk[nd].begin() + nkpts);
    int noccb_max = nocca_max;
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());
    int noccb_tot = 0;
    if (walker_type == COLLINEAR)
    {
      noccb_tot = std::accumulate(nelpk[nd].begin() + nkpts, nelpk[nd].begin() + 2 * nkpts, 0);
      noccb_max = *std::max_element(nelpk[nd].begin() + nkpts, nelpk[nd].begin() + 2 * nkpts);
    }
    assert(G.num_elements() == nwalk * (nocca_tot + noccb_tot) * nmo_tot);

    RealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
    SPComplexType one(1.0, 0.0);
    SPComplexType zero(0.0, 0.0);
    SPComplexType halfa(0.5 * a * scl, 0.0);
    SPComplexType minusimhalfa(0.0, -0.5 * a * scl);
    SPComplexType imhalfa(0.0, 0.5 * a * scl);

    {
      size_t i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(size_t(comm->rank()), size_t(v.size(0)), size_t(comm->size()));
      for (size_t i = i0; i < iN; ++i)
        ma::scal(c, v[i]);
    }
    comm->barrier();

    // Attempt #1
    size_t local_memory_needs = 2 * nu * nwalk + 2 * nchol_max * nwalk + nu * nocca_tot;
    size_t cnt                = 0;
    if (TMats.num_elements() < local_memory_needs)
      TMats.reextent({local_memory_needs, 1});
    SpMatrix_ref vlocal(TMats.origin() + cnt, {nchol_max, 2 * nwalk});
    cnt += vlocal.num_elements();
    auto Tua_ptr(TMats.origin() + cnt);
    cnt += (nu * nocca_max);
    std::fill_n(vlocal.origin(), vlocal.num_elements(), SPComplexType(0.0));
    SpMatrix_ref Fwu(TMats.origin() + cnt, {2 * nwalk, nu});
    cnt += Fwu.num_elements();
    boost::multi::array_cref<ComplexType, 3> G3Da(to_address(G.origin()), {nwalk, nocca_tot, nmo_tot});
    boost::multi::array_cref<ComplexType, 3> G3Db(to_address(G.origin()) + G3Da.num_elements() * (nspin - 1),
                                                  {nwalk, noccb_tot, nmo_tot});
    /*
Timer.reset("T0");
Timer.reset("T1");
Timer.reset("T2");
Timer.reset("T3");
Timer.reset("T4");
Timer.start("T0");
*/
    size_t nqk = 0;
    auto& cPua_nd(cPua[nd]);
    // for fine grained parallelization, split statically (Q,G) pairs to lead to "even" workload
    for (int Q = 0; Q < nkpts; ++Q)
    {
      if (Q % comm->size() == comm->rank())
      {
        int nG    = nGpk[Q];
        int nchol = ncholpQ[Q];
        std::fill_n(vlocal.origin(), vlocal.num_elements(), SPComplexType(0.0));
        auto&& vloc = vlocal.sliced(0, nchol);
        auto&& v1   = vlocal({0, nchol}, {0, nwalk});
        auto&& v2   = vlocal({0, nchol}, {nwalk, 2 * nwalk});
        for (int G = 0; G < nG; ++G)
        {
          std::fill_n(Fwu.origin(), Fwu.num_elements(), SPComplexType(0.0));
          for (int K = 0; K < nkpts; ++K)
          {
            if (QKToG[Q][K] != G)
              continue;
            int QK = QKToK2[Q][K];
            // (K,QK)
            int na1  = nelpk[nd][K];
            int na01 = std::accumulate(nelpk[nd].begin(), nelpk[nd].begin() + K, 0);
            int nk1  = nopk[QK];
            int nk01 = std::accumulate(nopk.begin(), nopk.begin() + QK, 0);
            auto&& Piu1_(Piu.sliced(nk01, nk01 + nk1));
            SpMatrix_ref Tua1(Tua_ptr, {nu, na1});
            // (QK,K)
            int na2  = nelpk[nd][QK];
            int na02 = std::accumulate(nelpk[nd].begin(), nelpk[nd].begin() + QK, 0);
            int nk2  = nopk[K];
            int nk02 = std::accumulate(nopk.begin(), nopk.begin() + K, 0);
            int nA   = cPua_nd.size(1);
            auto&& Piu2_(Piu.sliced(nk02, nk02 + nk2));
            SpMatrix_ref Tua2(Tua_ptr, {nu, na2});

            auto Fu1(Fwu.origin());
            for (int n = 0; n < nwalk; n++)
            {
              // Tua = sum_k T(Piu(k,u)) T(G[n](a,k))
              //Timer.start("T1");
              ma::product(ma::T(Piu1_), ma::T(G3Da[n]({na01, na01 + na1}, {nk01, nk01 + nk1})), Tua1);
              //Timer.stop("T1");
              //Timer.start("T3");
              // Fwu[w][u] = sum_a cPua(u,a) T(u,a)
              auto Tua1_a(Tua1.origin());
              auto cPua_nd_u(cPua_nd.origin() + na01);
              for (int u = 0; u < nu; u++, ++Fu1, cPua_nd_u += nA)
              {
                auto cPua_nd_a(cPua_nd_u);
                for (int ia = 0; ia < na1; ++ia, ++cPua_nd_a, ++Tua1_a)
                  *(Fu1) += (*cPua_nd_a) * (*Tua1_a);
              }
              //Timer.stop("T3");
            }
            // If LIK_n == conj(LKI_n), then v1(Q) = v2(-Q) and there is no need to calculate both components
            auto Fu2(Fwu.origin() + nwalk * nu);
            for (int n = 0; n < nwalk; n++)
            {
              //Timer.start("T1");
              ma::product(ma::T(Piu2_), ma::T(G3Da[n]({na02, na02 + na2}, {nk02, nk02 + nk2})), Tua2);
              //Timer.stop("T1");
              //Timer.start("T3");
              // Fwu[w][u] = sum_a cPua(u,a) T(u,a)
              auto Tua2_a(Tua2.origin());
              auto cPua_nd_u(cPua_nd.origin() + na02);
              for (int u = 0; u < nu; u++, ++Fu2, cPua_nd_u += nA)
              {
                auto cPua_nd_a(cPua_nd_u);
                for (int ia = 0; ia < na2; ++ia, ++cPua_nd_a, ++Tua2_a)
                  *(Fu2) += conj((*cPua_nd_a) * (*Tua2_a));
              }
              //Timer.stop("T3");
            }
          }
          //Timer.start("T2");
          // v1[n,w] = sum_u T(LQGun) * T(Fwu(w,u))
          ma::product(one, ma::T(LQGun[Q].sliced(G * nu, (G + 1) * nu)), ma::T(Fwu), one, vloc);
          //Timer.stop("T2");

          if (walker_type == COLLINEAR) {}
        }
        // conjugate v2 (use transform?)
        for (int i = 0; i < nchol; ++i)
        {
          auto v2_(v2[i]);
          for (int j = 0; j < nwalk; ++j)
            v2_[j] = conj(v2_[j]);
        }
        int nc0 = 2 * std::accumulate(ncholpQ.begin(), ncholpQ.begin() + Q, 0);
        /*
          for(int i=0; i<nchol; ++i) {
            // v+ = 0.5*a*(v1+v2)
            axpy(nwalk, halfa, v1[i].origin(), 1, v[nc0+i].origin(), 1);
            axpy(nwalk, halfa, v2[i].origin(), 1, v[nc0+i].origin(), 1);
          // v- = -0.5*a*i*(v1-v2)
            axpy(nwalk, minusimhalfa, v1[i].origin(), 1, v[nc0+nchol+i].origin(), 1);
            axpy(nwalk, imhalfa, v2[i].origin(), 1, v[nc0+nchol+i].origin(), 1);
          }
*/
      }
    }
    comm->barrier();
    /*
Timer.stop("T0");
app_log()<<" E time: "
<<Timer.total("T0") <<" "
<<Timer.total("T1") <<" "
<<Timer.total("T2") <<" "
<<Timer.total("T3") <<"\n";
*/
  }

  template<class Mat, class MatB>
  void generalizedFockMatrix(Mat&& G, MatB&& Fp, MatB&& Fm)
  {
    APP_ABORT(" Error: generalizedFockMatrix not implemented for this hamiltonian.\n");
  }

  bool distribution_over_cholesky_vectors() const { return true; }
  int number_of_ke_vectors() const { return std::accumulate(nGpk.begin(), nGpk.end(), 0) * rotPiu.size(1); }
  int local_number_of_cholesky_vectors() const { return 2 * std::accumulate(ncholpQ.begin(), ncholpQ.end(), 0); }
  int global_number_of_cholesky_vectors() const { return global_nCV; }

  // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
  bool transposed_G_for_vbias() const { return true; }
  bool transposed_G_for_E() const { return true; }
  // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
  bool transposed_vHS() const { return true; }

  bool fast_ph_energy() const { return false; }

protected:
  communicator* comm;

  WALKER_TYPES walker_type;

  int global_nCV;
  int local_nCV;

  // bare one body hamiltonian
  shmC3Tensor H1;

  // (potentially half rotated) one body hamiltonian
  shmCMatrix haj;

  // number of orbitals per k-point
  IVector nopk;

  // number of cholesky vectors per Q-point
  IVector ncholpQ;

  // position of (-K) in kp-list for every K
  IVector kminus;

  // number of G per Q-point
  std::vector<int> nGpk;
  std::vector<int> nG2pk;

  // number of electrons per k-point
  // nelpk[ndet][nspin*nkpts]
  shmIMatrix nelpk;

  // maps (Q,K) --> k2
  shmIMatrix QKToK2;

  // maps (Q,K) --> G
  shmIMatrix QKToG;

  /************************************************/
  // Used in the calculation of the energy
  // Coulomb matrix elements of interpolating vectors
  std::vector<shmSpMatrix> rotMuv;

  // Orbitals at interpolating points
  shmSpMatrix rotPiu;

  // Half-rotated Orbitals at interpolating points
  std::vector<shmSpMatrix> rotcPua;
  /************************************************/

  /************************************************/
  // Following 3 used in calculation of vbias and vHS
  std::vector<shmSpMatrix> LQGun;

  // Orbitals at interpolating points
  shmSpMatrix Piu;

  // Half-rotated Orbitals at interpolating points
  std::vector<shmSpMatrix> cPua;
  /************************************************/

  // Muv for energy
  //std::vector<shmSpMatrix> Muv;

  // one-body piece of Hamiltonian factorization
  shmC3Tensor vn0;

  ValueType E0;

  shmIMatrix KK2Q;

  // shared buffer space
  // using matrix since there are issues with vectors
  shmSpMatrix SM_TMats;
  SpMatrix TMats;


  std::vector<std::unique_ptr<shared_mutex>> mutex;

  myTimer Timer;

  void set_shm_buffer(size_t N)
  {
    if (SM_TMats.num_elements() < N)
      SM_TMats.reextent({N, 1});
  }

  boost::multi::array<ComplexType, 3> eloc;

  //Cholesky Tensor Lik[Q][nk][i][k][n]
  std::vector<shmSpMatrix> LQKikn;

  // half-transformed Cholesky tensor
  std::vector<shmSpMatrix> LQKank;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
