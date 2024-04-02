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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_THCOPS_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_THCOPS_HPP

#include <fstream>

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "Utilities/FairDivide.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "mpi3/shared_communicator.hpp"
#include "type_traits/complex_help.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"
#include "AFQMC/Wavefunctions/phmsd_helpers.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
// distribution:  size,  global,  offset
//   - rotMuv:    {rotnmu,grotnmu},{grotnmu,grotnmu},{rotnmu0,0}
//   - rotPiu:    {size_t(NMO),grotnmu},{size_t(NMO),grotnmu},{0,0}
//   - rotcPua    {grotnmu,nel_},{grotnmu,nel_},{0,0}
//   - Piu:       {size_t(NMO),nmu},{size_t(NMO),gnmu},{0,nmu0}
//   - Luv:       {nmu,gnmu},{gnmu,gnmu},{nmu0,0}
//   - cPua       {nmu,nel_},{gnmu,nel_},{nmu0,0}

class THCOps
{
  using communicator = boost::mpi3::shared_communicator;

  // allocators
  // device_allocator for local work space
  // localTG_allocator for shared work space
  // node_allocator for fixed arrays, e.g. Luv, Piu, ...

  template<class T>
  using device_alloc_type = DeviceBufferManager::template allocator_t<T>;
  template<class T>
  using shm_alloc_type = LocalTGBufferManager::template allocator_t<T>;

  // pointers
  using pointer    = typename device_alloc_type<ComplexType>::pointer;
  using sp_pointer = typename device_alloc_type<SPComplexType>::pointer;

  using const_pointer    = typename device_allocator<ComplexType>::const_pointer;
  using const_sp_pointer = typename device_allocator<SPComplexType>::const_pointer;

  template<class U, int N>
  using Array = boost::multi::static_array<U, N, device_alloc_type<U>>;
  template<class U, int N>
  using Array_ref = boost::multi::array_ref<U, N, typename device_alloc_type<U>::pointer>;
  template<class U, int N>
  using Array_cref = boost::multi::array_cref<U, N, typename device_allocator<U>::const_pointer>;
  // arrays on shared work space
  // remember that this is device memory when built with accelerator support
  template<class U, int N>
  using ShmArray = boost::multi::static_array<U, N, shm_alloc_type<U>>;

  // arrays on node allocator, for fixed arrays, e.g. Luv, Piu, ...
  // remember that this is device memory when built with accelerator support
  template<class U, int N>
  using nodeArray = boost::multi::array<U, N, node_allocator<U>>;

  // host array on shared memory
  using mpi3VMatrix   = boost::multi::array<ValueType, 2, shared_allocator<ValueType>>;
  using mpi3CMatrix   = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using mpi3SPVMatrix = boost::multi::array<SPValueType, 2, shared_allocator<SPValueType>>;
  using mpi3SPCMatrix = boost::multi::array<SPComplexType, 2, shared_allocator<SPComplexType>>;

public:
  static const HamiltonianTypes HamOpType = THC;
  HamiltonianTypes getHamType() const { return HamOpType; }

  /*
     * nup/ndown stands for number of active orbitals alpha/beta (instead of active electrons)
     */
  THCOps(communicator& c_,
         int nmo_,
         int naoa_,
         int naob_,
         WALKER_TYPES type,
         int nmu0_,
         int rotnmu0_,
         mpi3CMatrix&& hij_,
         mpi3CMatrix&& h1,
         mpi3SPVMatrix&& rotmuv_,
         mpi3SPVMatrix&& rotpiu_,
         std::vector<mpi3SPCMatrix>&& rotpau_,
         mpi3SPVMatrix&& luv_,
         mpi3SPVMatrix&& piu_,
         std::vector<mpi3SPCMatrix>&& pau_,
         mpi3CMatrix&& v0_,
         ValueType e0_,
         bool verbose = false)
      : comm(std::addressof(c_)),
        device_buffer_manager(),
        shm_buffer_manager(),
        NMO(nmo_),
        nup(naoa_),
        ndown(naob_),
        nelec{nup, ndown},
        nmu0(nmu0_),
        gnmu(0),
        rotnmu0(rotnmu0_),
        grotnmu(0),
        walker_type(type),
        hij(std::move(hij_)),
        haj(std::move(h1)),
        rotMuv(std::move(rotmuv_)),
        rotPiu(std::move(rotpiu_)),
        rotcPua(move_vector<nodeArray<SPComplexType, 2>>(std::move(rotpau_))),
        Luv(std::move(luv_)),
        Piu(std::move(piu_)),
        cPua(move_vector<nodeArray<SPComplexType, 2>>(std::move(pau_))),
        vn0(std::move(v0_)),
        E0(e0_)
  {
    gnmu    = std::get<1>(Luv.sizes());
    grotnmu = std::get<1>(rotMuv.sizes());
    if (haj.size() > 1)
      APP_ABORT(" Error: THC not yet implemented for multiple references.\n");
    assert(comm);
    // current partition over 'u' for L/Piu
    assert(Luv.size() == std::get<1>(Piu.sizes()));
    for (int i = 0; i < rotcPua.size(); i++)
    {
      // rot Ps are not yet distributed
      assert(rotcPua[i].size() == std::get<1>(rotPiu.sizes()));
      if (walker_type == CLOSED)
        assert(std::get<1>(rotcPua[i].sizes()) == nup);
      else if (walker_type == COLLINEAR)
        assert(std::get<1>(rotcPua[i].sizes()) == nup + ndown);
      else if (walker_type == NONCOLLINEAR)
        assert(std::get<1>(rotcPua[i].sizes()) == nup + ndown);
    }
    for (int i = 0; i < cPua.size(); i++)
    {
      assert(cPua[i].size() == Luv.size());
      if (walker_type == CLOSED)
        assert(std::get<1>(cPua[i].sizes()) == nup);
      else if (walker_type == COLLINEAR)
        assert(std::get<1>(cPua[i].sizes()) == nup + ndown);
      else if (walker_type == NONCOLLINEAR)
        assert(std::get<1>(cPua[i].sizes()) == nup + ndown);
    }
    if (walker_type == NONCOLLINEAR)
    {
      assert(Piu.size() == 2 * NMO);
      assert(rotPiu.size() == 2 * NMO);
    }
    else
    {
      assert(Piu.size() == NMO);
      assert(rotPiu.size() == NMO);
    }
  }

  ~THCOps() {}

  THCOps(THCOps const& other) = delete;
  THCOps& operator=(THCOps const& other) = delete;

  THCOps(THCOps&& other) = default;
  THCOps& operator=(THCOps&& other) = default;

  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(TaskGroup_& TG,
                                                                 boost::multi::array<ComplexType, 1> const& vMF)
  {
    using std::copy_n;
    using std::fill_n;
    int NMO = hij.size();
    // in non-collinear case with SO, keep SO matrix here and add it
    // for now, stay collinear

    ShmArray<ComplexType, 1> vMF_(vMF, shm_buffer_manager.get_generator().template get_allocator<ComplexType>());
    ShmArray<ComplexType, 1> P1D(iextensions<1u>{NMO * NMO}, ComplexType(0),
                                 shm_buffer_manager.get_generator().template get_allocator<ComplexType>());

    vHS(vMF_, P1D);
    if (TG.TG_Cores().size() > 1 && TG.TG_local().root())
      TG.TG_Cores().all_reduce_in_place_n(to_address(P1D.origin()), P1D.num_elements(), std::plus<>());
    TG.TG().barrier();

    boost::multi::array<ComplexType, 2> H1({NMO, NMO});
    copy_n(P1D.origin(), NMO * NMO, H1.origin());

    // add hij + vn0 and symmetrize
    using ma::conj;
    for (int i = 0; i < NMO; i++)
    {
      H1[i][i] += hij[i][i] + vn0[i][i];
      for (int j = i + 1; j < NMO; j++)
      {
        H1[i][j] += hij[i][j] + vn0[i][j];
        H1[j][i] += hij[j][i] + vn0[j][i];
#if defined(MIXED_PRECISION)
        if (std::abs(H1[i][j] - ma::conj(H1[j][i])) > 1e-5)
        {
#else
        if (std::abs(H1[i][j] - ma::conj(H1[j][i])) > 1e-6)
        {
#endif
          app_error() << " WARNING in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
          app_error() << i << " " << j << " " << H1[i][j] << " " << H1[j][i] << " "
                      << H1[i][j] - (hij[i][j] + vn0[i][j]) << " " << H1[j][i] - (hij[j][i] + vn0[j][i]) << " "
                      << hij[i][j] << " " << hij[j][i] << " " << vn0[i][j] << " " << vn0[j][i] << std::endl;
          //APP_ABORT("Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n");
        }
        H1[i][j] = 0.5 * (H1[i][j] + ma::conj(H1[j][i]));
        H1[j][i] = ma::conj(H1[i][j]);
      }
    }

    return H1;
  }

  template<class Mat, class MatB>
  void energy(Mat&& E, MatB const& G, int k, bool addH1 = true, bool addEJ = true, bool addEXX = true)
  {
    boost::multi::array<SPComplexType, 2>* Kr(nullptr);
    boost::multi::array<SPComplexType, 2>* Kl(nullptr);
    energy(E, G, k, Kl, Kr, addH1, addEJ, addEXX);
  }

  // Kl and Kr must be in shared memory for this to work correctly
  template<class Mat, class MatB, class MatC, class MatD>
  void energy(Mat&& E,
              MatB const& G,
              int k,
              MatC* Kl,
              MatD* Kr,
              bool addH1  = true,
              bool addEJ  = true,
              bool addEXX = true)
  {
    using std::copy_n;
    using std::fill_n;
    using GType = typename std::decay_t<typename MatB::element>;
    if (k > 0)
      APP_ABORT(" Error: THC not yet implemented for multiple references.\n");
    // G[nel][nmo]
    assert(std::get<0>(E.sizes()) == std::get<0>(G.sizes()));
    assert(std::get<1>(E.sizes()) == 3);
    int nwalk = G.size();
    int getKr = Kr != nullptr;
    int getKl = Kl != nullptr;

    // addH1
    fill_n(E.origin(), E.num_elements(), ComplexType(0.0));
    if (addH1)
    {
      ma::product(ComplexType(1.0), G, haj[k], ComplexType(0.0), E(E.extension(0), 0));
      for (int i = 0; i < nwalk; i++)
        E[i][0] += E0;
    }
    if (not(addEJ || addEXX))
      return;

    int nmo_  = rotPiu.size();
    int nu    = rotMuv.size();
    int nu0   = rotnmu0;
    int nv    = std::get<1>(rotMuv.sizes());
    int nel_  = std::get<1>(rotcPua[0].sizes());
    int nspin = (walker_type == COLLINEAR) ? 2 : 1;
    assert(std::get<1>(G.sizes()) == nel_ * nmo_);
    if (addEJ and getKl)
      assert(std::get<0>(Kl->sizes()) == nwalk && std::get<1>(Kl->sizes()) == nu);
    if (addEJ and getKr)
      assert(std::get<0>(Kr->sizes()) == nwalk && std::get<1>(Kr->sizes()) == nu);
    using ma::T;
    int u0, uN;
    std::tie(u0, uN) = FairDivideBoundary(comm->rank(), nu, comm->size());
    int v0, vN;
    std::tie(v0, vN) = FairDivideBoundary(comm->rank(), nv, comm->size());
    int k0, kN;
    std::tie(k0, kN) = FairDivideBoundary(comm->rank(), nmo_, comm->size());

    // calculate how many walkers can be done concurrently
    long mem_needs(0);
    if (not std::is_same<GType, SPComplexType>::value)
      mem_needs += G.num_elements();
    long Bytes = default_buffer_size_in_MB * 1024L * 1024L;
    Bytes -= mem_needs * long(sizeof(SPComplexType));
    Bytes /= long((nu * nv + nv + nv * nup) * sizeof(SPComplexType));
    int nwmax = std::min(nwalk, std::max(1, int(Bytes)));
    ShmArray<SPComplexType, 1> Gbuff(iextensions<1u>{mem_needs},
                                     shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());

    const_sp_pointer Gptr(nullptr);
    // setup origin of Gsp and copy_n_cast if necessary
    if (std::is_same<GType, SPComplexType>::value)
    {
      Gptr = pointer_cast<SPComplexType const>(make_device_ptr(G.origin()));
    }
    else
    {
      long i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(G.num_elements()), long(comm->size()));
      copy_n_cast(make_device_ptr(G.origin()) + i0, iN - i0, make_device_ptr(Gbuff.origin()) + i0);
      Gptr = make_device_ptr(Gbuff.origin());
    }
    Array_cref<SPComplexType, 2> Gsp(Gptr, G.extensions());

    // Guv[nspin][nu][nv]
    ShmArray<SPComplexType, 3> Guv({nwmax, nu, nv},
                                   shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    // Guu[u]: summed over spin
    ShmArray<SPComplexType, 2> Guu({nwmax, nv},
                                   shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());

    ShmArray<SPComplexType, 3> Tav({nwmax, nup, nv},
                                   shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());

    SPRealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
    int iw(0);
    while (iw < nwalk)
    {
      int nw = std::min(nwmax, nwalk - iw);
      fill_n(Guu.origin(), Guu.num_elements(), SPComplexType(0.0));
      for (int ispin = 0; ispin < nspin; ++ispin)
      {
        Guv_Guu(ispin, Gsp.sliced(iw, iw + nw), Guv, Guu, Tav, k);

        // Gwuv = Gwuv * rotMuv
        using ma::inplace_product;
        inplace_product(nw, nu, (vN - v0), make_device_ptr(rotMuv.origin()) + v0, nv,
                        make_device_ptr(Guv.origin()) + v0, nv);
        comm->barrier();

        // R[w,u][b] = sum_v Guv[w,u][v] * cPua[v][b]
        long i0, iN;
        std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(nw * nu), long(comm->size()));
        Array_ref<SPComplexType, 2> Rwub(make_device_ptr(Tav.origin()), {nw * nu, nelec[ispin]});
        Array_ref<SPComplexType, 2> Guv2D(make_device_ptr(Guv.origin()), {nw * nu, nv});
        ma::product(Guv2D.sliced(i0, iN), rotcPua[k]({0, nv}, {ispin * nup, nup + ispin * ndown}), Rwub.sliced(i0, iN));
        comm->barrier();

        //T[w][b][k] = sum_u R[w][u][b] * Piu[k][u]
        // need batching in this case
        Array_ref<SPComplexType, 3> Rwub3D(Rwub.origin(), {nw, nu, nelec[ispin]});
        Array_ref<SPComplexType, 3> Twbk(make_device_ptr(Guv.origin()), {nw, nelec[ispin], nmo_});
        Array_ref<SPComplexType, 2> Twbk2D(Twbk.origin(), {nw, nelec[ispin] * nmo_});
        std::vector<decltype(&(Rwub3D[0]))> vRwub;
        std::vector<decltype(&(rotPiu({0, 1}, {0, 1})))> vPku;
        std::vector<decltype(&(Twbk[0]({0, 1}, {0, 1})))> vTwbk;
        vRwub.reserve(nw);
        vPku.reserve(nw);
        vTwbk.reserve(nw);
        for (int w = 0; w < nw; ++w)
        {
          vRwub.emplace_back(&(Rwub3D[w]));
          vPku.emplace_back(&(rotPiu({k0, kN}, {nu0, nu0 + nu})));
          vTwbk.emplace_back(&(Twbk[w]({0, nelec[ispin]}, {k0, kN})));
        }
#if defined(QMC_COMPLEX)
        ma::BatchedProduct('T', 'T', vRwub, vPku, vTwbk);
#else
        // need to keep vPku on the left hand side in real build
        if (Guv.num_elements() >= 2 * Twbk.num_elements())
        {
          Array_ref<SPComplexType, 3> Twkb(Twbk.origin() + Twbk.num_elements(), {nw, nmo_, nelec[ispin]});
          std::vector<decltype(&(Twkb[0].sliced(0, 1)))> vTwkb;
          vTwkb.reserve(nw);
          for (int w = 0; w < nw; ++w)
            vTwkb.emplace_back(&(Twkb[w].sliced(k0, kN)));
          ma::BatchedProduct('N', 'N', vPku, vRwub, vTwkb);
          for (int w = 0; w < nw; ++w)
            ma::transpose(Twkb[w].sliced(k0, kN), Twbk[w]({0, nelec[ispin]}, {k0, kN}));
        }
        else
        {
          Array<SPComplexType, 3> Twkb({nw, (kN - k0), nelec[ispin]},
                                       device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
          std::vector<decltype(&(Twkb[0]))> vTwkb;
          vTwkb.reserve(nw);
          for (int w = 0; w < nw; ++w)
            vTwkb.emplace_back(&(Twkb[w]));
          ma::BatchedProduct('N', 'N', vPku, vRwub, vTwkb);
          for (int w = 0; w < nw; ++w)
            ma::transpose(Twkb[w], Twbk[w]({0, nelec[ispin]}, {k0, kN}));
        }
#endif
        comm->barrier();

        // E[w] = sum_bk T[w][bk] G[w][bk]
        // move to batched_dot!!!
        // or to batched_adotpby
        std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(nelec[ispin] * nmo_), long(comm->size()));
        using ma::adotpby;
        adotpby(SPComplexType(-0.5 * scl), Gsp({0, nw}, {ispin * nelec[0] * nmo_ + i0, ispin * nelec[0] * nmo_ + iN}),
                Twbk2D({0, nw}, {i0, iN}), ComplexType(0.0), E({iw, iw + nw}, 1));
        comm->barrier();
      }
      comm->barrier();
      if (addEJ)
      {
        Array<SPComplexType, 2> Twu({nw, (uN - u0)},
                                    device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
#if defined(QMC_COMPLEX)
        ma::product(Guu.sliced(0, nw), ma::T(rotMuv.sliced(u0, uN)), Twu);
#else
        // need to keep rotMuv on the left hand side in real build
        Array<SPComplexType, 2> Tuw({(uN - u0), nw},
                                    device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
        ShmArray<SPComplexType, 2> Gvw({nv, nw},
                                       shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());
        ma::transpose(Guu({0, nw}, {v0, vN}), Gvw.sliced(v0, vN));
        comm->barrier();
        ma::product(rotMuv.sliced(u0, uN), Gvw, Tuw);
        ma::transpose(Tuw, Twu);
#endif
        comm->barrier();
        using ma::adotpby;
        adotpby(SPComplexType(0.5 * scl * scl), Guu({0, nw}, {nu0 + u0, nu0 + uN}), Twu, ComplexType(0.0),
                E({iw, iw + nw}, 2));

        if (getKl)
          ma::add(SPComplexType(1.0), Guu({0, nw}, {nu0 + u0, nu0 + uN}), SPComplexType(0.0), Twu,
                  (*Kl)({iw, iw + nw}, {u0, uN}));
        if (getKr)
          ma::add(SPComplexType(1.0), Twu, SPComplexType(0.0), Twu, (*Kr)({iw, iw + nw}, {u0, uN}));
      }
      comm->barrier();
      iw += nw;
    }
    comm->barrier();
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
  {
    APP_ABORT(" Error: fast_energy not yet working");
    if (haj.size() != 1)
      APP_ABORT(" Error: Single reference implementation currently in THCOps::fast_energy.\n");
    if (walker_type != CLOSED)
      APP_ABORT(" Error: THCOps::fast_energy requires walker_type==CLOSED.\n");
    /*
       * E[nspins][maxn_unique_confg][nwalk][3]
       * Ov[nspins][maxn_unique_confg][nwalk]
       * GrefA[nwalk][nup][NMO]
       * GrefB[nwalk][ndown][NMO]
       * QQ0A[nwalk][nup][NAEA]
       * QQ0B[nwalk][nup][NAEA]
       */
    /*
      static_assert(std::decay<MatE>::type::dimensionality==4, "Wrong dimensionality");
      static_assert(std::decay<MatO>::type::dimensionality==3, "Wrong dimensionality");
      static_assert(std::decay<MatG>::type::dimensionality==3, "Wrong dimensionality");
      static_assert(std::decay<MatQ>::type::dimensionality==3, "Wrong dimensionality");
      //static_assert(std::decay<MatB>::type::dimensionality==3, "Wrong dimensionality");
      int nspin = E.size(0);
      int nrefs = haj.size();
      int nwalk = GrefA.size(0);
      int naoa_ = QQ0A.size(1);
      int naob_ = QQ0B.size(1);
      int nmo_ = rotPiu.size(0);
      int nu = rotMuv.size(0);
      int nu0 = rotnmu0; 
      int nv = rotMuv.size(1);
      int nel_ = rotcPua[0].size(1);
      // checking
      assert(E.size(2) == nwalk);
      assert(E.size(3) == 3);
      assert(Ov.size(0) == nspin);
      assert(Ov.size(1) == E.size(1));
      assert(Ov.size(2) == nwalk);
      assert(GrefA.size(1) == naoa_);
      assert(GrefA.size(2) == nmo_);
      assert(GrefB.size(0) == nwalk);
      assert(GrefB.size(1) == naob_);
      assert(GrefB.size(2) == nmo_);
      // limited to single reference now
      assert(rotcPua.size() == nrefs);
      assert(nel_ == naoa_);
      assert(nel_ == naob_);

      using ma::T;
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int k0,kN;
      std::tie(k0,kN) = FairDivideBoundary(comm->rank(),nel_,comm->size());
      // right now the algorithm uses 2 copies of matrices of size nuxnv in COLLINEAR case,
      // consider moving loop over spin to avoid storing the second copy which is not used
      // simultaneously
      size_t memory_needs = nu*nv + nv + nu  + nel_*(nv+2*nu+2*nel_);
      set_shmbuffer(memory_needs);
      size_t cnt=0;
      // if Alpha/Beta have different references, allocate the largest and
      // have distinct references for each
      // Guv[nu][nv]
      boost::multi::array_ref<ComplexType,2> Guv(to_address(SM_TMats.origin()),{nu,nv});
      cnt+=Guv.num_elements();
      // Gvv[v]: summed over spin
      boost::multi::array_ref<ComplexType,1> Gvv(to_address(SM_TMats.origin())+cnt,iextensions<1u>{nv});
      cnt+=Gvv.num_elements();
      // S[nel_][nv]
      boost::multi::array_ref<ComplexType,2> Scu(to_address(SM_TMats.origin())+cnt,{nel_,nv});
      cnt+=Scu.num_elements();
      // Qub[nu][nel_]:
      boost::multi::array_ref<ComplexType,2> Qub(to_address(SM_TMats.origin())+cnt,{nu,nel_});
      cnt+=Qub.num_elements();
      boost::multi::array_ref<ComplexType,1> Tuu(to_address(SM_TMats.origin())+cnt,iextensions<1u>{nu});
      cnt+=Tuu.num_elements();
      boost::multi::array_ref<ComplexType,2> Jcb(to_address(SM_TMats.origin())+cnt,{nel_,nel_});
      cnt+=Jcb.num_elements();
      boost::multi::array_ref<ComplexType,2> Xcb(to_address(SM_TMats.origin())+cnt,{nel_,nel_});
      cnt+=Xcb.num_elements();
      boost::multi::array_ref<ComplexType,2> Tub(to_address(SM_TMats.origin())+cnt,{nu,nel_});
      cnt+=Tub.num_elements();
      assert(cnt <= memory_needs);
      boost::multi::static_array<ComplexType,3,dev_buffer_type> eloc({2,nwalk,3}
                        device_buffer_manager.get_generator().template get_allocator<ComplexType>());
      std::fill_n(eloc.origin(),eloc.num_elements(),ComplexType(0.0));

      RealType scl = (walker_type==CLOSED?2.0:1.0);
      if(comm->root()) {
        std::fill_n(to_address(E.origin()),E.num_elements(),ComplexType(0.0));
        std::fill_n(to_address(Ov[0][1].origin()),nwalk*(Ov.size(1)-1),ComplexType(0.0));
        std::fill_n(to_address(Ov[1][1].origin()),nwalk*(Ov.size(1)-1),ComplexType(0.0));
        auto Ea = E[0][0];
        auto Eb = E[1][0];
        boost::multi::array_cref<ComplexType,2> G2DA(to_address(GrefA.origin()),
                                          {nwalk,GrefA[0].num_elements()});
        ma::product(ComplexType(1.0),G2DA,haj[0],ComplexType(0.0),Ea(Ea.extension(0),0));
        boost::multi::array_cref<ComplexType,2> G2DB(to_address(GrefA.origin()),
                                          {nwalk,GrefA[0].num_elements()});
        ma::product(ComplexType(1.0),G2DB,haj[0],ComplexType(0.0),Eb(Eb.extension(0),0));
        for(int i=0; i<nwalk; i++) {
            Ea[i][0] += E0;
            Eb[i][0] += E0;
        }
      }

      for(int wi=0; wi<nwalk; wi++) {

        { // Alpha
          auto Gw = GrefA[wi];
          boost::multi::array_cref<ComplexType,1> G1D(to_address(Gw.origin()),
                                                        iextensions<1u>{Gw.num_elements()});
          Guv_Guu2(Gw,Guv,Gvv,Scu,0);
          if(u0!=uN)
            ma::product(rotMuv.sliced(u0,uN),Gvv,
                      Tuu.sliced(u0,uN));
          auto Mptr = rotMuv[u0].origin();
          auto Gptr = to_address(Guv[u0].origin());
          for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
            (*Gptr) *= (*Mptr);
          if(u0!=uN)
            ma::product(Guv.sliced(u0,uN),rotcPua[0],
                      Qub.sliced(u0,uN));
          comm->barrier();
          if(k0!=kN)
            ma::product(Scu.sliced(k0,kN),Qub,
                      Xcb.sliced(k0,kN));
          // Tub = rotcPua.*Tu
          auto rPptr = rotcPua[0][nu0+u0].origin();
          auto Tuuptr = Tuu.origin()+u0;
          auto Tubptr = Tub[u0].origin();
          for(size_t u_=u0; u_<uN; ++u_, ++Tuuptr)
            for(size_t k=0; k<nel_; ++k, ++rPptr, ++Tubptr)
              (*Tubptr) = (*Tuuptr)*(*rPptr);
          comm->barrier();
          // Jcb = Scu*Tub
          if(k0!=kN)
            ma::product(Scu.sliced(k0,kN),Tub,
                      Jcb.sliced(k0,kN));
          for(int c=k0; c<kN; ++c)
            eloc[0][wi][1] += -0.5*scl*Xcb[c][c];
          for(int c=k0; c<kN; ++c)
            eloc[0][wi][2] += 0.5*scl*scl*Jcb[c][c];
          calculate_ph_energies(0,comm->rank(),comm->size(),
                                E[0],Ov[0],QQ0A,Qwork,
                                rotMuv,
                                abij,det_couplings);
        }

        { // Beta: Unnecessary in CLOSED walker type (on Walker)
          auto Gw = GrefB[wi];
          boost::multi::array_cref<ComplexType,1> G1D(to_address(Gw.origin()),
                                                        iextensions<1u>{Gw.num_elements()});
          Guv_Guu2(Gw,Guv,Gvv,Scu,0);
          if(u0!=uN)
            ma::product(rotMuv.sliced(u0,uN),Gvv,
                      Tuu.sliced(u0,uN));
          auto Mptr = rotMuv[u0].origin();
          auto Gptr = to_address(Guv[u0].origin());
          for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
            (*Gptr) *= (*Mptr);
          if(u0!=uN)
            ma::product(Guv.sliced(u0,uN),rotcPua[0],
                      Qub.sliced(u0,uN));
          comm->barrier();
          if(k0!=kN)
            ma::product(Scu.sliced(k0,kN),Qub,
                      Xcb.sliced(k0,kN));
          // Tub = rotcPua.*Tu
          auto rPptr = rotcPua[0][nu0+u0].origin();
          auto Tuuptr = Tuu.origin()+u0;
          auto Tubptr = Tub[u0].origin();
          for(size_t u_=u0; u_<uN; ++u_, ++Tuuptr)
            for(size_t k=0; k<nel_; ++k, ++rPptr, ++Tubptr)
              (*Tubptr) = (*Tuuptr)*(*rPptr);
          comm->barrier();
          // Jcb = Scu*Tub
          if(k0!=kN)
            ma::product(Scu.sliced(k0,kN),Tub,
                      Jcb.sliced(k0,kN));
          for(int c=k0; c<kN; ++c)
            eloc[1][wi][1] += -0.5*scl*Xcb[c][c];
          for(int c=k0; c<kN; ++c)
            eloc[1][wi][2] += 0.5*scl*scl*Jcb[c][c];
        }

      }
      comm->reduce_in_place_n(eloc.origin(),eloc.num_elements(),std::plus<>(),0);
      if(comm->root()) {
        // add Eref contributions to all configurations
        for(int nd=0; nd<E.size(1); ++nd) {
          auto Ea = E[0][nd];
          auto Eb = E[1][nd];
          for(int wi=0; wi<nwalk; wi++) {
            Ea[wi][1] += eloc[0][wi][1];
            Ea[wi][2] += eloc[0][wi][2];
            Eb[wi][1] += eloc[1][wi][1];
            Eb[wi][2] += eloc[1][wi][2];
          }
        }
      }
      comm->barrier();
*/
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vHS(MatA const& X, MatB&& v, double a = 1., double c = 0.)
  {
    using XType = typename std::decay_t<typename MatA::element>;
    using vType = typename std::decay<MatB>::type::element;
    boost::multi::array_ref<vType, 2, decltype(v.origin())> v_(v.origin(), {1, v.size()});
    boost::multi::array_ref<XType const, 2, decltype(X.origin())> X_(X.origin(), {X.size(), 1});
    vHS(X_, v_, a, c);
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    using ma::T;
    using XType = typename std::decay_t<typename MatA::element>;
    using vType = typename std::decay<MatB>::type::element;
    int nwalk   = std::get<1>(X.sizes());
#if defined(QMC_COMPLEX)
    int nchol = 2 * std::get<1>(Luv.sizes());
#else
    int nchol = std::get<1>(Luv.sizes());
#endif
    int nmo_ = std::get<0>(Piu.sizes());
    int nu   = std::get<1>(Piu.sizes());
    assert(std::get<0>(Luv.sizes()) == nu);
    assert(std::get<0>(X.sizes()) == nchol);
    assert(std::get<0>(v.sizes()) == nwalk);
    assert(std::get<1>(v.sizes()) == nmo_ * nmo_);

    size_t memory_needs = nu * nwalk;
    if (not std::is_same<XType, SPComplexType>::value)
      memory_needs += X.num_elements();
    if (not std::is_same<vType, SPComplexType>::value)
      memory_needs += v.num_elements();

    // calculate how many walkers can be done concurrently
    long Bytes = default_buffer_size_in_MB * 1024L * 1024L;
    // memory_needs = X, v, Tuw
    Bytes -= size_t(memory_needs * sizeof(SPComplexType)); // subtract other needs
    Bytes /= size_t(nmo_ * nu * sizeof(SPComplexType));
    int nwmax = std::min(nwalk, std::max(1, int(Bytes)));
    memory_needs += nwmax * nmo_ * nu;
    ShmArray<SPComplexType, 1> SM_TMats(iextensions<1u>(memory_needs),
                                        shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());

    size_t cnt(0);
    const_sp_pointer Xptr(nullptr);
    sp_pointer vptr(nullptr);
    // setup origin of vsp and copy_n_cast if necessary
    if (std::is_same<vType, SPComplexType>::value)
    {
      vptr = pointer_cast<SPComplexType>(make_device_ptr(v.origin()));
    }
    else
    {
      long i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(v.num_elements()), long(comm->size()));
      vptr             = make_device_ptr(SM_TMats.origin());
      cnt += size_t(v.num_elements());
      if (std::abs(c) > 1e-12)
        copy_n_cast(make_device_ptr(v.origin()) + i0, iN - i0, vptr + i0);
    }
    // setup origin of Xsp and copy_n_cast if necessary
    if (std::is_same<XType, SPComplexType>::value)
    {
      Xptr = pointer_cast<SPComplexType const>(make_device_ptr(X.origin()));
    }
    else
    {
      long i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(X.num_elements()), long(comm->size()));
      copy_n_cast(make_device_ptr(X.origin()) + i0, iN - i0, make_device_ptr(SM_TMats.origin()) + cnt + i0);
      Xptr = make_device_ptr(SM_TMats.origin()) + cnt;
      cnt += size_t(X.num_elements());
    }
    // setup array references
    Array_cref<SPComplexType, 2> Xsp(Xptr, X.extensions());
    Array_ref<SPComplexType, 2> vsp(vptr, v.extensions());

    int u0, uN;
    std::tie(u0, uN) = FairDivideBoundary(comm->rank(), nu, comm->size());
    Array_ref<SPComplexType, 2> Tuw(make_device_ptr(SM_TMats.origin()) + cnt, {nu, nwalk});
    // O[nwalk * nmu * nmu]
#if defined(QMC_COMPLEX)
    // reinterpret as RealType matrices with 2x the columns
    Array_ref<SPRealType, 2> Luv_R(pointer_cast<SPRealType>(make_device_ptr(Luv.origin())),
                                   {std::get<0>(Luv.sizes()), 2 * std::get<1>(Luv.sizes())});
    Array_cref<SPRealType, 2> X_R(pointer_cast<SPRealType const>(Xsp.origin()), {std::get<0>(Xsp.sizes()), 2 * std::get<1>(Xsp.sizes())});
    Array_ref<SPRealType, 2> Tuw_R(pointer_cast<SPRealType>(Tuw.origin()), {nu, 2 * nwalk});
    ma::product(Luv_R.sliced(u0, uN), X_R, Tuw_R.sliced(u0, uN));
#else
    ma::product(Luv.sliced(u0, uN), Xsp, Tuw.sliced(u0, uN));
#endif
    comm->barrier();
    int k0, kN;
    std::tie(k0, kN) = FairDivideBoundary(comm->rank(), nmo_, comm->size());
    int iw(0);
    while (iw < nwalk)
    {
      int nw = std::min(nwmax, nwalk - iw);
#if defined(QMC_COMPLEX)
      // Qwiu[w][i][u] = T[u][w] * conj(Piu[i][u])
      // v[w][i][k] = sum_u Qwiu[w][i][u] * Piu[k][u]
      Array_ref<SPComplexType, 2> Qwiu(Tuw.origin() + Tuw.num_elements(), {nw * nmo_, nu});
      using ma::element_wise_Aij_Bjk_Ckij;
      element_wise_Aij_Bjk_Ckij('C', (kN - k0), nu, nw, make_device_ptr(Piu[k0].origin()), nu,
                                make_device_ptr(Tuw.origin()) + iw, nwalk, make_device_ptr(Qwiu.origin()) + k0 * nu,
                                nmo_, nu);
      comm->barrier();
      // v[w][i][j] = sum_u Qwiu[w][i][u] * Piu[j][u]
      Array_ref<SPComplexType, 2> v_(vsp[iw].origin(), {nw * nmo_, nmo_});
      int wk0, wkN;
      std::tie(wk0, wkN) = FairDivideBoundary(comm->rank(), nw * nmo_, comm->size());
      ma::product(SPComplexType(a), Qwiu.sliced(wk0, wkN), T(Piu), SPComplexType(c), v_.sliced(wk0, wkN));
#else
      // Qwui[w][u][i] = Piu[i][u] * T[u][w]
      // v[w][i][j] = sum_u Piu[i][u] Qwui[w][u][j] // using batched blas
      Array_ref<SPComplexType, 3> Qwui(Tuw.origin() + Tuw.num_elements(), {nw, nu, nmo_});
      using ma::element_wise_Aij_Bjk_Ckji;
      element_wise_Aij_Bjk_Ckji((kN - k0), nu, nw, make_device_ptr(Piu[k0].origin()), nu,
                                make_device_ptr(Tuw.origin()) + iw, nwalk, make_device_ptr(Qwui.origin()) + k0, nmo_,
                                nu * nmo_);
      comm->barrier();
      Array_ref<SPComplexType, 3> v_(vsp[iw].origin(), {nw, nmo_, nmo_});
      std::vector<decltype(&(Piu.sliced(0, 1)))> vPiu;
      std::vector<decltype(&(Qwui[0]))> vQui;
      std::vector<decltype(&(v_[0].sliced(0, 1)))> vVij;
      vPiu.reserve(nw);
      vQui.reserve(nw);
      vVij.reserve(nw);
      for (int w = 0; w < nw; ++w)
      {
        vPiu.emplace_back(&(Piu.sliced(k0, kN)));
        vQui.emplace_back(&(Qwui[w]));
        vVij.emplace_back(&(v_[w].sliced(k0, kN)));
      }
      ma::BatchedProduct('N', 'N', SPRealType(a), vPiu, vQui, SPRealType(c), vVij);
#endif
      iw += nw;
      comm->barrier();
    }
    if (not std::is_same<vType, SPComplexType>::value)
    {
      long i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(v.num_elements()), long(comm->size()));
      copy_n_cast(vsp.origin() + i0, iN - i0, make_device_ptr(v.origin()) + i0);
    }
    comm->barrier();
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vbias(MatA const& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    using GType = typename std::decay_t<typename MatA::element>;
    using vType = typename std::decay<MatB>::type::element;
    boost::multi::array_ref<vType, 2, decltype(v.origin())> v_(v.origin(), {std::get<0>(v.sizes()), 1});
    boost::multi::array_ref<GType const, 2, decltype(G.origin())> G_(G.origin(), {1, std::get<0>(G.sizes())});
    vbias(G_, v_, a, c, k);
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vbias(MatA const& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    using GType = typename std::decay_t<typename MatA::element>;
    using vType = typename std::decay<MatB>::type::element;
    if (k > 0)
      APP_ABORT(" Error: THC not yet implemented for multiple references.\n");
    int nwalk = std::get<0>(G.sizes());
    int nmo_  = std::get<0>(Piu.sizes());
    int nu    = std::get<1>(Piu.sizes());
    int nel_  = std::get<1>(cPua[0].sizes());
#if defined(QMC_COMPLEX)
    int nchol = 2 * std::get<1>(Luv.sizes());
#else
    int nchol = std::get<1>(Luv.sizes());
#endif
    assert(std::get<1>(v.sizes()) == nwalk);
    assert(std::get<0>(v.sizes()) == nchol);
    using ma::T;
    int c0, cN;
    std::tie(c0, cN) = FairDivideBoundary(comm->rank(), nchol, comm->size());

    size_t memory_needs = nwalk * nu;
    if (not std::is_same<GType, SPComplexType>::value)
      memory_needs += G.num_elements();
    if (not std::is_same<vType, SPComplexType>::value)
      memory_needs += v.num_elements();
    ShmArray<SPComplexType, 1> SM_TMats(iextensions<1u>(memory_needs),
                                        shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    size_t cnt(0);
    const_sp_pointer Gptr(nullptr);
    sp_pointer vptr(nullptr);
    // setup origin of Gsp and copy_n_cast if necessary
    if (std::is_same<GType, SPComplexType>::value)
    {
      Gptr = pointer_cast<SPComplexType const>(make_device_ptr(G.origin()));
    }
    else
    {
      long i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(G.num_elements()), long(comm->size()));
      copy_n_cast(make_device_ptr(G.origin()) + i0, iN - i0, make_device_ptr(SM_TMats.origin()) + i0);
      cnt += size_t(G.num_elements());
      Gptr = make_device_ptr(SM_TMats.origin());
    }
    // setup origin of vsp and copy_n_cast if necessary
    if (std::is_same<vType, SPComplexType>::value)
    {
      vptr = pointer_cast<SPComplexType>(make_device_ptr(v.origin()));
    }
    else
    {
      long i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(long(comm->rank()), long(v.num_elements()), long(comm->size()));
      vptr             = make_device_ptr(SM_TMats.origin()) + cnt;
      cnt += size_t(v.num_elements());
      if (std::abs(c) > 1e-12)
        copy_n_cast(make_device_ptr(v.origin()) + i0, iN - i0, vptr + i0);
    }
    // setup array references
    Array_cref<SPComplexType, 2> Gsp(Gptr, G.extensions());
    Array_ref<SPComplexType, 2> vsp(vptr, v.extensions());

    if (haj.size() == 1)
    {
      Array_ref<SPComplexType, 2> Guu(make_device_ptr(SM_TMats.origin()) + cnt, {nu, nwalk});
      Guu_from_compact(Gsp, Guu);
#if defined(QMC_COMPLEX)
      // reinterpret as RealType matrices with 2x the columns
      Array_ref<SPRealType, 2> Luv_R(pointer_cast<SPRealType>(make_device_ptr(Luv.origin())),
                                     {std::get<0>(Luv.sizes()), 2 * std::get<1>(Luv.sizes())});
      Array_ref<SPRealType, 2> Guu_R(pointer_cast<SPRealType>(Guu.origin()), {nu, 2 * nwalk});
      Array_ref<SPRealType, 2> vsp_R(pointer_cast<SPRealType>(vsp.origin()), {std::get<0>(vsp.sizes()), 2 * std::get<1>(vsp.sizes())});
      ma::product(SPRealType(a), T(Luv_R(Luv_R.extension(0), {c0, cN})), Guu_R, SPRealType(c), vsp_R.sliced(c0, cN));
#else
      ma::product(SPRealType(a), T(Luv(Luv.extension(0), {c0, cN})), Guu, SPRealType(c), vsp.sliced(c0, cN));
#endif
    }
    else
    {
      Array_ref<SPComplexType, 2> Guu(make_device_ptr(SM_TMats.origin()) + cnt, {nu, nwalk});
      Guu_from_full(Gsp, Guu);
#if defined(QMC_COMPLEX)
      // reinterpret as RealType matrices with 2x the columns
      Array_ref<SPRealType, 2> Luv_R(pointer_cast<SPRealType>(make_device_ptr(Luv.origin())),
                                     {std::get<0>(Luv.sizes()), 2 * std::get<1>(Luv.sizes())});
      Array_ref<SPRealType, 2> Guu_R(pointer_cast<SPRealType>(Guu.origin()), {nu, 2 * nwalk});
      Array_ref<SPRealType, 2> vsp_R(pointer_cast<SPRealType>(vsp.origin()), {std::get<0>(vsp.sizes()), 2 * std::get<1>(vsp.sizes())});
      ma::product(SPRealType(a), T(Luv_R(Luv_R.extension(0), {c0, cN})), Guu_R, SPRealType(c), vsp_R.sliced(c0, cN));
#else
      ma::product(SPRealType(a), T(Luv(Luv.extension(0), {c0, cN})), Guu, SPRealType(c), vsp.sliced(c0, cN));
#endif
    }
    if (not std::is_same<vType, SPComplexType>::value)
    {
      copy_n_cast(make_device_ptr(vsp[c0].origin()), std::get<1>(vsp.sizes()) * (cN - c0), make_device_ptr(v[c0].origin()));
    }
    comm->barrier();
  }

  template<class Mat, class MatB>
  void generalizedFockMatrix(Mat&& G, MatB&& Fp, MatB&& Fm)
  {
    APP_ABORT(" Error: generalizedFockMatrix not implemented for this hamiltonian.\n");
  }

  bool distribution_over_cholesky_vectors() const { return false; }
  int number_of_ke_vectors() const { return std::get<0>(rotMuv.sizes()); }
#if defined(QMC_COMPLEX)
  int local_number_of_cholesky_vectors() const { return 2 * std::get<1>(Luv.sizes()); }
  int global_number_of_cholesky_vectors() const { return 2 * std::get<1>(Luv.sizes()); }
#else
  int local_number_of_cholesky_vectors() const { return std::get<1>(Luv.sizes()); }
  int global_number_of_cholesky_vectors() const { return std::get<1>(Luv.sizes()); }
#endif
  int global_origin_cholesky_vector() const { return 0; }

  // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
  bool transposed_G_for_vbias() const { return true; }
  bool transposed_G_for_E() const { return true; }
  // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
  bool transposed_vHS() const { return true; }

  bool fast_ph_energy() const { return false; }

  boost::multi::array<ComplexType, 2> getHSPotentials() { return boost::multi::array<ComplexType, 2>{}; }

protected:
  // Guu[nu][nwalk]
  template<class MatA, class MatB>
  void Guu_from_compact(MatA const& G, MatB&& Guu)
  {
    int nmo_ = int(std::get<0>(Piu.sizes()));
    int nu   = int(std::get<1>(Piu.sizes()));
    int nel_ = std::get<1>(cPua[0].sizes());
    int u0, uN;
    std::tie(u0, uN) = FairDivideBoundary(comm->rank(), nu, comm->size());
    int nw           = std::get<0>(G.sizes());

    assert(std::get<0>(G.sizes()) == std::get<1>(Guu.sizes()));
    assert(std::get<1>(G.sizes()) == nel_ * nmo_);
    assert(std::get<0>(Guu.sizes()) == nu);

    ComplexType a = (walker_type == CLOSED) ? ComplexType(2.0) : ComplexType(1.0);
    Array<SPComplexType, 2> T1({(uN - u0), nw * nel_},
                               device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    Array_cref<SPComplexType, 2> Gw(make_device_ptr(G.origin()), {nw * nel_, nmo_});
    comm->barrier();

    // transposing intermediary to make dot products faster in the next step
#if defined(QMC_COMPLEX)
    ma::product(ma::T(Piu({0, nmo_}, {u0, uN})), ma::T(Gw), T1);
#else
    {
      int k0, kN;
      std::tie(k0, kN) = FairDivideBoundary(comm->rank(), nmo_, comm->size());
      ShmArray<SPComplexType, 2> TGw({nmo_, nw * nel_},
                                     shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());
      ma::transpose(Gw(Gw.extension(0), {k0, kN}), TGw.sliced(k0, kN));
      comm->barrier();
      ma::product(ma::T(Piu({0, nmo_}, {u0, uN})), TGw, T1);
    }
#endif
    // Guu[u][w] = a * sum_n T1[u][w][n] * cPua[u][n]
    using ma::Auwn_Bun_Cuw;
    Auwn_Bun_Cuw(uN - u0, nw, nel_, SPComplexType(a), T1.origin(), make_device_ptr(cPua[0][u0].origin()),
                 make_device_ptr(Guu[u0].origin()));
    comm->barrier();
  }

  // Guu[nu][nwalk]
  template<class MatA, class MatB>
  void Guu_from_full(MatA const& G, MatB&& Guu)
  {
    using std::fill_n;
    int nmo_ = int(std::get<0>(Piu.sizes()));
    int nu   = int(std::get<1>(Piu.sizes()));
    int u0, uN;
    std::tie(u0, uN) = FairDivideBoundary(comm->rank(), nu, comm->size());
    int nwalk        = G.size();

    assert(std::get<0>(G.sizes()) == std::get<1>(Guu.sizes()));
    assert(std::get<0>(Guu.sizes()) == nu);
    assert(std::get<1>(G.sizes()) == nmo_ * nmo_);

    // calculate how many walkers can be done concurrently
    long Bytes = default_buffer_size_in_MB * 1024L * 1024L;
    Bytes /= size_t(nmo_ * nu * sizeof(SPComplexType));
    int nwmax = std::min(nwalk, std::max(1, int(Bytes)));

    ComplexType a = (walker_type == CLOSED) ? ComplexType(2.0) : ComplexType(1.0);
    Array<SPComplexType, 2> T1({nwmax * nmo_, (uN - u0)},
                               device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    comm->barrier();
    fill_n(Guu[u0].origin(), nwalk * (uN - u0), ComplexType(0.0));

    APP_ABORT(" Error: Finish Guu_from_full \n\n\n");
    int iw(0);
    while (iw < nwalk)
    {
      int nw = std::min(nwmax, nwalk - iw);
      Array_cref<SPComplexType, 2> Giw(make_device_ptr(G[iw].origin()), {nw * nmo_, nmo_});
      //        ma::product(Giw,Piu({0,nmo_},{u0,uN}),T1);
      // Guu[u+u0][w] = alpha * sum_i T[w][i][u] * P[i][u]
      using ma::Awiu_Biu_Cuw;
      Awiu_Biu_Cuw(uN - u0, nw, nmo_, SPComplexType(a), T1.origin(), make_device_ptr(Piu.origin()) + u0, nu,
                   make_device_ptr(Guu[u0].origin()) + iw, nwalk);
      iw += nw;
    }
    comm->barrier();
  }

  // since this is for energy, only compact is accepted
  // Computes Guv and Guu for a set of walkers
  // rotMuv is partitioned along 'u'
  // G[w][nel*nmo]
  // Guv[w][nu][nv]
  // Guu[w][v], accumulated on this routine, sum over spin is outside
  // Tav[w][nup][nv]
  template<class MatA, class MatB, class MatC, class MatD>
  void Guv_Guu(int ispin, MatA const& G, MatB&& Guv, MatC&& Guu, MatD&& Tav, int k)
  {
    static_assert(std::decay<MatA>::type::dimensionality == 2, "Wrong dimensionality");
    static_assert(std::decay<MatB>::type::dimensionality == 3, "Wrong dimensionality");
    static_assert(std::decay<MatC>::type::dimensionality == 2, "Wrong dimensionality");
    static_assert(std::decay<MatD>::type::dimensionality == 3, "Wrong dimensionality");
    int nmo_ = int(std::get<0>(rotPiu.sizes()));
    int nu   = int(std::get<0>(rotMuv.sizes())); // potentially distributed over nodes
    int nv   = int(std::get<1>(rotMuv.sizes())); // not distributed over nodes
    int nw   = int(G.size());
    assert(std::get<1>(rotPiu.sizes()) == nv);
    int v0, vN;
    std::tie(v0, vN) = FairDivideBoundary(comm->rank(), nv, comm->size());
    int k0, kN;
    std::tie(k0, kN) = FairDivideBoundary(comm->rank(), nmo_, comm->size());
    int nu0          = rotnmu0;
    SPComplexType zero(0.0, 0.0);

    // sync first
    comm->barrier();

    using const_array_ptr = boost::multi::array_ptr<SPComplexType, 2, const_sp_pointer>;
    using array_ptr       = boost::multi::array_ptr<SPComplexType, 2, sp_pointer>;
    auto Pua_ptr(&(rotcPua[k]({nu0, nu0 + nu}, {ispin * nup, nup + ispin * ndown})));

    std::vector<const_array_ptr> Gwaj;
    std::vector<decltype(&(rotPiu({0, 1}, {0, 1})))> Pjv;
    std::vector<decltype(&(Tav[0]({0, 1}, {0, 1})))> Twav;
    std::vector<decltype(Pua_ptr)> Pua;
    std::vector<decltype(&(Guv[0]({0, 1}, {0, 1})))> Gwuv;

    Gwaj.reserve(nw);
    Pjv.reserve(nw);
    Twav.reserve(nw);
    Pua.reserve(nw);
    Gwuv.reserve(nw);

    for (int iw = 0; iw < nw; ++iw)
    {
      Gwaj.emplace_back(make_device_ptr(G[iw].origin()) + ispin * nup * nmo_, iextensions<2u>{nelec[ispin], nmo_});
      Pjv.emplace_back(&(rotPiu({0, nmo_}, {v0, vN})));
      Twav.emplace_back(&(Tav[iw]({0, nelec[ispin]}, {v0, vN})));
      Pua.emplace_back(Pua_ptr);
      Gwuv.emplace_back(&(Guv[iw]({0, nu}, {v0, vN})));
    }
    // T[w][a][v] = sum_v G[w][a][j] * rotcPua[j][v]
#if defined(QMC_COMPLEX)
    ma::BatchedProduct('N', 'N', Gwaj, Pjv, Twav);
#else
    ShmArray<SPComplexType, 3> Gja({nw, nmo_, nelec[ispin]},
                                   shm_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    Array_ref<SPComplexType, 3> Tva(make_device_ptr(Guv.origin()), {nw, nv, nelec[ispin]});
    std::vector<array_ptr> Gwja;
    std::vector<decltype(&(Tva[0]({0, 1}, {0, 1})))> Twva;
    Twva.reserve(nw);
    Gwja.reserve(nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      Twva.emplace_back(&(Tva[iw]({v0, vN}, {0, nelec[ispin]})));
      Gwja.emplace_back(make_device_ptr(Gja[iw].origin()), iextensions<2u>{nmo_, nelec[ispin]});
      ma::transpose((*(Gwaj[iw]))({0, nelec[ispin]}, {k0, kN}), (*(Gwja[iw])).sliced(k0, kN));
    }
    comm->barrier();
    ma::BatchedProduct('T', 'N', Pjv, Gwja, Twva);
    for (int iw = 0; iw < nw; ++iw)
      ma::transpose(*Twva[iw], *Twav[iw]);
#endif
    comm->barrier();
    // G[w][u][v] = sum_a rotcPua[u][a] * T[w][a][v]
    ma::BatchedProduct('N', 'N', Pua, Twav, Gwuv);
    comm->barrier();

    // Gwv = Gwvv, in range v={nu0,nu0+nu}
    using ma::Aijk_Bkj_Cik;
    using ma::get_diagonal_strided;
    //  needs distribution
    if (comm->root())
    {
      get_diagonal_strided(Guv({0, nw}, {0, nu}, {nu0, nu0 + nu}), Guu({0, nw}, {nu0, nu0 + nu}));
      // dispatch these through ma_blas_extensions!!!
      // Gwv = sum_a Twav Pva
      if (nu0 > 0) // calculate Guu from u={0,nu0}
        Aijk_Bkj_Cik(nw, nelec[ispin], nu0, make_device_ptr(Tav.origin()), Tav.stride(1), Tav.stride(0),
                     make_device_ptr(rotcPua[k].origin()), rotcPua[k].stride(0), make_device_ptr(Guu.origin()), nv);
      if (nu0 + nu < nv) // calculate Guu from u={nu0+nu,nv}
        Aijk_Bkj_Cik(nw, nelec[ispin], nv - nu0 - nu, make_device_ptr(Tav.origin()) + nu0 + nu, Tav.stride(1),
                     Tav.stride(0), make_device_ptr(rotcPua[k][nu0 + nu].origin()), rotcPua[k].stride(0),
                     make_device_ptr(Guu.origin()) + nu0 + nu, nv);
    }
    comm->barrier();
  }

  /*
    // since this is for energy, only compact is accepted
    // Computes Guv and Guu for a single walker
    // As opposed to the other Guu routines,
    //  this routine expects G for the walker in matrix form
    // rotMuv is partitioned along 'u'
    // G[nel][nmo]
    // Guv[nu][nu]
    // Guu[u]: summed over spin
    // T1[nel_][nu]
    template<class MatA, class MatB, class MatC, class MatD>
    void Guv_Guu2(MatA const& G, MatB&& Guv, MatC&& Guu, MatD&& T1, int k) {

      static_assert(std::decay<MatA>::type::dimensionality == 2, "Wrong dimensionality");
      static_assert(std::decay<MatB>::type::dimensionality == 2, "Wrong dimensionality");
      static_assert(std::decay<MatC>::type::dimensionality == 1, "Wrong dimensionality");
      static_assert(std::decay<MatD>::type::dimensionality == 2, "Wrong dimensionality");
      int nmo_ = int(rotPiu.size(0));
      int nu = int(rotMuv.size(0));  // potentially distributed over nodes
      int nv = int(rotMuv.size(1));  // not distributed over nodes
      assert(rotPiu.size(1) == nv);
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int nu0 = rotnmu0; 
      ComplexType zero(0.0,0.0);

      assert(Guu.size(0) == nv);
      assert(Guv.size(0) == nu);
      assert(Guv.size(1) == nv);

      // sync first
      comm->barrier();
      int nel_ = (walker_type==CLOSED)?nup:(nup+ndown);
      assert(G.size(0) == size_t(nel_));
      assert(G.size(1) == size_t(nmo_));
      assert(T1.size(0) == size_t(nel_));
      assert(T1.size(1) == size_t(nv));

      ma::product(G,rotPiu({0,nmo_},{v0,vN}),
                  T1(T1.extension(0),{v0,vN}));
      // This operation might benefit from a 2-D work distribution
      ma::product(rotcPua[k].sliced(nu0,nu0+nu),
                  T1(T1.extension(0),{v0,vN}),
                  Guv(Guv.extension(0),{v0,vN}));
      for(int v=v0; v<vN; ++v)
        if( v < nu0 || v >= nu0+nu ) {
          Guu[v] = ma::dot(rotcPua[k][v],T1(T1.extension(0),v)); 
        } else
         Guu[v] = Guv[v-nu0][v];
      comm->barrier();
    }
*/
protected:
  communicator* comm;

  DeviceBufferManager device_buffer_manager;
  LocalTGBufferManager shm_buffer_manager;

  long default_buffer_size_in_MB = 4L * 1024L;

  int NMO, nup, ndown;
  int nelec[2];

  int nmu0, gnmu, rotnmu0, grotnmu;

  WALKER_TYPES walker_type;

  // bare one body hamiltonian
  mpi3CMatrix hij;

  // (potentially half rotated) one body hamiltonian
  nodeArray<ComplexType, 2> haj;

  /************************************************/
  // Used in the calculation of the energy
  // Coulomb matrix elements of interpolating vectors
  nodeArray<SPValueType, 2> rotMuv;

  // Orbitals at interpolating points
  nodeArray<SPValueType, 2> rotPiu;

  // Half-rotated Orbitals at interpolating points
  std::vector<nodeArray<SPComplexType, 2>> rotcPua;
  /************************************************/

  /************************************************/
  // Following 3 used in calculation of vbias and vHS
  // Cholesky factorization of Muv
  nodeArray<SPValueType, 2> Luv;

  // Orbitals at interpolating points
  nodeArray<SPValueType, 2> Piu;

  // Half-rotated Orbitals at interpolating points
  std::vector<nodeArray<SPComplexType, 2>> cPua;
  /************************************************/

  // one-body piece of Hamiltonian factorization
  mpi3CMatrix vn0;

  ValueType E0;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
