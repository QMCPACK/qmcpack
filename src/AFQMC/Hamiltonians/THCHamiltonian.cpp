#include <cstdlib>
#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <numeric>

#include "Configuration.h"
#include "type_traits/container_traits_multi.h"
#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "THCHamiltonian.h"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"
#include "AFQMC/HamiltonianOperations/THCOpsIO.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// Right now, cutvn, cutvn2, TGprop and TGwfn are completely ignored.
// Note: addCoulomb only has meaning on the sparse hamiltonians, not in THC
HamiltonianOperations THCHamiltonian::getHamiltonianOperations(bool pureSD,
                                                               bool addCoulomb,
                                                               WALKER_TYPES type,
                                                               std::vector<PsiT_Matrix>& PsiT,
                                                               double cutvn,
                                                               double cutv2,
                                                               TaskGroup_& TGprop,
                                                               TaskGroup_& TGwfn,
                                                               hdf_archive& hdf_restart)
{
  // hack until parallel hdf is in place
  bool write_hdf = false;
  if (TGwfn.Global().root())
    write_hdf = !hdf_restart.closed();
  //  if(TGwfn.Global().root()) write_hdf = (hdf_restart.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  if (type == COLLINEAR)
    assert(PsiT.size() % 2 == 0);
  int ndet = ((type != COLLINEAR) ? (PsiT.size()) : (PsiT.size() / 2));

  if (ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  // this communicator needs to be used for data structures that are distributed
  // over multiple nodes in a TG. When built with accelerator support, multiple
  // members of the TG will reside in the same node, so the Node() communicator will lead
  // to wrong results. For host only builds, this will just be a copy of Node.
  auto distNode(TG.Node().split(TGwfn.getLocalGroupNumber(), TG.Node().rank()));
  if (TGwfn.getLocalGroupNumber() != TGprop.getLocalGroupNumber())
  {
    // relax this later
    app_error() << " Error: nnodes in wavefunction must match value in Propagator. \n" << std::endl;
    APP_ABORT("");
  }

  size_t gnmu, grotnmu, nmu, rotnmu, nmu0, nmuN, rotnmu0, rotnmuN;
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
    dump.push("THC", false);
  }
  if (TG.Global().root())
  {
    std::vector<int> Idata(3);
    if (!dump.readEntry(Idata, "dims"))
    {
      app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                  << " Problems reading dims. \n";
      APP_ABORT("");
    }
    if (Idata[0] != NMO)
    {
      app_error() << " ERROR: NMO differs from value in integral file. \n";
      APP_ABORT(" Error: NMO differs from value in integral file. \n");
    }
    gnmu    = size_t(Idata[1]);
    grotnmu = size_t(Idata[2]);
  }
  TG.Global().broadcast_value(gnmu);
  TG.Global().broadcast_value(grotnmu);

  // setup partition, in general matrices are partitioned asize_t 'u'
  {
    int node_number            = TGwfn.getLocalGroupNumber();
    int nnodes_prt_TG          = TGwfn.getNGroupsPerTG();
    std::tie(rotnmu0, rotnmuN) = FairDivideBoundary(size_t(node_number), grotnmu, size_t(nnodes_prt_TG));
    rotnmu                     = rotnmuN - rotnmu0;

    node_number          = TGprop.getLocalGroupNumber();
    nnodes_prt_TG        = TGprop.getNGroupsPerTG();
    std::tie(nmu0, nmuN) = FairDivideBoundary(size_t(node_number), gnmu, size_t(nnodes_prt_TG));
    nmu                  = nmuN - nmu0;
  }

  using shmVMatrix    = boost::multi::array<ValueType, 2, shared_allocator<ValueType>>;
  using shmCMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using Vmatrix_ref   = boost::multi::array_ref<ValueType, 2>;
  using Cmatrix_ref   = boost::multi::array_ref<ComplexType, 2>;
  using shmSPVMatrix  = boost::multi::array<SPValueType, 2, shared_allocator<SPValueType>>;
  using shmSPCMatrix  = boost::multi::array<SPComplexType, 2, shared_allocator<SPComplexType>>;
  using SPVMatrix_ref = boost::multi::array_ref<SPValueType, 2>;
  using SPCMatrix_ref = boost::multi::array_ref<SPComplexType, 2>;

  // Until I figure something else, rotPiu and rotcPua are not distributed because a full copy is needed
  // distribution:  size,  global,  offset
  //   - rotMuv:    {rotnmu,grotnmu},{grotnmu,grotnmu},{rotnmu0,0}
  //   - rotPiu:    {size_t(NMO),grotnmu},{size_t(NMO),grotnmu},{0,0}
  //   - rotcPua    {grotnmu,nel_},{grotnmu,nel_},{0,0}
  //   - Piu:       {size_t(NMO),nmu},{size_t(NMO),gnmu},{0,nmu0}
  //   - Luv:       {nmu,gnmu},{gnmu,gnmu},{nmu0,0}
  //   - cPua       {nmu,nel_},{gnmu,nel_},{nmu0,0}
  //
  size_t nel_ = PsiT[0].size(0) + ((type == CLOSED) ? 0 : (PsiT[1].size(0)));
  shmSPVMatrix rotMuv({static_cast<shmSPVMatrix::size_type>(rotnmu), static_cast<shmSPVMatrix::size_type>(grotnmu)}, shared_allocator<SPValueType>{distNode});
  shmSPVMatrix rotPiu({NMO, static_cast<shmSPVMatrix::size_type>(grotnmu)}, shared_allocator<SPValueType>{distNode});
  std::vector<shmSPCMatrix> rotcPua;
  rotcPua.reserve(ndet);
  for (int i = 0; i < ndet; i++)
    rotcPua.emplace_back(shmSPCMatrix({static_cast<shmSPCMatrix::size_type>(grotnmu), static_cast<shmSPCMatrix::size_type>(nel_)}, shared_allocator<SPComplexType>{distNode}));
  shmSPVMatrix Piu({NMO, static_cast<shmSPVMatrix::size_type>(nmu)}, shared_allocator<SPValueType>{distNode});
  shmSPVMatrix Luv({static_cast<shmSPVMatrix::size_type>(nmu), static_cast<shmSPVMatrix::size_type>(gnmu)}, shared_allocator<SPValueType>{distNode});
  // right now only 1 reader. Use hyperslabs and parallel io later
  // read Half transformed first
  if (distNode.root())
  {
    using ma::conj;
    /***************************************/
    // read full matrix, not distributed for now
    if (!dump.readEntry(rotPiu, "HalfTransformedFullOrbitals"))
    {
      app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                  << " Problems reading HalfTransformedFullOrbitals. \n";
      APP_ABORT("");
    }
    /***************************************/
    hyperslab_proxy<shmSPVMatrix, 2> hslab(rotMuv, std::array<size_t, 2>{grotnmu, grotnmu},
                                           std::array<size_t, 2>{rotnmu, grotnmu}, std::array<size_t, 2>{rotnmu0, 0});
    if (!dump.readEntry(hslab, "HalfTransformedMuv"))
    {
      app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                  << " Problems reading HalfTransformedMuv. \n";
      APP_ABORT("");
    }
    /***************************************/
  }
  TG.global_barrier();

  if (distNode.root())
  {
    /***************************************/
    hyperslab_proxy<shmSPVMatrix, 2> hslab(Piu, std::array<size_t, 2>{size_t(NMO), gnmu},
                                           std::array<size_t, 2>{size_t(NMO), nmu}, std::array<size_t, 2>{0, nmu0});
    if (!dump.readEntry(hslab, "Orbitals"))
    {
      app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                  << " Problems reading Orbitals. \n";
      APP_ABORT("");
    }
    /***************************************/
    hyperslab_proxy<shmSPVMatrix, 2> hslab2(Luv, std::array<size_t, 2>{gnmu, gnmu}, std::array<size_t, 2>{nmu, gnmu},
                                            std::array<size_t, 2>{nmu0, 0});
    if (!dump.readEntry(hslab2, "Luv"))
    {
      app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                  << " Problems reading Luv. \n";
      APP_ABORT("");
    }
    /***************************************/
  }
  TG.global_barrier();

  shmCMatrix v0({NMO, NMO}, shared_allocator<ComplexType>{TG.Node()});
  if (TGprop.getNGroupsPerTG() > 1)
  {
    // MAM: doing this in single precision to be consistent with non-distributed case
    // very inefficient, find better work distribution for calculation of v0
    // that doesn't so much temporary space!!!
    boost::multi::array<SPValueType, 2> v0_({NMO, NMO});
    // TOO MUCH MEMORY, FIX FIX FIX!!!
    shmSPVMatrix Piu__({NMO, static_cast<shmSPVMatrix::size_type>(gnmu)}, shared_allocator<SPValueType>{TG.Node()});
    shmSPVMatrix Luv__({static_cast<shmSPVMatrix::size_type>(gnmu), static_cast<shmSPVMatrix::size_type>(gnmu)}, shared_allocator<SPValueType>{TG.Node()});
    if (TG.Node().root())
    {
      if (!dump.readEntry(Luv__, "Luv"))
      {
        app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                    << " Problems reading Orbitals. \n";
        APP_ABORT("");
      }
      if (!dump.readEntry(Piu__, "Orbitals"))
      {
        app_error() << " Error in THCHamiltonian::getHamiltonianOperations():"
                    << " Problems reading Orbitals. \n";
        APP_ABORT("");
      }
    }
    TG.Node().barrier();

    using ma::conj;
    using ma::H;
    using ma::T;
    size_t c0, cN, nc;
    std::tie(c0, cN) = FairDivideBoundary(size_t(TG.Global().rank()), gnmu, size_t(TG.Global().size()));
    nc               = cN - c0;
    boost::multi::array<SPValueType, 2> Tuv({static_cast<boost::multi::size_t>(gnmu), static_cast<boost::multi::size_t>(nc)});
    boost::multi::array<SPValueType, 2> Muv({static_cast<boost::multi::size_t>(gnmu), static_cast<boost::multi::size_t>(nc)});

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work
    ma::product(Luv__, H(Luv__.sliced(c0, cN)), Muv);

    // since generating v0 takes some effort and temporary space,
    // v0(i,l) = -0.5*sum_j <i,j|j,l>
    //         = -0.5 sum_j,u,v ma::conj(Piu(i,u)) ma::conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v ma::conj(Piu(i,u)) W(u,v) Piu(l,u), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) ma::conj(Piu(j,v))
    ma::product(H(Piu__), Piu__({0, long(NMO)}, {long(c0), long(cN)}), Tuv);
    auto itM = Muv.origin();
    auto itT = Tuv.origin();
    for (size_t i = 0; i < Muv.num_elements(); ++i, ++itT, ++itM)
      *(itT) = ma::conj(*itT) * (*itM);
    boost::multi::array<SPValueType, 2> T_({static_cast<boost::multi::size_t>(std::get<1>(Tuv.sizes())), NMO});
    ma::product(T(Tuv), H(Piu__), T_);
    ma::product(SPValueType(-0.5), T(T_), T(Piu__({0, long(NMO)}, {long(c0), long(cN)})), SPValueType(0.0), v0_);

    // reduce over Global
    TG.Global().all_reduce_in_place_n(v0_.origin(), v0_.num_elements(), std::plus<>());
    if (TG.Node().root())
    {
      copy_n_cast(v0_.origin(), NMO * NMO, to_address(v0.origin()));
#if defined(MIXED_PRECISION)
      // MAM: Since Muv gets large, might have problems with the check for hermicity below
      // fixing here
      for (int i = 0; i < NMO; i++)
        for (int j = i + 1; j < NMO; j++)
        {
          v0[i][j] = 0.5 * (v0[i][j] + ma::conj(v0[j][i]));
          v0[j][i] = ma::conj(v0[i][j]);
        }
#endif
    }
    TG.Node().barrier();
  }
  else
  {
    // very inefficient, find better work distribution for calculation of v0
    // that doesn't so much temporary space!!!
    boost::multi::array<SPValueType, 2> v0_({NMO, NMO});
    // very simple partitioning until something more sophisticated is in place!!!
    using ma::conj;
    using ma::H;
    using ma::T;
    size_t c0, cN, nc;
    std::tie(c0, cN) = FairDivideBoundary(size_t(TG.Global().rank()), gnmu, size_t(TG.Global().size()));
    nc               = cN - c0;
    boost::multi::array<SPValueType, 2> Tuv({static_cast<boost::multi::size_t>(gnmu), static_cast<boost::multi::size_t>(nc)});
    boost::multi::array<SPValueType, 2> Muv({static_cast<boost::multi::size_t>(gnmu), static_cast<boost::multi::size_t>(nc)});

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work
    ma::product(Luv, H(Luv.sliced(c0, cN)), Muv);

    // since generating v0 takes some effort and temporary space,
    // v0(i,l) = -0.5*sum_j <i,j|j,l>
    //         = -0.5 sum_j,u,v ma::conj(Piu(i,u)) ma::conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v ma::conj(Piu(i,u)) W(u,v) Piu(l,v), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) ma::conj(Piu(j,v))
    ma::product(H(Piu), Piu({0, long(NMO)}, {long(c0), long(cN)}), Tuv);
    auto itM = Muv.origin();
    auto itT = Tuv.origin();
    for (size_t i = 0; i < Muv.num_elements(); ++i, ++itT, ++itM)
      *(itT) = ma::conj(*itT) * (*itM);
    boost::multi::array<SPValueType, 2> T_({static_cast<boost::multi::size_t>(std::get<1>(Tuv.sizes())), NMO});
    ma::product(T(Tuv), H(Piu), T_);
    ma::product(SPValueType(-0.5), T(T_), T(Piu({0, long(NMO)}, {long(c0), long(cN)})), SPValueType(0.0), v0_);

    // reduce over Global
    TG.Global().all_reduce_in_place_n(v0_.origin(), v0_.num_elements(), std::plus<>());
    if (TG.Node().root())
    {
      copy_n_cast(v0_.origin(), NMO * NMO, to_address(v0.origin()));
#if defined(MIXED_PRECISION)
      // MAM: Since Muv gets large, might have problems with the check for hermicity below
      // fixing here
      for (int i = 0; i < NMO; i++)
        for (int j = i + 1; j < NMO; j++)
        {
          v0[i][j] = 0.5 * (v0[i][j] + ma::conj(v0[j][i]));
          v0[j][i] = ma::conj(v0[i][j]);
        }
#endif
    }
    TG.Node().barrier();
  }
  TG.global_barrier();

  long naea_ = PsiT[0].size(0);
  long naeb_ = ((type == COLLINEAR) ? PsiT.back().size(0) : 0);

  // half-rotated Pia
  std::vector<shmSPCMatrix> cPua;
  cPua.reserve(ndet);
  for (int i = 0; i < ndet; i++)
    cPua.emplace_back(shmSPCMatrix({static_cast<shmSPCMatrix::size_type>(nmu), static_cast<shmSPCMatrix::size_type>(nel_)}, shared_allocator<SPComplexType>{distNode}));
  if (distNode.root())
  {
    // simple
    using ma::H;
    if (type == COLLINEAR)
    {
      boost::multi::array<SPComplexType, 2> A({NMO, naea_});
      boost::multi::array<SPComplexType, 2> B({NMO, naeb_});
      for (int i = 0; i < ndet; i++)
      {
        // cPua = H(Piu) * ma::conj(A)
        ma::Matrix2MA('T', PsiT[2 * i], A);
        ma::product(H(Piu), A, cPua[i]({0, long(nmu)}, {0, long(naea_)}));
        ma::product(H(rotPiu), A, rotcPua[i]({0, long(grotnmu)}, {0, long(naea_)}));
        ma::Matrix2MA('T', PsiT[2 * i + 1], B);
        ma::product(H(Piu), B, cPua[i]({0, long(nmu)}, {naea_, long(nel_)}));
        ma::product(H(rotPiu), B, rotcPua[i]({0, long(grotnmu)}, {naea_, long(nel_)}));
      }
    }
    else
    {
      boost::multi::array<SPComplexType, 2> A({static_cast<boost::multi::size_t>(PsiT[0].size(1)), static_cast<boost::multi::size_t>(PsiT[0].size(0))});
      for (int i = 0; i < ndet; i++)
      {
        ma::Matrix2MA('T', PsiT[i], A);
        // cPua = H(Piu) * ma::conj(A)
        ma::product(H(Piu), A, cPua[i]);
        ma::product(H(rotPiu), A, rotcPua[i]);
      }
    }
  }
  TG.node_barrier();

  ValueType E0 = OneBodyHamiltonian::NuclearCoulombEnergy + OneBodyHamiltonian::FrozenCoreEnergy;

  shmCMatrix hij({ndet, (naea_ + naeb_) * NMO}, shared_allocator<ComplexType>{TG.Node()});
  shmCMatrix H1_({NMO, NMO}, shared_allocator<ComplexType>{TG.Node()});
  if (TG.Node().root())
  {
    // dense one body hamiltonian
    copy_n_cast((OneBodyHamiltonian::H1).origin(), NMO * NMO, to_address(H1_.origin()));
    int skp = ((type == COLLINEAR) ? 1 : 0);
    for (int n = 0, nd = 0; n < ndet; ++n, nd += (skp + 1))
    {
      check_wavefunction_consistency(type, &PsiT[nd], &PsiT[nd + skp], NMO, naea_, naeb_);
      auto hij_(rotateHij(type, &PsiT[nd], &PsiT[nd + skp], H1_));
      std::copy_n(hij_.origin(), hij_.num_elements(), to_address(hij[n].origin()));
    }
  }
  TG.Node().barrier();

  if (write_hdf)
    writeTHCOps(hdf_restart, type, NMO, naea_, naeb_, nmu0, rotnmu0, ndet, TGprop, TGwfn, H1_, rotPiu, rotMuv, Piu, Luv,
                v0, E0);

  if (distNode.root())
  {
    dump.pop();
    dump.pop();
    dump.close();
  }

  return HamiltonianOperations(THCOps(TGwfn.TG_local(), NMO, naea_, naeb_, type, nmu0, rotnmu0, std::move(H1_),
                                      std::move(hij), std::move(rotMuv), std::move(rotPiu), std::move(rotcPua),
                                      std::move(Luv), std::move(Piu), std::move(cPua), std::move(v0), E0));
}


} // namespace afqmc
} // namespace qmcplusplus
