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

#ifndef QMCPLUSPLUS_AFQMC_HAMHELPER_GENERATEHIJKL_HPP
#define QMCPLUSPLUS_AFQMC_HAMHELPER_GENERATEHIJKL_HPP

#include <cstdlib>
#include <algorithm>
#include <complex>
#include <map>
#include <vector>
#include <numeric>
#include <functional>

#include "Configuration.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Utilities/afqmc_TTI.hpp"
#include "AFQMC/Matrix/array_partition.hpp"
#include "AFQMC/Hamiltonians/createHamiltonian_Helper.hpp"

namespace qmcplusplus
{
namespace afqmc
{
namespace HamHelper
{
template<class Container, class Hijkl_op>
void generateHijkl(WALKER_TYPES walker_type,
                   bool addCoulomb,
                   TaskGroup_& TG,
                   Container& Vijkl,
                   IndexType NMO,
                   std::map<IndexType, std::pair<bool, IndexType>>& occ_a,
                   std::map<IndexType, std::pair<bool, IndexType>>& occ_b,
                   Hijkl_op& H,
                   bool reserve_to_fit_ = true,
                   bool global_reserve  = true,
                   RealType const cut   = 1e-6)
{
  using ma::conj;
  using std::abs;
  using std::sqrt;

  if (not addCoulomb)
    APP_ABORT("Error: addCoulomb=false not yet implemented in generateHijkl. \n");

  if (walker_type == NONCOLLINEAR)
  {
    APP_ABORT("Error: GHF density matrix only implemented with spinRestricted integrals. \n");
  }

  long npr = TG.getGlobalSize(), rk = TG.getGlobalRank();

  bool distribute_Ham = TG.getNumberOfTGs() > 1;
  if (distribute_Ham)
    APP_ABORT(" Distributed hamiltonian not yet implemented in generateHijkl. \n");

  ValueType zero = ValueType(0);

  std::vector<s4D<ValueType>> vs4D;
  vs4D.reserve(24);

  long cnter = 0;
  ValueType J1, J2, J3, J1a = zero, J2a = zero, J3a = zero;

  boost::multi::array<ValueType, 2> DiagHam({NMO, NMO});
  for (IndexType i = 0; i < NMO; i++)
    for (IndexType k = i; k < NMO; k++, cnter++)
    {
      if (cnter % npr != rk)
        continue;
      DiagHam[i][k] = H.H(i, k, k, i);
      DiagHam[k][i] = DiagHam[i][k];
#if defined(QMC_COMPLEX)
      if (imag(DiagHam[i][k]) > 1e-8)
      {
        app_error() << " Error: Found complex diagonal on hamiltonian. " << i << " " << k << " " << DiagHam[i][k]
                    << std::endl;
        APP_ABORT("Error: Found complex diagonal on hamiltonian.");
      }
      ComplexType Hik = H.H(i, k, i, k);
      if (imag(Hik) > cut)
      {
        app_error() << " Error: Problems with factorized Hamiltonian <ik|ik> should be real: " << i << " " << k << " "
                    << Hik << std::endl;
        APP_ABORT("Error: Problems with factorized hamiltonian.");
      }
#endif
    }
  TG.Global().all_reduce_in_place_n(DiagHam.origin(), DiagHam.num_elements(), std::plus<>());

  int occi, occj, occk, occl;
  if (reserve_to_fit_)
  {
    std::vector<std::size_t> sz_local(Vijkl.size(0));

    cnter = 0;
    // Approximation: (similar to non factorized case)
    //   - if <ij|kl>  <  cut, set to zero
    for (IndexType i = 0; i < NMO; i++)
    {
      for (IndexType j = i; j < NMO; j++, cnter++)
      {
        if (cnter % npr != rk)
          continue;
        occi = (occ_a[i].first || occ_b[i].first) ? 1 : 0;
        occj = (occ_a[j].first || occ_b[j].first) ? 1 : 0;
        for (IndexType k = j; k < NMO; k++)
        {
          occk = (occ_a[k].first || occ_b[k].first) ? 1 : 0;
          for (IndexType l = k; l < NMO; l++)
          {
            occl = (occ_a[l].first || occ_b[l].first) ? 1 : 0;
            if (occi + occj + occk + occl < 2)
              continue;

            // NOTE NOTE NOTE:
            // <ik|ik> can get a complex component due to truncation error, eliminate it here

            J1 = J2 = J3 = zero;
            // J1 = <ij|kl>
            // J1 < sqrt( <ik|ki> * <lj|jl> )
            if (sqrt(abs(DiagHam[i][k] * DiagHam[l][j])) > cut)
            {
              J1 = H.H(i, j, k, l);
#if defined(QMC_COMPLEX)
              if (i == k && j == l)
                J1 = ValueType(J1.real(), 0);
#endif
            }

            // J2 = <ij|lk>
            if (i == j || l == k)
            {
              J2 = J1;
            }
            else
            {
              if (sqrt(abs(DiagHam[i][l] * DiagHam[k][j])) > cut)
              {
                J2 = H.H(i, j, l, k);
#if defined(QMC_COMPLEX)
                if (i == l && j == k)
                  J2 = ValueType(J2.real(), 0);
#endif
              }
            }

            // J3 = <ik|jl>
            if (j == k)
            {
              J3 = J1;
            }
            else if (i == l)
            {
              J3 = ma::conj(J1);
            }
            else
            {
              if (sqrt(abs(DiagHam[i][j] * DiagHam[l][k])) > cut)
              {
                J3 = H.H(i, k, j, l);
#if defined(QMC_COMPLEX)
                if (i == j && k == l)
                  J3 = ValueType(J3.real(), 0);
#endif
              }
            }

#if defined(QMC_COMPLEX)
            J1a = J2a = J3a = zero;

            //  J2a = <ik|lj>
            if (l == j)
            {
              J2a = J3;
            }
            else if (i == k)
            {
              J2a = J3;
            }
            else if (k == j)
            {
              J2a = J2;
            }
            else if (i == l)
            {
              J2a = ma::conj(J2);
            }
            else
            {
              if (sqrt(abs(DiagHam[i][l] * DiagHam[j][k])) > cut)
                J2a = H.H(i, k, l, j);
            }

            //  J3a = <il|jk>
            if (l == j)
            {
              J3a = J2;
            }
            else if (i == k)
            {
              J3a = ma::conj(J2);
            }
            else if (k == l)
            {
              J3a = J3;
            }
            else if (i == j)
            {
              J3a = ma::conj(J3);
            }
            else
            {
              if (sqrt(abs(DiagHam[i][j] * DiagHam[k][l])) > cut)
                J3a = H.H(i, l, j, k);
            }

            //  For complex, there are 3 extra non-symmetric terms:
            //  J1a = <il|kj>
            if (k == l)
            {
              J1a = J2a;
            }
            else if (i == j)
            {
              J1a = ma::conj(J2a);
            }
            else if (j == k)
            {
              J1a = J3a;
            }
            else if (i == l)
            {
              J1a = J3a;
            }
            else if (l == j)
            {
              J1a = J1;
            }
            else if (i == k)
            {
              J1a = ma::conj(J1);
            }
            else
            {
              if (sqrt(abs(DiagHam[i][k] * DiagHam[j][l])) > cut)
                J1a = H.H(i, l, k, j);
            }
#endif

            if (abs(J1) < cut && abs(J2) < cut && abs(J3) < cut && abs(J1a) < cut && abs(J2a) < cut && abs(J3a) < cut)
              continue;

            vs4D.clear();
            if (walker_type == CLOSED)
            {
              find_all_contributions_to_hamiltonian_closed_shell(NMO, i, j, k, l, J1, J2, J3, J1a, J2a, J3a, cut, vs4D);
              count_allowed_terms(NMO, sz_local, vs4D, occ_a);
            }
            else if (walker_type == COLLINEAR)
            {
              APP_ABORT(" Error: Finish. \n");
              find_all_contributions_to_hamiltonian_collinear(NMO, i, j, k, l, J1, J2, J3, J1a, J2a, J3a, cut, vs4D);
              count_allowed_terms(sz_local, vs4D, occ_a, occ_b);
            }
            else if (walker_type == NONCOLLINEAR)
            {
              APP_ABORT(" Error: Finish. \n");
              find_all_contributions_to_hamiltonian_ghf(NMO, i, j, k, l, J1, J2, J3, J1a, J2a, J3a, cut, vs4D);
              count_allowed_terms(sz_local, vs4D, occ_a, occ_b);
            }
            else
            {
              APP_ABORT(" Error: Unknown walker type in generateHijkl.\n");
            }
          }
        }
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

    app_log() << "  Number of terms in Vijkl (generateHijkl): \n"
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

  } // reserve_to_fit_

  cnter = 0;
  // Approximation: (similar to non factorized case)
  //   - if <ij|kl>  <  cut, set to zero
  for (IndexType i = 0; i < NMO; i++)
  {
    for (IndexType j = i; j < NMO; j++, cnter++)
    {
      if (cnter % npr != rk)
        continue;
      occi = (occ_a[i].first || occ_b[i].first) ? 1 : 0;
      occj = (occ_a[j].first || occ_b[j].first) ? 1 : 0;
      for (IndexType k = j; k < NMO; k++)
      {
        occk = (occ_a[k].first || occ_b[k].first) ? 1 : 0;
        for (IndexType l = k; l < NMO; l++)
        {
          occl = (occ_a[l].first || occ_b[l].first) ? 1 : 0;
          if (occi + occj + occk + occl < 2)
            continue;

          J1 = J2 = J3 = zero;
          // J1 = <ij|kl>
          // J1 < sqrt( <ik|ki> * <lj|jl> )
          if (sqrt(abs(DiagHam[i][k] * DiagHam[l][j])) > cut)
          {
            J1 = H.H(i, j, k, l);
#if defined(QMC_COMPLEX)
            if (i == k && j == l)
              J1 = ValueType(J1.real(), 0);
#endif
          }

          // J2 = <ij|lk>
          if (i == j || l == k)
          {
            J2 = J1;
          }
          else
          {
            if (sqrt(abs(DiagHam[i][l] * DiagHam[k][j])) > cut)
            {
              J2 = H.H(i, j, l, k);
#if defined(QMC_COMPLEX)
              if (i == l && j == k)
                J2 = ValueType(J2.real(), 0);
#endif
            }
          }

          // J3 = <ik|jl>
          if (j == k)
          {
            J3 = J1;
          }
          else if (i == l)
          {
            J3 = ma::conj(J1);
          }
          else
          {
            if (sqrt(abs(DiagHam[i][j] * DiagHam[l][k])) > cut)
            {
              J3 = H.H(i, k, j, l);
#if defined(QMC_COMPLEX)
              if (i == j && k == l)
                J3 = ValueType(J3.real(), 0);
#endif
            }
          }

#if defined(QMC_COMPLEX)
          J1a = J2a = J3a = zero;

          //  J2a = <ik|lj>
          if (l == j)
          {
            J2a = J3;
          }
          else if (i == k)
          {
            J2a = J3;
          }
          else if (k == j)
          {
            J2a = J2;
          }
          else if (i == l)
          {
            J2a = ma::conj(J2);
          }
          else
          {
            if (sqrt(abs(DiagHam[i][l] * DiagHam[j][k])) > cut)
              J2a = H.H(i, k, l, j);
          }

          //  J3a = <il|jk>
          if (l == j)
          {
            J3a = J2;
          }
          else if (i == k)
          {
            J3a = ma::conj(J2);
          }
          else if (k == l)
          {
            J3a = J3;
          }
          else if (i == j)
          {
            J3a = ma::conj(J3);
          }
          else
          {
            if (sqrt(abs(DiagHam[i][j] * DiagHam[k][l])) > cut)
              J3a = H.H(i, l, j, k);
          }

          //  For complex, there are 3 extra non-symmetric terms:
          //  J1a = <il|kj>
          if (k == l)
          {
            J1a = J2a;
          }
          else if (i == j)
          {
            J1a = ma::conj(J2a);
          }
          else if (j == k)
          {
            J1a = J3a;
          }
          else if (i == l)
          {
            J1a = J3a;
          }
          else if (l == j)
          {
            J1a = J1;
          }
          else if (i == k)
          {
            J1a = ma::conj(J1);
          }
          else
          {
            if (sqrt(abs(DiagHam[i][k] * DiagHam[j][l])) > cut)
              J1a = H.H(i, l, k, j);
          }

#endif

          if (abs(J1) < cut && abs(J2) < cut && abs(J3) < cut && abs(J1a) < cut && abs(J2a) < cut && abs(J3a) < cut)
            continue;

          vs4D.clear();
          if (walker_type == CLOSED)
          {
            find_all_contributions_to_hamiltonian_closed_shell(NMO, i, j, k, l, J1, J2, J3, J1a, J2a, J3a, cut, vs4D);
            add_allowed_terms(NMO, vs4D, occ_a, Vijkl);
          }
          else if (walker_type == COLLINEAR)
          {
            APP_ABORT("Finish.\n");
            find_all_contributions_to_hamiltonian_collinear(NMO, i, j, k, l, J1, J2, J3, J1a, J2a, J3a, cut, vs4D);
            add_allowed_terms(NMO, vs4D, occ_a, occ_b, Vijkl, true, false);
          }
          else if (walker_type == NONCOLLINEAR)
          {
            APP_ABORT("Finish.\n");
            find_all_contributions_to_hamiltonian_ghf(NMO, i, j, k, l, J1, J2, J3, J1a, J2a, J3a, cut, vs4D);
          }
          else
          {
            APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
          }
        }
      }
    }
  }
  TG.global_barrier();
}

} // namespace HamHelper
} // namespace afqmc
} // namespace qmcplusplus
#endif
