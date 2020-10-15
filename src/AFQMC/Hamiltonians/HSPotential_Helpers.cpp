#include <cstdlib>
#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>
#include <numeric>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include "Configuration.h"

#include "AFQMC/config.h"
#include "HSPotential_Helpers.h"
#include "Utilities/FairDivide.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/Numerics/detail/utilities.hpp"

/********************************************************************
*  You get 2 potentials per Cholesky std::vector
*
*    vn(+-)_{i,k} = 0.5*( L^n_{i,k} +- ma::conj(L^n_{k,i}) )
********************************************************************/

namespace qmcplusplus
{
namespace afqmc
{
namespace HamHelper
{
namespace
{
//inline double const& ma::conj(double const& d){return d;}
//inline float const& ma::conj(float const& f){return f;}

void count_over_cholvec(double cut,
                        std::vector<std::size_t>& count,
                        int c0,
                        int c1,
                        SpVType_shm_csr_matrix::reference const& Lik,
                        SpVType_shm_csr_matrix::reference const& Lki)
{
  assert(c1 >= c0);
  if (c0 == c1)
    return;
  auto cik     = std::lower_bound(to_address(Lik.non_zero_indices2_data()),
                              to_address((Lik.non_zero_indices2_data() + Lik.num_non_zero_elements())), c0);
  auto cik_end = std::lower_bound(cik, to_address((Lik.non_zero_indices2_data() + Lik.num_non_zero_elements())), c1);
  auto vik     = to_address(Lik.non_zero_values_data()) + std::distance(to_address(Lik.non_zero_indices2_data()), cik);

  auto cki     = std::lower_bound(to_address(Lki.non_zero_indices2_data()),
                              to_address((Lki.non_zero_indices2_data() + Lki.num_non_zero_elements())), c0);
  auto cki_end = std::lower_bound(cki, to_address((Lki.non_zero_indices2_data() + Lki.num_non_zero_elements())), c1);
  auto vki     = to_address(Lki.non_zero_values_data()) + std::distance(to_address(Lki.non_zero_indices2_data()), cki);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  using ma::conj;
  using std::abs;
  using std::size_t;
  while (cik != cik_end && cki != cki_end)
  {
    if (*cik == *cki)
    { // both Lik and Lki have components on Chol Vec *cik==*cki
      if (abs(*vik + ma::conj(*vki)) > cut)
        count[2 * (*cik)] += size_t(2); // Lik + Lki* and Lki + Lik*
      if (abs(*vik - ma::conj(*vki)) > cut)
        count[2 * (*cik) + 1] += size_t(2); // Lik - Lki* and Lki - Lik*
      ++cik;
      ++vik;
      ++cki;
      ++vki;
    }
    else if (*cik < *cki)
    { // not on the same chol vector, only operate on the smallest
      if (abs(*vik) > cut)
      {
        count[2 * (*cik)] += size_t(2); // Lik + 0
#if defined(QMC_COMPLEX)
        count[2 * (*cik) + 1] += size_t(2); // Lik - 0
#endif
      }
      ++cik;
      ++vik;
    }
    else
    {
      if (abs(*vki) > cut)
      {
        count[2 * (*cki)] += size_t(2); // Lki + 0
#if defined(QMC_COMPLEX)
        count[2 * (*cki) + 1] += size_t(2); // Lki - 0
#endif
      }
      ++cki;
      ++vki;
    }
  }
  // either cik or cki are at end, check if any terms are missing
  while (cik != cik_end)
  {
    if (abs(*vik) > cut)
    {
      count[2 * (*cik)] += size_t(2); // Lik + 0
#if defined(QMC_COMPLEX)
      count[2 * (*cik) + 1] += size_t(2); // Lik - 0
#endif
    }
    ++cik;
    ++vik;
  }
  while (cki != cki_end)
  {
    if (abs(*vki) > cut)
    {
      count[2 * (*cki)] += size_t(2); // Lki + 0
#if defined(QMC_COMPLEX)
      count[2 * (*cki) + 1] += size_t(2); // Lki - 0
#endif
    }
    ++cki;
    ++vki;
  }
}

void count_over_cholvec(double cut,
                        std::vector<std::size_t>& count,
                        int c0,
                        int c1,
                        SpVType_shm_csr_matrix::reference const& Lii)
{
  assert(c1 >= c0);
  if (c0 == c1)
    return;
  auto ci     = std::lower_bound(to_address(Lii.non_zero_indices2_data()),
                             to_address((Lii.non_zero_indices2_data() + Lii.num_non_zero_elements())), c0);
  auto ci_end = std::lower_bound(ci, to_address((Lii.non_zero_indices2_data() + Lii.num_non_zero_elements())), c1);
  auto vi     = to_address(Lii.non_zero_values_data()) + std::distance(to_address(Lii.non_zero_indices2_data()), ci);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  using ma::conj;
  using std::abs;
  using std::size_t;
  while (ci != ci_end)
  {
    if (abs(*vi + ma::conj(*vi)) > cut)
      ++count[2 * (*ci)]; // Lii + Lii*
#if defined(QMC_COMPLEX)
    if (abs(*vi - ma::conj(*vi)) > cut)
      ++count[2 * (*ci) + 1]; // Lii - Lii*
#endif
    ++ci;
    ++vi;
  }
}

// In this case, c0/c1 refer to the ranges in the expanded list of CVs, e.g. 2*n/2*n+1
void count_nnz(double cut,
               std::size_t& nik,
               std::size_t& nki,
               int c0,
               int c1,
               SpVType_shm_csr_matrix::reference const& Lik,
               SpVType_shm_csr_matrix::reference const& Lki)
{
  using ma::conj;
  using std::abs;
  using std::size_t;
  assert(c1 >= c0);
  nik = nki = size_t(0);
  if (c0 == c1)
    return;
  auto cik = std::lower_bound(to_address(Lik.non_zero_indices2_data()),
                              to_address((Lik.non_zero_indices2_data() + Lik.num_non_zero_elements())), c0 / 2);
  auto cik_end =
      std::lower_bound(cik, to_address((Lik.non_zero_indices2_data() + Lik.num_non_zero_elements())), (c1 + 1) / 2);
  auto vik = to_address(Lik.non_zero_values_data()) + std::distance(to_address(Lik.non_zero_indices2_data()), cik);

  auto cki = std::lower_bound(to_address(Lki.non_zero_indices2_data()),
                              to_address((Lki.non_zero_indices2_data() + Lki.num_non_zero_elements())), c0 / 2);
  auto cki_end =
      std::lower_bound(cki, to_address((Lki.non_zero_indices2_data() + Lki.num_non_zero_elements())), (c1 + 1) / 2);
  auto vki = to_address(Lki.non_zero_values_data()) + std::distance(to_address(Lki.non_zero_indices2_data()), cki);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  while (cik != cik_end && cki != cki_end)
  {
    if (*cik == *cki)
    {                                         // both Lik and Lki have components on Chol Vec *cik==*cki
      if (abs(*vik + ma::conj(*vki)) > cut && // Lik + Lki* and Lki + Lik*
          2 * (*cik) >= c0 && 2 * (*cik) < c1)
      {
        nik++;
        nki++;
      }
#if defined(QMC_COMPLEX)
      if (abs(*vik - ma::conj(*vki)) > cut && // Lik - Lki* and Lki - Lik*
          2 * (*cik) + 1 >= c0 && 2 * (*cik) + 1 < c1)
      {
        nik++;
        nki++;
      }
#endif
      ++cik;
      ++vik;
      ++cki;
      ++vki;
    }
    else if (*cik < *cki)
    { // not on the same chol vector, only operate on the smallest
#if defined(QMC_COMPLEX)
      if (abs(*vik) > cut)
      {
        if (2 * (*cik) >= c0 && 2 * (*cik) < c1)
        {
          ++nik;
          ++nki;
        }
        if (2 * (*cik) + 1 >= c0 && 2 * (*cik) + 1 < c1)
        {
          ++nik;
          ++nki;
        }
      }
#else
      if (abs(*vik) > cut && 2 * (*cik) >= c0 && 2 * (*cik) < c1)
      {
        ++nik;
        ++nki;
      }
#endif
      ++cik;
      ++vik;
    }
    else
    {
#if defined(QMC_COMPLEX)
      if (abs(*vki) > cut)
      {
        if (2 * (*cki) >= c0 && 2 * (*cki) < c1)
        {
          ++nik;
          ++nki;
        }
        if (2 * (*cki) + 1 >= c0 && 2 * (*cki) + 1 < c1)
        {
          ++nik;
          ++nki;
        }
      }
#else
      if (abs(*vki) > cut && 2 * (*cki) >= c0 && 2 * (*cki) < c1)
      {
        ++nik;
        ++nki;
      }
#endif
      ++cki;
      ++vki;
    }
  }
  while (cik != cik_end)
  {
#if defined(QMC_COMPLEX)
    if (abs(*vik) > cut)
    {
      if (2 * (*cik) >= c0 && 2 * (*cik) < c1)
      {
        ++nik;
        ++nki;
      }
      if (2 * (*cik) + 1 >= c0 && 2 * (*cik) + 1 < c1)
      {
        ++nik;
        ++nki;
      }
    }
#else
    if (abs(*vik) > cut && 2 * (*cik) >= c0 && 2 * (*cik) < c1)
    {
      ++nik;
      ++nki;
    }
#endif
    ++cik;
    ++vik;
  }
  while (cki != cki_end)
  {
#if defined(QMC_COMPLEX)
    if (abs(*vki) > cut)
    {
      if (2 * (*cki) >= c0 && 2 * (*cki) < c1)
      {
        ++nik;
        ++nki;
      }
      if (2 * (*cki) + 1 >= c0 && 2 * (*cki) + 1 < c1)
      {
        ++nik;
        ++nki;
      }
    }
#else
    if (abs(*vki) > cut && 2 * (*cki) >= c0 && 2 * (*cki) < c1)
    {
      ++nik;
      ++nki;
    }
#endif
    ++cki;
    ++vki;
  }
}

void count_nnz(double cut, size_t& ni, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lii)
{
  using ma::conj;
  using std::abs;
  using std::size_t;
  assert(c1 >= c0);
  ni = size_t(0);
  if (c0 == c1)
    return;
  auto ci = std::lower_bound(to_address(Lii.non_zero_indices2_data()),
                             to_address((Lii.non_zero_indices2_data() + Lii.num_non_zero_elements())), c0 / 2);
  auto ci_end =
      std::lower_bound(ci, to_address((Lii.non_zero_indices2_data() + Lii.num_non_zero_elements())), (c1 + 1) / 2);
  auto vi = to_address(Lii.non_zero_values_data()) + std::distance(to_address(Lii.non_zero_indices2_data()), ci);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  while (ci != ci_end)
  {
    if (abs(*vi + ma::conj(*vi)) > cut && 2 * (*ci) >= c0 && 2 * (*ci) < c1)
      ++ni; // Lii + Lii*
#if defined(QMC_COMPLEX)
    if (abs(*vi - ma::conj(*vi)) > cut && 2 * (*ci) + 1 >= c0 && 2 * (*ci) + 1 < c1)
      ++ni; // Lii - Lii*
#endif
    ++ci;
    ++vi;
  }
}

void add_to_vn(SpVType_shm_csr_matrix& vn,
               std::vector<int> const& map_,
               double cut,
               int ik,
               int c0,
               int c1,
               SpVType_shm_csr_matrix::reference const& Lii)
{
  using ma::conj;
  using std::abs;
  using std::size_t;
  assert(c1 >= c0);
  if (c0 == c1)
    return;
  auto ci = std::lower_bound(to_address(Lii.non_zero_indices2_data()),
                             to_address((Lii.non_zero_indices2_data() + Lii.num_non_zero_elements())), c0 / 2);
  auto ci_end =
      std::lower_bound(ci, to_address((Lii.non_zero_indices2_data() + Lii.num_non_zero_elements())), (c1 + 1) / 2);
  auto vi = to_address(Lii.non_zero_values_data()) + std::distance(to_address(Lii.non_zero_indices2_data()), ci);

  int c_origin = vn.global_origin()[1];
  SPValueType half(0.5);
  SPComplexType im(0.0, 1.0);
  while (ci != ci_end)
  {
    if (abs(*vi + ma::conj(*vi)) > cut && 2 * (*ci) >= c0 && 2 * (*ci) < c1)
    {
      assert(map_[2 * (*ci)] >= 0);
      assert(map_[2 * (*ci)] - c_origin < vn.size(1));
      vn.emplace_back({ik, (map_[2 * (*ci)] - c_origin)},
                      static_cast<SPValueType>(half * (*vi + ma::conj(*vi)))); // Lii + Lii*
    }
#if defined(QMC_COMPLEX)
    if (abs(*vi - ma::conj(*vi)) > cut && 2 * (*ci) + 1 >= c0 && 2 * (*ci) + 1 < c1)
    {
      assert(map_[2 * (*ci) + 1] >= 0);
      assert(map_[2 * (*ci) + 1] - c_origin < vn.size(1));
      vn.emplace_back({ik, (map_[2 * (*ci) + 1] - c_origin)},
                      static_cast<SPValueType>(half * im * (*vi - ma::conj(*vi)))); // Lii - Lii*
    }
#endif
    ++ci;
    ++vi;
  }
}

void add_to_vn(SpVType_shm_csr_matrix& vn,
               std::vector<int> const& map_,
               double cut,
               int ik,
               int ki,
               int c0,
               int c1,
               SpVType_shm_csr_matrix::reference const& Lik,
               SpVType_shm_csr_matrix::reference const& Lki)
{
  using ma::conj;
  using std::abs;
  using std::size_t;
  assert(c1 >= c0);
  if (c0 == c1)
    return;
  auto cik = std::lower_bound(to_address(Lik.non_zero_indices2_data()),
                              to_address((Lik.non_zero_indices2_data() + Lik.num_non_zero_elements())), c0 / 2);
  auto cik_end =
      std::lower_bound(cik, to_address((Lik.non_zero_indices2_data() + Lik.num_non_zero_elements())), (c1 + 1) / 2);
  auto vik = to_address(Lik.non_zero_values_data()) + std::distance(to_address(Lik.non_zero_indices2_data()), cik);

  auto cki = std::lower_bound(to_address(Lki.non_zero_indices2_data()),
                              to_address((Lki.non_zero_indices2_data() + Lki.num_non_zero_elements())), c0 / 2);
  auto cki_end =
      std::lower_bound(cki, to_address((Lki.non_zero_indices2_data() + Lki.num_non_zero_elements())), (c1 + 1) / 2);
  auto vki = to_address(Lki.non_zero_values_data()) + std::distance(to_address(Lki.non_zero_indices2_data()), cki);

  SPComplexType im(0.0, 1.0);
  SPValueType half(0.5);
  int c_origin = vn.global_origin()[1];
  while (cik != cik_end && cki != cki_end)
  {
    if (*cik == *cki)
    { // both Lik and Lki have components on Chol Vec *cik==*cki
      if (abs(*vik + ma::conj(*vki)) > cut)
      { // Lik + Lki* and Lki + Lik*
        if (2 * (*cik) >= c0 && 2 * (*cik) < c1)
        {
          assert(map_[2 * (*cik)] >= 0);
          assert(map_[2 * (*cik)] - c_origin < vn.size(1));
          vn.emplace_back({ik, (map_[2 * (*cik)] - c_origin)},
                          static_cast<SPValueType>(half * (*vik + ma::conj(*vki)))); // Lik + Lki*
          vn.emplace_back({ki, (map_[2 * (*cki)] - c_origin)},
                          static_cast<SPValueType>(half * (*vki + ma::conj(*vik)))); // Lki + Lik*
        }
      }
#if defined(QMC_COMPLEX)
      if (abs(*vik - ma::conj(*vki)) > cut)
      { // Lik - Lki* and Lki - Lik*
        if (2 * (*cik) + 1 >= c0 && 2 * (*cik) + 1 < c1)
        {
          assert(map_[2 * (*cik) + 1] >= 0);
          assert(map_[2 * (*cik) + 1] - c_origin < vn.size(1));
          vn.emplace_back({ik, (map_[2 * (*cik) + 1] - c_origin)},
                          static_cast<SPValueType>(half * im * (*vik - ma::conj(*vki)))); // Lik - Lki*
          vn.emplace_back({ki, (map_[2 * (*cki) + 1] - c_origin)},
                          static_cast<SPValueType>(half * im * (*vki - ma::conj(*vik)))); // Lki - Lik*
        }
      }
#endif
      ++cik;
      ++vik;
      ++cki;
      ++vki;
    }
    else if (*cik < *cki)
    { // not on the same chol vector, only operate on the smallest
      if (abs(*vik) > cut)
      {
        if (2 * (*cik) >= c0 && 2 * (*cik) < c1)
        {
          assert(map_[2 * (*cik)] >= 0);
          assert(map_[2 * (*cik)] - c_origin < vn.size(1));
          vn.emplace_back({ik, (map_[2 * (*cik)] - c_origin)},
                          static_cast<SPValueType>(half * (*vik))); // Lik + 0
          vn.emplace_back({ki, (map_[2 * (*cik)] - c_origin)},
                          static_cast<SPValueType>(half * ma::conj(*vik))); // Lik + 0
        }
#if defined(QMC_COMPLEX)
        if (2 * (*cik) + 1 >= c0 && 2 * (*cik) + 1 < c1)
        {
          assert(map_[2 * (*cik) + 1] >= 0);
          assert(map_[2 * (*cik) + 1] - c_origin < vn.size(1));
          vn.emplace_back({ik, (map_[2 * (*cik) + 1] - c_origin)},
                          static_cast<SPValueType>(half * im * (*vik))); // Lik - 0
          vn.emplace_back({ki, (map_[2 * (*cik) + 1] - c_origin)},
                          static_cast<SPValueType>(-half * im * ma::conj(*vik))); // Lik - 0
        }
#endif
      }
      ++cik;
      ++vik;
    }
    else
    {
      if (abs(*vki) > cut)
      {
        if (2 * (*cki) >= c0 && 2 * (*cki) < c1)
        {
          assert(map_[2 * (*cki)] >= 0);
          assert(map_[2 * (*cki)] - c_origin < vn.size(1));
          vn.emplace_back({ik, (map_[2 * (*cki)] - c_origin)},
                          static_cast<SPValueType>(half * ma::conj(*vki))); // Lki + 0
          vn.emplace_back({ki, (map_[2 * (*cki)] - c_origin)},
                          static_cast<SPValueType>(half * (*vki))); // Lki + 0
        }
#if defined(QMC_COMPLEX)
        if (2 * (*cki) + 1 >= c0 && 2 * (*cki) + 1 < c1)
        {
          assert(map_[2 * (*cki) + 1] >= 0);
          assert(map_[2 * (*cki) + 1] - c_origin < vn.size(1));
          vn.emplace_back({ik, (map_[2 * (*cki) + 1] - c_origin)},
                          static_cast<SPValueType>(-half * im * ma::conj(*vki))); // Lki - 0
          vn.emplace_back({ki, (map_[2 * (*cki) + 1] - c_origin)},
                          static_cast<SPValueType>(half * im * (*vki))); // Lki - 0
        }
#endif
      }
      ++cki;
      ++vki;
    }
  }
  while (cik != cik_end)
  {
    if (abs(*vik) > cut)
    {
      if (2 * (*cik) >= c0 && 2 * (*cik) < c1)
      {
        assert(map_[2 * (*cik)] >= 0);
        assert(map_[2 * (*cik)] - c_origin < vn.size(1));
        vn.emplace_back({ik, (map_[2 * (*cik)] - c_origin)},
                        static_cast<SPValueType>(half * (*vik))); // Lik + 0
        vn.emplace_back({ki, (map_[2 * (*cik)] - c_origin)},
                        static_cast<SPValueType>(half * ma::conj(*vik))); // Lik + 0
      }
#if defined(QMC_COMPLEX)
      if (2 * (*cik) + 1 >= c0 && 2 * (*cik) + 1 < c1)
      {
        assert(map_[2 * (*cik) + 1] >= 0);
        assert(map_[2 * (*cik) + 1] - c_origin < vn.size(1));
        vn.emplace_back({ik, (map_[2 * (*cik) + 1] - c_origin)},
                        static_cast<SPValueType>(half * im * (*vik))); // Lik - 0
        vn.emplace_back({ki, (map_[2 * (*cik) + 1] - c_origin)},
                        static_cast<SPValueType>(-half * im * ma::conj(*vik))); // Lik - 0
      }
#endif
    }
    ++cik;
    ++vik;
  }
  while (cki != cki_end)
  {
    if (abs(*vki) > cut)
    {
      if (2 * (*cki) >= c0 && 2 * (*cki) < c1)
      {
        assert(map_[2 * (*cki)] >= 0);
        assert(map_[2 * (*cki)] - c_origin < vn.size(1));
        vn.emplace_back({ik, (map_[2 * (*cki)] - c_origin)},
                        static_cast<SPValueType>(half * ma::conj(*vki))); // Lki + 0
        vn.emplace_back({ki, (map_[2 * (*cki)] - c_origin)},
                        static_cast<SPValueType>(half * (*vki))); // Lki + 0
      }
#if defined(QMC_COMPLEX)
      if (2 * (*cki) + 1 >= c0 && 2 * (*cki) + 1 < c1)
      {
        assert(map_[2 * (*cki) + 1] >= 0);
        assert(map_[2 * (*cki) + 1] - c_origin < vn.size(1));
        vn.emplace_back({ik, (map_[2 * (*cki) + 1] - c_origin)},
                        static_cast<SPValueType>(-half * im * ma::conj(*vki))); // Lki - 0
        vn.emplace_back({ki, (map_[2 * (*cki) + 1] - c_origin)},
                        static_cast<SPValueType>(half * im * (*vki))); // Lki - 0
      }
#endif
    }
    ++cki;
    ++vki;
  }
}

} // namespace

std::vector<std::size_t> count_nnz_per_cholvec(double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO)
{
  if (TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: count_nnz_per_cholvec is not designed for distributed CholMat. \n");

  if (V2.size(0) != NMO * NMO)
    APP_ABORT(" Error in count_nnz_per_cholvec: V2.size(0) ! NMO*NMO \n");
  int ik0, ikN;
  // only upper triangular part since both ik/ki are processed simultaneously
  std::tie(ik0, ikN) = FairDivideBoundary(TG.getGlobalRank(), int(NMO * (NMO + 1) / 2), TG.getGlobalSize());

  if (cut < 1e-12)
    cut = 1e-12;
  std::vector<std::size_t> counts(2 * V2.size(1));
  int nvecs = V2.size(1);

  for (int i = 0, cnt = 0; i < NMO; i++)
    for (int k = i; k < NMO; k++, cnt++)
    {
      if (cnt < ik0)
        continue;
      if (cnt >= ikN)
        break;
      if (i == k)
        count_over_cholvec(cut, counts, 0, nvecs, V2[i * NMO + k]);
      else
        count_over_cholvec(cut, counts, 0, nvecs, V2[i * NMO + k], V2[k * NMO + i]);
    }
  TG.Global().all_reduce_in_place_n(counts.begin(), counts.size(), std::plus<>());

  return counts;
}

std::vector<std::size_t> count_nnz_per_ik(double cut,
                                          TaskGroup_& TG,
                                          SpVType_shm_csr_matrix& V2,
                                          int NMO,
                                          int cv0,
                                          int cvN)
{
  assert(cv0 >= 0 && cvN <= 2 * V2.size(1));
  if (TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: count_nnz_per_ik is not designed for distributed CholMat. \n");

  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if (V2.size(0) != NMO * NMO)
    APP_ABORT(" Error in count_nnz_per_cholvec: V2.size(0) ! NMO*NMO \n");
  int ik0, ikN;
  // only upper triangular part since both ik/ki are processed simultaneously
  // it is possible to ssubdivide this over equivalent nodes and then reduce over them
  std::tie(ik0, ikN) = FairDivideBoundary(coreid, int(NMO * (NMO + 1) / 2), ncores);

  if (cut < 1e-12)
    cut = 1e-12;
  std::vector<std::size_t> counts(V2.size(0));

  for (int i = 0, cnt = 0; i < NMO; i++)
    for (int k = i; k < NMO; k++, cnt++)
    {
      if (cnt < ik0)
        continue;
      if (cnt >= ikN)
        break;
      if (i == k)
        count_nnz(cut, counts[i * NMO + k], cv0, cvN, V2[i * NMO + k]);
      else
        count_nnz(cut, counts[i * NMO + k], counts[k * NMO + i], cv0, cvN, V2[i * NMO + k], V2[k * NMO + i]);
    }
  TG.Node().all_reduce_in_place_n(counts.begin(), counts.size(), std::plus<>());
  // if divided over equivalent nodes, reduce over the Core_Equivalent communicator
  //TG.EqvCore().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());

  return counts;
}

void generateHSPotential(SpVType_shm_csr_matrix& vn,
                         std::vector<int> const& map_,
                         double cut,
                         TaskGroup_& TG,
                         SpVType_shm_csr_matrix& V2,
                         int NMO,
                         int cv0,
                         int cvN)
{
  assert(cv0 >= 0 && cvN <= 2 * V2.size(1));
  if (TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: generateHSPotential is not designed for distributed CholMat. \n");

  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if (V2.size(0) != NMO * NMO)
    APP_ABORT(" Error in generateHSPotential: V2.size(0) ! NMO*NMO \n");
  int ik0, ikN;
  // only upper triangular part since both ik/ki are processed simultaneously
  // it is possible to ssubdivide this over equivalent nodes and then reduce over them
  std::tie(ik0, ikN) = FairDivideBoundary(coreid, int(NMO * (NMO + 1) / 2), ncores);

  if (cut < 1e-12)
    cut = 1e-12;

  for (int i = 0, cnt = 0; i < NMO; i++)
    for (int k = i; k < NMO; k++, cnt++)
    {
      if (cnt < ik0)
        continue;
      if (cnt >= ikN)
        break;
      if (i == k)
        add_to_vn(vn, map_, cut, i * NMO + k, cv0, cvN, V2[i * NMO + k]);
      else
        add_to_vn(vn, map_, cut, i * NMO + k, k * NMO + i, cv0, cvN, V2[i * NMO + k], V2[k * NMO + i]);
    }
  // if distributed over equivalent nodes, perform communication step here

  TG.node_barrier();
}

} // namespace HamHelper

} // namespace afqmc

} // namespace qmcplusplus
