#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<vector>
#include<numeric>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "Configuration.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HSPotential_Helpers.h"
#include <Utilities/FairDivide.h>
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"

/********************************************************************
*  You get 2 potentials per Cholesky std::vector
*
*    vn(+-)_{i,k} = 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )
********************************************************************/

namespace qmcplusplus
{

namespace afqmc
{

namespace HamHelper
{

namespace
{

inline double const& conj(double const& d){return d;}
inline float const& conj(float const& f){return f;}

void count_over_cholvec(double cut, std::vector<std::size_t>& count, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lik, SpVType_shm_csr_matrix::reference const& Lki)
{
  assert(c1>=c0);
  if(c0==c1) return;
  auto cik = std::lower_bound( std::addressof(*Lik.non_zero_indices2_data()),
                               std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                               c0);
  auto cik_end = std::lower_bound( cik,
                                   std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                                   c1);
  auto vik = std::addressof(*Lik.non_zero_values_data()) +
                    std::distance(std::addressof(*Lik.non_zero_indices2_data()),cik);

  auto cki = std::lower_bound( std::addressof(*Lki.non_zero_indices2_data()),
                               std::addressof(*(Lki.non_zero_indices2_data()
                                               +Lki.num_non_zero_elements())),
                               c0);
  auto cki_end = std::lower_bound( cki,
                                   std::addressof(*(Lki.non_zero_indices2_data()
                                               +Lki.num_non_zero_elements())),
                                   c1);
  auto vki = std::addressof(*Lki.non_zero_values_data()) +
                    std::distance(std::addressof(*Lki.non_zero_indices2_data()),cki);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  using std::conj;
  using std::abs;
  using std::size_t;
  while( cik!=cik_end && cki!=cki_end ) {
    if( *cik == *cki ) { // both Lik and Lki have components on Chol Vec *cik==*cki
      if(abs(*vik + conj(*vki)) > cut) count[2*(*cik)]+=size_t(2); // Lik + Lki* and Lki + Lik*
      if(abs(*vik - conj(*vki)) > cut) count[2*(*cik)+1]+=size_t(2); // Lik - Lki* and Lki - Lik*
      ++cik;
      ++vik;
      ++cki;
      ++vki;
    } else if( *cik < *cki )  { // not on the same chol vector, only operate on the smallest
      if(abs(*vik) > cut) {
        count[2*(*cik)]+=size_t(2);   // Lik + 0
#if defined(QMC_COMPLEX)
        count[2*(*cik)+1]+=size_t(2); // Lik - 0
#endif
      }
      ++cik;
      ++vik;
    } else {
      if(abs(*vki) > cut) {
        count[2*(*cki)]+=size_t(2);   // Lki + 0
#if defined(QMC_COMPLEX)
        count[2*(*cki)+1]+=size_t(2); // Lki - 0
#endif
      }
      ++cki;
      ++vki;
    }
  }
  // either cik or cki are at end, check if any terms are missing
  while( cik!=cik_end ) {
    if(abs(*vik) > cut) {
      count[2*(*cik)]+=size_t(2);   // Lik + 0
#if defined(QMC_COMPLEX)
      count[2*(*cik)+1]+=size_t(2); // Lik - 0
#endif
    }
    ++cik;
    ++vik;
  }
  while( cki!=cki_end ) {
    if(abs(*vki) > cut) {
      count[2*(*cki)]+=size_t(2);   // Lki + 0
#if defined(QMC_COMPLEX)
      count[2*(*cki)+1]+=size_t(2); // Lki - 0
#endif
    }
    ++cki;
    ++vki;
  }
}

void count_over_cholvec(double cut, std::vector<std::size_t>& count, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lii)
{
  assert(c1>=c0);
  if(c0==c1) return;
  auto ci = std::lower_bound( std::addressof(*Lii.non_zero_indices2_data()),
                              std::addressof(*(Lii.non_zero_indices2_data()
                                               +Lii.num_non_zero_elements())),
                              c0);
  auto ci_end = std::lower_bound( ci,
                                  std::addressof(*(Lii.non_zero_indices2_data()
                                               +Lii.num_non_zero_elements())),
                                  c1);
  auto vi = std::addressof(*Lii.non_zero_values_data()) +
                    std::distance(std::addressof(*Lii.non_zero_indices2_data()),ci);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  using std::conj;
  using std::abs;
  using std::size_t;
  while( ci!=ci_end ) {
    if(abs(*vi + conj(*vi)) > cut) ++count[2*(*ci)];   // Lii + Lii*
#if defined(QMC_COMPLEX)
    if(abs(*vi - conj(*vi)) > cut) ++count[2*(*ci)+1]; // Lii - Lii*
#endif
    ++ci;
    ++vi;
  }
}

// In this case, c0/c1 refer to the ranges in the expanded list of CVs, e.g. 2*n/2*n+1
void count_nnz(double cut, std::size_t& nik, std::size_t& nki, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lik, SpVType_shm_csr_matrix::reference const& Lki)
{
  using std::conj;
  using std::abs;
  using std::size_t;
  assert(c1>=c0);
  nik = nki = size_t(0);
  if(c0==c1) return;
  auto cik = std::lower_bound( std::addressof(*Lik.non_zero_indices2_data()),
                               std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                               c0/2);
  auto cik_end = std::lower_bound( cik,
                                   std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                                   (c1+1)/2);
  auto vik = std::addressof(*Lik.non_zero_values_data()) +
                    std::distance(std::addressof(*Lik.non_zero_indices2_data()),cik);

  auto cki = std::lower_bound( std::addressof(*Lki.non_zero_indices2_data()),
                               std::addressof(*(Lki.non_zero_indices2_data()
                                               +Lki.num_non_zero_elements())),
                               c0/2);
  auto cki_end = std::lower_bound( cki,
                                   std::addressof(*(Lki.non_zero_indices2_data()
                                               +Lki.num_non_zero_elements())),
                                   (c1+1)/2);
  auto vki = std::addressof(*Lki.non_zero_values_data()) +
                    std::distance(std::addressof(*Lki.non_zero_indices2_data()),cki);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  while( cik!=cik_end && cki!=cki_end ) {
    if( *cik == *cki ) { // both Lik and Lki have components on Chol Vec *cik==*cki
      if(abs(*vik + conj(*vki)) > cut && // Lik + Lki* and Lki + Lik*
          2*(*cik)>=c0 && 2*(*cik)<c1) { nik++; nki++; }
#if defined(QMC_COMPLEX)
      if(abs(*vik - conj(*vki)) > cut && // Lik - Lki* and Lki - Lik*
          2*(*cik)+1>=c0 && 2*(*cik)+1<c1) { nik++; nki++; }
#endif
      ++cik;
      ++vik;
      ++cki;
      ++vki;
    } else if( *cik < *cki )  { // not on the same chol vector, only operate on the smallest
#if defined(QMC_COMPLEX)
      if(abs(*vik) > cut) {
        if(2*(*cik)>=c0 && 2*(*cik)<c1) { ++nik; ++nki; }
        if(2*(*cik)+1>=c0 && 2*(*cik)+1<c1) { ++nik; ++nki; }
      }
#else
      if(abs(*vik) > cut && 2*(*cik)>=c0 && 2*(*cik)<c1) { ++nik; ++nki; }
#endif
      ++cik;
      ++vik;
    } else {
#if defined(QMC_COMPLEX)
      if(abs(*vki) > cut) {
        if(2*(*cki)>=c0 && 2*(*cki)<c1) { ++nik; ++nki; }
        if(2*(*cki)+1>=c0 && 2*(*cki)+1<c1) { ++nik; ++nki; }
      }
#else
      if(abs(*vki) > cut && 2*(*cki)>=c0 && 2*(*cki)<c1) { ++nik; ++nki; }
#endif
      ++cki;
      ++vki;
    }
  }
  while( cik!=cik_end ) {
#if defined(QMC_COMPLEX)
    if(abs(*vik) > cut) {
      if(2*(*cik)>=c0 && 2*(*cik)<c1) { ++nik; ++nki; }
      if(2*(*cik)+1>=c0 && 2*(*cik)+1<c1) { ++nik; ++nki; }
    }
#else
    if(abs(*vik) > cut && 2*(*cik)>=c0 && 2*(*cik)<c1) { ++nik; ++nki; }
#endif
    ++cik;
    ++vik;
  }
  while( cki!=cki_end ) {
#if defined(QMC_COMPLEX)
    if(abs(*vki) > cut) {
      if(2*(*cki)>=c0 && 2*(*cki)<c1) { ++nik; ++nki; }
      if(2*(*cki)+1>=c0 && 2*(*cki)+1<c1) { ++nik; ++nki; }
    }
#else
    if(abs(*vki) > cut && 2*(*cki)>=c0 && 2*(*cki)<c1) { ++nik; ++nki; }
#endif
    ++cki;
    ++vki;
  }
}

void count_nnz(double cut, size_t& ni, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lii)
{
  using std::conj;
  using std::abs;
  using std::size_t;
  assert(c1>=c0);
  ni=size_t(0);
  if(c0==c1) return;
  auto ci = std::lower_bound( std::addressof(*Lii.non_zero_indices2_data()),
                              std::addressof(*(Lii.non_zero_indices2_data()
                                               +Lii.num_non_zero_elements())),
                              c0/2);
  auto ci_end = std::lower_bound( ci,
                                  std::addressof(*(Lii.non_zero_indices2_data()
                                               +Lii.num_non_zero_elements())),
                                  (c1+1)/2);
  auto vi = std::addressof(*Lii.non_zero_values_data()) +
                    std::distance(std::addressof(*Lii.non_zero_indices2_data()),ci);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  while( ci!=ci_end ) {
    if(abs(*vi + conj(*vi)) > cut && 2*(*ci)>=c0 && 2*(*ci)<c1) ++ni; // Lii + Lii*
#if defined(QMC_COMPLEX)
    if(abs(*vi - conj(*vi)) > cut && 2*(*ci)+1>=c0 && 2*(*ci)+1<c1) ++ni; // Lii - Lii*
#endif
    ++ci;
    ++vi;
  }
}

void add_to_vn(SpVType_shm_csr_matrix& vn, std::vector<int> const& map_, double cut, int ik, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lii)
{
  using std::conj;
  using std::abs;
  using std::size_t;
  assert(c1>=c0);
  if(c0==c1) return;
  auto ci = std::lower_bound( std::addressof(*Lii.non_zero_indices2_data()),
                              std::addressof(*(Lii.non_zero_indices2_data()
                                               +Lii.num_non_zero_elements())),
                              c0/2);
  auto ci_end = std::lower_bound( ci,
                                  std::addressof(*(Lii.non_zero_indices2_data()
                                               +Lii.num_non_zero_elements())),
                                  (c1+1)/2);
  auto vi = std::addressof(*Lii.non_zero_values_data()) +
                    std::distance(std::addressof(*Lii.non_zero_indices2_data()),ci);

  int c_origin = vn.global_origin()[1];
  ComplexType im(0.0,1.0);
  while( ci!=ci_end ) {
    if(abs(*vi + conj(*vi)) > cut && 2*(*ci)>=c0 && 2*(*ci)<c1 ) {
        assert(map_[2*(*ci)] >= 0);
        assert(map_[2*(*ci)]-c_origin < vn.shape()[1]);
        vn.emplace_back({ik,(map_[2*(*ci)]-c_origin)},
                static_cast<SPValueType>(0.5*(*vi+conj(*vi)))); // Lii + Lii*
    }
#if defined(QMC_COMPLEX)
    if(abs(*vi - conj(*vi)) > cut && 2*(*ci)+1>=c0 && 2*(*ci)+1<c1) {
        assert(map_[2*(*ci)+1] >= 0);
        assert(map_[2*(*ci)+1]-c_origin < vn.shape()[1]);
        vn.emplace_back({ik,(map_[2*(*ci)+1]-c_origin)},
                static_cast<SPValueType>(0.5*im*(*vi-conj(*vi)))); // Lii - Lii*
    }
#endif
    ++ci;
    ++vi;
  }
}

void add_to_vn(SpVType_shm_csr_matrix& vn, std::vector<int> const& map_, double cut, int ik, int ki, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lik, SpVType_shm_csr_matrix::reference const& Lki)
{
  using std::conj;
  using std::abs;
  using std::size_t;
  assert(c1>=c0);
  if(c0==c1) return;
  auto cik = std::lower_bound( std::addressof(*Lik.non_zero_indices2_data()),
                               std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                               c0/2);
  auto cik_end = std::lower_bound( cik,
                                   std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                                   (c1+1)/2);
  auto vik = std::addressof(*Lik.non_zero_values_data()) +
                    std::distance(std::addressof(*Lik.non_zero_indices2_data()),cik);

  auto cki = std::lower_bound( std::addressof(*Lki.non_zero_indices2_data()),
                               std::addressof(*(Lki.non_zero_indices2_data()
                                               +Lki.num_non_zero_elements())),
                               c0/2);
  auto cki_end = std::lower_bound( cki,
                                   std::addressof(*(Lki.non_zero_indices2_data()
                                               +Lki.num_non_zero_elements())),
                                   (c1+1)/2);
  auto vki = std::addressof(*Lki.non_zero_values_data()) +
                    std::distance(std::addressof(*Lki.non_zero_indices2_data()),cki);

  ComplexType im(0.0,1.0);
  int c_origin = vn.global_origin()[1];
  while( cik!=cik_end && cki!=cki_end ) {
    if( *cik == *cki ) { // both Lik and Lki have components on Chol Vec *cik==*cki
      if(abs(*vik + conj(*vki)) > cut) { // Lik + Lki* and Lki + Lik*
          if(2*(*cik)>=c0 && 2*(*cik)<c1) {
            assert(map_[2*(*cik)] >= 0);
            assert(map_[2*(*cik)]-c_origin < vn.shape()[1]);
            vn.emplace_back({ik,(map_[2*(*cik)]-c_origin)},
                  static_cast<SPValueType>(0.5*(*vik+conj(*vki)))); // Lik + Lki*
            vn.emplace_back({ki,(map_[2*(*cki)]-c_origin)},
                  static_cast<SPValueType>(0.5*(*vki+conj(*vik)))); // Lki + Lik*
          }
      }
#if defined(QMC_COMPLEX)
      if(abs(*vik - conj(*vki)) > cut) { // Lik - Lki* and Lki - Lik*
          if(2*(*cik)+1>=c0 && 2*(*cik)+1<c1) {
            assert(map_[2*(*cik)+1] >= 0);
            assert(map_[2*(*cik)+1]-c_origin < vn.shape()[1]);
            vn.emplace_back({ik,(map_[2*(*cik)+1]-c_origin)},
                  static_cast<SPValueType>(0.5*im*(*vik-conj(*vki)))); // Lik - Lki*
            vn.emplace_back({ki,(map_[2*(*cki)+1]-c_origin)},
                  static_cast<SPValueType>(0.5*im*(*vki-conj(*vik)))); // Lki - Lik*
        }
      }
#endif
      ++cik;
      ++vik;
      ++cki;
      ++vki;
    } else if( *cik < *cki )  { // not on the same chol vector, only operate on the smallest
      if(abs(*vik) > cut) {
          if(2*(*cik)>=c0 && 2*(*cik)<c1) {
            assert(map_[2*(*cik)] >= 0);
            assert(map_[2*(*cik)]-c_origin < vn.shape()[1]);
            vn.emplace_back({ik,(map_[2*(*cik)]-c_origin)},
                  static_cast<SPValueType>(0.5*(*vik))); // Lik + 0
            vn.emplace_back({ki,(map_[2*(*cik)]-c_origin)},
                  static_cast<SPValueType>(0.5*conj(*vik))); // Lik + 0
          }
#if defined(QMC_COMPLEX)
          if(2*(*cik)+1>=c0 && 2*(*cik)+1<c1) {
            assert(map_[2*(*cik)+1] >= 0);
            assert(map_[2*(*cik)+1]-c_origin < vn.shape()[1]);
            vn.emplace_back({ik,(map_[2*(*cik)+1]-c_origin)},
                  static_cast<SPValueType>(0.5*im*(*vik))); // Lik - 0
            vn.emplace_back({ki,(map_[2*(*cik)+1]-c_origin)},
                  static_cast<SPValueType>(-0.5*im*conj(*vik))); // Lik - 0
          }
#endif
      }
      ++cik;
      ++vik;
    } else {
      if(abs(*vki) > cut) {
          if(2*(*cki)>=c0 && 2*(*cki)<c1) {
            assert(map_[2*(*cki)] >= 0);
            assert(map_[2*(*cki)]-c_origin < vn.shape()[1]);
            vn.emplace_back({ik,(map_[2*(*cki)]-c_origin)},
                  static_cast<SPValueType>(0.5*conj(*vki))); // Lki + 0
            vn.emplace_back({ki,(map_[2*(*cki)]-c_origin)},
                  static_cast<SPValueType>(0.5*(*vki))); // Lki + 0
          }
#if defined(QMC_COMPLEX)
          if(2*(*cki)+1>=c0 && 2*(*cki)+1<c1) {
            assert(map_[2*(*cki)+1] >= 0);
            assert(map_[2*(*cki)+1]-c_origin < vn.shape()[1]);
            vn.emplace_back({ik,(map_[2*(*cki)+1]-c_origin)},
                  static_cast<SPValueType>(-0.5*im*conj(*vki))); // Lki - 0
            vn.emplace_back({ki,(map_[2*(*cki)+1]-c_origin)},
                  static_cast<SPValueType>(0.5*im*(*vki))); // Lki - 0
          }
#endif
      }
      ++cki;
      ++vki;
    }
  }
  while( cik!=cik_end ) {
    if(abs(*vik) > cut) {
        if(2*(*cik)>=c0 && 2*(*cik)<c1) {
          assert(map_[2*(*cik)] >= 0);
          assert(map_[2*(*cik)]-c_origin < vn.shape()[1]);
          vn.emplace_back({ik,(map_[2*(*cik)]-c_origin)},
                static_cast<SPValueType>(0.5*(*vik))); // Lik + 0
          vn.emplace_back({ki,(map_[2*(*cik)]-c_origin)},
                static_cast<SPValueType>(0.5*conj(*vik))); // Lik + 0
        }
#if defined(QMC_COMPLEX)
        if(2*(*cik)+1>=c0 && 2*(*cik)+1<c1) {
          assert(map_[2*(*cik)+1] >= 0);
          assert(map_[2*(*cik)+1]-c_origin < vn.shape()[1]);
          vn.emplace_back({ik,(map_[2*(*cik)+1]-c_origin)},
                static_cast<SPValueType>(0.5*im*(*vik))); // Lik - 0
          vn.emplace_back({ki,(map_[2*(*cik)+1]-c_origin)},
                static_cast<SPValueType>(-0.5*im*conj(*vik))); // Lik - 0
        }
#endif
    }
    ++cik;
    ++vik;
  }
  while( cki!=cki_end ) {
    if(abs(*vki) > cut) {
        if(2*(*cki)>=c0 && 2*(*cki)<c1) {
          assert(map_[2*(*cki)] >= 0);
          assert(map_[2*(*cki)]-c_origin < vn.shape()[1]);
          vn.emplace_back({ik,(map_[2*(*cki)]-c_origin)},
                static_cast<SPValueType>(0.5*conj(*vki))); // Lki + 0
          vn.emplace_back({ki,(map_[2*(*cki)]-c_origin)},
                static_cast<SPValueType>(0.5*(*vki))); // Lki + 0
        }
#if defined(QMC_COMPLEX)
        if(2*(*cki)+1>=c0 && 2*(*cki)+1<c1) {
          assert(map_[2*(*cki)+1] >= 0);
          assert(map_[2*(*cki)+1]-c_origin < vn.shape()[1]);
          vn.emplace_back({ik,(map_[2*(*cki)+1]-c_origin)},
                static_cast<SPValueType>(-0.5*im*conj(*vki))); // Lki - 0
          vn.emplace_back({ki,(map_[2*(*cki)+1]-c_origin)},
                static_cast<SPValueType>(0.5*im*(*vki))); // Lki - 0
        }
#endif
    }
    ++cki;
    ++vki;
  }
}

}

std::vector<std::size_t> count_nnz_per_cholvec(double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO)
{
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: count_nnz_per_cholvec is not designed for distributed CholMat. \n");

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if(V2.shape()[0] != NMO*NMO)
    APP_ABORT(" Error in count_nnz_per_cholvec: V2.shape()[0] ! NMO*NMO \n");
  int ik0, ikN;
  // only upper triangular part since both ik/ki are processed simultaneously
  std::tie(ik0, ikN) = FairDivideBoundary(TG.getGlobalRank(),int(NMO*(NMO+1)/2),TG.getGlobalSize());

  if(cut < 1e-12) cut=1e-12;
  std::vector<std::size_t> counts(2*V2.shape()[1]);
  int nvecs = V2.shape()[1];

  for(int i=0, cnt=0; i<NMO; i++)
    for(int k=i; k<NMO; k++, cnt++) {
      if(cnt < ik0) continue;
      if(cnt >= ikN) break;
      if(i==k)
        count_over_cholvec(cut,counts,0,nvecs,V2[i*NMO+k]);
      else
        count_over_cholvec(cut,counts,0,nvecs,V2[i*NMO+k],V2[k*NMO+i]);
    }
  TG.Global().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());

  return counts;
}

std::vector<std::size_t> count_nnz_per_ik(double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO, int cv0, int cvN)
{
  assert(cv0 >= 0 && cvN <= 2*V2.shape()[1]);
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: count_nnz_per_ik is not designed for distributed CholMat. \n");

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if(V2.shape()[0] != NMO*NMO)
    APP_ABORT(" Error in count_nnz_per_cholvec: V2.shape()[0] ! NMO*NMO \n");
  int ik0, ikN;
  // only upper triangular part since both ik/ki are processed simultaneously
  // it is possible to ssubdivide this over equivalent nodes and then reduce over them
  std::tie(ik0, ikN) = FairDivideBoundary(coreid,int(NMO*(NMO+1)/2),ncores);

  if(cut < 1e-12) cut=1e-12;
  std::vector<std::size_t> counts(V2.shape()[0]);

  std::size_t nik,nki;
  for(int i=0, cnt=0; i<NMO; i++)
    for(int k=i; k<NMO; k++, cnt++) {
      if(cnt < ik0) continue;
      if(cnt >= ikN) break;
      if(i==k)
        count_nnz(cut,counts[i*NMO+k],cv0,cvN,V2[i*NMO+k]);
       else
        count_nnz(cut,counts[i*NMO+k],counts[k*NMO+i],cv0,cvN,V2[i*NMO+k],V2[k*NMO+i]);
    }
  TG.Node().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());
  // if divided over equivalent nodes, reduce over the Core_Equivalent communicator
  //TG.EqvCore().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());

  return counts;
}

void generateHSPotential(SpVType_shm_csr_matrix& vn, std::vector<int> const& map_, double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO, int cv0, int cvN)
{

  assert(cv0 >= 0 && cvN <= 2*V2.shape()[1]);
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: generateHSPotential is not designed for distributed CholMat. \n");

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if(V2.shape()[0] != NMO*NMO)
    APP_ABORT(" Error in generateHSPotential: V2.shape()[0] ! NMO*NMO \n");
  int ik0, ikN;
  // only upper triangular part since both ik/ki are processed simultaneously
  // it is possible to ssubdivide this over equivalent nodes and then reduce over them
  std::tie(ik0, ikN) = FairDivideBoundary(coreid,int(NMO*(NMO+1)/2),ncores);

  if(cut < 1e-12) cut=1e-12;

  for(int i=0, cnt=0; i<NMO; i++)
    for(int k=i; k<NMO; k++, cnt++) {
      if(cnt < ik0) continue;
      if(cnt >= ikN) break;
      if(i==k)
        add_to_vn(vn,map_,cut,i*NMO+k,cv0,cvN,V2[i*NMO+k]);
      else
        add_to_vn(vn,map_,cut,i*NMO+k,k*NMO+i,cv0,cvN,V2[i*NMO+k],V2[k*NMO+i]);
    }
  // if distributed over equivalent nodes, perform communication step here

  TG.node_barrier();
}

}

namespace HamHelperSymmetric
{

namespace
{

void count_over_cholvec(double cut, std::vector<std::size_t>& count, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lik)
{
  assert(c1>=c0);
  if(c0==c1) return;
  auto cik = std::lower_bound( std::addressof(*Lik.non_zero_indices2_data()),
                               std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                               c0);
  auto cik_end = std::lower_bound( cik,
                                   std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                                   c1);
  auto vik = std::addressof(*Lik.non_zero_values_data()) +
                    std::distance(std::addressof(*Lik.non_zero_indices2_data()),cik);

  using std::conj;
  using std::abs;
  using std::size_t;
  while( cik!=cik_end ) {
    if(abs(*vik) > cut) count[*cik]++;
    ++cik;
    ++vik;
  }
}

void count_nnz(double cut, size_t& ni, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lik)
{
  using std::conj;
  using std::abs;
  using std::size_t;
  assert(c1>=c0);
  ni=size_t(0);
  if(c0==c1) return;
  auto ci = std::lower_bound( std::addressof(*Lik.non_zero_indices2_data()),
                              std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                              c0);
  auto ci_end = std::lower_bound( ci,
                                  std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                                  c1);
  auto vi = std::addressof(*Lik.non_zero_values_data()) +
                    std::distance(std::addressof(*Lik.non_zero_indices2_data()),ci);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  while( ci!=ci_end ) {
    if(abs(*vi) > cut) ++ni;
    ++ci;
    ++vi;
  }
}

void add_to_vn(SpVType_shm_csr_matrix& vn, double cut, int ik, int c0, int c1, SpVType_shm_csr_matrix::reference const& Lik)
{
  using std::conj;
  using std::abs;
  using std::size_t;
  assert(c1>=c0);
  if(c0==c1) return;
  auto ci = std::lower_bound( std::addressof(*Lik.non_zero_indices2_data()),
                              std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                              c0);
  auto ci_end = std::lower_bound( ci,
                                  std::addressof(*(Lik.non_zero_indices2_data()
                                               +Lik.num_non_zero_elements())),
                                  c1);
  auto vi = std::addressof(*Lik.non_zero_values_data()) +
                    std::distance(std::addressof(*Lik.non_zero_indices2_data()),ci);

  // ignoring factor of 0.5 to keep it consistent with the old code for now
  int c_origin = vn.global_origin()[1];
  while( ci!=ci_end ) {
    if(abs(*vi) > cut) {
      assert((*ci)-c_origin < vn.shape()[1]);
      vn.emplace_back({ik,(*ci)-c_origin},static_cast<SPValueType>(*vi));
    }
    ++ci;
    ++vi;
  }
}

}

std::vector<std::size_t> count_nnz_per_cholvec(double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO)
{
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: count_nnz_per_cholvec is not designed for distributed CholMat. \n");

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if(V2.shape()[0] != NMO*NMO)
    APP_ABORT(" Error in count_nnz_per_cholvec: V2.shape()[0] ! NMO*NMO \n");
  int ik0, ikN;
  std::tie(ik0, ikN) = FairDivideBoundary(TG.getGlobalRank(),int(NMO*NMO),TG.getGlobalSize());

  if(cut < 1e-12) cut=1e-12;
  std::vector<std::size_t> counts(V2.shape()[1]);
  int nvecs = V2.shape()[1];

  for(int i=0, cnt=0; i<NMO; i++)
    for(int k=0; k<NMO; k++, cnt++) {
      if(cnt < ik0) continue;
      if(cnt == ikN) break;
      count_over_cholvec(cut,counts,0,nvecs,V2[cnt]);
    }
  TG.Global().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());

  return counts;
}

std::vector<std::size_t> count_nnz_per_ik(double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO, int cv0, int cvN)
{
  assert(cv0 >= 0 && cvN <= V2.shape()[1]);
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: count_nnz_per_ik is not designed for distributed CholMat. \n");

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if(V2.shape()[0] != NMO*NMO)
    APP_ABORT(" Error in count_nnz_per_cholvec: V2.shape()[0] ! NMO*NMO \n");
  int ik0, ikN;
  std::tie(ik0, ikN) = FairDivideBoundary(coreid,int(NMO*NMO),ncores);

  if(cut < 1e-12) cut=1e-12;
  std::vector<std::size_t> counts(V2.shape()[0]);

  std::size_t nik;
  for(int i=0, cnt=0; i<NMO; i++)
    for(int k=0; k<NMO; k++, cnt++) {
      if(cnt < ik0) continue;
      if(cnt == ikN) break;
      count_nnz(cut,counts[cnt],cv0,cvN,V2[cnt]);
    }
  TG.Node().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());
  // if divided over equivalent nodes, reduce over the Core_Equivalent communicator
  //TG.EqvCore().all_reduce_in_place_n(counts.begin(),counts.size(),std::plus<>());

  return counts;
}

void generateHSPotential(SpVType_shm_csr_matrix& vn, double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO, int cv0, int cvN)
{

  assert(cv0 >= 0 && cvN <= V2.shape()[1]);
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: generateHSPotential is not designed for distributed CholMat. \n");

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  if(V2.shape()[0] != NMO*NMO)
    APP_ABORT(" Error in generateHSPotential: V2.shape()[0] ! NMO*NMO \n");
  int ik0, ikN;
  std::tie(ik0, ikN) = FairDivideBoundary(coreid,int(NMO*NMO),ncores);

  if(cut < 1e-12) cut=1e-12;
  for(int i=0, cnt=0; i<NMO; i++)
    for(int k=0; k<NMO; k++, cnt++) {
      if(cnt < ik0) continue;
      if(cnt == ikN) break;
      add_to_vn(vn,cut,cnt,cv0,cvN,V2[cnt]);
    }
  // if distributed over equivalent nodes, perform communication step here
  TG.node_barrier();
}

}
/*

  bool FactorizedSparseHamiltonian_old::createHamiltonianForPureDeterminant(int walker_type, bool aa_only, std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b , std::vector<s1D<ValueType> >& hij, SPValueSMSpMat& Vijkl, const RealType cut)
  {
    if(walker_type==2) {
      APP_ABORT("Error: GHF density matrix only implemented with spinRestricted integrals. \n");
    }

    if(distribute_Ham) {
      APP_ABORT("Error: createHamiltonianForPureDeterminant not implemented with distributed Hamiltonian. \n");
    }

    // teporary until mpi3 is fully integrated
    TaskGroup_ TGham(TG,"DummyHS",1,TG.getNCoresPerTG());

    createHij(walker_type,NMO,occ_a,occ_b,H1,hij,cut);

    Timer.reset("Generic");
    Timer.start("Generic");

    ValueType zero = ValueType(0);

    std::vector<s4D<ValueType> > vs4D;
    vs4D.reserve(24);

    std::size_t cnt2=0;
    long cnter=0, npr = TG.getGlobalSize(), rk = TG.getGlobalRank();
    OrbitalType i,j,k,l,j1,k1,l1,j2,k2,l2;
    ValueType J1,J2,J3,J1a=zero,J2a=zero,J3a=zero,fct;

    DiagHam.resize(NMO,NMO);
    for(IndexType i=0; i<NMO; i++)
    for(IndexType k=i; k<NMO; k++, cnter++) {
      if( cnter%npr != rk ) continue;
      DiagHam(i,k) =  H(i,k,k,i);
      DiagHam(k,i) = DiagHam(i,k);
#if defined(QMC_COMPLEX)
      if(DiagHam(i,k).imag() > 1e-8) {
          app_error()<<" Error: Found complex diagonal on hamiltonian. " <<i <<" " <<k <<" " <<DiagHam(i,k) <<std::endl;
          APP_ABORT("Error: Found complex diagonal on hamiltonian.");
      }
      ComplexType Hik = H(i,k,i,k);
      if(Hik.imag() > cutoff_cholesky) {
          app_error()<<" Error: Problems with factorized Hamiltonian <ik|ik> should be real: "
                     <<i <<" " <<k <<" " <<Hik <<std::endl;
          APP_ABORT("Error: Problems with factorized hamiltonian.");
      }
#endif
    }
    {
      ValueMatrix dum(DiagHam);
      TG.Global().all_reduce(dum.begin(),dum.end(),DiagHam.begin());
    }

    int occi, occj, occk, occl;
    cnter=0;
    // Approximation: (similar to non factorized case)
    //   - if <ij|kl>  <  cutoff1bar, set to zero
    for(IndexType i=0; i<NMO; i++)  {
    for(IndexType j=i; j<NMO; j++,cnter++)  {
     if( cnter%npr != rk ) continue;
     occi = (occ_a[i]||occ_b[i])?1:0;
     occj = (occ_a[j]||occ_b[j])?1:0;
    for(IndexType k=j; k<NMO; k++)  {
     occk = (occ_a[k]||occ_b[k])?1:0;
    for(IndexType l=k; l<NMO; l++)  {

        occl = (occ_a[l]||occ_b[l])?1:0;
        if( occi + occj + occk + occl < 2 ) continue;

        // NOTE NOTE NOTE:
        // <ik|ik> can get a complex component due to truncation error, eliminate it here

        J1=J2=J3=zero;
        // J1 = <ij|kl>
        // J1 < sqrt( <ik|ki> * <lj|jl> )
        if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(l,j)) ) > cutoff1bar ) {
          J1 = H(i,j,k,l);
#if defined(QMC_COMPLEX)
          if(i==k && j==l) J1 = ValueType(J1.real(),0);
#endif
        }

        // J2 = <ij|lk>
        if(i==j || l==k) {
          J2=J1;
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(k,j)) ) > cutoff1bar ) {
            J2 = H(i,j,l,k);
#if defined(QMC_COMPLEX)
            if(i==l && j==k) J2 = ValueType(J2.real(),0);
#endif
          }
        }

        // J3 = <ik|jl>
        if(j==k) {
          J3=J1;
        } else if(i==l) {
          J3 = myconj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(l,k)) ) > cutoff1bar ) {
            J3 = H(i,k,j,l);
#if defined(QMC_COMPLEX)
            if(i==j && k==l) J3 = ValueType(J3.real(),0);
#endif
          }
        }

#if defined(QMC_COMPLEX)
        J1a=J2a=J3a=zero;

        //  J2a = <ik|lj>
        if(l==j) {
          J2a=J3;
        } else if(i==k) {
          J2a=J3;
        } else if(k==j) {
          J2a=J2;
        } else if(i==l) {
          J2a=std::conj(J2);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(j,k)) ) > cutoff1bar )
            J2a = H(i,k,l,j);
        }

        //  J3a = <il|jk>
        if(l==j) {
          J3a=J2;
        } else if(i==k) {
          J3a=std::conj(J2);
        } else if(k==l) {
          J3a=J3;
        } else if(i==j) {
          J3a=std::conj(J3);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(k,l)) ) > cutoff1bar )
            J3a = H(i,l,j,k);
        }

        //  For complex, there are 3 extra non-symmetric terms:
        //  J1a = <il|kj>
        if(k==l) {
          J1a=J2a;
        } else if(i==j) {
          J1a=std::conj(J2a);
        } else if(j==k) {
          J1a=J3a;
        } else if(i==l) {
          J1a=J3a;
        } else if(l==j) {
          J1a=J1;
        } else if(i==k) {
          J1a=std::conj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(j,l)) ) > cutoff1bar )
            J1a = H(i,l,k,j);
        }
#endif



        if( std::abs(J1)<cutoff1bar && std::abs(J2)<cutoff1bar && std::abs(J3)<cutoff1bar && std::abs(J1a)<cutoff1bar && std::abs(J2a)<cutoff1bar && std::abs(J3a)<cutoff1bar) continue;

        vs4D.clear();
        if(walker_type==0) {
          find_all_contributions_to_hamiltonian_closed_shell(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==1) {
          find_all_contributions_to_hamiltonian_collinear(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==2) {
          find_all_contributions_to_hamiltonian_ghf(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else {
          APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
        }
        cnt2+=count_allowed_terms(vs4D,occ_a,occ_b);

    }
    }
    }
    }

    cnt2 = TG.Global().all_reduce_value(cnt2);
    //number_of_terms = cnt2;
    app_log()<<" Number of terms in hamiltonian: " <<cnt2 <<std::endl;
    Vijkl.reserve(cnt2);
    cnt2=0;

    cnter=0;
    // Approximation: (similar to non factorized case)
    //   - if <ij|kl>  <  cutoff1bar, set to zero
    for(IndexType i=0; i<NMO; i++)  {
    for(IndexType j=i; j<NMO; j++,cnter++)  {
     if( cnter%npr != rk ) continue;
     occi = (occ_a[i]||occ_b[i])?1:0;
     occj = (occ_a[j]||occ_b[j])?1:0;
    for(IndexType k=j; k<NMO; k++)  {
     occk = (occ_a[k]||occ_b[k])?1:0;
    for(IndexType l=k; l<NMO; l++)  {

        occl = (occ_a[l]||occ_b[l])?1:0;
        if( occi + occj + occk + occl < 2 ) continue;

        J1=J2=J3=zero;
        // J1 = <ij|kl>
        // J1 < sqrt( <ik|ki> * <lj|jl> )
        if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(l,j)) ) > cutoff1bar ) {
          J1 = H(i,j,k,l);
#if defined(QMC_COMPLEX)
          if(i==k && j==l) J1 = ValueType(J1.real(),0);
#endif
        }

        // J2 = <ij|lk>
        if(i==j || l==k) {
          J2=J1;
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(k,j)) ) > cutoff1bar ) {
            J2 = H(i,j,l,k);
#if defined(QMC_COMPLEX)
            if(i==l && j==k) J2 = ValueType(J2.real(),0);
#endif
          }
        }

        // J3 = <ik|jl>
        if(j==k) {
          J3=J1;
        } else if(i==l) {
          J3 = myconj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(l,k)) ) > cutoff1bar ) {
            J3 = H(i,k,j,l);
#if defined(QMC_COMPLEX)
            if(i==j && k==l) J3 = ValueType(J3.real(),0);
#endif
          }
        }

#if defined(QMC_COMPLEX)
        J1a=J2a=J3a=zero;

        //  J2a = <ik|lj>
        if(l==j) {
          J2a=J3;
        } else if(i==k) {
          J2a=J3;
        } else if(k==j) {
          J2a=J2;
        } else if(i==l) {
          J2a=std::conj(J2);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(j,k)) ) > cutoff1bar )
            J2a = H(i,k,l,j);
        }

        //  J3a = <il|jk>
        if(l==j) {
          J3a=J2;
        } else if(i==k) {
          J3a=std::conj(J2);
        } else if(k==l) {
          J3a=J3;
        } else if(i==j) {
          J3a=std::conj(J3);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(k,l)) ) > cutoff1bar )
            J3a = H(i,l,j,k);
        }

        //  For complex, there are 3 extra non-symmetric terms:
        //  J1a = <il|kj>
        if(k==l) {
          J1a=J2a;
        } else if(i==j) {
          J1a=std::conj(J2a);
        } else if(j==k) {
          J1a=J3a;
        } else if(i==l) {
          J1a=J3a;
        } else if(l==j) {
          J1a=J1;
        } else if(i==k) {
          J1a=std::conj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(j,l)) ) > cutoff1bar )
            J1a = H(i,l,k,j);
        }

#endif

        if( std::abs(J1)<cutoff1bar && std::abs(J2)<cutoff1bar && std::abs(J3)<cutoff1bar && std::abs(J1a)<cutoff1bar && std::abs(J2a)<cutoff1bar && std::abs(J3a)<cutoff1bar) continue;

        vs4D.clear();
        if(walker_type==0) {
          find_all_contributions_to_hamiltonian_closed_shell(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==1) {
          find_all_contributions_to_hamiltonian_collinear(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==2) {
          find_all_contributions_to_hamiltonian_ghf(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else {
          APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
        }
        cnt2+=add_allowed_terms(NMO,vs4D,occ_a,occ_b, Vijkl, true, walker_type==2);
    }
    }
    }
    }
    TG.global_barrier();
    Timer.stop("Generic");
    app_log()<<"Time to generate full Hamiltonian from factorized form: " <<Timer.total("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");
    redistribute_sparse_matrix(TGham,Vijkl);
    Timer.stop("Generic");
    app_log()<<"Time to communicate Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");

    if(!Vijkl.remove_repeated_and_compress(TG.Node().impl_)) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }

    Timer.stop("Generic");
    app_log()<<"Time to remove_repeated_and_compress Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    return true;

  }
*/


}
}

