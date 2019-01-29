////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_KP_UTILITIES_HPP
#define AFQMC_KP_UTILITIES_HPP

#include<algorithm>

/*
 * Return true if PsiT is in block diagonal form, otherwise return false;
 * If it is block diagonal form, it also returns the number of states in each k-point block. 
 */
template<class Vector, class CSR, class Array>
bool get_nocc_per_kp(Vector const& nmo_per_kp, CSR const& PsiT, Array&& nocc_per_kp)
{
  int nkpts = nmo_per_kp.size();
  int N = PsiT.size(0);
  int M = PsiT.size(1);
  assert(nocc_per_kp.size() == nkpts);  

  std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
  std::vector<int> bounds(nkpts+1);
  bounds[0]=0;
  for(int k=0; k<nkpts; k++) 
    bounds[k+1] = bounds[k]+nmo_per_kp[k];
  int Q = 0;  
  for(int i=0; i<N; i++) {
    auto nt = PsiT.num_non_zero_elements(i);
    if(nt == 0) {
      std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
      return false;
    }
    auto col = PsiT.non_zero_indices2_data(i); 
    // check the kp index of the first non-zero column. Mut be either >= Q
    auto it = std::lower_bound(bounds.begin(),bounds.end(),*col+1)-1; 
    assert(it != bounds.end());
    int Q_ = std::distance(bounds.begin(),it);
    assert(Q_ >= 0 && Q_ < nkpts);
    if(Q_ < Q) {
      std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
      return false;
    }
    Q = Q_;
    for(int ni=0; ni<nt; ++ni, ++col) {
      if(*col < bounds[Q] || *col >= bounds[Q+1] ) {
        std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
        return false;
      }
    } 
    ++nocc_per_kp[Q];   
  }
  return true;
}

template<class Array, class Vector, class CSR>
Array get_PsiK(Vector const& nmo_per_kp, CSR const& PsiT, int K) 
{
  int nkpts = nmo_per_kp.size();
  int N = PsiT.size(0);
  int M = PsiT.size(1);

  int nel=0;
  std::vector<int> bounds(nkpts+1);
  bounds[0]=0;
  for(int k=0; k<nkpts; k++)
    bounds[k+1] = bounds[k]+nmo_per_kp[k];
  int Q = 0;
  for(int i=0; i<N; i++) {
    auto nt = PsiT.num_non_zero_elements(i);
    if(nt == 0) 
      APP_ABORT("Error: PsiT not in block-diagonal form in get_PsiK.\n");
    auto col = PsiT.non_zero_indices2_data(i);
    // check the kp index of the first non-zero column. Mut be either >= Q
    auto it = std::lower_bound(bounds.begin(),bounds.end(),*col+1)-1;
    assert(it != bounds.end());
    int Q_ = std::distance(bounds.begin(),it);
    assert(Q_ >= 0 && Q_ < nkpts);
    if(Q_ < Q) 
      APP_ABORT("Error: PsiT not in block-diagonal form in get_PsiK.\n");
    Q = Q_;
    for(int ni=0; ni<nt; ++ni, ++col) 
      if(*col < bounds[Q] || *col >= bounds[Q+1] ) 
        APP_ABORT("Error: PsiT not in block-diagonal form in get_PsiK.\n");
    if(Q==K) nel++;
  }
  using element = typename std::decay<Array>::type::element;  
  Array A({nel,nmo_per_kp[K]});
  using std::fill_n;
  fill_n(A.origin(),A.num_elements(),element(0));
  nel=0;
  for(int i=0; i<N; i++) {
    auto nt = PsiT.num_non_zero_elements(i);
    auto col = PsiT.non_zero_indices2_data(i);
    auto val = PsiT.non_zero_values_data(i);
    // check the kp index of the first non-zero column. Mut be either >= Q
    auto it = std::lower_bound(bounds.begin(),bounds.end(),*col+1)-1;
    int Q = std::distance(bounds.begin(),it);
    if(Q==K) {
      for(int ni=0; ni<nt; ++ni, ++col, ++val)
        A[nel][*col-bounds[K]] = static_cast<element>(*val);
      nel++;
    }
    if(Q > K) break;
  }  
  return A;
}



#endif
