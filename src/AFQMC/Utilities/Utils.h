#ifndef QMCPLUSPLUS_AFQMC_UTILS_H
#define QMCPLUSPLUS_AFQMC_UTILS_H

#include<iostream>
#include<cstdlib>
#include<vector>
#include<string>
#include<algorithm>
#include<AFQMC/config.0.h>
#include<AFQMC/config.h>
#include<Message/Communicate.h>

namespace qmcplusplus { 

// given a list of N integers, this routine attempts to find a partitioning of n continuous subsets  
// such that the sum of elements in each set is approximately homogeneous
// In other words, the routine will minimize the variance of the difference between the sums in each set
// The number of elements in bucket i are given by indx[i+1]-indx[i]. In other words, tasks from indx[i] through indx[i+1]
// are assigned to bucket i. 
template<typename IType>
void balance_partition_ordered_set(int N, IType* indx, std::vector<IType>& subsets); 

// careful 
// FIX FIX FIX
// this routine returns interchanged (i,j)/(k,l), so it is wrong due to a std::complex conjugation for std::complex matrix elements 
int cntExcitations(int NAEA, int NAEB, std::vector<IndexType>& DL, std::vector<IndexType>& DR, IndexType& n0, IndexType& n1, IndexType& n2, IndexType& n3, std::vector<IndexType>& occ, RealType& sg);

// use concepts???
template<class Iter, class Compare>
inline void parallel_inplace_merge(int np, int rk, Iter* beg, Iter* mid, Iter* end, MPI_Comm comm, Compare comp)
{

  if(np==1) {
    std::inplace_merge(beg,mid,end,comp);
    return;
  }

  MPI_Barrier(comm);

  Iter *p1, *p2;
  if( std::distance(beg,mid) >= std::distance(mid,end) ) {
    p1 = beg + std::distance(beg,mid)/2;
    auto it = std::lower_bound(mid,end,*p1,comp);
    p2 = &(*it);
  } else {
    p2 = mid + std::distance(mid,end)/2;
    auto it = std::lower_bound(beg,mid,*p2,comp);
    p1 = &(*it);
  }

  MPI_Barrier(comm);
  if(rk==0) std::rotate(p1,mid,p2);
  MPI_Barrier(comm);

  mid = p1 + std::distance(mid,p2);

  if(rk < np/2)
    parallel_inplace_merge<Iter,Compare>(np/2,rk,beg,p1,mid,comm,comp);
  else
    parallel_inplace_merge<Iter,Compare>(np/2,rk-np/2,mid,p2,end,comm,comp);
}

}

namespace std{

  void swap(std::tuple<int &, int &, qmcplusplus::RealType &> const& a, std::tuple<int &, int &, qmcplusplus::RealType &> const& b); 
  void swap(std::tuple<int &, int &, std::complex<qmcplusplus::RealType> &> const& a, std::tuple<int &, int &, std::complex<qmcplusplus::RealType> &> const& b); 

}
#endif
