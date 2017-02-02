#ifndef QMCPLUSPLUS_AFQMC_UTILS_H
#define QMCPLUSPLUS_AFQMC_UTILS_H

#include<iostream>
#include<cstdlib>
#include<vector>
#include<string>
#include<AFQMC/config.h>
#include<AFQMC/config.0.h>
#include<Message/Communicate.h>

namespace qmcplusplus { 

// given a list of N integers, this routine attempts to find a partitioning of n continuous subsets  
// such that the sum of elements in each set is approximately homogeneous
// In other words, the routine will minimize the variance of the difference between the sums in each set
// The number of elements in bucket i are given by indx[i+1]-indx[i]. In other words, tasks from indx[i] through indx[i+1]
// are assigned to bucket i. 
void balance_partition_ordered_set(int N, int* indx, std::vector<int>& subsets); 

// careful 
// FIX FIX FIX
// this routine returns interchanged (i,j)/(k,l), so it is wrong due to a std::complex conjugation for std::complex matrix elements 
int cntExcitations(int NAEA, int NAEB, std::vector<IndexType>& DL, std::vector<IndexType>& DR, IndexType& n0, IndexType& n1, IndexType& n2, IndexType& n3, std::vector<IndexType>& occ, RealType& sg);

}

namespace std{

void swap(std::tuple<int &, int &, qmcplusplus::ValueType &> const& a, std::tuple<int &, int &, qmcplusplus::ValueType &> const& b); 

}
#endif
