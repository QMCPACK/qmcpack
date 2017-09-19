
#include<AFQMC/Utilities/Utils.h>
#include <numeric>
#include<iostream>


namespace std{ 

  void swap(std::tuple<int &, int &, qmcplusplus::RealType &> const& a, std::tuple<int &, int &, qmcplusplus::RealType &> const& b) {
    using std::swap;
    swap(std::get<0>(a), std::get<0>(b));
    swap(std::get<1>(a), std::get<1>(b));
    swap(std::get<2>(a), std::get<2>(b));
  } 

  void swap(std::tuple<int &, int &, std::complex<qmcplusplus::RealType> &> const & a, std::tuple<int &, int &, std::complex<qmcplusplus::RealType> &> const& b) {
    using std::swap;
    swap(std::get<0>(a), std::get<0>(b));
    swap(std::get<1>(a), std::get<1>(b));
    swap(std::get<2>(a), std::get<2>(b));
  } 
 
}

namespace qmcplusplus { 

// put compiler guards for serial execution
/*
template<class T, class Compare>
void parallel_inplace_merge(int np, int rk, T beg, T mid, T end, MPI_Comm comm, Compare comp)
{

  if(np==1) {
    std::inplace_merge(beg,mid,end,comp);
    return;
  }

  MPI_Barrier(comm);

  T p1, p2;
  if( std::distance(beg,mid) >= std::distance(mid,end) ) {
    p1 = beg + std::distance(beg,mid)/2;
    auto it = std::lower_bound(mid,end,*p1,comp);
    p2 = it;
  } else {
    p2 = mid + std::distance(mid,end)/2;
    auto it = std::lower_bound(beg,mid,*p2,comp);
    p1 = it;
  }

  MPI_Barrier(comm);
  if(rk==0) std::rotate(p1,mid,p2);
  MPI_Barrier(comm);

  mid = p1 + std::distance(mid,p2);

  if(rk < np/2)
    parallel_inplace_merge(np/2,rk,beg,p1,mid,comm,comp);
  else
    parallel_inplace_merge(np/2,rk-np/2,mid,p2,end,comm,comp);
}
*/

// given a list of (N+1) integers, this routine attempts to find a partitioning of n continuous subsets  
// such that the sum of elements in each set is approximately homogeneous
// In other words, the routine will minimize the variance of the difference between the sums in each set
// The number of elements in bucket i are given by indx[i+1]-indx[i]. In other words, tasks from indx[i] through indx[i+1]
// are assigned to bucket i. There are N buckets 
template<typename IType>
void balance_partition_ordered_set(int N, IType* indx, std::vector<IType>& subsets) 
{
    int64_t avg=0;

    // finds optimal position for subsets[i] 
    auto step = [&] (int i) {
      IType i0 = subsets[i];
      subsets[i] = subsets[i-1]+1;
      int64_t vmin = std::abs(static_cast<int64_t>(*(indx+subsets[i]))
                            - static_cast<int64_t>(*(indx+subsets[i-1]))
                            - avg)
                   + std::abs(static_cast<int64_t>(*(indx+subsets[i+1]))
                            - static_cast<int64_t>(*(indx+subsets[i]))
                            - avg);
      for(int k=subsets[i-1]+2 ; k<subsets[i+1]; k++) {
        int64_t v = std::abs(static_cast<int64_t>(*(indx+k))
                           - static_cast<int64_t>(*(indx+subsets[i-1]))
                           - avg)
                  + std::abs(static_cast<int64_t>(*(indx+subsets[i+1]))
                           - static_cast<int64_t>(*(indx+k))
                           - avg);
        if( v < vmin ) {
          vmin=v;
          subsets[i] = k;
        }
      }
      return subsets[i]!=i0;
    };

    if(*(indx+N) == 0)
      APP_ABORT("Error in balance_partition_ordered_set(): empty hamiltonian. \n");

    IType nsets = subsets.size()-1;
    IType i0=0;
    IType iN = N;
    while( *(indx + i0) == *(indx + i0 + 1) ) i0++;
    while( *(indx + iN - 1) == *(indx + iN) ) iN--;
    int64_t avNpc = (iN-i0)/nsets;
    int64_t extra = (iN-i0)%nsets;
    for(IType i=0; i<nsets; i++)
      subsets[i]=( i<extra )?(i0+i*(avNpc+1)):(i0+i*avNpc+extra);
    subsets[nsets]=iN;

    for(IType i=0; i<nsets; i++)
      avg += static_cast<int64_t>(*(indx+subsets[i+1])) - static_cast<int64_t>(*(indx+subsets[i]));
    avg /= nsets;
    bool changed;
    do {
      changed=false;
      for(IType i=1; i<nsets; i++)
        changed |= step(i);
    } while( changed );

}

// careful 
// FIX FIX FIX
// this routine returns interchanged (i,j)/(k,l), so it is wrong due to a std::complex conjugation for std::complex matrix elements 
int cntExcitations(int NAEA, int NAEB, std::vector<IndexType>& DL, std::vector<IndexType>& DR, IndexType& n0, IndexType& n1, IndexType& n2, IndexType& n3, std::vector<IndexType>& occ, RealType& sg)
{
  std::vector<IndexType>::iterator itR = DR.begin();
  std::vector<IndexType>::iterator itL = DL.begin();
  sg = 0.0;
  int cnt=0,pos=0,ind[20],cnt2=0,nq=0,cnt3=0;
  bool found;
  int dummy = 1000000;
  n0=n1=n2=n3=dummy;

  for(int i=0; i<NAEA; i++) {
   found=false;
   for(int j=0; j<NAEA; j++)
     if(*(itL+i) == *(itR+j)) {
       found = true;
       occ[cnt2++] = *(itL+i);
       *(itL+i) = dummy;
       *(itR+j) = dummy;
       break;
     }
   if(!found) {
     if(cnt<2) ind[cnt]=i;
     cnt++;
     if(cnt > 2) {
       sg=0.0;
       return 2*cnt;
     }
   }
  }
  for(int i=NAEA; i<NAEA+NAEB; i++) {
   found=false;
   for(int j=NAEA; j<NAEA+NAEB; j++)
     if(*(itL+i) == *(itR+j)) {
       found = true;
       occ[cnt2++] = *(itL+i);
       *(itL+i) = dummy;
       *(itR+j) = dummy;
       break;
     }
   if(!found) {
     if(cnt<2) ind[cnt]=i;
     cnt++;
     if(cnt > 2) {
       sg=0.0;
       return 2*cnt;
     }
   }
  }
  if(cnt == 1) {
    n1=static_cast<IndexType>(*(itL+ind[0]));
    for(int i=0; i<NAEA+NAEB; i++)
     if(*(itR+i) != dummy) {   // there should be only one
       nq = ind[0]-i;
       n0=static_cast<IndexType>(*(itR+i));
       break;
     }
    sg = nq%2==0?1.0:-1.0;
  } else if(cnt == 2)  {
    int iq1=-1,iq2=-1;
    n2=static_cast<IndexType>(*(itL+ind[0]));
    n3=static_cast<IndexType>(*(itL+ind[1]));
    for(int i=0; i<NAEA+NAEB; i++)
     if(*(itR+i) != dummy) {  
       n0=static_cast<IndexType>(*(itR+i));
       iq1=i;
       break;
     }
    for(int i=iq1+1; i<NAEA+NAEB; i++)
     if(*(itR+i) != dummy) {   // there should be only one
       n1=static_cast<IndexType>(*(itR+i));
       iq2=i;
       break;
     }
    if(iq1<0 || iq2<0)
      APP_ABORT("Error in: cntExcitations.\n");
    nq = ind[0]-iq1+ind[1]-iq2;
    sg = nq%2==0?1.0:-1.0;
  } else
    sg=0.0;
  return 2*cnt;
}

template
void balance_partition_ordered_set(int N, uint32_t* indx, std::vector<uint32_t>& subsets); 
template
void balance_partition_ordered_set(int N, int* indx, std::vector<int>& subsets); 
template
void balance_partition_ordered_set(int N, long* indx, std::vector<long>& subsets); 
template
void balance_partition_ordered_set(int N, std::size_t* indx, std::vector<std::size_t>& subsets); 

}

