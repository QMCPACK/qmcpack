
#include<AFQMC/Utilities/Utils.h>
#include <numeric>
#include<iostream>


namespace std{ 

  void swap(std::tuple<int &, int &, qmcplusplus::ValueType &> const& a, std::tuple<int &, int &, qmcplusplus::ValueType &> const& b) {
    using std::swap;
    swap(std::get<0>(a), std::get<0>(b));
    swap(std::get<1>(a), std::get<1>(b));
    swap(std::get<2>(a), std::get<2>(b));
  } 

  void swap(std::tuple<int &, int &, std::complex<qmcplusplus::ValueType> &> const& a, std::tuple<int &, int &, std::complex<qmcplusplus::ValueType> &> const& b) {
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
void balance_partition_ordered_set(int N, int* indx, std::vector<int>& subsets) 
{

    int nsets = subsets.size()-1;
    std::vector<int64_t> sums(nsets);
    int64_t total,avg;
    int64_t var,old_var;
    auto get_var = [&] () {
      int64_t var_=0;
      for(std::vector<int64_t>::iterator it=sums.begin();
          it<sums.end(); it++)
        var_ += (*it-avg)*(*it-avg);
      return var_;
    };

    // attempts to move subsets[i] to the right or the left if it reduces var 
    auto single_step_boundary = [&] (int i) {

      int64_t dv1 = 0, dv2=0;
      if( subsets[i]-1==subsets[i-1] ) {
        // can't move left    
        if(subsets[i]+1==subsets[i+1]) return; // can't move right either 
        // try moving right
        int64_t sm1_;
        int64_t sm2_;
        sm1_ = *(indx+subsets[i]+1) - *(indx+subsets[i-1]);
        sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]+1);
        dv1 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
          - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
        while( dv1 <= 0 ) {
          var += dv1;
          subsets[i]++;
          sums[i-1] = sm1_;
          sums[i] = sm2_;
          sm1_ = *(indx+subsets[i]+1) - *(indx+subsets[i-1]);
          sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]+1);
          dv1 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
            - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
        }

      } else {
        if(subsets[i]+1==subsets[i+1]) {
          // can only move left

          int64_t sm1_;
          int64_t sm2_;
          sm1_ = *(indx+subsets[i]-1) - *(indx+subsets[i-1]);
          sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]-1);
          dv1 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
            - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
          while( dv1 <= 0 ) {
            var += dv1;
            subsets[i]--;
            sums[i-1] = sm1_;
            sums[i] = sm2_;
            sm1_ = *(indx+subsets[i]-1) - *(indx+subsets[i-1]);
            sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]-1);
            dv1 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
              - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
          }

        } else {
          // can move either way
          int64_t osm1 = sums[i-1], osm2 = sums[i], oset = subsets[i];
          int64_t dvtot1 = 0, dvtot2 = 0;

          // try moving left 
          int64_t sm1_;
          int64_t sm2_;
          sm1_ = *(indx+subsets[i]-1) - *(indx+subsets[i-1]);
          sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]-1);
          dv1 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
            - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
          while( dv1 <= 0 ) {
            dvtot1 += dv1;
            subsets[i]--;
            sums[i-1] = sm1_;
            sums[i] = sm2_;
            sm1_ = *(indx+subsets[i]-1) - *(indx+subsets[i-1]);
            sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]-1);
            dv1 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
              - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
          }
          //store
          int64_t lsm1 = sums[i-1], lsm2 = sums[i], lset = subsets[i];
          // restore
          sums[i-1]=osm1;
          sums[i]=osm2;
          subsets[i]=oset;

          // try moving right
          sm1_ = *(indx+subsets[i]+1) - *(indx+subsets[i-1]);
          sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]+1);
          dv2 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
            - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
          while( dv2 <= 0 ) {
            dvtot2 += dv2;
            subsets[i]++;
            sums[i-1] = sm1_;
            sums[i] = sm2_;
            sm1_ = *(indx+subsets[i]+1) - *(indx+subsets[i-1]);
            sm2_ = *(indx+subsets[i+1]) - *(indx+subsets[i]+1);
            dv2 = (sm1_-avg)*(sm1_-avg) + (sm2_-avg)*(sm2_-avg)
              - (sums[i]-avg)*(sums[i]-avg) - (sums[i-1]-avg)*(sums[i-1]-avg) ;
          }
          if(dvtot1 < dvtot2) {
            sums[i-1]=lsm1;
            sums[i]=lsm2;
            subsets[i]=lset;
            var += dvtot1;
          } else
            var += dvtot2;
          return;
        }
      }
    };

    if(*(indx+N) == 0)
      APP_ABORT("Error in PureSingleDeterminant::split_Ham_rows(): empty hamiltonian. \n");

    // stupid algorithm right now
    int i0=0;
    int iN = N;
    while( *(indx + i0) == *(indx + i0 + 1) ) i0++;
    while( *(indx + iN - 1) == *(indx + iN) ) iN--;
    int64_t avNpc = (iN-i0)/nsets;
    int64_t extra = (iN-i0)%nsets;
    for(int i=0; i<nsets; i++)
      subsets[i]=( i<extra )?(i0+i*(avNpc+1)):(i0+i*avNpc+extra);
    subsets[nsets]=iN;

    for(int i=0; i<nsets; i++)
      sums[i] = *(indx+subsets[i+1]) - *(indx+subsets[i]);

    total = std::accumulate(sums.begin(),sums.end(),0);
    avg = total/nsets;

    //var=get_var();
    do {
      //old_var = var;
      var=0;
      for(int i=1; i<nsets; i++) 
        single_step_boundary(i);

    } while( var < 0  );
    //} while( std::abs(old_var-var) > 0  );

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

}

