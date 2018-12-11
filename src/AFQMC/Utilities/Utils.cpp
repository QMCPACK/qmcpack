
#include<AFQMC/Utilities/Utils.h>
#include <numeric>
#include <stack>
#include<iostream>
#include<complex>

namespace qmcplusplus { 

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

