
#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cassert>
#include<complex>
#include<map>

#include "AFQMC/config.0.h"
#include "AFQMC/Wavefunctions/WavefunctionHelper.h" 

// Helper functions for slater determinant routines 

namespace qmcplusplus
{

// taking from FCIQMC code. I'm sure there's a nicer way to do this using STL
int cmpDets(int NAEA, int NAEB, int* n, double &sg, std::vector<IndexType>::iterator sdet1, std::vector<IndexType>::iterator sdet2, std::vector<IndexType>& work )
{

  sg=0.0;
  int cnt=0,pos=0,ind[20],cnt2=0,nq=0;
  bool found;
  IndexType dummy = 30000;
  for(int i=0; i<4; i++) n[i]=-1;

  work.resize(NAEA+NAEB); 
  std::copy(sdet2,sdet2+NAEA+NAEB,work.begin()); 

  for(int i=0; i<NAEA; i++) {
   found=false;
   for(int j=0; j<NAEA; j++)
     if(*(sdet1+i) == *(sdet2+j)) {
       found = true;
       work[j]=dummy;
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
     if(*(sdet1+i) == *(sdet2+j)) {
       found = true;
       work[j]=dummy;
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
    n[1]=static_cast<int>(  *(sdet1+ind[0]) );
    for(int i=0; i<NAEA+NAEB; i++) {
     if(work[i] != dummy) {   // there should be only one
       nq = ind[0]-i;
       n[0]=static_cast<int>(work[i]);
       break;
     }
    }
    sg = nq%2==0?1.0:-1.0;
  } else if(cnt == 2)  {
    int iq1=-1,iq2=-1;
    n[2]=static_cast<int>(  *(sdet1+ind[0]) );
    n[3]=static_cast<int>(  *(sdet1+ind[1]) );
    for(int i=0; i<NAEA+NAEB; i++)
     if(work[i] != dummy) {   // there should be only one
       n[0]=static_cast<int>(work[i]);
       iq1=i;
       break;
     }
    for(int i=iq1+1; i<NAEA+NAEB; i++)
     if(work[i] != dummy) {   // there should be only one
       n[1]=static_cast<int>(work[i]);
       iq2=i;
       break;
     }
    if(iq1 < 0 || iq2 < 0) {
std::cout<<"Problems in cmpDet: \n"
 <<"det1: ";
for(int i=0; i<NAEA+NAEB; i++) std::cout<<*(sdet1+i) <<" ";  
std::cout<<"\ndet2: ";
for(int i=0; i<NAEA+NAEB; i++) std::cout<<*(sdet2+i) <<" ";  
std::cout<<std::endl;
std::cout.flush();
    }
    assert(iq1>=0 && iq2>=0);
    nq = ind[0]-iq1+ind[1]-iq2;
    sg = nq%2==0?1.0:-1.0;
  } else
    sg=0.0;
  return 2*cnt;
}

}
