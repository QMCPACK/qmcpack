//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/GlobalWalkerControl.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
using namespace qmcplusplus;

/** partition ntot elements among npart
 * @param ntot total number of elements
 * @param npart number of partitions
 * @param adist distribution offset 
 *
 * adist[ip-1]-adist[ip] is the number of elements of ip partition
 * This method makes the zero-th node equal to or less than 1.
 */
template<class IV>
inline void FairPartition(int ntot, int npart, IV& adist) {
  int bat=ntot/npart;
  int residue = npart-ntot%npart;
  adist[0] = 0;
  for(int i=0; i<npart; i++) {
    if(i<residue)
      adist[i+1] = adist[i] + bat;
    else
      adist[i+1] = adist[i] + bat+1;
  }
}

/** default constructor
 *
 * set SwapMode
 */
GlobalWalkerControl::GlobalWalkerControl() {
  SwapMode=1;
  NumContexts=OHMMS::Controller->ncontexts();
  MyContext=OHMMS::Controller->mycontext();
  NumPerNode.resize(NumContexts);
  OffSet.resize(NumContexts+1);
  FairOffSet.resize(NumContexts+1);
  NumSwaps=0;

#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  ofstream fout(fname);
#endif
}

int 
GlobalWalkerControl::branch(int iter, MCWalkerConfiguration& W, RealType trigger) {

  std::fill(NumPerNode.begin(),NumPerNode.end(),0);

  sortWalkers(W);

  NumPerNode[MyContext] = NumWalkers;

  int nw = copyWalkers(W);

  //wait until everynode comes here
  OHMMS::Controller->barrier();
  gsum(NumPerNode,0);

  Cur_min=Nmax; Cur_max=0; Cur_pop=0;
  for(int i=0; i<NumContexts; i++) {
    Cur_pop+= NumPerNode[i];
    Cur_min = std::min(Cur_min,NumPerNode[i]);
    Cur_max = std::max(Cur_max,NumPerNode[i]);
  }

  int max_diff = std::max(Cur_max*NumContexts-Cur_pop,Cur_pop-Cur_min*NumContexts);
  double diff_pop = static_cast<double>(max_diff)/static_cast<double>(Cur_pop);

  if(diff_pop > trigger) { swapWalkersMap(W); }

  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end) {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }

  return Cur_pop;
}

/** swap Walkers with Irecv/Send 
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional. 
 */
void GlobalWalkerControl::swapWalkersAsync(MCWalkerConfiguration& W) {
  NumSwaps++;

  OffSet[0]=0;
  for(int i=0; i<NumContexts; i++) OffSet[i+1]=OffSet[i]+NumPerNode[i];
  FairPartition(Cur_pop,NumContexts,FairOffSet);

  int toLeft=FairOffSet[MyContext]-OffSet[MyContext];
  int toRight=FairOffSet[MyContext+1]-OffSet[MyContext+1];

  Walker_t& wRef(*W[0]);

  int num_deleted=0;
  int last=NumPerNode[MyContext]-1;

  OOMPI_Packed* recvLeftBuffer=0;
  OOMPI_Packed* recvRightBuffer=0;
  OOMPI_Request recvLeftRequest;
  OOMPI_Request recvRightRequest;

  //initiate Irecv
  if(toLeft<0) { // recv from node-1 
    int dn=-toLeft;
    recvLeftBuffer = new OOMPI_Packed(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    recvLeftRequest = OOMPI_COMM_WORLD[MyContext-1].Irecv(*recvLeftBuffer,MPI_ANY_TAG);
  } 
  if(toRight>0) { //recv from node+1
    recvRightBuffer = new OOMPI_Packed(toRight*wRef.byteSize(),OOMPI_COMM_WORLD);
    recvRightRequest = OOMPI_COMM_WORLD[MyContext+1].Irecv(*recvRightBuffer,MPI_ANY_TAG);
  }
  
  //send the data
  if(toLeft>0) { //send to node-1
    int dn=toLeft;
    num_deleted+=dn;
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    OOMPI_COMM_WORLD[MyContext-1].Send(sendBuffer);
  }
  if(toRight<0) { // send to node+1
    int dn=-toRight;
    num_deleted+=dn;
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    OOMPI_COMM_WORLD[MyContext+1].Send(sendBuffer);
  } 

  while(num_deleted>0) {
    W.pop_back(); --num_deleted;
  }

  if(recvLeftBuffer) {
    recvLeftRequest.Wait();
    int dn=-toLeft;
    while(dn) {
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(*recvLeftBuffer);
      W.push_back(awalker);
      --dn;
    }
    delete recvLeftBuffer;
    recvLeftBuffer=0;
  }

  if(recvRightBuffer) {
    recvRightRequest.Wait();
    int dn=toRight;
    while(dn) {
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(*recvRightBuffer);
      W.push_back(awalker);
      --dn;
    }
    delete recvRightBuffer;
    recvRightBuffer=0;
  }
}

/** swap Walkers with Recv/Send 
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional. 
 */
void GlobalWalkerControl::swapWalkersBlocked(MCWalkerConfiguration& W) {
  NumSwaps++;

  OffSet[0]=0;
  for(int i=0; i<NumContexts; i++) OffSet[i+1]=OffSet[i]+NumPerNode[i];
  FairPartition(Cur_pop,NumContexts,FairOffSet);

  int toLeft=FairOffSet[MyContext]-OffSet[MyContext];
  int toRight=FairOffSet[MyContext+1]-OffSet[MyContext+1];

  Walker_t& wRef(*W[0]);

  vector<Walker_t*> newW;

  int num_deleted=0;
  int last=NumPerNode[MyContext]-1;
  //scehdule irecv
  if(toLeft<0) { // recv from node-1 
    int dn=-toLeft;
    OOMPI_Packed recvBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_COMM_WORLD[MyContext-1].Recv(recvBuffer);
    while(dn) {
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(recvBuffer);
      newW.push_back(awalker);
      --dn;
    }
  } else if(toLeft>0) { //send to node-1
    int dn=toLeft;
    num_deleted+=dn;
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    OOMPI_COMM_WORLD[MyContext-1].Send(sendBuffer);
  }

  if(toRight<0) { // send to node+1
    int dn=-toRight;
    num_deleted+=dn;
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    OOMPI_COMM_WORLD[MyContext+1].Send(sendBuffer);
  } else if(toRight>0) { //recv from node+1
    OOMPI_Packed recvBuffer(toRight*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_COMM_WORLD[MyContext+1].Recv(recvBuffer);
    int dn=toRight;
    while(dn) {
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(recvBuffer);
      newW.push_back(awalker);
      --dn;
    }
  }

  while(num_deleted>0) {
    W.pop_back(); --num_deleted;
  }

  if(newW.size()) W.insert(W.end(),newW.begin(),newW.end());
}

/** swap walkers using (low,high) ordered pairs
 *
 * The communication occur only between the (low,high) pairs.
 * This does not guarantee a perfect balance swapWalkersBlocked and swapWalkersAsync
 * try to achieve. However, the number of messages and their sizes are less than
 * other methods.
 */
void GlobalWalkerControl::swapWalkersMap(MCWalkerConfiguration& W) {

  NumSwaps++;
  multimap<int,int> nw_map;
  for(int i=0; i<NumContexts; i++) {
    nw_map.insert(pair<int,int>(NumPerNode[i],i));
  }
  // multimap key is sorted with ascending order 
  multimap<int,int>::iterator it(nw_map.begin());
  multimap<int,int>::reverse_iterator it_b(nw_map.end());
  bool notpaired=true;
  int target_context=-1;
  int half=NumContexts/2;
  int item=0;
  bool minorcontext;
  while(notpaired &&item<half) {
    int i=(*it).second;
    int j=(*it_b).second;
    if(i == MyContext) {
      target_context=j;
      notpaired=false;
      minorcontext=true;
    } else if(j == MyContext) {
      target_context= i;
      notpaired=false;
      minorcontext=false;
    } 
    ++it; ++it_b; ++item;
  }

  int nw_tot=NumPerNode[MyContext]+NumPerNode[target_context];
  int nw_L=nw_tot/2;
  int nw_R=nw_tot-nw_L;
  int dnw(0);
  if(minorcontext) {
    dnw=nw_L-NumPerNode[MyContext];//how many to get
  } else {
    dnw=NumPerNode[MyContext]-nw_R;//how many to send
  }

  if(dnw) {//something to swap
    if(minorcontext) {//open recv buffer
      Walker_t& wRef(*W[0]);
      OOMPI_Packed recvBuffer(dnw*wRef.byteSize(),OOMPI_COMM_WORLD);
      OOMPI_COMM_WORLD[target_context].Recv(recvBuffer);
      //create walkers
      int last = W.getActiveWalkers();
      while(dnw) {
        Walker_t *awalker= new Walker_t(wRef);
        awalker->getMessage(recvBuffer);
        W.push_back(awalker);
        --dnw; ++last;
      }
    } else {
      Walker_t& wRef(*W[0]);
      OOMPI_Packed sendBuffer(dnw*wRef.byteSize(),OOMPI_COMM_WORLD);
      int last=W.getActiveWalkers()-1;
      while(dnw) {
        W[last]->putMessage(sendBuffer);
        --dnw; --last;
      }
      OOMPI_COMM_WORLD[target_context].Send(sendBuffer);
      //last=WalkerList.size()-1;
      //while(dnw_save) {
      //  delete WalkerList[last]; 
      //  WalkerList.pop_back();
      //  --dnw_save; --last;
      //}
      //destroyWalkers(WalkerList.begin()+nsub[MyContext], WalkerList.end());
      W.destroyWalkers(W.begin()+nw_R, W.end());
    }
  }

  /* not used yet
  struct lessNode {
    inline bool operator()(const pair<int,int>& a,
        const pair<int,int>& b) const {
      return a.second < b.second;
    }
  };
  //typedef pair<int,int> mytype;
  //vector<mytype> id2n(NumContexts);
  //for(int i=0; i<NumContexts; i++) {
  //  id2n=pair<int,int>(i,NumPerNode[i]);
  //}
  //
  //std::sort(id2n.begin(),id2n.end(),lessNode);
  */
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

