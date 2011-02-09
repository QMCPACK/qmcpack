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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/WalkerControlMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/NewTimer.h"
#include "Utilities/Timer.h"
using namespace qmcplusplus;

//#if defined(PRINT_DEBUG)
//#define DMC_BRANCH_START(NOW) NOW
//#define DMC_BRANCH_STOP(TID,TM) TID=TM
//#define DMC_BRANCH_DUMP(IT,T1,T2,T3)  \
//  OhmmsInfo::Debug->getStream()  << "BRANCH " \
//<< setw(8) << IT \
//<< " SORT" << setw(16) << T1 \
//<< " GSUM" << setw(16) << T2 \
//<< " SWAP" << setw(16) << T3 << std::endl
//#else
#define DMC_BRANCH_START(NOW) 
#define DMC_BRANCH_STOP(TID,TM) 
#define DMC_BRANCH_DUMP(IT,T1,T2,T3) 
//#endif


//#define MCWALKERSET_MPI_DEBUG

/** default constructor
 *
 * set SwapMode
 */
WalkerControlMPI::WalkerControlMPI(Communicate* c): WalkerControlBase(c) 
{
  SwapMode=1;
  Cur_min=0;
  Cur_max=0; 
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  ofstream fout(fname);
#endif

  myTimers.push_back(new NewTimer("WalkerControlMPI::branch")); //timer for the branch
  myTimers.push_back(new NewTimer("WalkerControlMPI::pre-loadbalance")); //timer for the branch
  myTimers.push_back(new NewTimer("WalkerControlMPI::loadbalance")); //timer for the branch
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[2]);
}

int 
WalkerControlMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger) 
{

  DMC_BRANCH_START(Timer localTimer);
  TinyVector<RealType,3> bTime(0.0);

  myTimers[0]->start();

  myTimers[1]->start();
  std::fill(curData.begin(),curData.end(),0);
  //std::fill(NumPerNode.begin(),NumPerNode.end(),0);
  sortWalkers(W);

  //use NumWalkersSent from the previous exchange
  curData[SENTWALKERS_INDEX]=NumWalkersSent;
  //update the number of walkers for this node
  curData[LE_MAX+MyContext]=NumWalkers;

  DMC_BRANCH_STOP(bTime[0],localTimer.elapsed());
  DMC_BRANCH_START(localTimer.restart());

  int nw = copyWalkers(W);
  myComm->allreduce(curData);
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;

  DMC_BRANCH_STOP(bTime[1],localTimer.elapsed());
  DMC_BRANCH_START(localTimer.restart());

  ////update the samples and weights
  //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  //RealType wgtInv(1.0/curData[WEIGHT_INDEX]);
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];

  Cur_pop=0;
  for(int i=0, j=LE_MAX; i<NumContexts; i++,j++) {
    Cur_pop+= NumPerNode[i]=static_cast<int>(curData[j]);
  }
  myTimers[1]->stop();

  myTimers[2]->start();
  swapWalkersSimple(W); 
  myTimers[2]->stop();

  //Do not need to use a trigger.
  //Cur_min=Nmax; 
  //Cur_max=0; 
  //Cur_pop=0;
  //for(int i=0, j=LE_MAX; i<NumContexts; i++,j++) {
  //  Cur_pop+= NumPerNode[i]=static_cast<int>(curData[j]);
  //  Cur_min = std::min(Cur_min,NumPerNode[i]);
  //  Cur_max = std::max(Cur_max,NumPerNode[i]);
  //}
  //int max_diff = std::max(Cur_max*NumContexts-Cur_pop,Cur_pop-Cur_min*NumContexts);
  //double diff_pop = static_cast<double>(max_diff)/static_cast<double>(Cur_pop);
  //if(diff_pop > trigger) { 
  //  swapWalkersSimple(W); 
  //  //swapWalkersMap(W); 
  //}

  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end) {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }

  //update the global number of walkers and offsets
  W.setGlobalNumWalkers(Cur_pop);
  W.setWalkerOffsets(FairOffSet);

  DMC_BRANCH_STOP(bTime[2],localTimer.elapsed());

  DMC_BRANCH_DUMP(iter,bTime[0],bTime[1],bTime[2]);

  myTimers[0]->stop();
  return Cur_pop;
}

/** swap Walkers with Recv/Send 
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional. 
 */
void WalkerControlMPI::swapWalkersSimple(MCWalkerConfiguration& W) {

  FairDivideLow(Cur_pop,NumContexts,FairOffSet);
  vector<int> minus, plus;
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++) {
    int dn=NumPerNode[ip]-(FairOffSet[ip+1]-FairOffSet[ip]);
    if(ip == MyContext) deltaN=dn;
    if(dn>0) {
      plus.insert(plus.end(),dn,ip); 
    } else if(dn<0){
      minus.insert(minus.end(),-dn,ip); 
    }
  }

  Walker_t& wRef(*W[0]);
  vector<Walker_t*> newW;
  vector<Walker_t*> oldW; 

#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  ofstream fout(fname, ios::app);
  //fout << NumSwaps << " " << Cur_pop << " ";
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerNode[ic] << " ";
  //fout << " | ";
  //for(int ic=0; ic<NumContexts; ic++) fout << FairOffSet[ic+1]-FairOffSet[ic] << " ";
  //fout << " | ";
  for(int ic=0; ic<plus.size(); ic++) {
    fout << plus[ic] << " ";
  }
  fout << " | ";
  for(int ic=0; ic<minus.size(); ic++) {
    fout << minus[ic] << " ";
  }
  fout << endl;
#endif

  int nswap=std::min(plus.size(), minus.size());
  int last=W.getActiveWalkers()-1;
  int nsend=0;
  for(int ic=0; ic<nswap; ic++) {
    if(plus[ic]==MyContext) {
      //OOMPI_Packed sendBuffer(wRef.byteSize(),OOMPI_COMM_WORLD);
      OOMPI_Packed sendBuffer(wRef.byteSize(),myComm->getComm());
      W[last]->putMessage(sendBuffer);
      //OOMPI_COMM_WORLD[minus[ic]].Send(sendBuffer);
      myComm->getComm()[minus[ic]].Send(sendBuffer);
      --last; ++nsend;
    }
    if(minus[ic]==MyContext) {
      //OOMPI_Packed recvBuffer(wRef.byteSize(),OOMPI_COMM_WORLD);
      OOMPI_Packed recvBuffer(wRef.byteSize(),myComm->getComm());
      //OOMPI_COMM_WORLD[plus[ic]].Recv(recvBuffer);
      myComm->getComm()[plus[ic]].Recv(recvBuffer);
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(recvBuffer);
      newW.push_back(awalker);
    }
  }

  //save the number of walkers sent
  NumWalkersSent=nsend;

  if(nsend) {
    nsend=NumPerNode[MyContext]-nsend;
    W.destroyWalkers(W.begin()+nsend, W.end());
  }


  //add walkers from other node
  if(newW.size()) W.insert(W.end(),newW.begin(),newW.end());
}

/** swap Walkers with Irecv/Send 
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional. 
 */
void WalkerControlMPI::swapWalkersAsync(MCWalkerConfiguration& W) {

  OffSet[0]=0;
  for(int i=0; i<NumContexts; i++) OffSet[i+1]=OffSet[i]+NumPerNode[i];
  FairDivideLow(Cur_pop,NumContexts,FairOffSet);

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
    //recvLeftBuffer = new OOMPI_Packed(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    //recvLeftRequest = OOMPI_COMM_WORLD[MyContext-1].Irecv(*recvLeftBuffer,MPI_ANY_TAG);
    recvLeftBuffer = new OOMPI_Packed(dn*wRef.byteSize(),myComm->getComm());
    recvLeftRequest = myComm->getComm()[MyContext-1].Irecv(*recvLeftBuffer,MPI_ANY_TAG);
  } 
  if(toRight>0) { //recv from node+1
    //recvRightBuffer = new OOMPI_Packed(toRight*wRef.byteSize(),OOMPI_COMM_WORLD);
    //recvRightRequest = OOMPI_COMM_WORLD[MyContext+1].Irecv(*recvRightBuffer,MPI_ANY_TAG);
    recvRightBuffer = new OOMPI_Packed(toRight*wRef.byteSize(),myComm->getComm());
    recvRightRequest = myComm->getComm()[MyContext+1].Irecv(*recvRightBuffer,MPI_ANY_TAG);
  }
  
  //send the data
  if(toLeft>0) { //send to node-1
    int dn=toLeft;
    num_deleted+=dn;
    //OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),myComm->getComm());
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    //OOMPI_COMM_WORLD[MyContext-1].Send(sendBuffer);
    myComm->getComm()[MyContext-1].Send(sendBuffer);
  }
  if(toRight<0) { // send to node+1
    int dn=-toRight;
    num_deleted+=dn;
    //OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),myComm->getComm());
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    //OOMPI_COMM_WORLD[MyContext+1].Send(sendBuffer);
    myComm->getComm()[MyContext+1].Send(sendBuffer);
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
void WalkerControlMPI::swapWalkersBlocked(MCWalkerConfiguration& W) {

  OffSet[0]=0;
  for(int i=0; i<NumContexts; i++) OffSet[i+1]=OffSet[i]+NumPerNode[i];
  FairDivideLow(Cur_pop,NumContexts,FairOffSet);

  int toLeft=FairOffSet[MyContext]-OffSet[MyContext];
  int toRight=FairOffSet[MyContext+1]-OffSet[MyContext+1];

  Walker_t& wRef(*W[0]);

  vector<Walker_t*> newW;

  int num_deleted=0;
  int last=NumPerNode[MyContext]-1;
  //scehdule irecv
  if(toLeft<0) { // recv from node-1 
    int dn=-toLeft;
    //OOMPI_Packed recvBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    //OOMPI_COMM_WORLD[MyContext-1].Recv(recvBuffer);
    OOMPI_Packed recvBuffer(dn*wRef.byteSize(),myComm->getComm());
    myComm->getComm()[MyContext-1].Recv(recvBuffer);
    while(dn) {
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(recvBuffer);
      newW.push_back(awalker);
      --dn;
    }
  } else if(toLeft>0) { //send to node-1
    int dn=toLeft;
    num_deleted+=dn;
    //OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),myComm->getComm());
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    //OOMPI_COMM_WORLD[MyContext-1].Send(sendBuffer);
    myComm->getComm()[MyContext-1].Send(sendBuffer);
  }

  if(toRight<0) { // send to node+1
    int dn=-toRight;
    num_deleted+=dn;
    //OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),myComm->getComm());
    while(dn) {
      W[last]->putMessage(sendBuffer); --dn; --last;
    }
    //OOMPI_COMM_WORLD[MyContext+1].Send(sendBuffer);
    myComm->getComm()[MyContext+1].Send(sendBuffer);
  } else if(toRight>0) { //recv from node+1
    //OOMPI_Packed recvBuffer(toRight*wRef.byteSize(),OOMPI_COMM_WORLD);
    //OOMPI_COMM_WORLD[MyContext+1].Recv(recvBuffer);
    OOMPI_Packed recvBuffer(toRight*wRef.byteSize(),myComm->getComm());
    myComm->getComm()[MyContext+1].Recv(recvBuffer);
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
void WalkerControlMPI::swapWalkersMap(MCWalkerConfiguration& W) {

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
      //OOMPI_Packed recvBuffer(dnw*wRef.byteSize(),OOMPI_COMM_WORLD);
      //OOMPI_COMM_WORLD[target_context].Recv(recvBuffer);
      OOMPI_Packed recvBuffer(dnw*wRef.byteSize(),myComm->getComm());
      myComm->getComm()[target_context].Recv(recvBuffer);
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
      //OOMPI_Packed sendBuffer(dnw*wRef.byteSize(),OOMPI_COMM_WORLD);
      OOMPI_Packed sendBuffer(dnw*wRef.byteSize(),myComm->getComm());
      int last=W.getActiveWalkers()-1;
      while(dnw) {
        W[last]->putMessage(sendBuffer);
        --dnw; --last;
      }
      //OOMPI_COMM_WORLD[target_context].Send(sendBuffer);
      myComm->getComm()[target_context].Send(sendBuffer);
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
 * $RCSfile: WalkerControlMPI.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

