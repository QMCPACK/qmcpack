//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include<cassert>
#include <memory>
#include <mpi.h>
#include<AFQMC/config.0.h>
#include <Utilities/UtilityFunctions.h>

namespace qmcplusplus
{

namespace afqmc
{

/** swap Walkers with Recv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
template<class WlkBucket, 
         class IVec = std::vector<int>
         >
// eventually generalize MPI_Comm to a MPI wrapper
inline int swapWalkersSimple(WlkBucket& wlk, IVec& CurrNumPerNode, IVec& NewNumPerNode, MPI_Comm comm)
{
  int NumContexts, MyContext; 
  MPI_Comm_size(comm,&NumContexts);
  MPI_Comm_rank(comm,&MyContext);
  assert(CurrNumPerNode.size() >= NumContexts);  
  assert(NewNumPerNode.size() >= NumContexts);  
  auto Cur_pop = std::accumulate(CurrNumPerNode.begin(), CurrNumPerNode.begin()+NumContexts, 0);
  std::vector<int> minus, plus;
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    int dn=CurrNumPerNode[ip]-NewNumPerNode[ip];
    if(ip == MyContext)
      deltaN=dn;
    if(dn>0)
    {
      plus.insert(plus.end(),dn,ip);
    }
    else
      if(dn<0)
      {
        minus.insert(minus.end(),-dn,ip);
      }
  }
  int nswap=std::min(plus.size(), minus.size());
  int nsend=0;
  int wlk_size = wlk.single_walker_size();
  // 1 walker at a time
  std::vector<ComplexType> buff(wlk_size);  
  for(int ic=0; ic<nswap; ic++)
  {
    if(plus[ic]==MyContext)
    {
      wlk.pop_walkers(1,buff);
      MPI_Send(buff.data(),2*buff.size(),MPI_DOUBLE,minus[ic],plus[ic]+999,comm);
      ++nsend;
    }
    if(minus[ic]==MyContext)
    {
      MPI_Status st;
      MPI_Recv(buff.data(),2*buff.size(),MPI_DOUBLE,plus[ic],plus[ic]+999,comm,&st);
      wlk.push_walkers(1,buff);
    }
  }
  return nswap; 
}

/** swap Walkers with Irecv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
template<class WlkBucket, 
         class IVec = std::vector<int>
         >
// eventually generalize MPI_Comm to a MPI wrapper
inline int swapWalkersAsync(WlkBucket& wlk, IVec& CurrNumPerNode, IVec& NewNumPerNode, MPI_Comm comm)
{
  int NumContexts, MyContext;
  MPI_Comm_size(comm,&NumContexts);
  MPI_Comm_rank(comm,&MyContext);
  assert(CurrNumPerNode.size() >= NumContexts);
  assert(NewNumPerNode.size() >= NumContexts);
  auto Cur_pop = std::accumulate(CurrNumPerNode.begin(), CurrNumPerNode.begin()+NumContexts, 0);
  std::vector<int> minus, plus;
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    int dn=CurrNumPerNode[ip]-NewNumPerNode[ip];
    if(ip == MyContext)
      deltaN=dn;
    if(dn>0)
    {
      plus.insert(plus.end(),dn,ip);
    }
    else
      if(dn<0)
      {
        minus.insert(minus.end(),-dn,ip);
      }
  }
  int nswap=std::min(plus.size(), minus.size());
  int nsend=0;
  int wlk_size = wlk.single_walker_size();
  int countSend = 1;
  std::vector<ComplexType*> buffers;
  std::vector<MPI_Request> requests;
  std::vector<int> recvCounts;
  for(int ic=0; ic<nswap; ic++)
  {
    if(plus[ic]==MyContext)
    {
      if((ic < nswap - 1) && (plus[ic] == plus[ic+1]) && (minus[ic] == minus[ic+1]))
      {
        countSend++;
      }
      else
      {
        ComplexType* bf = new ComplexType[countSend*wlk_size];   
        buffers.push_back(bf);
        wlk.pop_walkers(countSend,bf);
        requests.push_back(MPI_Request());
        MPI_Isend(bf,2*countSend*wlk_size,MPI_DOUBLE,minus[ic],plus[ic]+1999,comm,std::addressof(requests.back()));
        nsend += countSend;
        countSend = 1;
      }
    }
    if(minus[ic]==MyContext)
    {
      if((ic < nswap - 1) && (plus[ic] == plus[ic+1]) && (minus[ic] == minus[ic+1]))
      {
        countSend++;
      }
      else
      {
        ComplexType* bf = new ComplexType[countSend*wlk_size];
        buffers.push_back(bf);
        requests.push_back(MPI_Request());
        recvCounts.push_back(countSend);
        MPI_Irecv(bf,2*countSend*wlk_size,MPI_DOUBLE,plus[ic],plus[ic]+1999,comm,std::addressof(requests.back()));
        countSend = 1;
      }
    }
  }
  if(deltaN < 0) {
    // receiving nodes
    MPI_Status st;
    for(int ip = 0; ip < requests.size(); ++ip)
    {
      MPI_Wait(std::addressof(requests[ip]),std::addressof(st));
      wlk.push_walkers(recvCounts[ip],buffers[ip]);
      delete[] buffers[ip];
    }
  } else {
    // sending nodes
    MPI_Status st;
    for(int ip = 0; ip < requests.size(); ++ip)
    {
      MPI_Wait(std::addressof(requests[ip]),std::addressof(st));
      delete[] buffers[ip]; 
    }
  }
  return nswap;
}

/** swap Walkers with Recv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
/*
void WalkerControlMPI::swapWalkersBlocked(MCWalkerConfiguration& W)
{
  OffSet[0]=0;
  for(int i=0; i<NumContexts; i++)
    OffSet[i+1]=OffSet[i]+NumPerNode[i];
  FairDivideLow(Cur_pop,NumContexts,FairOffSet);
  int toLeft=FairOffSet[MyContext]-OffSet[MyContext];
  int toRight=FairOffSet[MyContext+1]-OffSet[MyContext+1];
  Walker_t& wRef(*W[0]);
  std::vector<Walker_t*> newW;
  int num_deleted=0;
  int last=NumPerNode[MyContext]-1;
  //scehdule irecv
  if(toLeft<0)
    // recv from node-1
  {
    int dn=-toLeft;
    //OOMPI_Packed recvBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    //OOMPI_COMM_WORLD[MyContext-1].Recv(recvBuffer);
    OOMPI_Packed recvBuffer(dn*wRef.byteSize(),myComm->getComm());
    myComm->getComm()[MyContext-1].Recv(recvBuffer);
    while(dn)
    {
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(recvBuffer);
      newW.push_back(awalker);
      --dn;
    }
  }
  else
    if(toLeft>0)
      //send to node-1
    {
      int dn=toLeft;
      num_deleted+=dn;
      //OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
      OOMPI_Packed sendBuffer(dn*wRef.byteSize(),myComm->getComm());
      while(dn)
      {
        W[last]->putMessage(sendBuffer);
        --dn;
        --last;
      }
      //OOMPI_COMM_WORLD[MyContext-1].Send(sendBuffer);
      myComm->getComm()[MyContext-1].Send(sendBuffer);
    }
  if(toRight<0)
    // send to node+1
  {
    int dn=-toRight;
    num_deleted+=dn;
    //OOMPI_Packed sendBuffer(dn*wRef.byteSize(),OOMPI_COMM_WORLD);
    OOMPI_Packed sendBuffer(dn*wRef.byteSize(),myComm->getComm());
    while(dn)
    {
      W[last]->putMessage(sendBuffer);
      --dn;
      --last;
    }
    //OOMPI_COMM_WORLD[MyContext+1].Send(sendBuffer);
    myComm->getComm()[MyContext+1].Send(sendBuffer);
  }
  else
    if(toRight>0)
      //recv from node+1
    {
      //OOMPI_Packed recvBuffer(toRight*wRef.byteSize(),OOMPI_COMM_WORLD);
      //OOMPI_COMM_WORLD[MyContext+1].Recv(recvBuffer);
      OOMPI_Packed recvBuffer(toRight*wRef.byteSize(),myComm->getComm());
      myComm->getComm()[MyContext+1].Recv(recvBuffer);
      int dn=toRight;
      while(dn)
      {
        Walker_t *awalker= new Walker_t(wRef);
        awalker->getMessage(recvBuffer);
        newW.push_back(awalker);
        --dn;
      }
    }
  while(num_deleted>0)
  {
    W.pop_back();
    --num_deleted;
  }
  if(newW.size())
    W.insert(W.end(),newW.begin(),newW.end());
}
*/

}

}
