//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <QMCDrivers/DMC/WalkerControlMPI.h>
#include <qmc_common.h>
#include <Utilities/IteratorUtility.h>
#include <Utilities/UtilityFunctions.h>
#include <Utilities/NewTimer.h>
#include <Utilities/Timer.h>

namespace qmcplusplus
{

//#define MCWALKERSET_MPI_DEBUG

enum DMC_MPI_Timers
{
  DMC_MPI_branch,
  DMC_MPI_imbalance,
  DMC_MPI_prebalance,
  DMC_MPI_copyWalkers,
  DMC_MPI_allreduce,
  DMC_MPI_loadbalance
};

TimerNameList_t<DMC_MPI_Timers> DMCMPITimerNames =
{
  {DMC_MPI_branch, "WalkerControlMPI::branch"},
  {DMC_MPI_imbalance, "WalkerControlMPI::imbalance"},
  {DMC_MPI_prebalance, "WalkerControlMPI::pre-loadbalance"},
  {DMC_MPI_copyWalkers, "WalkerControlMPI::copyWalkers"},
  {DMC_MPI_allreduce, "WalkerControlMPI::allreduce"},
  {DMC_MPI_loadbalance, "WalkerControlMPI::loadbalance"}
};

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
  std::ofstream fout(fname);
#endif
  setup_timers(myTimers, DMCMPITimerNames, timer_level_medium);
}

/** perform branch and swap walkers as required
 *
 *  it takes 4 steps:
 *    1. shortWalkers marks good and bad walkers.
 *    2. allreduce and make the decision of load balancing.
 *    3. send/recv walkers. Receiving side recycle bad walkers' memory first.
 *    4. copyWalkers generate walker copies of good walkers.
 */
int WalkerControlMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  myTimers[DMC_MPI_branch]->start();
  myTimers[DMC_MPI_prebalance]->start();
  std::fill(curData.begin(),curData.end(),0);
  sortWalkers(W);
  //use NumWalkersSent from the previous exchange
  curData[SENTWALKERS_INDEX]=NumWalkersSent;
  //update the number of walkers for this node
  curData[LE_MAX+MyContext]=NumWalkers;
  //myTimers[DMC_MPI_imbalance]->start();
  //myComm->barrier();
  //myTimers[DMC_MPI_imbalance]->stop();
  myTimers[DMC_MPI_allreduce]->start();
  myComm->allreduce(curData);
  myTimers[DMC_MPI_allreduce]->stop();
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  Cur_pop=0;
  for(int i=0, j=LE_MAX; i<NumContexts; i++,j++)
  {
    Cur_pop+= NumPerNode[i]=static_cast<int>(curData[j]);
  }
  myTimers[DMC_MPI_prebalance]->stop();
  if(qmc_common.async_swap)
  {
    myTimers[DMC_MPI_copyWalkers]->start();
    copyWalkers(W);
    //myComm->barrier();
    myTimers[DMC_MPI_copyWalkers]->stop();
    myTimers[DMC_MPI_loadbalance]->start();
    swapWalkersAsync(W);
    myTimers[DMC_MPI_loadbalance]->stop();
  }
  else
  {
    myTimers[DMC_MPI_loadbalance]->start();
    swapWalkersSimple(W);
    //myComm->barrier();
    myTimers[DMC_MPI_loadbalance]->stop();
    myTimers[DMC_MPI_copyWalkers]->start();
    copyWalkers(W);
    //myComm->barrier();
    myTimers[DMC_MPI_copyWalkers]->stop();
  }
  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
  //update the global number of walkers and offsets
  W.setGlobalNumWalkers(Cur_pop);
  W.setWalkerOffsets(FairOffSet);

  // walker count safety check
  //if(W.getActiveWalkers()!=FairOffSet[MyContext+1]-FairOffSet[MyContext])
  //  std::cout << " Strange on rank " << MyContext
  //            << " should be " << FairOffSet[MyContext+1]-FairOffSet[MyContext]
  //            << " walker but actually " << W.getActiveWalkers() << std::endl;
  myTimers[DMC_MPI_branch]->stop();
  return Cur_pop;
}

// determine new walker population on each node
void determineNewWalkerPopulation(int Cur_pop, int NumContexts, int MyContext, const std::vector<int> &NumPerNode, std::vector<int> &FairOffSet, std::vector<int> &minus, std::vector<int> &plus)
{
  // Cur_pop - in - current population
  // NumContexts -in - number of MPI processes
  // MyContext - in - my MPI rank
  // NumPerNode - in - current walkers per node
  // FairOffSet - out - new walker offset
  // minus - out - number of walkers to be removed from each node
  // plus -  out - number of walkers to be added to each node

  FairDivideLow(Cur_pop,NumContexts,FairOffSet);
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    int dn=NumPerNode[ip]-(FairOffSet[ip+1]-FairOffSet[ip]);
    if(ip == MyContext)
      deltaN=dn;
    if(dn>0)
    {
      plus.insert(plus.end(),dn,ip);
    }
    else if(dn<0)
    {
      minus.insert(minus.end(),-dn,ip);
    }
  }
}

/** swap Walkers with Recv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
void WalkerControlMPI::swapWalkersSimple(MCWalkerConfiguration& W)
{
  std::vector<int> minus, plus;
  determineNewWalkerPopulation(Cur_pop, NumContexts, MyContext, NumPerNode, FairOffSet, minus, plus);

  Walker_t& wRef(*good_w[0]);
  std::vector<Walker_t*> newW;
  std::vector<int> ncopy_newW;
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  std::ofstream fout(fname, std::ios::app);
  //fout << NumSwaps << " " << Cur_pop << " ";
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerNode[ic] << " ";
  //fout << " | ";
  //for(int ic=0; ic<NumContexts; ic++) fout << FairOffSet[ic+1]-FairOffSet[ic] << " ";
  //fout << " | ";
  for(int ic=0; ic<plus.size(); ic++)
  {
    fout << plus[ic] << " ";
  }
  fout << " | ";
  for(int ic=0; ic<minus.size(); ic++)
  {
    fout << minus[ic] << " ";
  }
  fout << std::endl;
#endif
  int nswap=std::min(plus.size(), minus.size());

  // sort good walkers by the number of copies
  assert(good_w.size()==ncopy_w.size());
  std::vector<std::pair<int,int> > ncopy_pairs;
  for(int iw=0; iw<ncopy_w.size(); iw++)
    ncopy_pairs.push_back(std::make_pair(ncopy_w[iw],iw));
  std::sort(ncopy_pairs.begin(), ncopy_pairs.end());

  int nsend=0;
  for(int ic=0; ic<nswap; ic++)
  {
    if(plus[ic]==MyContext)
    {
      // always send the last good walker
      Walker_t* &awalker = good_w[ncopy_pairs.back().second];

      // count the possible copies in one send
      auto &nsentcopy = awalker->NumSentCopies;
      nsentcopy = 0;

      for(int id=ic+1; id<nswap; id++)
        if(plus[ic]==plus[id]&&minus[ic]==minus[id]&&ncopy_pairs.back().first>0)
        { // increment copy counter
          ncopy_pairs.back().first--;
          nsentcopy++;
        }
        else
        { // not enough copies to send or not the same send/recv pattern
          break;
        }

      // pack data and send
      size_t byteSize = awalker->byteSize();
      awalker->updateBuffer();
      OOMPI_Message sendBuffer(awalker->DataSet.data(), byteSize);
      myComm->getComm()[minus[ic]].Send(sendBuffer);
      // update counter and cursor
      ++nsend;
      ic+=nsentcopy;

      // update copy counter
      if(ncopy_pairs.back().first>0)
      {
        ncopy_pairs.back().first--;
        std::sort(ncopy_pairs.begin(), ncopy_pairs.end());
      }
      else
      {
        ncopy_pairs.pop_back();
        bad_w.push_back(awalker);
      }
    }
    if(minus[ic]==MyContext)
    {
      Walker_t* awalker;
      if(bad_w.empty())
      {
        awalker=new Walker_t(wRef);
      }
      else
      {
        awalker=bad_w.back();
        bad_w.pop_back();
      }

      // receive and unpack data
      size_t byteSize = awalker->byteSize();
      OOMPI_Message recvBuffer(awalker->DataSet.data(), byteSize);
      myComm->getComm()[plus[ic]].Recv(recvBuffer);
      awalker->copyFromBuffer();
      newW.push_back(awalker);
      ncopy_newW.push_back(awalker->NumSentCopies);

      // update counter and cursor
      ic+=awalker->NumSentCopies;
    }
  }
  //save the number of walkers sent
  NumWalkersSent=nsend;
  // rebuild good_w and ncopy_w
  std::vector<Walker_t*> good_w_temp(good_w);
  good_w.resize(ncopy_pairs.size());
  ncopy_w.resize(ncopy_pairs.size());
  for(int iw=0; iw<ncopy_pairs.size(); iw++)
  {
    good_w[iw]=good_w_temp[ncopy_pairs[iw].second];
    ncopy_w[iw]=ncopy_pairs[iw].first;
  }
  //add walkers from other node
  if(newW.size())
  {
    good_w.insert(good_w.end(),newW.begin(),newW.end());
    ncopy_w.insert(ncopy_w.end(),ncopy_newW.begin(),ncopy_newW.end());
  }
}

/** swap Walkers with Irecv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
void WalkerControlMPI::swapWalkersAsync(MCWalkerConfiguration& W)
{
  std::vector<int> minus, plus;
  determineNewWalkerPopulation(Cur_pop, NumContexts, MyContext, NumPerNode, FairOffSet, minus, plus);
  Walker_t& wRef(*W[0]);
  std::vector<Walker_t*> newW;
  std::vector<Walker_t*> oldW;
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  std::ofstream fout(fname, std::ios::app);
  //fout << NumSwaps << " " << Cur_pop << " ";
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerNode[ic] << " ";
  //fout << " | ";
  //for(int ic=0; ic<NumContexts; ic++) fout << FairOffSet[ic+1]-FairOffSet[ic] << " ";
  //fout << " | ";
  for(int ic=0; ic<plus.size(); ic++)
  {
    fout << plus[ic] << " ";
  }
  fout << " | ";
  for(int ic=0; ic<minus.size(); ic++)
  {
    fout << minus[ic] << " ";
  }
  fout << std::endl;
#endif
  int nswap=std::min(plus.size(), minus.size());
  int last=W.getActiveWalkers()-1;
  int nsend=0;
  int countSend = 1;
  OOMPI_Packed ** sendBuffers = new OOMPI_Packed*[NumContexts];
  OOMPI_Packed ** recvBuffers = new OOMPI_Packed*[NumContexts];
  std::vector<OOMPI_Request> requests(NumContexts);
  std::vector<int> sendCounts(NumContexts,0);
  for(int ip=0; ip < NumContexts; ++ip)
    sendBuffers[ip] = 0;
  for(int ip=0; ip < NumContexts; ++ip)
    recvBuffers[ip] = 0;
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
        //OOMPI_Packed sendBuffer(wRef.byteSize(),OOMPI_COMM_WORLD);
        sendBuffers[minus[ic]] = new OOMPI_Packed(countSend * wRef.byteSize(),myComm->getComm());
        for(int cs = 0; cs < countSend; ++cs)
        {
          W[last]->putMessage(*(sendBuffers[minus[ic]]));
          --last;
        }
        //OOMPI_COMM_WORLD[minus[ic]].Send(sendBuffer);
        requests[minus[ic]] = myComm->getComm()[minus[ic]].Isend(*(sendBuffers[minus[ic]]), plus[ic]);
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
        //OOMPI_Packed recvBuffer(wRef.byteSize(),OOMPI_COMM_WORLD);
        recvBuffers[plus[ic]] = new OOMPI_Packed(countSend * wRef.byteSize(),myComm->getComm());
        //OOMPI_COMM_WORLD[plus[ic]].Recv(recvBuffer);
        requests[plus[ic]] = myComm->getComm()[plus[ic]].Irecv(*(recvBuffers[plus[ic]]), plus[ic]);
        sendCounts[plus[ic]] = countSend;
        countSend = 1;
      }
    }
  }
  for(int ip = 0; ip < NumContexts; ++ip)
  {
    if(recvBuffers[ip])
    {
      requests[ip].Wait();
      for(int cs = 0; cs < sendCounts[ip]; ++cs)
      {
        Walker_t *awalker= new Walker_t(wRef);
        awalker->getMessage(*(recvBuffers[ip]));
        newW.push_back(awalker);
      }
      delete recvBuffers[ip];
    }
  }
  for(int ip = 0; ip < NumContexts; ++ip)
  {
    if(sendBuffers[ip])
    {
      requests[ip].Wait();
      delete sendBuffers[ip];
    }
  }
  delete[] sendBuffers;
  delete[] recvBuffers;
  //save the number of walkers sent
  NumWalkersSent=nsend;
  if(nsend)
  {
    nsend=NumPerNode[MyContext]-nsend;
    W.destroyWalkers(W.begin()+nsend, W.end());
  }
  //add walkers from other node
  if(newW.size())
    W.insert(W.end(),newW.begin(),newW.end());
}


}

