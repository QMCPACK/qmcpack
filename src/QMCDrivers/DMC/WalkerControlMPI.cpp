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
  DMC_MPI_loadbalance,
  DMC_MPI_send,
  DMC_MPI_recv,
};

TimerNameList_t<DMC_MPI_Timers> DMCMPITimerNames =
{
  {DMC_MPI_branch, "WalkerControlMPI::branch"},
  {DMC_MPI_imbalance, "WalkerControlMPI::imbalance"},
  {DMC_MPI_prebalance, "WalkerControlMPI::pre-loadbalance"},
  {DMC_MPI_copyWalkers, "WalkerControlMPI::copyWalkers"},
  {DMC_MPI_allreduce, "WalkerControlMPI::allreduce"},
  {DMC_MPI_loadbalance, "WalkerControlMPI::loadbalance"},
  {DMC_MPI_send, "WalkerControlMPI::send"},
  {DMC_MPI_recv, "WalkerControlMPI::recv"}
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

/** Perform branch and swap walkers as required
 *
 *  It takes 4 steps:
 *    1. sortWalkers marks good and bad walkers.
 *    2. allreduce and make the decision of load balancing.
 *    3. send/recv walkers. Receiving side recycles bad walkers' memory first.
 *    4. copyWalkers generates copies of good walkers.
 *  In order to minimize the memory footprint fluctuation
 *  the walker copying is placed as the last step.
 *  In order to reduce the time for allocating walker memory,
 *  this algorithm does not destroy the bad walkers in step 1.
 *  All the bad walkers are recycled as much as possible in step 3/4.
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
  myTimers[DMC_MPI_loadbalance]->start();
  swapWalkersSimple(W);
  myTimers[DMC_MPI_loadbalance]->stop();
  myTimers[DMC_MPI_copyWalkers]->start();
  copyWalkers(W);
  myTimers[DMC_MPI_copyWalkers]->stop();
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

/** swap Walkers with Recv/Send or Irecv/Isend
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * Each MPI rank can only send or receive or be silent.
 * The communication is one-dimensional and very local.
 * If multiple copies of a walker need to be sent to the target rank,
 * only one walker is sent and the number of copies is encoded in the message.
 * The blocking send/recv may become serialized and worsen load imbalance.
 * Non blocking send/recv algorithm avoids serialization completely.
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
  if(plus.size()!=minus.size())
  {
    app_error() << "Walker send/recv pattern doesn't match. "
                << "The send size " << plus.size()
                << " is not equal to the receive size " << minus.size()
                << " ." << std::endl;
    APP_ABORT("WalkerControlMPI::swapWalkersSimple");
  }
  int nswap=plus.size();
  // sort good walkers by the number of copies
  assert(good_w.size()==ncopy_w.size());
  std::vector<std::pair<int,int> > ncopy_pairs, nrecv_pairs;
  for(int iw=0; iw<ncopy_w.size(); iw++)
    ncopy_pairs.push_back(std::make_pair(ncopy_w[iw],iw));
  std::sort(ncopy_pairs.begin(), ncopy_pairs.end());

  int nsend=0;
  std::vector<OOMPI_Request> requests;
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
        { // not enough copies to send or not the same send/recv pair
          break;
        }

      // pack data and send
      size_t byteSize = awalker->byteSize();
      awalker->updateBuffer();
      OOMPI_Message sendBuffer(awalker->DataSet.data(), byteSize);
      if(use_nonblocking)
        requests.push_back(myComm->getComm()[minus[ic]].Isend(sendBuffer));
      else
      {
        myTimers[DMC_MPI_send]->start();
        myComm->getComm()[minus[ic]].Send(sendBuffer);
        myTimers[DMC_MPI_send]->stop();
      }
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
      // count receive pairs, (source,copy)
      nrecv_pairs.push_back(std::make_pair(plus[ic],0));
      for(int id=ic+1; id<nswap; id++)
        if(plus[ic]==plus[id]&&minus[ic]==minus[id])
          nrecv_pairs.back().second++;
        else
          break;
      // update cursor
      ic+=nrecv_pairs.back().second;
    }
  }

  if(nsend>0)
  {
    if(use_nonblocking)
    {
      // wait all the isend
      for(int im=0; im<requests.size(); im++)
      {
        myTimers[DMC_MPI_send]->start();
        requests[im].Wait();
        myTimers[DMC_MPI_send]->stop();
      }
      requests.clear();
    }
  }
  else
  {
    struct job {
      OOMPI_Request request;
      int walkerID;
      int queueID;
      job(const OOMPI_Request &req, int wid, int qid): request(req), walkerID(wid), queueID(qid) {};
    };
    std::vector<job> job_list;
    std::vector<bool> queue_status(nrecv_pairs.size(),true);

    bool completed=false;
    while(!completed)
    {
      // receive data
      for(int ic=0; ic<nrecv_pairs.size(); ic++)
        if(queue_status[ic]&&nrecv_pairs[ic].second>=0)
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

          size_t byteSize = awalker->byteSize();
          myTimers[DMC_MPI_recv]->start();
          OOMPI_Message recvBuffer(awalker->DataSet.data(), byteSize);
          if(use_nonblocking)
          {
            job_list.push_back(job(myComm->getComm()[nrecv_pairs[ic].first].Irecv(recvBuffer),
                               newW.size(),ic));
            queue_status[ic]=false;
          }
          else
          {
            myComm->getComm()[nrecv_pairs[ic].first].Recv(recvBuffer);
            job_list.push_back(job(OOMPI_Request(),newW.size(),ic));
          }
          myTimers[DMC_MPI_recv]->stop();

          newW.push_back(awalker);
          ncopy_newW.push_back(0);
        }

      if(use_nonblocking)
      {
        OOMPI_Status status;
        for(auto jobit=job_list.begin(); jobit!=job_list.end(); jobit++)
          if(jobit->request.Test(status))
          {
            auto &awalker=newW[jobit->walkerID];
            // unpack data
            awalker->copyFromBuffer();
            ncopy_newW[jobit->walkerID]=awalker->NumSentCopies;
            // update counter
            nrecv_pairs[jobit->queueID].second-=(awalker->NumSentCopies+1);
            queue_status[jobit->queueID]=true;
            job_list.erase(jobit);
            break;
          }
      }
      else
      {
        for(auto jobit=job_list.begin(); jobit!=job_list.end(); jobit++)
        {
          auto &awalker=newW[jobit->walkerID];
          // unpack data
          awalker->copyFromBuffer();
          ncopy_newW[jobit->walkerID]=awalker->NumSentCopies;
          // update counter
          nrecv_pairs[jobit->queueID].second-=(awalker->NumSentCopies+1);
        }
        job_list.clear();
      }

      // check the completion of queues
      completed=true;
      for(int ic=0; ic<nrecv_pairs.size(); ic++)
        completed = completed && (nrecv_pairs[ic].second==-1);
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

}

