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
    
    
#include<tuple>
#include<cassert>
#include <memory>
#include <stack>
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

/** 
 * Implements the paired branching algorithm on a popultion of walkers.
 *   - buff: array of walker info (weight,num).
 */ 
template<class Random 
         >
inline void SerialPairBranching(std::vector<char>& buff, Random& rng, double max_c, double min_c)
{
  typedef std::tuple<double,int,int>  tp; 
  typedef std::vector<tp>::iterator  tp_it; 
// slow for now, not efficient!!!
  int nw = buff.size()/(sizeof(double)+sizeof(int));
  std::vector<tp> wlks(nw);
  double e;
  for(int i=0,sz=0,ni=0; i<nw; i++) {
    std::memcpy(&e,buff.data()+sz,sizeof(double));  
    sz+=sizeof(double)+sizeof(int);
    std::get<0>(wlks[i]) = e;
    std::get<1>(wlks[i]) = 1;
    std::get<2>(wlks[i]) = i;
  }

  std::sort( wlks.begin(), wlks.end(),  
             [] (const tp& a, const tp& b) {
               return std::get<0>(a) < std::get<0>(b);
             }
  );

  tp_it it_s = wlks.begin();
  tp_it it_l = wlks.end()-1;
 
  while( it_s < it_l ) {

    if( std::abs(std::get<0>(*it_s)) < min_c || std::abs(std::get<0>(*it_l)) > max_c) {
      double w12 = std::get<0>(*it_s) + std::get<0>(*it_l); 
      if(rng() < std::get<0>(*it_l)/w12) {
        std::get<0>(*it_l) = 0.5*w12;
        std::get<0>(*it_s) = 0.0;
        std::get<1>(*it_l) = 2;      
        std::get<1>(*it_s) = 0;      
      } else {
        std::get<0>(*it_s) = 0.5*w12;
        std::get<0>(*it_l) = 0.0;
        std::get<1>(*it_s) = 2;      
        std::get<1>(*it_l) = 0;      
      }
      it_s++;
      it_l--;
    } else 
      break;    

  }

  int sz = sizeof(double) + sizeof(int);
  int nnew=0;
  for(int i=0,ni=0,ic=0; i<nw; i++) {
    e = std::get<0>(wlks[i]);
    ni = std::get<1>(wlks[i]);
    ic = std::get<2>(wlks[i]);
    nnew += ni;
    std::memcpy(buff.data()+ic*sz,&e,sizeof(double));
    std::memcpy(buff.data()+ic*sz+sizeof(double),&ni,sizeof(int));
  }
  assert(nw==nnew);
 
}

/** Pair branching with with Irecv/Send
 *    THIS ROUTINE DOESN'T WORK !!! ONLY FOR TESTING PURPOSES !!!
 */
template<class WlkBucket, 
         class Random // = std::uniform_real_distribution<double>,
         >
inline int NaivePairBranching(WlkBucket& wlk, Random& rng, MPI_Comm comm)
{
  int NumContexts, MyContext;
  MPI_Comm_size(comm,&NumContexts);
  MPI_Comm_rank(comm,&MyContext);

  // local pair branching 
  std::vector<int> windx; // this are the positions of walkers that are candidates for global pair branching
  int sg = wlk.local_pair_branching(windx);  // sg: +0 for large weights, -1 for small weights
  int nw2b = windx.size()*sg;
//std::cout<<MyContext <<" " <<nw2b <<std::endl;

  std::vector<int> WCnt(NumContexts);
  MPI_Allgather(&nw2b, 1, MPI_INT, WCnt.data(), 1, MPI_INT, comm);

  std::vector<int> minus, plus;
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    int dn=WCnt[ip];
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
  int nleft = plus.size()-minus.size();
  if(nswap == 0)
    return 0;

  int nsend=0;
  int wlk_size = wlk.single_walker_size();
  int countSend = 1, wlkPos=0;
  std::vector<ComplexType*> buffers;
  std::vector<MPI_Request> requests;
  std::vector<double> wgts_out(2*windx.size()); 
  std::vector<double> wgts_in(2*windx.size()); 
  requests.reserve(2*windx.size());
  buffers.reserve(2*windx.size());
  for(int i=0; i<windx.size(); i++) 
  {
    wgts_out[2*i] = wlk.getWeight(windx[i]).real();
    wgts_out[2*i+1] = rng(); 
  }
  // send weigts and random numbers (on plus side)
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
        requests.push_back(MPI_Request());
        MPI_Isend(wgts_out.data()+2*wlkPos,2*countSend,MPI_DOUBLE,minus[ic],plus[ic]+1999,comm,std::addressof(requests.back()));
        requests.push_back(MPI_Request());
        MPI_Irecv(wgts_in.data()+2*wlkPos,2*countSend,MPI_DOUBLE,minus[ic],plus[ic]+999,comm,std::addressof(requests.back()));
        wlkPos+=countSend;
        countSend = 1;
        nsend++;
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
        requests.push_back(MPI_Request());
        MPI_Isend(wgts_out.data()+2*wlkPos,2*countSend,MPI_DOUBLE,plus[ic],plus[ic]+999,comm,std::addressof(requests.back()));
        requests.push_back(MPI_Request());
        MPI_Irecv(wgts_in.data()+2*wlkPos,2*countSend,MPI_DOUBLE,plus[ic],plus[ic]+1999,comm,std::addressof(requests.back()));
        wlkPos+=countSend;
        countSend = 1;
        nsend++;
      }
    }
  }  

  MPI_Status st;
  std::vector<int> exch(nsend);
  std::stack<int> branch, remv;
  std::stack<int> gremv;
  wlkPos=0; 
  nsend=0;
  countSend = 1;
  // now perform branching and exchange data 
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
        // wait for Irecv
        MPI_Wait(std::addressof(requests[2*nsend+1]),std::addressof(st)); 
        // count exchanges
        assert(branch.size()==0);
        assert(remv.size()==0);
        for(int i=0; i<countSend; i++) {
          double w12 = wgts_out[2*(wlkPos+i)] + wgts_in[2*(wlkPos+i)]; 
          if( wgts_out[2*(wlkPos+i)+1] < wgts_out[2*(wlkPos+i)]/w12 ) {
            // branch large weight 
            branch.push(wlkPos+i);
          } else {
            // branch small weight 
            remv.push(wlkPos+i);
          }
        }
        int nlocal = std::min(branch.size(),remv.size());
        for(int i=0; i<nlocal; i++) {
          int ib = branch.top();
          int ir = remv.top();
          double w12 = 0.5*(wgts_out[2*ib] + wgts_in[2*ib]);
          wlk.pair_branch(w12,windx[ib],windx[ir]);
          branch.pop();
          remv.pop();
        } 
        if(branch.size() > 0) {
          exch[nsend] = branch.size();

          ComplexType* bf = new ComplexType[exch[nsend]*wlk_size];
          buffers.push_back(bf);
          int bsz = branch.size();
          for(int i=0; i<bsz; i++) {
            int ib = branch.top();
            branch.pop();
            double w12 = 0.5*(wgts_out[2*ib] + wgts_in[2*ib]);
            wlk.setWeight(windx[ib],w12);
            wlk.copy_to_buffer(windx[ib],bf+i*wlk_size);
          }
          MPI_Isend(bf,2*exch[nsend]*wlk_size,MPI_DOUBLE,minus[ic],plus[ic]+2999,comm,std::addressof(requests[2*nsend+1]));
        } else if(remv.size() > 0) {
          exch[nsend] = remv.size();

          for(int i=0; i<exch[nsend]; i++) {
            gremv.push( windx[remv.top()] );
            remv.pop();
          } 
 
          ComplexType* bf = new ComplexType[exch[nsend]*wlk_size];
          buffers.push_back(bf);
          MPI_Irecv(bf,2*exch[nsend]*wlk_size,MPI_DOUBLE,minus[ic],plus[ic]+3999,comm,std::addressof(requests[2*nsend+1]));

          exch[nsend] *= -1;
        } else 
          exch[nsend] = 0;

        nsend++;
        wlkPos+=countSend;
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
        // wait for Irecv
        MPI_Wait(std::addressof(requests[2*nsend+1]),std::addressof(st)); 
        // count exchanges
        assert(branch.size()==0);
        assert(remv.size()==0);
        for(int i=0; i<countSend; i++) {
          double w12 = wgts_out[2*(wlkPos+i)] + wgts_in[2*(wlkPos+i)]; 
          if( wgts_in[2*(wlkPos+i)+1] < wgts_in[2*(wlkPos+i)]/w12 ) {
            // branch large weight 
            remv.push(wlkPos+i);
          } else {
            // branch small weight 
            branch.push(wlkPos+i);
          }
        }
        int nlocal = std::min(branch.size(),remv.size());
        for(int i=0; i<nlocal; i++) {
          int ib = branch.top();
          int ir = remv.top();
          branch.pop();
          remv.pop();
          double w12 = 0.5*(wgts_out[2*ib] + wgts_in[2*ib]);
          wlk.pair_branch(w12,windx[ib],windx[ir]);
        } 
        if(branch.size() > 0) {
          exch[nsend] = branch.size();

          ComplexType* bf = new ComplexType[exch[nsend]*wlk_size];
          buffers.push_back(bf);
          int bsz = branch.size();
          for(int i=0; i<bsz; i++) {
            int ib = branch.top();
            branch.pop();
            double w12 = 0.5*(wgts_out[2*ib] + wgts_in[2*ib]);
            wlk.setWeight(windx[ib],w12);
            wlk.copy_to_buffer(windx[ib],bf+i*wlk_size);
          }
          MPI_Isend(bf,2*exch[nsend]*wlk_size,MPI_DOUBLE,plus[ic],plus[ic]+3999,comm,std::addressof(requests[2*nsend+1]));
        } else if(remv.size() > 0) {
          exch[nsend] = remv.size();

          for(int i=0; i<exch[nsend]; i++) { 
            gremv.push( windx[remv.top()] );
            remv.pop();
          } 

          ComplexType* bf = new ComplexType[exch[nsend]*wlk_size];
          buffers.push_back(bf);
          MPI_Irecv(bf,2*exch[nsend]*wlk_size,MPI_DOUBLE,plus[ic],plus[ic]+2999,comm,std::addressof(requests[2*nsend+1]));

          exch[nsend] *= -1;
        } else 
          exch[nsend] = 0;

        nsend++;
        wlkPos+=countSend;
        countSend = 1;
      }
    }
  }

  // now receive messages and Wait for communications
  nsend=0;
  for(int ic=0,bfcnt=0; ic<exch.size(); ic++)
  {
    // wait for original ISend
    MPI_Wait(std::addressof(requests[2*ic]),std::addressof(st)); 
    if(exch[ic] > 0) {
      // just wait for ISend and delete buffer
      MPI_Wait(std::addressof(requests[2*ic+1]),std::addressof(st)); 
      delete[] buffers[bfcnt++]; 
    } else if(exch[ic] < 0) {
      MPI_Wait(std::addressof(requests[2*ic+1]),std::addressof(st)); 
      int nw = -1*exch[ic];
      for(int i=0; i<nw; i++) {
        int ir = gremv.top();
        gremv.pop();
        wlk.copy_from_buffer(ir,buffers[bfcnt]+i*wlk_size);
      }
      delete[] buffers[bfcnt++];
    }
  }

  // stack of deleted walkers should be empty
  assert(gremv.size()==0);
  
  return nleft;

}


}

}

