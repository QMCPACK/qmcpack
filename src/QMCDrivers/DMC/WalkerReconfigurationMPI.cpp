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
#include "QMCDrivers/DMC/WalkerReconfigurationMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
WalkerReconfigurationMPI::WalkerReconfigurationMPI(): TotalWalkers(0) {
  SwapMode=1;
  NumContexts=OHMMS::Controller->ncontexts();
  MyContext=OHMMS::Controller->mycontext();
  UnitZeta=Random();

  MPI_Bcast(&UnitZeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  app_log() << "  First weight [0,1) for reconfiguration =" << UnitZeta << endl;

  //ostringstream o;
  //o << "check." << MyContext << ".dat";
  //ofstream fout(o.str().c_str());
  //fout << "UnitZeta " << UnitZeta << endl;
}

int 
WalkerReconfigurationMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger) {

  OHMMS::Controller->barrier();

  int nwkept = swapWalkers(W);

  OHMMS::Controller->barrier();

  gsum(nwkept,0);

  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end) {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }

  return nwkept;
}

int WalkerReconfigurationMPI::swapWalkers(MCWalkerConfiguration& W) {
  //ostringstream o;
  //o << "check." << MyContext << ".dat";
  //ofstream fout(o.str().c_str(),ios::app);

  int nw=W.getActiveWalkers();
  if(TotalWalkers ==0) {
    FirstWalker=nw*MyContext;
    LastWalker=FirstWalker+nw;
    TotalWalkers = nw*NumContexts;
    nwInv = 1.0/static_cast<RealType>(TotalWalkers);
    DeltaStep = UnitZeta*nwInv;

    ncopy_w.resize(nw);
    wConf.resize(nw);
    wSum.resize(NumContexts);
    wOffset.resize(NumContexts+1);
  }

  std::fill(wSum.begin(),wSum.end(),0.0);

  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  RealType wtot=0.0;
  int iw=0;
  while(it != it_end) {
    wtot+=wConf[iw++]=(*it)->Weight;
    ++it;
  }

  wSum[MyContext]=wtot;
  gsum(wSum,0);

  //wOffset[ip] is the partial sum update to ip
  wOffset[0]=0;
  for(int ip=0; ip<NumContexts; ip++) wOffset[ip+1]=wOffset[ip]+wSum[ip];

  //wtot is the total weight
  wtot=wOffset[NumContexts];

  //find the lower and upper bound of index
  int minIndex=(wOffset[MyContext]/wtot-DeltaStep)*static_cast<RealType>(TotalWalkers)-1;
  int maxIndex=(wOffset[MyContext+1]/wtot-DeltaStep)*static_cast<RealType>(TotalWalkers)+1;
  int nb=maxIndex-minIndex+1;
  vector<RealType> Zeta(nb);

  for(int i=minIndex, ii=0; i<maxIndex; i++,ii++) {
    Zeta[ii]= wtot*(DeltaStep+static_cast<RealType>(i)*nwInv);
  }

  RealType wCur=wOffset[MyContext];
  int ind=0;
  while(Zeta[ind]<wCur) {ind++;} 

  //surviving walkers
  int icdiff=0;
  for(iw=0; iw<nw; iw++) {
    RealType tryp=wCur+fabs(wConf[iw]);
    int ni=0;
    while(Zeta[ind]<tryp && Zeta[ind] >= wCur) {
      ind++;
      ni++;
    }
    wCur+=fabs(wConf[iw]);
    if(ni) {
      icdiff++;
    } 
    ncopy_w[iw]=ni;
  }

  vector<int> plus, minus;
  for(iw=0;iw<nw; iw++) {
    int m=ncopy_w[iw];
    if(m>1) {// add the index of this walker to plus, duplicate m-1 times
      plus.insert(plus.end(),m-1,iw);
    } else if(m==0) { // add the walker index to be killed/overwritten
      minus.push_back(iw);
    }
  }

  //copy within the local node
  int lower=std::min(plus.size(),minus.size()); 
  while(lower>0) {
    --lower;
    W[minus[lower]]->assign(*(W[plus[lower]]));
    minus.pop_back();
    plus.pop_back();
  }


  //dN[ip] extra/missing walkers
  vector<int> dN(NumContexts,0);
  if(plus.size()) { dN[MyContext]=plus.size();}
  if(minus.size()) { dN[MyContext]=-minus.size();}

  //collect the data
  gsum(dN,0);

  vector<int> minusN, plusN;
  for(int ip=0; ip<NumContexts; ip++) {
    if(dN[ip]>0) 
      plusN.insert(plusN.end(),dN[ip],ip);
    else if(dN[ip]<0)
      minusN.insert(minusN.end(),-dN[ip],ip);
  }

  //fout << " | p="<<plusN.size() << " ";
  //std::copy(plusN.begin(), plusN.end(), ostream_iterator<int>(fout," "));
  //fout << " | m="<< minusN.size() << " " ;
  //std::copy(minusN.begin(), minusN.end(), ostream_iterator<int>(fout," "));

  int wbuffer_size=W[0]->byteSize();

  //minusN.size() == plusN.size()
  //a node either send walker or recv walker but never both
  int nswap=plusN.size();
  int last = abs(dN[MyContext])-1;
  int ic=0;
  while(ic<nswap && last>=0) {
    if(plusN[ic]==MyContext) {
      OOMPI_Packed sendBuffer(wbuffer_size,OOMPI_COMM_WORLD);
      W[plus[last]]->putMessage(sendBuffer);
      OOMPI_COMM_WORLD[minusN[ic]].Send(sendBuffer);
      --last; 
    } 
    if(minusN[ic]==MyContext) {
      OOMPI_Packed recvBuffer(wbuffer_size,OOMPI_COMM_WORLD);
      OOMPI_COMM_WORLD[plusN[ic]].Recv(recvBuffer);
      W[minus[last]]->getMessage(recvBuffer);
      --last;
    }
    ++ic;
  }

  //collect surviving walkers
  return icdiff;
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

