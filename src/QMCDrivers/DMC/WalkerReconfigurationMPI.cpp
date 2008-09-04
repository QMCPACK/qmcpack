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
#include "QMCDrivers/DMC/WalkerReconfigurationMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
WalkerReconfigurationMPI::WalkerReconfigurationMPI(Communicate* c): 
WalkerControlBase(c), TotalWalkers(0) 
{
  SwapMode=1;
  UnitZeta=Random();
  myComm->bcast(UnitZeta);
  app_log() << "  First weight [0,1) for reconfiguration =" << UnitZeta << endl;
}

int 
WalkerReconfigurationMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger) {

  int nwkept = swapWalkers(W);

  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;

  //RealType wgtInv(1.0/curData[WEIGHT_INDEX]);
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += nwkept;
  ////accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];

  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end) 
  {
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
  if(TotalWalkers ==0) 
  {
    FirstWalker=nw*MyContext;
    LastWalker=FirstWalker+nw;
    TotalWalkers = nw*NumContexts;
    nwInv = 1.0/static_cast<RealType>(TotalWalkers);
    DeltaStep = UnitZeta*nwInv;

    ncopy_w.resize(nw);
    wConf.resize(nw);
    //wSum.resize(NumContexts);
    wOffset.resize(NumContexts+1);
    dN.resize(NumContexts+1);
  }

  //std::fill(wSum.begin(),wSum.end(),0.0);

  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  int iw=0;
  RealType esum=0.0,e2sum=0.0,wtot=0.0,ecum=0.0;
  RealType r2_accepted=0.0,r2_proposed=0.0;
  while(it != it_end) 
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    RealType wgt((*it)->Weight);
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    wtot += wgt;
    ecum += e;
    wConf[iw++]=wgt;
    ++it;
  }
  //wSum[MyContext]=wtot;
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;

  std::fill(curData.begin()+LE_MAX,curData.end(),0.0);
  curData[LE_MAX+MyContext]=wtot;

  //collect everything
  myComm->allreduce(curData);

  //update EnsembleProperty
  W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];

  //wOffset[ip] is the partial sum update to ip
  wOffset[0]=0;
  //for(int ip=0; ip<NumContexts; ip++) wOffset[ip+1]=wOffset[ip]+wSum[ip];
  for(int ip=0,jp=LE_MAX; ip<NumContexts; ip++,jp++) 
    wOffset[ip+1]=wOffset[ip]+curData[jp];

  wtot=wOffset[NumContexts]; //wtot is the total weight

  //find the lower and upper bound of index
  int minIndex=static_cast<int>((wOffset[MyContext]/wtot-DeltaStep)*static_cast<RealType>(TotalWalkers))-1;
  int maxIndex=static_cast<int>((wOffset[MyContext+1]/wtot-DeltaStep)*static_cast<RealType>(TotalWalkers))+1;
  int nb=maxIndex-minIndex+1;
  vector<RealType> Zeta(nb);

  for(int i=minIndex, ii=0; i<maxIndex; i++,ii++) 
  {
    Zeta[ii]= wtot*(DeltaStep+static_cast<RealType>(i)*nwInv);
  }

  RealType wCur=wOffset[MyContext];
  int ind=0;
  while(Zeta[ind]<wCur) {ind++;} 

  //surviving walkers
  int icdiff=0;
  for(iw=0; iw<nw; iw++) {
    RealType tryp=wCur+abs(wConf[iw]);
    int ni=0;
    while(Zeta[ind]<tryp && Zeta[ind] >= wCur) {
      ind++;
      ni++;
    }
    wCur+=abs(wConf[iw]);
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
  while(lower>0) 
  {
    --lower;
    int im=minus[lower], ip=plus[lower];
    W[im]->makeCopy(*(W[ip]));
    W[im]->ParentID=W[ip]->ID;
    W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
    minus.pop_back();
    plus.pop_back();
  }


  //dN[ip] extra/missing walkers
  //dN[NumContexts] contains the number of surviving walkers
  std::fill(dN.begin(),dN.end(),0);
  dN[NumContexts]=icdiff;
  if(plus.size()) { dN[MyContext]=plus.size();}
  if(minus.size()) { dN[MyContext]=-minus.size();}

  //collect the data
  myComm->allreduce(dN);

  if(plus.size()) sendWalkers(W,plus);
  if(minus.size()) recvWalkers(W,minus);

  //vector<int> minusN, plusN;
  //bool tosend=false, torecv=false;
  //for(int ip=0; ip<NumContexts; ip++) {
  //  if(dN[ip]>0) {
  //    plusN.insert(plusN.end(),dN[ip],ip);
  //    tosend=true;
  //  } else if(dN[ip]<0) {
  //    minusN.insert(minusN.end(),-dN[ip],ip);
  //    torecv=true;
  //  }
  //}

  //int wbuffer_size=W[0]->byteSize();

  //int nswap=plusN.size();
  //int last = abs(dN[MyContext])-1;
  //int ic=0;
  //while(ic<nswap && last>=0) {
  //  if(plusN[ic]==MyContext) {
  //    OOMPI_Packed sendBuffer(wbuffer_size,OOMPI_COMM_WORLD);
  //    W[plus[last]]->putMessage(sendBuffer);
  //    OOMPI_COMM_WORLD[minusN[ic]].Send(sendBuffer);
  //    --last; 
  //  } 
  //  if(minusN[ic]==MyContext) {
  //    OOMPI_Packed recvBuffer(wbuffer_size,OOMPI_COMM_WORLD);
  //    OOMPI_COMM_WORLD[plusN[ic]].Recv(recvBuffer);
  //    W[minus[last]]->getMessage(recvBuffer);
  //    --last;
  //  }
  //  ++ic;
  //}

  //collect surviving walkers
  return dN[NumContexts];
}

void WalkerReconfigurationMPI::sendWalkers(MCWalkerConfiguration& W,
    const vector<IndexType>& plus) {
  vector<int> minusN, plusN;
  for(int ip=0; ip<NumContexts; ip++) {
    if(dN[ip]>0) {
      plusN.insert(plusN.end(),dN[ip],ip);
    } else if(dN[ip]<0) {
      minusN.insert(minusN.end(),-dN[ip],ip);
    }
  }

  int wbuffer_size=W[0]->byteSize();

  int nswap=plusN.size();
  int last = abs(dN[MyContext])-1;
  int ic=0;
  while(ic<nswap && last>=0) {
    if(plusN[ic]==MyContext) {
      //OOMPI_Packed sendBuffer(wbuffer_size,OOMPI_COMM_WORLD);
      OOMPI_Packed sendBuffer(wbuffer_size,myComm->getComm());
      W[plus[last]]->putMessage(sendBuffer);
      //OOMPI_COMM_WORLD[minusN[ic]].Send(sendBuffer);
      myComm->getComm()[minusN[ic]].Send(sendBuffer);
      --last; 
    } 
    ++ic;
  }
}

void WalkerReconfigurationMPI::recvWalkers(MCWalkerConfiguration& W,
    const vector<IndexType>& minus) {
  vector<IndexType> minusN, plusN;
  for(int ip=0; ip<NumContexts; ip++) {
    if(dN[ip]>0) {
      plusN.insert(plusN.end(),dN[ip],ip);
    } else if(dN[ip]<0) {
      minusN.insert(minusN.end(),-dN[ip],ip);
    }
  }

  int wbuffer_size=W[0]->byteSize();

  int nswap=plusN.size();
  int last = abs(dN[MyContext])-1;
  int ic=0;
  while(ic<nswap && last>=0) {
    if(minusN[ic]==MyContext) {
      //OOMPI_Packed recvBuffer(wbuffer_size,OOMPI_COMM_WORLD);
      //OOMPI_COMM_WORLD[plusN[ic]].Recv(recvBuffer);
      OOMPI_Packed recvBuffer(wbuffer_size,myComm->getComm());
      myComm->getComm()[plusN[ic]].Recv(recvBuffer);
      int im=minus[last];
      W[im]->getMessage(recvBuffer);
      W[im]->ParentID=W[im]->ID;
      W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
      --last;
    }
    ++ic;
  }
}


/***************************************************************************
 * $RCSfile: WalkerReconfigurationMPI.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
