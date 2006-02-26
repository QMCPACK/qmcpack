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
#include "QMCDrivers/DMC/WalkerReconfiguration.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
WalkerReconfiguration::WalkerReconfiguration() {
  SwapMode=1;
  NumContexts=OHMMS::Controller->ncontexts();
  MyContext=OHMMS::Controller->mycontext();

  UnitZeta=Random();

  //ofstream fout("check.dat");
}

int WalkerReconfiguration::getIndexPermutation(MCWalkerConfiguration& W) {

  if(Zeta.empty()) {
    Zeta.resize(W.getActiveWalkers()+1);
    IndexCopy.resize(W.getActiveWalkers());
    wConf.resize(W.getActiveWalkers());
  }

  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  RealType wtot=0.0;
  int nw=0;
  while(it != it_end) {
    wtot += wConf[nw]=(*it)->Weight;
    ++nw;++it;
  }

  RealType nwInv=1.0/static_cast<RealType>(nw);
  RealType dstep=UnitZeta*nwInv;

  for(int iw=0; iw<nw;iw++) {
    Zeta[iw]=wtot*(dstep+static_cast<RealType>(iw)*nwInv);
  }
  Zeta[nw]=wtot+1.0;

  //for(int iw=0; iw<nw; iw++) {
  //  fout << iw << " " << Zeta[iw+1]-Zeta[iw] << " " << wConf[iw] << endl;
  //}

  //assign negative
  //std::fill(IndexCopy.begin(),IndexCopy.end(),-1);

  int ind=0;
  RealType wCur=0.0;
  //surviving walkers
  int icdiff=0;
  it=W.begin();
  vector<int> ipip(nw,0);
  for(int iw=0; iw<nw; iw++) {
    RealType tryp=wCur+fabs(wConf[iw]);
    int ni=0;
    while(Zeta[ind]<tryp && Zeta[ind] >= wCur) {
      //IndexCopy[ind]=iw;
      ind++;
      ni++;
    }
    wCur+=fabs(wConf[iw]);
    if(ni) {
      icdiff++;
    } 
    ipip[iw]=ni;
  }

  //ofstream fout("check.dat", ios::app);
  //fout << wtot << " " << icdiff << endl;

  vector<int> plus,minus;
  for(int iw=0; iw<nw; iw++) {
    int m=ipip[iw];
    if(m>1) 
      plus.insert(plus.end(),m-1,iw);
    else if(m==0) 
      minus.push_back(iw);
  }

  for(int i=0; i<plus.size(); i++) {
    W[minus[i]]->assign(*(W[plus[i]]));
  }
  //int killed = shuffleIndex(nw);
  //fout << "# Total weight " << wtot << " " << killed <<  endl;
  //cout << "<<<< CopyIndex " << endl;
  //std::copy(IndexCopy.begin(), IndexCopy.end(), ostream_iterator<int>(cout, " "));
  //cout << endl << "<<<<<<" << endl;

  //for(int iw=0; iw<nw; iw++) {
  //  if(IndexCopy[iw] != iw) {
  //    W[iw]->assign(*(W[IndexCopy[iw]]));
  //  }
  //}

  return icdiff;
}

int WalkerReconfiguration::shuffleIndex(int nw) {
  vector<int> ipip(nw,0);
  for(int iw=0; iw<nw; iw++) ipip[IndexCopy[iw]]+=1;

  vector<int> indz;
  for(int iw=0; iw<nw; iw++) {
    if(ipip[iw]==0) {
      indz.push_back(iw);
    }
  }

  int ikilled=0;
  for(int iw=0; iw<nw; iw++) {
    if(ipip[iw] != 0) {
      IndexCopy[iw]=iw;
      for(int i=1;i<ipip[iw]; i++) {
        IndexCopy[indz[ikilled++]]=iw;
      }
    }
  }

  return indz.size();
}

int 
WalkerReconfiguration::branch(int iter, MCWalkerConfiguration& W, RealType trigger) {

  int nwkept = getIndexPermutation(W);

  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end) {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }

  return nwkept;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

