//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/WalkerReconfiguration.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
WalkerReconfiguration::WalkerReconfiguration(Communicate* c) :WalkerControlBase(c)
{
  SwapMode=1;
  //ofstream fout("check.dat");
}

int WalkerReconfiguration::getIndexPermutation(MCWalkerConfiguration& W)
{
  int nw(W.getActiveWalkers());
  if(Zeta.size()!=nw)
  {
    Zeta.resize(nw+1);
    IndexCopy.resize(nw);
    wConf.resize(nw);
  }
  //accumulate the energies
  RealType esum=0.0,e2sum=0.0,wtot=0.0,ecum=0.0;
  MCWalkerConfiguration::iterator it(W.begin());
  RealType r2_accepted=0.0,r2_proposed=0.0;
  for(int iw=0; iw<nw; iw++)
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    RealType wgt((*it)->Weight);
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    ecum += e;
    wtot += wConf[iw]=wgt;
    ++it;
  }
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;
  RealType nwInv=1.0/static_cast<RealType>(nw);
  UnitZeta=Random();
  RealType dstep=UnitZeta*nwInv;
  for(int iw=0; iw<nw; iw++)
  {
    Zeta[iw]=wtot*(dstep+static_cast<RealType>(iw)*nwInv);
  }
  Zeta[nw]=wtot+1.0;
  //for(int iw=0; iw<nw; iw++) {
  //  fout << iw << " " << Zeta[iw+1]-Zeta[iw] << " " << wConf[iw] << std::endl;
  //}
  //assign negative
  //std::fill(IndexCopy.begin(),IndexCopy.end(),-1);
  int ind=0;
  RealType wCur=0.0;
  //surviving walkers
  int icdiff=0;
  it=W.begin();
  std::vector<int> ipip(nw,0);
  for(int iw=0; iw<nw; iw++)
  {
    RealType tryp=wCur+std::abs(wConf[iw]);
    int ni=0;
    while(Zeta[ind]<tryp && Zeta[ind] >= wCur)
    {
      //IndexCopy[ind]=iw;
      ind++;
      ni++;
    }
    wCur+=std::abs(wConf[iw]);
    if(ni)
    {
      icdiff++;
    }
    ipip[iw]=ni;
  }
  //ofstream fout("check.dat", std::ios::app);
  //fout << wtot << " " << icdiff << std::endl;
  std::vector<int> plus,minus;
  for(int iw=0; iw<nw; iw++)
  {
    int m=ipip[iw];
    if(m>1)
      plus.insert(plus.end(),m-1,iw);
    else
      if(m==0)
        minus.push_back(iw);
  }
  curData[FNSIZE_INDEX]=nw-minus.size();
  curData[RNONESIZE_INDEX]=minus.size();
  for(int i=0; i<plus.size(); i++)
  {
    int im=minus[i],ip=plus[i];
    W[im]->makeCopy(*(W[ip]));
    W[im]->ParentID=W[ip]->ID;
    W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
  }
  //int killed = shuffleIndex(nw);
  //fout << "# Total weight " << wtot << " " << killed <<  std::endl;
  //cout << "<<<< CopyIndex " << std::endl;
  //std::copy(IndexCopy.begin(), IndexCopy.end(), std::ostream_iterator<int>(std::cout, " "));
  //cout << std::endl << "<<<<<<" << std::endl;
  //for(int iw=0; iw<nw; iw++) {
  //  if(IndexCopy[iw] != iw) {
  //    W[iw]->assign(*(W[IndexCopy[iw]]));
  //  }
  //}
  return icdiff;
}

int WalkerReconfiguration::shuffleIndex(int nw)
{
  std::vector<int> ipip(nw,0);
  for(int iw=0; iw<nw; iw++)
    ipip[IndexCopy[iw]]+=1;
  std::vector<int> indz;
  for(int iw=0; iw<nw; iw++)
  {
    if(ipip[iw]==0)
    {
      indz.push_back(iw);
    }
  }
  int ikilled=0;
  for(int iw=0; iw<nw; iw++)
  {
    if(ipip[iw] != 0)
    {
      IndexCopy[iw]=iw;
      for(int i=1; i<ipip[iw]; i++)
      {
        IndexCopy[indz[ikilled++]]=iw;
      }
    }
  }
  return indz.size();
}

int
WalkerReconfiguration::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  int nwkept = getIndexPermutation(W);
  //update EnsembleProperty
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
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
  //curData[WALKERSIZE_INDEX]=nwkept;
  return nwkept;
}


