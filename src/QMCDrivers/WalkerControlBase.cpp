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
#include "QMCDrivers/WalkerControlBase.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{

WalkerControlBase::WalkerControlBase(Communicate* c, bool rn)
  : MPIObjectBase(c), SwapMode(0), Nmin(1), Nmax(10)
  , MaxCopy(2), NumWalkersCreated(0), NumWalkersSent(0)
  , targetEnergyBound(10), targetVar(2), targetSigma(10)
  , dmcStream(0), WriteRN(rn)
{
  NumContexts=myComm->size();
  MyContext=myComm->rank();
  curData.resize(LE_MAX+NumContexts);
  NumPerNode.resize(NumContexts);
  OffSet.resize(NumContexts+1);
  FairOffSet.resize(NumContexts+1);
  accumData.resize(LE_MAX);
}

WalkerControlBase::~WalkerControlBase()
{
  if(dmcStream)
    delete dmcStream;
}


//disable it: everything is done by a constructor
//void WalkerControlBase::setCommunicator(Communicate* c)
//{
//  initCommunicator(c);
//  NumContexts=myComm->size();
//  MyContext=myComm->rank();
//  curData.resize(LE_MAX+NumContexts);
//  NumPerNode.resize(NumContexts);
//  OffSet.resize(NumContexts+1);
//  FairOffSet.resize(NumContexts+1);
//  accumData.resize(LE_MAX);
//}

void WalkerControlBase::start()
{
  if(MyContext == 0)
  {
    string hname(myComm->getName());
    if (WriteRN)
      hname.append(".rn.dat");
    else
      hname.append(".dmc.dat");
    if(hname != dmcFname)
    {
      if(dmcStream)
        delete dmcStream;
      dmcStream= new ofstream(hname.c_str());
      //oa = new boost::archive::binary_oarchive (*dmcStream);
      dmcStream->setf(ios::scientific, ios::floatfield);
      dmcStream->precision(10);
      (*dmcStream) << "# Index "
                   << setw(20) << "LocalEnergy"
                   << setw(20) << "Variance"
                   << setw(20) << "Weight"
                   << setw(20) << "NumOfWalkers"
                   << setw(20) << "AvgSentWalkers";  //add the number of walkers
      if (WriteRN)
      {
        (*dmcStream) << setw(20) << "RNWalkers"
                     << setw(20) << "AlternateEnergy";
      }
      (*dmcStream)   << setw(20) << "TrialEnergy"
                     << setw(20) << "DiffEff";
//         if (WriteRN)
      (*dmcStream)  << setw(20) << "LivingFraction";
      (*dmcStream) << endl;
      dmcFname=hname;
    }
  }
}

void WalkerControlBase::setWalkerID(MCWalkerConfiguration& walkers)
{
  start(); //do the normal start
  MCWalkerConfiguration::iterator wit(walkers.begin());
  MCWalkerConfiguration::iterator wit_end(walkers.end());
  for(; wit != wit_end; ++wit)
  {
    if((*wit)->ID==0)
    {
      (*wit)->ID=(++NumWalkersCreated)*NumContexts+MyContext;
      (*wit)->ParentID=(*wit)->ID;
    }
  }
}

void WalkerControlBase::measureProperties(int iter)
{
  //taking average over the walkers
  RealType wgtInv(1.0/curData[WEIGHT_INDEX]);
  RealType eavg=curData[ENERGY_INDEX]*wgtInv;
  EnsembleProperty.Energy=eavg;
  EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  EnsembleProperty.Variance=(curData[ENERGY_SQ_INDEX]*wgtInv-eavg*eavg);
  EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  EnsembleProperty.R2Accepted=curData[R2ACCEPTED_INDEX];
  EnsembleProperty.R2Proposed=curData[R2PROPOSED_INDEX];
  EnsembleProperty.LivingFraction= static_cast<RealType>(curData[FNSIZE_INDEX]) / static_cast<RealType>(curData[FNSIZE_INDEX]+curData[RNONESIZE_INDEX]);
  EnsembleProperty.AlternateEnergy=curData[B_ENERGY_INDEX]/curData[B_WGT_INDEX];
  EnsembleProperty.RNSamples=curData[RNSIZE_INDEX];
  if(dmcStream)
  {
    //boost::archive::text_oarchive oa(*dmcStream);
    //(*oa) & iter  & eavg_cur & wgt_cur & Etrial  & pop_old;
    (*dmcStream) << setw(10) << iter
                 << setw(20) << EnsembleProperty.Energy
                 << setw(20) << EnsembleProperty.Variance
                 << setw(20) << EnsembleProperty.Weight
                 << setw(20) << EnsembleProperty.NumSamples
                 << setw(20) << curData[SENTWALKERS_INDEX]/static_cast<double>(NumContexts);
    if (WriteRN)
    {
      (*dmcStream) << setw(20) << EnsembleProperty.RNSamples
                   << setw(20) << EnsembleProperty.AlternateEnergy;
    }
    (*dmcStream)
        << setw(20) << trialEnergy
        << setw(20) << EnsembleProperty.R2Accepted/EnsembleProperty.R2Proposed;
//       if (WriteRN) (*dmcStream)
    (*dmcStream) << setw(20) << EnsembleProperty.LivingFraction;
    (*dmcStream)  << endl;
  }
}

void WalkerControlBase::reset()
{
  std::fill(accumData.begin(),accumData.end(),0.0);
  std::fill(curData.begin(),curData.end(),0.0);
}

int WalkerControlBase::doNotBranch(int iter, MCWalkerConfiguration& W)
{
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  RealType esum=0.0,e2sum=0.0,wsum=0.0,ecum=0.0, w2sum=0.0, besum=0.0, bwgtsum=0.0;
  RealType r2_accepted=0.0,r2_proposed=0.0;
  int nrn(0),ncr(0),nfn(0),ngoodfn(0), nc(0);
  for(; it!=it_end; ++it)
  {
    bool inFN=(((*it)->ReleasedNodeAge)==0);
    nc= std::min(static_cast<int>((*it)->Multiplicity),MaxCopy);
    if(WriteRN)
    {
      if ((*it)->ReleasedNodeAge==1)
        ncr+=1;
      else
        if ((*it)->ReleasedNodeAge==0)
        {
          nfn+=1;
          ngoodfn+=nc;
        }
      r2_accepted+=(*it)->Properties(R2ACCEPTED);
      r2_proposed+=(*it)->Properties(R2PROPOSED);
      RealType e((*it)->Properties(LOCALENERGY));
      RealType bfe((*it)->Properties(ALTERNATEENERGY));
      RealType wgt=((*it)->Weight);
      RealType rnwgt=((*it)->ReleasedNodeWeight);
      esum += wgt*rnwgt*e;
      e2sum += wgt*rnwgt*e*e;
      wsum += rnwgt*wgt;
      w2sum += rnwgt*rnwgt*wgt*wgt;
      ecum += e;
      besum += bfe*wgt;
      bwgtsum += wgt;
    }
    else
    {
      if (nc>0)
        nfn++;
      else
        ncr++;
      r2_accepted+=(*it)->Properties(R2ACCEPTED);
      r2_proposed+=(*it)->Properties(R2PROPOSED);
      RealType e((*it)->Properties(LOCALENERGY));
      RealType wgt=((*it)->Weight);
      esum += wgt*e;
      e2sum += wgt*e*e;
      wsum += wgt;
      w2sum += wgt*wgt;
      ecum += e;
    }
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(),curData.end(),0);
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=W.getActiveWalkers();
  curData[WEIGHT_INDEX]=wsum;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;
  curData[FNSIZE_INDEX]=static_cast<RealType>(nfn);
  curData[RNONESIZE_INDEX]=static_cast<RealType>(ncr);
  curData[RNSIZE_INDEX]=nrn;
  curData[B_ENERGY_INDEX]=besum;
  curData[B_WGT_INDEX]=bwgtsum;
  myComm->allreduce(curData);
  measureProperties(iter);
  trialEnergy=EnsembleProperty.Energy;
  W.EnsembleProperty=EnsembleProperty;
  return W.getActiveWalkers();
}

int WalkerControlBase::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  int prefill_numwalkers = sortWalkers(W);
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  //un-biased variance but we use the saimple one
  //W.EnsembleProperty.Variance=(e2sum*wsum-esum*esum)/(wsum*wsum-w2sum);
  ////add to the accumData for block average: REMOVE THIS
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];
  int nw_tot = copyWalkers(W);
  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
  //set the global number of walkers
  W.setGlobalNumWalkers(nw_tot);
  return nw_tot;
}

void Write2XYZ(MCWalkerConfiguration& W)
{
  ofstream fout("bad.xyz");
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  int nptcls(W.getTotalNum());
  while(it != it_end)
  {
    fout << nptcls << endl
         << "# E = " << (*it)->Properties(LOCALENERGY)
         << " Wgt= " << (*it)->Weight << endl;
    for(int i=0; i<nptcls; i++)
      fout << "H " << (*it)->R[i] << endl;
    ++it;
  }
}


/** evaluate curData and mark the bad/good walkers
 */
int WalkerControlBase::sortWalkers(MCWalkerConfiguration& W)
{
  MCWalkerConfiguration::iterator it(W.begin());
  vector<Walker_t*> bad,good_rn;
  vector<int> ncopy_rn;
  NumWalkers=0;
  MCWalkerConfiguration::iterator it_end(W.end());
  RealType esum=0.0,e2sum=0.0,wsum=0.0,ecum=0.0, w2sum=0.0, besum=0.0, bwgtsum=0.0;
  RealType r2_accepted=0.0,r2_proposed=0.0;
  int nfn(0),nrn(0),ngoodfn(0),ncr(0),nc(0);
  while(it != it_end)
  {
    bool inFN=(((*it)->ReleasedNodeAge)==0);
    nc= std::min(static_cast<int>((*it)->Multiplicity),MaxCopy);
    if(WriteRN)
    {
      if ((*it)->ReleasedNodeAge==1)
        ncr+=1;
      else
        if ((*it)->ReleasedNodeAge==0)
        {
          nfn+=1;
          ngoodfn+=nc;
        }
      r2_accepted+=(*it)->Properties(R2ACCEPTED);
      r2_proposed+=(*it)->Properties(R2PROPOSED);
      RealType e((*it)->Properties(LOCALENERGY));
      RealType bfe((*it)->Properties(ALTERNATEENERGY));
      RealType wgt=((*it)->Weight);
      RealType rnwgt=((*it)->ReleasedNodeWeight);
      esum += wgt*rnwgt*e;
      e2sum += wgt*rnwgt*e*e;
      wsum += rnwgt*wgt;
      w2sum += rnwgt*rnwgt*wgt*wgt;
      ecum += e;
      besum += bfe*wgt;
      bwgtsum += wgt;
    }
    else
    {
      if (nc>0)
        nfn++;
      else
        ncr++;
      r2_accepted+=(*it)->Properties(R2ACCEPTED);
      r2_proposed+=(*it)->Properties(R2PROPOSED);
      RealType e((*it)->Properties(LOCALENERGY));
      RealType wgt=((*it)->Weight);
      esum += wgt*e;
      e2sum += wgt*e*e;
      wsum += wgt;
      w2sum += wgt*wgt;
      ecum += e;
    }
    if((nc) && (inFN))
    {
      NumWalkers += nc;
      good_w.push_back(*it);
      ncopy_w.push_back(nc-1);
    }
    else
      if (nc)
      {
        NumWalkers += nc;
        nrn+=nc;
        good_rn.push_back(*it);
        ncopy_rn.push_back(nc-1);
      }
      else
      {
        bad.push_back(*it);
      }
    ++it;
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(),curData.end(),0);
  //update curData
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=W.getActiveWalkers();
  curData[WEIGHT_INDEX]=wsum;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;
  curData[FNSIZE_INDEX]=static_cast<RealType>(good_w.size());
  curData[RNONESIZE_INDEX]=static_cast<RealType>(ncr);
  curData[RNSIZE_INDEX]=nrn;
  curData[B_ENERGY_INDEX]=besum;
  curData[B_WGT_INDEX]=bwgtsum;
  ////this should be move
  //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  //W.EnsembleProperty.Energy=(esum/=wsum);
  //W.EnsembleProperty.Variance=(e2sum/wsum-esum*esum);
  //W.EnsembleProperty.Variance=(e2sum*wsum-esum*esum)/(wsum*wsum-w2sum);
  //remove bad walkers empty the container
  for(int i=0; i<bad.size(); i++)
    delete bad[i];
  if (!WriteRN)
  {
    if(good_w.empty())
    {
      app_error() << "All the walkers have died. Abort. " << endl;
      APP_ABORT("WalkerControlBase::sortWalkers");
    }
    int sizeofgood = good_w.size();
    //check if the projected number of walkers is too small or too large
    if(NumWalkers>Nmax)
    {
      int nsub=0;
      int nsub_target=(NumWalkers-nrn)-static_cast<int>(0.9*Nmax);
      int i=0;
      while(i< sizeofgood && nsub<nsub_target)
      {
        if(ncopy_w[i])
        {
          ncopy_w[i]--;
          nsub++;
        }
        ++i;
      }
      NumWalkers -= nsub;
    }
    else
      if(NumWalkers < Nmin)
      {
        int nadd=0;
        int nadd_target = static_cast<int>(Nmin*1.1)-(NumWalkers-nrn);
        if(nadd_target> sizeofgood)
        {
          app_warning() << "The number of walkers is running low. Requested walkers "
                        << nadd_target << " good walkers = " << sizeofgood << endl;
        }
        int i=0;
        while(i< sizeofgood && nadd<nadd_target)
        {
          ncopy_w[i]++;
          ++nadd;
          ++i;
        }
        NumWalkers +=  nadd;
      }
  }
  else
  {
    it=good_rn.begin();
    it_end=good_rn.end();
    int indy(0);
    while(it!=it_end)
    {
      good_w.push_back(*it);
      ncopy_w.push_back(ncopy_rn[indy]);
      it++,indy++;
    }
  }
  return NumWalkers;
}

int WalkerControlBase::copyWalkers(MCWalkerConfiguration& W)
{
  //clear the WalkerList to populate them with the good walkers
  W.clear();
  W.insert(W.begin(), good_w.begin(), good_w.end());
  int cur_walker = good_w.size();
  for(int i=0; i<good_w.size(); i++)
    //,ie+=ncols) {
  {
    for(int j=0; j<ncopy_w[i]; j++, cur_walker++)
    {
      Walker_t* awalker=new Walker_t(*(good_w[i]));
      awalker->ID=(++NumWalkersCreated)*NumContexts+MyContext;
      awalker->ParentID=good_w[i]->ParentID;
      W.push_back(awalker);
    }
  }
  //clear good_w and ncopy_w for the next branch
  good_w.clear();
  ncopy_w.clear();
  return W.getActiveWalkers();
}

bool WalkerControlBase::put(xmlNodePtr cur)
{
  ParameterSet params;
  params.add(targetEnergyBound,"energyBound","double");
  params.add(targetSigma,"sigmaBound","double");
  params.add(MaxCopy,"maxCopy","int");
  bool success=params.put(cur);
  app_log() << "  WalkerControlBase parameters " << endl;
  //app_log() << "    energyBound = " << targetEnergyBound << endl;
  //app_log() << "    sigmaBound = " << targetSigma << endl;
  app_log() << "    maxCopy = " << MaxCopy << endl;
  return true;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

