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

namespace qmcplusplus {

  WalkerControlBase::WalkerControlBase(Communicate* c): 
  SwapMode(0), Nmin(1), Nmax(10), MaxCopy(5), 
  targetEnergyBound(20), targetVar(2), targetSigma(10), dmcStream(0)
  {
    setCommunicator(c);

  }

  WalkerControlBase::~WalkerControlBase()
  {
    if(dmcStream) delete dmcStream;
  }

 
  void WalkerControlBase::setCommunicator(Communicate* c)
  {
    if(c) 
      myComm=c;
    else
      myComm = OHMMS::Controller;

    NumContexts=myComm->ncontexts();
    MyContext=myComm->mycontext();

    curData.resize(LE_MAX+NumContexts);
    NumPerNode.resize(NumContexts);
    OffSet.resize(NumContexts+1);
    FairOffSet.resize(NumContexts+1);
    accumData.resize(LE_MAX);
  }

  void WalkerControlBase::start()
  {
    if(MyContext == 0)
    {
      if(dmcStream) delete dmcStream;
      string hname(myComm->getName());
      hname.append(".dmc.dat");
      dmcStream= new ofstream(hname.c_str());
      //oa = new boost::archive::binary_oarchive (*dmcStream);
      dmcStream->setf(ios::scientific, ios::floatfield);
      dmcStream->precision(10);
      (*dmcStream) << setw(10) << "# Index "
        << setw(14) << "LocalEnergy"
        << setw(20) << "Variance" 
        << setw(22) << "Weight "
        << setw(20) << "NumOfWalkers" 
        << setw(20) << "TrialEnergy " 
        << endl;
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

    if(dmcStream) 
    {
      //boost::archive::text_oarchive oa(*dmcStream);
      //(*oa) & iter  & eavg_cur & wgt_cur & Etrial  & pop_old;
      (*dmcStream) << setw(10) << iter 
        << setw(20) << EnsembleProperty.Energy
        << setw(20) << EnsembleProperty.Variance
        << setw(20) << EnsembleProperty.Weight
        << setw(20) << EnsembleProperty.NumSamples
        << setw(20) << trialEnergy 
        << endl;
    }
  }

  void WalkerControlBase::reset() 
  {
    std::fill(accumData.begin(),accumData.end(),0.0);
  }

  int WalkerControlBase::branch(int iter, MCWalkerConfiguration& W, RealType trigger) {

    sortWalkers(W);

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
    while(it != it_end) {
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
    while(it != it_end) {
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
  void WalkerControlBase::sortWalkers(MCWalkerConfiguration& W) {

    MCWalkerConfiguration::iterator it(W.begin());

    vector<Walker_t*> bad;
    NumWalkers=0;
    MCWalkerConfiguration::iterator it_end(W.end());
    RealType esum=0.0,e2sum=0.0,wsum=0.0,ecum=0.0, w2sum=0.0;
    //RealType sigma=std::max(5.0*targetVar,targetEnergyBound);
    //RealType ebar= targetAvg;
    while(it != it_end) 
    {
      RealType e((*it)->Properties(LOCALENERGY));
      int nc= std::min(static_cast<int>((*it)->Multiplicity),MaxCopy);
      RealType wgt((*it)->Weight);

      esum += wgt*e;
      e2sum += wgt*e*e;
      wsum += wgt;
      w2sum += wgt*wgt;
      ecum += e;

      if(nc) 
      {
        NumWalkers += nc;
        good_w.push_back(*it);
        ncopy_w.push_back(nc-1);
      } else {
        bad.push_back(*it);
      }
      ++it;
    }

    //temp is an array to perform reduction operations
    std::fill(curData.begin(),curData.end(),0);

    //evaluate variance of this block
    //curVar=(e2sum-esum*esum/wsum)/wsum;
    //THIS IS NOT USED ANYMORE:BELOW
    //if(curVar>sigma) {
    //  app_error() << "Unphysical block variance is detected. Stop simulations." << endl;
    //  Write2XYZ(W);
    //  //Does not work some reason
    //  //OHMMS::Controller->abort();
    //#if defined(HAVE_MPI)
    //      OOMPI_COMM_WORLD.Abort();
    //#else
    //      abort();
    //#endif
    //    }
    //THIS IS NOT USED ANYMORE:ABOVE

    //update curData
    curData[ENERGY_INDEX]=esum;
    curData[ENERGY_SQ_INDEX]=e2sum;
    curData[WALKERSIZE_INDEX]=W.getActiveWalkers();
    curData[WEIGHT_INDEX]=wsum;
    curData[EREF_INDEX]=ecum;
    
    ////this should be move
    //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
    //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
    //W.EnsembleProperty.Energy=(esum/=wsum);
    //W.EnsembleProperty.Variance=(e2sum/wsum-esum*esum);
    //W.EnsembleProperty.Variance=(e2sum*wsum-esum*esum)/(wsum*wsum-w2sum);

    //remove bad walkers empty the container
    for(int i=0; i<bad.size(); i++) delete bad[i];

    if(good_w.empty()) {
      app_error() << "All the walkers have died. Abort. " << endl;
      OHMMS::Controller->abort();
    }

    int sizeofgood = good_w.size();

    //check if the projected number of walkers is too small or too large
    if(NumWalkers>Nmax) {
      int nsub=0;
      int nsub_target=NumWalkers-static_cast<int>(0.9*Nmax);
      int i=0;
      while(i< sizeofgood && nsub<nsub_target) {
        if(ncopy_w[i]) {ncopy_w[i]--; nsub++;}
        ++i;
      }
      NumWalkers -= nsub;
    } else  if(NumWalkers < Nmin) {
      int nadd=0;
      int nadd_target = static_cast<int>(Nmin*1.1)-NumWalkers;
      if(nadd_target> sizeofgood) {
        app_warning() << "The number of walkers is running low. Requested walkers " 
          << nadd_target << " good walkers = " << sizeofgood << endl;
      }
      int i=0;
      while(i< sizeofgood && nadd<nadd_target) {
        ncopy_w[i]++; ++nadd;++i;
      }
      NumWalkers +=  nadd;
    }
  }

  int WalkerControlBase::copyWalkers(MCWalkerConfiguration& W) {
    //clear the WalkerList to populate them with the good walkers
    W.clear();
    W.insert(W.begin(), good_w.begin(), good_w.end());

    int cur_walker = good_w.size();
    for(int i=0; i<good_w.size(); i++) { //,ie+=ncols) {
      for(int j=0; j<ncopy_w[i]; j++, cur_walker++) {
        W.push_back(new Walker_t(*(good_w[i])));
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
    app_log() << "    energyBound = " << targetEnergyBound << endl;
    app_log() << "    sigmaBound = " << targetSigma << endl;
    app_log() << "    maxCopy = " << MaxCopy << endl;
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

