//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/WalkerControlBase.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/ParameterSet.h"
#include "type_traits/template_types.hpp"

#include <numeric>

namespace qmcplusplus
{
WalkerControlBase::WalkerControlBase(Communicate* c, bool rn)
    : MPIObjectBase(c),
      n_min_(1),
      n_max_(10),
      MaxCopy(2),
      target_sigma_(10),
      dmcStream(0),
      NumWalkersCreated(0),
      SwapMode(0),
      write_release_nodes_(rn)
{
  method_    = -1; //assign invalid method
  num_contexts_ = myComm->size();
  MyContext   = myComm->rank();
  curData.resize(LE_MAX + num_contexts_);
  NumPerNode.resize(num_contexts_);
  OffSet.resize(num_contexts_ + 1);
  FairOffSet.resize(num_contexts_ + 1);
  accumData.resize(LE_MAX);
}

WalkerControlBase::~WalkerControlBase()
{
  if (dmcStream)
    delete dmcStream;
}


//disable it: everything is done by a constructor
//void WalkerControlBase::setCommunicator(Communicate* c)
//{
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
  if (MyContext == 0)
  {
    std::string hname(myComm->getName());
    if (write_release_nodes_)
      hname.append(".rn.dat");
    else
      hname.append(".dmc.dat");
    if (hname != dmcFname)
    {
      if (dmcStream)
        delete dmcStream;
      dmcStream = new std::ofstream(hname.c_str());
      //oa = new boost::archive::binary_oarchive (*dmcStream);
      dmcStream->setf(std::ios::scientific, std::ios::floatfield);
      dmcStream->precision(10);
      (*dmcStream) << "# Index " << std::setw(20) << "LocalEnergy" << std::setw(20) << "Variance" << std::setw(20)
                   << "Weight" << std::setw(20) << "NumOfWalkers" << std::setw(20)
                   << "AvgSentWalkers"; //add the number of walkers
      if (write_release_nodes_)
      {
        (*dmcStream) << std::setw(20) << "RNWalkers" << std::setw(20) << "AlternateEnergy";
      }
      (*dmcStream) << std::setw(20) << "TrialEnergy" << std::setw(20) << "DiffEff";
      //         if (WriteRN)
      (*dmcStream) << std::setw(20) << "LivingFraction";
      (*dmcStream) << std::endl;
      dmcFname = hname;
    }
  }
}

void WalkerControlBase::setWalkerID(MCWalkerConfiguration& walkers)
{
  start(); //do the normal start
  MCWalkerConfiguration::iterator wit(walkers.begin());
  MCWalkerConfiguration::iterator wit_end(walkers.end());
  for (; wit != wit_end; ++wit)
  {
    if ((*wit)->ID == 0)
    {
      (*wit)->ID       = (++NumWalkersCreated) * num_contexts_ + MyContext;
      (*wit)->ParentID = (*wit)->ID;
    }
  }
}

void WalkerControlBase::setWalkerID(MCPopulation& population)
{
  start(); //do the normal start
}

/** Depends on alot of state
 *
 *  Does not depend on state refactored to  PopulationAdjustment directly
 */
void WalkerControlBase::measureProperties(int iter)
{
  //taking average over the walkers
  FullPrecRealType wgtInv(1.0 / curData[WEIGHT_INDEX]);
  FullPrecRealType eavg                   = curData[ENERGY_INDEX] * wgtInv;
  ensemble_property_.Energy         = eavg;
  ensemble_property_.Weight         = curData[WEIGHT_INDEX];
  ensemble_property_.Variance       = (curData[ENERGY_SQ_INDEX] * wgtInv - eavg * eavg);
  ensemble_property_.NumSamples     = curData[WALKERSIZE_INDEX];
  ensemble_property_.R2Accepted     = curData[R2ACCEPTED_INDEX];
  ensemble_property_.R2Proposed     = curData[R2PROPOSED_INDEX];
  ensemble_property_.LivingFraction = static_cast<FullPrecRealType>(curData[FNSIZE_INDEX]) /
      static_cast<FullPrecRealType>(curData[FNSIZE_INDEX] + curData[RNONESIZE_INDEX]);
  ensemble_property_.AlternateEnergy = curData[B_ENERGY_INDEX] / curData[B_WGT_INDEX];
  ensemble_property_.RNSamples       = curData[RNSIZE_INDEX];
  if (dmcStream)
  {
    //boost::archive::text_oarchive oa(*dmcStream);
    //(*oa) & iter  & eavg_cur & wgt_cur & Etrial  & pop_old;
    (*dmcStream) << std::setw(10) << iter << std::setw(20) << ensemble_property_.Energy << std::setw(20)
                 << ensemble_property_.Variance << std::setw(20) << ensemble_property_.Weight << std::setw(20)
                 << ensemble_property_.NumSamples << std::setw(20)
                 << curData[SENTWALKERS_INDEX] / static_cast<double>(num_contexts_);
    if (write_release_nodes_)
    {
      (*dmcStream) << std::setw(20) << ensemble_property_.RNSamples << std::setw(20) << ensemble_property_.AlternateEnergy;
    }
    (*dmcStream) << std::setw(20) << trialEnergy << std::setw(20)
                 << ensemble_property_.R2Accepted / ensemble_property_.R2Proposed;
    //       if (WriteRN) (*dmcStream)
    (*dmcStream) << std::setw(20) << ensemble_property_.LivingFraction;
    (*dmcStream) << std::endl;
  }
}

void WalkerControlBase::reset()
{
  std::fill(accumData.begin(), accumData.end(), 0.0);
  std::fill(curData.begin(), curData.end(), 0.0);
}

int WalkerControlBase::doNotBranch(int iter, MCWalkerConfiguration& W)
{
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0, ecum = 0.0, w2sum = 0.0, besum = 0.0, bwgtsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int nrn(0), ncr(0), nfn(0), ngoodfn(0), nc(0);
  for (; it != it_end; ++it)
  {
    bool inFN = (((*it)->ReleasedNodeAge) == 0);
    nc        = std::min(static_cast<int>((*it)->Multiplicity), MaxCopy);
    if (write_release_nodes_)
    {
      if ((*it)->ReleasedNodeAge == 1)
        ncr += 1;
      else if ((*it)->ReleasedNodeAge == 0)
      {
        nfn += 1;
        ngoodfn += nc;
      }
      r2_accepted += (*it)->Properties(R2ACCEPTED);
      r2_proposed += (*it)->Properties(R2PROPOSED);
      FullPrecRealType e((*it)->Properties(LOCALENERGY));
      FullPrecRealType bfe((*it)->Properties(ALTERNATEENERGY));
      FullPrecRealType wgt   = ((*it)->Weight);
      FullPrecRealType rnwgt = ((*it)->ReleasedNodeWeight);
      esum += wgt * rnwgt * e;
      e2sum += wgt * rnwgt * e * e;
      wsum += rnwgt * wgt;
      w2sum += rnwgt * rnwgt * wgt * wgt;
      ecum += e;
      besum += bfe * wgt;
      bwgtsum += wgt;
    }
    else
    {
      if (nc > 0)
        nfn++;
      else
        ncr++;
      r2_accepted += (*it)->Properties(R2ACCEPTED);
      r2_proposed += (*it)->Properties(R2PROPOSED);
      FullPrecRealType e((*it)->Properties(LOCALENERGY));
      // This is a trick to estimate the number of walkers
      // after the first iterration branching.
      //RealType wgt=((*it)->Weight);
      FullPrecRealType wgt = FullPrecRealType(nc);
      esum += wgt * e;
      e2sum += wgt * e * e;
      wsum += wgt;
      w2sum += wgt * wgt;
      ecum += e;
    }
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(), curData.end(), 0);
  curData[ENERGY_INDEX]      = esum;
  curData[ENERGY_SQ_INDEX]   = e2sum;
  curData[WALKERSIZE_INDEX]  = W.getActiveWalkers();
  curData[WEIGHT_INDEX]      = wsum;
  curData[EREF_INDEX]        = ecum;
  curData[R2ACCEPTED_INDEX]  = r2_accepted;
  curData[R2PROPOSED_INDEX]  = r2_proposed;
  curData[FNSIZE_INDEX]      = static_cast<FullPrecRealType>(nfn);
  curData[RNONESIZE_INDEX]   = static_cast<FullPrecRealType>(ncr);
  curData[RNSIZE_INDEX]      = nrn;
  curData[B_ENERGY_INDEX]    = besum;
  curData[B_WGT_INDEX]       = bwgtsum;

  myComm->allreduce(curData);
  measureProperties(iter);
  trialEnergy        = ensemble_property_.Energy;
  W.EnsembleProperty = ensemble_property_;
  //return W.getActiveWalkers();
  return int(curData[WEIGHT_INDEX]);
}

int WalkerControlBase::doNotBranch(int iter, MCPopulation& pop)
{
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(pop.get_walkers()));
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0, ecum = 0.0, w2sum = 0.0, besum = 0.0, bwgtsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int nrn(0), ncr(0), nfn(0), ngoodfn(0), nc(0);
  for(MCPWalker& walker : walkers)
  {
    bool inFN = ((walker.ReleasedNodeAge) == 0);
    nc        = std::min(static_cast<int>(walker.Multiplicity), MaxCopy);
    if (write_release_nodes_)
    {
      if (walker.ReleasedNodeAge == 1)
        ncr += 1;
      else if (walker.ReleasedNodeAge == 0)
      {
        nfn += 1;
        ngoodfn += nc;
      }
      r2_accepted += walker.Properties(R2ACCEPTED);
      r2_proposed += walker.Properties(R2PROPOSED);
      FullPrecRealType e(walker.Properties(LOCALENERGY));
      FullPrecRealType bfe(walker.Properties(ALTERNATEENERGY));
      FullPrecRealType wgt   = (walker.Weight);
      FullPrecRealType rnwgt = (walker.ReleasedNodeWeight);
      esum += wgt * rnwgt * e;
      e2sum += wgt * rnwgt * e * e;
      wsum += rnwgt * wgt;
      w2sum += rnwgt * rnwgt * wgt * wgt;
      ecum += e;
      besum += bfe * wgt;
      bwgtsum += wgt;
    }
    else
    {
      if (nc > 0)
        nfn++;
      else
        ncr++;
      r2_accepted += walker.Properties(R2ACCEPTED);
      r2_proposed += walker.Properties(R2PROPOSED);
      FullPrecRealType e(walker.Properties(LOCALENERGY));
      // This is a trick to estimate the number of walkers
      // after the first iterration branching.
      //RealType wgt=(walker.Weight);
      FullPrecRealType wgt = FullPrecRealType(nc);
      esum += wgt * e;
      e2sum += wgt * e * e;
      wsum += wgt;
      w2sum += wgt * wgt;
      ecum += e;
    }
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(), curData.end(), 0);
  curData[ENERGY_INDEX]      = esum;
  curData[ENERGY_SQ_INDEX]   = e2sum;
  curData[WALKERSIZE_INDEX]  = pop.get_num_global_walkers();
  curData[WEIGHT_INDEX]      = wsum;
  curData[EREF_INDEX]        = ecum;
  curData[R2ACCEPTED_INDEX]  = r2_accepted;
  curData[R2PROPOSED_INDEX]  = r2_proposed;
  curData[FNSIZE_INDEX]      = static_cast<FullPrecRealType>(nfn);
  curData[RNONESIZE_INDEX]   = static_cast<FullPrecRealType>(ncr);
  curData[RNSIZE_INDEX]      = nrn;
  curData[B_ENERGY_INDEX]    = besum;
  curData[B_WGT_INDEX]       = bwgtsum;

  myComm->allreduce(curData);
  measureProperties(iter);
  trialEnergy        = ensemble_property_.Energy;
  pop.set_ensemble_property(ensemble_property_);
  //return W.getActiveWalkers();
  return int(curData[WEIGHT_INDEX]);
}

int WalkerControlBase::branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger)
{
  NumPerNode[0] = sortWalkers(W);
  measureProperties(iter);
  W.EnsembleProperty = ensemble_property_;
  //un-biased variance but we use the saimple one
  //W.EnsembleProperty.Variance=(e2sum*wsum-esum*esum)/(wsum*wsum-w2sum);
  ////add to the accumData for block average: REMOVE THIS
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];
  int current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);
  applyNmaxNmin(current_population);
  int nw_tot = copyWalkers(W);
  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  while (it != it_end)
  {
    (*it)->Weight       = 1.0;
    (*it)->Multiplicity = 1.0;
    ++it;
  }
  //set the global number of walkers
  W.setGlobalNumWalkers(nw_tot);
  // Update offsets in non-MPI case, needed to ensure checkpoint outputs the correct
  // number of configurations.
  if (W.WalkerOffsets.size() == 2)
  {
    W.WalkerOffsets[1] = nw_tot;
  }
  return nw_tot;
}

int WalkerControlBase::branch(int iter, MCPopulation& pop, FullPrecRealType trigger)
{
  PopulationAdjustment adjust(calcPopulationAdjustment(pop));
  // adjust.num_new_walkers now caries what was stashed in state variable NumPerNode[0]
  pop.set_num_global_walkers(adjustPopulation(pop, adjust));
  measureProperties(iter);
  // More ParticleSet is a set of particles or a collection of particles and sort of both.
  
  pop.set_ensemble_property(ensemble_property_);

  applyNmaxNmin(pop.get_num_global_walkers());

  while( ! adjust.bad_walkers.empty() )
  {
    pop.killWalker(adjust.bad_walkers.back());
    adjust.bad_walkers.pop_back();
  }
  for( int iw = 0; iw < adjust.good_walkers.size(); ++iw )
  {
    for(int i_copies = 0; i_copies < adjust.copies_to_make[iw]; ++i_copies )
    {
      auto& walker = pop.spawnWalker();
      walker = adjust.good_walkers[iw];
      // IF these are really unique ID's they should be UUID's or something
      // old algorithm seems to reuse them in a way that I'm not sure avoids
      // duplicates even at a particular time.
      walker.ID = pop.get_num_local_walkers();
      walker.ParentID = adjust.good_walkers[iw].get().ParentID;
      walker.Weight  = 1.0;
      walker.Multiplicity = 1.0;
    }
  }
  return pop.get_num_local_walkers();
}

// int nw_tot = copyWalkers(W);
  // //set Weight and Multiplicity to default values
  // MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  // while (it != it_end)
  // {
  //   (*it)->Weight       = 1.0;
  //   (*it)->Multiplicity = 1.0;
  //   ++it;
  // }
  // //set the global number of walkers
  // W.setGlobalNumWalkers(nw_tot);
  // // Update offsets in non-MPI case, needed to ensure checkpoint outputs the correct
  // // number of configurations.
  // if (W.WalkerOffsets.size() == 2)
  // {
  //   W.WalkerOffsets[1] = nw_tot;
  // }
  // return nw_tot;


void WalkerControlBase::Write2XYZ(MCWalkerConfiguration& W)
{
  std::ofstream fout("bad.xyz");
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  int nptcls(W.getTotalNum());
  while (it != it_end)
  {
    fout << nptcls << std::endl << "# E = " << (*it)->Properties(LOCALENERGY) << " Wgt= " << (*it)->Weight << std::endl;
    for (int i = 0; i < nptcls; i++)
      fout << "H " << (*it)->R[i] << std::endl;
    ++it;
  }
}


/** evaluate curData and mark the bad/good walkers.
 *
 *  Each good walker has a counter registering the
 *  number of copies needed to be generated from this walker.
 *  Bad walkers will be either recycled or removed later.
 */
int WalkerControlBase::sortWalkers(MCWalkerConfiguration& W)
{
  MCWalkerConfiguration::iterator it(W.begin());
  std::vector<Walker_t*> good_rn;
  std::vector<int> ncopy_rn;
  NumWalkers = 0;
  MCWalkerConfiguration::iterator it_end(W.end());
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0, ecum = 0.0, w2sum = 0.0, besum = 0.0, bwgtsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int nfn(0), nrn(0), ngoodfn(0), ncr(0), nc(0);
  while (it != it_end)
  {
    bool inFN = (((*it)->ReleasedNodeAge) == 0);
    nc        = std::min(static_cast<int>((*it)->Multiplicity), MaxCopy);
    if (write_release_nodes_)
    {
      if ((*it)->ReleasedNodeAge == 1)
        ncr += 1;
      else if ((*it)->ReleasedNodeAge == 0)
      {
        nfn += 1;
        ngoodfn += nc;
      }
      r2_accepted += (*it)->Properties(R2ACCEPTED);
      r2_proposed += (*it)->Properties(R2PROPOSED);
      FullPrecRealType local_energy((*it)->Properties(LOCALENERGY));
      FullPrecRealType alternate_energy((*it)->Properties(ALTERNATEENERGY));
      FullPrecRealType wgt   = ((*it)->Weight);
      FullPrecRealType rnwgt = ((*it)->ReleasedNodeWeight);
      esum += wgt * rnwgt * local_energy;
      e2sum += wgt * rnwgt * local_energy * local_energy;
      wsum += rnwgt * wgt;
      w2sum += rnwgt * rnwgt * wgt * wgt;
      ecum += local_energy;
      besum += alternate_energy* wgt;
      bwgtsum += wgt;
    }
    else
    {
      if (nc > 0)
        nfn++;
      else
        ncr++;
      r2_accepted += (*it)->Properties(R2ACCEPTED);
      r2_proposed += (*it)->Properties(R2PROPOSED);
      FullPrecRealType e((*it)->Properties(LOCALENERGY));
      FullPrecRealType wgt = ((*it)->Weight);
      esum += wgt * e;
      e2sum += wgt * e * e;
      wsum += wgt;
      w2sum += wgt * wgt;
      ecum += e;
    }
    if ((nc) && (inFN))
    {
      NumWalkers += nc;
      good_w.push_back(*it);
      ncopy_w.push_back(nc - 1);
    }
    else if (nc)
    {
      NumWalkers += nc;
      nrn += nc;
      good_rn.push_back(*it);
      ncopy_rn.push_back(nc - 1);
    }
    else
    {
      bad_w.push_back(*it);
    }
    ++it;
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(), curData.end(), 0);
  //update curData
  curData[ENERGY_INDEX]     = esum;
  curData[ENERGY_SQ_INDEX]  = e2sum;
  curData[WALKERSIZE_INDEX] = W.getActiveWalkers();
  curData[WEIGHT_INDEX]     = wsum;
  curData[EREF_INDEX]       = ecum;
  curData[R2ACCEPTED_INDEX] = r2_accepted;
  curData[R2PROPOSED_INDEX] = r2_proposed;
  // This is really an integral type but is implicitly converted
  curData[FNSIZE_INDEX]     = good_w.size();
  // This is really an integral type but is implicitly converted
  curData[RNONESIZE_INDEX]  = ncr;
  // This is really an integral type but is implicitly converted
  curData[RNSIZE_INDEX]     = nrn;
  curData[B_ENERGY_INDEX]   = besum;
  curData[B_WGT_INDEX]      = bwgtsum;
  ////this should be move
  //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  //W.EnsembleProperty.Energy=(esum/=wsum);
  //W.EnsembleProperty.Variance=(e2sum/wsum-esum*esum);
  //W.EnsembleProperty.Variance=(e2sum*wsum-esum*esum)/(wsum*wsum-w2sum);
  if (write_release_nodes_)
  {
    it     = good_rn.begin();
    it_end = good_rn.end();
    int indy(0);
    while (it != it_end)
    {
      good_w.push_back(*it);
      ncopy_w.push_back(ncopy_rn[indy]);
      it++, indy++;
    }
  }
  return NumWalkers;
}

/** Unified Driver version: evaluate curData and mark the bad/good walkers.
 *
 *  Somehow I don't think this belongs here
 *  
 *  Each good walker has a counter registering the
 *  number of copies needed to be generated from this walker.
 *  Bad walkers will be either recycled or removed later.
 */
WalkerControlBase::PopulationAdjustment WalkerControlBase::calcPopulationAdjustment(MCPopulation& pop)
{
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(pop.get_walkers()));
  RefVector<MCPWalker> good_walkers;
  std::vector<int> copies_to_make_good_walkers;

  NumPerNode[0] = pop.get_num_local_walkers();
  
  PopulationAdjustment adjustment;
  adjustment.num_walkers = 0;
  // So many sums and counters
  FullPrecRealType esum = 0.0, e2sum = 0.0, wsum = 0.0, ecum = 0.0, w2sum = 0.0, besum = 0.0, bwgtsum = 0.0;
  FullPrecRealType r2_accepted = 0.0, r2_proposed = 0.0;
  int nfn(0), nrn(0), ngoodfn(0), ncr(0);
  for(MCPWalker& walker : walkers)
  {
    bool inFN = (walker.ReleasedNodeAge == 0);

    // Written as a lambda for emphasis
    auto calcNumberWalkerCopies = [this](int multiplicity)->int {
            return std::min(multiplicity, MaxCopy);
                                  };
    int nc = calcNumberWalkerCopies(walker.Multiplicity);

    if (write_release_nodes_)
    {
      if (walker.ReleasedNodeAge == 1)
        ncr += 1;
      else if (walker.ReleasedNodeAge == 0)
      {
        nfn += 1;
        ngoodfn += nc;
      }
      r2_accepted += walker.Properties(R2ACCEPTED);
      r2_proposed += walker.Properties(R2PROPOSED);
      FullPrecRealType local_energy(walker.Properties(LOCALENERGY));
      FullPrecRealType alternate_energy(walker.Properties(ALTERNATEENERGY));
      FullPrecRealType wgt   = walker.Weight;
      FullPrecRealType rnwgt = walker.ReleasedNodeWeight;
      esum += wgt * rnwgt * local_energy;
      e2sum += wgt * rnwgt * local_energy * local_energy;
      wsum += rnwgt * wgt;
      w2sum += rnwgt * rnwgt * wgt * wgt;
      ecum += local_energy;
      besum += alternate_energy * wgt;
      bwgtsum += wgt;
    }
    else
    {
      if (nc > 0)
        nfn++;
      else
        ncr++;
      r2_accepted += walker.Properties(R2ACCEPTED);
      r2_proposed += walker.Properties(R2PROPOSED);
      FullPrecRealType local_energy(walker.Properties(LOCALENERGY));
      FullPrecRealType wgt = walker.Weight;
      esum += wgt * local_energy;
      e2sum += wgt * local_energy * local_energy;
      wsum += wgt;
      w2sum += wgt * wgt;
      ecum += local_energy;
    }
    if ((nc) && (inFN))
    {
      adjustment.num_walkers += nc;
      adjustment.good_walkers.push_back(walker);
      adjustment.copies_to_make.push_back(nc - 1);
    }
    else if (nc)
    {
      // This seems odd since if write_release_nodes_ == false
      // num_new_walkers != sum(copies_to_make)
      // qmcpack style would be to use this as an ersatz flag somewhere
      adjustment.num_walkers += nc;
      nrn += nc;
      // Nothing is every done with this except put its size in
      // curData[FNSIZE_INDEX]
      good_walkers.push_back(walker);
      copies_to_make_good_walkers.push_back(nc - 1);
    }
    else
    {
      adjustment.bad_walkers.push_back(walker);
    }
  }
  //temp is an array to perform reduction operations
  std::fill(curData.begin(), curData.end(), 0);
  //update curData
  curData[ENERGY_INDEX]     = esum;
  curData[ENERGY_SQ_INDEX]  = e2sum;
  curData[WALKERSIZE_INDEX] = pop.get_walkers().size();
  curData[WEIGHT_INDEX]     = wsum;
  curData[EREF_INDEX]       = ecum;
  curData[R2ACCEPTED_INDEX] = r2_accepted;
  curData[R2PROPOSED_INDEX] = r2_proposed;
  curData[FNSIZE_INDEX]     = good_walkers.size();
  curData[RNONESIZE_INDEX]  = ncr;
  curData[RNSIZE_INDEX]     = nrn;
  curData[B_ENERGY_INDEX]   = besum;
  curData[B_WGT_INDEX]      = bwgtsum;
  if (write_release_nodes_)
  {
    auto addWalkersWithReleaseNodeAge = [&adjustment](MCPWalker& walker, int copies) {
      // so now if release nodes is on we push all the good nodes to good_walkers_temp,
      // if inFN you'll have two copies of each.
      // I'm just going to preserve this logic but there has to be a simpler way to express
      // whatever the point of this is.
      adjustment.good_walkers.push_back(walker);
      adjustment.copies_to_make.push_back(copies);
                                        };                                          
    for(int iw = 0; iw < good_walkers.size(); ++iw)
      addWalkersWithReleaseNodeAge(good_walkers[iw],copies_to_make_good_walkers[iw]);
  }
  return adjustment;
}

/** Not yet compatible with MCPopulation
 */
int WalkerControlBase::applyNmaxNmin(int current_population)
{

  // limit Nmax
  int current_max = (current_population + num_contexts_ - 1) / num_contexts_;
  if (current_max > n_max_)
  {
    app_warning() << "Exceeding Max Walkers per MPI rank : " << n_max_ << ". Ceiling is applied" << std::endl;
    int nsub = current_population - n_max_ * num_contexts_;
    for (int inode = 0; inode < num_contexts_; inode++)
      if (NumPerNode[inode] > n_max_)
      {
        int n_remove = std::min(nsub, NumPerNode[inode] - n_max_);
        NumPerNode[inode] -= n_remove;
        nsub -= n_remove;

        if (inode == MyContext)
        {
          for (int iw = 0; iw < ncopy_w.size(); iw++)
          {
            int n_remove_walker = std::min(ncopy_w[iw], n_remove);
            ncopy_w[iw] -= n_remove_walker;
            n_remove -= n_remove_walker;
            if (n_remove == 0)
              break;
          }

          if (n_remove > 0)
          {
            app_warning() << "Removing copies of good walkers is not enough. "
                          << "Removing good walkers." << std::endl;
            do
            {
              bad_w.push_back(good_w.back());
              good_w.pop_back();
              ncopy_w.pop_back();
              --n_remove;
            } while (n_remove > 0 && !good_w.empty());
          }

          if (n_remove)
            APP_ABORT("WalkerControlBase::applyNmaxNmin not able to remove sufficient walkers on a node!");
          if (std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) != NumPerNode[inode])
            APP_ABORT("WalkerControlBase::applyNmaxNmin walker removal mismatch!");
        }

        if (nsub == 0)
          break;
      }

    if (nsub)
      APP_ABORT("WalkerControlBase::applyNmaxNmin not able to remove sufficient walkers overall!");
  }

  // limit Nmin
  if (current_population / num_contexts_ < n_min_)
  {
    app_warning() << "The number of walkers is running lower than Min Walkers per MPI rank : " << n_min_
                  << ". Floor is applied" << std::endl;
    int nadd = n_min_ * num_contexts_ - current_population;
    for (int inode = 0; inode < num_contexts_; inode++)
      if (NumPerNode[inode] > 0 && NumPerNode[inode] < n_min_)
      {
        int n_insert = std::min(nadd, n_min_ - NumPerNode[inode]);
        NumPerNode[inode] += n_insert;
        nadd -= n_insert;

        if (inode == MyContext)
        {
          int n_avg_insert_per_walker = (n_insert + ncopy_w.size() - 1) / ncopy_w.size();
          for (int iw = 0; iw < ncopy_w.size(); iw++)
          {
            int n_insert_walker = std::min(n_avg_insert_per_walker, n_insert);
            ncopy_w[iw] += n_insert_walker;
            n_insert -= n_insert_walker;
            if (n_insert == 0)
              break;
          }

          if (std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) != NumPerNode[inode])
            APP_ABORT("WalkerControlBase::applyNmaxNmin walker insertion mismatch!");
        }

        if (nadd == 0)
          break;
      }

    if (nadd)
      app_warning() << "WalkerControlBase::applyNmaxNmin not able to add sufficient walkers overall!" << std::endl;
  }

  // check current population
  current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);
  // at least one walker after load-balancing
  if (current_population / num_contexts_ == 0)
  {
    app_error() << "Some MPI ranks have no walkers after load balancing. This should not happen."
                << "Improve the trial wavefunction or adjust the simulation parameters." << std::endl;
    APP_ABORT("WalkerControlBase::sortWalkers");
  }

  return current_population;
}

/**  Not really sorting so its been rename something more sensible.
 *
 *   Not actually finished.  needs to use adjust's walker vectors
 */
int WalkerControlBase::adjustPopulation(MCPopulation& population, PopulationAdjustment& adjust)
{
  // In the unified driver design each ranks MCPopulation knows this and it is not stored a bunch of other places
  int current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);

  // limit Nmax
  // TODO:  this seems to be the wrong pace to do this.
  int current_max = (current_population + num_contexts_ - 1) / num_contexts_;
  if (current_max > n_max_)
  {
    app_warning() << "Exceeding Max Walkers per MPI rank : " << n_max_ << ". Ceiling is applied" << std::endl;
    int nsub = current_population - n_max_ * num_contexts_;
    for (int inode = 0; inode < num_contexts_; inode++)
      if (NumPerNode[inode] > n_max_)
      {
        int n_remove = std::min(nsub, NumPerNode[inode] - n_max_);
        NumPerNode[inode] -= n_remove;
        nsub -= n_remove;

        if (inode == MyContext)
        {
          for (int iw = 0; iw < ncopy_w.size(); iw++)
          {
            int n_remove_walker = std::min(ncopy_w[iw], n_remove);
            ncopy_w[iw] -= n_remove_walker;
            n_remove -= n_remove_walker;
            if (n_remove == 0)
              break;
          }

          if (n_remove > 0)
          {
            app_warning() << "Removing copies of good walkers is not enough. "
                          << "Removing good walkers." << std::endl;
            do
            {
              bad_w.push_back(good_w.back());
              good_w.pop_back();
              ncopy_w.pop_back();
              --n_remove;
            } while (n_remove > 0 && !good_w.empty());
          }

          if (n_remove)
            APP_ABORT("WalkerControlBase::applyNmaxNmin not able to remove sufficient walkers on a node!");
          if (std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) != NumPerNode[inode])
            APP_ABORT("WalkerControlBase::applyNmaxNmin walker removal mismatch!");
        }

        if (nsub == 0)
          break;
      }

    if (nsub)
      APP_ABORT("WalkerControlBase::applyNmaxNmin not able to remove sufficient walkers overall!");
  }

  // limit Nmin
  if (current_population / num_contexts_ < n_min_)
  {
    app_warning() << "The number of walkers is running lower than Min Walkers per MPI rank : " << n_min_
                  << ". Floor is applied" << std::endl;
    int nadd = n_min_ * num_contexts_ - current_population;
    for (int inode = 0; inode < num_contexts_; inode++)
      if (NumPerNode[inode] > 0 && NumPerNode[inode] < n_min_)
      {
        int n_insert = std::min(nadd, n_min_ - NumPerNode[inode]);
        NumPerNode[inode] += n_insert;
        nadd -= n_insert;

        if (inode == MyContext)
        {
          int n_avg_insert_per_walker = (n_insert + ncopy_w.size() - 1) / ncopy_w.size();
          for (int iw = 0; iw < ncopy_w.size(); iw++)
          {
            int n_insert_walker = std::min(n_avg_insert_per_walker, n_insert);
            ncopy_w[iw] += n_insert_walker;
            n_insert -= n_insert_walker;
            if (n_insert == 0)
              break;
          }

          if (std::accumulate(ncopy_w.begin(), ncopy_w.end(), ncopy_w.size()) != NumPerNode[inode])
            APP_ABORT("WalkerControlBase::applyNmaxNmin walker insertion mismatch!");
        }

        if (nadd == 0)
          break;
      }

    if (nadd)
      app_warning() << "WalkerControlBase::applyNmaxNmin not able to add sufficient walkers overall!" << std::endl;
  }

  // check current population
  current_population = std::accumulate(NumPerNode.begin(), NumPerNode.end(), 0);
  // at least one walker after load-balancing
  if (current_population / num_contexts_ == 0)
  {
    app_error() << "Some MPI ranks have no walkers after load balancing. This should not happen."
                << "Improve the trial wavefunction or adjust the simulation parameters." << std::endl;
    APP_ABORT("WalkerControlBase::sortWalkers");
  }

  return current_population;
}


/** copy good walkers to W
 *
 *  Good walkers are copied based on the registered number of copies
 *  Bad walkers are recycled to avoid memory allocation and deallocation.
 */
int WalkerControlBase::copyWalkers(MCWalkerConfiguration& W)
{
  // save current good walker size.
  const int size_good_w = good_w.size();
  std::vector<int> copy_list;
  for (int i = 0; i < size_good_w; i++)
  {
    for (int j = 0; j < ncopy_w[i]; j++)
    {
      if (bad_w.empty())
      {
        good_w.push_back(nullptr);
      }
      else
      {
        good_w.push_back(bad_w.back());
        bad_w.pop_back();
      }
      copy_list.push_back(i);
    }
  }

#pragma omp parallel for
  for (int i = size_good_w; i < good_w.size(); i++)
  {
    auto& wRef    = good_w[copy_list[i - size_good_w]];
    auto& awalker = good_w[i];
    if (awalker == nullptr)
      awalker = new Walker_t(*wRef);
    else
      *awalker = *wRef;
    // not fully sure this is correct or even used
    awalker->ID       = (i - size_good_w) * num_contexts_ + MyContext;
    awalker->ParentID = wRef->ParentID;
  }

  //clear the WalkerList to populate them with the good walkers
  W.clear();
  W.insert(W.begin(), good_w.begin(), good_w.end());

  //remove bad walkers if there is any left
  for (int i = 0; i < bad_w.size(); i++)
    delete bad_w[i];

  //clear good_w and ncopy_w for the next branch
  good_w.clear();
  bad_w.clear();
  ncopy_w.clear();
  return W.getActiveWalkers();
}

bool WalkerControlBase::put(xmlNodePtr cur)
{
  int nw_target = 0, nw_max = 0;
  std::string nonblocking = "yes";
  ParameterSet params;
  params.add(target_sigma_, "sigmaBound", "double");
  params.add(MaxCopy, "maxCopy", "int");
  params.add(nw_target, "targetwalkers", "int");
  params.add(nw_max, "max_walkers", "int");
  params.add(nonblocking, "use_nonblocking", "string");

  bool success = params.put(cur);

  // validating input
  if (nonblocking == "yes")
  {
    use_nonblocking = true;
  }
  else if (nonblocking == "no")
  {
    use_nonblocking = false;
  }
  else
  {
    APP_ABORT("WalkerControlBase::put unknown use_nonblocking option " + nonblocking);
  }

  setMinMax(nw_target, nw_max);

  app_log() << "  WalkerControlBase parameters " << std::endl;
  //app_log() << "    energyBound = " << targetEnergyBound << std::endl;
  //app_log() << "    sigmaBound = " << targetSigma << std::endl;
  app_log() << "    maxCopy = " << MaxCopy << std::endl;
  app_log() << "    Max Walkers per MPI rank " << n_max_ << std::endl;
  app_log() << "    Min Walkers per MPI rank " << n_min_ << std::endl;
  app_log() << "    Using " << (use_nonblocking ? "non-" : "") << "blocking send/recv" << std::endl;
  return true;
}

void WalkerControlBase::setMinMax(int nw_in, int nmax_in)
{
  if (nw_in > 0)
  {
    int npernode = nw_in / num_contexts_;
    if (method_)
    {
      n_max_ = npernode;
      n_min_ = npernode;
    }
    else
    {
      n_max_ = MaxCopy * npernode + 1;
      n_min_ = npernode / 5 + 1;
      if (nmax_in > 0)
        n_max_ = nmax_in;
    }
  }
}

} // namespace qmcplusplus
