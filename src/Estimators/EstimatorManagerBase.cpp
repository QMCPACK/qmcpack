//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Particle/MCWalkerConfiguration.h"
#include "Estimators/EstimatorManagerBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommUtilities.h"
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "Estimators/RMCLocalEnergyEstimator.h"
#include "Estimators/CollectablesEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
#include "Estimators/CSEnergyEstimator.h"
//leave it for serialization debug
//#define DEBUG_ESTIMATOR_ARCHIVE

namespace qmcplusplus
{

/** enumeration for EstimatorManagerBase.Options
 */
enum {COLLECT=0,
      MANAGE,
      RECORD,
      POSTIRECV,
      APPEND
     };

//initialize the name of the primary estimator
EstimatorManagerBase::EstimatorManagerBase(Communicate* c)
  : RecordCount(0),h_file(-1), FieldWidth(20)
  , MainEstimatorName("LocalEnergy"), Archive(0), DebugArchive(0)
  , myComm(0), MainEstimator(0), Collectables(0)
  , max4ascii(8)
{
  setCommunicator(c);
}

EstimatorManagerBase::EstimatorManagerBase(EstimatorManagerBase& em)
  : RecordCount(0),h_file(-1), FieldWidth(20)
  , MainEstimatorName(em.MainEstimatorName), Options(em.Options), Archive(0), DebugArchive(0)
  , myComm(0), MainEstimator(0), Collectables(0)
  , EstimatorMap(em.EstimatorMap), max4ascii(em.max4ascii)
{
  //inherit communicator
  setCommunicator(em.myComm);
  for(int i=0; i<em.Estimators.size(); i++)
    Estimators.push_back(em.Estimators[i]->clone());
  MainEstimator=Estimators[EstimatorMap[MainEstimatorName]];
  if(em.Collectables)
    Collectables=em.Collectables->clone();
}

EstimatorManagerBase::~EstimatorManagerBase()
{
  delete_iter(Estimators.begin(), Estimators.end());
  delete_iter(RemoteData.begin(), RemoteData.end());
  delete_iter(h5desc.begin(), h5desc.end());
  if(Collectables)
    delete Collectables;
}

void EstimatorManagerBase::setCommunicator(Communicate* c)
{
  if(myComm && myComm == c)
    return;
  myComm = c ? c:OHMMS::Controller;
  //set the default options
  Options.set(COLLECT,myComm->size()>1);
  Options.set(MANAGE,myComm->rank() == 0);
  if(RemoteData.empty())
  {
    RemoteData.push_back(new BufferType);
    RemoteData.push_back(new BufferType);
  }
}

/** set CollectSum
 * @param collect if true, global sum is done over the values
 */
void EstimatorManagerBase::setCollectionMode(bool collect)
{
  if(!myComm)
    setCommunicator(0);
  Options.set(COLLECT,(myComm->size() == 1)? false:collect);
  //force to be false for serial runs
  //CollectSum = (myComm->size() == 1)? false:collect;
}

/** reset names of the properties
 *
 * The number of estimators and their order can vary from the previous state.
 * Clear properties before setting up a new BlockAverage data list.
 */
void EstimatorManagerBase::reset()
{
  weightInd = BlockProperties.add("BlockWeight");
  cpuInd = BlockProperties.add("BlockCPU");
  acceptInd = BlockProperties.add("AcceptRatio");
  BlockAverages.clear();//cleaup the records
  for(int i=0; i<Estimators.size(); i++)
    Estimators[i]->add2Record(BlockAverages);
  max4ascii+=BlockAverages.size();
  if(Collectables)
    Collectables->add2Record(BlockAverages);
}

void EstimatorManagerBase::resetTargetParticleSet(ParticleSet& p)
{
}

void EstimatorManagerBase::addHeader(std::ostream& o)
{
  o.setf(std::ios::scientific, std::ios::floatfield);
  o.setf(std::ios::left,std::ios::adjustfield);
  o.precision(10);
  for (int i=0; i<BlockAverages.size(); i++)
    FieldWidth = std::max(FieldWidth, BlockAverages.Names[i].size()+2);
  for (int i=0; i<BlockProperties.size(); i++)
    FieldWidth = std::max(FieldWidth, BlockProperties.Names[i].size()+2);
  int maxobjs=std::min(BlockAverages.size(),max4ascii);
  o << "#   index    ";
  for(int i=0; i<maxobjs; i++)
    o << std::setw(FieldWidth) << BlockAverages.Names[i];
  for(int i=0; i<BlockProperties.size(); i++)
    o << std::setw(FieldWidth) << BlockProperties.Names[i];
  o << std::endl;
  o.setf(std::ios::right,std::ios::adjustfield);
}

void EstimatorManagerBase::start(int blocks, bool record)
{
  for(int i=0; i<Estimators.size(); i++)
    Estimators[i]->setNumberOfBlocks(blocks);
  reset();
  RecordCount=0;
  energyAccumulator.clear();
  varAccumulator.clear();
  int nc=(Collectables)?Collectables->size():0;
  BlockAverages.setValues(0.0);
  AverageCache.resize(BlockAverages.size()+nc);
  SquaredAverageCache.resize(BlockAverages.size()+nc);
  PropertyCache.resize(BlockProperties.size());
  //count the buffer size for message
  BufferSize=2*AverageCache.size()+PropertyCache.size();
  int sources=2;
  //allocate buffer for data collection
  if(RemoteData.empty())
    for(int i=0; i<sources; ++i)
      RemoteData.push_back(new BufferType(BufferSize));
  else
    for(int i=0; i<RemoteData.size(); ++i)
      RemoteData[i]->resize(BufferSize);
#if defined(DEBUG_ESTIMATOR_ARCHIVE)
  if(record && DebugArchive ==0)
  {
    char fname[128];
    sprintf(fname,"%s.p%03d.scalar.dat",myComm->getName().c_str(), myComm->rank());
    DebugArchive = new std::ofstream(fname);
    addHeader(*DebugArchive);
  }
#endif
  //set Options[RECORD] to enable/disable output
  Options.set(RECORD,record&&Options[MANAGE]);
  if(Options[RECORD])
  {
    if(Archive)
      delete Archive;
    std::string fname(myComm->getName());
    fname.append(".scalar.dat");
    Archive = new std::ofstream(fname.c_str());
    addHeader(*Archive);
    if(h5desc.size())
    {
      delete_iter(h5desc.begin(),h5desc.end());
      h5desc.clear();
    }
    fname=myComm->getName()+".stat.h5";
    h_file= H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    for(int i=0; i<Estimators.size(); i++)
      Estimators[i]->registerObservables(h5desc,h_file);
    if(Collectables)
      Collectables->registerObservables(h5desc,h_file);
  }
}

void EstimatorManagerBase::stop(const std::vector<EstimatorManagerBase*> est)
{
  int num_threads=est.size();
  //normalize by the number of threads per node
  RealType tnorm=1.0/static_cast<RealType>(num_threads);
  //add averages and divide them by the number of threads
  AverageCache=est[0]->AverageCache;
  for(int i=1; i<num_threads; i++)
    AverageCache+=est[i]->AverageCache;
  AverageCache*=tnorm;
  SquaredAverageCache=est[0]->SquaredAverageCache;
  for(int i=1; i<num_threads; i++)
    SquaredAverageCache+=est[i]->SquaredAverageCache;
  SquaredAverageCache*=tnorm;
  //add properties and divide them by the number of threads except for the weight
  PropertyCache=est[0]->PropertyCache;
  for(int i=1; i<num_threads; i++)
    PropertyCache+=est[i]->PropertyCache;
  for(int i=1; i<PropertyCache.size(); i++)
    PropertyCache[i] *= tnorm;
  stop();
}

/** Stop a run
 *
 * Collect data in Cache and print out data into hdf5 and ascii file.
 * This should not be called in a OpenMP parallel region or should
 * be guarded by master/single.
 * Keep the ascii output for now
 */
void EstimatorManagerBase::stop()
{
  //close any open files
  if(Archive)
  {
    delete Archive;
    Archive=0;
  } 
  if (h_file != -1)
  {
    H5Fclose(h_file);
    h_file = -1;
  }
}


void EstimatorManagerBase::startBlock(int steps)
{
  MyTimer.restart();
  BlockWeight=0.0;
}

/** take statistics of a block
 * @param accept acceptance rate of this block
 * @param collectall if true, need to gather data over MPI tasks
 */
void EstimatorManagerBase::stopBlock(RealType accept, bool collectall)
{
  //take block averages and update properties per block
  PropertyCache[weightInd]=BlockWeight;
  PropertyCache[cpuInd] = MyTimer.elapsed();
  PropertyCache[acceptInd] = accept;
  for(int i=0; i<Estimators.size(); i++)
    Estimators[i]->takeBlockAverage(AverageCache.begin(),SquaredAverageCache.begin());
  if(Collectables)
  {
    Collectables->takeBlockAverage(AverageCache.begin(),SquaredAverageCache.begin());
  }
  if(collectall)
    collectBlockAverages(1);
}

void EstimatorManagerBase::stopBlock(const std::vector<EstimatorManagerBase*>& est)
{
  //normalized it by the thread
  int num_threads=est.size();
  RealType tnorm=1.0/num_threads;
  AverageCache=est[0]->AverageCache;
  for(int i=1; i<num_threads; i++)
    AverageCache +=est[i]->AverageCache;
  AverageCache *= tnorm;
  SquaredAverageCache=est[0]->SquaredAverageCache;
  for(int i=1; i<num_threads; i++)
    SquaredAverageCache +=est[i]->SquaredAverageCache;
  SquaredAverageCache *= tnorm;
  PropertyCache=est[0]->PropertyCache;
  for(int i=1; i<num_threads; i++)
    PropertyCache+=est[i]->PropertyCache;
  for(int i=1; i<PropertyCache.size(); i++)
    PropertyCache[i] *= tnorm;
  //for(int i=0; i<num_threads; ++i)
  //varAccumulator(est[i]->varAccumulator.mean());
  collectBlockAverages(num_threads);
}


void EstimatorManagerBase::collectBlockAverages(int num_threads)
{
  if(Options[COLLECT])
  {
    //copy cached data to RemoteData[0]
    int n1=AverageCache.size();
    int n2=n1+AverageCache.size();
    int n3=n2+PropertyCache.size();
    {
      BufferType::iterator cur(RemoteData[0]->begin());
      copy(AverageCache.begin(),AverageCache.end(),cur);
      copy(SquaredAverageCache.begin(),SquaredAverageCache.end(),cur+n1);
      copy(PropertyCache.begin(),PropertyCache.end(),cur+n2);
    }
    myComm->reduce(*RemoteData[0]);
    if(Options[MANAGE])
    {
      BufferType::iterator cur(RemoteData[0]->begin());
      copy(cur,cur+n1, AverageCache.begin());
      copy(cur+n1,cur+n2, SquaredAverageCache.begin());
      copy(cur+n2,cur+n3, PropertyCache.begin());
      RealType nth=1.0/static_cast<RealType>(myComm->size());
      AverageCache *= nth;
      SquaredAverageCache *= nth;
      //do not weight weightInd
      for(int i=1; i<PropertyCache.size(); i++)
        PropertyCache[i] *= nth;
    }
  }
  //add the block average to summarize
  energyAccumulator(AverageCache[0]);
  varAccumulator(SquaredAverageCache[0]-AverageCache[0]*AverageCache[0]);
  if(Archive)
  {
    *Archive << std::setw(10) << RecordCount;
    int maxobjs=std::min(BlockAverages.size(),max4ascii);
    for(int j=0; j<maxobjs; j++)
      *Archive << std::setw(FieldWidth) << AverageCache[j];
    for(int j=0; j<PropertyCache.size(); j++)
      *Archive << std::setw(FieldWidth) << PropertyCache[j];
    *Archive << std::endl;
    for(int o=0; o<h5desc.size(); ++o)
      h5desc[o]->write(AverageCache.data(),SquaredAverageCache.data());
    H5Fflush(h_file,H5F_SCOPE_LOCAL);
  }
  RecordCount++;
}

/** accumulate Local energies and collectables
 * @param W ensemble
 */
void EstimatorManagerBase::accumulate(MCWalkerConfiguration& W)
{
  BlockWeight += W.getActiveWalkers();
  RealType norm=1.0/W.getGlobalNumWalkers();
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->accumulate(W,W.begin(),W.end(),norm);
  if(Collectables)//collectables are normalized by QMC drivers
    Collectables->accumulate_all(W.Collectables,1.0);
}

void EstimatorManagerBase::accumulate(MCWalkerConfiguration& W
                                  , MCWalkerConfiguration::iterator it, MCWalkerConfiguration::iterator it_end)
{
  BlockWeight += it_end-it;
  RealType norm=1.0/W.getGlobalNumWalkers();
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->accumulate(W,it,it_end,norm);
  if(Collectables)
    Collectables->accumulate_all(W.Collectables,1.0);
}

void EstimatorManagerBase::getEnergyAndWeight(RealType& e, RealType& w, RealType& var)
{
  if(Options[COLLECT])//need to broadcast the value
  {
    RealType tmp[3];
    tmp[0]= energyAccumulator.result();
    tmp[1]= energyAccumulator.count();
    tmp[2]= varAccumulator.mean();
    myComm->bcast(tmp,3);
    e=tmp[0];
    w=tmp[1];
    var=tmp[2];
  }
  else
  {
    e= energyAccumulator.result();
    w= energyAccumulator.count();
    var= varAccumulator.mean();
  }
}

void EstimatorManagerBase::getCurrentStatistics(MCWalkerConfiguration& W
    , RealType& eavg, RealType& var)
{
  LocalEnergyOnlyEstimator energynow;
  energynow.clear();
  energynow.accumulate(W,W.begin(),W.end(),1.0);
  std::vector<RealType> tmp(3);
  tmp[0]= energynow.scalars[0].result();
  tmp[1]= energynow.scalars[0].result2();
  tmp[2]= energynow.scalars[0].count();
  myComm->allreduce(tmp);
  eavg=tmp[0]/tmp[2];
  var=tmp[1]/tmp[2]-eavg*eavg;
}

EstimatorManagerBase::EstimatorType* EstimatorManagerBase::getMainEstimator()
{
  if(MainEstimator==0)
    add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
  return MainEstimator;
}

EstimatorManagerBase::EstimatorType* EstimatorManagerBase::getEstimator(const std::string& a)
{
  std::map<std::string,int>::iterator it = EstimatorMap.find(a);
  if(it == EstimatorMap.end())
    return 0;
  else
    return Estimators[(*it).second];
}

/** This should be moved to branch engine */
bool EstimatorManagerBase::put(MCWalkerConfiguration& W, QMCHamiltonian& H, xmlNodePtr cur)
{
  std::vector<std::string> extra;
  cur = cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "estimator")
    {
      std::string est_name(MainEstimatorName);
      std::string use_hdf5("yes");
      OhmmsAttributeSet hAttrib;
      hAttrib.add(est_name, "name");
      hAttrib.add(use_hdf5, "hdf5");
      hAttrib.put(cur);
      if( (est_name == MainEstimatorName) || (est_name=="elocal") )
      {
        max4ascii=H.sizeOfObservables()+3;
        add(new LocalEnergyEstimator(H,use_hdf5=="yes"),MainEstimatorName);
      }
      else
        if (est_name=="RMC")
        {
          int nobs(20);
          OhmmsAttributeSet hAttrib;
          hAttrib.add(nobs, "nobs");
          hAttrib.put(cur);
          max4ascii=nobs*H.sizeOfObservables()+3;
          add(new RMCLocalEnergyEstimator(H,nobs),MainEstimatorName);
        }
        else if ( est_name=="CSLocalEnergy" )
        {
          OhmmsAttributeSet hAttrib;
          int nPsi=1;
          hAttrib.add(nPsi, "nPsi");
          hAttrib.put(cur);	      
	      add(new CSEnergyEstimator(H,nPsi), MainEstimatorName);
	       app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
	    }
        else  
          extra.push_back(est_name);
    }
    cur = cur->next;
  }
  if(Estimators.empty())
  {
    app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
    max4ascii=H.sizeOfObservables()+3;
    add(new LocalEnergyEstimator(H,true),MainEstimatorName);
  }
  //Collectables is special and should not be added to Estimators
  if(Collectables == 0 && H.sizeOfCollectables())
  {
    app_log() << "  Using CollectablesEstimator for collectables, e.g. sk, gofr, density " << std::endl;
    Collectables=new CollectablesEstimator(H);
  }
  return true;
}

int EstimatorManagerBase::add(EstimatorType* newestimator, const std::string& aname)
{
  std::map<std::string,int>::iterator it = EstimatorMap.find(aname);
  int n =  Estimators.size();
  if(it == EstimatorMap.end())
  {
    Estimators.push_back(newestimator);
    EstimatorMap[aname] = n;
  }
  else
  {
    n=(*it).second;
    app_log() << "  EstimatorManagerBase::add replace " << aname << " estimator." << std::endl;
    delete Estimators[n];
    Estimators[n]=newestimator;
  }
  //check the name and set the MainEstimator
  if(aname == MainEstimatorName)
    MainEstimator=newestimator;
  return n;
}

int EstimatorManagerBase::addObservable(const char* aname)
{
  int mine = BlockAverages.add(aname);
  int add = TotalAverages.add(aname);
  if(mine < Block2Total.size())
    Block2Total[mine] = add;
  else
    Block2Total.push_back(add);
  return mine;
}

void EstimatorManagerBase::getData(int i, std::vector<RealType>& values)
{
  int entries = TotalAveragesData.rows();
  values.resize(entries);
  for (int a=0; a<entries; a++)
    values[a] = TotalAveragesData(a,Block2Total[i]);
}

//void EstimatorManagerBase::updateRefEnergy()
//{
//  CumEnergy[0]+=1.0;
//  RealType et=AverageCache(RecordCount-1,0);
//  CumEnergy[1]+=et;
//  CumEnergy[2]+=et*et;
//  CumEnergy[3]+=std::sqrt(MainEstimator->d_variance);

//  RealType wgtnorm=1.0/CumEnergy[0];
//  RefEnergy[0]=CumEnergy[1]*wgtnorm;
//  RefEnergy[1]=CumEnergy[2]*wgtnorm-RefEnergy[0]*RefEnergy[0];
//  if(CumEnergy[0]>1)
//    RefEnergy[2]=std::sqrt(RefEnergy[1]*wgtnorm/(CumEnergy[0]-1.0));
//  RefEnergy[3]=CumEnergy[3]*wgtnorm;//average block variance
//}

}

