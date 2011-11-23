////////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
#include "Particle/MCWalkerConfiguration.h"
#include "Estimators/EstimatorManager.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommUtilities.h"
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/ReleasedNodeEnergyEstimator.h"
#include "Estimators/AlternateReleasedNodeEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "Estimators/WFMCOnlyEstimator.h"
#include "Estimators/LocalEnergyEstimatorHDF.h"
#include "Estimators/FWEstimator.h"
#include "Estimators/CollectablesEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
//#define QMC_ASYNC_COLLECT
//leave it for serialization debug
//#define DEBUG_ESTIMATOR_ARCHIVE

namespace qmcplusplus {

  /** enumeration for EstimatorManager.Options
   */
  enum {COLLECT=0, 
  MANAGE,
  RECORD, 
  POSTIRECV, 
  APPEND};

  //initialize the name of the primary estimator
  EstimatorManager::EstimatorManager(Communicate* c)
    : RecordCount(0),h_file(-1), FieldWidth(20)
      , MainEstimatorName("LocalEnergy"), Archive(0), DebugArchive(0)
      , myComm(0), MainEstimator(0), Collectables(0)
      , max4ascii(8), pendingRequests(0)
  { 
    setCommunicator(c);
  }

  EstimatorManager::EstimatorManager(EstimatorManager& em)
    : RecordCount(0),h_file(-1), FieldWidth(20)
      , MainEstimatorName(em.MainEstimatorName), Options(em.Options), Archive(0), DebugArchive(0)
      , myComm(0), MainEstimator(0), Collectables(0) 
      , EstimatorMap(em.EstimatorMap), max4ascii(em.max4ascii), pendingRequests(0)
  {
    //inherit communicator
    setCommunicator(em.myComm);
    for(int i=0; i<em.Estimators.size(); i++) 
      Estimators.push_back(em.Estimators[i]->clone());
    MainEstimator=Estimators[EstimatorMap[MainEstimatorName]];
    if(em.Collectables) Collectables=em.Collectables->clone();
  }

  EstimatorManager::~EstimatorManager()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
    delete_iter(RemoteData.begin(), RemoteData.end());
    delete_iter(h5desc.begin(), h5desc.end());
    if(Collectables) delete Collectables;
  }

  void EstimatorManager::setCommunicator(Communicate* c) 
  {
    if(myComm && myComm == c) return;
    myComm = c ? c:OHMMS::Controller;

    //set the default options
    Options.set(COLLECT,myComm->size()>1);
    Options.set(MANAGE,myComm->rank() == 0);

    myRequest.resize(2);
    if(Options[COLLECT] && Options[MANAGE]) myRequest.resize(myComm->size()-1);

    if(RemoteData.empty())
    {
#if defined(QMC_ASYNC_COLLECT)
      for(int i=0; i<myComm->size(); i++) RemoteData.push_back(new BufferType);
#else
      RemoteData.push_back(new BufferType);
      RemoteData.push_back(new BufferType);
#endif
    }
  }

  /** set CollectSum
   * @param collect if true, global sum is done over the values
   */
  void EstimatorManager::setCollectionMode(bool collect) 
  {
    if(!myComm) setCommunicator(0);
    Options.set(COLLECT,(myComm->size() == 1)? false:collect);
    //force to be false for serial runs
    //CollectSum = (myComm->size() == 1)? false:collect;
  }

  /** reset names of the properties 
   *
   * The number of estimators and their order can vary from the previous state.
   * Clear properties before setting up a new BlockAverage data list.
   */
  void EstimatorManager::reset()
  {

    weightInd = BlockProperties.add("BlockWeight");
    cpuInd = BlockProperties.add("BlockCPU");
    acceptInd = BlockProperties.add("AcceptRatio");

    BlockAverages.clear();//cleaup the records
    for(int i=0; i<Estimators.size(); i++) 
      Estimators[i]->add2Record(BlockAverages);

    if(Collectables) Collectables->add2Record(BlockAverages);
  }

  void EstimatorManager::resetTargetParticleSet(ParticleSet& p)
  {
  }

  void EstimatorManager::addHeader(ostream& o)
  {
    o.setf(ios::scientific, ios::floatfield);
    o.setf(ios::left,ios::adjustfield);
    o.precision(10);
    for (int i=0; i<BlockAverages.size(); i++)
      FieldWidth = std::max(FieldWidth, BlockAverages.Names[i].size()+2);
    for (int i=0; i<BlockProperties.size(); i++)
      FieldWidth = std::max(FieldWidth, BlockProperties.Names[i].size()+2);

    int maxobjs=std::min(BlockAverages.size(),max4ascii);
    o << "#   index    ";
    for(int i=0; i<maxobjs; i++) o << setw(FieldWidth) << BlockAverages.Names[i];
    for(int i=0; i<BlockProperties.size(); i++) o << setw(FieldWidth) << BlockProperties.Names[i];
    o << endl;
    o.setf(ios::right,ios::adjustfield);
  }

  void EstimatorManager::start(int blocks, bool record)
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

#if defined(QMC_ASYNC_COLLECT)
    int sources=myComm->size();
#else
    int sources=2;
#endif

    //allocate buffer for data collection
    if(RemoteData.empty())
      for(int i=0; i<sources; ++i) RemoteData.push_back(new BufferType(BufferSize));
    else
      for(int i=0; i<RemoteData.size(); ++i) RemoteData[i]->resize(BufferSize);

#if defined(DEBUG_ESTIMATOR_ARCHIVE)
    if(record && DebugArchive ==0)
    {
      char fname[128];
      sprintf(fname,"%s.p%03d.scalar.dat",myComm->getName().c_str(), myComm->rank());
      DebugArchive = new ofstream(fname);
      addHeader(*DebugArchive);
    }
#endif

    //set Options[RECORD] to enable/disable output 
    Options.set(RECORD,record&&Options[MANAGE]);

    if(Options[RECORD])
    {
      if(Archive) delete Archive;
      string fname(myComm->getName());
      fname.append(".scalar.dat");
      Archive = new ofstream(fname.c_str());
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

      if(Collectables) Collectables->registerObservables(h5desc,h_file);

#if defined(QMC_ASYNC_COLLECT)
      if(Options[COLLECT])
      {//issue a irecv
        pendingRequests=0;
        for(int i=1,is=0; i<myComm->size(); i++,is++) 
        {
          myRequest[is]=myComm->irecv(i,i,*RemoteData[i]);//request only has size-1
          pendingRequests++;
        }
      }
#endif
    }
  }

  void EstimatorManager::stop(const vector<EstimatorManager*> est)
  {

    int num_threads=est.size();
    //normalize by the number of threads per node
    RealType tnorm=1.0/static_cast<RealType>(num_threads);

    //add averages and divide them by the number of threads
    AverageCache=est[0]->AverageCache;
    for(int i=1; i<num_threads; i++) AverageCache+=est[i]->AverageCache;
    AverageCache*=tnorm;

    SquaredAverageCache=est[0]->SquaredAverageCache;
    for(int i=1; i<num_threads; i++) SquaredAverageCache+=est[i]->SquaredAverageCache;
    SquaredAverageCache*=tnorm;

    //add properties and divide them by the number of threads except for the weight
    PropertyCache=est[0]->PropertyCache;
    for(int i=1; i<num_threads; i++) PropertyCache+=est[i]->PropertyCache;
    for(int i=1; i<PropertyCache.size(); i++) PropertyCache[i] *= tnorm;

    stop();
  }

  /** Stop a run
   *
   * Collect data in Cache and print out data into hdf5 and ascii file.
   * This should not be called in a OpenMP parallel region or should
   * be guarded by master/single.
   * Keep the ascii output for now
   */
  void EstimatorManager::stop() 
  {
    //clean up pending messages
    if(pendingRequests) 
    {
      cancel(myRequest);
      pendingRequests=0;
    }

    //close any open files
    if(Archive) 
    {
      delete Archive; Archive=0;
    }
  }

  void EstimatorManager::startBlock(int steps)
  { 
    MyTimer.restart();
    BlockWeight=0.0;
  }

  /** take statistics of a block
   * @param accept acceptance rate of this block
   * @param collectall if true, need to gather data over MPI tasks
   */
  void EstimatorManager::stopBlock(RealType accept, bool collectall)
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

    if(collectall) collectBlockAverages(1);
  }

  void EstimatorManager::stopBlock(const vector<EstimatorManager*>& est)
  {

    //normalized it by the thread
    int num_threads=est.size();
    RealType tnorm=1.0/num_threads;

    AverageCache=est[0]->AverageCache;
    for(int i=1; i<num_threads; i++) AverageCache +=est[i]->AverageCache;
    AverageCache *= tnorm;

    SquaredAverageCache=est[0]->SquaredAverageCache;
    for(int i=1; i<num_threads; i++) SquaredAverageCache +=est[i]->SquaredAverageCache;
    SquaredAverageCache *= tnorm;

    PropertyCache=est[0]->PropertyCache;
    for(int i=1; i<num_threads; i++) PropertyCache+=est[i]->PropertyCache;
    for(int i=1; i<PropertyCache.size(); i++) PropertyCache[i] *= tnorm;

    //for(int i=0; i<num_threads; ++i)
      //varAccumulator(est[i]->varAccumulator.mean()); 

    collectBlockAverages(num_threads);
  }


  void EstimatorManager::collectBlockAverages(int num_threads)
  {
    if(Options[COLLECT])
    { //copy cached data to RemoteData[0]
      int n1=AverageCache.size();
      int n2=n1+AverageCache.size();
      int n3=n2+PropertyCache.size();
      {
        BufferType::iterator cur(RemoteData[0]->begin());
        std::copy(AverageCache.begin(),AverageCache.end(),cur);
        std::copy(SquaredAverageCache.begin(),SquaredAverageCache.end(),cur+n1);
        std::copy(PropertyCache.begin(),PropertyCache.end(),cur+n2);
      }

#if defined(QMC_ASYNC_COLLECT)
      if(Options[MANAGE]) 
      { //wait all the message but we can choose to wait one-by-one with a timer
        wait_all(myRequest.size(),&myRequest[0]);
        for(int is=1; is<myComm->size(); is++) 
          accumulate_elements(RemoteData[is]->begin(),RemoteData[is]->end(), RemoteData[0]->begin());
      } 
      else //not a master, pack and send the data
        myRequest[0]=myComm->isend(0,myComm->rank(),*RemoteData[0]);
#else
      myComm->reduce(*RemoteData[0]);
#endif
      if(Options[MANAGE])
      {
        BufferType::iterator cur(RemoteData[0]->begin());
        std::copy(cur,cur+n1, AverageCache.begin());
        std::copy(cur+n1,cur+n2, SquaredAverageCache.begin());
        std::copy(cur+n2,cur+n3, PropertyCache.begin());

        RealType nth=1.0/static_cast<RealType>(myComm->size());
        AverageCache *= nth;
        SquaredAverageCache *= nth;
        //do not weight weightInd
        for(int i=1; i<PropertyCache.size(); i++) PropertyCache[i] *= nth;
      }
    }

    //add the block average to summarize
    energyAccumulator(AverageCache[0]);
    varAccumulator(SquaredAverageCache[0]-AverageCache[0]*AverageCache[0]);

    if(Archive)
    {
      *Archive << setw(10) << RecordCount;
      int maxobjs=std::min(BlockAverages.size(),max4ascii);
      for(int j=0; j<maxobjs; j++) *Archive << setw(FieldWidth) << AverageCache[j];
      for(int j=0; j<PropertyCache.size(); j++) *Archive << setw(FieldWidth) << PropertyCache[j];
      *Archive << endl;

      for(int o=0; o<h5desc.size(); ++o)
        h5desc[o]->write(AverageCache.data(),SquaredAverageCache.data());
      H5Fflush(h_file,H5F_SCOPE_LOCAL);
    }
    RecordCount++;
  }

  /** accumulate Local energies and collectables
   * @param W ensemble
   */
  void EstimatorManager::accumulate(MCWalkerConfiguration& W)
  {
    BlockWeight += W.getActiveWalkers();
    RealType norm=1.0/W.getGlobalNumWalkers();

    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate(W,W.begin(),W.end(),norm);

    if(Collectables)//collectables are normalized by QMC drivers
      Collectables->accumulate_all(W.Collectables,1.0);
  }

  void EstimatorManager::accumulate(MCWalkerConfiguration& W 
     , MCWalkerConfiguration::iterator it, MCWalkerConfiguration::iterator it_end)
  {
    BlockWeight += it_end-it;
    RealType norm=1.0/W.getGlobalNumWalkers();

    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate(W,it,it_end,norm);

    if(Collectables)
      Collectables->accumulate_all(W.Collectables,1.0);
  }

  void EstimatorManager::accumulate( HDF5_FW_observables& OBS, HDF5_FW_weights& WGTS, vector<int>& Dims )
  {
    BlockWeight=1;
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate_fw(OBS,WGTS,Dims);
  }


  void EstimatorManager::getEnergyAndWeight(RealType& e, RealType& w, RealType& var) 
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

  void EstimatorManager::getCurrentStatistics(MCWalkerConfiguration& W
      , RealType& eavg, RealType& var)
  {
    LocalEnergyOnlyEstimator energynow;
    energynow.clear();
    energynow.accumulate(W,W.begin(),W.end(),1.0);
    vector<RealType> tmp(3);
    tmp[0]= energynow.scalars[0].result();
    tmp[1]= energynow.scalars[0].result2();
    tmp[2]= energynow.scalars[0].count();
    myComm->allreduce(tmp);
    eavg=tmp[0]/tmp[2];
    var=tmp[1]/tmp[2]-eavg*eavg;
  }

  EstimatorManager::EstimatorType* EstimatorManager::getMainEstimator() 
  {
    if(MainEstimator==0) 
      add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
    return MainEstimator;
  }

  EstimatorManager::EstimatorType* EstimatorManager::getEstimator(const string& a) 
  {
    std::map<string,int>::iterator it = EstimatorMap.find(a);
    if(it == EstimatorMap.end()) 
      return 0;
    else 
      return Estimators[(*it).second];
  }

  /** This should be moved to branch engine */
  bool EstimatorManager::put(MCWalkerConfiguration& W, QMCHamiltonian& H, xmlNodePtr cur) 
  {
    vector<string> extra;
    cur = cur->children;
    while(cur != NULL) 
    {
      string cname((const char*)(cur->name));
      if(cname == "estimator") 
      {
        string est_name(MainEstimatorName);
        string use_hdf5("yes");
        OhmmsAttributeSet hAttrib;
        hAttrib.add(est_name, "name");
        hAttrib.add(use_hdf5, "hdf5");
        hAttrib.put(cur);
        if( (est_name == MainEstimatorName) || (est_name=="elocal") )
        {
          if(use_hdf5 == "yes")
          {
            max4ascii=H.size()+3;//write only physical energies
            add(new LocalEnergyEstimatorHDF(H),MainEstimatorName);
          }
          else
          {//fall back to the ascii file
            max4ascii=H.sizeOfObservables()+3;
            add(new LocalEnergyEstimator(H),MainEstimatorName);
          }
        }
        else if (est_name=="WFMConly")
        {
          max4ascii=H.sizeOfObservables()+10;
          app_log() << "  Using WFMConly for the MainEstimator " << endl;
          add(new WFMCOnlyEstimator(H),MainEstimatorName);
          est_name=MainEstimatorName;
        }
        else if (est_name=="releasednode")
        { 
          int Smax(100);
          int primary(1);
          OhmmsAttributeSet hAttrib;
          hAttrib.add(Smax, "Smax");
          hAttrib.add(primary, "primary");
          hAttrib.put(cur);
        
          max4ascii=H.sizeOfObservables()+ 4 + 3*(Smax+1);
          app_log() << "  Using ReleasedNode for the MainEstimator with Smax="<<Smax<<" and max4ascii="<<max4ascii << endl;
          if (primary==2) add(new ReleasedNodeEnergyEstimator(H,Smax),MainEstimatorName);
          else add(new AlternateReleasedNodeEnergyEstimator(H,Smax),MainEstimatorName);
          est_name=MainEstimatorName;
        }
        else if (est_name=="forwardwalking")
        {
          max4ascii=2*H.sizeOfObservables()+4;
          app_log() << "  Doing forwardwalking on hdf5 " << endl;
          add(new ForwardWalkingEstimator(H),MainEstimatorName);
          est_name=MainEstimatorName;
        }        else 
          extra.push_back(est_name);
      } 
      cur = cur->next;
    }

    if(Estimators.empty()) 
    {
      app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << endl;
      max4ascii=H.sizeOfObservables()+3;
      add(new LocalEnergyEstimator(H),MainEstimatorName);
      //add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
    } 

    //Collectables is special and should not be added to Estimators
    if(Collectables == 0 && H.sizeOfCollectables())
    {
      app_log() << "  Using CollectablesEstimator for collectables, e.g. sk, gofr, density " << endl;
      Collectables=new CollectablesEstimator(H);
    }

    return true;
  }

  int EstimatorManager::add(EstimatorType* newestimator, const string& aname) 
  { 

    std::map<string,int>::iterator it = EstimatorMap.find(aname);
    int n =  Estimators.size();
    if(it == EstimatorMap.end()) 
    {
      Estimators.push_back(newestimator);
      EstimatorMap[aname] = n;
    } 
    else 
    {
      n=(*it).second;
      app_log() << "  EstimatorManager::add replace " << aname << " estimator." << endl;
      delete Estimators[n]; 
      Estimators[n]=newestimator;
    }

    //check the name and set the MainEstimator
    if(aname == MainEstimatorName) MainEstimator=newestimator;
    return n;
  }

  int EstimatorManager::addObservable(const char* aname) 
  {
    int mine = BlockAverages.add(aname);
    int add = TotalAverages.add(aname);
    if(mine < Block2Total.size()) 
      Block2Total[mine] = add;
    else 
      Block2Total.push_back(add);
    return mine;
  }

  void EstimatorManager::getData(int i, vector<RealType>& values)
  {
    int entries = TotalAveragesData.rows();
    values.resize(entries);
    for (int a=0; a<entries; a++)
      values[a] = TotalAveragesData(a,Block2Total[i]);
  }

  //void EstimatorManager::updateRefEnergy()
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

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
