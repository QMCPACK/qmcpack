//////////////////////////////////////////////////////////////////
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
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "Estimators/CompositeEstimators.h"
#include "Estimators/GofREstimator.h"
#include "Estimators/SkEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#include "Estimators/PolarizationEstimator.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"

namespace qmcplusplus {

  //initialize the name of the primary estimator
  EstimatorManager::EstimatorManager(Communicate* c): 
    MainEstimatorName("elocal"),
  Manager(false), CollectSum(true), AppendRecord(false), Collected(false), 
  ThreadCount(1), h_file(-1), h_obs(-1), myComm(0), //H(h),
  MainEstimator(0), CompEstimators(0)
  { 
    setCommunicator(c);
    //Block2Total.resize(0); 
  }

  EstimatorManager::EstimatorManager(EstimatorManager& em):
    MainEstimatorName("elocal"),
  Manager(false), AppendRecord(false), Collected(false), 
  ThreadCount(1),h_file(-1), h_obs(-1), 
  MainEstimator(0), CompEstimators(0),
  EstimatorMap(em.EstimatorMap)
  {
    CollectSum=em.CollectSum, 
    //inherit communicator
    setCommunicator(em.myComm);
    for(int i=0; i<em.Estimators.size(); i++) 
    {
      Estimators.push_back(em.Estimators[i]->clone());
    }
    MainEstimator=Estimators[EstimatorMap[MainEstimatorName]];
  }

  EstimatorManager::~EstimatorManager()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
    if(CompEstimators) delete CompEstimators;
  }

  void EstimatorManager::setCommunicator(Communicate* c) 
  {
    if(myComm && myComm == c) return;
    myComm = c ? c:OHMMS::Controller;
    Manager = myComm->master();
    CollectSum = (myComm->ncontexts()>1);
  }

  /** set CollectSum
   * @param collect if true, global sum is done over the values
   */
  void EstimatorManager::setCollectionMode(bool collect) 
  {
    if(!myComm) setCommunicator(0);
    //force to be false for serial runs
    CollectSum = (myComm->ncontexts() == 1)? false:collect;
    //for(int i=0; i< Estimators.size(); i++) Estimators[i]->CollectSum = CollectSum;
  }

  /** reset names of the properties 
   *
   * The number of estimators and their order can vary from the previous state.
   * Clear properties before setting up a new BlockAverage data list.
   */
  void EstimatorManager::reset()
  {

    if(Estimators.empty()) 
      add(new LocalEnergyOnlyEstimator(),MainEstimatorName);

    BlockAverages.clear();
    //no need for RemoteData. We don't collect anything during a run
    BufferType tmp;
    for(int i=0; i<Estimators.size(); i++) 
    {
      Estimators[i]->add2Record(BlockAverages,tmp);
      Estimators[i]->reset();
    }

    //weightInd = BlockProperties.add("WeightSum");
    cpuInd = BlockProperties.add("BlockCPU");
    acceptInd = BlockProperties.add("AcceptRatio");
  }

  void EstimatorManager::start(int blocks, bool record)
  {
    reset();

    Collected= false;
    RecordCount=0;
    BlockAverages.setValues(0.0);

    //resize does not reset the value
    TotalWeight.resize(blocks);

    AverageCache.resize(blocks,BlockAverages.size());
    AverageCache=0.0;

    PropertyCache.resize(blocks,BlockProperties.size());

    if(record)
    {
      if(h_obs>-1)  H5Gclose(h_obs);
      if(h_file>-1)  H5Fclose(h_file);

      //open obsevables and allocate spaces 
      string h5file=RootName+".config.h5";
      h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);

      herr_t status = H5Eset_auto(NULL, NULL);
      status = H5Gget_objinfo (h_file, "observables", 0, NULL);
      if(status ==0) //observables are already there. remove them
      {
        h_obs = H5Gunlink(h_file,"observables");
      }

      h_obs = H5Gcreate(h_file,"observables",0);
      HDFAttribIO<int> i(RecordCount);
      i.write(h_obs,"count");
      HDFAttribIO<int> t(ThreadCount);
      t.write(h_obs,"threads");
      HDFAttribIO<Matrix<RealType> > m(AverageCache);
      m.write(h_obs,"scalars");

      ostringstream o;
      for(int i=0; i<BlockAverages.size()-1; i++) o << BlockAverages.Name[i] << ":";
      o<<BlockAverages.Name.back();
      string banner(o.str());
      HDFAttribIO<string> so(banner);
      so.write(h_obs,"scalar_ids");

      if(CompEstimators) CompEstimators->open(h_obs);

      //H5Gclose(h_obs);
      //H5Fclose(h_file);
    }
  }

  void EstimatorManager::stop(const vector<EstimatorManager*> est)
  {
    RecordCount=est[0]->RecordCount;

    TotalWeight=est[0]->TotalWeight;
    for(int i=1; i<ThreadCount; i++)
      TotalWeight+=est[i]->TotalWeight;

    AverageCache=est[0]->AverageCache;
    for(int i=1; i<ThreadCount; i++)
      AverageCache+=est[i]->AverageCache;

    //normalize by the number of threads per node
    RealType wgt=1.0/static_cast<RealType>(ThreadCount);
    AverageCache*=wgt;

    PropertyCache=est[0]->PropertyCache;
    for(int i=1; i<ThreadCount; i++)
      PropertyCache+=est[i]->PropertyCache;

    PropertyCache*=wgt;

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
    if(CollectSum)
    {
      //save the local data to scalars_node so that we can trace it back to the node
      if(h_obs>-1)
      {
        HDFAttribIO<Matrix<RealType> > m(AverageCache);
        m.write(h_obs,"scalars_node");

        HDFAttribIO<Matrix<RealType> > p(PropertyCache);
        p.write(h_obs,"run_properties_node");
      }

      //we will use serialization
      int n1=AverageCache.size(), n2=PropertyCache.size(), n3=TotalWeight.size();
      vector<RealType> c(n1+n2+n3);
      std::copy(AverageCache.begin(),AverageCache.end(),c.begin());
      std::copy(PropertyCache.begin(),PropertyCache.end(),c.begin()+n1);
      std::copy(TotalWeight.begin(),TotalWeight.end(),c.begin()+n1+n2);

      //use allreduce
      myComm->allreduce(c);

      std::copy(c.begin(),c.begin()+n1,AverageCache.begin());
      std::copy(c.begin()+n1,c.begin()+n1+n2,PropertyCache.begin());
      std::copy(c.begin()+n1+n2,c.end(),TotalWeight.begin());

      RealType norm=1.0/static_cast<RealType>(myComm->ncontexts());
      AverageCache *=norm;
      PropertyCache *=norm;
    }

    //close the group
    if(h_obs>-1)
    {
      //reset to 1 to indicate that scalars do not have the count
      HDFAttribIO<int> i(RecordCount,true);
      i.write(h_obs,"count");
      HDFAttribIO<int> t(ThreadCount,true);
      t.write(h_obs,"threads");
      HDFAttribIO<Matrix<RealType> > m(AverageCache,true);
      m.write(h_obs,"scalars");

      //save the run-time properties
      ostringstream o;
      for(int i=0; i<BlockProperties.size()-1; i++) o << BlockProperties.Name[i] << ":";
      o<<BlockProperties.Name.back();
      string banner(o.str());
      HDFAttribIO<string> so(banner);
      so.write(h_obs,"run_properties_ids");

      HDFAttribIO<Matrix<RealType> > p(PropertyCache);
      p.write(h_obs,"run_properties");

      HDFAttribIO<Vector<RealType> > w(TotalWeight);
      p.write(h_obs,"weights");

      if(CompEstimators) CompEstimators->close();

      H5Gclose(h_obs); 
      h_obs=-1;
    }

    if(h_file>-1)
    {
      H5Fclose(h_file); 
      h_file=-1;
    }

    if(Manager) //write to a ascii file. Plan to disable it entirely.
    {
      string fname(RootName);
      fname.append(".scalar.dat");
      ofstream fout(fname.c_str());
      fout.setf(ios::scientific, ios::floatfield);
      fout.setf(ios::left,ios::adjustfield);
      fout << "#   index    ";
      for(int i=0; i<BlockAverages.size(); i++) fout << setw(16) << BlockAverages.Name[i];
      fout << setw(16) << "WeightSum";
      for(int i=0; i<BlockProperties.size(); i++) fout << setw(16) << BlockProperties.Name[i];
      fout << endl;
      fout.setf(ios::right,ios::adjustfield);
      const RealType* restrict rptr=AverageCache.data();
      const RealType* restrict pptr=PropertyCache.data();
      int pc=PropertyCache.cols();
      int nc=AverageCache.cols();
      for(int i=0; i<RecordCount; i++)
      {
        fout << setw(10) << i;
        for(int j=0; j<nc; j++) fout << setw(16) << *rptr++;
        fout << setw(16) << TotalWeight[i];
        for(int j=0; j<pc; j++) fout << setw(16) << *pptr++;
        fout << endl;
      }
    }

    Collected=true;
  }

  void EstimatorManager::startBlock(int steps)
  { 
    MyTimer.restart();
    //if(CompEstimators) CompEstimators->startBlock(steps);
  }

  void EstimatorManager::stopBlock(RealType accept)
  {
    TotalWeight[RecordCount]=Estimators[0]->d_wgt;
    //PropertyCache(RecordCount,weightInd) = Estimators[0]->d_wgt;
    PropertyCache(RecordCount,cpuInd)    = MyTimer.elapsed();
    PropertyCache(RecordCount,acceptInd) = accept;

    for(int i=0; i<Estimators.size(); i++) 
    {
      Estimators[i]->takeBlockAverage(AverageCache[RecordCount]);
    }
    RecordCount++;

    //group is closed. Do not save it to hdf
    if(h_obs<-1) return;

    HDFAttribIO<int> i(RecordCount,true);
    i.write(h_obs,"count");
    HDFAttribIO<Matrix<RealType> > m(AverageCache,true);
    m.write(h_obs,"scalars");

    if(CompEstimators) CompEstimators->stopBlock();
  }

  void EstimatorManager::stopBlock(const vector<EstimatorManager*> est)
  {
    ThreadCount=est.size();
    //RecordCount=est[0]->RecordCount;
    //AverageCache=est[0]->AverageCache;
    //for(int i=1; i<est.size(); i++)
    //{
    //  AverageCache+=est[i]->AverageCache;
    //}
    RecordCount=est[0]->RecordCount;
    for(int i=0; i<ThreadCount; i++)
    {
      int rc=est[i]->RecordCount-1;
      accumulate_elements(est[i]->PropertyCache[rc],est[i]->PropertyCache[rc+1],PropertyCache[rc]);
      accumulate_elements(est[i]->AverageCache[rc], est[i]->AverageCache[rc+1],AverageCache[rc]);
      //FastAccumulate(est[i]->PropertyCache[rc],est[i]->PropertyCache[rc+1],PropertyCache[rc]);
      //FastAccumulate(est[i]->AverageCache[rc], est[i]->AverageCache[rc+1],AverageCache[rc]);
    }

    if(h_obs<-1) return;

    HDFAttribIO<int> i(RecordCount,true);
    i.write(h_obs,"count");
    HDFAttribIO<int> t(ThreadCount,true);
    t.write(h_obs,"threads");
    HDFAttribIO<Matrix<RealType> > m(AverageCache,true);
    m.write(h_obs,"scalars");
  }

  void EstimatorManager::accumulate(MCWalkerConfiguration& W)
  {
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->accumulate(W.begin(),W.end());
    if(CompEstimators) CompEstimators->accumulate(W);
  }

  void EstimatorManager::accumulate(MCWalkerConfiguration::iterator it,
      MCWalkerConfiguration::iterator it_end)
  {
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->accumulate(it,it_end);
  }

  void EstimatorManager::accumulate(ParticleSet& P, 
      MCWalkerConfiguration::Walker_t& awalker) 
  {
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate(P,awalker);
  }

  void 
    EstimatorManager::getEnergyAndWeight(RealType& e, RealType& w) 
  {
    int nc=AverageCache.cols();
    if(nc==0) return;
    EPSum[0]=0.0;
    EPSum[1]=RecordCount;
    const RealType* restrict eptr=AverageCache.data();
    for(int i=0; i<RecordCount; i++, eptr+=nc) EPSum[0] += *eptr;
//#if defined(HAVE_MPI)
//    MPI_Bcast(EPSum.begin(),2,MPI_DOUBLE,0,myComm->getMPI());
//#endif
    e=EPSum[0];
    w=EPSum[1];
  }

  EstimatorManager::EstimatorType* 
    EstimatorManager::getMainEstimator() 
    {
      if(MainEstimator==0) 
        add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
      return MainEstimator;
    }

  EstimatorManager::EstimatorType* 
    EstimatorManager::getEstimator(const string& a) 
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
        xmlChar* att=xmlGetProp(cur,(const xmlChar*)"name");
        if(att) 
        {
          string aname((const char*)att);
          if(aname == "LocalEnergy") 
            add(new LocalEnergyEstimator(H),MainEstimatorName);
          else 
            extra.push_back(aname);
        }
      } 
      cur = cur->next;
    }

    if(Estimators.empty()) 
    {
      app_log() << "  Using a default LocalEnergyOnlyEstimator for the MainEstimator " << endl;
      add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
    } 

    string SkName("sk");
    string GofRName("gofr");
    for(int i=0; i< extra.size(); i++)
    {
      if(extra[i] == SkName && W.Lattice.SuperCellEnum)
      {
        if(CompEstimators == 0) CompEstimators = new CompositeEstimatorSet;
        if(CompEstimators->missing(SkName))
        {
          app_log() << "  EstimatorManager::add " << SkName << endl;
          CompEstimators->add(new SkEstimator(W),SkName);
        }
      } else if(extra[i] == GofRName)
      {
        if(CompEstimators == 0) CompEstimators = new CompositeEstimatorSet;
        if(CompEstimators->missing(GofRName))
        {
          app_log() << "  EstimatorManager::add " << GofRName << endl;
          CompEstimators->add(new GofREstimator(W),GofRName);
        }
      }
    }
    //add extra
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
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
