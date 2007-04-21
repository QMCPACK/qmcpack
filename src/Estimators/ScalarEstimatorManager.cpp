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
#include "Estimators/ScalarEstimatorManager.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#include "Estimators/PolarizationEstimator.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"

namespace qmcplusplus {

  template<typename IT>
    inline void FastAccumulate(IT first, IT last, IT result)
    {
      while(first != last) *result++ += *first++;
    }

  //initialize the name of the primary estimator
  ScalarEstimatorManager::ScalarEstimatorManager(QMCHamiltonian& h, 
      Communicate* c): 
    MainEstimatorName("elocal"),
  Manager(false), CollectSum(true), AppendRecord(false), Collected(false), 
  ThreadCount(1), h_file(-1), h_obs(-1), myComm(0), H(h),MainEstimator(0)
  { 
    //Block2Total.resize(0); 
  }

  ScalarEstimatorManager::ScalarEstimatorManager(ScalarEstimatorManager& em, 
      QMCHamiltonian& h):
    MainEstimatorName("elocal"),
  Manager(false), AppendRecord(false), Collected(false), 
  ThreadCount(1),h_file(-1), h_obs(-1), H(h), 
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

  ScalarEstimatorManager::~ScalarEstimatorManager()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
  }

  void ScalarEstimatorManager::setCommunicator(Communicate* c) 
  {
    if(myComm && myComm == c) return;
    myComm = c ? c:OHMMS::Controller;
    Manager = myComm->master();
    CollectSum = (myComm->ncontexts()>1);
  }

  /** set CollectSum
   * @param collect if true, global sum is done over the values
   */
  void ScalarEstimatorManager::setCollectionMode(bool collect) 
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
  void ScalarEstimatorManager::reset()
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

    weightInd = BlockProperties.add("WeightSum");
    cpuInd = BlockProperties.add("BlockCPU");
    acceptInd = BlockProperties.add("AcceptRatio");
  }

  void ScalarEstimatorManager::start(int blocks, bool record)
  {
    reset();

    Collected= false;
    RecordCount=0;
    BlockAverages.setValues(0.0);
    AverageCache.resize(blocks,BlockAverages.size());
    PropertyCache.resize(blocks,BlockProperties.size());

    AverageCache=0.0;

    if(record)
    {
      //open obsevables and allocate spaces 
      string h5file=RootName+".config.h5";
      h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
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

      //H5Gclose(h_obs);
      H5Fclose(h_file);
    }
  }

  void ScalarEstimatorManager::stop(const vector<ScalarEstimatorManager*> est)
  {
    RecordCount=est[0]->RecordCount;
    AverageCache=est[0]->AverageCache;
    for(int i=1; i<ThreadCount; i++)
    {
      AverageCache+=est[i]->AverageCache;
    }

    RealType wgt=1.0/static_cast<RealType>(ThreadCount);
    AverageCache*=wgt;

    //reset to 1 to indicate that scalars do not have the count
    HDFAttribIO<int> i(RecordCount,true);
    i.write(h_obs,"count");
    ThreadCount=1;
    HDFAttribIO<int> t(ThreadCount,true);
    t.write(h_obs,"threads");
    HDFAttribIO<Matrix<RealType> > m(AverageCache,true);
    m.write(h_obs,"scalars");

    stop();
  }

  /** Stop a run
   *
   * Collect data in Cache and print out data into hdf5 and ascii file.
   * This should not be called in a OpenMP parallel region or should
   * be guarded by master/single.
   * Keep the ascii output for now
   */
  void ScalarEstimatorManager::stop() 
  {
    if(CollectSum)
    {
      Matrix<RealType> c(AverageCache);
      myComm->reduce(c.data(),AverageCache.data(),AverageCache.size());
      RealType norm=1.0/static_cast<RealType>(myComm->ncontexts());
      AverageCache *= norm;
      if(Manager)
      {//write scalars_gsum 
        HDFAttribIO<Matrix<RealType> > m(AverageCache);
        m.write(h_obs,"scalars_gsum");
      }
    }

    //close the group
    if(h_obs>-1)
    {
      H5Gclose(h_obs); 
      h_obs=-1;
    }

    if(Manager) 
    {
      int nc=AverageCache.cols();
      EPSum[0]=0.0;
      EPSum[1]=RecordCount;
      const RealType* restrict eptr=AverageCache.data();
      for(int i=0; i<RecordCount; i++, eptr+=nc) EPSum[0] += *eptr;

      //This will go away with hdf5 record and utility
      string fname(RootName);
      fname.append(".scalar.dat");
      ofstream fout(fname.c_str());
      fout.setf(ios::scientific, ios::floatfield);
      const RealType* restrict rptr=AverageCache.data();
      const RealType* restrict pptr=PropertyCache.data();
      int pc=PropertyCache.cols();
      fout << "#   index ";
      for(int i=0; i<BlockAverages.size(); i++) fout << setw(16) << BlockAverages.Name[i];
      for(int i=0; i<BlockProperties.size(); i++) fout << setw(16) << BlockProperties.Name[i];
      fout << endl;
      fout.setf(ios::right,ios::adjustfield);
      for(int i=0; i<RecordCount; i++)
      {
        fout << setw(10) << i;
        for(int j=0; j<nc; j++) fout << setw(16) << *rptr++;
        for(int j=0; j<pc; j++) fout << setw(16) << *pptr++;
        fout << endl;
      }
    }

    Collected=true;
  }

  void ScalarEstimatorManager::stopBlock(RealType accept)
  {
    PropertyCache(RecordCount,weightInd) = Estimators[0]->d_wgt;
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
  }

  void ScalarEstimatorManager::stopBlock(const vector<ScalarEstimatorManager*> est)
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
      FastAccumulate(est[i]->AverageCache[rc], est[i]->AverageCache[rc+1],AverageCache[rc]);
    }

    if(h_obs<-1) return;

    PropertyCache=est[0]->PropertyCache;
    HDFAttribIO<int> i(RecordCount,true);
    i.write(h_obs,"count");
    HDFAttribIO<int> t(ThreadCount,true);
    t.write(h_obs,"threads");
    HDFAttribIO<Matrix<RealType> > m(AverageCache,true);
    m.write(h_obs,"scalars");
  }

  void ScalarEstimatorManager::accumulate(MCWalkerConfiguration::iterator it,
      MCWalkerConfiguration::iterator it_end)
  {
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->accumulate(it,it_end);
  }

  void ScalarEstimatorManager::accumulate(ParticleSet& P, 
      MCWalkerConfiguration::Walker_t& awalker) 
  {
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate(P,awalker);
  }

  void 
    ScalarEstimatorManager::getEnergyAndWeight(RealType& e, RealType& w) 
  {
#if defined(HAVE_MPI)
    MPI_Bcast(EPSum.begin(),2,MPI_DOUBLE,0,myComm->getMPI());
#endif
    e=EPSum[0];
    w=EPSum[1];
  }

  ScalarEstimatorManager::EstimatorType* 
    ScalarEstimatorManager::getMainEstimator() 
    {
      if(MainEstimator==0) 
        add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
      return MainEstimator;
    }

  ScalarEstimatorManager::EstimatorType* 
    ScalarEstimatorManager::getEstimator(const string& a) 
    {
      std::map<string,int>::iterator it = EstimatorMap.find(a);
      if(it == EstimatorMap.end()) 
        return 0;
      else 
        return Estimators[(*it).second];
    }

  /** This should be moved to branch engine */
  bool ScalarEstimatorManager::put(xmlNodePtr cur) 
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
            extra.push_back(cname);
        }
      } 
      cur = cur->next;
    }

    if(Estimators.empty()) 
    {
      app_log() << "  Using a default LocalEnergyOnlyEstimator for the MainEstimator " << endl;
      add(new LocalEnergyOnlyEstimator(),MainEstimatorName);
    } 

    //add extra
    return true;
  }

  int ScalarEstimatorManager::add(EstimatorType* newestimator, const string& aname) 
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
      app_log() << "  ScalarEstimatorManager::add replace " << aname << " estimator." << endl;
      delete Estimators[n]; 
      Estimators[n]=newestimator;
    }

    //check the name and set the MainEstimator
    if(aname == MainEstimatorName) MainEstimator=newestimator;
    return n;
  }

  int ScalarEstimatorManager::addObservable(const char* aname) 
  {
    int mine = BlockAverages.add(aname);
    int add = TotalAverages.add(aname);
    if(mine < Block2Total.size()) 
      Block2Total[mine] = add;
    else 
      Block2Total.push_back(add);
    return mine;
  }

  void ScalarEstimatorManager::getData(int i, vector<RealType>& values)
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
