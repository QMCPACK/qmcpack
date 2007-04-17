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

  //initialize the name of the primary estimator

  ScalarEstimatorManager::ScalarEstimatorManager(QMCHamiltonian& h, Communicate* c): 
    MainEstimatorName("elocal"),
  Manager(false), CollectSum(true), AppendRecord(false), Collected(false), 
  h_file(-1), h_obs(-1), myComm(0), OutStream(0), H(h),MainEstimator(0)
  { 
    //Block2Total.resize(0); 
  }

  ScalarEstimatorManager::~ScalarEstimatorManager()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
    if(OutStream) delete OutStream;
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

  ////reset when a new driver is used
  //void ScalarEstimatorManager::reset(const string& aname, bool append)
  //{
  //  RootName = aname;
  //  AppendRecord=append;
  //}

  void ScalarEstimatorManager::reset()
  {

    if(Estimators.empty()) 
      add(new LocalEnergyOnlyEstimator(),MainEstimatorName);

    //no need for RemoteData. We don't collect anything during a run
    BufferType tmp;
    for(int i=0; i<Estimators.size(); i++) 
      Estimators[i]->add2Record(BlockAverages,tmp);

    MyIndex[WEIGHT_INDEX] = BlockAverages.add("WeightSum");
    MyIndex[BLOCK_CPU_INDEX] = BlockAverages.add("BlockCPU");
    MyIndex[ACCEPT_RATIO_INDEX] = BlockAverages.add("AcceptRatio");
  }

  void ScalarEstimatorManager::start(int blocks)
  {
    Collected= false;
    RecordCount=0;
    BlockAverages.setValues(0.0);
    Cache.resize(blocks,BlockAverages.size());

    //open obsevables and allocate spaces 
    string h5file=RootName+".config.h5";
    h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    h_obs = H5Gcreate(h_file,"observables",0);
    HDFAttribIO<int> i(RecordCount);
    i.write(h_obs,"count");
    HDFAttribIO<Matrix<RealType> > m(Cache);
    m.write(h_obs,"scalars");

    ostringstream o;
    for(int i=0; i<BlockAverages.size()-1; i++) o << BlockAverages.Name[i] << ":";
    o<<BlockAverages.Name.back();
    string t(o.str());
    HDFAttribIO<string> so(t);
    so.write(h_obs,"scalar_ids");

    //H5Gclose(h_obs);
    H5Fclose(h_file);
  }

  /** Stop a run
   *
   * Collect data in Cache and print out
   * Keep the ascii output for now
   */
  void ScalarEstimatorManager::stop() 
  {
    if(CollectSum)
    {
      Matrix<RealType> c(Cache);
      myComm->reduce(c.data(),Cache.data(),Cache.size());
      RealType norm=1.0/static_cast<RealType>(myComm->ncontexts());
      Cache *= norm;
      if(Manager)
      {//write scalars_gsum 
        HDFAttribIO<Matrix<RealType> > m(Cache);
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
      int nc=Cache.cols();
      EPSum[0]=0.0;
      EPSum[1]=RecordCount;
      const RealType* restrict eptr=Cache.data()+BlockAverages.add("LocalEnergy");
      for(int i=0; i<RecordCount; i++, eptr+=nc) EPSum[0] += *eptr;

      //This will go away with hdf5 record and utility
      string fname(RootName);
      fname.append(".scalar.dat");
      ofstream fout(fname.c_str());
      fout.setf(ios::scientific, ios::floatfield);
      const RealType* restrict rptr=Cache.data();
      fout << "#   index ";
      for(int i=0; i<BlockAverages.size(); i++) fout << setw(16) << BlockAverages.Name[i];
      fout << endl;
      fout.setf(ios::right,ios::adjustfield);
      for(int i=0; i<RecordCount; i++)
      {
        fout << setw(10) << i;
        for(int j=0; j<nc; j++) fout << setw(16) << *rptr++;
        fout << endl;
      }
    }

    Collected=true;
  }


  void ScalarEstimatorManager::stopBlock(RealType accept)
  {
    BlockAverages[MyIndex[WEIGHT_INDEX]]=MyData[WEIGHT_INDEX];
    BlockAverages[MyIndex[BLOCK_CPU_INDEX]] = MyTimer.elapsed();
    BlockAverages[MyIndex[ACCEPT_RATIO_INDEX]] = accept;
    RealType wgtinv = 1.0/MyData[WEIGHT_INDEX];
    for(int i=0; i<Estimators.size(); i++) 
      Estimators[i]->report(BlockAverages,wgtinv);
    MyData=0.0;
    std::copy(BlockAverages.begin(),BlockAverages.end(),Cache[RecordCount]);
    RecordCount++;

    HDFAttribIO<int> i(RecordCount,true);
    i.write(h_obs,"count");
    HDFAttribIO<Matrix<RealType> > m(Cache,true);
    m.write(h_obs,"scalars");
  }

  /** accumulate data for all the estimators
  */
  void ScalarEstimatorManager::accumulate(MCWalkerConfiguration& W) 
  {
    RealType wgt=0.0;
    for(int i=0; i< Estimators.size(); i++) 
      wgt += Estimators[i]->accumulate(W.begin(),W.end());
    //increment BinSize
    MyData[WEIGHT_INDEX]+=wgt;
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
