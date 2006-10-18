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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
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
#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#include "Estimators/PolarizationEstimator.h"
#include "Utilities/IteratorUtility.h"
//#define QMC_ASYNC_COLLECT

using namespace qmcplusplus;

ScalarEstimatorManager::ScalarEstimatorManager(QMCHamiltonian& h): 
 FileManager(true), CollectSum(false), 
 AccumulateBlocks(false), ManagerNode(false),
 firstReport(true), 
 Period(1000), NodeWeight(1.0), H(h), RootName("estimator"), 
  OutStream(0){ 
    Block2Total.resize(0); 
    myNodeID=OHMMS::Controller->mycontext();
    numNodes=OHMMS::Controller->getNumNodes();
    myGroupID=0;
    if(numNodes ==1) {
      RemoteData.push_back(new BufferType);
      myRequest.resize(1);
      myStatus.resize(1);
    } else {
      if(myNodeID) {
        RemoteData.push_back(new BufferType);
        RemoteData.push_back(new BufferType);
        myRequest.resize(1);
        myStatus.resize(1);
      } else {
        for(int i=0; i<numNodes;i++) RemoteData.push_back(new BufferType);
        myRequest.resize(numNodes-1);
        myStatus.resize(numNodes-1);
      }
    }
  }

ScalarEstimatorManager::~ScalarEstimatorManager(){ 
  delete_iter(Estimators.begin(), Estimators.end());
  delete_iter(RemoteData.begin(), RemoteData.end());
  if(OutStream) delete OutStream;
}

int
ScalarEstimatorManager::add(EstimatorType* newestimator, const string& aname) { 

  std::map<string,int>::iterator it = EstimatorMap.find(aname);
  if(it == EstimatorMap.end()) {
    int n =  Estimators.size();
    Estimators.push_back(newestimator);
    EstimatorMap[aname] = n;
    return n;
  } else {
    return (*it).second;
  }
}

/** reset the internal data of all the estimators for new averages
 */
void 
ScalarEstimatorManager::reset() {
  BinSize=0;
  MyData=0.0;
  for(int i=0; i< Estimators.size(); i++) Estimators[i]->reset();
}

/**  print the header to the output file
 *
 * When ColllectSum is true, the manager node initiates asychronous recv.
 */
void 
ScalarEstimatorManager::reportHeader(bool append) {
  EPSum=0.0;
  if(FileManager) {
    if(!append)  {
      *OutStream << "#    index     ";
      for(int i=0; i<BlockAverages.size(); i++) 
        (*OutStream) << setw(16) << BlockAverages.Name[i];
      (*OutStream) << endl;
    }
    OutStream->setf(ios::right,ios::adjustfield);
  }

#if defined(HAVE_MPI) && defined(QMC_ASYNC_COLLECT)
  if(ManagerNode) { //initiate the recv
    for(int i=1,is=0; i<numNodes; i++,is++) {
      MPI_Irecv(RemoteData[i]->begin(),msgBufferSize,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&(myRequest[is]));
    }
  }
#endif
}


void ScalarEstimatorManager::getEnergyAndWeight(RealType& e, RealType& w) {
  //EPSum is always valid for the manager node
#if defined(HAVE_MPI)
  MPI_Bcast(EPSum.begin(),2,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  e=EPSum[0];
  w=EPSum[1];
}

/** closes the stream to the output file
 */
void 
ScalarEstimatorManager::finalize() {
  if(OutStream) delete OutStream;
  OutStream = 0;
#if defined(HAVE_MPI) && defined(QMC_ASYNC_COLLECT)
  if(ManagerNode) {//cancel irecv initiazed by flush
    for(int is=0; is<numNodes-1; is++) MPI_Cancel(&myRequest[is]);
  }
#endif
}

//void 
//ScalarEstimatorManager::finalize() {
//  if(OutStream) delete OutStream;
//  OutStream = 0;
//#if defined(HAVE_MPI) 
//#if defined(QMC_ASYNC_COLLECT)
//  if(ManagerNode) {//cancel irecv initiazed by flush
//    for(int is=0; is<numNodes-1; is++) MPI_Cancel(&myRequest[is]);
//  }
//#endif
//  MPI_Bcast(EPSum.begin(),2,MPI_DOUBLE,0,MPI_COMM_WORLD);
//#endif
//}
///** closes the stream to the output file
// */
//void 
//ScalarEstimatorManager::finalize(SimpleFixedNodeBranch& branchEngine) {
//  if(OutStream) delete OutStream;
//  OutStream = 0;
//#if defined(HAVE_MPI) 
//#if defined(QMC_ASYNC_COLLECT)
//  if(ManagerNode) {//cancel irecv initiazed by flush
//    for(int is=0; is<numNodes-1; is++) MPI_Cancel(&myRequest[is]);
//  }
//#endif
//  MPI_Bcast(EPSum.begin(),2,MPI_DOUBLE,0,MPI_COMM_WORLD);
//#endif
//  branchEngine.setTrialEnergy(EPSum[0],EPSum[1]);
//}
//


/** accumulate data for all the estimators
 */
void 
ScalarEstimatorManager::accumulate(MCWalkerConfiguration& W) {

  RealType wgt=0.0;
  for(int i=0; i< Estimators.size(); i++) 
    wgt += Estimators[i]->accumulate(W.begin(),W.end());
  //increment BinSize
  BinSize++;
  MyData[WEIGHT_INDEX]+=wgt;
}


/**  compute the averages for all the estimators and reset
 *
 * @todo MPI calls will be rmoved.
 */
void ScalarEstimatorManager::flush(){

#if defined(HAVE_MPI) 
  if(CollectSum) {
    RemoteData[0]->rewind();
    RemoteData[0]->put(MyData.begin(),MyData.end());
    for(int i=0; i<Estimators.size(); i++) {
      Estimators[i]->copy2Buffer(*RemoteData[0]);
    }

#if defined(QMC_ASYNC_COLLECT)
    if(myNodeID) {
      MPI_Isend(RemoteData[0]->data(),msgBufferSize, MPI_DOUBLE,0,myNodeID,MPI_COMM_WORLD,&(myRequest[0]));
    } else {
      MPI_Waitall(numNodes-1,&(myRequest[0]),&(myStatus[0]));
      for(int is=1; is<numNodes; is++) {
        BufferType::iterator tit(RemoteData[0]->begin());
        BufferType::const_iterator rit(RemoteData[is]->begin()), rit_end(RemoteData[is]->end());
        while(rit != rit_end) (*tit++) += (*rit++);
      }
    }

    RemoteData[0]->rewind();
    RemoteData[0]->get(MyData.begin(),MyData.end());
    RealType wgtinv = static_cast<RealType>(numNodes)/MyData[WEIGHT_INDEX];
    for(int i=0; i<Estimators.size(); i++) 
      Estimators[i]->report(BlockAverages,wgtinv,*RemoteData[0]);
#else
    //using reduce function
    MPI_Reduce(RemoteData[0]->data(),RemoteData[1]->data(),RemoteData[0]->size(),
        MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    RemoteData[1]->rewind();
    RemoteData[1]->get(MyData.begin(),MyData.end());
    RealType wgtinv = static_cast<RealType>(numNodes)/MyData[WEIGHT_INDEX];
    for(int i=0; i<Estimators.size(); i++) 
      Estimators[i]->report(BlockAverages,wgtinv,*RemoteData[1]);
#endif
  } else {//CollectSum == false
#endif
    RealType wgtinv = 1.0/MyData[WEIGHT_INDEX];
    for(int i=0; i<Estimators.size(); i++) 
      Estimators[i]->report(BlockAverages,wgtinv,*RemoteData[0]);

#if defined(HAVE_MPI) 
  }
#endif

    
  BlockAverages[MyIndex[WEIGHT_INDEX]]=MyData[WEIGHT_INDEX];
  BlockAverages[MyIndex[BLOCK_CPU_INDEX]] = MyData[BLOCK_CPU_INDEX]*NodeWeight;
  BlockAverages[MyIndex[ACCEPT_RATIO_INDEX]] = MyData[ACCEPT_RATIO_INDEX]*NodeWeight;

  BinSize=0;
  MyData=0.0;
  EPSum[0]+= Estimators[0]->average();
  EPSum[1]+= 1.0;

#if defined(HAVE_MPI) && defined(QMC_ASYNC_COLLECT)
  if(ManagerNode) {
    for(int i=1,is=0; i<numNodes; i++,is++) {
      MPI_Irecv(RemoteData[i]->begin(),msgBufferSize,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&(myRequest[is]));
    }
  }
#endif
}


/** print the averages for all the estimators to a file
 * @param iter the interval 
 */
void ScalarEstimatorManager::report(int iter){
  if(AccumulateBlocks) {
    if(firstReport){
      index = 0;
      TotalAveragesData.resize(0,TotalAverages.size());
      firstReport = false;
    }
    TotalAveragesData.add(1);
    for(int i=0; i<Block2Total.size(); i++) {
      TotalAveragesData(index,Block2Total[i])=BlockAverages[i];
    }
    index++;
  }// else {
  if(FileManager) {
    (*OutStream) << setw(10) << iter;
    for(int i=0; i<BlockAverages.size();i++)
      (*OutStream) << setw(16) << BlockAverages[i];
    (*OutStream) << endl;
  }
  //}
}


int ScalarEstimatorManager::addObservable(const char* aname) {
  int mine = BlockAverages.add(aname);
  int add = TotalAverages.add(aname);
  if(mine < Block2Total.size()) {
    Block2Total[mine] = add;
  } else {
    Block2Total.push_back(add);
  }
  return mine;
}

void ScalarEstimatorManager::getData(int i, vector<RealType>& values){
  int entries = TotalAveragesData.rows();
  values.resize(entries);
  for (int a=0; a<entries; a++){	
    values[a] = TotalAveragesData(a,Block2Total[i]);
  }
}

/** combines the functionality of flush and report
 * @param iter the interval
 */
void ScalarEstimatorManager::flushreport(int iter){
  flush();
  if(FileManager) {
    (*OutStream) << setw(10) << iter;
    for(int i=0; i<BlockAverages.size();i++) (*OutStream) << setw(16) << BlockAverages[i];
    (*OutStream) << endl;
  }
}

/** set CollectSum
 * @param collect if true, global sum is done over the values
 *
 * FileManager is set when collect is true so that only the first node writes.
 */
void ScalarEstimatorManager::setCollectionMode(bool collect) {
  CollectSum=collect;
  for(int i=0; i< Estimators.size(); i++) Estimators[i]->CollectSum = collect;
  if(collect) {
    FileManager = OHMMS::Controller->master();
    NodeWeight = 1.0/static_cast<RealType>(OHMMS::Controller->ncontexts());
  } else {
    FileManager = true;
    NodeWeight = 1.0;
  }
  ManagerNode = (CollectSum && myNodeID==0);
}


void 
ScalarEstimatorManager::resetReportSettings(const string& aname, bool append) {

  //at least have local energy
  if(Estimators.empty()) {
    add(new LocalEnergyEstimator(H),"elocal");
  } 
  
  RemoteData[0]->clear();
  RemoteData[0]->rewind();
  RemoteData[0]->add(MyData.begin(),MyData.end());
  //update the weight index table and set up the buffer
  for(int i=0; i<Estimators.size(); i++) {
    Estimators[i]->add2Record(BlockAverages,*RemoteData[0]);
  }
  msgBufferSize=RemoteData[0]->current();

  //copy extra remote data
  for(int i=1;i<RemoteData.size(); i++) (*RemoteData[i]) = (*RemoteData[0]);

  MyIndex[WEIGHT_INDEX] = BlockAverages.add("WeightSum");
  MyIndex[BLOCK_CPU_INDEX] = BlockAverages.add("BlockCPU");
  MyIndex[ACCEPT_RATIO_INDEX] = BlockAverages.add("AcceptRatio");

  RootName = aname;

  if(FileManager) {
    string fname(aname);
    fname.append(".scalar.dat");
    if(OutStream) {delete OutStream; OutStream=0;}

    if(append) 
      OutStream = new ofstream(fname.c_str(), ios::app);
    else
      OutStream = new ofstream(fname.c_str());

    OutStream->setf(ios::scientific, ios::floatfield);
    OutStream->setf(ios::left,ios::adjustfield);
  }

  BlockAverages.setValues(0.0);
  BinSize=0;

}


ScalarEstimatorManager::EstimatorType* 
ScalarEstimatorManager::getMainEstimator() {
  return getEstimator("elocal");
}

ScalarEstimatorManager::EstimatorType* 
ScalarEstimatorManager::getEstimator(const string& a) {
  std::map<string,int>::iterator it = EstimatorMap.find(a);
  if(it == EstimatorMap.end()) {
    return 0;
  } else {
    return Estimators[(*it).second];
  }
}

bool ScalarEstimatorManager::put(xmlNodePtr cur) {

  vector<string> extra;

  cur = cur->children;

  //count the number of Hamiltonian copies
  int nsystem=0;
  while(cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname == "qmcsystem") {
      nsystem++;
    } else if(cname == "estimator") {
      xmlChar* att=xmlGetProp(cur,(const xmlChar*)"name");
      if(att) {
	string aname((const char*)att);
	if(aname == "LocalEnergy") {
          //att=xmlGetProp(cur,(const xmlChar*)"size");
	  //int ncopy(1);
          //if(att) {ncopy=atoi((const char*)att);}
          //add(new LocalEnergyEstimator<RealType>(H,ncopy),"elocal");
          add(new LocalEnergyEstimator(H),"elocal");
	} else { 
	  extra.push_back(cname);
	}
      }
    } 
    cur = cur->next;
  }

  if(Estimators.empty()) {
    //if(nsystem == 0) nsystem=1;
    add(new LocalEnergyEstimator(H),"elocal");
  } 

  //for(int i=0; i<extra.size(); i++) {
  //  if(extra[i] == "Polarization"){
  //    add(new PolarizationEstimator<RealType>(),extra[i]);
  //  }
  //}

  //add extra
  return true;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
