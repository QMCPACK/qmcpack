//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim  and Jordan Vincent
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
#include "Estimators/PolarizationEstimator.h"

using namespace ohmmsqmc;

ScalarEstimatorManager::ScalarEstimatorManager(QMCHamiltonian& h): 
  Stride(1000), WeightSum(0.0), H(h), RootName("estimator"), 
  OutStream(NULL) { }

ScalarEstimatorManager::~ScalarEstimatorManager(){ 
  Estimators.erase(Estimators.begin(), Estimators.end());
  if(OutStream) delete OutStream;
}

int
ScalarEstimatorManager::add(EstimatorType* newestimator, const string& aname) { 

  std::map<string,int>::iterator it = EstimatorMap.find(aname);
  if(it == EstimatorMap.end()) {
    int n =  Estimators.size();
    Estimators.push_back(newestimator);
    EstimatorMap[aname] = n;
    //newestimator->add2Record(BlockAverages);
    return n;
  } else {
    cout << "Already added estimator " << aname << endl;
    return (*it).second;
  }
}

/**
   @brief reset the internal data of all the estimators for new averages
*/
void 
ScalarEstimatorManager::reset() {
  WeightSum = 0.0;
  for(int i=0; i< Estimators.size(); i++) Estimators[i]->reset();

}

/**
   @brief accumulate data for all the estimators
*/
void 
ScalarEstimatorManager::accumulate(const MCWalkerConfiguration& W) {

  for(MCWalkerConfiguration::const_iterator it = W.begin(); 
      it != W.end(); it++){    
    RealType wgt = (*it)->Properties(WEIGHT);
    WeightSum += wgt;
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate(**it,wgt);
  }
}

/**
   @brief compute the averages for all the estimators and reset
*/
void ScalarEstimatorManager::flush(){
  RealType wgtinv = 1.0/WeightSum;
  for(int i=0; i<Estimators.size(); i++) 
    Estimators[i]->report(BlockAverages,wgtinv);

  BlockAverages[WeightIndex] = WeightSum;
  WeightSum = 0.0;
}

void ScalarEstimatorManager::setCollectionMode(bool collect) {
  for(int i=0; i< Estimators.size(); i++) Estimators[i]->CollectSum = collect;
}

/**
   @param iter the interval 
   @brief print the averages for all the estimators to a file
*/
void ScalarEstimatorManager::report(int iter){
  (*OutStream) << setw(10) << iter;
  for(int i=0; i<BlockAverages.size();i++)
    (*OutStream) << setw(16) << BlockAverages[i];
  (*OutStream) << endl;
}

/**
   @param iter the interval
   @brief combines the functionality of flush and report
*/
void ScalarEstimatorManager::flushreport(int iter){

  gsum(WeightSum,0);
  RealType wgtinv = 1.0/WeightSum;
  for(int i=0; i<Estimators.size(); i++)  Estimators[i]->report(BlockAverages,wgtinv);
  
  BlockAverages[WeightIndex] = WeightSum;
  WeightSum = 0.0;

  (*OutStream) << setw(10) << iter;
  for(int i=0; i<BlockAverages.size();i++) (*OutStream) << setw(16) << BlockAverages[i];
  (*OutStream) << endl;
}

void 
ScalarEstimatorManager::resetReportSettings(const string& aname) {

  //at least have local energy
  if(Estimators.empty()) {
    add(new LocalEnergyEstimator<RealType>(H),"elocal");
  } 

  //update the weight index
  for(int i=0; i<Estimators.size(); i++) Estimators[i]->add2Record(BlockAverages);
  WeightIndex = BlockAverages.add("WeightSum");
  if(aname != RootName) {
    RootName = aname;
    string fname(aname);
    fname.append(".scalar.dat");
    if(OutStream) delete OutStream;
    OutStream = new ofstream(fname.c_str());
    OutStream->setf(ios::scientific, ios::floatfield);
    OutStream->setf(ios::left,ios::adjustfield);
  }

  BlockAverages.setValues(0.0);
}

/**
   @brief print the header to the output file
*/
void 
ScalarEstimatorManager::reportHeader() {
  *OutStream << "#    index     ";
  for(int i=0; i<BlockAverages.size(); i++) 
    (*OutStream) << setw(16) << BlockAverages.Name[i];
  (*OutStream) << endl;
  OutStream->setf(ios::right,ios::adjustfield);
}

/**
   @brief closes the stream to the output file
*/
void 
ScalarEstimatorManager::finalize() {
  if(OutStream) delete OutStream;
  OutStream = NULL;
}

/**
   @param istride an inverval of reportting/flushing the cummulative quantities
   @brief set the stride of all the estimators to istride
*/
void 
ScalarEstimatorManager::setStride(int istride) {
  Stride = istride;
}

ScalarEstimatorManager::EstimatorType* 
ScalarEstimatorManager::getEstimator(const string& a) {
  std::map<string,int>::iterator it = EstimatorMap.find(a);
  if(it == EstimatorMap.end()) {
    return NULL;
  } else {
    return Estimators[(*it).second];
  }
}

bool ScalarEstimatorManager::put(xmlNodePtr cur) {

  vector<string> extra;
  cur = cur->children;
  while(cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname == "estimator") {
      xmlChar* att=xmlGetProp(cur,(const xmlChar*)"name");
      if(att) {
	string aname((const char*)att);
	if(aname == "LocalEnergy") {
          att=xmlGetProp(cur,(const xmlChar*)"size");
	  int ncopy(1);
          if(att) {ncopy=atoi((const char*)att);}
          add(new LocalEnergyEstimator<RealType>(H,ncopy),"elocal");
	} else { 
	  extra.push_back(cname);
	}
      }
    } 
    cur = cur->next;
  }

  if(Estimators.empty()) {
    add(new LocalEnergyEstimator<RealType>(H),"elocal");
  } 

  for(int i=0; i<extra.size(); i++) {
    if(extra[i] == "Polarization"){
      add(new PolarizationEstimator<RealType>(),extra[i]);
    }
  }

  //add extra
  return true;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
