//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by:  Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/ForwardWalking/FWSingle.h"
#include "QMCDrivers/ForwardWalkingStructure.h"

namespace qmcplusplus
{

/// Constructor.
FWSingle::FWSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool)
  : QMCDriver(w,psi,h,ppool), weightFreq(1), weightLength(5), fname(""), verbose(0)
  , WIDstring("WalkerID"), PIDstring("ParentID"),gensTransferred(1),startStep(0)
{
  RootName = "FW";
  QMCType ="FWSingle";
  xmlrootName="";
  m_param.add(xmlrootName,"rootname","string");
  m_param.add(weightLength,"numbersteps","int");
  m_param.add(weightFreq,"skipsteps","int");
  m_param.add(verbose,"verbose","int");
  m_param.add(startStep,"ignore","int");
}

bool FWSingle::run()
{
  Estimators->start(weightLength,1);
  fillIDMatrix();
  //we do this once because we only want to link parents to parents if we need to
  //     if (verbose>1) app_log()<<" getting weights for generation "<<gensTransferred<< std::endl;
  std::vector<std::vector<vector<int> > > WeightHistory;
  WeightHistory.push_back(Weights);
  for(int ill=1; ill<weightLength; ill++)
  {
    transferParentsOneGeneration();
    FWOneStep();
    WeightHistory.push_back(Weights);
  }
  if (verbose>0)
    app_log()<<" Done Computing Weights"<< std::endl;
  int nprops = H.sizeOfObservables();//local energy, local potnetial and all hamiltonian elements
  int FirstHamiltonian = H.startIndex();
  std::vector<std::vector<vector<RealType> > > savedValues;
  int nelectrons = W[0]->R.size();
  int nfloats=OHMMS_DIM*nelectrons;
  ForwardWalkingData fwer;
  fwer.resize(nelectrons);
  //      MCWalkerConfiguration* savedW = new MCWalkerConfiguration(W);
  for(int step=0; step<numSteps; step++)
  {
    std::vector<float> ALLcoordinates;
    readInFloat(step,ALLcoordinates);
    std::vector<float> SINGLEcoordinate(nfloats);
    std::vector<float>::iterator fgg(ALLcoordinates.begin()), fgg2(ALLcoordinates.begin()+nfloats);
    W.resetCollectables();
    std::vector<std::vector<RealType> > stepObservables;
    for(int wstep=0; wstep<walkersPerBlock[step]; wstep++)
    {
      copy( fgg,fgg2,SINGLEcoordinate.begin());
      fwer.fromFloat(SINGLEcoordinate);
      W.R=fwer.Pos;
      fgg+=nfloats;
      fgg2+=nfloats;
      W.update();
      RealType logpsi(Psi.evaluateLog(W));
      RealType eloc=H.evaluate( W );
      //             (*W[0]).resetProperty(logpsi,1,eloc);
      H.auxHevaluate(W);
      H.saveProperty(W.getPropertyBase());
      std::vector<RealType> walkerObservables(nprops+2,0);
      walkerObservables[0]= eloc;
      walkerObservables[1]= H.getLocalPotential();
      const RealType* restrict ePtr = W.getPropertyBase();
      for(int i=0; i<nprops; i++)
        walkerObservables[i+2] = ePtr[FirstHamiltonian+i] ;
      stepObservables.push_back(walkerObservables);
    }
    savedValues.push_back(stepObservables);
  }
  for(int ill=0; ill<weightLength; ill++)
  {
    Estimators->startBlock(1);
    Estimators->accumulate(savedValues,WeightHistory[ill],getNumberOfSamples(ill));
    Estimators->stopBlock(getNumberOfSamples(ill));
  }
  Estimators->stop();
  return true;
}

int FWSingle::getNumberOfSamples(int omittedSteps)
{
  int returnValue(0);
  for(int i=startStep; i<(numSteps-omittedSteps); i++)
    returnValue +=walkersPerBlock[i];
  return returnValue;
}

void FWSingle::readInLong(int step, std::string IDstring, std::vector<long>& data_out)
{
  int RANK=1;
  hid_t       dataset;
  hid_t       filespace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[1];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  herr_t      status, status_n;
  int         rank, rank_chunk;
  hsize_t hi, hj;
  c_file = H5Fopen(fname.str().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  std::stringstream gname("");
  gname<<"Block_"<< step;
  hid_t d_file = H5Gopen(c_file,gname.str().c_str());
  dataset = H5Dopen(d_file, IDstring.c_str());
  filespace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(filespace);
  status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  //     if (verbose>1) printf("dataset ",IDstring.c_str(),"rank %d, dimensions %lu x %lu\n", rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));
  data_out.resize(dims[0]);
  cparms = H5Dget_create_plist(dataset);
  if (H5D_CHUNKED == H5Pget_layout(cparms))
  {
    rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    //       if (verbose>1) printf("chunk rank %d, dimensions %lu \n", rank_chunk, (unsigned long)(chunk_dims[0]) );
  }
  memspace = H5Screate_simple(RANK,dims,NULL);
  status = H5Dread(dataset, H5T_NATIVE_LONG, memspace, filespace, H5P_DEFAULT, &(data_out[0]));
  //     if(verbose>2)
  //     {
  //       printf("\n");
  //       app_log()<<IDstring.c_str()<<" Dataset: \n"<< std::endl;
  //       for (int j = 0; j < dims[0]; j++) app_log()<<data_out[j]<<" ";
  //       app_log()<< std::endl;
  //     }
  H5Pclose(cparms);
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Gclose(d_file);
  H5Fclose(c_file);
}

void FWSingle::readInFloat(int step, std::vector<float>& data_out)
{
  int RANK=1;
  hid_t       dataset;
  hid_t       filespace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[2];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  hsize_t     col_dims[1];
  hsize_t     count[2];
  hsize_t     offset[2];
  herr_t      status, status_n;
  int         rank, rank_chunk;
  hsize_t hi, hj;
  c_file = H5Fopen(fname.str().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  std::stringstream gname("");
  gname<<"Block_"<< step;
  hid_t d_file = H5Gopen(c_file,gname.str().c_str());
  dataset = H5Dopen(d_file, "Positions");
  filespace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(filespace);
  status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  //     if (verbose>1) printf("dataset rank %d, dimensions %lu x %lu\n", rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));
  data_out.resize(dims[0]);
  cparms = H5Dget_create_plist(dataset);
  if (H5D_CHUNKED == H5Pget_layout(cparms))
  {
    rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    //       if (verbose>1) printf("chunk rank %d, dimensions %lu \n", rank_chunk, (unsigned long)(chunk_dims[0]) );
  }
  memspace = H5Screate_simple(RANK,dims,NULL);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, &(data_out[0]));
  if(verbose>2)
  {
    printf("\n");
    printf("Dataset: \n");
    for (int j = 0; j < dims[0]; j++)
      app_log()<<data_out[j]<<" ";
    app_log()<< std::endl;
  }
  H5Pclose(cparms);
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Gclose(d_file);
  H5Fclose(c_file);
}

void FWSingle::fillIDMatrix()
{
  if (verbose>0)
    app_log()<<" There are "<<numSteps<<" steps"<< std::endl;
  IDs.resize(numSteps);
  PIDs.resize(numSteps);
  Weights.resize(numSteps);
  std::vector<std::vector<long> >::iterator stepIDIterator(IDs.begin());
  std::vector<std::vector<long> >::iterator stepPIDIterator(PIDs.begin());
  int st(0);
  do
  {
    readInLong(st ,WIDstring ,*(stepIDIterator));
    readInLong(st ,PIDstring ,*(stepPIDIterator));
    walkersPerBlock.push_back( (*stepIDIterator).size() );
    stepIDIterator++;
    stepPIDIterator++;
    st++;
    if (verbose>1)
      app_log()<<"step:"<<st<< std::endl;
  }
  while (st<numSteps);
  //     Weights.resize( IDs.size());
  for(int i=0; i<numSteps; i++)
    Weights[i].resize(IDs[i].size(),1);
  realPIDs = PIDs;
  realIDs = IDs;
}

void FWSingle::fillWalkerPositionsandWeights(int nstep)
{
  //     int needed = walkersPerBlock[nstep] - W.getActiveWalkers();
  //     if (needed>0) W.createWalkers(needed);
  //     else if (needed<0) W.destroyWalkers(-1*needed-1);
  std::vector<float> ALLcoordinates;
  readInFloat(nstep,ALLcoordinates);
  int nelectrons = W[0]->R.size();
  int nfloats=OHMMS_DIM*nelectrons;
  std::vector<float> SINGLEcoordinate(nfloats);
  ForwardWalkingData fwer;
  fwer.resize(nelectrons);
  std::vector<float>::iterator fgg(ALLcoordinates.begin()), fgg2(ALLcoordinates.begin()+nfloats);
  for(int nw=0; nw<W.getActiveWalkers(); nw++)
  {
    copy( fgg,fgg2,SINGLEcoordinate.begin());
    fwer.fromFloat(SINGLEcoordinate);
    W[nw]->R=fwer.Pos;
    W[nw]->Weight = 1.0;
    fgg+=nfloats;
    fgg2+=nfloats;
  }
}

void FWSingle::resetWeights()
{
  if (verbose>2)
    app_log()<<" Resetting Weights"<< std::endl;
  Weights.clear();
  Weights.resize(numSteps);
  for(int i=0; i<numSteps; i++)
    Weights[i].resize(IDs[i].size(),0);
}

void FWSingle::FWOneStep()
{
  //create an ordered version of the ParentIDs to make the weight calculation faster
  std::vector<std::vector<long> > orderedPIDs=(PIDs);
  std::vector<std::vector<long> >::iterator it(orderedPIDs.begin());
  do
  {
    std::sort((*it).begin(),(*it).end());
    it++;
  }
  while(it<orderedPIDs.end());
  if (verbose>2)
    app_log()<<" Done Sorting IDs"<< std::endl;
  resetWeights();
  std::vector<std::vector<long> >::iterator stepIDIterator(IDs.begin());
  std::vector<std::vector<long> >::iterator stepPIDIterator(orderedPIDs.begin() + gensTransferred);
  std::vector<std::vector<int> >::iterator stepWeightsIterator(Weights.begin());
  //we start comparing the next generations ParentIDs with the current generations IDs
  int i=0;
  do
  {
    if (verbose>2)
      app_log()<<"  calculating weights for gen:"<<gensTransferred<<" step:"<<i<<"/"<<orderedPIDs.size()<< std::endl;
    //       if (verbose>2) app_log()<<"Nsamples ="<<(*stepWeightsIteratoetWeights).size()<< std::endl;
    std::vector<long>::iterator IDit( (*stepIDIterator).begin()     );
    std::vector<long>::iterator PIDit( (*stepPIDIterator).begin()   );
    std::vector<int>::iterator  Wit( (*stepWeightsIterator).begin() );
    if (verbose>2)
      app_log()<<"ID size:"<<(*stepIDIterator).size()<<" PID size:"<<(*stepPIDIterator).size()<<" Weight size:"<<(*stepWeightsIterator).size()<< std::endl;
    do
    {
      if ((*PIDit)==(*IDit))
      {
        (*Wit)++;
        PIDit++;
      }
      else
      {
        IDit++;
        Wit++;
        if (IDit==(*stepIDIterator).end())
        {
          IDit=(*stepIDIterator).begin();
          Wit=(*stepWeightsIterator).begin();
        }
      }
    }
    while(PIDit<(*stepPIDIterator).end());
    //       if (verbose>2) { printIDs((*stepIDIterator));printIDs((*stepPIDIterator));}
    //       if (verbose>2) printInts((*stepWeightsIterator));
    stepIDIterator++;
    stepPIDIterator++;
    stepWeightsIterator++;
    i++;
  }
  while(stepPIDIterator<orderedPIDs.end());
}

void FWSingle::printIDs(std::vector<long> vi)
{
  for (int j=0; j<vi.size(); j++)
    app_log()<<vi[j]<<" ";
  app_log()<< std::endl;
}
void FWSingle::printInts(std::vector<int> vi)
{
  for (int j=0; j<vi.size(); j++)
    app_log()<<vi[j]<<" ";
  app_log()<< std::endl;
}

void FWSingle::transferParentsOneGeneration( )
{
  std::vector<std::vector<long> >::reverse_iterator stepIDIterator(IDs.rbegin());
  std::vector<std::vector<long> >::reverse_iterator stepPIDIterator(PIDs.rbegin()), nextStepPIDIterator(realPIDs.rbegin());
  stepIDIterator+=gensTransferred;
  nextStepPIDIterator+=gensTransferred;
  int i(0);
  do
  {
    std::vector<long>::iterator hereID( (*stepIDIterator).begin() ) ;
    std::vector<long>::iterator nextStepPID( (*nextStepPIDIterator).begin() );
    std::vector<long>::iterator herePID( (*stepPIDIterator).begin() );
    if (verbose>2)
      app_log()<<"  calculating Parent IDs for gen:"<<gensTransferred<<" step:"<<i<<"/"<<PIDs.size()-gensTransferred<< std::endl;
    if (verbose>2)
    {
      printIDs((*nextStepPIDIterator));
      printIDs((*stepIDIterator));
      printIDs((*stepPIDIterator));
    }
    do
    {
      if ((*herePID)==(*hereID))
      {
        (*herePID)=(*nextStepPID);
        herePID++;
      }
      else
      {
        hereID++;
        nextStepPID++;
        if (hereID==(*stepIDIterator).end())
        {
          hereID=(*stepIDIterator).begin();
          nextStepPID=(*nextStepPIDIterator).begin();
          //             if (verbose>2) app_log()<<"resetting to beginning of parents"<< std::endl;
        }
      }
    }
    while(herePID<(*stepPIDIterator).end());
    stepIDIterator++;
    nextStepPIDIterator++;
    stepPIDIterator++;
    i++;
  }
  while(stepIDIterator<IDs.rend());
  gensTransferred++; //number of gens backward to compare parents
  if (verbose>2)
    app_log()<<"  Finished generation block"<< std::endl;
}


bool
FWSingle::put(xmlNodePtr q)
{
  fname<<xmlrootName<<".storeConfig.h5";
  c_file = H5Fopen(fname.str().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  hsize_t numGrps = 0;
  H5Gget_num_objs(c_file, &numGrps);
  numSteps = static_cast<int> (numGrps)-3;
  app_log()<<"  Total number of steps in input file "<<numSteps<< std::endl;
  if (weightFreq<1)
    weightFreq=1;
  int numberDataPoints = weightLength/weightFreq;
  //     pointsToCalculate.resize(numberDataPoints);
  //     for(int i=0;i<numberDataPoints;i++) pointsToCalculate[i]=i*weightFreq;
  app_log()<<"  "<<numberDataPoints<<" observables will be calculated each "<<weightFreq<<" steps"<< std::endl;
  app_log()<<"  Config Generations skipped for thermalization: "<<startStep<< std::endl;//<<" steps. At: ";
  //     for(int i=0;i<numberDataPoints;i++) app_log()<<pointsToCalculate[i]<<" ";
  app_log()<< std::endl;
  if (H5Fclose(c_file)>-1)
    c_file=-1;
  return true;
}
}

