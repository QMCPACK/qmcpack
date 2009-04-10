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
#include "QMCDrivers/ForwardWalking/FWSingle.h"
#include "QMCDrivers/ForwardWalkingStructure.h"

namespace qmcplusplus { 

  /// Constructor.
  FWSingle::FWSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), weightFreq(1), weightLength(5), fname(""), verbose(0), WID("WalkerID"), PID("ParentID"),gensTransferred(1)
  { 
    RootName = "FW";
    QMCType ="FWSingle";
    xmlrootName="";
    m_param.add(xmlrootName,"rootname","string");
    m_param.add(weightLength,"numbersteps","int");
    m_param.add(weightFreq,"skipsteps","int");
    m_param.add(verbose,"verbose","int");
  }
  
  bool FWSingle::run() {
    Estimators->start(weightLength,1);
    fillIDMatrix();
    vector<vector<vector<int> > > WeightHistory(weightLength);
    //we do this once because we only want to link parents to parents if we need to
    if (verbose>1) app_log()<<" getting weights for generation "<<gensTransferred<<endl;
    this->FWOneStep();
    WeightHistory.push_back(Weights);
    for(int ill=1;ill<weightLength;ill++)
    {
      Estimators->startBlock(1);
      if (verbose>1) app_log()<<endl<<" getting weights for generation "<<gensTransferred<<" "<<ill<<endl<<endl;
      this->transferParentsOneGeneration();
      this->FWOneStep();
      WeightHistory.push_back(Weights);
      //now the weights are correct. Next we need to load the coordinates into the configurations, evaluate the Hamiltonian,
      //weight them according to each entry in weighthistory, and accumulate statistics.
//       Estimators->startBlock(numSteps-ill);
      for(int step=0;step<(numSteps-ill);step++)
      {
        fillWalkerPositionsandWeights(step);
        for(int wstep=0;wstep<walkersPerBlock[step];wstep++)
        {
          W.resetCollectables();//do I need to do this?
          W.R = W[wstep]->R;
          W.update();
          return_t eloc=H.evaluate( W );
          (*W[wstep]).resetProperty(1,1,eloc);
          H.auxHevaluate(W,*(W[wstep]));
          H.saveProperty((*W[wstep]).getPropertyBase());
        }
        Estimators->accumulate(W);
      }
      Estimators->stopBlock(1);
//       recordBlock(ill);
    }

    
    Estimators->stop();
    return true;
  }
  
  void FWSingle::readInLong(int step, string IDstring, vector<long>& data_out)
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
    stringstream gname("");
        
    gname<<"Block_"<< step;
    hid_t d_file = H5Gopen(c_file,gname.str().c_str());
    dataset = H5Dopen(d_file, IDstring.c_str());
    filespace = H5Dget_space(dataset); 
    rank = H5Sget_simple_extent_ndims(filespace);
    status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);
    if (verbose>1) printf("dataset ",IDstring.c_str(),"rank %d, dimensions %lu x %lu\n", rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));
    data_out.resize(dims[0]);
    cparms = H5Dget_create_plist(dataset);
    if (H5D_CHUNKED == H5Pget_layout(cparms))  {
      rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
      if (verbose>1) printf("chunk rank %d, dimensions %lu \n", rank_chunk, (unsigned long)(chunk_dims[0]) );
    }
    memspace = H5Screate_simple(RANK,dims,NULL);
    
    status = H5Dread(dataset, H5T_NATIVE_LONG, memspace, filespace, H5P_DEFAULT, &(data_out[0]));
    if(verbose>2)
    {
      printf("\n");
      app_log()<<IDstring.c_str()<<" Dataset: \n"<<endl;
      for (int j = 0; j < dims[0]; j++) app_log()<<data_out[j]<<" ";
      app_log()<<endl;
    }
    H5Pclose(cparms);
    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);

    H5Gclose(d_file);
    H5Fclose(c_file);
  }
  
  void FWSingle::readInFloat(int step, vector<float>& data_out)
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
    stringstream gname("");
        
    gname<<"Block_"<< step;
    hid_t d_file = H5Gopen(c_file,gname.str().c_str());
    dataset = H5Dopen(d_file, "Positions");
    filespace = H5Dget_space(dataset); 
    rank = H5Sget_simple_extent_ndims(filespace);
    status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);
    if (verbose>1) printf("dataset rank %d, dimensions %lu x %lu\n", rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));
    data_out.resize(dims[0]);
    cparms = H5Dget_create_plist(dataset);
    if (H5D_CHUNKED == H5Pget_layout(cparms))  {
      rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
      if (verbose>1) printf("chunk rank %d, dimensions %lu \n", rank_chunk, (unsigned long)(chunk_dims[0]) );
    }
    memspace = H5Screate_simple(RANK,dims,NULL);
    
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, &(data_out[0]));
    if(verbose>2)
    {
      printf("\n");
      printf("Dataset: \n");
      for (int j = 0; j < dims[0]; j++) app_log()<<data_out[j]<<" ";
      app_log()<<endl;
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
    IDs.resize(numSteps); PIDs.resize(numSteps); Weights.resize(numSteps);
    vector<vector<long> >::iterator stepIDIterator(IDs.begin());
    vector<vector<long> >::iterator stepPIDIterator(PIDs.begin());
    int st(0);
    do
    {
      this->readInLong(st ,WID ,*(stepIDIterator));
      this->readInLong(st ,PID ,*(stepPIDIterator));
      walkersPerBlock.push_back( (*stepIDIterator).size() );
      stepIDIterator++; stepPIDIterator++; st++;
    } while (st<numSteps);
    this->resetWeights();
    realPIDs = PIDs;
  }
  
  void FWSingle::fillWalkerPositionsandWeights(int nstep)
  {
    int needed = walkersPerBlock[nstep] - W.getActiveWalkers();
    if (needed>0) W.createWalkers(needed);
    else if (needed<0) W.destroyWalkers(-1*needed);
    vector<float> ALLcoordinates;
    this->readInFloat(nstep,ALLcoordinates);
    int nelectrons = W[0]->R.size();
    int nfloats=OHMMS_DIM*nelectrons;
    vector<float> SINGLEcoordinate(nfloats);
    ForwardWalkingData* fwer = new ForwardWalkingData(nelectrons);
    for(int nw=0;nw<W.getActiveWalkers();nw++)
    {
      std::copy( &ALLcoordinates[nw*nfloats],&ALLcoordinates[(nw+1)*nfloats],&(SINGLEcoordinate[0]));
      fwer->fromFloat(SINGLEcoordinate);
      W[nw]->R=fwer->Pos;
      W[nw]->Weight = 1.0*Weights[nstep][nw];
    }
  }
  
  void FWSingle::resetWeights()
  {
    if (verbose>2) app_log()<<" Resetting Weights"<<endl;
    Weights.resize( IDs.size());
    for(int i=0;i<IDs.size();i++) Weights[i].resize(IDs[i].size(),0);
  }

  void FWSingle::FWOneStep()
  {
    //create an ordered version of the ParentIDs to make the weight calculation faster
    vector<vector<long> > orderedPIDs(PIDs);
    vector<vector<long> >::iterator it(orderedPIDs.begin());
    do
    {
      std::sort((*it).begin(),(*it).end());
      it++;
    } while(it<orderedPIDs.end());
    
    if (verbose>2) app_log()<<" Done Sorting IDs"<<endl;
    this->resetWeights();
    vector<vector<long> >::iterator stepIDIterator(IDs.begin());
    vector<vector<long> >::iterator stepPIDIterator(orderedPIDs.begin() + gensTransferred), last_stepPIDIterator(orderedPIDs.end());
    vector<vector<int> >::iterator stepWeightsIterator(Weights.begin());
    //we start comparing the next generations ParentIDs with the current generations IDs
    int i=0;
    do
    {
      if (verbose>2) app_log()<<"  calculating weights for gen:"<<gensTransferred<<" step:"<<i<<"/"<<orderedPIDs.size()<<endl;
      if (verbose>2) app_log()<<"Nsamples ="<<(*stepWeightsIterator).size()<<endl;
      
      vector<long>::iterator IDit( (*stepIDIterator).begin()     ), last_IDit( (*stepIDIterator).end());
      vector<long>::iterator PIDit( (*stepPIDIterator).begin()   ), last_PIDit( (*stepPIDIterator).end());
      vector<int>::iterator  Wit( (*stepWeightsIterator).begin() );
      if (verbose>2) app_log()<<"ID size:"<<(*stepIDIterator).size()<<" PID size:"<<(*stepPIDIterator).size()<<" Weight size:"<<(*stepWeightsIterator).size()<<endl;
      
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
        }
      }while(PIDit<last_PIDit);
      if (verbose>2) { printIDs((*stepIDIterator));printIDs((*stepPIDIterator));}
      if (verbose>2) printInts((*stepWeightsIterator));
      stepIDIterator++; stepPIDIterator++; stepWeightsIterator++; i++;
      
    }while(stepPIDIterator<orderedPIDs.end());
    
  }
  
  void FWSingle::printIDs(vector<long> vi)
  {
    for (int j=0;j<vi.size();j++) app_log()<<vi[j]<<" ";
    app_log()<<endl;
  }
    void FWSingle::printInts(vector<int> vi)
  {
    for (int j=0;j<vi.size();j++) app_log()<<vi[j]<<" ";
    app_log()<<endl;
  }
  
  void FWSingle::transferParentsOneGeneration( )
  {
//     set parentIDs such that they reflect the parents one more generations ago
//     must start from rear and work forward.
//     vector<vector<long> > orderedPIDs(PIDs);
//     vector<vector<long> >::iterator it(orderedPIDs.begin()),last_it(orderedPIDs.end());
//     do
//     {
//       std::sort((*it).begin(),(*it).end());
//       it++;
//     } while(it!=last_it);
//     vector<vector<long> >::iterator stepOrderedPIDIterator(orderedPIDs.begin());


    vector<vector<long> >::reverse_iterator stepIDIterator(IDs.rbegin()+gensTransferred);
    vector<vector<long> >::reverse_iterator stepPIDIterator(PIDs.rbegin()), nextStepPIDIterator(realPIDs.rbegin()+gensTransferred);
    int i(0);
    do
    {
      if (verbose>2) {printIDs((*nextStepPIDIterator));printIDs((*stepIDIterator));printIDs((*stepPIDIterator));}
      vector<long>::iterator hereID( (*stepIDIterator).begin() ) ;
      vector<long>::iterator nextStepPID( (*nextStepPIDIterator).begin() );
      vector<long>::iterator herePID( (*stepPIDIterator).begin() );
      if (verbose>2) app_log()<<"  calculating Parent IDs for gen:"<<gensTransferred<<" step:"<<i<<"/"<<PIDs.size()<<endl;
      do
      {
        if ((*herePID)==(*hereID)) 
        {
          (*herePID)=(*nextStepPID);
          herePID++;
        }
        else
        {
          hereID++; nextStepPID++;
          if (hereID==(*stepIDIterator).end()) 
          {
            hereID=(*stepIDIterator).begin();
            nextStepPID=(*nextStepPIDIterator).begin();
//             if (verbose>2) app_log()<<"resetting to beginning of parents"<<endl;
          }
        }
      }while(herePID<(*stepPIDIterator).end());
      stepIDIterator++;nextStepPIDIterator++;stepPIDIterator++;i++;
    }while(stepIDIterator<IDs.rend());
    gensTransferred++; //number of gens backward to compare parents
  }


  bool 
  FWSingle::put(xmlNodePtr q){

    fname<<xmlrootName<<".storeConfig.h5";
    c_file = H5Fopen(fname.str().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    hsize_t numGrps = 0;
    H5Gget_num_objs(c_file, &numGrps);
    numSteps = static_cast<int> (numGrps)-3;
    app_log()<<"Total number of steps in input file "<<numSteps<<endl;
    if (weightFreq<1) weightFreq=1;
    int numberDataPoints = weightLength/weightFreq;
    pointsToCalculate.resize(numberDataPoints);
    for(int i=0;i<numberDataPoints;i++) pointsToCalculate[i]=i*weightFreq;
    app_log()<<"  Observables will be calculated each "<<weightFreq<<" steps. At: ";
    for(int i=0;i<numberDataPoints;i++) app_log()<<pointsToCalculate[i]<<" ";
    app_log()<<endl;
    H5Fclose(c_file);
    
    return true;
  }
}

/***************************************************************************
 * $RCSfile: VMCParticleByParticle.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCParticleByParticle.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
