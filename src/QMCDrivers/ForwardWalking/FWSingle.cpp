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

namespace qmcplusplus { 

  /// Constructor.
  FWSingle::FWSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), weightFreq(1), weightLength(5), fname(""), verbose(0), WID("WalkerID"), PID("ParentID")
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

    
    calculateWeights();
    
    
    return finalize(0);
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
    if (verbose>1) printf("dataset rank %d, dimensions %lu x %lu\n", rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));
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
  
  void FWSingle::calculateWeights()
  {
    vector<long> lastRow;
    this->readInLong(numSteps-1 ,WID ,lastRow);
    int lastWalkerID=lastRow[lastRow.size()-1];
    if (verbose>1) app_log()<<"  Total number of walkers throughout entire run is "<<lastWalkerID<<endl;
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
