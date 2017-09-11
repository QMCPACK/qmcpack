
#include <fstream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <random>

#include "Configuration.h"
#include "AFQMC/config.h"
#include <Message/MPIObjectBase.h>
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/libxmldefs.h"
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Wavefunctions/WavefunctionBase.h"
#include "AFQMC/Wavefunctions/PureSingleDeterminant.h"
//#include "AFQMC/Walkers/SlaterDetWalker.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Propagators/phaseless_ImpSamp_ForceBias.h"
#include "AFQMC/Matrix/hdf5_readers.hpp"
#include "AFQMC/Matrix/array_partition.hpp"

#include"AFQMC/Numerics/SparseMatrixOperations.h"
#include"AFQMC/Numerics/DenseMatrixOperations.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/UtilityFunctions.h"
#include "AFQMC/Utilities/Utils.h"

namespace qmcplusplus
{


/*
 * Notes:
 * 1. The indexing on Spvn will be defined by walker_type:
 *      - walker_type=0 --> ik = [0,M2)
 *      - walker_type=1 --> ik = [0,2*M2)
 *      - walker_type=2 --> ik = [0,4*M2)
 *    Wavefunctions should adhere to this convention.
 * 2. One body operators in Wavefunction will now assume the same structure as Spvn for now.
 *    For off-diagonal observables in non-GHF case, define specific interfaces later.  
 */

bool phaseless_ImpSamp_ForceBias::parse(xmlNodePtr cur)
{

    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    std::string sub("yes");
    std::string use("yes");
    std::string constrain("yes");
    std::string save_mem("no");
    std::string hyb("no");
    std::string impsam("yes");
    ParameterSet m_param;
    m_param.add(sub,"substractMF","std::string");
    m_param.add(use,"useCholesky","std::string");
    m_param.add(cutoff,"cutoff_propg","double");
    m_param.add(cutoff,"cutoff","double");
    m_param.add(save_mem,"save_memory","std::string");
    m_param.add(vbias_bound,"vbias_bound","double");
    m_param.add(constrain,"apply_constrain","std::string");
    m_param.add(impsam,"importance_sampling","std::string");
    m_param.add(hyb,"hybrid_method","std::string");
    m_param.add(hyb,"hybrid","std::string");
    m_param.add(nnodes_per_TG,"nnodes_per_TG","int");
    m_param.add(nnodes_per_TG,"nnodes","int");
    m_param.add(nnodes_per_TG,"nodes","int");
    // hdf
    m_param.add(hdf_read_tag,"hdf_read_tag","std::string");
    m_param.add(hdf_read_file,"hdf_read_file","std::string");
    m_param.add(hdf_write_tag,"hdf_write_tag","std::string");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");
    // fix bias propagation
    m_param.add(fix_bias,"fix_bias","int");

    std::string par("yes");
    m_param.add(par,"paral_fac","std::string");
    m_param.add(par,"parallel_fac","std::string");
    m_param.add(par,"parallel_factorization","std::string");

    std::string par2("yes");
    m_param.add(par2,"parallel_propagation","std::string");

    std::string spr("no");
    m_param.add(spr,"dense","std::string");
    m_param.add(spr,"dense_propagator","std::string");

    m_param.add(walkerBlock,"walkerBlock","int");
    m_param.add(walkerBlock,"block","int");
    m_param.add(walkerBlock,"Block","int");

    m_param.add(n_reading_cores,"num_io_cores","int");

/* to be completed later

    m_param.add(useMFbias,"useMFbias","std::string");
    m_param.add(useHybrid,"useHybrid","std::string");

    // phaseless (require substractMF=no, apply_constrain=no, useMFbias=true, useHybrid=true)
    m_param.add(free_projection,"free_projection","std::string");

*/

    m_param.put(cur);

    substractMF=true;
    use_eig=false;
    apply_constrain=true;
    std::transform(sub.begin(),sub.end(),sub.begin(),(int (*)(int)) tolower);
    if(sub == "no" || sub == "no") substractMF = false;
    std::transform(use.begin(),use.end(),use.begin(),(int (*)(int)) tolower);
    if(use == "no" || use == "false") use_eig = true;
    std::transform(constrain.begin(),constrain.end(),constrain.begin(),(int (*)(int)) tolower);
    if(constrain == "no" || constrain == "false") apply_constrain = false;
    std::transform(save_mem.begin(),save_mem.end(),save_mem.begin(),(int (*)(int)) tolower);
    if(save_mem == "yes" || save_mem == "true") save_memory = true;  
    std::transform(par.begin(),par.end(),par.begin(),(int (*)(int)) tolower);
    if(par == "no" || par == "false") parallel_factorization = false;  
    std::transform(impsam.begin(),impsam.end(),impsam.begin(),(int (*)(int)) tolower);
    if(impsam == "false" || impsam == "no") imp_sampl = false;  
    std::transform(hyb.begin(),hyb.end(),hyb.begin(),(int (*)(int)) tolower);
    if(hyb == "yes" || hyb == "true") hybrid_method = true;  
    std::transform(par2.begin(),par2.end(),par2.begin(),(int (*)(int)) tolower);
    if(par2 == "no" || par2 == "false") parallelPropagation = false;  
    std::transform(spr.begin(),spr.end(),spr.begin(),(int (*)(int)) tolower);
    if(spr == "yes" || spr == "true") sparsePropagator = false;  

    app_log()<<"\n\n --------------- Parsing Propagator input ------------------ \n\n";

    if(substractMF) 
      app_log()<<" Using mean-field substraction in propagator. \n";

    if(parallel_factorization)
      app_log()<<" Calculating factorization of 2-body hamiltonian in parallel. \n";
    else
      app_log()<<" Calculating factorization of 2-body hamiltonian in serial. \n";

    if(parallelPropagation)
      app_log()<<" Using algorithm for parallel propagation (regardless of nnodes/ncores). \n";  
    else
      app_log()<<" Using algorithm for serial propagation. \n";  

    if(nnodes_per_TG > 1 && !parallel_factorization) {
      parallel_factorization=true;
      app_log()<<" Distributed propagator (nnodes_per_TG(propagator)>1) requires a parallel factorization. Setting parallel_factorization to true. \n"; 
      app_log()<<" Calculating factorization of 2-body hamiltonian in parallel. \n";
    }

    if(use_eig)
      app_log()<<" Calculating factorization of 2 body interaction with direct diagonalization.\n";
    else 
      app_log()<<" Calculating factorization of 2 body interaction with Cholesky method.\n"; 

    if(!imp_sampl) {
      app_log()<<" AFQMC propagation WITHOUT Importance Sampling! " <<"\n"; 
      app_log()<<" CAREFUL Incomplete implementation of propagation WITHOUT Importance Sampling! " <<"\n"; 
      app_log()<<" CAREFUL Weights are wrong and energies are incorrect !!! " <<"\n"; 
    }

    if(hybrid_method)
      app_log()<<" Using hybrid method to calculate the weights during the propagation." <<"\n";

    app_log()<<" Running Propagator with " <<nnodes_per_TG <<" nodes per task group. \n"; 

    if(fix_bias>0 && imp_sampl)
      app_log()<<" Using propagation scheme where bias potential is updated every " <<fix_bias <<" sub steps. \n";

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }


    app_log()<<std::endl;

    return true;

}

bool phaseless_ImpSamp_ForceBias::hdf_write_transposed(hdf_archive& dump,const std::string& tag)
{

  if(!sparsePropagator) {
    app_error()<<" WARNING: Propagator restart file can not be written yet in dense form. \n FILE NOT WRITTEN. \n";
    return true;    
  }

  // only heads of nodes in TG_number==0 communicate 
  int rk,npr;
  if(rank()==0) {

    if(TG.getTGNumber()!=0 || TG.getCoreRank()!=0 || TG.getTGRank()!=0)
      APP_ABORT("Error in phaseless_ImpSamp_ForceBias::hdf_write(): Global root is not head of TG=0.\n\n");


    std::string path = "/Propagators/phaseless_ImpSamp_ForceBias"; 
    if(tag != std::string("")) path += std::string("/")+tag; 
    if(dump.is_group( path )) {
      app_error()<<" ERROR: H5Group /Propagators/phaseless_ImpSamp_ForceBias/{tag} already exists in restart file. This is a bug and should not happen. Contact a developer.\n";
      return false;
    }

    // scale Spvn by 1/std::sqrt(dt) and write 
    RealType scale = 1.0/std::sqrt(dt); 
    Spvn *= scale; 

    // ranks of other heads in TG
    std::vector<int> ranks;
    int pos=0;
    if(nnodes_per_TG>1) TG.getRanksOfRoots(ranks,pos);
    else {
      ranks.push_back(0);
    }
    npr=ranks.size();

    std::vector<long> Ldata(4);
    long ntot=Spvn.size(), ntot_=Spvn.size();
    if(nnodes_per_TG>1) 
      MPI_Reduce(&ntot_, &ntot, 1, MPI_LONG, MPI_SUM, 0, TG.getTGCOMM()); 
    Ldata[0]=ntot;  // this needs to be summed over relevant nodes
    Ldata[1]=Spvn.rows();  // same here
    Ldata[2]=nCholVecs;
    Ldata[3]=NMO;

    dump.push("Propagators");
    dump.push("phaseless_ImpSamp_ForceBias");
    if(tag != std::string("")) dump.push(tag);
    dump.write(Ldata,"Spvn_dims");

    dump.write(vn0,"Spvn_vn0");


    // transpose matrix for easier access to Cholesky vectors 
    Spvn.transpose();   

    std::vector<int> Idata(4);
    Idata.resize(nCholVecs);
    for(int i=0; i<nCholVecs; i++)
      Idata[i] = *(Spvn.rowIndex_begin()+i+1) - *(Spvn.rowIndex_begin()+i);
    std::vector<int> sz_local;
    sz_local = Idata; 
    if(nnodes_per_TG>1) myComm->gsum(Idata,TG.getTGCOMM());

    dump.write(Idata,std::string("Spvn_nterms_per_vector"));

    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<SPValueType> vvec;
    ivec.reserve(nmax);
    vvec.reserve(nmax);

    MPI_Status st;  // do I need to delete these objects 
    std::vector<int> nvpp(1);
    nvpp[0]=nCholVecs;
    if(nnodes_per_TG > 1) {
      nvpp.resize(nnodes_per_TG);
      int cnt=0;
      for(int i=0; i<nnodes_per_TG; i++) {
        cnt+=nCholVec_per_node[i];
        nvpp[i]=cnt;
      }
    }

    int orig=0; 
    for(ValueSMSpMat::intType i=0; i<nCholVecs; i++) {

      if(nnodes_per_TG>1 && i == nvpp[orig]) orig++;
      if(orig==0) { 
        if(sz_local[i] == 0) 
          APP_ABORT("Error in phaseless_ImpSamp_ForceBias::hdf_write:: Problems with Spvn - sz_local[i]==0. \n\n");
        ivec.resize(Idata[i]);
        vvec.resize(Idata[i]);
        std::copy(Spvn.cols_begin()+(*(Spvn.rowIndex_begin()+i)),Spvn.cols_begin()+(*(Spvn.rowIndex_begin()+i+1)),ivec.begin());
        std::copy(Spvn.vals_begin()+(*(Spvn.rowIndex_begin()+i)),Spvn.vals_begin()+(*(Spvn.rowIndex_begin()+i+1)),vvec.begin());

      } else {
        if(sz_local[i] != 0) 
          APP_ABORT("Error in phaseless_ImpSamp_ForceBias::hdf_write:: Problems with Spvn - sz_local[i]!=0. \n\n");
        ivec.resize(Idata[i]);
        vvec.resize(Idata[i]);
        myComm->recv(ivec.data(),ivec.size(),ranks[orig],i,TG.getTGCOMM(),&st);
        myComm->recv(vvec.data(),vvec.size(),ranks[orig],i,TG.getTGCOMM(),&st);
      }

      dump.write(ivec,std::string("Spvn_index_")+std::to_string(i)); 
      dump.write(vvec,std::string("Spvn_vals_")+std::to_string(i)); 

    }

    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();

    Spvn.transpose();

    // scale back by std::sqrt(dt) 
    scale = std::sqrt(dt); 
    Spvn *= scale; 

    dump.flush();  

  } else if(nnodes_per_TG>1 && TG.getTGNumber()==0 && TG.getCoreRank()==0) {

    RealType scale = 1.0/std::sqrt(dt);
    Spvn *= scale;

    // ranks of other heads in TG
    std::vector<int> ranks;
    int pos=0;
    TG.getRanksOfRoots(ranks,pos);
    npr=ranks.size();

    long ntot=Spvn.size(), ntot_=Spvn.size();
    MPI_Reduce(&ntot_, &ntot, 1, MPI_LONG, MPI_SUM, 0, TG.getTGCOMM()); 

    Spvn.transpose();

    std::vector<int> Idata(nCholVecs);
    for(int i=0; i<nCholVecs; i++)
      Idata[i] = *(Spvn.rowIndex_begin()+i+1) - *(Spvn.rowIndex_begin()+i);
    std::vector<int> sz_local;
    sz_local = Idata;
    myComm->gsum(Idata,TG.getTGCOMM());

    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<SPValueType> vvec;
    ivec.reserve(nmax);
    vvec.reserve(nmax);

    for(ValueSMSpMat::intType i=cvec0; i<cvecN; i++) {

      ivec.resize(Idata[i]);
      vvec.resize(Idata[i]);
      std::copy(Spvn.cols_begin()+(*(Spvn.rowIndex_begin()+i)),Spvn.cols_begin()+(*(Spvn.rowIndex_begin()+i+1)),ivec.begin());
      std::copy(Spvn.vals_begin()+(*(Spvn.rowIndex_begin()+i)),Spvn.vals_begin()+(*(Spvn.rowIndex_begin()+i+1)),vvec.begin());

      myComm->send(ivec.data(),ivec.size(),0,i,TG.getTGCOMM());
      myComm->send(vvec.data(),vvec.size(),0,i,TG.getTGCOMM());
    }
    Spvn.transpose();

    scale = std::sqrt(dt);
    Spvn *= scale;

  } else if(nnodes_per_TG>1 && TG.getTGNumber()==0) {

    long ntot=0, ntot_=0;
    if(nnodes_per_TG>1)
      MPI_Reduce(&ntot_, &ntot, 1, MPI_LONG, MPI_SUM, 0, TG.getTGCOMM()); 

    std::vector<int> Idata;
    Idata.resize(nCholVecs);
    std::fill(Idata.begin(),Idata.end(),0);
    myComm->gsum(Idata,TG.getTGCOMM());

  } 

  return true;
 
}

bool phaseless_ImpSamp_ForceBias::hdf_write(hdf_archive& dump,const std::string& tag)
{

  if(!sparsePropagator) {
    app_error()<<" WARNING: Propagator restart file can not be written yet in dense form. \n FILE NOT WRITTEN. \n";
    return true;    
  }

  // only heads of nodes in TG_number==0 communicate 
  if(rank()==0) {

    if(TG.getTGNumber()!=0 || TG.getCoreRank()!=0 || TG.getTGRank()!=0)
      APP_ABORT("Error in phaseless_ImpSamp_ForceBias::hdf_write(): Global root is not head of TG=0.\n\n");

    std::string path = "/Propagators/phaseless_ImpSamp_ForceBias"; 
    if(tag != std::string("")) path += std::string("/")+tag; 
    if(dump.is_group( path )) {
      app_error()<<" ERROR: H5Group /Propagators/phaseless_ImpSamp_ForceBias/{tag} already exists in restart file. This is a bug and should not happen. Contact a developer.\n";
      return false;
    }

    dump.push("Propagators");
    dump.push("phaseless_ImpSamp_ForceBias");
    if(tag != std::string("")) dump.push(tag);
    dump.write(vn0,"Spvn_vn0");

    std::vector<SPComplexType> vvec;

    // write propg (not to be used in restart, only for miniapp)
    if(spinRestricted) 
      vvec.resize(NMO*NMO);
    else
      vvec.resize(2*NMO*NMO);

    for(int n=0; n<Propg_H1.size(); n++) { 
      int i, j;
      SPComplexType v;
      std::tie(i,j,v) = Propg_H1[n];
      vvec[i*NMO+j] = v;
    }
    dump.write(vvec,"Spvn_propg1");

    // scale Spvn by 1/std::sqrt(dt) and write 
    RealType scale = 1.0/std::sqrt(dt); 
    Spvn *= scale; 

    // write matrix
    long nblocks,ntot;
    std::tie(nblocks,ntot) = afqmc::write_hdf5_SpMat(Spvn,dump,std::string("Spvn"),intgs_per_block,TG);

    std::vector<long> Idata(5);
    Idata[0]=ntot;  
    Idata[1]=Spvn.rows(); 
    Idata[2]=nCholVecs;
    Idata[3]=NMO;
    Idata[4]=nblocks;
    dump.write(Idata,"Spvn_dims");

    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();

    // scale back by std::sqrt(dt) 
    scale = std::sqrt(dt); 
    Spvn *= scale; 

    dump.flush();  

  } else if(nnodes_per_TG>1 && TG.getTGNumber()==0 && TG.getCoreRank()==0) {

    RealType scale = 1.0/std::sqrt(dt);
    Spvn *= scale;

    int nblocks,ntot;
    std::tie(nblocks,ntot) = afqmc::write_hdf5_SpMat(Spvn,dump,std::string("Spvn"),intgs_per_block,TG);

    scale = std::sqrt(dt);
    Spvn *= scale;

  } else if(nnodes_per_TG>1 && TG.getTGNumber()==0) {

    int nblocks,ntot;
    std::tie(nblocks,ntot) = afqmc::write_hdf5_SpMat(Spvn,dump,std::string("Spvn"),intgs_per_block,TG);

  } 
  // just for safety
  myComm->barrier();

  return true;
 
}

bool phaseless_ImpSamp_ForceBias::hdf_read(hdf_archive& dump,const std::string& tag)
{
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID() , coreid = TG.getCoreID();
  int nread = (n_reading_cores<=0)?(TG.getTotalCores()):(n_reading_cores);
  std::vector<long> dims(5);
  int rk,npr; 
  if(n_reading_cores <=0 || coreid < n_reading_cores) {
    // these cores access hdf file
    if(!dump.push("Propagators",false)) return false;
    if(!dump.push("phaseless_ImpSamp_ForceBias",false)) return false;
    if(tag != std::string("")) 
      if(!dump.push(tag,false)) return false;
  }

  // only 1 core reads small datasets
  if( rank() == 0 ) {
    if(!dump.read(dims,"Spvn_dims")) return false;

    if(int(dims[3]) != NMO) {
      app_error()<<" Error in phaseless_ImpSamp_ForceBias::hdf_read. NMO is not consistent between hdf5 file and current run. \n";
      return false; 
    }  

    // read basic info
    MPI_Bcast(dims.data(), dims.size(), MPI_LONG, 0, myComm->getMPI());
    long ntot = dims[0];
    int nrows = int(dims[1]);  
    nCholVecs = int(dims[2]); 
    int nblk = int(dims[4]);

    // read vn0
    if(!dump.read(vn0,"Spvn_vn0")) return false;
    myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,myComm->getMPI());

    afqmc::simple_matrix_partition<afqmc::TaskGroup,IndexType,RealType> split(nrows,nCholVecs,cutoff);
    std::vector<IndexType> counts;
    // count dimensions of sparse matrix
    afqmc::count_entries_hdf5_SpMat(dump,split,std::string("Spvn"),nblk,false,counts,TG,true,nread);

    std::vector<IndexType> sets;
    split.partition(TG,false,counts,sets);

    app_log()<<" Partitioning of Cholesky Vectors: ";
    for(int i=0; i<=nnodes_per_TG; i++)
      app_log()<<sets[i] <<" ";
    app_log()<<std::endl;
    app_log()<<" Number of terms in each partitioning: ";
    for(int i=0; i<nnodes_per_TG; i++)
      app_log()<<accumulate(counts.begin()+sets[i],counts.begin()+sets[i+1],0) <<" ";
    app_log()<<std::endl;

    MPI_Bcast(sets.data(), sets.size(), MPI_INT, 0, myComm->getMPI());
    cvec0 = sets[TG.getLocalNodeNumber()];
    cvecN = sets[TG.getLocalNodeNumber()+1];
    for(int i=0; i<nnodes_per_TG; i++)
      nCholVec_per_node[i] = sets[i+1]-sets[i];
    GlobalSpvnSize = std::accumulate(counts.begin(),counts.end(),0);

    // resize Spvn
    int sz = std::accumulate(counts.begin()+cvec0,counts.begin()+cvecN,0);
    MPI_Bcast(&sz, 1, MPI_INT, 0, TG.getNodeCommLocal());

    Spvn.setDims(nrows,nCholVecs);
    Spvn.allocate(sz);

    // read Spvn
    afqmc::read_hdf5_SpMat(Spvn,split,dump,std::string("Spvn"),nblk,TG,true,nread);

    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();


  } else { 
    // these cores only read sparse matrix 
    
    MPI_Bcast(dims.data(), dims.size(), MPI_LONG, 0, myComm->getMPI());
    int nrows = int(dims[1]);
    nCholVecs = int(dims[2]);
    int nblk = int(dims[4]);

    myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,myComm->getMPI());

    std::vector<IndexType> counts;
    afqmc::simple_matrix_partition<afqmc::TaskGroup,IndexType,RealType> split(nrows,nCholVecs,cutoff);
    // count dimensions of sparse matrix
    afqmc::count_entries_hdf5_SpMat(dump,split,std::string("Spvn"),nblk,false,counts,TG,true,nread);

    std::vector<IndexType> sets(nnodes_per_TG+1);
    MPI_Bcast(sets.data(), sets.size(), MPI_INT, 0, myComm->getMPI());
    cvec0 = sets[TG.getLocalNodeNumber()];
    cvecN = sets[TG.getLocalNodeNumber()+1];
    for(int i=0; i<nnodes_per_TG; i++)
      nCholVec_per_node[i] = sets[i+1]-sets[i];

    // resize Spvn
    int sz;
    if( coreid==0 ) 
      sz = std::accumulate(counts.begin()+cvec0,counts.begin()+cvecN,0);
    MPI_Bcast(&sz, 1, MPI_INT, 0, TG.getNodeCommLocal());

    Spvn.setDims(nrows,nCholVecs);
    Spvn.allocate(sz);

    if(n_reading_cores <=0 || coreid < n_reading_cores) {
      split.partition(TG,false,counts,sets);

      // read Spvn
      afqmc::read_hdf5_SpMat(Spvn,split,dump,std::string("Spvn"),nblk,TG,true,nread);

      if(tag != std::string("")) dump.pop();
      dump.pop();
      dump.pop();
    } else {
      // read Spvn
      afqmc::read_hdf5_SpMat(Spvn,split,dump,std::string("Spvn"),nblk,TG,true,nread);
    }
  } 

  // compression is done in read_hdf5_SpMat
  //Spvn.compress(TG.getNodeCommLocal());

  if(coreid == 0) {
    // scale back by std::sqrt(dt) 
    RealType scale = std::sqrt(dt); 
    Spvn *= scale; 
  }

  myComm->barrier();

  return true;
}

// do not remove, in case you find machines with bad I/O bottlenecks, e.g. BGQ
bool phaseless_ImpSamp_ForceBias::hdf_read_transposed(hdf_archive& dump,const std::string& tag)
{

  std::vector<int> Idata(4);
  int rk,npr; 
  // in this case, head_of_nodes communicate.
  // the partitioning between nodes in TG is done here
  // all vectors are sent to all heads and the decision to keep the data is made locally 
  if(rank() == 0) {

    if(!dump.push("Propagators",false)) return false;
    if(!dump.push("phaseless_ImpSamp_ForceBias",false)) return false;
    if(tag != std::string("")) 
      if(!dump.push(tag,false)) return false;
    if(!dump.read(Idata,"Spvn_dims")) return false;

    if(Idata[3] != NMO) {
      app_error()<<" Error in phaseless_ImpSamp_ForceBias::hdf_read. NMO is not consistent between hdf5 file and current run. \n";
      return false; 
    }  

    MPI_Comm_rank(MPI_COMM_HEAD_OF_NODES,&rk);
    MPI_Comm_size(MPI_COMM_HEAD_OF_NODES,&npr);

    myComm->bcast(Idata);

    int ntot = Idata[0];
    int nrows = Idata[1];  
    nCholVecs = Idata[2]; 

    if(!dump.read(vn0,"Spvn_vn0")) return false;
    myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,myComm->getMPI());

    Idata.resize(nCholVecs);
    if(!dump.read(Idata,"Spvn_nterms_per_vector")) return false;
    myComm->bcast(Idata,MPI_COMM_HEAD_OF_NODES);    

    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<ValueType> vvec;
    ivec.reserve(nmax);
    vvec.reserve(nmax);     

    // distribute cholesky vectors among nodes in TG
    if(npr > nCholVecs)
      APP_ABORT("Error in hdf_read(): nnodes_per_TG > nCholVecs. \n\n\n");

    std::vector<int> nv;
    int node_number = TG.getLocalNodeNumber();  
    if(nnodes_per_TG > 1) {
      nv.resize(nCholVecs+1);
      nv[0]=0;
      for(int i=0,cnt=0; i<nCholVecs; i++) {
        cnt+=Idata[i];
        nv[i+1]=cnt;
      } 
      std::vector<int> sets(nnodes_per_TG+1);
      balance_partition_ordered_set(nCholVecs,nv.data(),sets);
      cvec0 = sets[node_number];
      cvecN = sets[node_number+1];
      app_log()<<" Partitioning of Cholesky Vectors: ";
      for(int i=0; i<=nnodes_per_TG; i++)
        app_log()<<sets[i] <<" ";
      app_log()<<std::endl;
      app_log()<<" Number of terms in each partitioning: "; 
      for(int i=0; i<nnodes_per_TG; i++)
        app_log()<<accumulate(Idata.begin()+sets[i],Idata.begin()+sets[i+1],0) <<" "; 
      app_log()<<std::endl;
      myComm->bcast(sets,MPI_COMM_HEAD_OF_NODES);
      for(int i=0; i<=nnodes_per_TG; i++)
        nCholVec_per_node[i] = sets[i+1]-sets[i]; 
      myComm->bcast(nCholVec_per_node);
      GlobalSpvnSize = std::accumulate(Idata.begin(),Idata.end(),0);
    } else {
      cvec0 = 0;
      cvecN = nCholVecs; 
      nCholVec_per_node[0]=nCholVecs;
      GlobalSpvnSize = std::accumulate(Idata.begin(),Idata.end(),0);
    }

    int sz = std::accumulate(Idata.begin()+cvec0,Idata.begin()+cvecN,0);
    myComm->bcast(&sz,1,TG.getNodeCommLocal());

    Spvn.setDims(nrows,nCholVecs);
    Spvn.resize(sz);

    Spvn.share(&cvec0,1,true);
    Spvn.share(&cvecN,1,true);

    SPValueSMSpMat::iterator itv = Spvn.vals_begin();
    ValueSMSpMat::int_iterator itc = Spvn.cols_begin();
    ValueSMSpMat::int_iterator itr = Spvn.rows_begin();
/*
#if MPI_VERSION >= 3
#define HAVE_MPI_IALLGATHER 
#define HAVE_MPI_IGATHER 
#define HAVE_MPI_IBCAST 
#endif
#ifdef HAVE_MPI_IBCAST
    // should I start assuming c++14??? Much easier there! 
    typedef std::tuple< MPI_Request*, std::vector<IndexType>*, MPI_Request*, std::vector<ValueType>* >  commData; 

    auto nextit = [] (const std::vector<commData>::iterator& it, const std::vector<commData>::iterator& itbegin, const std::vector<commData>::iterator& itend) { return ((it+1)==itend)?itbegin:(it+1); };

    std::vector<IndexType> ivec2;
    std::vector<ValueType> vvec2;
    ivec2.reserve(ivec.size());
    vvec2.reserve(vvec.size());

    MPI_Request req1, req2, req3, req4;
    MPI_Status st1, st2;
    std::vector<commData> data(2);
    data[0] = std::make_tuple( &req1, &ivec, &req2, &vvec );    
    data[1] = std::make_tuple( &req3, &ivec2, &req4, &vvec2 );    
    std::vector<commData>::iterator it = data.begin(), itend = data.end(), itbeg = data.begin();

    std::get<1>(*it)->resize(Idata[0]);
    std::get<3>(*it)->resize(Idata[0]);
    if(!dump.read(*(std::get<1>(*it)),std::string("Spvn_index_")+std::to_string(0))) return false;
    if(!dump.read(*(std::get<3>(*it)),std::string("Spvn_vals_")+std::to_string(0))) return false;       
    MPI_Ibcast( std::get<1>(*it)->data() , std::get<1>(*it)->size(),MPI_INT,0,MPI_COMM_HEAD_OF_NODES,std::get<0>(*it));
#if defined(QMC_COMPLEX)
    MPI_Ibcast( std::get<3>(*it)->data() , 2*std::get<3>(*it)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it));
#else
    MPI_Ibcast( std::get<3>(*it)->data() , std::get<3>(*it)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it));
#endif
    for(ValueSMSpMat::intType i=0; i<nCholVecs; i++) {

      std::vector<commData>::iterator it2 = nextit(it,itbeg,itend);
      // initiate communication
      if(i+1 < nCholVecs) {
        std::get<1>(*it2)->resize(Idata[i+1]);
        std::get<1>(*it2)->resize(Idata[i+1]);
        if(!dump.read(*(std::get<1>(*it2)),std::string("Spvn_index_")+std::to_string(i+1))) return false;
        if(!dump.read(*(std::get<3>(*it2)),std::string("Spvn_vals_")+std::to_string(i+1))) return false;       
        MPI_Ibcast( std::get<1>(*it2)->data() , std::get<1>(*it2)->size(),MPI_INT,0,MPI_COMM_HEAD_OF_NODES,std::get<0>(*it2));
#if defined(QMC_COMPLEX)
        MPI_Ibcast( std::get<3>(*it2)->data() , 2*std::get<3>(*it2)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it2));
#else
        MPI_Ibcast( std::get<3>(*it2)->data() , std::get<3>(*it2)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it2));
#endif
      }
      
      // wait for current message
      MPI_Wait(std::get<0>(*it),&st1);
      MPI_Wait(std::get<2>(*it),&st2);

      if(nnodes_per_TG == 1 || (i>=cvec0 && i < cvecN)) {
        for(int k=0; k<std::get<1>(*it)->size(); k++) {
          *(itv++) = static_cast<SPValueType>((*std::get<3>(*it))[k]);
          *(itr++) = (*std::get<1>(*it))[k];
          *(itc++) = i;
        }
      }

      // push iterator forward
      it = it2;

    }
#else
*/
    for(ValueSMSpMat::intType i=0; i<nCholVecs; i++) {  

      ivec.resize(Idata[i]); 
      vvec.resize(Idata[i]); 
      if(!dump.read(ivec,std::string("Spvn_index_")+std::to_string(i))) return false;
      if(!dump.read(vvec,std::string("Spvn_vals_")+std::to_string(i))) return false;           
      myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
#if defined(QMC_COMPLEX)
      myComm->bcast<RealType>(reinterpret_cast<RealType*>(vvec.data()),2*vvec.size(),MPI_COMM_HEAD_OF_NODES);
#else
      myComm->bcast<ValueType>(vvec.data(),vvec.size(),MPI_COMM_HEAD_OF_NODES);
#endif

      if(nnodes_per_TG == 1 || (i>=cvec0 && i < cvecN)) {
        for(int k=0; k<ivec.size(); k++) {
          *(itv++) = static_cast<SPValueType>(vvec[k]); 
          *(itr++) = ivec[k]; 
          *(itc++) = i; 
        }
      }

    }
//#endif

    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();

    Spvn.compress(TG.getNodeCommLocal());

    // scale back by std::sqrt(dt) 
    RealType scale = std::sqrt(dt); 
    Spvn *= scale; 
  
  } else if(head_of_nodes) {

    MPI_Comm_rank(MPI_COMM_HEAD_OF_NODES,&rk);
    MPI_Comm_size(MPI_COMM_HEAD_OF_NODES,&npr);

    myComm->bcast(Idata);
    int ntot = Idata[0];
    int nrows = Idata[1];
    nCholVecs = Idata[2];

    myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,myComm->getMPI());

    Idata.resize(nCholVecs);
    myComm->bcast(Idata,MPI_COMM_HEAD_OF_NODES);
    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<ValueType> vvec;
    ivec.reserve(nmax);
    vvec.reserve(nmax);     

    cvec0 = 0;
    cvecN = nCholVecs; 
    int node_number = TG.getLocalNodeNumber();
    if(nnodes_per_TG > 1) {
      std::vector<int> sets(nnodes_per_TG+1);
      myComm->bcast(sets,MPI_COMM_HEAD_OF_NODES);
      cvec0 = sets[node_number];
      cvecN = sets[node_number+1];
      myComm->bcast(nCholVec_per_node);
    } else {
      nCholVec_per_node[0]=nCholVecs;
    }
    int sz = std::accumulate(Idata.begin()+cvec0,Idata.begin()+cvecN,0);
    myComm->bcast(&sz,1,TG.getNodeCommLocal());   
    
    Spvn.setDims(nrows,nCholVecs);
    Spvn.resize(sz);

    Spvn.share(&cvec0,1,true);
    Spvn.share(&cvecN,1,true);

    SPValueSMSpMat::iterator itv = Spvn.vals_begin();
    ValueSMSpMat::int_iterator itc = Spvn.cols_begin();
    ValueSMSpMat::int_iterator itr = Spvn.rows_begin();
/*
#if MPI_VERSION >= 3
#define HAVE_MPI_IALLGATHER 
#define HAVE_MPI_IGATHER 
#define HAVE_MPI_IBCAST 
#endif
#ifdef HAVE_MPI_IBCAST
    // should I start assuming c++14??? Much easier there! 
    typedef std::tuple< MPI_Request*, std::vector<IndexType>*, MPI_Request*, std::vector<ValueType>* >  commData; 

    auto nextit = [] (const std::vector<commData>::iterator& it, const std::vector<commData>::iterator& itbegin, const std::vector<commData>::iterator& itend) { return ((it+1)==itend)?itbegin:(it+1); };

    std::vector<IndexType> ivec2;
    std::vector<ValueType> vvec2;
    ivec2.reserve(ivec.size());
    vvec2.reserve(vvec.size());

    MPI_Request req1, req2, req3, req4;
    MPI_Status st1, st2;
    std::vector<commData> data(2);
    data[0] = std::make_tuple( &req1, &ivec, &req2, &vvec );    
    data[1] = std::make_tuple( &req3, &ivec2, &req4, &vvec2 );    
    std::vector<commData>::iterator it = data.begin(), itend = data.end(), itbeg = data.begin();

    std::get<1>(*it)->resize(Idata[0]);
    std::get<3>(*it)->resize(Idata[0]);
    if(!dump.read(*(std::get<1>(*it)),std::string("Spvn_index_")+std::to_string(0))) return false;
    if(!dump.read(*(std::get<3>(*it)),std::string("Spvn_vals_")+std::to_string(0))) return false;       
    MPI_Ibcast( std::get<1>(*it)->data() , std::get<1>(*it)->size(),MPI_INT,0,MPI_COMM_HEAD_OF_NODES,std::get<0>(*it));
#if defined(QMC_COMPLEX)
    MPI_Ibcast( std::get<3>(*it)->data() , 2*std::get<3>(*it)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it));
#else
    MPI_Ibcast( std::get<3>(*it)->data() , std::get<3>(*it)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it));
#endif
    for(ValueSMSpMat::intType i=0; i<nCholVecs; i++) {

      std::vector<commData>::iterator it2 = nextit(it,itbeg,itend);
      // initiate communication
      if(i+1 < nCholVecs) {
        std::get<1>(*it2)->resize(Idata[i+1]);
        std::get<1>(*it2)->resize(Idata[i+1]);
        MPI_Ibcast( std::get<1>(*it2)->data() , std::get<1>(*it2)->size(),MPI_INT,0,MPI_COMM_HEAD_OF_NODES,std::get<0>(*it2));
#if defined(QMC_COMPLEX)
        MPI_Ibcast( std::get<3>(*it2)->data() , 2*std::get<3>(*it2)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it2));
#else
        MPI_Ibcast( std::get<3>(*it2)->data() , std::get<3>(*it2)->size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES,std::get<2>(*it2));
#endif
      }
      
      // wait for current message
      MPI_Wait(std::get<0>(*it),&st1);
      MPI_Wait(std::get<2>(*it),&st2);

      if(nnodes_per_TG == 1 || (i>=cvec0 && i < cvecN)) {
        for(int k=0; k<std::get<1>(*it)->size(); k++) {
          *(itv++) = static_cast<SPValueType>((*std::get<3>(*it))[k]);
          *(itr++) = (*std::get<1>(*it))[k];
          *(itc++) = i;
        }
      }

      // push iterator forward
      it = it2;

    }
#else
*/
    for(ValueSMSpMat::intType i=0; i<nCholVecs; i++) {

      ivec.resize(Idata[i]);
      vvec.resize(Idata[i]);

      myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
#if defined(QMC_COMPLEX)
      myComm->bcast<RealType>(reinterpret_cast<RealType*>(vvec.data()),2*vvec.size(),MPI_COMM_HEAD_OF_NODES);
#else
      myComm->bcast<ValueType>(vvec.data(),vvec.size(),MPI_COMM_HEAD_OF_NODES);
#endif

      if(nnodes_per_TG == 1 || (i>=cvec0 && i < cvecN)) {
        for(int k=0; k<ivec.size(); k++) {
          *(itv++) = static_cast<SPValueType>(vvec[k]); 
          *(itr++) = ivec[k];
          *(itc++) = i;      
        }
      }

    }
//#endif

    Spvn.compress(TG.getNodeCommLocal());

    // scale back by std::sqrt(dt) 
    RealType scale = std::sqrt(dt); 
    Spvn *= scale; 

  } else {

    myComm->bcast(Idata);
    int ntot = Idata[0];
    int nrows = Idata[1];
    nCholVecs = Idata[2];

    myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,myComm->getMPI());

    if(nnodes_per_TG > 1) 
      myComm->bcast(nCholVec_per_node);
    else 
      nCholVec_per_node[0]=nCholVecs;

    int sz;
    myComm->bcast(&sz,1,TG.getNodeCommLocal());   

    Spvn.setDims(nrows,nCholVecs);
    Spvn.resize(sz);

    Spvn.share(&cvec0,1,false);
    Spvn.share(&cvecN,1,false);

    Spvn.compress(TG.getNodeCommLocal());

  } 

  myComm->barrier();

  return true;
}



bool phaseless_ImpSamp_ForceBias::setup(std::vector<int>& TGdata, SPComplexSMVector* v, HamiltonianBase* ham,WavefunctionHandler* w, RealType dt_, hdf_archive& dump_read, const std::string& hdf_restart_tag, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm node_heads_comm)
{

  dt = dt_;

  wfn=w;
  sizeOfG = 0;
  vHS_size = NMO*NMO;
  if(!spinRestricted) vHS_size+=NMO*NMO;

  ncores_per_TG=TGdata[4];
  if(!TG.quick_setup(ncores_per_TG,nnodes_per_TG,TGdata[0],TGdata[1],TGdata[2],TGdata[3]))
    return false;
  TG.setBuffer(v);
  core_rank = TG.getCoreRank();
  TG.setNodeCommLocal(node_comm);
  TG.setTGCommLocal(tg_comm);
  MPI_COMM_HEAD_OF_NODES = node_heads_comm;
  TG.setHeadOfNodesComm(MPI_COMM_HEAD_OF_NODES);
  head_of_nodes = (TG.getCoreID()==0);

  walker_per_node.resize(nnodes_per_TG);

  parallelPropagation = (parallelPropagation||(nnodes_per_TG>1 || ncores_per_TG>1));
  distributeSpvn = nnodes_per_TG>1;

  // this is the number of ComplexType needed to store 
  // all the information necessary to calculate 
  // the local contribution to vHS. This includes Green functions, overlaps, etc.
  if(parallelPropagation)
    sizeOfG = wfn->sizeOfInfoForDistributedPropagation("ImportanceSampling");

  bool read_Spvn_from_file=false;
  Propg_H1_indx.resize(4); 
  ComplexMatrix Hadd(2*NMO,NMO);
  for(int i=0; i<2*NMO; i++) 
   for(int j=0; j<NMO; j++) 
    Hadd(i,j)=0;
  
  if(sparsePropagator)
    Spvn.setup(head_of_nodes,"Spvn",TG.getNodeCommLocal());
  else
    Dvn.setup(head_of_nodes,"Dvn",TG.getNodeCommLocal());
  if(spinRestricted)  
    vn0.resize(NMO,NMO);
  else
    vn0.resize(2*NMO,NMO);
    
//  std::string str_ = std::string("debug.") + std::to_string(rank());
//  out_debug.open(str_.c_str());

  nCholVec_per_node.resize(nnodes_per_TG);

  if(hdf_read_file!=std::string("")) {

    if(!sparsePropagator)
      APP_ABORT("Error: Propagator restart not implemented with dense propagator yet. \n");

    hdf_archive dump(myComm);
    int coreid = TG.getCoreID();
    if( n_reading_cores <=0 || coreid < n_reading_cores ) {  // only those that will read
      if(dump.open(hdf_read_file,H5F_ACC_RDONLY)) { 
        read_Spvn_from_file = hdf_read(dump,hdf_read_tag);
        dump.close();
        if(read_Spvn_from_file) 
          app_log()<<"Successfully read HS potentials from file: " <<hdf_read_file <<"\n";
      } else {
        APP_ABORT("Error in phaseless_ImpSamp_ForceBias::setup(): problems opening hdf5 file. \n");
      }
    } else
      read_Spvn_from_file = hdf_read(dump,hdf_read_tag);

  } 

  if(!read_Spvn_from_file && !parallel_factorization) {

    if(nnodes_per_TG > 1) { 
      app_error()<<" Error in phaseless_ImpSamp_ForceBias::setup(): nnodes_per_TG > 1 only implemented with parallel_factorization=true \n" <<std::endl;
      return false;
    }  

    if(rank()==0) {

      app_log()<<" Calculating HS potentials from scratch. \n";

      Timer.reset("Generic1");
      Timer.start("Generic1");

      // calculates Hubbard-Stratonovich potentials (vn) 
      if(use_eig) {
        if(nnodes_per_TG > 1) { 
          app_error()<<" Error: nnodes_per_TG > 1 not implemented with use_eig=true \n" <<std::endl;
          return false;
        }  
        ham->calculateHSPotentials_Diagonalization(cutoff,dt,vn0,Spvn,Dvn,TG,nCholVec_per_node,sparsePropagator,false);
      } else
        ham->calculateHSPotentials(cutoff, dt, vn0,Spvn,Dvn,TG,nCholVec_per_node,sparsePropagator,false);

      Timer.stop("Generic1");
      app_log()<<" -- Time to calculate HS potentials: " <<Timer.average("Generic1") <<"\n";

      myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,MPI_COMM_HEAD_OF_NODES);
      if(sparsePropagator) {
        std::vector<int> ni(3);
        ni[0]=Spvn.size();
        ni[1]=Spvn.rows();
        ni[2]=Spvn.cols();
        myComm->bcast(ni);        

        nCholVecs = Spvn.cols();
        cvec0 = 0;
        cvecN = nCholVecs;

        myComm->bcast<char>(reinterpret_cast<char*>(Spvn.values()),sizeof(SPValueType)*Spvn.size(),0,MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.row_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.column_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.row_index(),Spvn.rows()+1,MPI_COMM_HEAD_OF_NODES);

        myComm->barrier();
        myComm->barrier();

        GlobalSpvnSize = Spvn.size();
      } else {
        std::vector<int> ni(3);
        ni[0]=Dvn.size();
        ni[1]=((spinRestricted)?(NMO*NMO):(2*NMO*NMO));
        ni[2]=nCholVec_per_node[0];
        myComm->bcast(ni);

        nCholVecs = nCholVec_per_node[0];
        cvec0 = 0;
        cvecN = nCholVecs;

        myComm->bcast<char>(reinterpret_cast<char*>(Dvn.values()),sizeof(SPValueType)*Dvn.size(),0,MPI_COMM_HEAD_OF_NODES);

        myComm->barrier();
        myComm->barrier();

        GlobalSpvnSize = Dvn.size(); 
      } 

    } else {

      std::vector<int> ni(3);
      myComm->bcast(ni);

      if(sparsePropagator) {
        Spvn.setDims(ni[1],ni[2]);
    
        nCholVecs = ni[2];
        cvec0 = 0;
        cvecN = nCholVecs;

        if(head_of_nodes) {
          Spvn.allocate_serial(ni[0]);
          Spvn.resize_serial(ni[0]);

          myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<char>(reinterpret_cast<char*>(Spvn.values()),sizeof(SPValueType)*Spvn.size(),0,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<int>(Spvn.row_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<int>(Spvn.column_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<int>(Spvn.row_index(),Spvn.rows()+1,MPI_COMM_HEAD_OF_NODES);
        }
        Spvn.setCompressed();

        myComm->barrier();
        if(!head_of_nodes) Spvn.initializeChildren();
        myComm->barrier();

      } else {  // !sparsePropagator

        nCholVecs = ni[2];
        cvec0 = 0;
        cvecN = nCholVecs;

        if(head_of_nodes) {
          Dvn.allocate_serial(ni[0]);
          Dvn.resize_serial(ni[0]);

          myComm->bcast<char>(reinterpret_cast<char*>(vn0.data()),sizeof(ValueType)*vn0.size(),0,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<char>(reinterpret_cast<char*>(Dvn.values()),sizeof(SPValueType)*Dvn.size(),0,MPI_COMM_HEAD_OF_NODES);
        }

        myComm->barrier();
        if(!head_of_nodes) Dvn.initializeChildren();
        myComm->barrier();

      }
    }

  } else if(!read_Spvn_from_file) {
    // calculating Spvn in parallel

    app_log()<<" Calculating HS potentials from scratch. \n";

    Timer.reset("Generic1");
    Timer.start("Generic1");

    // calculates Hubbard-Stratonovich potentials (vn)
    if(use_eig) {
      if(nnodes_per_TG > 1) { 
        app_error()<<" Error: nnodes_per_TG > 1 not implemented with use_eig=true \n" <<std::endl;
        return false;
      }  
      ham->calculateHSPotentials_Diagonalization(cutoff,dt,vn0,Spvn,Dvn,TG,nCholVec_per_node,sparsePropagator,true);
    } else
      ham->calculateHSPotentials(cutoff, dt, vn0,Spvn,Dvn,TG,nCholVec_per_node,sparsePropagator,true);

    Timer.stop("Generic1");
    app_log()<<" -- Time to calculate HS potentials: " <<Timer.average("Generic1") <<"\n";

    nCholVecs = std::accumulate(nCholVec_per_node.begin(),nCholVec_per_node.end(),0); 
    cvec0 = 0;
    cvecN = nCholVecs;
    if(nnodes_per_TG>1) {
      if(sparsePropagator)
        GlobalSpvnSize=Spvn.size();
      int sz_ = (core_rank==0)?GlobalSpvnSize:0;
      MPI_Allreduce(&sz_,&GlobalSpvnSize,1,MPI_INT,MPI_SUM,TG.getTGCOMM());  
      int cnt=0,node_number = TG.getLocalNodeNumber();
      if(node_number==0)
        cvec0=0;
      else
        cvec0 = std::accumulate(nCholVec_per_node.begin(),nCholVec_per_node.begin()+node_number,0);
      cvecN = cvec0+nCholVec_per_node[node_number];
    } else {
      // collect this from head processors in TG=0
      if(sparsePropagator)
        GlobalSpvnSize=Spvn.size();
    }
  } 

  if(!sparsePropagator)
    Dvn.setDims( (spinRestricted?NMO*NMO:2*NMO*NMO) , cvecN-cvec0 );   

  if(sparsePropagator) {
    app_log()<<"Memory used by HS potential: " <<(Spvn.capacity()*sizeof(SPValueType)+(2*Spvn.capacity()+NMO*NMO+1)*sizeof(int)+8000)/1024.0/1024.0 <<" MB " <<std::endl;
    RealType vnmax1=0,vnmax2=0;
    for(int i=0; i<Spvn.size(); i++) {
      if( std::abs( *(Spvn.values()+i) ) > vnmax1 )
        vnmax1 = std::abs(*(Spvn.values()+i)); 
#if defined(QMC_COMPLEX)
      if( (Spvn.values()+i)->imag()  > vnmax2 )
        vnmax2 = (Spvn.values()+i)->imag(); 
#endif
    }
    app_log()<<" Largest term in Vn: " <<vnmax1 <<" " <<vnmax2 <<std::endl; 
  } else {
    app_log()<<"Memory used by HS potential: " <<(Dvn.size()*sizeof(SPValueType)+8000)/1024.0/1024.0 <<" MB " <<std::endl;
    RealType vnmax1=0,vnmax2=0;
    for(int i=0; i<Dvn.size(); i++) {
      if( std::abs( *(Dvn.values()+i) ) > vnmax1 )
        vnmax1 = std::abs(*(Dvn.values()+i));
#if defined(QMC_COMPLEX)
      if( (Dvn.values()+i)->imag()  > vnmax2 )
        vnmax2 = (Dvn.values()+i)->imag();
#endif
    }
    app_log()<<" Largest term in Vn: " <<vnmax1 <<" " <<vnmax2 <<std::endl;
  }

  hybrid_weight.reserve(1000);
  MFfactor.reserve(1000);
  vMF.resize(nCholVecs);  // keeping storage for all vectors for now, changing this takes modifications to Spvn
  SPvHS.resize(2*NMO,NMO);
  for(int i=0; i<vMF.size(); i++) vMF[i]=0;

  // vMF is actually sqrt(dt)*<vn>, since Spvn is also scaled by sqrt(dt)
  if(substractMF) { 
    if(sparsePropagator) 
      wfn->calculateMeanFieldMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,Spvn,vMF);
    else
      wfn->calculateMeanFieldMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,Dvn,vMF);
    if(nnodes_per_TG>1) {
      if(rank()==0) {
        MPI_Status st;
        std::vector<ComplexType> v_(vMF.size());
        std::vector<int> ranks;
        int pos=0;
        TG.getRanksOfRoots(ranks,pos);
        for(int i=1; i<nnodes_per_TG; i++) {
          myComm->recv(v_.data(),v_.size(),ranks[i],i,TG.getTGCOMM(),&st);
          for(int k=0; k<vMF.size(); k++) vMF[k] += v_[k];
        }
      } else if(TG.getTGNumber()==0 && TG.getCoreRank()==0) {
        std::vector<int> ranks;
        int pos=0;
        TG.getRanksOfRoots(ranks,pos);        
        myComm->send(vMF.data(),vMF.size(),0,pos,TG.getTGCOMM());  
      }
      myComm->bcast(vMF.data(),vMF.size(),0,MPI_COMM_WORLD);
    }
  } 

  // CALCULATE AND COMMUNICATE 1-Body PROPAGATOR
  // NOTE: Propagator sorted in a specific format.
  // DO NOT MODIFY THE ORDERING !!!!   
  if(rank()==0) { 

    // check that vMF is real
    for(int i=0; i<vMF.size(); i++) {
      if( std::abs(vMF[i].imag()) > 1e-12 ) {
        app_error()<<" Error: Found complex vMF. \n";
        app_error()<<i <<" " <<vMF[i] <<std::endl;
        APP_ABORT(" Error: Complex vMF. \n\n\n");
      }
    } 

    // to correct for the extra 2 factor of sqrt(dt) from each vn and <vn>_MF term
    if(substractMF) {
      SPValueType one = SPValueType(static_cast<SPValueType>(1.0/dt));
      SPValueType zero = SPValueType(0.0);
      if(sparsePropagator)
        SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn.values(),Spvn.column_data(),Spvn.row_index(),vMF.data(),zero,SPvHS.data());
      else
        DenseMatrixOperators::product_Ax( Dvn.rows(), Dvn.cols(), one, Dvn.values(), Dvn.cols(), vMF.data(), zero, SPvHS.data() ); 
      std::copy(SPvHS.begin(),SPvHS.end(),Hadd.begin());
        
      if(nnodes_per_TG>1) {
        MPI_Status st;
        std::vector<ComplexType> v_(Hadd.size());
        std::vector<int> ranks;
        int pos=0;
        TG.getRanksOfRoots(ranks,pos);
        for(int i=1; i<nnodes_per_TG; i++) {
          myComm->recv(v_.data(),v_.size(),ranks[i],i,TG.getTGCOMM(),&st);
          for(int k=0; k<Hadd.size(); k++) Hadd(k) += v_[k];
        }
      }
    }

    
    if(spinRestricted)   
      for(int k=0; k<NMO*NMO; k++) Hadd(k) += vn0(k);
    else
      for(int k=0; k<2*NMO*NMO; k++) Hadd(k) += vn0(k);
 
    Timer.reset("Generic1");
    Timer.start("Generic1");    
    ham->calculateOneBodyPropagator(cutoff, dt, Hadd, Propg_H1);
    Timer.stop("Generic1");
    app_log()<<" -- Time to calculate one body propagator: " <<Timer.average("Generic1") <<"\n";

    int ni=Propg_H1.size();
    myComm->bcast<int>(&ni,1); 
    myComm->bcast<char>(reinterpret_cast<char*>(Propg_H1.data()), sizeof(s2D<ComplexType>)*Propg_H1.size());

  } else {

    if(substractMF && nnodes_per_TG>1 && TG.getTGNumber()==0 && TG.getCoreRank()==0) {
      SPValueType one = SPValueType(static_cast<SPValueType>(1.0/dt));
      SPValueType zero = SPValueType(0.0);
      if(sparsePropagator)
        SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn.values(),Spvn.column_data(),Spvn.row_index(),vMF.data(),zero,SPvHS.data());
      else
        DenseMatrixOperators::product_Ax( Dvn.rows(), Dvn.cols(), one, Dvn.values(), Dvn.cols(), vMF.data(), zero, SPvHS.data() );
      std::copy(SPvHS.begin(),SPvHS.end(),Hadd.begin());
      std::vector<int> ranks;
      int pos=0;
      TG.getRanksOfRoots(ranks,pos);
      myComm->send(Hadd[0],Hadd.size(),0,pos,TG.getTGCOMM());
    }

    int ni;
    myComm->bcast<int>(&ni,1);
    Propg_H1.resize(ni); 
    myComm->bcast<char>(reinterpret_cast<char*>(Propg_H1.data()), sizeof(s2D<ComplexType>)*Propg_H1.size());

  }

  // write restart if desired
  if(rank()==0) {
    if(hdf_write_file!=std::string("")) {
      hdf_archive dump(myComm);
      if(dump.create(hdf_write_file)) {
        if(!hdf_write(dump,hdf_write_tag)) {
          app_error()<<" Problems writing hdf5 file in phaseless_ImpSamp_ForceBias::setup(). \n";
          return false;
        }
        dump.close();
      } else {
        app_error()<<" Problems opening hdf5 file in phaseless_ImpSamp_ForceBias::setup(). \n";
        APP_ABORT(" Problems opening hdf5 file in phaseless_ImpSamp_ForceBias::setup(). \n");
        return false;
      }
    } 
  } else if(hdf_write_file!=std::string("")) { 
    hdf_archive dump(myComm);
    hdf_write(dump,std::string(""));
  }

  app_log()<<" -- Number of terms in one body propagator: " <<Propg_H1.size() <<"\n";
  app_log()<<" -- Total number of terms in Cholesky vectors: " <<GlobalSpvnSize <<std::endl;
  if(sparsePropagator) 
    app_log()<<" -- Sparsity of Cholesky Matrix, # of Cholesky vectors: " <<double(GlobalSpvnSize)/double(NMO*NMO*Spvn.cols()) <<" " <<Spvn.cols() <<std::endl;


  // setup matrices
  vHS.resize(2*NMO,NMO);
  PHS.resize(2*NMO,NMO);
  vbias.resize(nCholVecs);
  sigma.resize(nCholVecs);
  CV0.resize(nCholVecs);
  for(int i=0; i<vbias.size(); i++) vbias[i]=0;
  for(int i=0; i<sigma.size(); i++) sigma[i]=0;

  // FIX THIS
  ik0=pik0=0;
  ikN = NMO*NMO;
  if(!spinRestricted) ikN *= 2;
  if(ncores_per_TG > 1) {
    int nextra,nvpc = NMO*NMO;
    if(!spinRestricted) nvpc *= 2; 
    nextra = nvpc%ncores_per_TG;
    nvpc /= ncores_per_TG; 
    if(core_rank < nextra ) {
      ik0 = core_rank*(nvpc+1);
      ikN = ik0 + nvpc+1;
    } else {
      ik0 = core_rank*nvpc + nextra;
      ikN = ik0 + nvpc;        
    }
    if(sparsePropagator) {
      pik0 = *(Spvn.row_index()+ik0);
      pikN = *(Spvn.row_index()+ikN);
    } 
  }
  if(parallelPropagation) {
    // used to store quantities local to the TG on a node
    local_buffer.setup(TG.getCoreRank()==0,std::string("localBuffer_")+std::to_string(TG.getTGNumber()),TG.getTGCommLocal());
  }

  // setup temporary storage
  T1.resize(NMO,NAEA);
  T2.resize(NMO,NAEA);

  S1.resize(NMO,NAEA);
  S2.resize(NMO,NAEA);

  // local energy bounds
  dEloc = std::sqrt(2.0/dt);  
  
  if(spinRestricted) {
    Propg_H1_indx[0]=0;
    Propg_H1_indx[1]=Propg_H1.size();
    Propg_H1_indx[2]=0;
    Propg_H1_indx[3]=Propg_H1.size();
  } else {
    Propg_H1_indx[0]=0;
    for(int i=1; i<Propg_H1.size(); i++)
      if( std::get<0>(Propg_H1[i]) < std::get<0>(Propg_H1[i-1]) ) { 
        Propg_H1_indx[1]=i;
        Propg_H1_indx[2]=i;
        Propg_H1_indx[3]=Propg_H1.size()-i;
        break;
      } 
  }

  // these become meaningful only if Cholesky matrix is transposed
  cvec0_loc=0;
  cvecN_loc=cvecN-cvec0;
  if(!save_memory && imp_sampl) {
    app_log()<<" Generating transposed matrix of Cholesky vectors." <<std::endl;
    Timer.reset("Generic1");
    Timer.start("Generic1");
    // distribution is done over cholesky vectors if the propagator is transposed
    if(ncores_per_TG>1) 
      std::tie(cvec0_loc, cvecN_loc) = FairDivideBoundary(core_rank,int(cvecN-cvec0),ncores_per_TG);

    if( sparsePropagator ) { 
      wfn->generateTransposedOneBodyOperator(spinRestricted,"ImportanceSampling",Spvn,SpvnT);
    } else {
      wfn->generateTransposedOneBodyOperator(spinRestricted,"ImportanceSampling",Dvn,DvnT);
    }
    myComm->barrier();
    Timer.stop("Generic1");
    app_log()<<" -- Time to calculate transposed Cholesky matrix: " <<Timer.average("Generic1") <<"\n";
  }

  return true;
}

void phaseless_ImpSamp_ForceBias::Propagate(int steps, int& steps_total, WalkerHandlerBase* wlk, RealType& E1) 
{

  if(fix_bias > 0) {

    int nblk = steps/fix_bias;
    int nextra = steps%fix_bias;
    if(parallelPropagation) {

      for(int i=0; i<nblk; i++, steps_total+=fix_bias) 
        dist_propagation_multiple_steps(fix_bias,wlk,E1);
      if(nextra>0) {
        steps_total+=nextra;
        dist_propagation_multiple_steps(nextra,wlk,E1);
      }
    } else {

      APP_ABORT("Error: fix_bias > 0 only implemented with parallel_propagation = yes. \n\n\n");
      for(int i=0; i<nblk; i++, steps_total+=fix_bias) 
        serial_propagation_multiple_steps(fix_bias,wlk,E1);
      if(nextra>0) {
        steps_total+=nextra;
        serial_propagation_multiple_steps(nextra,wlk,E1);
      }
    }

  } else {

    if(parallelPropagation) {

      for(int i=0; i<steps; i++, steps_total++) 
        dist_propagation_single_step(wlk,E1);

    } else {

      for(int i=0; i<steps; i++, steps_total++) 
        serial_propagation_single_step(wlk,E1);

    }

  }

}

// right now using local energy form of important sampling
void phaseless_ImpSamp_ForceBias::serial_propagation_single_step(WalkerHandlerBase* wset, RealType& Eshift)
{

  const SPComplexType im = SPComplexType(0.0,1.0);
  const SPComplexType halfim = SPComplexType(0.0,0.5);
  int nw = wset->numWalkers(true); // true includes dead walkers in list
  MFfactor.resize(nw);
  hybrid_weight.resize(nw);

  for(int i=0; i<nw; i++) {
    if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue; 

    ComplexType* Sdet = wset->getSM(i);
    wset->setCurrToOld(i);

    S1=ComplexType(0.0); 
    S2=ComplexType(0.0); 

    // propagate forward half a timestep with mean-field propagator
    Timer.start("Propagate::product_SD");
    SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],Sdet,NAEA,S1.data(),NAEA);
    // this will be a problem with DistWalkerHandler in closed_shell case
    SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],Sdet+NAEA*NMO,NAEA,S2.data(),NAEA);
    Timer.stop("Propagate::product_SD");

    // propagate forward full timestep with potential (HS) propagator 

    // 1. sample gaussian field
    Timer.start("Propagate::sampleGaussianFields");  
    sampleGaussianFields(); 
    Timer.stop("Propagate::sampleGaussianFields");  
  
    // 2. calculate force-bias potential. Careful, this gets both alpha and beta  
    Timer.start("Propagate::calculateMixedMatrixElementOfOneBodyOperators");
    if(imp_sampl) {
      SPComplexType* dummy=NULL;
      if(sparsePropagator)
        wfn->calculateMixedMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,Sdet,dummy,Spvn,SpvnT,vbias,!save_memory,true);
      else 
        wfn->calculateMixedMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,Sdet,dummy,Dvn,DvnT,vbias,!save_memory,true);
      apply_bound_vbias();
    }
    Timer.stop("Propagate::calculateMixedMatrixElementOfOneBodyOperators");

    // 3. generate and apply HS propagator (or a good approx of it)
    //  to the S1 and S2 Slater Determinants. The projected determinants are 
    //  returned in S1 and S2. 
    Timer.start("Propagate::applyHSPropagator");
    MFfactor[i]=0.0;
    applyHSPropagator(S1,S2,MFfactor[i],6);  
    Timer.stop("Propagate::applyHSPropagator");

    if(hybrid_method && imp_sampl) {
      SPComplexType tmp = SPComplexType(0.0,0.0);
      for(int ii=0; ii<sigma.size(); ii++) 
        tmp += im*(vMF[ii]-vbias[ii])* ( sigma[ii] - halfim*(vMF[ii]-vbias[ii])  ) ;
      hybrid_weight[i] = tmp;  
    }

    std::fill(Sdet,Sdet+2*NMO*NAEA,ComplexType(0,0));
    // propagate forward half a timestep with mean-field propagator
    Timer.start("Propagate::product_SD");
    SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],S1.data(),NAEA,Sdet,NAEA);
    SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],S2.data(),NAEA,Sdet+NAEA*NMO,NAEA);
    Timer.stop("Propagate::product_SD");

  }  // finished propagating walkers

  Timer.start("Propagate::overlaps_and_or_eloc");
  // calculate overlaps and local energy for all walkers
  if(hybrid_method) { 
    if(imp_sampl)
      wfn->evaluateOverlap("ImportanceSampling",0,wset);
  } else {
      wfn->evaluateLocalEnergyAndOverlap("ImportanceSampling",0,wset);
  }
  Timer.stop("Propagate::overlaps_and_or_eloc");

  if(parallelPropagation) 
    TG.local_barrier();

  // now only head core in TG works
  if(TG.getCoreRank() != 0) return;

  // calculate new weight of walker
  RealType scale = 1.0;
  for(int i=0; i<nw; i++) {
    ComplexType eloc, oldeloc, w0, oa, ob, ooa, oob;
    ComplexType* dummy = wset->getWalker(i,w0,eloc,oa,ob);
    if(!wset->isAlive(i) || std::abs(w0) <= 1e-6) continue;
    wset->getOldWalker(i,oldeloc,ooa,oob);

    scale = 1.0;
    if(hybrid_method) {
      ComplexType ratioOverlaps = ComplexType(1.0,0.0);
      if(imp_sampl) 
        ratioOverlaps = oa*ob/(ooa*oob); 
      //if( (!std::isfinite(ratioOverlaps.real()) || std::abs(oa*ob) < 1e-8) && apply_constrain && imp_sampl ) { 
      if( !std::isfinite(ratioOverlaps.real()) && apply_constrain && imp_sampl ) { 
        scale = 0.0;
        eloc = oldeloc; 
      } else {  
        //scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps*MFfactor[i]) ))):(1.0));
        scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps)-MFfactor[i].imag() ) )):1.0);
        if(imp_sampl) 
          eloc = ( MFfactor[i] - hybrid_weight[i] - std::log(ratioOverlaps) )/dt;
        else 
          eloc = MFfactor[i]/dt;
      }
    } else {
      ComplexType ratioOverlaps = oa*ob/(ooa*oob);
      //if( (!std::isfinite( (ratioOverlaps*MFfactor[i]).real()) || std::abs(oa*ob) < 1e-8) && apply_constrain ) { 
      if( !std::isfinite( (ratioOverlaps*MFfactor[i]).real()) && apply_constrain ) { 
        scale = 0.0;
        eloc = oldeloc; 
      } else  
        scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps)-MFfactor[i].imag() ) )):1.0);
    }

    // temporary fix when FTZ is not set
    if( (!std::isfinite(eloc.real())) || (std::abs(eloc.real())<std::numeric_limits<RealType>::min()) ) {
      scale=0.0;
      eloc=oldeloc;
    } else { 
      eloc = apply_bound_eloc(eloc,Eshift);
    }

    w0 *= ComplexType(scale*std::exp( -dt*( 0.5*( eloc.real() + oldeloc.real() ) - Eshift )),0.0);

    wset->setWalker(i,w0,eloc);

  }  // loop over walkers

}

// since sparse and dense propagation have different optimal settings, 
// treat the layout of TG.buff and CV0 as variables. TG.buff determines
// the marix operation for vbias and CV0 determines the operation for vHS.
// These two are assembled independently, so it can be done.
// The actual assembly of each of these matrices may not be optimal, 
// but for now implement with the efficientcy of the MM in mind.
void phaseless_ImpSamp_ForceBias::dist_propagation_single_step(WalkerHandlerBase* wset, RealType& Eshift)
{
  // structure in TG [hybrid_w, MFfactor, G(1:{2*}NMO*NMO), sigma(1:nCholVecs), vHS(1:{2*}NMO*NMO)] 
  //  1. You either need to send sigma or communicate back vbias, choosing the former 
  //     In principle, you can replace sigma by vbias on the return communication if you need to accumulate vbias

  Timer.start("Propagate::setup");
  const SPComplexType im = SPComplexType(0.0,1.0);
  const SPComplexType halfim = SPComplexType(0.0,0.5);
  int nw = wset->numWalkers(true); // true includes dead walkers in list
  int nw0 = wset->numWalkers(false); 
  MFfactor.resize(nw);
  hybrid_weight.resize(nw);

  int sz = 2 + sizeOfG + nCholVecs + vHS_size; 
  int cnt=0;
  walker_per_node[0] = nw0;
  ComplexType factor;
  int nwglobal = nw0*sz;

  transposed_walker_buffer = (walkerBlock!=1);
  int wstride=1, ax=sz;
  if(transposed_walker_buffer) {
    wstride = nw0;                    
    ax = 1;
  }

  // this routine implies communication inside TG to resize buffers, it also implies a barrier within local TG  
  TG.resize_buffer(nwglobal);
  nwglobal /= sz;
  //local_buffer.resize(nCholVecs*nwglobal);
  local_buffer.resize((cvecN-cvec0)*nwglobal);
  Timer.stop("Propagate::setup");

  Timer.start("Propagate::evalG");  
  // add G to buff 
  wfn->evaluateOneBodyMixedDensityMatrix("ImportanceSampling",wset,(TG.commBuff),sz,2,transposed_walker_buffer,save_memory);
  Timer.stop("Propagate::evalG");  

  // add sigma to buff: This does not need to be communicated, unless you want to keep the option open to
  // store a copy later on 
  Timer.start("Propagate::sampleGaussianFields");
  if(transposed_walker_buffer) {
    if(core_rank==0) {  // doesn't seem useful to split this up
      std::fill( (TG.commBuff)->begin(), (TG.commBuff)->begin()+2*nw0, ComplexType(0,0));
      std::fill( (TG.commBuff)->begin()+nw0*(2+sizeOfG+nCholVecs), (TG.commBuff)->begin()+nw0*sz, ComplexType(0,0));
    }
    cnt=0;
    int ne = (nw0*nCholVecs)%ncores_per_TG, ntpc = (nw0*nCholVecs)/ncores_per_TG;
    int nb = core_rank*ntpc + ((core_rank<ne)?core_rank:ne) ,nt = ((core_rank<ne)?ntpc+1:ntpc);
    sampleGaussianFields( (TG.commBuff)->values()+nw0*(2+sizeOfG)+nb, nt);
  } else {
    cnt=0;
    for(int i=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
      if(cnt%ncores_per_TG != core_rank) {
        cnt++;
        continue;
      }
      sampleGaussianFields( (TG.commBuff)->values()+cnt*sz+2+sizeOfG, nCholVecs); 
      // zero out vHS, and hybrid_w, MFfactor 
      *((TG.commBuff)->values()+cnt*sz) = 0;  
      *((TG.commBuff)->values()+cnt*sz+1) = 0;  
      std::fill( (TG.commBuff)->begin()+cnt*sz+2+sizeOfG+nCholVecs, (TG.commBuff)->begin()+(cnt+1)*sz, ComplexType(0,0)); 
      cnt++;
    }
  }
  Timer.stop("Propagate::sampleGaussianFields");
  Timer.start("Propagate::idle");
  TG.local_barrier();
  Timer.add_time("Propagate::idle");

  Timer.start("Propagate::addvHS");  
  // calculate vHS 
  int currnw=nw0; 
  for(int tgi = 0; tgi<nnodes_per_TG; tgi++) {
    // calculates vbias and accumulates vHS
    addvHS( TG.commBuff ,currnw,sz, wset);   

    // rotate date
    if(distributeSpvn) TG.rotate_buffer(currnw,sz);   
  } 
  Timer.stop("Propagate::addvHS");

  Timer.start("Propagate::idle");
  // synchronize
  TG.local_barrier();
  Timer.add_time("Propagate::idle");

  // propagate walkers now
  if(transposed_walker_buffer) 
    wstride = currnw; 
  cnt=0; 

  Timer.start("Propagate::propg");
    for(int i=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue; 
      
      hybrid_weight[i] = *( (TG.commBuff)->begin() + cnt*ax ); 
      MFfactor[i] = *( (TG.commBuff)->begin() + cnt*ax + wstride ); 
      //MFfactor[i] = std::exp(-MFfactor[i]); 

      if(cnt%ncores_per_TG != core_rank) {
        cnt++;
        continue;
      }

      ComplexType* Sdet = wset->getSM(i);
      wset->setCurrToOld(i);

      S1=ComplexType(0.0); 
      S2=ComplexType(0.0); 

      // propagate forward half a timestep with mean-field propagator
      Timer.start("Propagate::product_SD");
      SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],Sdet,NAEA,S1.data(),NAEA);
      // this will be a problem with DistWalkerHandler in closed_shell case
      SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],Sdet+NAEA*NMO,NAEA,S2.data(),NAEA);
      Timer.stop("Propagate::product_SD");

      Timer.start("Propagate::applyHSPropagator");
      // propagate forward full timestep with potential (HS) propagator 
      // copy vHS from TG.commBuff
#if defined(AFQMC_SP)      
      SPComplexType *spptr = (TG.commBuff)->values()+ ax*cnt + wstride*(2+sizeOfG+nCholVecs);
      ComplexType *ptr = vHS.data(); 
      for(int i=0; i<vHS_size; i++, spptr+=wstride, ptr++) 
        *ptr = static_cast<ComplexType>(*spptr);
#else
      zcopy (vHS_size, (TG.commBuff)->values()+ ax*cnt + wstride*(2+sizeOfG+nCholVecs) , wstride, vHS.data(), 1); 
#endif

      // 3. generate and apply HS propagator (or a good approx of it)
      //  to the S1 and S2 Slater Determinants. The projected determinants are 
      //  returned in S1 and S2. 
      applyHSPropagator(S1,S2,factor,6,false);  
      Timer.stop("Propagate::applyHSPropagator");

      std::fill(Sdet,Sdet+2*NMO*NAEA,ComplexType(0,0));
      // propagate forward half a timestep with mean-field propagator
      Timer.start("Propagate::product_SD");
      SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],S1.data(),NAEA,Sdet,NAEA);
      SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],S2.data(),NAEA,Sdet+NAEA*NMO,NAEA);
      Timer.stop("Propagate::product_SD");

      cnt++;
    }  // finished propagating walkers
  Timer.stop("Propagate::propg");
 
  Timer.start("Propagate::idle");
  TG.local_barrier();
  Timer.add_time("Propagate::idle");

  // in case I implement load balance features in parallelPropagation
  nw = wset->numWalkers(true);

  Timer.start("Propagate::overlaps_and_or_eloc");
  // calculate overlaps and local energy for all walkers
  if(hybrid_method) { 
    if(imp_sampl)
      wfn->evaluateOverlap("ImportanceSampling",0,wset);
  } else {
      wfn->evaluateLocalEnergyAndOverlap("ImportanceSampling",0,wset);
  }
  Timer.stop("Propagate::overlaps_and_or_eloc");

  Timer.start("Propagate::idle");
  if(parallelPropagation) 
    TG.local_barrier();
  Timer.stop("Propagate::idle");

  // now only head core in TG works
  if(TG.getCoreRank() != 0) return;

  Timer.start("Propagate::finish_step");
  // calculate new weight of walker
  RealType scale = 1.0;
  for(int i=0; i<nw; i++) {
    ComplexType eloc, oldeloc, w0, oa, ob, ooa, oob;
    ComplexType* dummy = wset->getWalker(i,w0,eloc,oa,ob);
    if(!wset->isAlive(i) || std::abs(w0) <= 1e-6) continue;
    wset->getOldWalker(i,oldeloc,ooa,oob);

    scale = 1.0;
    if(hybrid_method) {
      ComplexType ratioOverlaps = ComplexType(1.0,0.0);
      if(imp_sampl) 
        ratioOverlaps = oa*ob/(ooa*oob); 
      //if( (!std::isfinite(ratioOverlaps.real()) || std::abs(oa*ob) < 1e-8) && apply_constrain && imp_sampl ) { 
      if( !std::isfinite(ratioOverlaps.real()) && apply_constrain && imp_sampl ) { 
        scale = 0.0;
        eloc = oldeloc; 
      } else {  
        //scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps*MFfactor[i]) ))):(1.0));
        scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps)-MFfactor[i].imag() ) )):1.0);
        if(imp_sampl) 
          eloc = ( MFfactor[i] - hybrid_weight[i] - std::log(ratioOverlaps) )/dt;
        else 
          eloc = MFfactor[i]/dt;
      }
    } else {
      ComplexType ratioOverlaps = oa*ob/(ooa*oob);
      //if( (!std::isfinite( (ratioOverlaps*MFfactor[i]).real()) || std::abs(oa*ob) < 1e-8) && apply_constrain ) { 
      if( !std::isfinite( (ratioOverlaps*MFfactor[i]).real()) && apply_constrain ) { 
        scale = 0.0;
        eloc = oldeloc; 
      } else  
        scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps)-MFfactor[i].imag() ) )):1.0);
    }

    // temporary fix when FTZ is not set
    if( (!std::isfinite(eloc.real())) || (std::abs(eloc.real())<std::numeric_limits<RealType>::min()) ) {
      scale=0.0;
      eloc=oldeloc;
    } else { 
      eloc = apply_bound_eloc(eloc,Eshift);
    }

    w0 *= ComplexType(scale*std::exp( -dt*( 0.5*( eloc.real() + oldeloc.real() ) - Eshift )),0.0);

    wset->setWalker(i,w0,eloc);

  }  // loop over walkers
  Timer.stop("Propagate::finish_step");

}

void phaseless_ImpSamp_ForceBias::dist_propagation_multiple_steps(int n, WalkerHandlerBase* wset, RealType& Eshift)
{

  // structure in TG [hybrid_w(1:n), MFfactor(1:n), G(1:{2*}NMO*NMO), vHS(1:{2*}NMO*NMO, 1:n)] 
  
  Timer.start("Propagate::setup");
  const SPComplexType im = SPComplexType(0.0,1.0);
  const SPComplexType halfim = SPComplexType(0.0,0.5);
  int nw = wset->numWalkers(true); // true includes dead walkers in list
  int nw0 = wset->numWalkers(false); 
  MFfactor.resize(nw);
  hybrid_weight.resize(nw);

  // notice that now I calculate n versions of vHS, hybrid_w, MFfactor
  // random numbers are generated locally
  int sz = 2*n + sizeOfG + vHS_size*n; 
  int cnt=0;
  walker_per_node[0] = nw0;
  ComplexType factor;
  int nwglobal = nw0*sz;

  transposed_walker_buffer = (walkerBlock!=1);
  int wstride=1, ax=sz;
  if(transposed_walker_buffer) {
    wstride = nw0;                    
    ax = 1;
  }

  // this routine implies communication inside TG to resize buffers, it also implies a barrier within local TG  
  TG.resize_buffer(nwglobal);
  nwglobal /= sz;
  // in this case, local_buffer holds:  [ vbias(1:(cvecN-cvec0)*nw), sigma(1,:) ]
  if(transposed_walker_buffer) 
    local_buffer.resize((cvecN-cvec0)*nwglobal*(1+n));
  else 
    local_buffer.resize((cvecN-cvec0)*(nwglobal+n));
  Timer.add_time("Propagate::setup");

  Timer.start("Propagate::evalG");  
  // add G to buff 
  wfn->evaluateOneBodyMixedDensityMatrix("ImportanceSampling",wset,(TG.commBuff),sz,2*n,transposed_walker_buffer,save_memory);
  Timer.stop("Propagate::evalG");  

  // add sigma to buff: This does not need to be communicated, unless you want to keep the option open to
  // store a copy later on 
  Timer.start("Propagate::setup");
  if(transposed_walker_buffer) {
    if(core_rank==0) {  // doesn't seem useful to split this up
      std::fill( (TG.commBuff)->begin(), (TG.commBuff)->begin()+2*n*nw0, ComplexType(0,0));
      std::fill( (TG.commBuff)->begin()+nw0*(2*n+sizeOfG), (TG.commBuff)->begin()+nw0*sz, ComplexType(0,0));
    }
  } else {
    cnt=0;
    for(int i=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
      if(cnt%ncores_per_TG != core_rank) {
        cnt++;
        continue;
      }
      std::fill( (TG.commBuff)->begin()+cnt*sz, (TG.commBuff)->begin()+(cnt+1)*sz, ComplexType(0,0));
      cnt++;
    }
  }
  Timer.stop("Propagate::setup");

  Timer.start("Propagate::idle");  
  TG.local_barrier();
  Timer.add_time("Propagate::idle");  

  Timer.start("Propagate::addvHS");  
  // calculate vHS 
  int currnw=nw0; 
  for(int tgi = 0; tgi<nnodes_per_TG; tgi++) {
    // calculates vbias and accumulates vHS
    addvHS_multiple(n, TG.commBuff ,currnw,sz, wset);   

    // rotate date
    if(distributeSpvn) TG.rotate_buffer(currnw,sz);   
  } 
  Timer.stop("Propagate::addvHS");  

  Timer.start("Propagate::idle");  
  // synchronize
  TG.local_barrier();
  Timer.add_time("Propagate::idle");  

  // in case I implement load balance features in parallelPropagation
  nw = wset->numWalkers(true);

  // propagate walkers now
  if(transposed_walker_buffer) 
    wstride = currnw; 
  for(int ni=0; ni<n; ni++) {
    cnt=0; 
    Timer.start("Propagate::propg");  
    for(int i=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue; 
      
      if(transposed_walker_buffer) { 
        hybrid_weight[i] = *( (TG.commBuff)->begin() + ni*wstride + cnt ); 
        MFfactor[i] = *( (TG.commBuff)->begin() + (n+ni)*wstride + cnt ); 
      } else {
        hybrid_weight[i] = *( (TG.commBuff)->begin() + cnt*sz + ni ); 
        MFfactor[i] = *( (TG.commBuff)->begin() + cnt*sz + n + ni ); 
      }

      if(cnt%ncores_per_TG != core_rank) {
        cnt++;
        continue;
      }

      ComplexType* Sdet = wset->getSM(i);
      wset->setCurrToOld(i);

      S1=ComplexType(0.0); 
      S2=ComplexType(0.0); 

      // propagate forward half a timestep with mean-field propagator
      Timer.start("Propagate::product_SD");
      SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],Sdet,NAEA,S1.data(),NAEA);
      // this will be a problem with DistWalkerHandler in closed_shell case
      SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],Sdet+NAEA*NMO,NAEA,S2.data(),NAEA);
      Timer.stop("Propagate::product_SD");

      // propagate forward full timestep with potential (HS) propagator 
      // copy vHS from TG.commBuff
#if defined(AFQMC_SP)      
      SPComplexType *spptr = (TG.commBuff)->values()+ ax*cnt + wstride*(2+sizeOfG+nCholVecs);
      ComplexType *ptr = vHS.data(); 
      for(int i=0; i<vHS_size; i++, spptr+=wstride, ptr++) 
        *ptr = static_cast<ComplexType>(*spptr);
#else
      if(transposed_walker_buffer) 
        zcopy (vHS_size, (TG.commBuff)->values()+ ni*wstride + cnt + wstride*(2*n+sizeOfG) , wstride*n, vHS.data(), 1); 
      else
        zcopy (vHS_size, (TG.commBuff)->values()+ sz*cnt + (2*n+sizeOfG+ni) , n, vHS.data(), 1); 
#endif

      // 3. generate and apply HS propagator (or a good approx of it)
      //  to the S1 and S2 Slater Determinants. The projected determinants are 
      //  returned in S1 and S2. 
      Timer.start("Propagate::applyHSPropagator");
      applyHSPropagator(S1,S2,factor,6,false);  
      Timer.stop("Propagate::applyHSPropagator");

      std::fill(Sdet,Sdet+2*NMO*NAEA,ComplexType(0,0));
      // propagate forward half a timestep with mean-field propagator
      Timer.start("Propagate::product_SD");
      SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],S1.data(),NAEA,Sdet,NAEA);
      SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],S2.data(),NAEA,Sdet+NAEA*NMO,NAEA);
      Timer.stop("Propagate::product_SD");

      cnt++;
    }  // finished propagating walkers
    Timer.stop("Propagate::propg");  

    Timer.start("Propagate::idle");
    TG.local_barrier();
    Timer.add_time("Propagate::idle");

    Timer.start("Propagate::overlaps_and_or_eloc");
    // calculate overlaps and local energy for all walkers
    if(hybrid_method) { 
      if(imp_sampl)
        wfn->evaluateOverlap("ImportanceSampling",0,wset);
    } else {
        wfn->evaluateLocalEnergyAndOverlap("ImportanceSampling",0,wset);
    }
    Timer.stop("Propagate::overlaps_and_or_eloc");

    Timer.start("Propagate::idle");
    TG.local_barrier();
    Timer.stop("Propagate::idle");

    Timer.start("Propagate::finish_step");
    if(TG.getCoreRank() == 0) {
      // calculate new weight of walker
      RealType scale = 1.0;
      for(int i=0; i<nw; i++) {
        ComplexType eloc, oldeloc, w0, oa, ob, ooa, oob;
        ComplexType* dummy = wset->getWalker(i,w0,eloc,oa,ob);
        if(!wset->isAlive(i) || std::abs(w0) <= 1e-6) continue;
        wset->getOldWalker(i,oldeloc,ooa,oob);
 
        scale = 1.0;
        if(hybrid_method) {
          ComplexType ratioOverlaps = ComplexType(1.0,0.0);
          if(imp_sampl) 
            ratioOverlaps = oa*ob/(ooa*oob); 
          //if( (!std::isfinite(ratioOverlaps.real()) || std::abs(oa*ob) < 1e-8) && apply_constrain && imp_sampl ) { 
          if( !std::isfinite(ratioOverlaps.real()) && apply_constrain && imp_sampl ) { 
            scale = 0.0;
            eloc = oldeloc; 
          } else {  
            //scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps*MFfactor[i]) ))):(1.0));
            scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps)-MFfactor[i].imag() ) )):1.0);
            if(imp_sampl) 
              eloc = ( MFfactor[i] - hybrid_weight[i] - std::log(ratioOverlaps) )/dt;
            else 
              eloc = MFfactor[i]/dt;
          }
        } else {
          ComplexType ratioOverlaps = oa*ob/(ooa*oob);
          //if( (!std::isfinite( (ratioOverlaps*MFfactor[i]).real()) || std::abs(oa*ob) < 1e-8) && apply_constrain ) { 
          if( !std::isfinite( (ratioOverlaps*MFfactor[i]).real()) && apply_constrain ) { 
            scale = 0.0;
            eloc = oldeloc; 
          } else  
            scale = (apply_constrain?(std::max(0.0,std::cos( std::arg(ratioOverlaps)-MFfactor[i].imag() ) )):1.0);
        }

        // temporary fix when FTZ is not set
        if( (!std::isfinite(eloc.real())) || (std::abs(eloc.real())<std::numeric_limits<RealType>::min()) ) {
          scale=0.0;
          eloc=oldeloc;
        } else { 
          eloc = apply_bound_eloc(eloc,Eshift);
        }

        w0 *= ComplexType(scale*std::exp( -dt*( 0.5*( eloc.real() + oldeloc.real() ) - Eshift )),0.0);

        wset->setWalker(i,w0,eloc);

      }  // loop over walkers
    }
    TG.local_barrier();
    Timer.stop("Propagate::finish_step");
    
  }  // loop over steps 

}


void phaseless_ImpSamp_ForceBias::serial_propagation_multiple_steps(int n, WalkerHandlerBase* wset, RealType& Eshift) {
}

// structure in TG [hybrid_w(1:n), MFfactor(1:n), G(1:{2*}NMO*NMO), vHS(1:{2*}NMO*NMO,1:n)]
// NOTE: If walkerBlock>1, buff has been transposed
// local_buffer is used as temporary storage for vbias in the same format as the TG buffer (transposed or not)
// if walkerBlock > 1, the MM algorithm is used and vbias and CV0 are resized. vbias is used as temporaty storage for both vbias and vHS for all walkers
void phaseless_ImpSamp_ForceBias::addvHS_multiple(int nstep, SPComplexSMVector *buff, int nw, int sz, WalkerHandlerBase* wset)
{

  const SPComplexType im = SPComplexType(0.0,1.0);
  const SPComplexType halfim = SPComplexType(0.0,0.5);
  SPComplexType one = SPComplexType(1.0,0.0);
  SPComplexType zero = SPComplexType(0.0,0.0);
  SPValueType vone = SPValueType(1.0);
  SPValueType vzero = SPValueType(0.0);
  int nt = NMO*NMO;
  if(!spinRestricted) nt *= 2; 

  Timer.start("Propagate::addvHS::setup");
  if(TG.getCoreRank()==0) 
    std::fill(local_buffer.begin(),local_buffer.begin()+nw*(cvecN-cvec0),SPComplexType(0.0,0.0));

  TG.local_barrier();

  if(walkerBlock == 1) {
    sigma.resize( static_cast<int>(std::ceil( nstep*(cvecN-cvec0)/ncores_per_TG )) );
    MFs.resize(nstep);
    HWs.resize(nstep);
    if(sparsePropagator && !save_memory) {
      vbias.resize(nstep*std::max(int(ikN-ik0),nCholVecs));
    } else {
      vbias.resize(nstep*std::max(int(ikN-ik0),cvecN-cvec0));
    }
  } else {
    sigma.resize( static_cast<int>(std::ceil( nw*nstep*(cvecN-cvec0)/ncores_per_TG )) );
    MFs.resize(nw*nstep);
    HWs.resize(nw*nstep);
    if(sparsePropagator && !save_memory) {
      vbias.resize(nw*nstep*std::max(int(ikN-ik0),nCholVecs));
    } else {
      vbias.resize(nw*nstep*std::max(int(ikN-ik0),cvecN-cvec0));
    }
  }
  Timer.stop("Propagate::addvHS::setup");

  // Calculate vbias for all walkers.
  //  Careful here:
  //  Dense and Sparse modes operate slightly differently, because in sparse format the value of the cholesky 
  //  vectors on Spvn go from [cvec-:cvecN] and in dense mode they go from [0:cvecN-cvec0]
  SPComplexType *pos = buff->values();
  int vni0=ik0, vniN=ikN;
  if(!save_memory) {
    // sparse uses CV index based on global value
    // dense uses CV index based on local segment
    vni0 = sparsePropagator?(cvec0+cvec0_loc):(cvec0_loc);
    vniN = sparsePropagator?(cvec0+cvecN_loc):(cvecN_loc);
  }
  int shft = 0;
  if(sparsePropagator) shft = cvec0;
  if(imp_sampl) {
    if(walkerBlock == 1) {
      for(int wlk=0; wlk<nw; wlk++, pos+=sz) {

        // calculate vbias for local sector in 'ik' space. This routine assumes a buffer created by wfn->evaluateOneBodyMixedDensityMatrix 
        Timer.start("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");

        if(sparsePropagator)
          wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2*nstep,vni0,vniN,Spvn,SpvnT,vbias,1,1,!save_memory,false);    
        else
          wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2*nstep,vni0,vniN,Dvn,DvnT,vbias,1,1,!save_memory,false);    

        Timer.stop("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");
        Timer.start("Propagate::addvHS::shm_copy");
        if(!save_memory) {
          // if transposed vn (!save_memory), cores work on non-overlapping segments of [cvec0:cvecN]. No need for sync
          zcopy(cvecN_loc-cvec0_loc,vbias.data(),1,local_buffer.values()+wlk*(cvecN-cvec0)+cvec0_loc,1);  
        } else { 
          boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
          BLAS::axpy(cvecN-cvec0,one,vbias.data()+shft,1,local_buffer.values()+wlk*(cvecN-cvec0),1);  
        } 
        Timer.stop("Propagate::addvHS::shm_copy");
      }
    } else {

      Timer.start("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");

      if(sparsePropagator)
        wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2*nw*nstep,vni0,vni0,Spvn,SpvnT,vbias,((walkerBlock<=0)?nw:walkerBlock),nw,!save_memory,false);
      else 
        wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2*nw*nstep,vni0,vniN,Dvn,DvnT,vbias,((walkerBlock<=0)?nw:walkerBlock),nw,!save_memory,false);

      Timer.stop("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");

      Timer.start("Propagate::addvHS::shm_copy");
      if(!save_memory) {
        // if transposed vn (!save_memory), cores work on non-overlapping segments of [cvec0:cvecN]. No need for sync
        zcopy(nw*(cvecN_loc-cvec0_loc),vbias.data(),1,local_buffer.values()+nw*cvec0_loc,1);  
      } else {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
        BLAS::axpy(nw*(cvecN-cvec0),one,vbias.data()+shft*nw,1,local_buffer.values(),1);
      }
      Timer.stop("Propagate::addvHS::shm_copy");
    }
  }

  Timer.start("Propagate::addvHS::idle");
  TG.local_barrier();
  Timer.stop("Propagate::addvHS::idle");

  if(walkerBlock==1) { 

    pos = buff->values();
    for(int wlk=0; wlk<nw; wlk++, pos+=sz) {

      Timer.start("Propagate::addvHS::build_CV0");
      // calculate enough random numbers
      sampleGaussianFields(); 
     // this is where the sigma matrix begins
      SPComplexType* Bmat = local_buffer.values() + nw*(cvecN-cvec0) + shft*nstep; 
      SPRealType* sg = sigma.data(); 
      SPComplexType* vb = local_buffer.values()+wlk*(cvecN-cvec0);
      std::fill(MFs.begin(),MFs.end(),0);
      std::fill(HWs.begin(),HWs.end(),0);
      for(int i=cvec0,ip=0; i<cvecN; i++,vb++) { 
        SPComplexType vbi = apply_bound_vbias(*vb);
        for(int n=0; n<nstep; n++, ip++) {
          if( ip%ncores_per_TG != core_rank ) {
            Bmat++;
            continue;
          }
          *Bmat = *(sg) + im*( vbi - vMF[i] );
          HWs[n] += static_cast<ComplexType>(im*(vMF[i]-vbi)*( *(sg) - halfim*(vMF[i]-vbi) ));
          MFs[n] += static_cast<ComplexType>(*Bmat)*static_cast<ComplexType>(im*vMF[i]); 
          sg++;   Bmat++;
        }
      }
      TG.local_barrier();
      Timer.stop("Propagate::addvHS::build_CV0");
      Timer.start("Propagate::addvHS::build_vHS");
      // vHS
      if(sparsePropagator)
        SparseMatrixOperators::product_SpMatM(int(ikN-ik0), nstep, Spvn.cols(), vone, Spvn.values() + pik0, Spvn.column_data() + pik0, Spvn.row_index()+ik0, local_buffer.values() + nw*(cvecN-cvec0), nstep, vzero, vbias.data(), nstep);
      else
        DenseMatrixOperators::product(int(ikN-ik0), nstep, Dvn.cols(), vone, Dvn.values() + ik0*Dvn.cols(), Dvn.cols(), local_buffer.values() + nw*(cvecN-cvec0), nstep, vzero, vbias.data(), nstep);

      Timer.stop("Propagate::addvHS::build_vHS");
      Timer.start("Propagate::addvHS::shm_copy");

      {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
        BLAS::axpy(nstep,one,HWs.data(),1,pos,1);
        BLAS::axpy(nstep,one,MFs.data(),1,pos+nstep,1);
      }
      // Removed this from lock
      BLAS::axpy(nstep*int(ikN-ik0),one,vbias.data(),1,pos+2*nstep+sizeOfG+ik0*nstep,1);  
      Timer.stop("Propagate::addvHS::shm_copy");

    } 

  } else {

    Timer.start("Propagate::addvHS::build_CV0");
    // in this case, doing all walkers and all steps simultaneously

    // build CV0
    sampleGaussianFields();
    SPRealType* sg =  sigma.data(); 
    SPComplexType* Bmat = local_buffer.values() + nw*(cvecN-cvec0) + shft*nstep*nw;
    SPComplexType* vb = local_buffer.values();
    std::fill(MFs.begin(),MFs.end(),0);
    std::fill(HWs.begin(),HWs.end(),0);
    long ip=0;
    for(int i=cvec0; i<cvecN; i++) {
      SPComplexType vmf_ = im*vMF[i]; 
      for(int n=0; n<nstep; n++, ip++) { 
        for(int wlk=0; wlk<nw; wlk++) {
          if( ip%ncores_per_TG != core_rank ) {
            Bmat++;
            continue;
          }
          SPComplexType vbi = im*apply_bound_vbias(*(vb + i*nw + wlk));
          *Bmat = *(sg) + ( vbi - vmf_ );
          MFs[n*nw + wlk] += static_cast<ComplexType>(*Bmat)*static_cast<ComplexType>(vmf_);
          HWs[n*nw + wlk] += static_cast<ComplexType>((vmf_-vbi)*( *(sg) - static_cast<SPValueType>(0.5)*(vmf_-vbi) ));      
          sg++; Bmat++;  
        }
      }
    } 
    Timer.stop("Propagate::addvHS::build_CV0");
    Timer.start("Propagate::addvHS::build_vHS");
    // vHS: using vbias as temporary storage
    if(sparsePropagator)
      SparseMatrixOperators::product_SpMatM(int(ikN-ik0), nw*nstep, Spvn.cols(), vone, Spvn.values() + pik0, Spvn.column_data() + pik0, Spvn.row_index()+ik0, local_buffer.values() + nw*(cvecN-cvec0), nw*nstep, vzero, vbias.data(), nw*nstep);
    else
      DenseMatrixOperators::product(int(ikN-ik0), nw*nstep, Dvn.cols(), vone, Dvn.values() + ik0*Dvn.cols(), Dvn.cols(), local_buffer.values() + nw*(cvecN-cvec0), nw*nstep, vzero, vbias.data(), nw*nstep);

    Timer.stop("Propagate::addvHS::build_vHS");
    Timer.start("Propagate::addvHS::shm_copy");

    // These are rigurously non-overlapping memory regions, but they are adjacent.
    // Do I still need to control concurrent access??? Probably NOT!!!
    {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
#if defined(AFQMC_SP) 
      SPComplexType *spptr = buff->values();
      // hybrid_weight
      ComplexType *ptr = HWs.data();
      for(int i=0; i<nw*nstep; i++, spptr++, ptr++)
        *spptr += static_cast<SPComplexType>(*ptr);
      // MFfactor
      ptr = MFs.data();
      for(int i=0; i<nw*nstep; i++, spptr++, ptr++)
        *spptr += static_cast<SPComplexType>(*ptr);
#else
      // hybrid_weight
      BLAS::axpy(nw*nstep,one,HWs.data(),1,buff->values(),1);
      // MFfactor
      BLAS::axpy(nw*nstep,one,MFs.data(),1,buff->values()+nw*nstep,1);
#endif
    }
    // remmoved from lock
    BLAS::axpy(nw*nstep*(ikN-ik0),one,vbias.data(),1,buff->values()+nw*(2*nstep+sizeOfG+nstep*ik0),1);  

    Timer.stop("Propagate::addvHS::shm_copy");
  }

  Timer.start("Propagate::addvHS::idle");
  TG.local_barrier();
  Timer.stop("Propagate::addvHS::idle");

}

// structure in TG [hybrid_w, MFfactor, G(1:{2*}NMO*NMO), sigma(1:#CholVects), vHS(1:{2*}NMO*NMO)]
// NOTE: If walkerBlock>1, buff has been transposed
// local_buffer is used as temporary storage for vbias in the same format as the TG buffer (transposed or not)
// if walkerBlock > 1, the MM algorithm is used and vbias and CV0 are resized. vbias is used as temporaty storage for both vbias and vHS for all walkers
void phaseless_ImpSamp_ForceBias::addvHS(SPComplexSMVector *buff, int nw, int sz, WalkerHandlerBase* wset)
{

  Timer.start("Propagate::addvHS::setup");
  const SPComplexType im = SPComplexType(0.0,1.0);
  const SPComplexType halfim = SPComplexType(0.0,0.5);
  SPComplexType one = SPComplexType(1.0,0.0);
  SPComplexType zero = SPComplexType(0.0,0.0);
  SPValueType vone = SPValueType(1.0);
  SPValueType vzero = SPValueType(0.0);
  int nt = NMO*NMO;
  if(!spinRestricted) nt *= 2; 

  if(TG.getCoreRank()==0) 
    std::fill(local_buffer.begin(),local_buffer.end(),SPComplexType(0.0,0.0));

  TG.local_barrier();

  if(walkerBlock != 1) {
    MFs.resize(nw);
    HWs.resize(nw);
    if(sparsePropagator) 
      CV0.resize(nw*nCholVecs);
    else 
      CV0.resize(nw*(cvecN-cvec0));
    if (save_memory)
      if(sparsePropagator)
        vbias.resize(nw*std::max(int(ikN-ik0),nCholVecs));   
      else
        vbias.resize(nw*std::max(int(ikN-ik0),cvecN-cvec0));   
    else 
      vbias.resize(nw*std::max(int(ikN-ik0),cvecN_loc-cvec0_loc));   
  }
  for(int i=0; i<CV0.size(); i++) CV0[i]=zero; 
  Timer.stop("Propagate::addvHS::setup");

  // Calculate vbias for all walkers.
  //  Careful here:
  //  Dense and Sparse modes operate slightly differently, because in sparse format the value of the cholesky 
  //  vectors on Spvn go from [cvec-:cvecN] and in dense mode they go from [0:cvecN-cvec0]
  SPComplexType *pos = buff->values();
  int vni0=ik0, vniN=ikN;
  if(!save_memory) {
    // sparse uses CV index based on global value
    // dense uses CV index based on local segment
    vni0 = sparsePropagator?(cvec0+cvec0_loc):(cvec0_loc);
    vniN = sparsePropagator?(cvec0+cvecN_loc):(cvecN_loc);
  }
  int shft = 0;
  if(sparsePropagator) shft = cvec0;
  if(imp_sampl) {
    if(walkerBlock == 1) {
      for(int wlk=0; wlk<nw; wlk++, pos+=sz) {
        // calculate vbias for local sector in 'ik' space. This routine assumes a buffer created by wfn->evaluateOneBodyMixedDensityMatrix 
        Timer.start("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");

        if(sparsePropagator)
          wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2,vni0,vniN,Spvn,SpvnT,vbias,1,1,!save_memory,false);    
        else
          wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2,vni0,vniN,Dvn,DvnT,vbias,1,1,!save_memory,false);    

        Timer.stop("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");
        Timer.start("Propagate::addvHS::shm_copy");
        if(!save_memory) {
          // if transposed vn (!save_memory), cores work on non-overlapping segments of [cvec0:cvecN]. No need for sync
          zcopy(cvecN_loc-cvec0_loc,vbias.data(),1,local_buffer.values()+wlk*(cvecN-cvec0)+cvec0_loc,1);
        } else {
          boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
          BLAS::axpy(cvecN-cvec0,one,vbias.data()+shft,1,local_buffer.values()+wlk*(cvecN-cvec0),1);  
        } 
        Timer.stop("Propagate::addvHS::shm_copy");
      }
    } else {

      Timer.start("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");

      if(sparsePropagator)
        wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2*nw,vni0,vniN,Spvn,SpvnT,vbias,((walkerBlock<=0)?nw:walkerBlock),nw,!save_memory,false);
      else 
        wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2*nw,vni0,vniN,Dvn,DvnT,vbias,((walkerBlock<=0)?nw:walkerBlock),nw,!save_memory,false);

      Timer.stop("Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer");

      Timer.start("Propagate::addvHS::shm_copy");
      if(!save_memory) {
        // if transposed vn (!save_memory), cores work on non-overlapping segments of [cvec0:cvecN]. No need for sync
        zcopy(nw*(cvecN_loc-cvec0_loc),vbias.data(),1,local_buffer.values()+nw*cvec0_loc,1);
      } else {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
        BLAS::axpy(nw*(cvecN-cvec0),one,vbias.data()+shft*nw,1,local_buffer.values(),1);
      }
      Timer.stop("Propagate::addvHS::shm_copy");

    }
  }

  Timer.start("Propagate::addvHS::idle");
  TG.local_barrier();
  Timer.stop("Propagate::addvHS::idle");

  if(walkerBlock==1) { 

    pos = buff->values();
    for(int wlk=0; wlk<nw; wlk++, pos+=sz) {

      Timer.start("Propagate::addvHS::build_CV0");
      SPComplexType* sg = pos+2+sizeOfG+cvec0;
      SPComplexType* vb = local_buffer.values()+wlk*(cvecN-cvec0);
      ComplexType mf=zero, hw=zero; 
      for(int i=cvec0, j=shft; i<cvecN; i++, sg++,vb++, j++) { 
        SPComplexType vbi = apply_bound_vbias(*vb);
        CV0[j] = *(sg) + im*( vbi - vMF[i] );
        mf += static_cast<ComplexType>(CV0[j])*static_cast<ComplexType>(im*vMF[i]); 
        hw += static_cast<ComplexType>(im*(vMF[i]-vbi)*( *(sg) - halfim*(vMF[i]-vbi) ));
      }
      Timer.stop("Propagate::addvHS::build_CV0");
      Timer.start("Propagate::addvHS::build_vHS");
      // vHS
      if(sparsePropagator)
        SparseMatrixOperators::product_SpMatV(int(ikN-ik0), Spvn.cols(), vone, Spvn.values() + pik0, Spvn.column_data() + pik0, Spvn.row_index()+ik0, CV0.data(), vzero, SPvHS.data()+ik0);
      else
        DenseMatrixOperators::product_Ax(int(ikN-ik0), Dvn.cols(), vone, Dvn.values() + ik0*Dvn.cols(), Dvn.cols(), CV0.data(), vzero, SPvHS.data()+ik0); 

      Timer.stop("Propagate::addvHS::build_vHS");
      Timer.start("Propagate::addvHS::shm_copy");

      // These are rigurously non-overlapping memory regions, but they are adjacent.
      // Do I still need to control concurrent access??? Probably NOT!!!
      if(TG.getCoreRank()==0)
      {
        //boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
        // hybrid_weight
        *pos += hw;
        // MFfactor
        *(pos+1) += mf;
      } 
      BLAS::axpy(int(ikN-ik0),one,SPvHS.data()+ik0,1,pos+2+sizeOfG+nCholVecs+ik0,1);  
      Timer.stop("Propagate::addvHS::shm_copy");

    } 

  } else {

    Timer.start("Propagate::addvHS::build_CV0");

    // build CV0
    SPComplexType* sg = buff->values()+nw*(2+sizeOfG+cvec0);
    SPComplexType* vb = local_buffer.values();
    std::fill(MFs.begin(),MFs.end(),0);
    std::fill(HWs.begin(),HWs.end(),0);
    for(int i=cvec0, j=shft*nw; i<cvecN; i++) {
      SPComplexType vmf_ = im*vMF[i]; 
      for(int wlk=0; wlk<nw; wlk++,j++,sg++,vb++) {
        SPComplexType vbi = im*apply_bound_vbias(*vb);
        CV0[j] = *(sg) + ( vbi - vmf_ );
        MFs[wlk] += CV0[j]*vmf_;
        HWs[wlk] += (vmf_-vbi)*( *(sg) - static_cast<SPValueType>(0.5)*(vmf_-vbi) );      
      }
    } 
    Timer.stop("Propagate::addvHS::build_CV0");

    Timer.start("Propagate::addvHS::build_vHS");
    // vHS: using vbias as temporary storage
    if(sparsePropagator)
      SparseMatrixOperators::product_SpMatM(int(ikN-ik0), nw, Spvn.cols(), vone, Spvn.values() + pik0, Spvn.column_data() + pik0, Spvn.row_index()+ik0, CV0.data(), nw, vzero, vbias.data(), nw);
    else
      DenseMatrixOperators::product(int(ikN-ik0), nw, Dvn.cols(), vone, Dvn.values() + ik0*Dvn.cols(), Dvn.cols(), CV0.data(), nw, vzero, vbias.data(), nw);

    Timer.stop("Propagate::addvHS::build_vHS");

    Timer.start("Propagate::addvHS::shm_copy");

    // These are rigurously non-overlapping memory regions, but they are adjacent.
    // Do I still need to control concurrent access??? Probably NOT!!!
    if(TG.getCoreRank()==0)
    {
    //  boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
#if defined(AFQMC_SP) 
      SPComplexType *spptr = buff->values();
      // hybrid_weight
      ComplexType *ptr = HWs.data();
      for(int i=0; i<nw; i++, spptr++, ptr++)
        *spptr += static_cast<SPComplexType>(*ptr);
      // MFfactor
      ptr = MFs.data();
      for(int i=0; i<nw; i++, spptr++, ptr++)
        *spptr += static_cast<SPComplexType>(*ptr);
#else
      // hybrid_weight
      BLAS::axpy(nw,one,HWs.data(),1,buff->values(),1);
      // MFfactor
      BLAS::axpy(nw,one,MFs.data(),1,buff->values()+nw,1);
#endif
    }
    BLAS::axpy(nw*(ikN-ik0),one,vbias.data(),1,buff->values()+nw*(2+sizeOfG+nCholVecs+ik0),1);  

    Timer.stop("Propagate::addvHS::shm_copy");
  }

  Timer.start("Propagate::addvHS::idle");
  TG.local_barrier();
  Timer.stop("Propagate::addvHS::idle");

}

void phaseless_ImpSamp_ForceBias::applyHSPropagator(ComplexMatrix& M1, ComplexMatrix& M2, ComplexType& factor,  int order, bool calculatevHS)
{

  if(order < 0) order = Order_Taylor_Expansion;

  const SPComplexType im = SPComplexType(0.0,1.0);
  ComplexType zero = ComplexType(0.0,0.0); 

  if(calculatevHS) {
    factor = zero;

    Timer.start("Propagate::addvHS::build_vHS");
    for(int i=0; i<sigma.size(); i++) {
      CV0[i] = sigma[i] + im*(vbias[i]-vMF[i]);
      factor += CV0[i]*im*vMF[i];
    }

    SPComplexType *vHSptr = SPvHS.data(); 
#ifndef AFQMC_SP
    vHSptr = vHS.data(); 
#endif

    if(sparsePropagator)
      SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),SPValueType(1),Spvn.values(),Spvn.column_data(),Spvn.row_index(),CV0.data(),SPValueType(0),vHSptr);
    else
      DenseMatrixOperators::product_Ax(Dvn.rows(),Dvn.cols(),SPValueType(1),Dvn.values(),Dvn.cols(),CV0.data(),SPValueType(0),vHSptr);

#if defined(AFQMC_SP)
    // if working with single precision, copy to vHS
    SPComplexType *spptr = SPvHS.data(); 
    ComplexType *ptr = vHS.data(); 
    for(int i=0; i<vHS_size; i++, spptr++, ptr++)
      *ptr = static_cast<ComplexType>(*spptr);
#endif

    Timer.stop("Propagate::addvHS::build_vHS");
    //factor = std::exp(-factor);  // not exponentiated anymore
  }

  // calculate exp(vHS)*S through a Taylor expansion of exp(vHS)
  Timer.start("Propagate::apply_expvHS_Ohmms");
  T1=M1;
  for(int n=1; n<=order; n++) {
    ComplexType fact = ComplexType(0.0,1.0)*static_cast<ComplexType>(1.0/static_cast<double>(n)); 
    DenseMatrixOperators::product(NMO,NAEA,NMO,fact,vHS.data(),NMO,T1.data(),NAEA,zero,T2.data(),NAEA);
    T1  = T2;
    M1 += T1;
  }

  //if(closed_shell) {
  //  Timer.stop("Propagate::apply_expvHS_Ohmms");
  //  M2=M1;
  //  return;
  //}   

//  if M1 == M2 on entry, no need to do this in spinRestricted case
  int disp = ((spinRestricted)?0:NMO*NMO);
  T1 = M2;
  for(int n=1; n<=order; n++) {
    ComplexType fact = ComplexType(0.0,1.0)*static_cast<ComplexType>(1.0/static_cast<double>(n));
    DenseMatrixOperators::product(NMO,NAEB,NMO,fact,vHS.data()+disp,NMO,T1.data(),NAEA,zero,T2.data(),NAEA);
    T1  = T2;
    M2 += T1;
  } 
  Timer.stop("Propagate::apply_expvHS_Ohmms");
}

void phaseless_ImpSamp_ForceBias::sampleGaussianFields()
{

  int n = sigma.size();
  for (int i=0; i+1<n; i+=2)
  {
    SPRealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    SPRealType mag = std::sqrt(-2.0*std::log(temp1));
    sigma[i]  =mag*std::cos(6.283185306*temp2);
    sigma[i+1]=mag*std::sin(6.283185306*temp2);
  }
  if (n%2==1)
  {
    SPRealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    sigma[n-1]=std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2);
  }

}

void phaseless_ImpSamp_ForceBias::sampleGaussianFields(SPComplexType* v, int n)
{ 
  for (int i=0; i+1<n; i+=2)
  { 
    SPRealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    SPRealType mag = std::sqrt(-2.0*std::log(temp1));
    *(v++)  = SPComplexType(mag*std::cos(6.283185306*temp2),0);
    *(v++)  = SPComplexType(mag*std::sin(6.283185306*temp2),0);
  }
  if (n%2==1)
  { 
    SPRealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    *v = SPComplexType(std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2),0);
  }

}

void phaseless_ImpSamp_ForceBias::benchmark(std::string& blist, int maxnW, int delnW, int nrepeat, WalkerHandlerBase* wlk)
{

  SPValueType vone = SPValueType(1.0);
  SPValueType vzero = SPValueType(0.0);
  SPComplexType cone = SPComplexType(1.0);
  SPComplexType czero = SPComplexType(0.0);

  // time different components of the calculation as a function of the number of walkers
  // 1. Propagate
  if(blist.find("propg")!=std::string::npos) {
    std::ofstream out;
    if(myComm->rank() == 0) 
      out.open(myComm->getName()+".propagate.dat");

    std::vector<std::string> tags(17); 
    tags[0]="Benchmark_tot";
    tags[1]="Benchmark_wait"; 
    tags[2]="Propagate::setup";
    tags[3]="Propagate::evalG"; 
    tags[4]="Propagate::sampleGaussianFields"; 
    tags[5]="Propagate::addvHS"; 
    tags[6]="Propagate::propg"; 
    tags[7]="Propagate::overlaps_and_or_eloc";
    tags[8]="Propagate::finish_step";
    tags[9]="Propagate::idle";
    tags[10]="Propagate::addvHS::setup";
    tags[11]="Propagate::addvHS::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer";
    tags[12]="Propagate::addvHS::build_CV0"; 
    tags[13]="Propagate::addvHS::shm_copy"; 
    tags[14]="Propagate::addvHS::build_vHS"; 
    tags[15]="Propagate::addvHS::idle";
    tags[16]="PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer::setup";
    out<<"# #wlk  #steps "; 
    for( std::string& str: tags) out<<str <<" ";
    out<<std::endl;
 
    int dummy,nstp = (fix_bias>0)?fix_bias:1; 
    double Es=wlk->getEloc(0).real();
    // call it once to set timers
    wlk->resetNumberWalkers(1);
    Propagate(nstp,dummy,wlk,Es);
    Timer.add("Benchmark_tot");
    Timer.add("Benchmark_wait");
    int nw=1;
    while(nw <= maxnW) {

      wlk->resetNumberWalkers(nw);
      Timer.reset_all(); 
      double t1, t2;
      myComm->barrier();
      for(int i=0; i<nrepeat; i++) {
        Timer.start("Benchmark_tot");
        Propagate(nstp,dummy,wlk,Es);
        Timer.stop("Benchmark_tot");
        Timer.start("Benchmark_wait");
        myComm->barrier();
        Timer.stop("Benchmark_wait");
      }     
      if(myComm->rank() == 0) { 
        out<<nw <<" " <<nstp <<" ";
        for( std::string& str: tags) out<<Timer.total(str)/double(nrepeat) <<" ";
        out<<std::endl;
      }

      if(delnW <=0) nw *= 2;
      else nw += delnW;
    }

    if(myComm->rank() == 0) out.close(); 

  } 

  // 2. MatM: if sparse, do both.
  //          if dense, do dense only
  if(blist.find("spmm")!=std::string::npos && sparsePropagator) {
    std::ofstream out;
    if(myComm->rank() == 0)
      out.open(myComm->getName()+".spmm.dat");

    std::vector<std::string> tags(3);
    tags[0]="SpMatM_CMat*Vbias";
    tags[1]="SpMatTM_Cmat*G";
    tags[2]="SpMatM_TCmat*G_compact";
    out<<"# #wlk "; 
    for( std::string& str: tags) out<<str <<" ";
    out<<std::endl;

    for( std::string& str: tags) Timer.reset(str);
    
    int nw=1;
    while(nw <= maxnW) {

      ComplexMatrix G(NMO*NMO,nw);
      ComplexMatrix vb(nCholVecs,nw);
      Timer.reset_all();
      double t1, t2;
      for(int i=0; i<nrepeat; i++) {
        if(nw > 1) {
          Timer.start("SpMatM_CMat*Vbias");
          SparseMatrixOperators::product_SpMatM(Spvn.rows(), nw, Spvn.cols(), vone, Spvn.values(), Spvn.column_data(), Spvn.row_index(), vb.data(), nw, vzero, G.data(), nw);
          Timer.stop("SpMatM_CMat*Vbias");
          Timer.start("SpMatTM_Cmat*G");
          SparseMatrixOperators::product_SpMatTM(Spvn.rows(), nw, Spvn.cols(), vone, Spvn.values(), Spvn.column_data(), Spvn.row_index(), G.data(), nw, vzero, vb.data(), nw);
          Timer.stop("SpMatTM_Cmat*G");
          if(!save_memory) {
            Timer.start("SpMatM_TCmat*G_compact");
            SparseMatrixOperators::product_SpMatM(SpvnT.rows(), nw, SpvnT.cols(), cone, SpvnT.values(), SpvnT.column_data(), SpvnT.row_index(), G.data(), nw, czero, vb.data(), nw);
            Timer.stop("SpMatM_TCmat*G_compact"); 
          }
        } else {
          Timer.start("SpMatM_CMat*Vbias");
          SparseMatrixOperators::product_SpMatV(Spvn.rows(), Spvn.cols(), vone, Spvn.values(), Spvn.column_data(), Spvn.row_index(), vb.data(), vzero, G.data());
          Timer.stop("SpMatM_CMat*Vbias");
          Timer.start("SpMatTM_Cmat*G");
          SparseMatrixOperators::product_SpMatTV(Spvn.rows(), Spvn.cols(), vone, Spvn.values(), Spvn.column_data(), Spvn.row_index(), G.data(), vzero, vb.data());
          Timer.stop("SpMatTM_Cmat*G");
          if(!save_memory) {
            Timer.start("SpMatM_TCmat*G_compact");
            SparseMatrixOperators::product_SpMatV(SpvnT.rows(), SpvnT.cols(), cone, SpvnT.values(), SpvnT.column_data(), SpvnT.row_index(), G.data(), czero, vb.data());
            Timer.stop("SpMatM_TCmat*G_compact"); 
          }
        }
      }
      if(myComm->rank() == 0) {
        out<<nw <<" " ;
        for( std::string& str: tags) out<<Timer.total(str)/double(nrepeat) <<" ";
        out<<std::endl;
      }

      if(delnW <=0) nw *= 2;
      else nw += delnW;
    }

    if(myComm->rank() == 0) out.close(); 
    myComm->barrier();
 
  }  
  
  // 4. evalG
  
  // 5. eloc
  // 6. distMM

}

}
