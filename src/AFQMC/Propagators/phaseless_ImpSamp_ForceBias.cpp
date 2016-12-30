
#include<algorithm>

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

#include"AFQMC/Numerics/SparseMatrixOperations.h"
#include"AFQMC/Numerics/DenseMatrixOperations.h"
#include "Utilities/RandomGenerator.h"
#include "AFQMC/Utilities/Utils.h"

namespace qmcplusplus
{

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
    m_param.add(test_library,"test","int");
    m_param.add(sub,"substractMF","std::string");
    m_param.add(use,"useCholesky","std::string");
    m_param.add(cutoff,"cutoff_propg","double");
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

    std::string par("no");
    m_param.add(par,"paral_fac","std::string");
    m_param.add(par,"parallel_fac","std::string");
    m_param.add(par,"parallel_factorization","std::string");

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

    if(substractMF) 
      app_log()<<"Using mean-field substraction in propagator. \n";

    if(parallel_factorization)
      app_log()<<"Calculating factorization of 2-body hamiltonian in parallel. \n";

    if(nnodes_per_TG > 1 && !parallel_factorization) {
      parallel_factorization=true;
      app_log()<<"Distributed propagator (nnodes_per_TG(propagator)>1) requires a parallel factorization. Setting parallel_factorization to true. \n"; 
      app_log()<<"Calculating factorization of 2-body hamiltonian in parallel. \n";
    }

    if(use_eig)
      app_log()<<"Calculating factorization of 2 body interaction with direct diagonalization.\n";
    else 
      app_log()<<"Calculating factorization of 2 body interaction with Cholesky method.\n"; 

    if(!imp_sampl) {
      app_log()<<"AFQMC propagation WITHOUT Importance Sampling! " <<std::endl; 
      app_log()<<"CAREFUL Incomplete implementation of propagation WITHOUT Importance Sampling! " <<std::endl; 
      app_log()<<"CAREFUL Weights are wrong and energies are incorrect !!! " <<std::endl; 
    }

    if(hybrid_method)
      app_log()<<"Using hybrid method to calculate the weights during the propagation." <<std::endl;

    app_log()<<" Running Propagator with " <<nnodes_per_TG <<" nodes per task group. \n"; 

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }

    return true;

}

bool phaseless_ImpSamp_ForceBias::hdf_write(hdf_archive& dump,const std::string& tag)
{

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

    std::vector<int> Idata(4);
    int ntot=Spvn.size(), ntot_=Spvn.size();
    if(nnodes_per_TG>1) 
      MPI_Reduce(&ntot_, &ntot, 1, MPI_INT, MPI_SUM, 0, TG.getTGCOMM()); 
    Idata[0]=ntot;  // this needs to be summed over relevant nodes
    Idata[1]=Spvn.rows();  // same here
    Idata[2]=nCholVecs;
    Idata[3]=NMO;

    dump.push("Propagators");
    dump.push("phaseless_ImpSamp_ForceBias");
    if(tag != std::string("")) dump.push(tag);
    dump.write(Idata,"Spvn_dims");

    // transpose matrix for easier access to Cholesky vectors 
    Spvn.transpose();   

    Idata.resize(nCholVecs);
    for(int i=0; i<nCholVecs; i++)
      Idata[i] = *(Spvn.rowIndex_begin()+i+1) - *(Spvn.rowIndex_begin()+i);
    std::vector<int> sz_local;
    sz_local = Idata; 
    if(nnodes_per_TG>1) myComm->gsum(Idata,TG.getTGCOMM());

    dump.write(Idata,std::string("Spvn_nterms_per_vector"));

    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<ComplexType> vvec;
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
    for(ComplexSMSpMat::indxType i=0; i<nCholVecs; i++) {

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

    return true;

  } else if(nnodes_per_TG>1 && TG.getTGNumber()==0 && TG.getCoreRank()==0) {

    RealType scale = 1.0/std::sqrt(dt);
    Spvn *= scale;

    // ranks of other heads in TG
    std::vector<int> ranks;
    int pos=0;
    TG.getRanksOfRoots(ranks,pos);
    npr=ranks.size();

    int ntot=Spvn.size(), ntot_=Spvn.size();
    MPI_Reduce(&ntot_, &ntot, 1, MPI_INT, MPI_SUM, 0, TG.getTGCOMM()); 

    Spvn.transpose();

    std::vector<int> Idata(nCholVecs);
    for(int i=0; i<nCholVecs; i++)
      Idata[i] = *(Spvn.rowIndex_begin()+i+1) - *(Spvn.rowIndex_begin()+i);
    std::vector<int> sz_local;
    sz_local = Idata;
    myComm->gsum(Idata,TG.getTGCOMM());

    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<ComplexType> vvec;
    ivec.reserve(nmax);
    vvec.reserve(nmax);

    for(ComplexSMSpMat::indxType i=cvec0; i<cvecN; i++) {

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

    return true;

  } else if(nnodes_per_TG>1 && TG.getTGNumber()==0) {

    int ntot=0, ntot_=0;
    if(nnodes_per_TG>1)
      MPI_Reduce(&ntot_, &ntot, 1, MPI_INT, MPI_SUM, 0, TG.getTGCOMM()); 

    std::vector<int> Idata;
    Idata.resize(nCholVecs);
    std::fill(Idata.begin(),Idata.end(),0);
    myComm->gsum(Idata,TG.getTGCOMM());

  } else {
    return true;
  } 
 
  return false;
}

bool phaseless_ImpSamp_ForceBias::hdf_read(hdf_archive& dump,const std::string& tag)
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

    Idata.resize(nCholVecs);
    if(!dump.read(Idata,"Spvn_nterms_per_vector")) return false;
    myComm->bcast(Idata,MPI_COMM_HEAD_OF_NODES);    

    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<ComplexType> vvec;
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

    ComplexSMSpMat::iterator itv = Spvn.vals_begin();
    ComplexSMSpMat::int_iterator itc = Spvn.cols_begin();
    ComplexSMSpMat::int_iterator itr = Spvn.rows_begin();
    for(ComplexSMSpMat::indxType i=0; i<nCholVecs; i++) {  

      ivec.resize(Idata[i]); 
      vvec.resize(Idata[i]); 
      if(!dump.read(ivec,std::string("Spvn_index_")+std::to_string(i))) return false;
      if(!dump.read(vvec,std::string("Spvn_vals_")+std::to_string(i))) return false;           

      myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<RealType>(reinterpret_cast<RealType*>(vvec.data()),2*vvec.size(),MPI_COMM_HEAD_OF_NODES);

      if(nnodes_per_TG == 1 || (i>=cvec0 && i < cvecN)) {
        for(int k=0; k<ivec.size(); k++) {
          *(itv++) = vvec[k]; 
          *(itr++) = ivec[k]; 
          *(itc++) = i; 
        }
      }

    }

    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();

    Spvn.compress();

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

    Idata.resize(nCholVecs);
    myComm->bcast(Idata,MPI_COMM_HEAD_OF_NODES);
    int nmax = *std::max_element(Idata.begin(),Idata.end());
    std::vector<IndexType> ivec;
    std::vector<ComplexType> vvec;
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

    ComplexSMSpMat::iterator itv = Spvn.vals_begin();
    ComplexSMSpMat::int_iterator itc = Spvn.cols_begin();
    ComplexSMSpMat::int_iterator itr = Spvn.rows_begin();
    for(ComplexSMSpMat::indxType i=0; i<nCholVecs; i++) {

      ivec.resize(Idata[i]);
      vvec.resize(Idata[i]);

      myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<RealType>(reinterpret_cast<RealType*>(vvec.data()),2*vvec.size(),MPI_COMM_HEAD_OF_NODES);

      if(nnodes_per_TG == 1 || (i>=cvec0 && i < cvecN)) {
        for(int k=0; k<ivec.size(); k++) {
          *(itv++) = vvec[k]; 
          *(itr++) = ivec[k];
          *(itc++) = i;      
        }
      }

    }

    Spvn.compress();

    // scale back by std::sqrt(dt) 
    RealType scale = std::sqrt(dt); 
    Spvn *= scale; 

  } else {

    myComm->bcast(Idata);
    int ntot = Idata[0];
    int nrows = Idata[1];
    nCholVecs = Idata[2];

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

    Spvn.setCompressed();

  } 

  myComm->barrier();

  return true;
}



bool phaseless_ImpSamp_ForceBias::setup(std::vector<int>& TGdata, ComplexSMVector* v, HamiltonianBase* ham,WavefunctionHandler* w, RealType dt_, hdf_archive& dump_read, const std::string& hdf_restart_tag, MPI_Comm tg_comm, MPI_Comm node_comm)
{
  dt = dt_;

  wfn=w;
  sizeOfG = 0;

  ncores_per_TG=TGdata[4];
  if(!TG.quick_setup(ncores_per_TG,nnodes_per_TG,TGdata[0],TGdata[1],TGdata[2],TGdata[3]))
    return false;
  TG.setBuffer(v);
  core_rank = TG.getCoreRank();
  TG.setNodeCommLocal(node_comm);
  TG.setTGCommLocal(tg_comm);

  walker_per_node.resize(nnodes_per_TG);

  parallelPropagation = (nnodes_per_TG>1 || ncores_per_TG>1);
  distributeSpvn = nnodes_per_TG>1;

  if(parallelPropagation)
    sizeOfG = wfn->sizeOfInfoForDistributedPropagation("ImportanceSampling");

  bool read_Spvn_from_file=false;
  Propg_H1_indx.resize(4); 
  ComplexMatrix Hadd(2*NMO,NMO);
  for(int i=0; i<2*NMO; i++) 
   for(int j=0; j<NMO; j++) 
    Hadd(i,j)=0;

  Spvn.setup(head_of_nodes,"Spvn",TG.getNodeCommLocal());

// FIX FIX FIXL must define
// nCholVecs 
// cvec0, cvecN

  nCholVec_per_node.resize(nnodes_per_TG);
 
  // Only master tries to read 
  if(myComm->rank() == 0) {

    if(hdf_read_file!=std::string("")) {

      hdf_archive dump(myComm);
      if(dump.open(hdf_read_file,H5F_ACC_RDONLY,false)) { 
        read_Spvn_from_file = hdf_read(dump,hdf_read_tag);
        dump.close();
        if(read_Spvn_from_file) 
          app_log()<<"Successfully read HS potentials from file: " <<hdf_read_file <<"\n";
      } else {
        APP_ABORT("Error in phaseless_ImpSamp_ForceBias::setup(): problems opening hdf5 file. \n");
      }

    } 

  } else {

    if(hdf_read_file!=std::string("")) { 
      hdf_archive dump(myComm);
      read_Spvn_from_file = hdf_read(dump,std::string(""));
    }

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
        ham->calculateHSPotentials_Diagonalization(cutoff,dt,Spvn,TG,nCholVec_per_node,false);
      } else
        ham->calculateHSPotentials(cutoff, dt, Spvn,TG,nCholVec_per_node,false);

      Timer.stop("Generic1");
      app_log()<<" -- Time to calculate HS potentials: " <<Timer.average("Generic1") <<"\n";

      std::vector<int> ni(3);
      ni[0]=Spvn.size();
      ni[1]=Spvn.rows();
      ni[2]=Spvn.cols();
      myComm->bcast(ni);        

      nCholVecs = Spvn.cols();
      cvec0 = 0;
      cvecN = nCholVecs;

      myComm->bcast<RealType>(reinterpret_cast<RealType*>(Spvn.values()),2*Spvn.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<int>(Spvn.row_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<int>(Spvn.column_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<int>(Spvn.row_index(),Spvn.rows()+1,MPI_COMM_HEAD_OF_NODES);

      myComm->barrier();
      myComm->barrier();

      GlobalSpvnSize = Spvn.size();

    } else {

      std::vector<int> ni(3);
      myComm->bcast(ni);
      Spvn.setDims(ni[1],ni[2]);
    
      nCholVecs = ni[2];
      cvec0 = 0;
      cvecN = nCholVecs;

      if(head_of_nodes) {
        Spvn.allocate_serial(ni[0]);
        Spvn.resize_serial(ni[0]);

        myComm->bcast<RealType>(reinterpret_cast<RealType*>(Spvn.values()),2*Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.row_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.column_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.row_index(),Spvn.rows()+1,MPI_COMM_HEAD_OF_NODES);
      }
      Spvn.setCompressed();

      myComm->barrier();
      if(!head_of_nodes) Spvn.initializeChildren();
      myComm->barrier();
 
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
      ham->calculateHSPotentials_Diagonalization(cutoff,dt,Spvn,TG,nCholVec_per_node,true);
    } else
      ham->calculateHSPotentials(cutoff, dt, Spvn,TG,nCholVec_per_node,true);

    Timer.stop("Generic1");
    app_log()<<" -- Time to calculate HS potentials: " <<Timer.average("Generic1") <<"\n";

    nCholVecs = Spvn.cols();
    cvec0 = 0;
    cvecN = nCholVecs;
    if(nnodes_per_TG>1) {
      int cnt=0,node_number = TG.getLocalNodeNumber();
      if(node_number==0)
        cvec0=0;
      else
        cvec0 = std::accumulate(nCholVec_per_node.begin(),nCholVec_per_node.begin()+node_number,0);
      cvecN = cvec0+nCholVec_per_node[node_number];
    }
    // collect this from head processors in TG=0
    GlobalSpvnSize=Spvn.size();
  } 

  int nt = Spvn.capacity();
  app_log()<<"Memory used by HS potential: " <<(Spvn.capacity()*sizeof(ValueType)+(2*Spvn.capacity()+NMO*NMO+1)*sizeof(int)+8000)/1024.0/1024.0 <<" MB " <<std::endl;
  RealType vnmax1=0,vnmax2=0;
  for(int i=0; i<Spvn.size(); i++) {
    if( std::abs( *(Spvn.values()+i) ) > vnmax1 )
      vnmax1 = std::abs(*(Spvn.values()+i)); 
    if( (Spvn.values()+i)->imag()  > vnmax2 )
      vnmax2 = (Spvn.values()+i)->imag(); 
  }
  app_log()<<" Largest term in Vn: " <<vnmax1 <<" " <<vnmax2 <<std::endl; 

  hybrid_weight.reserve(1000);
  MFfactor.reserve(1000);
  vMF.resize(Spvn.cols());
  for(int i=0; i<vMF.size(); i++) vMF[i]=0;

  // vMF is actually i*sqrt(dt)*<vn>, since Spvn is also scaled by i*sqrt(dt)
  if(substractMF) { 
    wfn->calculateMeanFieldMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,Spvn,vMF);
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

    // to correct for the extra 2 factor of i*sqrt(dt) from each vn and <vn>_MF term
    if(substractMF) {
      ComplexType one = ComplexType(-1.0/dt,0.0);
      ComplexType zero = ComplexType(0.0,0.0);
      SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn,vMF.data(),zero,Hadd.data());
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
      ComplexType one = ComplexType(-1.0/dt,0.0);
      ComplexType zero = ComplexType(0.0,0.0);
      SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn,vMF.data(),zero,Hadd.data());
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

  // setup matrices
  vHS.resize(2*NMO,NMO);
  PHS.resize(2*NMO,NMO);
  vbias.resize(Spvn.cols());
  sigma.resize(Spvn.cols());
  CV0.resize(Spvn.cols());
  for(int i=0; i<vbias.size(); i++) vbias[i]=0;
  for(int i=0; i<sigma.size(); i++) sigma[i]=0;

  // FIX THIS
  ik0=pik0=0;
  ikN = NMO*NMO;
  if(!spinRestricted) ikN *= 2;
  if(nnodes_per_TG > 1) {
    app_log()<<"  WARNING: Disabling save_memory in propagator because nnodes_per_TG. FIX FIX FIX later. \n";
    save_memory = true; // until I fix distributed SpvnT
  }
  if(ncores_per_TG > 1) {
    app_log()<<"  WARNING: Disabling save_memory in propagator because ncores_per_TG. FIX FIX FIX later. \n";
    app_log()<<"  WARNING: Simple uniform distribution of Spvn over cores. Might not be efficient. FIX FIX FIX \n"; 
    save_memory = true; // until I fix SpvnT
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
    pik0 = *(Spvn.row_index()+ik0);
    pikN = *(Spvn.row_index()+ikN);
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

  if(test_library) test_linear_algebra();

  Spvn_for_onebody = &Spvn;
  if(!save_memory && imp_sampl) {

// NEED TO REDEFINE SpvnT if ncores_per_TG > 1  

    SpvnT.setup(head_of_nodes,"SpvnT",TG.getNodeCommLocal());
    int NMO2 = NMO*NMO;
    if(!spinRestricted) NMO2*=2;
    SpvnT.setDims(Spvn.cols(),NMO2);
    if(head_of_nodes) {
      ComplexSMSpMat::const_int_iterator itik = Spvn.rowIndex_begin();
      int cnt1=0;
      for(int i=0; i<NMO; i++)
       for(int k=0; k<NMO; k++, ++itik) {
         if( *itik == *(itik+1) ) continue;
         if(wfn->isOccupAlpha("ImportanceSampling",i)) cnt1+=((*(itik+1))-(*(itik)));
       }
      if(!spinRestricted) {
        for(int i=0; i<NMO; i++)
         for(int k=0; k<NMO; k++, ++itik) {
           if( *itik == *(itik+1) ) continue;
           if(wfn->isOccupBeta("ImportanceSampling",i+NMO)) cnt1+=((*(itik+1))-(*(itik)));
         }
      }
      SpvnT.allocate_serial(cnt1); 
      itik = Spvn.rowIndex_begin();
      const int* cols = Spvn.column_data();
      const ComplexType* vals = Spvn.values(); 
      for(int i=0; i<NMO; i++)
       for(int k=0; k<NMO; k++, ++itik) {
         if( *itik == *(itik+1) ) continue;
         if(wfn->isOccupAlpha("ImportanceSampling",i)) {
           for(int p=(*itik); p<(*(itik+1)); p++) 
             SpvnT.add( *(cols+p) , i*NMO+k , *(vals+p) ); 
         }
       }
      if(!spinRestricted) { 
        for(int i=0; i<NMO; i++)
         for(int k=0; k<NMO; k++, ++itik) {
           if( *itik == *(itik+1) ) continue;
           if(wfn->isOccupBeta("ImportanceSampling",i+NMO)) { 
             for(int p=(*itik); p<(*(itik+1)); p++) 
               SpvnT.add( *(cols+p) , NMO*NMO+i*NMO+k , *(vals+p) ); 
           }
         }
      }
      SpvnT.compress();
    }  
    myComm->barrier();
    if(!head_of_nodes) SpvnT.initializeChildren();
    myComm->barrier();
    Spvn_for_onebody = &SpvnT;
  }

  return true;
}

// right now using local energy form of important sampling
void phaseless_ImpSamp_ForceBias::Propagate(int n, WalkerHandlerBase* wset, RealType& Eshift, const RealType Eav)
{

// split over cholesky vectors Spvn. That way you only need 1 communiation ievent per propagation
// during parallel execution, you'll still need to communicate both vbias and vhs for each walker
// if you are using the hybrid algorithm since you need vbias for hybrid_weight.
//
  
  int nw = wset->numWalkers(true); // true includes dead walkers in list
  MFfactor.resize(nw);
  hybrid_weight.resize(nw);
  if(parallelPropagation) {
  //if(true) {
    // hybrid_weight, MFfactor, and Sdet for each walker is updated 
    dist_Propagate(wset);
    TG.local_barrier();
  } else {
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
        wfn->calculateMixedMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,Sdet,*Spvn_for_onebody,vbias,!save_memory,true);
        apply_bound_vbias();
      }
      Timer.stop("Propagate::calculateMixedMatrixElementOfOneBodyOperators");

      // 3. generate and apply HS propagator (or a good approx of it)
      //  to the S1 and S2 Slater Determinants. The projected determinants are 
      //  returned in S1 and S2. 
      Timer.start("Propagate::applyHSPropagator");
      MFfactor[i]=1.0;
      applyHSPropagator(S1,S2,MFfactor[i],6);  
      Timer.stop("Propagate::applyHSPropagator");

      if(hybrid_method && imp_sampl) {
        ComplexType tmp = ComplexType(0.0,0.0);
        for(int ii=0; ii<sigma.size(); ii++) 
          tmp += (vMF[ii]-vbias[ii])* ( sigma[ii] -0.5*(vMF[ii]-vbias[ii])  ) ;
        hybrid_weight[i] = tmp;  
      }

      std::fill(Sdet,Sdet+2*NMO*NAEA,ComplexType(0,0));
      // propagate forward half a timestep with mean-field propagator
      Timer.start("Propagate::product_SD");
      SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],S1.data(),NAEA,Sdet,NAEA);
      SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],S2.data(),NAEA,Sdet+NAEA*NMO,NAEA);
      Timer.stop("Propagate::product_SD");

    }  // finished propagating walkers
  }

  // in case I implement load balance features in parallelPropagation
  nw = wset->numWalkers(true);

  Timer.start("Propagate::overlaps_and_or_eloc");
  // calculate overlaps and local energy for all walkers
  if(hybrid_method) { 
    if(imp_sampl)
      wfn->evaluateOverlap("ImportanceSampling",n,wset);
  } else {
    wfn->evaluateLocalEnergyAndOverlap("ImportanceSampling",n,wset);
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
      ComplexType factor = ComplexType(1.0,0.0);
      if(imp_sampl) {
        ratioOverlaps = oa*ob/(ooa*oob); 
        factor *= std::exp(hybrid_weight[i]) * ratioOverlaps; 
      }  
      if( (std::isnan(ratioOverlaps.real()) || std::abs(oa*ob) < 1e-8) && apply_constrain && imp_sampl ) { 
        scale = 0.0;
        eloc = oldeloc; 
      } else {  
        scale = apply_constrain?(std::max(0.0,std::cos(std::arg(ratioOverlaps*MFfactor[i])))):(1.0);
        eloc = -std::log( MFfactor[i]*factor)/dt; 
      }
    } else {
      ComplexType ratioOverlaps = MFfactor[i]*oa*ob/(ooa*oob);
      if( (std::isnan(ratioOverlaps.real()) || std::abs(oa*ob) < 1e-8) && apply_constrain ) { 
        scale = 0.0;
        eloc = oldeloc; 
      } else  
        scale = apply_constrain?(std::max(0.0,std::cos(std::arg(ratioOverlaps)))):(1.0);
    }
    eloc = apply_bound_eloc(eloc,Eav);

    w0 *= ComplexType(scale*std::exp( -dt*( 0.5*( eloc.real() + oldeloc.real() ) - Eshift )),0.0);
//app_log()<<myComm->rank() <<" " <<i <<" " <<eloc.real() <<"  " <<oldeloc.real() <<" " <<w0 <<" " <<oa <<" " <<ob <<" " <<ooa <<" " <<oob <<"  " <<scale <<" " <<MFfactor[i] <<std::endl;  
    wset->setWalker(i,w0,eloc);

  }  // loop over walkers

}

void phaseless_ImpSamp_ForceBias::dist_Propagate(WalkerHandlerBase* wset)
{
  // structure in TG [hybrid_w, MFfactor, G(1:{2*}NMO*NMO), sigma(1:nCholVecs), vHS(1:{2*}NMO*NMO)] 
  //  1. You either need to send sigma or communicate back vbias, choosing the former 
  //     In principle, you can replace sigma by vbias on the return communication if you need to accumulate vbias
  
  register int sz = sizeOfG + nCholVecs + NMO*NMO + 2; 
  int nw = wset->numWalkers(true); 
  int nw0 = wset->numWalkers(false); 
  if(!spinRestricted) sz+=NMO*NMO;
  int cnt=0;
  walker_per_node[0] = nw0;
  ComplexType factor;
  int nwglobal = nw0*sz;

  // this routine implies communication inside TG to resize buffers, it also implies a barrier within local TG   
  TG.resize_buffer(nwglobal);
  nwglobal /= sz;
  local_buffer.resize(nCholVecs*nwglobal);

  // add G to buff 
  wfn->evaluateOneBodyMixedDensityMatrix("ImportanceSampling",wset,(TG.commBuff),sz,2,true);
  // add sigma to buff: This does not need to be communicated, unless you want to keep the option open to
  // store a copy later on 
  cnt=0;
  for(int i=0; i<nw; i++) {
    if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
    if(cnt%ncores_per_TG != core_rank) {
      cnt++;
      continue;
    }
    Timer.start("Propagate::sampleGaussianFields");  
    sampleGaussianFields( (TG.commBuff)->values()+cnt*sz+2+sizeOfG, nCholVecs); 
    Timer.stop("Propagate::sampleGaussianFields");  
    // zero out vHS, and hybrid_w, MFfactor 
    *((TG.commBuff)->values()+cnt*sz) = 0;  
    *((TG.commBuff)->values()+cnt*sz+1) = 0;  
    std::fill( (TG.commBuff)->begin()+cnt*sz+2+sizeOfG+nCholVecs, (TG.commBuff)->begin()+(cnt+1)*sz, ComplexType(0,0)); 
    cnt++;
  }
  TG.local_barrier();
  
  // calculate vHS 
  int currnw=nw0; 
  for(int tgi = 0; tgi<nnodes_per_TG; tgi++) {
    // calculates vbias and accumulates vHS
    addvHS( TG.commBuff ,currnw,sz, wset);   

    // rotate date
    if(distributeSpvn) TG.rotate_buffer(currnw,sz);   
  } 

  // synchronize
  TG.local_barrier();

  // propagate walkers now
  cnt=0; 
    for(int i=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue; 
      
      hybrid_weight[i] = *( (TG.commBuff)->begin() + cnt*sz ); 
      MFfactor[i] = *( (TG.commBuff)->begin() + cnt*sz + 1 ); 
      MFfactor[i] = std::exp(-MFfactor[i]); 

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
      // NOTE: In principle, I can define a new applyHSPropagator that uses vHS and vbias from provided pointers
      // instead of using local arrays. I'm doing it this way because I think that it is less efficient to have 
      // the memory being accessed concurrently by all threads in shared memory.
      // But of course this involves an extra std::copy.
      // WHICH ONE IS MORE EFFICIENT??? 
      std::copy((TG.commBuff)->begin()+cnt*sz+2+sizeOfG+nCholVecs,(TG.commBuff)->begin()+(cnt+1)*sz,vHS.begin()); 

      // 3. generate and apply HS propagator (or a good approx of it)
      //  to the S1 and S2 Slater Determinants. The projected determinants are 
      //  returned in S1 and S2. 
      Timer.start("Propagate::applyHSPropagator");
      applyHSPropagator(S1,S2,factor,6,false);  
//debug
//applyHSPropagator(S1,S2,factor,6,true);  
      Timer.stop("Propagate::applyHSPropagator");

      std::fill(Sdet,Sdet+2*NMO*NAEA,ComplexType(0,0));
      // propagate forward half a timestep with mean-field propagator
      Timer.start("Propagate::product_SD");
      SparseMatrixOperators::product_SD(NAEA,Propg_H1.data()+Propg_H1_indx[0],Propg_H1_indx[1],S1.data(),NAEA,Sdet,NAEA);
      SparseMatrixOperators::product_SD(NAEB,Propg_H1.data()+Propg_H1_indx[2],Propg_H1_indx[3],S2.data(),NAEA,Sdet+NAEA*NMO,NAEA);
      Timer.stop("Propagate::product_SD");

      cnt++;
    }  // finished propagating walkers
}

// structure in TG [hybrid_w, MFfactor, G(1:{2*}NMO*NMO), sigma(1:#CholVects), vHS(1:{2*}NMO*NMO)]
void phaseless_ImpSamp_ForceBias::addvHS(ComplexSMVector *buff, int nw, int sz, WalkerHandlerBase* wset)
{

  ComplexType one = ComplexType(1.0,0.0);
  ComplexType minusone = ComplexType(-1.0,0.0);
  ComplexType zero = ComplexType(0.0,0.0);
  int nt = NMO*NMO;
  if(!spinRestricted) nt *= 2; 

  for(int i=0; i<CV0.size(); i++) CV0[i]=zero; 
  if(TG.getCoreRank()==0) 
    std::fill(local_buffer.begin(),local_buffer.end(),zero);

  TG.local_barrier();

  // Calculate vbias for all walkers.
  ComplexType *pos = buff->values();
  if(imp_sampl) {
    for(int wlk=0; wlk<nw; wlk++, pos+=sz) {

      // calculate vbias for local sector in 'ik' space. This routine assumes a buffer created by wfn->evaluateOneBodyMixedDensityMatrix 
      Timer.start("Propagate::calculateMixedMatrixElementOfOneBodyOperators");

      wfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(spinRestricted,"ImportanceSampling",-1,pos+2,ik0,ikN,pik0,*Spvn_for_onebody,vbias,!save_memory,false);    

      Timer.stop("Propagate::calculateMixedMatrixElementOfOneBodyOperators");
      Timer.start("Propagate::addvHS::shm_copy");
      {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
        zaxpy(cvecN-cvec0,one,vbias.data()+cvec0,1,local_buffer.values()+wlk*nCholVecs+cvec0,1);  
      } 
      Timer.stop("Propagate::addvHS::shm_copy");
    }
  }

  TG.local_barrier();

//debug debug debug
/*
for(int wlk=0; wlk<nw; wlk++) {
app_log()<<" GF (diag): \n";
for(int ii=0; ii<NAEA; ii++)
  app_log()<<ii <<" " <<*( buff->begin() + wlk*sz + 4 + ii*NMO + ii) <<std::endl;
 std::vector<ComplexType> vbias_(vbias.size());
 wfn->calculateMixedMatrixElementOfOneBodyOperators(spinRestricted,"ImportanceSampling",-1,wset->getSM(wlk),*Spvn_for_onebody,vbias_,!save_memory,true);
 RealType dif = 0;
 for(int ii=0; ii<vbias.size(); ii++)
   dif += std::abs(vbias_[ii]- *(local_buffer.values()+wlk*nCholVecs+cvec0+ii) );
 std::cout<<rank() <<" " <<wlk <<" " <<dif <<std::endl; 
 for(int ii=0; ii<vbias.size(); ii++)
   app_log()<<ii <<" " <<vbias_[ii] <<" " <<*(local_buffer.values()+wlk*nCholVecs+cvec0+ii) <<" " <<abs(vbias_[ii]- *(local_buffer.values()+wlk*nCholVecs+cvec0+ii)) <<std::endl;
 app_log()<<std::endl;
}
*/
//for(int ii=0; ii<vbias.size(); ii++)
//  vbias[ii] = vbias_[ii];

//for(int ii=0; ii<sigma.size(); ii++)
//  sigma[ii] = (buff->values()+2+sizeOfG+cvec0+ii)->real(); 
  
  pos = buff->values();
  for(int wlk=0; wlk<nw; wlk++, pos+=sz) {

    // accumulate vHS
    //for(ComplexMatrix::iterator it=vHS.begin(); it!=vHS.end(); it++) *it=zero;  // no need to zero out

    Timer.start("Propagate::build_vHS");
    ComplexType* sg = pos+2+sizeOfG+cvec0;
    ComplexType* vb = local_buffer.values()+wlk*nCholVecs+cvec0;
    ComplexType mf=zero, hw=zero; 
    for(int i=cvec0; i<cvecN; i++, sg++,vb++) { 
      ComplexType vbi = apply_bound_vbias(*vb);
      CV0[i] = *(sg) + ( vbi - vMF[i] );
      mf += CV0[i]*vMF[i]; 
      hw += (vMF[i]-vbi)*( *(sg) -0.5*(vMF[i]-vbi) );
    }
    // vHS
    SparseMatrixOperators::product_SpMatV(int(ikN-ik0), Spvn.cols(), one, Spvn.values() + pik0, Spvn.column_data() + pik0, Spvn.row_index()+ik0, CV0.data(), zero, vHS.data()+ik0);
    Timer.stop("Propagate::build_vHS");
    Timer.start("Propagate::addvHS::shm_copy");

    // These are rigurously non-overlapping memory regions, but they are adjacent.
    // Do I still need to control concurrent access??? Probably NOT!!!
    if(TG.getCoreRank()==0)
    {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
      // hybrid_weight
      *pos += hw;
      // MFfactor
      *(pos+1) += mf;
      zaxpy(int(ikN-ik0),one,vHS.data()+ik0,1,pos+2+sizeOfG+nCholVecs+ik0,1);  
    } else {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buff->getMutex()));
      zaxpy(int(ikN-ik0),one,vHS.data()+ik0,1,pos+2+sizeOfG+nCholVecs+ik0,1);  
    } 
    Timer.stop("Propagate::addvHS::shm_copy");

  } 

  TG.local_barrier();

}

void phaseless_ImpSamp_ForceBias::applyHSPropagator(ComplexMatrix& M1, ComplexMatrix& M2, ComplexType& factor,  int order, bool calculatevHS)
{

  if(order < 0) order = Order_Taylor_Expansion;

  ComplexType one = ComplexType(1.0,0.0); 
  ComplexType minusone = ComplexType(-1.0,0.0); 
  ComplexType zero = ComplexType(0.0,0.0); 


  //if(calculatevHS) {
  if(calculatevHS) {
    factor = zero;
//    for(ComplexMatrix::iterator it=vHS.begin(); it!=vHS.end(); it++) *it=zero;

    Timer.start("Propagate::build_vHS");
    for(int i=0; i<sigma.size(); i++) {
      CV0[i] = sigma[i] + (vbias[i]-vMF[i]);
      factor += CV0[i]*vMF[i];
    }
    SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn,CV0.data(),zero,vHS.data());
    Timer.stop("Propagate::build_vHS");
    factor = std::exp(-factor);
  }

  // calculate exp(vHS)*S through a Taylor expansion of exp(vHS)
  Timer.start("Propagate::apply_expvHS_Ohmms");
  T1=M1;
  for(int n=1; n<=order; n++) {
    ComplexType fact = static_cast<ComplexType>(1.0/static_cast<double>(n)); 
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
  int disp = (spinRestricted)?0:NMO*NMO;
  T1 = M2;
  for(int n=1; n<=order; n++) {
    ComplexType fact = static_cast<ComplexType>(1.0/static_cast<double>(n));
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
    RealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    RealType mag = std::sqrt(-2.0*std::log(temp1));
    sigma[i]  =mag*std::cos(6.283185306*temp2);
    sigma[i+1]=mag*std::sin(6.283185306*temp2);
  }
  if (n%2==1)
  {
    RealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    sigma[n-1]=std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2);
  }

}

void phaseless_ImpSamp_ForceBias::sampleGaussianFields(ComplexType* v, int n)
{ 
  for (int i=0; i+1<n; i+=2)
  { 
    RealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    RealType mag = std::sqrt(-2.0*std::log(temp1));
    *(v++)  = ComplexType(mag*std::cos(6.283185306*temp2),0);
    *(v++)  = ComplexType(mag*std::sin(6.283185306*temp2),0);
  }
  if (n%2==1)
  { 
    RealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    *v = ComplexType(std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2),0);
  }

}

// this must be called after the setup has been made 
void phaseless_ImpSamp_ForceBias::test_linear_algebra()
{

  if(myComm->rank() != 0) return;
  if(!spinRestricted) return;

  // sparse versus dense Spvn matrices

  ComplexType one = ComplexType(1.0,0.0);
  ComplexType minusone = ComplexType(-1.0,0.0);
  ComplexType zero = ComplexType(0.0,0.0);

  RealSMSpMat SpvnReal;
  SpvnReal.setup(true,"Test_Test_Test",MPI_COMM_WORLD);
  SpvnReal.setDims(Spvn.rows(),Spvn.cols());
  SpvnReal.allocate_serial(Spvn.size());
  SpvnReal.resize_serial(Spvn.size());
  
  std::copy( Spvn.cols_begin(), Spvn.cols_end(), SpvnReal.cols_begin());
  std::copy( Spvn.rows_begin(), Spvn.rows_end(), SpvnReal.rows_begin());
  ComplexSMSpMat::iterator itv=Spvn.vals_begin();
  RealSMSpMat::iterator ritv=SpvnReal.vals_begin();
  for(int i=0; i<Spvn.size(); i++)
  {
    *(ritv++) = (itv++)->imag(); 
  } 
  SpvnReal.compress();

  SMSparseMatrix<float> SpvnSP;
  SpvnSP.setup(true,"Test_Test_Test_2",MPI_COMM_WORLD);
  SpvnSP.setDims(Spvn.rows(),Spvn.cols());
  SpvnSP.allocate_serial(Spvn.size());
  SpvnSP.resize_serial(Spvn.size());

  std::copy( Spvn.cols_begin(), Spvn.cols_end(), SpvnSP.cols_begin());
  std::copy( Spvn.rows_begin(), Spvn.rows_end(), SpvnSP.rows_begin());
  itv=Spvn.vals_begin();
  SMSparseMatrix<float>::iterator spitv=SpvnSP.vals_begin();
  for(int i=0; i<Spvn.size(); i++)
  {
    *(spitv++) = static_cast<float>((itv++)->imag());
  }
  SpvnSP.compress();

  for(ComplexMatrix::iterator it=vHS.begin(); it!=vHS.end(); it++) *it=zero;
  sampleGaussianFields();
  for(int i=0; i<sigma.size(); i++) 
    CV0[i] = sigma[i];

  // for possible setup/delays 
  for(int i=0; i<10; i++)
    SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn,CV0.data(),zero,vHS.data());
  int ntimes = 20;
  Timer.reset("Propagate::Sparse::build_vHS");
  Timer.start("Propagate::Sparse::build_vHS");
  for(int i=0; i<ntimes; i++)
    SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn,CV0.data(),zero,vHS.data());
  Timer.stop("Propagate::Sparse::build_vHS"); 

  app_log()<<"Sparsity value: " <<static_cast<double>(Spvn.size())/NMO/NMO/Spvn.cols() <<std::endl;
  app_log()<<"Time for std::complex SpMV product:" <<Timer.total("Propagate::Sparse::build_vHS")/ntimes <<std::endl;

  char trans1 = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  double oned = 1, zerod=0;

  ComplexMatrix vHS2(NMO,NMO);
  for(int i=0; i<10; i++)
    SparseMatrixOperators::product_SpMatM(SpvnReal.rows(),2,SpvnReal.cols(),1,SpvnReal,reinterpret_cast<double*>(CV0.data()),2,0,reinterpret_cast<double*>(vHS2.data()),2);
  Timer.reset("Propagate::Sparse::build_vHS");
  Timer.start("Propagate::Sparse::build_vHS");
  for(int i=0; i<ntimes; i++)
    SparseMatrixOperators::product_SpMatM(SpvnReal.rows(),2,SpvnReal.cols(),1,SpvnReal,reinterpret_cast<double*>(CV0.data()),2,0,reinterpret_cast<double*>(vHS2.data()),2);
  Timer.stop("Propagate::Sparse::build_vHS");

  app_log()<<"Time for real SpMM product:" <<Timer.total("Propagate::Sparse::build_vHS")/ntimes <<std::endl;

  std::cout<<" Looking for difference in vHS: \n";   
  for(int i=0; i<NMO; i++)
   for(int j=0; j<NMO; j++)
    if(std::abs(vHS(i,j).imag()-vHS2(i,j)) > 1e-8)
     std::cout<<i <<" " <<j <<" " <<vHS(i,j).imag() <<" " <<vHS2(i,j) <<std::endl;
  std::cout<<std::endl;

  Matrix<float> vHS3(NMO,NMO);
  std::vector<std::complex<float> > CV1(CV0.size());
  for(int i=0; i<CV0.size(); i++) CV1[i] = static_cast<std::complex<float> >(CV0[i]);
  for(int i=0; i<10; i++)
    SparseMatrixOperators::product_SpMatM(SpvnSP.rows(),2,SpvnSP.cols(),1,SpvnSP,reinterpret_cast<float*>(CV1.data()),2,0,reinterpret_cast<float*>(vHS3.data()),2);
  Timer.reset("Propagate::Sparse::build_vHS");
  Timer.start("Propagate::Sparse::build_vHS");
  for(int i=0; i<ntimes; i++)
    SparseMatrixOperators::product_SpMatM(SpvnSP.rows(),2,SpvnSP.cols(),1,SpvnSP,reinterpret_cast<float*>(CV1.data()),2,0,reinterpret_cast<float*>(vHS3.data()),2);
  Timer.stop("Propagate::Sparse::build_vHS");

  app_log()<<"Time for SP real SpMM product:" <<Timer.total("Propagate::Sparse::build_vHS")/ntimes <<std::endl;

  const char trans = 'T';
  int nrows = NMO*NMO, ncols=Spvn.cols();
  int inc=1;
  ComplexMatrix Dvn(nrows,ncols);
  // stupid but easy
  for(int i=0; i<Spvn.size(); i++)
    Dvn(*(Spvn.row_data()+i),*(Spvn.column_data()+i)) = *(Spvn.values()+i); 

  // setup time
  for(int i=0; i<10; i++)
    zgemv(trans,Dvn.cols(),Dvn.rows(),one,Dvn.data(),Dvn.cols(),CV0.data(),inc,zero,vHS2.data(),inc);

  Timer.reset("Propagate::Dense::build_vHS");
  Timer.start("Propagate::Dense::build_vHS");
  for(int i=0; i<ntimes; i++)
    zgemv(trans,Dvn.cols(),Dvn.rows(),one,Dvn.data(),Dvn.cols(),CV0.data(),inc,zero,vHS2.data(),inc);
  Timer.stop("Propagate::Dense::build_vHS");

  app_log()<<"Time for dense zgemv product:" <<Timer.total("Propagate::Dense::build_vHS")/ntimes <<std::endl;

  RealMatrix Rvn(nrows,ncols);
  Matrix<float> Fvn(nrows,ncols);
  for(int i=0; i<Spvn.size(); i++)
    Rvn(*(Spvn.row_data()+i),*(Spvn.column_data()+i)) = (Spvn.values()+i)->imag();
  for(int i=0; i<Spvn.size(); i++)
    Fvn(*(Spvn.row_data()+i),*(Spvn.column_data()+i)) = static_cast<float>((Spvn.values()+i)->imag());

  // setup time: Rvn*CV0 = vHS
  for(int i=0; i<10; i++)
    dgemm('N','N', 2, Rvn.rows(), Rvn.cols(), 1, reinterpret_cast<double*>(CV0.data()), 2, Rvn.data(), Rvn.cols(), 0, reinterpret_cast<double*>(vHS.data()), 2); 
  Timer.reset("Propagate::Dense::build_vHS");
  Timer.start("Propagate::Dense::build_vHS");
  for(int i=0; i<ntimes; i++)
    dgemm('N','N', 2, Rvn.rows(), Rvn.cols(), 1, reinterpret_cast<double*>(CV0.data()), 2, Rvn.data(), Rvn.cols(), 0, reinterpret_cast<double*>(vHS.data()), 2); 
  Timer.stop("Propagate::Dense::build_vHS");

  app_log()<<"Time for dense dgemm product:" <<Timer.total("Propagate::Dense::build_vHS")/ntimes <<std::endl;

  for(int i=0; i<10; i++)
    sgemm('N','N', 2, Fvn.rows(), Fvn.cols(), 1, reinterpret_cast<float*>(CV1.data()), 2, Fvn.data(), Fvn.cols(), 0, reinterpret_cast<float*>(vHS3.data()), 2);
  Timer.reset("Propagate::Dense::build_vHS");
  Timer.start("Propagate::Dense::build_vHS");
  for(int i=0; i<ntimes; i++)
    sgemm('N','N', 2, Fvn.rows(), Fvn.cols(), 1, reinterpret_cast<float*>(CV1.data()), 2, Fvn.data(), Fvn.cols(), 0, reinterpret_cast<float*>(vHS3.data()), 2);
  Timer.stop("Propagate::Dense::build_vHS");

  app_log()<<"Time for SP dense dgemm product:" <<Timer.total("Propagate::Dense::build_vHS")/ntimes <<std::endl;


  APP_ABORT("TESTING TESTING TESTING. \n\n\n");

}

}
