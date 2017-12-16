
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
#include "AFQMC/Walkers/SlaterDetWalker.h"
#include "AFQMC/Propagators/VMCPropagator.h"
#include "AFQMC/Hamiltonians/ProjectorBase.h"
#include "AFQMC/Hamiltonians/DDProjector.h"
#include "AFQMC/Hamiltonians/CCProjector.h"

#include"AFQMC/Numerics/SparseMatrixOperations.h"
#include"AFQMC/Numerics/DenseMatrixOperations.h"
#include "Utilities/RandomGenerator.h"

namespace qmcplusplus
{

bool VMCPropagator::parse(xmlNodePtr cur)
{

    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    ParameterSet m_param;
    // hdf
    m_param.add(hdf_read_tag,"hdf_read_tag","std::string");
    m_param.add(hdf_read_file,"hdf_read_file","std::string");
    m_param.add(hdf_write_tag,"hdf_write_tag","std::string");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");

    m_param.put(cur);

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="projector") {
        std::string type("DD");
        xmlNodePtr curRoot=cur;
        OhmmsAttributeSet oAttrib;
        oAttrib.add(type,"type");
        oAttrib.put(cur);
        if(type=="DD" || type=="Gutzwiller" || type=="dd" || type == "gutzwiller" || type=="DDProjector" )
        { 
          proj0 = new DDProjector(myComm);
          proj0->parse(cur);
        } else if(type=="CC" || type=="cc" || type=="cluster" || type=="Cluster" )
        {
          proj0 = new CCProjector(myComm);
          proj0->parse(cur);
        } else {
          std::cerr<<"Unknown projector type: " <<type <<std::endl;
        }
      }
      cur = cur->next;
    }

    return true;

}

bool VMCPropagator::hdf_write(hdf_archive& dump,const std::string& tag)
{

  std::string path = "/Propagators/VMCPropagator"; 
  if(tag != std::string("")) path += std::string("/")+tag; 
  if(dump.is_group( path )) {
    app_error()<<" ERROR: H5Group /Propagators/VMCPropagator/{tag} already exists in restart file. This is a bug and should not happen. Contact a developer.\n";
    return false;
  }

  std::vector<int> Idata(4);
  Idata[0]=Spvn.size();
  Idata[1]=Spvn.rows();
  Idata[2]=Spvn.cols();
  Idata[3]=NMO;

  dump.push("Propagators");
  dump.push("VMCPropagator");
  if(tag != std::string("")) dump.push(tag);
  dump.write(Idata,"Spvn_dims");
  dump.write(*(Spvn.getVals()),"Spvn_vals");
  dump.write(*(Spvn.getRows()),"Spvn_rows");
  dump.write(*(Spvn.getCols()),"Spvn_cols");
  dump.write(*(Spvn.getRowIndex()),"Spvn_rowIndex");
  if(tag != std::string("")) dump.pop();
  dump.pop();
  dump.pop();

  dump.flush();  

  return true;
}

bool VMCPropagator::hdf_read(hdf_archive& dump,const std::string& tag)
{
  std::vector<int> Idata(4);

  if(!dump.push("Propagators",false)) return false;
  if(!dump.push("VMCPropagator",false)) return false;
  if(tag != std::string("")) 
    if(!dump.push(tag,false)) return false;
  if(!dump.read(Idata,"Spvn_dims")) return false;

  if(Idata[3] != NMO) {
    app_error()<<" Error in VMCPropagator::hdf_read. NMO is not consistent between hdf5 file and current run. \n";
    return false; 
  }  
  Spvn.setDims(Idata[1],Idata[2]);
  Spvn.allocate_serial(Idata[0]);
  Spvn.resize_serial(Idata[0]);

  if(!dump.read(*(Spvn.getVals()),"Spvn_vals")) return false;
  if(!dump.read(*(Spvn.getRows()),"Spvn_rows")) return false;
  if(!dump.read(*(Spvn.getCols()),"Spvn_cols")) return false;
  if(!dump.read(*(Spvn.getRowIndex()),"Spvn_rowIndex")) return false;
  if(tag != std::string("")) dump.pop();
  dump.pop();
  dump.pop();

  // check that everything went fine

  return true;
}



bool VMCPropagator::setup(std::vector<int>& TGdata, SPComplexSMVector *v, HamiltonianBase* ham,WavefunctionHandler* w, RealType dt_, hdf_archive& dump_read, const std::string& hdf_restart_tag,MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm node_heads_comm)
{

  MPI_COMM_HEAD_OF_NODES = node_heads_comm;
  head_of_nodes = (TGdata[1]==0);

  dt = dt_;

  proj0->copyInfo(*this);

  bool read_Spvn_from_file=false;
  Spvn.setup(head_of_nodes,"Spvn",MPI_COMM_WORLD); // THIS IS WRONG FIX FIX FIX

  proj0->init(ham);

  // Only master tries to read 
  if(myComm->rank() == 0) {

    if(hdf_read_file!=std::string("")) {

      hdf_archive dump(myComm);
      if(dump.open(hdf_read_file,H5F_ACC_RDONLY)) { 
        read_Spvn_from_file = hdf_read(dump,hdf_read_tag);
        dump.close();
        if(read_Spvn_from_file) 
          app_log()<<"Successfully read HS potentials from file: " <<hdf_read_file <<"\n";
      }

    } else {
      if(dump_read.file_id != hdf_archive::is_closed) { 
        read_Spvn_from_file = hdf_read(dump_read,hdf_restart_tag); 

        if(read_Spvn_from_file) 
          app_log()<<"Successfully read HS potentials from restart file. \n";
        else 
          while(dump_read.top() != hdf_archive::is_closed ) { dump_read.pop(); } 
      }
    }

  }

  int success = read_Spvn_from_file?0:1; 
  myComm->bcast<int>(&success,1); 

  if(read_Spvn_from_file || !parallel_factorization) {

    if(rank()==0) {

      if(!read_Spvn_from_file) {

        app_log()<<" Calculating HS potentials from scratch. \n";

        Timer.reset("Generic1");
        Timer.start("Generic1");

        // calculates Hubbard-Stratonovich potentials (vn) 
        //proj0->calculateHSPotentials(Spvn);
	
        Timer.stop("Generic1");
        app_log()<<" -- Time to calculate HS potentials: " <<Timer.average("Generic1") <<"\n";
      }

      std::vector<int> ni(3);
      ni[0]=Spvn.size();
      ni[1]=Spvn.rows();
      ni[2]=Spvn.cols();
      myComm->bcast(ni);        

      // do this later through a proper MPI object
      myComm->bcast<char>(reinterpret_cast<char*>(Spvn.values()),sizeof(ComplexType)*Spvn.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<int>(Spvn.row_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<int>(Spvn.column_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
      myComm->bcast<int>(Spvn.row_index(),Spvn.rows()+1,MPI_COMM_HEAD_OF_NODES);

      myComm->barrier();
      myComm->barrier();

    } else {

      std::vector<int> ni(3);
      myComm->bcast(ni);
      Spvn.setDims(ni[1],ni[2]);
      if(head_of_nodes) {

        Spvn.allocate_serial(ni[0]);
        Spvn.resize_serial(ni[0]);

        myComm->bcast<char>(reinterpret_cast<char*>(Spvn.values()),sizeof(ComplexType)*Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.row_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.column_data(),Spvn.size(),MPI_COMM_HEAD_OF_NODES);
        myComm->bcast<int>(Spvn.row_index(),Spvn.rows()+1,MPI_COMM_HEAD_OF_NODES);

      }
      Spvn.setCompressed();

      myComm->barrier();
      if(!head_of_nodes) Spvn.initializeChildren();
      myComm->barrier();
 
    }

  } else {
    // calculating Spvn in parallel
    APP_ABORT("Parallel calculation of Spvn not implemented\n");
  } 

  // write restart if desired
  if(rank()==0) {
    if(hdf_write_file!=std::string("")) {
      hdf_archive dump(myComm);
      if(dump.create(hdf_write_file)) {
        if(!hdf_write(dump,hdf_write_tag)) {
          app_error()<<" Problems writing hdf5 file in VMCPropagator::setup(). \n";
          return false;
        }
        dump.close();
      } else {
        app_error()<<" Problems opening hdf5 file in VMCPropagator::setup(). \n";
        return false;
      }
    } 
  }

  app_log()<<" -- Total number of terms in Cholesky vectors: " <<Spvn.size() <<std::endl; 

  // setup matrices
  vHS.resize(2*NMO,NMO);
  sigmaL.resize(Spvn.cols());
  sigmaR.resize(Spvn.cols());
  for(int i=0; i<sigmaR.size(); i++) sigmaR[i]=0;
  for(int i=0; i<sigmaL.size(); i++) sigmaL[i]=0;

  // setup temporary storage
  T1.resize(NMO,NAEA);
  T2.resize(NMO,NAEA);

  SL.resize(2*NMO,NAEA);
  SR.resize(2*NMO,NAEA);

  return true;
}

// right now using local energy form of important sampling
void VMCPropagator::Propagate(int steps, int& steps_total, WalkerHandlerBase*, RealType& E1)
{
/*
  int sz = NMO*NAEA;
  std::copy(w.SlaterMat.begin(),w.SlaterMat.begin()+2*sz,SL.begin());
  std::copy(w.SlaterMat.begin()+2*sz,w.SlaterMat.begin()+4*sz,SR.begin());

  // 1. sample gaussian field
  Timer.start("Propagate::sampleGaussianFields");  
  sampleGaussianFields(sigmaL); 
  sampleGaussianFields(sigmaR); 
  Timer.stop("Propagate::sampleGaussianFields");  
  
  // 3. generate and apply HS propagator (or a good approx of it)
  //  to the S1 and S2 Slater Determinants. The projected determinants are 
  //  returned in S1 and S2. 
  Timer.start("Propagate::applyHSPropagator");
  applyHSPropagator(SL.data(),sigmaL,6);  
  applyHSPropagator(SR.data(),sigmaR,6);  
  Timer.stop("Propagate::applyHSPropagator");

  Timer.start("Propagate::ovlp");
  ComplexType wgt;
  SDetOps->green_function(SL,SR,wgt,T1,false);
  Timer.stop("Propagate::ovlp");

  if( (*rng)() < std::abs(wgt)/std::abs(w.weight) ) {
    accept++;
    std::copy(SL.begin(),SL.end(),w.SlaterMat.begin());
    std::copy(SR.begin(),SR.end(),w.SlaterMat.begin()+2*sz);
    w.weight = wgt; 
  } 
*/
}

void VMCPropagator::applyHSPropagator(ComplexType* SD, std::vector<ComplexType>& sigma,  int order)
{

  if(order < 0) order = Order_Taylor_Expansion;

  ComplexType one = ComplexType(1.0,0.0); 
  ComplexType minusone = ComplexType(-1.0,0.0); 
  ComplexType zero = ComplexType(0.0,0.0); 

  for(ComplexMatrix::iterator it=vHS.begin(); it!=vHS.end(); it++) *it=zero;

  Timer.start("build_vHS");
  SparseMatrixOperators::product_SpMatV(Spvn.rows(),Spvn.cols(),one,Spvn.values(),Spvn.column_data(),Spvn.row_index(),sigma.data(),zero,vHS.data());
  Timer.stop("build_vHS");

  // calculate exp(vHS)*S through a Taylor expansion of exp(vHS)
  Timer.start("apply_expvHS_Ohmms");
  int sz=NMO*NAEA;
  std::copy(SD,SD+sz,T1.begin());
  for(int n=1; n<=order; n++) {
    ComplexType fact = static_cast<ComplexType>(1.0/static_cast<double>(n)); 
    DenseMatrixOperators::product(NMO,NAEA,NMO,fact,vHS.data(),NMO,T1.data(),NAEA,zero,T2.data(),NAEA);
    T1  = T2;
    ComplexType* itSD = SD;
    for(ComplexMatrix::iterator it=T1.begin(); it!=T1.end(); it++, itSD++)
      *itSD += *it;
  }

//  if M1 == M2 on entry, no need to do this in spinRestricted case
  std::copy(SD+sz,SD+2*sz,T1.begin());
  for(int n=1; n<=order; n++) {
    ComplexType fact = static_cast<ComplexType>(1.0/static_cast<double>(n));
    DenseMatrixOperators::product(NMO,NAEB,NMO,fact,vHS.data()+NMO*NMO,NMO,T1.data(),NAEA,zero,T2.data(),NAEA);
    T1  = T2;
    ComplexType* itSD = SD+sz;
    for(ComplexMatrix::iterator it=T1.begin(); it!=T1.end(); it++, itSD++)
      *itSD += *it;
  } 
  Timer.stop("apply_expvHS_Ohmms");
}

void VMCPropagator::sampleGaussianFields(std::vector<ComplexType>& sigma)
{
  int n = sigma.size();
  for (int i=0; i+1<n; i+=2)
  {
    RealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    sigma[i]  =ComplexType(std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2),0.0);
    sigma[i+1]=ComplexType(std::sqrt(-2.0*std::log(temp1))*std::sin(6.283185306*temp2),0.0);
  }
  if (n%2==1)
  {
    RealType temp1=1-0.9999999999*(*rng)(), temp2=(*rng)();
    sigma[n-1]=ComplexType(std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2),0.0);
  }

}

}
