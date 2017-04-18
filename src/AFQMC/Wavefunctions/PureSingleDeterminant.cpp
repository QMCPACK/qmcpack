#include<cstdlib>
#include<complex>
#include<iostream>
#include<fstream>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/SimpleParser.h"
#include "Message/CommOperators.h"
#include "Configuration.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Hamiltonians/SparseGeneralHamiltonian.h"
#include "AFQMC/Wavefunctions/WavefunctionBase.h"
#include "AFQMC/Wavefunctions/PureSingleDeterminant.h"
#include "AFQMC/Utilities/Utils.h"

#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

//#include "mkl.h"
//#include "mkl_service.h"

namespace qmcplusplus
{

bool PureSingleDeterminant::parse(xmlNodePtr cur)
{
    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur;
   
    std::string type("");
    init_type = std::string("");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.add(type,"type");
    oAttrib.add(init_type,"init");
    oAttrib.put(cur);

    std::transform(type.begin(),type.end(),type.begin(),(int (*)(int)) tolower);
    std::transform(init_type.begin(),init_type.end(),init_type.begin(),(int (*)(int)) tolower);

    if(type != "puresd") {
      app_error()<<" ERROR: Problems in PureSingleDeterminant::parse: type should be PureSD. \n"  <<std::endl; 
      return false;
    }

    std::string str("false");
    filename = std::string("none");
    filetype = std::string("ascii");
    ParameterSet m_param;
    m_param.add(filename,"filename","std::string");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(nnodes_per_TG,"nnodes_per_TG","int");
    m_param.add(nnodes_per_TG,"nodes","int");
    m_param.add(nnodes_per_TG,"nnodes","int");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");
    m_param.add(cutoff,"cutoff","double");
    m_param.put(cur);

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }

    return true;
}

bool PureSingleDeterminant::initFromAscii(std::string fileName)
{

  occup_alpha.clear();
  occup_beta.clear();

  int rnk=0;
#if defined(USE_MPI)
  rnk = rank();
#endif

  std::ifstream in;
  in.open(fileName.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
     return false;
  }

  // this class should know nothing about core states
  NCA = NCB = 0;
  bool fullMOMat = false;
  bool Cstyle = true;

  // Read header, but do not overwrite variables that are >= 0 (initialized in xml input). 
  std::vector<std::string> words;
  getwords(words,in);       
  do {
    if(words.size() == 0) 
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
    for(std::vector<std::string>::iterator it=words.begin(); it!=words.end(); it++) {
      if(*it == "&FCI") {
        // do nothing   
      } else if(*it == "NORB") {
        if( it+1 == words.end() ) {
          app_error()<<"Format error in ASCII integral file. NORB \n";
          return false;
        }
        if(NMO < 0) NMO = atoi((++it)->c_str());   
        else it++;
      } else if(*it == "NAEA") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NAEA \n";
          return false;
        }
        if(NAEA < 0) NAEA = atoi((++it)->c_str());   
        else it++;
      } else if(*it == "NAEB") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NAEB \n";
          return false;
        }
        if(NAEB < 0) NAEB = atoi((++it)->c_str());   
        else it++;
      } else if(*it == "NCA" || *it == "NCB") {
        it++;
      } else if(*it == "NELEC") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NETOT \n";
          return false;
        }
        if(NETOT < 0) NETOT = atoi((++it)->c_str());
        else it++;
      } else if(*it == "Type" || *it == "TYPE" || *it == "type") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NAEB \n";
          return false;
        }
        if( *(it+1) == "rotated" || *(it+1) == "rotate" ) rotated_hamiltonian=true;
        else rotated_hamiltonian=false;
        it++;
      } else if(*it == "UHF" || *it == "GHF") {
        if( it+1 == words.end() ) {
          app_error()<<"Format error in ASCII integral file. UHF/GHF \n";
          return false;
        }
        wfntype = atoi((++it)->c_str());
        switch(wfntype) {
          case 0:
          {
            app_log()<<"Using a RHF-type trial wave-function in lSlaterDeterminant. \n";
            break;
          }
          case 1:
          {
            app_log()<<"Using a UHF-type trial wave-function in PureSingleDeterminant. \n";
            break;
          }
          case 2:
          {
            app_log()<<"Using a GHF-type trial wave-function in PureSingleDeterminant. \n";
            app_error()<<" GHF type not implemented. \n";
            return false;
            break;
          }
          default:
          {
            app_error()<<"Unknown wave-function type in PureSingleDeterminant: " <<wfntype <<std::endl;
            return false;
          }
        }
      } else if(*it == "OCCUP_ALPHA") {  
        if( NAEA < 0 ) {
          app_error()<<" NAEA  must be defined before OCCUP_ALPHA in ASCII integral file.\n";
          return false;
        }
        if( it+NAEA == words.end() ) {
          app_error()<<"Format error in ASCII integral file. OCCUP_ALPHA \n";
          return false;
        }
        occup_alpha.resize(NAEA);
        it++;
        for(int i=0; i<NAEA; i++,it++) occup_alpha[i] = atoi(it->c_str())-1;
        std::sort(occup_alpha.begin(),occup_alpha.end());
      } else if(*it == "OCCUP_BETA") {
        if( NAEB < 0 ) {
          app_error()<<"NAEB  must be defined before OCCUP_ALPHA in ASCII integral file.\n";
          return false;
        }
        if( it+NAEB == words.end() ) {
          app_error()<<"Format error in ASCII integral file. OCCUP_BETA \n";
          return false;
        }
        occup_beta.resize(NAEB);
        it++;
        for(int i=0; i<NAEB; i++,it++) occup_beta[i] = atoi(it->c_str())-1+NMO;
        std::sort(occup_beta.begin(),occup_beta.end());
      } else if(*it == "OCCUP") {
        if( NAEB < 0 || NAEA < 0 || NAEA != NAEB || NCA!=NCB ) {
          app_error()<<"OCCUP std::string in ASCII integral file requires NCA=NCB,NAEA=NAEB,NAEA>0,NAEB>0. \n" << std::endl;
          return false;
        }
        if( words.size() < NAEA+1 ) {
          app_error()<<"Format error in ASCII integral file. OCCUP \n" <<std::endl;
          return false;
        }
        occup_alpha.resize(NAEA);
        occup_beta.resize(NAEB);
        for(int i=0; i<NAEA; i++) occup_alpha[i] = atoi((++it)->c_str())-1;
        for(int i=0; i<NAEB; i++) occup_beta[i] = occup_alpha[i]+NMO;
        std::sort(occup_beta.begin(),occup_beta.end());
        std::sort(occup_alpha.begin(),occup_alpha.end());
      } else if(*it == "FullMO" || *it == "FULLMO") {
        fullMOMat = true;
        app_log()<<" Expecting full MO matrix in PureSingleDeterminant.\n";
      } else if(*it == "CMajor") {
        Cstyle = false;
        app_log()<<" Expecting MO matrix in Column-major format in PureSingleDeterminant.\n";
      } else {
//        app_log()<<"Ignoring unknown tag in ASCII integral file: " <<*it <<std::endl;
      }
    }
    getwords(words,in);
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
  } while((words[0].find(std::string("/"))==std::string::npos && words[0].find(std::string("&END"))==std::string::npos));

  if(rotated_hamiltonian) {

    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);
    for(int i=0; i<NAEA; i++) occup_alpha[i] = i;
    for(int i=0; i<NAEB; i++) occup_beta[i] = i+NMO;

    int ncols = NAEA;
    //int nrows = NMO;
    if(wfntype == 2)
      ncols = NAEA+NAEB;

    OrbMat.resize(2*NMO,ncols);
    ComplexType dummy;
    int nread;

    if (wfntype == 0 ) {

       if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
       } else {
        nread = fullMOMat?NMO:NAEA;
        for(int j=0; j<nread; j++)
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
       }

    } else if(wfntype == 1) {

    
      if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
        nread = fullMOMat?NMO:NAEB;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
      } else {
        nread = fullMOMat?NMO:NAEA;
        for(int j=0; j<nread; j++) 
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
        nread = fullMOMat?NMO:NAEB;
        for(int j=0; j<nread; j++) 
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
      }

    } else if(wfntype == 2) {

      if(Cstyle) {

       nread = fullMOMat?NMO:NAEA;
       int nread2 = fullMOMat?NMO:NAEB;
       for(int i=0; i<NMO; i++) {
        for(int j=0; j<nread; j++) {
          in>>dummy;
          if(j<NAEA) OrbMat(i,j) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
        for(int j=0; j<nread2; j++) {
          in>>dummy;
          if(j<NAEB) OrbMat(i,j+NAEA) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
       }
       for(int i=0; i<NMO; i++) {
        for(int j=0; j<nread; j++) {
          in>>dummy;
          if(j<NAEA) OrbMat(i+NMO,j) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
        for(int j=0; j<nread2; j++) {
          in>>dummy;
          if(j<NAEB) OrbMat(i+NMO,j+NAEA) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
       }

     } else {


/*
      nread = fullMOMat?NMO:NAEA;  
      for(int j=0; j<nread; j++) {
        for(int i=0; i<2*NMO; i++) {
         in>>dummy; 
         if(j<NAEA) OrbMat(i,j) = dummy;
         if(in.fail()) {
           app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
           in.close();
           return false;
         }
       }
      } 
      nread = fullMOMat?NMO:NAEB;  
      for(int j=0; j<nread; j++) {
        for(int i=0; i<2*NMO; i++) {
         in>>dummy;
         if(j<NAEB) OrbMat(i,j+NAEA) = dummy;
         if(in.fail()) {
           app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
           in.close();
           return false;
         } 
       }
      }  
 */ 

      nread = fullMOMat?2*NMO:NAEA+NAEB;
      for(int j=0; j<nread; j++) {
        for(int i=0; i<2*NMO; i++) {
         in>>dummy;
         if(j<NAEA+NAEB) OrbMat(i,j) = dummy;
         if(in.fail()) {
           app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
           in.close();
           return false;
         }
       }
      }

     } // Cstyle
   } // wfntype 
  } // rotate

  in.close();

  return setup_local();

}

bool PureSingleDeterminant::setup_local()
{

  if(init_type == "ground" && !rotated_hamiltonian && !readHamFromFile ) {
    occup_alpha.clear();
    occup_beta.clear();
    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);
    for(int i=0; i<NAEA; i++) occup_alpha[i]=i;
    for(int i=0; i<NAEB; i++) occup_beta[i]=NMO+i;
  } else {
    if((occup_alpha.size() != NAEA ) || (occup_beta.size() != NAEB)) {
      app_error()<<"ERROR in PureSingleDeterminant::initFromAscii. occup_alpha/occup_beta wrongly defined. " <<std::endl;  
      return false;
    }
  } 

  std::sort(occup_beta.begin(),occup_beta.end());
  std::sort(occup_alpha.begin(),occup_alpha.end());

  isOcc_alpha.clear();
  isOcc_beta.clear();
  for(IndexType i=0; i<2*NMO; i++) isOcc_alpha[i]=false;
  for(IndexType i=0; i<2*NMO; i++) isOcc_beta[i]=false;
  for(int i=0; i<occup_alpha.size(); i++) isOcc_alpha[ occup_alpha[i] ]=true;
  for(int i=0; i<occup_beta.size(); i++) isOcc_beta[ occup_beta[i] ]=true;

  //trial_density_matrix.resize(NAEA+NAEB,NMO);
  trial_density_matrix.resize(2*NMO,NMO);

  //mixed_density_matrix.resize(NAEA+NAEB,NMO);
  mixed_density_matrix.resize(2*NMO,NMO);

  // not used right now 
  overlap_inv.resize(NAEA+NAEB,NAEA); 


  // temporary storage
  S0.resize(NAEA,NAEA); 
  S1.resize(NAEB,NAEB); 
  SS0.resize(2*NMO,NAEA); 
  V0.resize(2*NMO*NMO); 

  //SpHij.setDims(1,2*NMO*NMO);
  if(!readHamFromFile) SMSpHijkl.setDims(2*NMO*NMO,2*NMO*NMO);
  if(!readHamFromFile) SMSpHijkl.setup(head_of_nodes,name+std::string("SMSpHijkl"),TG.getNodeCommLocal());

  Cwork.resize(2*NMO);
  pivot.resize(2*NMO);

  Iwork.resize(NAEA+NAEB);

  return true;
}

bool PureSingleDeterminant::setup(HamPtr cur) {
  return getHamiltonian(cur);
}

bool PureSingleDeterminant::getHamiltonian(HamPtr h)
{
  ham0 = h;
  sHam = dynamic_cast<SparseGeneralHamiltonian*>(h);
  if(!sHam) { 
    app_error()<<" Error in PureSingleDeterminant::getHamiltonian. \n"
               <<" Hamiltonian associated with PureSingleDeterminant must of the \n"
               <<" type SparseGeneralHamiltonian. \n";
    return false;
  }

  spinRestricted = sHam->RHF();
  if(NCA==NCB && NAEA == NAEB && occup_alpha.size() == occup_beta.size() && spinRestricted && wfntype==0 ) {
    closed_shell = true;
    for(int i=0; i<occup_alpha.size(); i++)
      if( occup_alpha[i]+NMO != occup_beta[i] )
        closed_shell = false;
  } else
    closed_shell = false;

  if(closed_shell) {
    app_log()<<"Found closed shell system. " <<std::endl;
  }

  NuclearCoulombEnergy = static_cast<ValueType>(sHam->NuclearCoulombEnergy);
  if(!readHamFromFile) {

    if(rotated_hamiltonian) {
      app_log()<<" PureSingleDeterminant - Creating Hamiltonian for Rotated Determinant. \n"; 
      // no RHF/GHF walker yet
      int wtype = closed_shell?0:1;   
      if(!sHam->createHamiltonianForGeneralDeterminant(wtype,OrbMat,haj,SMSpHijkl,cutoff)) {
        app_error()<<"Error in createHamiltonianForGeneralDeterminant. \n";
        return false;
      }
    } else {
      app_log()<<" PureSingleDeterminant - Creating Hamiltonian for Pure Determinant. \n"; 
      if(!sHam->createHamiltonianForPureDeterminant(isOcc_alpha,isOcc_beta,hij,SMSpHijkl,cutoff,closed_shell)) {
        app_error()<<"Error in createHamiltonianForPureDeterminant. \n";
        return false;
      }
    }
  }

  // is this correct if hamiltonian actually implements closed_shell?
  // check that SMSpHijkl.rows is consistent with use below
  // FIX FIX FIX
  split_Ham_rows(SMSpHijkl.rows(),SMSpHijkl.rowIndex_begin(),ik0,ikN);
  pik0 = *(SMSpHijkl.row_index()+ik0);

  if(closed_shell)
    local_buff.resize(NMO*NMO+3); 
  else
    local_buff.resize(2*NMO*NMO+3); 

  app_log()<<std::endl <<"*********************************************************************: \n"
           <<" PureSingleDeterminant: \n"
           <<"     Number of terms and memory usage of hij:    " <<(rotated_hamiltonian?haj.size():hij.size()) <<"  " <<(rotated_hamiltonian?haj.size():hij.size())*sizeof(s1D<ValueType>)/1.0e6 <<"  MB. " <<std::endl
           <<"     Number of terms and memory usage of Vijkl:  " <<SMSpHijkl.size() <<"  " <<SMSpHijkl.size()*sizeof(s2D<ValueType>)/1.0e6 <<"  MB. " <<std::endl; 

    ComplexType e1,e2,o1,o2;
    HF.resize(2*NMO,NAEA);
    for(int i=0; i<NAEA; i++) HF(occup_alpha[i],i)=ComplexType(1.0,0.0);
    for(int i=0; i<NAEB; i++) HF(occup_beta[i],i)=ComplexType(1.0,0.0);
    evaluateLocalEnergy(HF.data(),e1,e2,o1,o2);

  app_log()<<" Ehf:     " <<std::setprecision(12) <<e1+e2  <<"  \n " //<<ea+eb <<std::endl
           <<" Ekin:     " <<std::setprecision(12) <<e1    <<"  \n " //<<ea <<std::endl
           <<" Epot:     " <<std::setprecision(12) <<e2    <<"  \n " // <<eb <<std::endl
           <<"*********************************************************************: \n" <<std::endl <<std::endl;

#ifdef AFQMC_TIMER
    Timer.reset("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix");
    Timer.reset("PureSingleDeterminant:evaluateLocalEnergy");
#endif

  if(rank() == 0) bool wrote = hdf_write();

  return true;

}

bool PureSingleDeterminant::hdf_write()
{

  if(hdf_write_file == std::string("")) return true;
  if(readHamFromFile) return true;

  hdf_archive dump(myComm);
  if(!dump.create(hdf_write_file)) {
    app_error()<<" Error opening restart file in PureSingleDeterminant. \n"; 
    return false;
  }

  std::string path = "/Wavefunctions/PureSingleDeterminant";
  if(dump.is_group( path )) {
    app_error()<<" ERROR: H5Group /Wavefunctions/PureSingleDeterminant already exists in restart file. Not over-writing data in file. \n";
    return false;
  }

  dump.push("Wavefunctions");
  dump.push("PureSingleDeterminant");

  std::vector<int> Idata(8);
  Idata[0]=hij.size();
  Idata[1]=SMSpHijkl.size();
  Idata[2]=SMSpHijkl.rows();
  Idata[3]=SMSpHijkl.cols();
  Idata[4]=NMO;
  Idata[5]=NAEA;
  Idata[6]=NAEB;
  Idata[7]=spinRestricted?(0):(1);
  dump.write(Idata,"dims");

  Idata.resize(NAEA+NAEB);
  for(int i=0; i<NAEA; i++) Idata[i] = occup_alpha[i];
  for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) Idata[i] = occup_beta[j];
  dump.write(Idata,"occups");

  std::vector<IndexType> ivec;
  ivec.resize(hij.size());
  for(int i=0; i<hij.size(); i++)
    std::tie (ivec[i],std::ignore) = hij[i];
  dump.write(ivec,"hij_indx");

  std::vector<ValueType> vvec;
  vvec.resize(hij.size());
  for(int i=0; i<hij.size(); i++)
    std::tie (std::ignore,vvec[i]) = hij[i];
  dump.write(vvec,"hij");

  dump.write(*(SMSpHijkl.getVals()),"SpHijkl_vals");
  dump.write(*(SMSpHijkl.getCols()),"SpHijkl_cols");
  dump.write(*(SMSpHijkl.getRowIndex()),"SpHijkl_rowIndex");

  dump.pop();
  dump.pop();

  dump.flush();
  dump.close();

  return true;
}


bool PureSingleDeterminant::initFromHDF5(hdf_archive& read,const std::string& tag)
{
  SMSpHijkl.setup(head_of_nodes,"SMSpHijkl",TG.getNodeCommLocal());
  if(head_of_nodes) {

    std::string path = "/Wavefunctions/PureSingleDeterminant";
    if(!read.is_group( path )) {
      app_error()<<" ERROR: H5Group /Wavefunctions/PureSingleDeterminant does not exist. \n"; 
      return false;
    }

    if(!read.push("Wavefunctions")) return false;
    if(!read.push("PureSingleDeterminant")) return false;

    std::vector<int> Idata(8);
    if(!read.read(Idata,"dims")) return false;
    hij.resize(Idata[0]);

    NCA = NCB = 0;
    if(NMO < 0) NMO = Idata[4];
    if(NAEA < 0) NAEA = Idata[5];
    if(NAEB < 0) NAEB = Idata[6];
    if(Idata[4] != NMO) {
      app_error()<<" ERROR: NMO differs from value in wfn file. \n";
      return false;
    }
    if(Idata[5] != NAEA) {
      app_error()<<" ERROR: NAEA differs from value in wfn file. \n";
      return false;
    }
    if(Idata[6] != NAEB) {
      app_error()<<" ERROR: NAEB differs from value in wfn file. \n";
      return false;
    }
    spinRestricted = (Idata[7]==0)?(true):(false);

    SMSpHijkl.setDims(Idata[2],Idata[3]);
    if(!SMSpHijkl.allocate_serial(Idata[1])) return false;
    SMSpHijkl.resize_serial(Idata[1]);

    myComm->bcast(Idata);   

    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);
    Idata.resize(NAEA+NAEB);
    if(!read.read(Idata,"occups")) return false;
    for(int i=0; i<NAEA; i++) occup_alpha[i] = Idata[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) occup_beta[j] = Idata[i];

    myComm->bcast(Idata);   

    std::vector<int> ivec;
    ivec.resize(hij.size());
    if(!read.read(ivec,"hij_indx")) return false;
    for(int i=0; i<hij.size(); i++)
      std::get<0>(hij[i]) = ivec[i]; 
    myComm->bcast(ivec);

    std::vector<ValueType> vvec;
    vvec.resize(hij.size());
    if(!read.read(vvec,"hij")) return false;
    for(int i=0; i<hij.size(); i++)
      std::get<1>(hij[i]) = vvec[i]; 
    myComm->bcast(vvec);

    if(!read.read(*(SMSpHijkl.getVals()),"SpHijkl_vals")) return false;
    if(!read.read(*(SMSpHijkl.getCols()),"SpHijkl_cols")) return false;
    if(!read.read(*(SMSpHijkl.getRowIndex()),"SpHijkl_rowIndex")) return false;
    SMSpHijkl.setRowsFromRowIndex();

    myComm->barrier();
    myComm->barrier();

    read.pop();
    read.pop();

  } else {

    std::vector<int> Idata(8);
    myComm->bcast(Idata);

    hij.resize(Idata[0]);
    SMSpHijkl.setDims(Idata[2],Idata[3]);
    NCA = NCB = 0;
    if(NMO < 0) NMO = Idata[4];
    if(NAEA < 0) NAEA = Idata[5];
    if(NAEB < 0) NAEB = Idata[6];
    spinRestricted = (Idata[7]==0)?(true):(false);

    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);
    Idata.resize(NAEA+NAEB);
    myComm->bcast(Idata);
    for(int i=0; i<NAEA; i++) occup_alpha[i] = Idata[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) occup_beta[j] = Idata[i];     

    std::vector<int> ivec;
    ivec.resize(hij.size());
    myComm->bcast(ivec);
    for(int i=0; i<hij.size(); i++)
      std::get<0>(hij[i]) = ivec[i];

    std::vector<ValueType> vvec;
    vvec.resize(hij.size());
    myComm->bcast(vvec);
    for(int i=0; i<hij.size(); i++)
      std::get<1>(hij[i]) = vvec[i];

    myComm->barrier();
    if(!SMSpHijkl.initializeChildren()) return false;
    myComm->barrier();

  }
  readHamFromFile=true;
  return setup_local();
}

// approximate the k-set balance partition problem for ordered sets. k0/kN are the partition of the local core
void PureSingleDeterminant::split_Ham_rows(IndexType N, ComplexSMSpMat::int_iterator indx, IndexType& k0, IndexType& kN)
{
  if(ncores_per_TG == 1) {
    k0=0;
    kN=N;
    return;
  } 

  std::vector<ComplexSMSpMat::indxType> subsets(ncores_per_TG+1); 
  if(core_rank==0) {

    balance_partition_ordered_set(N,&(*indx),subsets); 

    if(TG.getTGNumber()==0) {
      std::cout<<" Ham split over cores in TG: \n";
      std::cout<<" Index: \n"; 
      for(int i=0; i<ncores_per_TG+1; i++)
        std::cout<<subsets[i] <<" ";
      std::cout<<std::endl;
      std::cout<<" Number of terms per core: \n"; 
      for(int i=0; i<ncores_per_TG; i++)
        std::cout<<*(indx+subsets[i+1]) - *(indx+subsets[i]) <<" ";
      std::cout<<std::endl; 
//cout<<"Dump: \n";
//for(int i=0; i<N; i++)
//cout<<i <<" " <<*(indx+i) <<std::endl;
    }

  } 

  // since all cores on a node share the same Hamiltonian regardless of TG
  myComm->bcast(subsets,TG.getNodeCommLocal());
  k0 = subsets[core_rank];
  kN = subsets[core_rank+1];

}

void PureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool full)
{

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);

  } else {

    // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
    // lots of simplification because A is an identity matrix with possibly exchanged columns, 
    // look at derivations for information

    // copy rows corresponding to occupied orbitals to S0
    //   S0(i,:) = SlaterMat(ik,:), where ik is the ith occupied orbital    
    ComplexMatrix::iterator itS0 = S0.begin(); 
    for(VIndexit it=occup_alpha.begin(); it!=occup_alpha.end(); it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEA;
      for(; itSM!=itSMend; ++itSM, ++itS0) 
        *itS0 = *itSM;
    }
  }

  ovl_alpha = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

  // SS0 = SlaterMat * S0
  DenseMatrixOperators::product(NMO,NAEA,NAEA,SlaterMat,NAEA,S0.data(),NAEA,SS0.data(),NAEA);

  for(ComplexMatrix::iterator it=mixed_density_matrix.begin(); it!=mixed_density_matrix.end(); it++) *it=zero;

 
  if(rotated_hamiltonian && full) {

    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,mixed_density_matrix.data(),NMO);
    DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data(),NMO);

  } else {

    // copy to mixed_density_matrix 
    for(int i=0; i<NAEA; i++)
     for(int j=0; j<NMO; j++)
      mixed_density_matrix(occup_alpha[i],j) = SS0(j,i);

  }

  if(closed_shell) {
    ovl_beta=ovl_alpha;
    // once this fully works, you will not need the beta sector at all
    std::copy(mixed_density_matrix.begin(),mixed_density_matrix.begin()+NMO*NMO,mixed_density_matrix.begin()+NMO*NMO); 
    return;
  }

  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B    
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NAEA*NMO,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);
   
  } else {
  
    // repeat for beta
    ComplexMatrix::iterator itS1 = S1.begin(); 
    for(VIndexit it=occup_beta.begin(); it!=occup_beta.end(); it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEB;
      for(; itSM!=itSMend; ++itSM, ++itS1) 
        *itS1 = *itSM;
    }
  }

  ovl_beta = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

  // SS0(beta) = SlaterMat(beta) * S1
  DenseMatrixOperators::product(NMO,NAEB,NAEB,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,SS0.data()+NAEA*NMO,NAEA);

  if(rotated_hamiltonian && full) {

    // G(beta) = SS0*transpose(conjg(Dbeta)) 
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,mixed_density_matrix.data()+NMO*NMO,NMO);
    DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data()+NMO*NMO,NMO);

  } else {

    // copy to mixed_density_matrix 
    for(int i=0; i<NAEB; i++)
     for(int j=0; j<NMO; j++)
      mixed_density_matrix(occup_beta[i],j) = SS0(j+NMO,i);

  }

} 

void PureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, ComplexType* dm, bool full)
{

  // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
  // lots of simplification because A is an identity matrix with possibly exchanged columns, 
  // look at derivations for information

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);

  } else {

    // copy rows corresponding to occupied orbitals to S0
    //   S0(i,:) = SlaterMat(ik,:), where ik is the ith occupied orbital    
    ComplexMatrix::iterator itS0 = S0.begin(); 
    for(VIndexit it=occup_alpha.begin(); it!=occup_alpha.end(); it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEA;
      for(; itSM!=itSMend; ++itSM, ++itS0) 
        *itS0 = *itSM;
    }
  }

  ovl_alpha = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());
  
  // SS0 = SlaterMat * S0
  DenseMatrixOperators::product(NMO,NAEA,NAEA,SlaterMat,NAEA,S0.data(),NAEA,SS0.data(),NAEA);


  if(rotated_hamiltonian && full) {

    int nt = closed_shell?NMO:2*NMO;
    for(int i=0; i<nt*NMO; i++) *(dm+i)=zero;
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,dm,NMO);
    DenseMatrixOperators::transpose(NMO,dm,NMO);

  } else {

    int nt = closed_shell?NAEA:(NAEA+NAEB);
    for(int i=0; i<nt*NMO; i++) *(dm+i)=zero;
    // copy to mixed_density_matrix 
    for(int i=0; i<NAEA; i++)
     for(int j=0; j<NMO; j++)
      *(dm++) = SS0(j,i); 
      //*(dm+occup_alpha[i]*NMO+j) = SS0(j,i); 

  }

  if(closed_shell) {
    ovl_beta=ovl_alpha;
    return;
  }

  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B    
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NAEA*NMO,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);
    
  } else {
    // repeat for beta
    ComplexMatrix::iterator itS1 = S1.begin(); 
    for(VIndexit it=occup_beta.begin(); it!=occup_beta.end(); it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEB;
      for(; itSM!=itSMend; ++itSM, ++itS1) 
        *itS1 = *itSM;
    }
  }
  ovl_beta = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

  // SS0(beta) = SlaterMat(beta) * S1
  DenseMatrixOperators::product(NMO,NAEB,NAEB,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,SS0.data()+NAEA*NMO,NAEA);

  if(rotated_hamiltonian && full) {

    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,dm+NMO*NMO,NMO);
    DenseMatrixOperators::transpose(NMO,dm+NMO*NMO,NMO);

  } else {

    // copy to mixed_density_matrix 
    for(int i=0; i<NAEB; i++)
     for(int j=0; j<NMO; j++)
      *(dm++) = SS0(j,i); 
      //*(dm+occup_beta[i]*NMO+j) = SS0(j,i); 

  }
}
    
// on exit, commBuff contains the green function information for all walkers in appropriate place   
// buffer: (this routine assumes the buffer is large enough) 
//    wlksz: size of the data of a walker
//    gfoffset: location of the gf block within a walker
//
void PureSingleDeterminant::evaluateOneBodyMixedDensityMatrix(WalkerHandlerBase* wset, ComplexSMVector* buf, int wlksz, int gfoffset, bool full)
{

  // nwalkers > ncores_per_TG: serial computation with round robin distribution
  // nwalkers < ncores_per_TG: parallel computation with round robin distribution

  ComplexType oa,ob;
  int nw = wset->numWalkers(true), cnt=0;
// Careful here, only returns full matrix if (rotated_hamiltonian&&full), otherwise only NAEA*NMO
// Also careful with closed_shell=true, only returns alpha sector
  int sz = 2 + ((full&&rotated_hamiltonian)?NMO*NMO:NAEA*NMO);
  if(!closed_shell) sz += (full&&rotated_hamiltonian)?NMO*NMO:NAEB*NMO; 
  // in Buff:  {..., oa, ob, DM, ...} for each walker 
  for(int i=0; i<nw; i++) {
    if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;  
    if(cnt%ncores_per_TG == core_rank ) {
      local_evaluateOneBodyMixedDensityMatrix(wset->getSM(i),local_buff[0],local_buff[1],local_buff.data()+2,full);
// Do I need a lock?
      std::copy(local_buff.begin(),local_buff.begin()+sz, buf->begin()+cnt*wlksz+gfoffset);
      //{
      //  boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buffer->getMutex()));
      //  std::copy(local_buff.begin(),local_buff.begin()+sz, buf->begin()+cnt*wlksz+gfoffset);
      //}
    }
    ++cnt;
  }
  
} 

void PureSingleDeterminant::local_evaluateOneBodyTrialDensityMatrix(bool full)
{

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  for(int i=0; i<2*NMO; i++)
   for(int j=0; j<NMO; j++)
    trial_density_matrix(i,j) = zero; 

  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,OrbMat.data(),NAEA,zero,S0.data(),NAEA);

    ComplexType ovl = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());
  
    // SS0 = SlaterMat * S0
    DenseMatrixOperators::product(NMO,NAEA,NAEA,OrbMat.data(),NAEA,S0.data(),NAEA,SS0.data(),NAEA);

    if(full) {

      // G(alpha) = SS0*transpose(conjg(Dalpha))
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,trial_density_matrix.data(),NMO);

      DenseMatrixOperators::transpose(NMO,trial_density_matrix.data(),NMO);

    } else {

      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        trial_density_matrix(i,j) = SS0(j,i);
 
    }

    if(closed_shell) {
      std::copy(trial_density_matrix.begin(),trial_density_matrix.begin()+NMO*NMO,trial_density_matrix.begin()+NMO*NMO);
      return;
    }

    // S0 = transpose(conjg(A))*B    
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,S1.data(),NAEB);
    
    ovl = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

    // SS0(beta) = SlaterMat(beta) * S1
    DenseMatrixOperators::product(NMO,NAEB,NAEB,OrbMat.data()+NAEA*NMO,NAEA,S1.data(),NAEB,SS0.data()+NAEA*NMO,NAEA);

    if(full) {

      // G(beta) = SS0*transpose(conjg(Dbeta))
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,trial_density_matrix.data()+NMO*NMO,NMO);
      DenseMatrixOperators::transpose(NMO,trial_density_matrix.data()+NMO*NMO,NMO);

    } else {

      // copy to mixed_density_matrix 
      // Careful here!!!!! Only writting NAEB terms 
      for(int i=0; i<NAEB; i++)
       for(int j=0; j<NMO; j++)
        trial_density_matrix(i+NMO,j) = SS0(j+NMO,i);
   
    }

  } else {
    for(int i=0; i<NAEA; i++)
      trial_density_matrix(occup_alpha[i],occup_alpha[i]) = one;

    for(int i=0; i<NAEB; i++)
      trial_density_matrix(occup_beta[i],occup_beta[i]-NMO) = one;
  }
} 
  

  void PureSingleDeterminant::evaluateOverlap(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

    const ComplexType one = ComplexType(1.0);
    const ComplexType zero = ComplexType(0.0); 

    if(rotated_hamiltonian) {
    
      // S0 = transpose(conjg(A))*B
      DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);
           
    } else {

      ComplexMatrix::iterator itS0 = S0.begin(); 
      for(VIndexit it=occup_alpha.begin(); it!=occup_alpha.end(); it++) {
        const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
        const ComplexType* itSMend = itSM+NAEA;
        for(; itSM!=itSMend; ++itSM, ++itS0)
          *itS0 = *itSM;
      }
    }

    ovl_alpha = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

    if(closed_shell) {
      ovl_beta=ovl_alpha;
      return;
    }
  
    if(rotated_hamiltonian) {
    
      // S0 = transpose(conjg(A))*B
      DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NAEA*NMO,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);
           
    } else {
      ComplexMatrix::iterator itS1 = S1.begin(); 
      for(VIndexit it=occup_beta.begin(); it!=occup_beta.end(); it++) {
        const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
        const ComplexType* itSMend = itSM+NAEB;
        for(; itSM!=itSMend; ++itSM, ++itS1)                                                                                 
          *itS1 = *itSM;
      }
    }
    ovl_beta = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());
  
  }

  void PureSingleDeterminant::dist_evaluateOverlap(WalkerHandlerBase* wset , bool first, const int n ) {
    // testing
    //if(TG.getCoreRank()==0) serial_evaluateOverlap(wset,first,n);

    // right now use round-robin with serial evaluations of overlaps 

    int nw = wset->numWalkers(true);
    if(nw==0) return;
    ComplexType ovlp_a,ovlp_b;
    for(int i=0, cnt=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
      if(cnt%ncores_per_TG == core_rank) {
        evaluateOverlap(wset->getSM(i),ovlp_a,ovlp_b,n);
        if(first)
          wset->setOvlp(i,ovlp_a,ovlp_b);
        else
          wset->setOvlp2(i,ovlp_a,ovlp_b);
      } 
      cnt++;
    }

  }

  void PureSingleDeterminant::dist_evaluateLocalEnergy(WalkerHandlerBase* wset , bool first, const int n ) {
    // testing
    //if(TG.getCoreRank()==0) serial_evaluateLocalEnergy(wset,first,n);
 
    //return;

    if(distribute_Ham)
      APP_ABORT(" Error: PureSingleDeterminant::dist_evaluateLocalEnergy with nnodes_per_TG > 1 not implemented yet");    

  // structure in TG [eloc, oa, ob, G(1:{2*}NMO*NMO)] 

    int sz = 3 + NAEA*NMO;
    if(!closed_shell) sz += NAEB*NMO;
    int nw0 = wset->numWalkers(false);
    int cnt=0;
    int nwglobal = nw0*sz;
    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    int nr1=1, nc1=2*NMO*NMO;
    ComplexType oa,ob,epot,ekin = zero;    
    ComplexSMVector* buffer = TG.commBuff;
 
    TG.resize_buffer(nwglobal);
    nwglobal /= sz;

    // calculate [oa,ob,G] for all walkers
    evaluateOneBodyMixedDensityMatrix(wset,buffer,sz,1,false);

    if(nwglobal > local_buff.size())
      local_buff.resize(nwglobal);

    if(core_rank==0) {
      ComplexType *val = buffer->values();
      for(int i=0; i<nw0; i++,val+=sz) *val = 0;
    } 

    // synchronize
    TG.local_barrier(); 

    // add 1-body term, round-robin
    for(int wlk=0; wlk<nw0; wlk++) { 
      ekin=zero;
      mixed_density_matrix = zero;
      ComplexType* ptr = buffer->values() + wlk*sz + 3;
      ComplexMatrix::iterator itG = mixed_density_matrix.begin();
      for(int i=0; i<NAEA; i++, ptr+=NMO )  // can use std::copy here or axpy 
       std::copy( ptr, ptr+NMO, mixed_density_matrix[occup_alpha[i]]); 
      if(!closed_shell) {
       for(int i=0; i<NAEB; i++, ptr+=NMO)
        std::copy( ptr, ptr+NMO, mixed_density_matrix[occup_beta[i]]);
      }
      if( wlk%ncores_per_TG == core_rank ) { 
        itG = mixed_density_matrix.begin();
        if(rotated_hamiltonian) {
          std::vector<s1D<ComplexType> >::iterator  end1 = haj.end();
          for(std::vector<s1D<ComplexType> >::iterator it = haj.begin(); it != end1; it++)
            ekin += *(itG + std::get<0>(*it)) * std::get<1>(*it);
        } else {
          s1Dit end1 = hij.end();
          for(s1Dit it = hij.begin(); it != end1; it++)
            ekin += *(itG + std::get<0>(*it)) * std::get<1>(*it);
        }       
        ekin+=NuclearCoulombEnergy;
      }

      SparseMatrixOperators::product_SpMatV(int(ikN-ik0),SMSpHijkl.cols(),one,SMSpHijkl.values() + pik0, SMSpHijkl.column_data() + pik0, SMSpHijkl.row_index() + ik0,  mixed_density_matrix.data(),zero,V0.data()+ik0);
      itG = mixed_density_matrix.begin()+ik0;
      ComplexVector::iterator itV = V0.begin()+ik0;
      epot = zero;
      for(int i=ik0; i<ikN; i++,++itG,++itV) epot += (*itV) * (*itG);
      local_buff[wlk] = ekin + 0.5*epot;
    }

    {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buffer->getMutex()));
      zaxpy (nw0, one, local_buff.data(), 1, buffer->values(), sz);
    }
  
    // synchronize
    TG.local_barrier(); 

    int currnw=nw0;
    for(int tgi = 0; tgi<nnodes_per_TG; tgi++) {

      if(tgi > 0) {
        // adds local contribution to 2-body 
        for(int wlk=0; wlk<currnw; wlk++) { 
          mixed_density_matrix = zero;
          ComplexType* ptr = buffer->values() + wlk*sz + 3;
          for(int i=0; i<NAEA; i++, ptr+=NMO )  // can use std::copy here or axpy 
           std::copy( ptr, ptr+NMO, mixed_density_matrix[occup_alpha[i]]);
          if(!closed_shell) {
           for(int i=0; i<NAEB; i++, ptr+=NMO)
            std::copy( ptr, ptr+NMO, mixed_density_matrix[occup_beta[i]]);
          }

          SparseMatrixOperators::product_SpMatV(int(ikN-ik0),SMSpHijkl.cols(),one,SMSpHijkl.values() + pik0, SMSpHijkl.column_data() + pik0, SMSpHijkl.row_index() + ik0,  mixed_density_matrix.data(),zero,V0.data()+ik0);
          ComplexMatrix::iterator itG = mixed_density_matrix.begin()+ik0;
          ComplexVector::iterator itV = V0.begin()+ik0;
          epot = zero;
          for(int i=ik0; i<ikN; i++,++itG,++itV) epot += (*itV) * (*itG);
          local_buff[wlk] = ekin + 0.5*epot;
          //local_buff[wlk] = ekin + 0.5*zdotu(int(ikN-ik0),mixed_density_matrix.data()+ik0,1,V0.data()+ik0,1);
          //local_buff[wlk] = 0; 
        }
        {
          boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buffer->getMutex()));
          zaxpy (currnw, one, local_buff.data(), 1, buffer->values(), sz);
        }
      }

      // rotate buffer 
      if(distribute_Ham) TG.rotate_buffer(currnw,sz);
    }

    // synchronize: rotate_buffer calls barrier, so no need to do it again 
    if(!distribute_Ham) TG.local_barrier(); 

    if(core_rank != 0) return;
    
    int nw = wset->numWalkers(true);
    ComplexType* ptr = buffer->values();
    for(int i=0,cnt=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
      if(first)
        wset->setWalker(i,*ptr,*(ptr+1),*(ptr+2));
      else {
        wset->setEloc2(i,*ptr);
        wset->setOvlp2(i,*(ptr+1),*(ptr+2));
      }
      ptr+=sz;
      cnt++;
    }
  }

  void PureSingleDeterminant::evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

#if _DEBUG_AFQMC_
   // assert(SlaterMat.rows() != NMO && SlaterMat.cols() != NMO )
 #endif
#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix"); 
#endif
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,ovl_alpha,ovl_beta,false);
#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix"); 
#endif

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:evaluateLocalEnergy"); 
#endif

    ekin = 0;
    ComplexMatrix::iterator itG = mixed_density_matrix.begin();
    if(rotated_hamiltonian) {
      std::vector<s1D<ComplexType> >::iterator  end1 = haj.end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj.begin(); it != end1; it++) 
        ekin += *(itG + std::get<0>(*it)) * std::get<1>(*it);
    } else {
      s1Dit end1 = hij.end();
      for(s1Dit it = hij.begin(); it != end1; it++) 
        ekin += *(itG + std::get<0>(*it)) * std::get<1>(*it);
    }

// NOTES:
// Expand this routine to enable the following options:
//   1. single precision matrix
//       - in this case copy mixed_density_matrix to a single precision vector and call the SP routine.
//       - this should be automatically detected when ValueType != METype. (to be introduced)
//   2. matrix stored in real type, detected when isComplex(METype)==false
//       - do this by calling a SpMatM routine, where you interpret mixed_density_matrix as a 2 column matrix     
//       - combine with DP->SP option  

    epot = 0; 
    int nr1=1, nc1=2*NMO*NMO;
    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(nc1,nc1,one,SMSpHijkl,mixed_density_matrix.data(),zero,V0.data()); 
    itG = mixed_density_matrix.begin();
    ComplexVector::iterator itV = V0.begin();
    for(int i=0; i<nc1; i++,++itG,++itV) epot += (*itV) * (*itG); 
    epot = 0.5*epot+NuclearCoulombEnergy;   
  
#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:evaluateLocalEnergy"); 
#endif
    
  }

  void PureSingleDeterminant::evaluateLocalEnergy(bool addBetaBeta, RealType dt, const ComplexType* SlaterMat, const ComplexSMSpMat& Spvn, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta,  bool transposed, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix"); 
#endif
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,ovl_alpha,ovl_beta,false);
#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix"); 
#endif

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:evaluateLocalEnergy_factHam"); 
#endif

    ekin = 0;
    ComplexMatrix::iterator itG = mixed_density_matrix.begin();
    if(rotated_hamiltonian) {
      std::vector<s1D<ComplexType> >::iterator  end1 = haj.end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj.begin(); it != end1; it++)
        ekin += *(itG + std::get<0>(*it)) * std::get<1>(*it);
    } else {
      s1Dit end1 = hij.end();
      for(s1Dit it = hij.begin(); it != end1; it++)                                 
        ekin += *(itG + std::get<0>(*it)) * std::get<1>(*it);
    }

    epot = 0; 
    ComplexType epot_exch = 0; 
    // Coulomb:   sum_n ( sum_ik v^n_ik * G_ik )^2
    if(transposed)
     T1.resize(Spvn.rows());
    else
     T1.resize(Spvn.cols());
    ComplexType one = ComplexType(1.0,0.0);
    ComplexType one_ = ComplexType(1.0,0.0);
    ComplexType two = ComplexType(2.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(closed_shell) one_ = ComplexType(2.0,0.0);
    if(transposed) {
      SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one_,Spvn,mixed_density_matrix.data(),zero,T1.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,mixed_density_matrix.data()+NMO*NMO,one,T1.data());
    } else {
      SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one_,Spvn,mixed_density_matrix.data(),zero,T1.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,mixed_density_matrix.data()+NMO*NMO,one,T1.data());
    }
    for(std::vector<ComplexType>::iterator it = T1.begin(); it<T1.end(); ++it) epot += (*it)*(*it); 

    // Exchange:  - sum_n  sum_ijkl v^n_ik v^n_lj Gil Gjk   
    itG = mixed_density_matrix.begin();
    register int NMO2 = NMO*NMO;
    if(transposed) {

      APP_ABORT("Error: evaluateLocalEnergy disabled with transposed Spvn. Very Slow \n\n\n");

      // Spvn( n, ik )
      for(ComplexSMSpMat::const_int_iterator itri = Spvn.rowIndex_begin(); (itri+1) != Spvn.rowIndex_end(); itri++) { 
        register int i1 = *itri; 
        register int i3 = *(itri+1); 
        if( i1 == i3 ) continue;
        ComplexSMSpMat::const_int_iterator it1 = Spvn.cols_begin()+i1;
        ComplexSMSpMat::const_iterator itv1 = Spvn.vals_begin()+i1;
        register int i2 = i3; // look for end of alpha-alpha block 
        if( addBetaBeta && !closed_shell) {
          for(int ik = i1; ik < i2; ++ik, ++it1, ++itv1) { 
            ComplexSMSpMat::const_int_iterator it2 = Spvn.cols_begin()+ik;
            ComplexSMSpMat::const_iterator itv2 = Spvn.vals_begin()+ik;
            div_t divresult = div(*it1,NMO);
            register int i_ = divresult.quot; 
            register int k_ = divresult.rem; 
            // diagonal term first
            {  
              divresult = div(*it2,NMO);
              register int l_ = divresult.quot; 
              register int j_ = divresult.rem; 
              epot_exch -= (*itv1) * (*itv2) * 
                  (   ( *(itG + i_*NMO + l_ ) ) 
                    * ( *(itG + j_*NMO + k_ ) ) 
                    + ( *(itG + NMO2 + i_*NMO + l_ ) ) 
                    * ( *(itG + NMO2 + j_*NMO + k_ ) ) ); 
              ++it2;
              ++itv2;
            }   
            for(int lj = ik+1; lj < i2; ++lj, ++it2, ++itv2) { 
              divresult = div(*it2,NMO);
              register int l_ = divresult.quot; 
              register int j_ = divresult.rem; 
              epot_exch -= two * (*itv1) * (*itv2) * 
                  (   ( *(itG + i_*NMO + l_ ) ) 
                    * ( *(itG + j_*NMO + k_ ) ) 
                    + ( *(itG + NMO2 + i_*NMO + l_ ) ) 
                    * ( *(itG + NMO2 + j_*NMO + k_ ) ) ); 
            }
          }
        } else {

          if( !closed_shell ) {
            for(int ik = i1; ik < i3; ++ik, ++it1) { 
              if( (*it1)/NMO >= NMO) {
                i2 = ik;
                break;
              }
            }
            it1 = Spvn.cols_begin()+i1;
          } 
          // now i2 is the first beta-beta element for the current vector 
          for(int ik = i1; ik < i2; ++ik, ++it1, ++itv1) { 
            ComplexSMSpMat::const_int_iterator it2 = Spvn.cols_begin()+ik;
            ComplexSMSpMat::const_iterator itv2 = Spvn.vals_begin()+ik;
            div_t divresult = div(*it1,NMO);
            register int i_ = divresult.quot; 
            register int k_ = divresult.rem; 
            // diagonal term first
            {  
              divresult = div(*it2,NMO);
              register int l_ = divresult.quot; 
              register int j_ = divresult.rem; 
              epot_exch -= (*itv1) * (*itv2) 
                    * ( *(itG + i_*NMO + l_ ) ) 
                    * ( *(itG + j_*NMO + k_ ) );
              ++it2;
              ++itv2;
            }   
            for(int lj = ik+1; lj < i2; ++lj, ++it2, ++itv2) { 
              divresult = div(*it2,NMO);
              register int l_ = divresult.quot; 
              register int j_ = divresult.rem; 
              epot_exch -= two * (*itv1) * (*itv2) 
                    * ( *(itG + i_*NMO + l_ ) ) 
                    * ( *(itG + j_*NMO + k_ ) );
            }
          }

          // if closed_shell, i2==i3 
          // now do remaining beta-beta sector
          for(int ik = i2; ik < i3; ++ik, ++it1, ++itv1) { 
            ComplexSMSpMat::const_int_iterator it2 = Spvn.cols_begin()+ik;
            ComplexSMSpMat::const_iterator itv2 = Spvn.vals_begin()+ik;
            div_t divresult = div(*it1,NMO);
            register int i_ = divresult.quot; 
            register int k_ = divresult.rem; 
            // diagonal term first
            {  
              divresult = div(*it2,NMO);
              register int l_ = divresult.quot - NMO; 
              register int j_ = divresult.rem + NMO; 
              epot_exch -= (*itv1) * (*itv2) 
                    * ( *(itG + i_*NMO + l_ ) ) 
                    * ( *(itG + j_*NMO + k_ ) );
              ++it2;
              ++itv2;
            }   
            for(int lj = ik+1; lj < i2; ++lj, ++it2, ++itv2) { 
              divresult = div(*it2,NMO);
              register int l_ = divresult.quot - NMO; 
              register int j_ = divresult.rem + NMO; 
              epot_exch -= two * (*itv1) * (*itv2) 
                    * ( *(itG + i_*NMO + l_ ) ) 
                    * ( *(itG + j_*NMO + k_ ) );
            }
          }
        }
      }      
 
    } else {

      if(!closed_shell) {
        APP_ABORT("Error: evaluateLocalEnergy not implemented with Spvn for non closed-shelli yet. \n\n\n");
      }
  
      // currently needs care in use
      if(setup_vn_occ_indx) {
        int cnt1 = 0;
        int cnt2 = 0;
        if( std::distance(Spvn.rowIndex_begin(),Spvn.rowIndex_end()) != NMO*NMO+1 ) { 
          app_error()<<" Error: Spvn number of rows: " <<std::distance(Spvn.rowIndex_begin(),Spvn.rowIndex_end()) <<std::endl;
          APP_ABORT("Error: Problems with number of rows in evaluateLocalEnergy with Spvn. \n ");  
        }
        ComplexSMSpMat::const_int_iterator itik = Spvn.rowIndex_begin();
        for(int i=0; i<NMO; i++)
         for(int k=0; k<NMO; k++, ++itik) {
           if( *itik == *(itik+1) ) continue;
           if(isOcc_alpha[i]) cnt1++; 
           if(isOcc_alpha[k]) cnt2++; 
         }
        vn_occ_ik.reserve(cnt1); 
        vn_occ_lj.reserve(cnt2); 
        itik = Spvn.rowIndex_begin();
        for(int i=0; i<NMO; i++)
         for(int k=0; k<NMO; k++, ++itik) {
           if( *itik == *(itik+1) ) continue;
           if(isOcc_alpha[i]) vn_occ_ik.push_back(std::forward_as_tuple( std::distance(Spvn.rowIndex_begin(),itik),i,k)); 
           if(isOcc_alpha[k]) vn_occ_lj.push_back(std::forward_as_tuple( std::distance(Spvn.rowIndex_begin(),itik),i,k)); 
         }
        setup_vn_occ_indx = false;
      } 

      const int* cols = Spvn.column_data();
      const int* indx = Spvn.row_index(); 
      const ComplexType* vals = Spvn.values();
      // Spvn( ik, n )
      for(std::vector<std::tuple<int,int,int>>::iterator it1 = vn_occ_ik.begin(); it1!=vn_occ_ik.end(); ++it1) {
        for(std::vector<std::tuple<int,int,int>>::iterator it2 = vn_occ_lj.begin(); it2!=vn_occ_lj.end(); ++it2) {
          int ik = std::get<0>(*it1);
          int lj = std::get<0>(*it2);
          ComplexType val = SparseMatrixOperators::product_SpVSpV<ComplexType>( indx[ik+1]-indx[ik],cols+indx[ik],vals+indx[ik],indx[lj+1]-indx[lj],cols+indx[lj],vals+indx[lj]);  
          epot_exch -= val 
               * ( *(itG + std::get<1>(*it1)*NMO + std::get<1>(*it2) ) )
               * ( *(itG + std::get<2>(*it2)*NMO + std::get<2>(*it1) ) );          
        }
      } 

    }
    if(closed_shell) epot_exch*=two;
    epot = (epot+epot_exch)*ComplexType(-0.5/dt,0.0)+NuclearCoulombEnergy;   
  
#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:evaluateLocalEnergy_factHam"); 
#endif
    
  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators"); 
#endif
    ComplexType o1,o2;
    const ComplexType *GF=SlaterMat;
    if(needsG) {
     GF = mixed_density_matrix.data(); 
     local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);
    }

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(closed_shell) one = ComplexType(2.0,0.0); 
    if(transposed) {
      SparseMatrixOperators::product_SpMatV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF+NMO*NMO,one,v.data());
    } else { 
      SparseMatrixOperators::product_SpMatTV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators"); 
#endif

  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif
    ComplexType o1,o2;
    const ComplexType *GF=SlaterMat;
    if(needsG) {
     GF = mixed_density_matrix.data();
     local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);
    }

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(closed_shell) one = ComplexType(2.0,0.0); 
    if(transposed) {
      SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF+NMO*NMO,one,v.data());
    } else {
      SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,GF+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int i0, int iN, int pi0, ComplexSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {
    APP_ABORT(" Error: Routine not implemented: PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(ComplexSpMat). \n\n\n");
  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int i0, int iN, int pi0, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {
    const ComplexType *GF = buff+2;
    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(closed_shell) one = ComplexType(2.0,0.0);
    if(!rotated_hamiltonian) {  
      mixed_density_matrix = zero; 
      // copy to mixed_density_matrix 
      buff+=2; 
      for(int i=0; i<NAEA; i++)  // can use std::copy here or axpy 
       for(int j=0; j<NMO; j++)
        mixed_density_matrix(occup_alpha[i],j) = *(buff++);
      if(!closed_shell) {
       for(int i=0; i<NAEB; i++)
        for(int j=0; j<NMO; j++)
         mixed_density_matrix(occup_beta[i],j) = *(buff++);
      }
      GF = mixed_density_matrix.data();
    }

    if(transposed) {
      SparseMatrixOperators::product_SpMatV( int(iN-i0), Spvn.cols(), one, Spvn.values() + pi0, Spvn.column_data() + pi0, Spvn.row_index()+i0, GF+i0, zero, v.data());    
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatV( int(iN-i0), Spvn.cols(), one, Spvn.values() + pi0, Spvn.column_data() + pi0, Spvn.row_index()+i0, GF+NMO*NMO+i0, one, v.data());    
    } else {
      SparseMatrixOperators::product_SpMatTV( int(iN-i0), Spvn.cols(), one, Spvn.values() + pi0, Spvn.column_data() + pi0, Spvn.row_index()+i0, GF+i0, zero, v.data());    
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV( int(iN-i0), Spvn.cols(), one, Spvn.values() + pi0, Spvn.column_data() + pi0, Spvn.row_index()+i0, GF+NMO*NMO+i0, one, v.data());    
    }
  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSpMat& Spvn, std::vector<ComplexType>& v, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators"); 
#endif
    ComplexType o1,o2;
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);


    ComplexMatrix::iterator itG = mixed_density_matrix.begin();
    for(int i=0; i<vn_indx.size()-1; i++) {
      v[i] = static_cast<ComplexType>(0.0); 
      for(int n = vn_indx[i]; n<vn_indx[i+1]; n++) {
        IndexType ik = Index2Mat(std::get<0>(vn[n]),std::get<2>(vn[n])); 
        IndexType jl = Index2Mat(std::get<1>(vn[n]),std::get<3>(vn[n])); 
        IndexType il = Index2Mat(std::get<0>(vn[n]),std::get<3>(vn[n])); 
        IndexType jk = Index2Mat(std::get<1>(vn[n]),std::get<2>(vn[n])); 
        v[i] += (*(itG + ik)) * (*(itG + jl) - *(itG + il)) * (*(itG + jk)) * std::get<4>(vn[n]);
      }
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators"); 
#endif

  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators");
#endif
    ComplexType o1,o2;
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);


    ComplexMatrix::iterator itG = mixed_density_matrix.begin();
    for(int i=0; i<vn_indx.size()-1; i++) {
      v[i] = static_cast<ComplexType>(0.0);
      for(int n = vn_indx[i]; n<vn_indx[i+1]; n++) {
        IndexType ik = Index2Mat(std::get<0>(vn[n]),std::get<2>(vn[n]));
        IndexType jl = Index2Mat(std::get<1>(vn[n]),std::get<3>(vn[n]));
        IndexType il = Index2Mat(std::get<0>(vn[n]),std::get<3>(vn[n]));
        IndexType jk = Index2Mat(std::get<1>(vn[n]),std::get<2>(vn[n]));
        v[i] += (*(itG + ik)) * (*(itG + jl) - *(itG + il)) * (*(itG + jk)) * std::get<4>(vn[n]);
      }
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators");
#endif

  }

  void PureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSpMat& Spvn, std::vector<ComplexType>& v, const int n)
  {

    if(trialDensityMatrix_needsupdate) {
      trialDensityMatrix_needsupdate = false;
      ComplexType o1,o2;
      local_evaluateOneBodyTrialDensityMatrix(true);
    }

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(closed_shell) one = ComplexType(2.0,0.0);
    SparseMatrixOperators::product_SpMatTV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,trial_density_matrix.data(),zero,v.data());
    if(addBetaBeta && !closed_shell)
      SparseMatrixOperators::product_SpMatTV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,trial_density_matrix.data()+NMO*NMO,one,v.data());

  }

  void PureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, const int n)
  {

    if(trialDensityMatrix_needsupdate) {
      trialDensityMatrix_needsupdate = false;
      ComplexType o1,o2;
      local_evaluateOneBodyTrialDensityMatrix(true);
    }

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(closed_shell) one = ComplexType(2.0,0.0);
    SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,trial_density_matrix.data(),zero,v.data());
    if(addBetaBeta && !closed_shell)
      SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,trial_density_matrix.data()+NMO*NMO,one,v.data());

  }

}

      
