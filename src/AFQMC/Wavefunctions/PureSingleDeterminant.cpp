#include<cstdlib>
#include<complex>
#include<iostream>
#include<fstream>
#include <string.h>
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
#include "Numerics/OhmmsBlas.h"
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

    app_log()<<"\n\n --------------- Parsing PureSD input ------------------ \n\n";

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
    m_param.add(write_trial_density_matrix,"trial_density_matrix","std::string");
    m_param.add(cutoff,"cutoff","double");
    m_param.add(initialDet,"initialDetType","int");
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
        wfn_type = atoi((++it)->c_str());
        switch(wfn_type) {
          case 0:
          {
            app_log()<<"Using a RHF-type trial wave-function in SlaterDeterminant. \n";
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
            app_error()<<"Unknown wave-function type in PureSingleDeterminant: " <<wfn_type <<std::endl;
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
    if(wfn_type == 2)
      ncols = NAEA+NAEB;

    OrbMat.resize(2*NMO,ncols);
    ComplexType dummy;
    int nread;

    if (wfn_type == 0 ) {

       if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant.  \n";
              app_error()<<i <<" " <<j <<std::endl;
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
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              return false;
            }
          }
       }

    } else if(wfn_type == 1) {

    
      if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (alpha) \n";
              app_error()<<i <<" " <<j <<std::endl;
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
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (beta) \n";
              app_error()<<i <<" " <<j <<std::endl;
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
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (alpha) \n";
              app_error()<<i <<" " <<j <<std::endl;
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
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (beta) \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              return false;
            }
          }
      }

    } else if(wfn_type == 2) {

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
   } // wfn_type 
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
  SPtrial_density_matrix.resize(2*NMO,NMO);

  temp_density_matrix.resize(NMO,NMO);
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
  if(!readHamFromFile) SMSpHabkl.setDims(2*NMO*NMO,2*NMO*NMO);
  if(!readHamFromFile) SMSpHabkl.setup(head_of_nodes,name+std::string("SMSpHabkl"),TG.getNodeCommLocal());

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
  if(NCA==NCB && NAEA == NAEB && occup_alpha.size() == occup_beta.size() && spinRestricted && wfn_type==0 ) {
    closed_shell = true;
    for(int i=0; i<occup_alpha.size(); i++)
      if( occup_alpha[i]+NMO != occup_beta[i] ) 
        closed_shell = false;
  } else 
    closed_shell = false;

  // walker_type will come from WalkerHandler, but for now only ROHF/UHF is implemented correctly.
  walker_type=1; 
  dm_type = std::max(wfn_type,1); 
  if(closed_shell) {
    app_log()<<" Found closed shell system. " <<std::endl;
    dm_type=0;
  } else {
    app_log()<<" System is not closed shell. " <<std::endl;
  } 

  if(useFacHam) {
    APP_ABORT(" Error: Use of factorized hamiltonian is not implemented in PureSD. \n\n\n");
  }

  NuclearCoulombEnergy = static_cast<ValueType>(sHam->NuclearCoulombEnergy);
  if(!readHamFromFile) {

    if(rotated_hamiltonian) {
      app_log()<<" PureSingleDeterminant - Creating Hamiltonian for Rotated Determinant. \n"; 
      // no RHF/GHF walker yet
      if(!sHam->createHamiltonianForGeneralDeterminant(dm_type,OrbMat,haj,SMSpHabkl,cutoff)) {
        app_error()<<"Error in createHamiltonianForGeneralDeterminant. \n";
        return false;
      }
    } else {
      app_log()<<" PureSingleDeterminant - Creating Hamiltonian for Pure Determinant. \n"; 
      if(!sHam->createHamiltonianForPureDeterminant(dm_type,useFacHam,isOcc_alpha,isOcc_beta,hij,SMSpHijkl,cutoff)) {
        app_error()<<"Error in createHamiltonianForPureDeterminant. \n";
        return false;
      }
    }
  }

  // is this correct if hamiltonian actually implements closed_shell?
  // check that SMSpHijkl.rows is consistent with use below
  // FIX FIX FIX
  if(rotated_hamiltonian) {
    split_Ham_rows(SMSpHabkl.rows(),SMSpHabkl.rowIndex_begin(),ik0,ikN);
    pik0 = *(SMSpHabkl.row_index()+ik0);
  } else {
    split_Ham_rows(SMSpHijkl.rows(),SMSpHijkl.rowIndex_begin(),ik0,ikN);
    pik0 = *(SMSpHijkl.row_index()+ik0);
  }

  if(closed_shell)
    local_buff.resize(NMO*NMO+3); 
  else
    local_buff.resize(2*NMO*NMO+3); 

  app_log()<<std::endl <<"*********************************************************************: \n"
           <<"  PureSingleDeterminant: \n"
           <<"     Number of terms and memory usage of hij:    " <<(rotated_hamiltonian?haj.size():hij.size()) <<"  " <<(rotated_hamiltonian?haj.size():hij.size())*sizeof(s1D<ValueType>)/1.0e6 <<"  MB. " <<std::endl
           <<"     Number of terms and memory usage of Vijkl:  " <<(rotated_hamiltonian?SMSpHabkl.size():SMSpHijkl.size()) <<"  " <<(rotated_hamiltonian?SMSpHabkl.size()*sizeof(s2D<SPComplexType>):SMSpHijkl.size()*sizeof(s2D<SPValueType>))/1.0e6 <<"  MB. " <<std::endl;

    ComplexType e1,e2,o1,o2;
    HF.resize(2*NMO,NAEA);
    HF = ComplexType(0.0,0.0);
    if(rotated_hamiltonian) {
        std::copy(OrbMat.data(),OrbMat.data()+NAEA*NMO,HF.data());
        if(wfn_type==0 || initialDet==0)
            std::copy(OrbMat.data(),OrbMat.data()+NAEA*NMO,HF.data()+NAEA*NMO);
        else
            std::copy(OrbMat.data()+NAEA*NMO,OrbMat.data()+2*NAEA*NMO,HF.data()+NAEA*NMO);
    } else {
        for(int i=0; i<NAEA; i++) HF(occup_alpha[i],i)=ComplexType(1.0,0.0);
        for(int i=0; i<NAEB; i++) HF(occup_beta[i],i)=ComplexType(1.0,0.0);
    }
    evaluateLocalEnergy(HF.data(),e1,e2,o1,o2);

  app_log()<<"  Ehf:      " <<std::setprecision(12) <<e1+e2  <<"  \n" //<<ea+eb <<std::endl
           <<"  Ekin:     " <<std::setprecision(12) <<e1    <<"  \n" //<<ea <<std::endl
           <<"  Epot:     " <<std::setprecision(12) <<e2    <<"  \n" // <<eb <<std::endl
           <<"*********************************************************************: \n" <<std::endl <<std::endl;

  if(write_trial_density_matrix != "" && rank() == 0) {
    local_evaluateOneBodyTrialDensityMatrix(true);
    std::ofstream out(write_trial_density_matrix.c_str());
    out<<"# trial density matrix: NMO, NAEA, NAEB:" <<NMO <<" " <<NAEA <<" " <<NAEB <<"\n";
    for(int i=0; i<2*NMO; i++) {
      for(int j=0; j<NMO; j++)
        out<<trial_density_matrix(i,j) <<" ";
      out<<"\n";
    }
    out.flush();
    out.close();
  }
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
  if(rotated_hamiltonian) { 
   Idata[0]=haj.size();
   Idata[1]=SMSpHabkl.size();
   Idata[2]=SMSpHabkl.rows();
   Idata[3]=SMSpHabkl.cols();
  } else {
   Idata[0]=hij.size();
   Idata[1]=SMSpHijkl.size();
   Idata[2]=SMSpHijkl.rows();
   Idata[3]=SMSpHijkl.cols();
  }
  Idata[4]=NMO;
  Idata[5]=NAEA;
  Idata[6]=NAEB;
  Idata[7]=spinRestricted?(0):(1);
  dump.write(Idata,"dims");

  Idata.resize(NAEA+NAEB);
  for(int i=0; i<NAEA; i++) Idata[i] = occup_alpha[i];
  for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) Idata[i] = occup_beta[j];
  dump.write(Idata,"occups");

  // write wavefunction (mainly for miniapp) 
  if(rotated_hamiltonian) {
    std::vector<ComplexType> vvec(OrbMat.size());
    std::copy(OrbMat.begin(),OrbMat.end(),vvec.begin());
    dump.write(vvec,"Wavefun");
  }

  std::vector<IndexType> ivec;
  if(rotated_hamiltonian) {
    ivec.resize(haj.size());
    for(int i=0; i<haj.size(); i++)
      std::tie (ivec[i],std::ignore) = haj[i];
  } else {
    ivec.resize(hij.size());
    for(int i=0; i<hij.size(); i++)
      std::tie (ivec[i],std::ignore) = hij[i];
  }
  dump.write(ivec,"hij_indx");

  if(rotated_hamiltonian) {
   std::vector<ComplexType> vvec;
   vvec.resize(haj.size());
   for(int i=0; i<haj.size(); i++)
     std::tie (std::ignore,vvec[i]) = haj[i];
    dump.write(vvec,"hij");
  } else {
   std::vector<ValueType> vvec;
   vvec.resize(hij.size());
   for(int i=0; i<hij.size(); i++)
     std::tie (std::ignore,vvec[i]) = hij[i];
    dump.write(vvec,"hij");
  }

  if(rotated_hamiltonian) {
    dump.write(*(SMSpHabkl.getVals()),"SpHijkl_vals");
    dump.write(*(SMSpHabkl.getCols()),"SpHijkl_cols");
    dump.write(*(SMSpHabkl.getRowIndex()),"SpHijkl_rowIndex");
  } else {
    dump.write(*(SMSpHijkl.getVals()),"SpHijkl_vals");
    dump.write(*(SMSpHijkl.getCols()),"SpHijkl_cols");
    dump.write(*(SMSpHijkl.getRowIndex()),"SpHijkl_rowIndex");
  }

  dump.pop();
  dump.pop();

  dump.flush();
  dump.close();

  return true;
}


bool PureSingleDeterminant::initFromHDF5(hdf_archive& read,const std::string& tag)
{
  SMSpHijkl.setup(head_of_nodes,"SMSpHijkl",TG.getNodeCommLocal());
  SMSpHabkl.setup(head_of_nodes,"SMSpHabkl",TG.getNodeCommLocal());
  if(head_of_nodes) {

    std::string path = "/Wavefunctions/PureSingleDeterminant";
    if(!read.is_group( path )) {
      app_error()<<" ERROR: H5Group /Wavefunctions/PureSingleDeterminant does not exist. \n"; 
      return false;
    }

    if(!read.push("Wavefunctions")) return false;
    if(!read.push("PureSingleDeterminant")) return false;

    std::vector<int> Idata;
    Idata.reserve(9);
    Idata.resize(8);
    if(!read.read(Idata,"dims")) return false;

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

    std::vector<ComplexType> vvec0;
    if(read.read(vvec0,"Wavefun")) {
      int nc_ = NAEA;
      if(wfn_type==2) nc_ = NAEA+NAEB;
      if(vvec0.size() != 2*NMO*nc_) return false;
      OrbMat.resize(2*NMO,nc_);
      std::copy(vvec0.begin(),vvec0.end(),OrbMat.begin());
      Idata.push_back(1);
    } else {
      Idata.push_back(0);
      rotated_hamiltonian = false;
    }

    myComm->bcast(Idata);   
    if(rotated_hamiltonian) 
      myComm->bcast(vvec0);  

    if(rotated_hamiltonian) 
      haj.resize(Idata[0]);
    else
      hij.resize(Idata[0]);

    if(rotated_hamiltonian) {
      SMSpHabkl.setDims(Idata[2],Idata[3]);
      if(!SMSpHabkl.allocate_serial(Idata[1])) return false;
      SMSpHabkl.resize_serial(Idata[1]);
    } else {
      SMSpHijkl.setDims(Idata[2],Idata[3]);
      if(!SMSpHijkl.allocate_serial(Idata[1])) return false;
      SMSpHijkl.resize_serial(Idata[1]);
    }


    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);
    Idata.resize(NAEA+NAEB);
    if(!read.read(Idata,"occups")) return false;
    for(int i=0; i<NAEA; i++) occup_alpha[i] = Idata[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) occup_beta[j] = Idata[i];
    myComm->bcast(Idata);   


    std::vector<int> ivec;
    if(rotated_hamiltonian) {
     ivec.resize(haj.size());
     if(!read.read(ivec,"hij_indx")) return false;
     for(int i=0; i<haj.size(); i++)
       std::get<0>(haj[i]) = ivec[i]; 
     myComm->bcast(ivec);
    } else {
     ivec.resize(hij.size());
     if(!read.read(ivec,"hij_indx")) return false;
     for(int i=0; i<hij.size(); i++)
       std::get<0>(hij[i]) = ivec[i]; 
     myComm->bcast(ivec);
    }

    if(rotated_hamiltonian) {
      std::vector<ComplexType> vvec;
      vvec.resize(haj.size());
      if(!read.read(vvec,"hij")) return false;
      for(int i=0; i<haj.size(); i++)
        std::get<1>(haj[i]) = vvec[i]; 
      myComm->bcast(vvec);

    } else {
      std::vector<ValueType> vvec;
      vvec.resize(hij.size());
      if(!read.read(vvec,"hij")) return false;
      for(int i=0; i<hij.size(); i++)
        std::get<1>(hij[i]) = vvec[i]; 
      myComm->bcast(vvec);
    }

    if(rotated_hamiltonian) {
     if(!read.read(*(SMSpHabkl.getVals()),"SpHijkl_vals")) return false;
     if(!read.read(*(SMSpHabkl.getCols()),"SpHijkl_cols")) return false;
     if(!read.read(*(SMSpHabkl.getRowIndex()),"SpHijkl_rowIndex")) return false;
     SMSpHabkl.setRowsFromRowIndex();
    } else {
     if(!read.read(*(SMSpHijkl.getVals()),"SpHijkl_vals")) return false;
     if(!read.read(*(SMSpHijkl.getCols()),"SpHijkl_cols")) return false;
     if(!read.read(*(SMSpHijkl.getRowIndex()),"SpHijkl_rowIndex")) return false;
     SMSpHijkl.setRowsFromRowIndex();
    }

    myComm->barrier();
    myComm->barrier();

    read.pop();
    read.pop();

  } else {

    std::vector<int> Idata(9);
    myComm->bcast(Idata);

    rotated_hamiltonian = (Idata[8]!=0);
    if(rotated_hamiltonian) 
    {
      std::vector<ComplexType> vvec(2*NMO*((wfn_type==2)?(NAEA+NAEB):NAEA));
      myComm->bcast(vvec);
      OrbMat.resize(2*NMO,(wfn_type==2)?(NAEA+NAEB):NAEA);
      std::copy(vvec.begin(),vvec.end(),OrbMat.begin());
    }

    if(rotated_hamiltonian) {
     haj.resize(Idata[0]);
     SMSpHabkl.setDims(Idata[2],Idata[3]);
    } else { 
     hij.resize(Idata[0]);
     SMSpHijkl.setDims(Idata[2],Idata[3]);
    }
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
    if(rotated_hamiltonian) {
     ivec.resize(haj.size());
     myComm->bcast(ivec);
     for(int i=0; i<haj.size(); i++)
       std::get<0>(haj[i]) = ivec[i];
    } else {
     ivec.resize(hij.size());
     myComm->bcast(ivec);
     for(int i=0; i<hij.size(); i++)
       std::get<0>(hij[i]) = ivec[i];
    }

    if(rotated_hamiltonian) {
     std::vector<ComplexType> vvec;
     vvec.resize(haj.size());
     myComm->bcast(vvec);
     for(int i=0; i<haj.size(); i++)
       std::get<1>(haj[i]) = vvec[i];
    } else {
     std::vector<ValueType> vvec;
     vvec.resize(hij.size());
     myComm->bcast(vvec);
     for(int i=0; i<hij.size(); i++)
       std::get<1>(hij[i]) = vvec[i];
    }

    {
      std::vector<ComplexType> vvec(2*NMO*((wfn_type==2)?(NAEA+NAEB):NAEA));
      myComm->bcast(vvec);
      OrbMat.resize(2*NMO,(wfn_type==2)?(NAEA+NAEB):NAEA);
      std::copy(vvec.begin(),vvec.end(),OrbMat.begin());
    } 

    myComm->barrier();
    if(rotated_hamiltonian) {
     if(!SMSpHabkl.initializeChildren()) return false;
    } else {
     if(!SMSpHijkl.initializeChildren()) return false;
    }
    myComm->barrier();

  }
  readHamFromFile=true;
  return setup_local();
}

// approximate the k-set balance partition problem for ordered sets. k0/kN are the partition of the local core
void PureSingleDeterminant::split_Ham_rows(IndexType N, SPValueSMSpMat::int_iterator indx, IndexType& k0, IndexType& kN)
{
  if(ncores_per_TG == 1) {
    k0=0;
    kN=N;
    return;
  } 

  std::vector<SPValueSMSpMat::intType> subsets(ncores_per_TG+1); 
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

// this routine builds mixed_density_matrix
// temp_density_matrix is only used temporary
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
  DenseMatrixOperators::product(NMO,NAEA,NAEA,one,SlaterMat,NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);

  if(rotated_hamiltonian && full) {

#if defined(AFQMC_SP)
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,temp_density_matrix.data(),NMO);
    std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, mixed_density_matrix.begin());
#else
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,mixed_density_matrix.data(),NMO);
#endif
    DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data(),NMO);

  } else {
    std::fill(mixed_density_matrix.begin(), mixed_density_matrix.begin()+NMO*NMO, 0);
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
  DenseMatrixOperators::product(NMO,NAEB,NAEB,one,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

  if(rotated_hamiltonian && full) {

    // G(beta) = SS0*transpose(conjg(Dbeta)) 
#if defined(AFQMC_SP)
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,temp_density_matrix.data(),NMO);
    std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, mixed_density_matrix.begin()+NMO*NMO);
#else
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,mixed_density_matrix.data()+NMO*NMO,NMO);
#endif
    DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data()+NMO*NMO,NMO);

  } else {

    std::fill(mixed_density_matrix.begin()+NMO*NMO, mixed_density_matrix.end(), 0);
    // copy to mixed_density_matrix 
    for(int i=0; i<NAEB; i++)
     for(int j=0; j<NMO; j++)
      mixed_density_matrix(occup_beta[i],j) = SS0(j+NMO,i);

  }

} 

void PureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, SPComplexType* dm, bool full)
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
  DenseMatrixOperators::product(NMO,NAEA,NAEA,one,SlaterMat,NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);


  if(rotated_hamiltonian && full) {

#if defined(AFQMC_SP)
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,temp_density_matrix.data(),NMO);
    std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, dm); 
#else
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,dm,NMO);
#endif
    DenseMatrixOperators::transpose(NMO,dm,NMO);

  } else {

    // copy to mixed_density_matrix 
    for(int i=0; i<NAEA; i++)
     for(int j=0; j<NMO; j++)
      *(dm++) = SS0(j,i); 

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
  DenseMatrixOperators::product(NMO,NAEB,NAEB,one,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

  if(rotated_hamiltonian && full) {

#if defined(AFQMC_SP)
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,temp_density_matrix.data(),NMO);
    std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, dm+NMO*NMO);
#else
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,dm+NMO*NMO,NMO);
#endif    
    DenseMatrixOperators::transpose(NMO,dm+NMO*NMO,NMO);

  } else {

    // copy to mixed_density_matrix 
    for(int i=0; i<NAEB; i++)
     for(int j=0; j<NMO; j++)
      *(dm++) = SS0(j+NMO,i); 

  }
}
    
// on exit, commBuff contains the green function information for all walkers in appropriate place   
// buffer: (this routine assumes the buffer is large enough) 
//    wlksz: size of the data of a walker
//    gfoffset: location of the gf block within a walker
//    transposed: if true: walker index along columns
//                if false: walker index along rows (continuous walker data)
void PureSingleDeterminant::evaluateOneBodyMixedDensityMatrix(WalkerHandlerBase* wset, SPComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full)
{

  // nwalkers > ncores_per_TG: serial computation with round robin distribution
  // nwalkers < ncores_per_TG: parallel computation with round robin distribution

  ComplexType oa,ob;
  int nw = wset->numWalkers(true), cnt=0;
  int nw0 = wset->numWalkers(false);
// Careful here, only returns full matrix if (rotated_hamiltonian&&full), otherwise only NAEA*NMO
// Also careful with closed_shell=true, only returns alpha sector
  int sz = 2 + ((full&&rotated_hamiltonian)?NMO*NMO:NAEA*NMO);
  if(!closed_shell) sz += (full&&rotated_hamiltonian)?NMO*NMO:NAEB*NMO; 
  int wstride = transposed?nw0:1;
  int ax = transposed?1:wlksz; 
  // in Buff:  {..., oa, ob, DM, ...} for each walker 
  for(int i=0; i<nw; i++) {
    if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;  
    if(cnt%ncores_per_TG == core_rank ) {
      local_evaluateOneBodyMixedDensityMatrix(wset->getSM(i),oa,ob,local_buff.data()+2,full);
      local_buff[0]=oa;
      local_buff[1]=ob;
// Do I really need a lock?
      {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buf->getMutex()));
        BLAS::copy(sz,local_buff.data(),1,buf->values() + ax*cnt + wstride*gfoffset ,wstride);  
      }
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
    DenseMatrixOperators::product(NMO,NAEA,NAEA,one,OrbMat.data(),NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);

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
    DenseMatrixOperators::product(NMO,NAEB,NAEB,one,OrbMat.data()+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

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
    SPComplexType one = SPComplexType(1.0,0.0);
    SPComplexType zero = SPComplexType(0.0,0.0);
    int nr1=1, nc1=2*NMO*NMO;
    ComplexType oa,ob,epot,ekin = zero;    
    SPComplexSMVector* buffer = TG.commBuff;
 
    TG.resize_buffer(nwglobal);
    nwglobal /= sz;

// FIX FIX FIX: For rotated hamiltonian, inconsistent convention for location and index of beta blocks between
// hamiltonian and green function!!!
    // calculate [oa,ob,G] for all walkers
    evaluateOneBodyMixedDensityMatrix(wset,buffer,sz,1,false,false);

    if(nwglobal > local_buff.size())
      local_buff.resize(nwglobal);

    if(core_rank==0) {
      SPComplexType *val = buffer->values();
      for(int i=0; i<nw0; i++,val+=sz) *val = 0;
    } 

    // synchronize
    TG.local_barrier(); 

    // add 1-body term, round-robin
    for(int wlk=0; wlk<nw0; wlk++) { 
      ekin=ComplexType(0.0,0.0);
      mixed_density_matrix = zero;
      SPComplexType* ptr = buffer->values() + wlk*sz + 3;
      SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
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
            ekin += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
        } else {
          s1Dit end1 = hij.end();
          for(s1Dit it = hij.begin(); it != end1; it++)
            ekin += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
        }       
        ekin+=NuclearCoulombEnergy;
      }

      if(rotated_hamiltonian) 
        SparseMatrixOperators::product_SpMatV(int(ikN-ik0),SMSpHabkl.cols(),one,SMSpHabkl.values() + pik0, SMSpHabkl.column_data() + pik0, SMSpHabkl.row_index() + ik0,  mixed_density_matrix.data(),zero,V0.data()+ik0);
      else
        SparseMatrixOperators::product_SpMatV(int(ikN-ik0),SMSpHijkl.cols(),SPValueType(1),SMSpHijkl.values() + pik0, SMSpHijkl.column_data() + pik0, SMSpHijkl.row_index() + ik0,  mixed_density_matrix.data(),SPValueType(0),V0.data()+ik0);
      itG = mixed_density_matrix.begin()+ik0;
      SPComplexVector::iterator itV = V0.begin()+ik0;
      epot = ComplexType(0,0);
      for(int i=ik0; i<ikN; i++,++itG,++itV) epot += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
      local_buff[wlk] = ekin + 0.5*epot;
    }

    {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buffer->getMutex()));
      BLAS::axpy (nw0, one, local_buff.data(), 1, buffer->values(), sz);
    }
  
    // synchronize
    TG.local_barrier(); 

    int currnw=nw0;
    for(int tgi = 0; tgi<nnodes_per_TG; tgi++) {

      if(tgi > 0) {
        // adds local contribution to 2-body 
        for(int wlk=0; wlk<currnw; wlk++) { 
          mixed_density_matrix = zero;
          SPComplexType* ptr = buffer->values() + wlk*sz + 3;
          for(int i=0; i<NAEA; i++, ptr+=NMO )  // can use std::copy here or axpy 
           std::copy( ptr, ptr+NMO, mixed_density_matrix[occup_alpha[i]]);
          if(!closed_shell) {
           for(int i=0; i<NAEB; i++, ptr+=NMO)
            std::copy( ptr, ptr+NMO, mixed_density_matrix[occup_beta[i]]);
          }

          if(rotated_hamiltonian)
            SparseMatrixOperators::product_SpMatV(int(ikN-ik0),SMSpHabkl.cols(),one,SMSpHabkl.values() + pik0, SMSpHabkl.column_data() + pik0, SMSpHabkl.row_index() + ik0,  mixed_density_matrix.data(),zero,V0.data()+ik0);
          else
            SparseMatrixOperators::product_SpMatV(int(ikN-ik0),SMSpHijkl.cols(),SPValueType(1),SMSpHijkl.values() + pik0, SMSpHijkl.column_data() + pik0, SMSpHijkl.row_index() + ik0,  mixed_density_matrix.data(),SPValueType(0),V0.data()+ik0);
          SPComplexMatrix::iterator itG = mixed_density_matrix.begin()+ik0;
          SPComplexVector::iterator itV = V0.begin()+ik0;
          epot = ComplexType(0,0);
          for(int i=ik0; i<ikN; i++,++itG,++itV) epot += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
          local_buff[wlk] = ekin + 0.5*epot;
          //local_buff[wlk] = ekin + 0.5*zdotu(int(ikN-ik0),mixed_density_matrix.data()+ik0,1,V0.data()+ik0,1);
          //local_buff[wlk] = 0; 
        }
        {
          boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buffer->getMutex()));
          BLAS::axpy (currnw, one, local_buff.data(), 1, buffer->values(), sz);
        }
      }

      // rotate buffer 
      if(distribute_Ham && nnodes_per_TG>1) TG.rotate_buffer(currnw,sz);
    }

    // synchronize: rotate_buffer calls barrier, so no need to do it again 
    if(!distribute_Ham) TG.local_barrier(); 

    if(core_rank != 0) return;
    
    int nw = wset->numWalkers(true);
    SPComplexType* ptr = buffer->values();
    for(int i=0,cnt=0; i<nw; i++) {
      if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
      ComplexType eloc=*ptr;
      ComplexType oa=*(ptr+1);
      ComplexType ob=*(ptr+2);
      if(first)
        wset->setWalker(i,eloc,oa,ob);
      else {
        wset->setEloc2(i,eloc);
        wset->setOvlp2(i,oa,ob);
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
    SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
    if(rotated_hamiltonian) {
      std::vector<s1D<ComplexType> >::iterator  end1 = haj.end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj.begin(); it != end1; it++) 
       ekin += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
    } else {
      s1Dit end1 = hij.end();
      for(s1Dit it = hij.begin(); it != end1; it++) { 
        ekin += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
      }  
    }

    epot = 0; 
    int nr1=1, nc1=2*NMO*NMO;
    SPComplexType one = SPComplexType(1.0,0.0);
    SPComplexType zero = SPComplexType(0.0,0.0);
    if(rotated_hamiltonian) 
      SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHabkl.values(),SMSpHabkl.column_data(),SMSpHabkl.row_index(),mixed_density_matrix.data(),zero,V0.data()); 
    else
      SparseMatrixOperators::product_SpMatV(nc1,nc1,SPValueType(1),SMSpHijkl.values(),SMSpHijkl.column_data(),SMSpHijkl.row_index(),mixed_density_matrix.data(),SPValueType(0),V0.data()); 
    itG = mixed_density_matrix.begin();
    SPComplexVector::iterator itV = V0.begin();
    for(int i=0; i<nc1; i++,++itG,++itV) epot += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG); 
    epot = 0.5*epot+NuclearCoulombEnergy;   

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:evaluateLocalEnergy"); 
#endif
    
  }

  // for PureSD, I half-rotate the transposed operator to the basis of PsiT.
  // It is then always in compressed form. 
  void PureSingleDeterminant::generateTransposedOneBodyOperator(bool addBetaBeta, SPValueSMVector& Dvn, SPComplexSMVector& DvnT) {
    if(!addBetaBeta)
        APP_ABORT("  Error: PureSingleDeterminant::generateTransposedOneBodyOperator not implemented with UHF integrals.");  
    
    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    int NMO2=NMO*NMO;
    int ncol = closed_shell?NAEA:(NAEA+NAEB);

    // resize matrix
    DvnT.setup(head_of_nodes,"DvnT",TG.getNodeCommLocal());
    DvnT.setDims(Dvn.cols(),ncol*NMO);
    DvnT.resize(Dvn.cols()*ncol*NMO);    

    if(!rotated_hamiltonian) { 
      OrbMat.resize(2*NMO,NAEA);
      std::fill(OrbMat.begin(), OrbMat.end(), zero);
      for(int i=0; i<NAEA; i++)
        OrbMat(occup_alpha[i],i) = one; 
      if(!closed_shell) 
        for(int i=0; i<NAEB; i++)
          OrbMat(occup_beta[i],i) = one; 
    }

    int npr,rnk;
    MPI_Comm_rank(TG.getNodeCommLocal(),&rnk);
    MPI_Comm_size(TG.getNodeCommLocal(),&npr);

    ComplexMatrix CVn(NMO,NMO); 
    ComplexMatrix CVnrot(ncol,NMO); 
    for(int i=0, iend=Dvn.cols(); i<iend; i++) {

      if(i%npr != rnk) continue;

      std::fill(CVn.begin(), CVn.end(), zero);
      // extract Cholesky vector i
#if !defined(QMC_COMPLEX) || defined(AFQMC_SP)
      for(int ik=0; ik<NMO2; ik++)
        CVn(ik) = static_cast<ComplexType>(*(Dvn.values()+ik*iend+i));
#else      
      zcopy(NMO2, Dvn.values()+i, Dvn.cols(), CVn.data(), 1);
#endif

      // rotate the matrix: CVnrot = OrbMat^H * CVn
      DenseMatrixOperators::product_AhB(NAEA,CVn.cols(),CVn.rows(), one, OrbMat.data(), OrbMat.cols(), CVn.data(), CVn.cols(), zero, CVnrot.data(), CVnrot.cols() );
      if(!closed_shell) 
        DenseMatrixOperators::product_AhB(NAEB,CVn.cols(),CVn.rows(), one, OrbMat.data()+NAEA*NMO, OrbMat.cols(), CVn.data(), CVn.cols(), zero, CVnrot.data()+NAEA+NMO, CVnrot.cols() );
   
      // insert it in DvnT
#if !defined(QMC_COMPLEX) || defined(AFQMC_SP)
      for(int ik=0, ikend=ncol*NMO; ik<ikend; ik++)
        *(DvnT.values()+i*DvnT.cols()+ik) = static_cast<SPComplexType>(CVnrot(ik));
#else      
       zcopy(ncol*NMO, CVnrot.data() , 1, DvnT.values()+i*DvnT.cols(),1);
#endif

    }
    MPI_Barrier(TG.getNodeCommLocal());

  }
 
  void PureSingleDeterminant::generateTransposedOneBodyOperator(bool addBetaBeta, SPValueSMSpMat& Spvn, SPComplexSMSpMat& SpvnT) {

    if(!addBetaBeta)
        APP_ABORT("  Error: PureSingleDeterminant::generateTransposedOneBodyOperator not implemented with UHF integrals.");

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    int NMO2=NMO*NMO;
    int ncol = closed_shell?NAEA:(NAEA+NAEB);

    // resize matrix
    SpvnT.setup(head_of_nodes,"SpvnT",TG.getNodeCommLocal());
    SpvnT.setDims(Spvn.cols(),ncol*NMO);

    int npr,rnk;
    MPI_Comm_rank(TG.getNodeCommLocal(),&rnk);
    MPI_Comm_size(TG.getNodeCommLocal(),&npr);

    ComplexMatrix CVn, CVnrot;

    std::size_t sz=0;
    SPValueSMSpMat::intPtr row = Spvn.row_data();   
    SPValueSMSpMat::intPtr row_index = Spvn.row_index();   
    SPValueSMSpMat::intPtr col = Spvn.column_data();   
    SPValueSMSpMat::pointer val = Spvn.values();   

    if(rotated_hamiltonian) {

      CVn.resize(NMO,NMO);
      CVnrot.resize(ncol,NMO);

      app_log()<<" Temporarily transposing Cholesky matrix. \n";
      Spvn.transpose(TG.getNodeCommLocal());

      col = Spvn.column_data();
      val = Spvn.values();
      row_index = Spvn.row_index();
      for(int i=0, iend=Spvn.rows(); i<iend; i++, row_index++) {

        int nterms = *(row_index+1) - (*row_index);
        if(nterms==0) continue;

        if(i%npr != rnk) {
          col+=nterms;
          val+=nterms;
          continue;
        }
 
        std::fill(CVn.begin(), CVn.end(), zero);
        // extract Cholesky vector i
        for(int nt=0; nt<nterms; nt++, col++, val++) 
          CVn(*col) = static_cast<ComplexType>(*val);

        // rotate the matrix: CVnrot = OrbMat^H * CVn
        DenseMatrixOperators::product_AhB(NAEA,CVn.cols(),CVn.rows(), one, OrbMat.data(), OrbMat.cols(), CVn.data(), CVn.cols(), zero, CVnrot.data(), CVnrot.cols() );
        if(!closed_shell) 
          DenseMatrixOperators::product_AhB(NAEB,CVn.cols(),CVn.rows(), one, OrbMat.data()+NAEA*NMO, OrbMat.cols(), CVn.data(), CVn.cols(), zero, CVnrot.data()+NAEA*NMO, CVnrot.cols() );
 
      // insert it in DvnT
       for(int ik=0, ikend=ncol*NMO; ik<ikend; ik++)
        if(std::abs(CVnrot(ik)) > 1e-6)  // fixed cutoff for now
          sz++;
      }
    
    } else {

      int nterms = Spvn.size();
      int ntpc = nterms/npr, nex = nterms%npr;
      int p0 = rnk*ntpc + std::min(rnk,nex);
      int p1 = p0 + ntpc + ((rnk<nex)?1:0);
       
      col = Spvn.column_data()+p0;
      val = Spvn.values()+p0;
      row = Spvn.row_data()+p0;
      SPValueSMSpMat::pointer vend = Spvn.values()+p1;
      for(; val!= vend; val++, col++, row++) {
        IndexType i = (*row)/NMO;
        IndexType k = (*row)%NMO;
        if(isOcc_alpha[i])
          sz++;
        if(!closed_shell && isOcc_beta[i+NMO])
          sz++;
      }

    }

    std::size_t sz_=sz;
    MPI_Allreduce(&sz_,&sz,1,MPI_UNSIGNED_LONG,MPI_SUM,TG.getNodeCommLocal());
    SpvnT.allocate(sz);    

    using itype = SPComplexSMSpMat::intType;
    using vtype = SPComplexSMSpMat::value_type;

    int nmax = 10000;
    std::vector< std::tuple<itype,itype,vtype> > ints;
    ints.reserve(nmax); 

    if(rotated_hamiltonian) {

      col = Spvn.column_data();   
      val = Spvn.values();   
      row_index = Spvn.row_index();   
      for(int i=0, iend=Spvn.rows(); i<iend; i++, row_index++) {

        int nterms = *(row_index+1) - (*row_index);
        if(nterms==0) continue;

        if(i%npr != rnk) {
          col+=nterms;
          val+=nterms;
          continue;
        }

        std::fill(CVn.begin(), CVn.end(), zero);
        // extract Cholesky vector i
        for(int nt=0; nt<nterms; nt++, col++, val++) 
          CVn(*col) = static_cast<ComplexType>(*val);

        // rotate the matrix: CVnrot = OrbMat^H * CVn
        DenseMatrixOperators::product_AhB(NAEA,CVn.cols(),CVn.rows(), one, OrbMat.data(), OrbMat.cols(), CVn.data(), CVn.cols(), zero, CVnrot.data(), CVnrot.cols() );
        if(!closed_shell) 
          DenseMatrixOperators::product_AhB(NAEB,CVn.cols(),CVn.rows(), one, OrbMat.data()+NAEA*NMO, OrbMat.cols(), CVn.data(), CVn.cols(), zero, CVnrot.data()+NAEA*NMO, CVnrot.cols() );
   
        // insert it in DvnT
        for(int ik=0, ikend=ncol*NMO; ik<ikend; ik++)
          if(std::abs(CVnrot(ik)) > 1e-6)  {// fixed cutoff for now
            ints.push_back(std::make_tuple(i,ik,static_cast<SPComplexType>(CVnrot(ik))));
            if(ints.size() == nmax) {
              SpvnT.add(ints,true);
              ints.clear();
            }
          }
      }
      if(ints.size() > 0) {
        SpvnT.add(ints,true);
        ints.clear();
      }

    } else {

      int nterms = Spvn.size();
      int ntpc = nterms/npr, nex = nterms%npr;
      int p0 = rnk*ntpc + std::min(rnk,nex);
      int p1 = p0 + ntpc + ((rnk<nex)?1:0);

      std::vector<SPComplexSMSpMat::intType> mymap(2*NMO,-1);
      for(int i=0; i<NAEA; i++)
        mymap[occup_alpha[i]] = i;  
      for(int i=0; i<NAEB; i++)
        mymap[occup_beta[i]] = NAEA+i;  

      col = Spvn.column_data()+p0;
      val = Spvn.values()+p0;
      row = Spvn.row_data()+p0;
      SPValueSMSpMat::pointer vend = Spvn.values()+p1;
      for(; val!= vend; val++, col++, row++) {
        itype i = (*row)/NMO;
        itype k = (*row)%NMO;
        if(isOcc_alpha[i]) {
          assert(mymap[i] >= 0);
          itype ik = mymap[i]*NMO+k;
          ints.push_back(std::make_tuple(*col,ik,static_cast<SPComplexType>(*val)));
          if(ints.size() == nmax) {
            SpvnT.add(ints,true);
            ints.clear();
          }
        }
        if(!closed_shell && isOcc_beta[i+NMO]) {
          assert(mymap[i+NMO] >= 0);
          itype ik = mymap[i+NMO]*NMO+k;
          ints.push_back(std::make_tuple(*col,ik,static_cast<SPComplexType>(*val)));
          if(ints.size() == nmax) {
            SpvnT.add(ints,true);
            ints.clear();
          }
        }
      }
      if(ints.size() > 0) {
        SpvnT.add(ints,true);
        ints.clear();
      }

    }


    app_log()<<" Compressing transposed Cholesky matrix. \n";
    SpvnT.compress(TG.getNodeCommLocal());
    if(rotated_hamiltonian) {
      app_log()<<" Transposing Cholesky matrix back to original form. \n";
      Spvn.transpose(TG.getNodeCommLocal());
    }

  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMSpMat& vn, SPComplexSMSpMat& vnT, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif
    ComplexType o1,o2;
    const SPComplexType *GF=GG;
    if(needsG) {
      GF = mixed_density_matrix.data();
      if(transposed) // transposed implies compact storage 
       local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,mixed_density_matrix.data(),false);
      else 
       local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);
    }

    if(transposed) {
      SPComplexType one = SPComplexType(1.0);
      const SPComplexType zero = SPComplexType(0.0);
      if(closed_shell) one = SPComplexType(2.0);      
      // since transposed implies compact storage, both components are done in one call
      // if closed_shell, only 1 component is there and factor of 2 is accounted for through "one"
      SparseMatrixOperators::product_SpMatV(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.column_data(),vnT.row_index(),GF,zero,v.data());
    } else {
      SPValueType one = SPValueType(1.0);
      const SPValueType zero = SPValueType(0.0);
      if(closed_shell) one = SPValueType(2.0);      
      SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),GF+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMVector& vn, SPComplexSMVector& vnT, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif
    ComplexType o1,o2;
    const SPComplexType *GF=GG;
    if(needsG) {
      GF = mixed_density_matrix.data();
      if(transposed) // transposed implies compact storage 
       local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,mixed_density_matrix.data(),false);
      else 
       local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);
    }

    if(transposed) {
      SPComplexType one = SPComplexType(1.0);
      const SPComplexType zero = SPComplexType(0.0);
      if(closed_shell) one = SPComplexType(2.0);
      DenseMatrixOperators::product_Ax(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.cols(),GF,zero,v.data());
    } else {
      SPValueType one = SPValueType(1.0);
      const SPValueType zero = SPValueType(0.0);
      if(closed_shell) one = SPValueType(2.0);
      DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one,vn.values(),vn.cols(),GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one,vn.values(),vn.cols(),GF+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  // buf is a matrix of walker information which includes the green's function. This matrix
  // was created above by "". The number of columns is walkerBlock and the number of rows is the size of G (since walkers go along columns, the rest of the buffer is irrelevant)
  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int i0, int iN, SPValueSMSpMat& vn, SPComplexSMSpMat& vnT, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n)
  {
    if(transposed)
      assert(v.size() >= nW*(iN-i0));
    else
      assert(v.size() >= nW*vn.cols());
    const SPComplexType *GF = buff+2*nW;
    int GFi0 = i0;
    int GFi02 = NMO*NMO+i0;
    SPValueType one = SPValueType(1.0);
    const SPValueType zero = SPValueType(0.0);
    if(closed_shell) one = SPValueType(2.0);      
    const SPComplexType czero = SPComplexType(0.0,0.0);
    SPComplexType cone = SPComplexType(1.0);
    if(closed_shell) cone = SPComplexType(2.0);      

    Timer.start("PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer::setup");
    if(!rotated_hamiltonian && !transposed) {  
      // check if there is any work to be done
      // for spinRestricted (addBetaBeta==true)
      //   1: i0/NMO <= occup_alpha[NAEA-1]   
      //   2: iN/NMO >= occup_alpha[0]
      // for addBetaBeta==false
      //   1: i0/NMO <= occup_beta[NAEB-1]   
      //   2: iN/NMO >= occup_beta[0]
      // Notice that if addBetaBeta==true, then i0/iN are within 0:NMO*NMO, otherwise they are from 0:2*NMO*NMO 
      //
      int j0 = i0/NMO;
      int l0 = i0%NMO;
      int jN = iN/NMO;
      int lN = iN%NMO;

      // I have 3 scenarios: 1: addBetaBeta==false, 2: addBetaBeta==true && closed_shell==true, 3: addBetaBeta==true && closed_shell==false
      int mult=1; 
      bool ovlpA = (j0 <= occup_alpha[NAEA-1]) && (jN >= occup_alpha[0]);  
      bool ovlpB; // assume no overlap 
      if(addBetaBeta) {   // {j0,l0,jN,lN} within 0:NMO*NMO
        if(closed_shell)  // no need to build beta sector
          ovlpB=false;  
        else {  // build beta sector and define GFi02
          mult=2;  
          // in this case, must compare shifted indexes 
          ovlpB = (j0 <= occup_beta[NAEB-1]-NMO) && (jN >= occup_beta[0]-NMO);
        }
      } else  // {j0,l0,jN,lN} within 0:2*NMO*NMO
        ovlpB = (j0 <= occup_beta[NAEB-1]) && (jN >= occup_beta[0]);

      if( (!ovlpA && !ovlpB) ) { //nothing to do
        std::fill(v.begin(),v.begin()+nW*vn.cols(),czero);  // v might be bigger, zero out only expected region 
        return;
      }

      // slightly larger resize to make routine simpler
      cGF.resize(nW*mult*(jN-j0+1)*NMO);
      // no need to generate GF matrix outside {i0:iN}
      std::fill(cGF.begin(),cGF.end(),czero); 
      // copy to cGF 
      buff+=2*nW; 
      // copy entire orbital that maps to [j0,jN], shifted back by j0*NMO
      // GFi0 is the location of i0 in the reduced matrix 
      GFi0=l0;
      // GFi02 is the location of i0+NMO*NMO in the reduced matrix
      GFi02=nW*(jN-j0+1)*NMO+l0;

      for(int i=0; i<NAEA; i++, buff+=nW*NMO) { 
        if(occup_alpha[i] >= j0 && occup_alpha[i] <= jN)  
          std::copy(buff,buff+nW*NMO,cGF.begin()+(occup_alpha[i]-j0)*NMO*nW); 
      }

      if(!closed_shell) {  
        if(addBetaBeta) {
          for(int i=0; i<NAEB; i++,buff+=nW*NMO)
            if(occup_beta[i]-NMO >= j0 && occup_beta[i]-NMO <= jN)
              std::copy(buff,buff+nW*NMO,cGF.begin()+nW*(jN-j0+1)*NMO+(occup_beta[i]-NMO-j0)*NMO*nW);
        } else {
          for(int i=0; i<NAEB; i++,buff+=nW*NMO)
            if(occup_beta[i] >= j0 && occup_beta[i] <= jN)
              std::copy(buff,buff+nW*NMO,cGF.begin()+(occup_beta[i]-j0)*NMO*nW);
        }
      }
      GF = cGF.data();

    }
    Timer.stop("PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer::setup");

    std::size_t p_;
    if(transposed)
      p_ = static_cast<std::size_t>(*(vnT.row_index()+i0));
    else 
      p_ = static_cast<std::size_t>(*(vn.row_index()+i0));
    // walkerBlock implies MatV
    if(walkerBlock==1) {

      if(transposed) {      
        // transposed implies compact DM
        SparseMatrixOperators::product_SpMatV( int(iN-i0), vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF, czero, v.data());
      } else {
        SparseMatrixOperators::product_SpMatTV( int(iN-i0), vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+GFi0, zero, v.data());
        if(addBetaBeta && !closed_shell)
          SparseMatrixOperators::product_SpMatTV( int(iN-i0), vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+GFi02, one, v.data());
      }

    } else {
    // walkerBlock>1, so use MatM  
    // almost certainly, it is highly favourable  to set walkerBlock to nW
   
      int nblk = nW/walkerBlock;
      int nextra = nW%walkerBlock;      

      for(int i=0; i<nblk; i++) {

        if(transposed) {
          // transposed implies compact DM
          SparseMatrixOperators::product_SpMatM( int(iN-i0), walkerBlock , vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+i*walkerBlock, nW, czero, v.data()+i*walkerBlock, nW);
        } else {
          SparseMatrixOperators::product_SpMatTM( int(iN-i0), walkerBlock , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+GFi0*nW+i*walkerBlock, nW, zero, v.data()+i*walkerBlock, nW);
          if(addBetaBeta && !closed_shell)
            SparseMatrixOperators::product_SpMatTM( int(iN-i0), walkerBlock , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+GFi02*nW+i*walkerBlock, nW, one, v.data()+i*walkerBlock, nW);
        } 
 
      }
      // perform remaining walkers. Should I pad this to the closest power of 2???
      if(nextra > 0) {
        if(transposed) {
          // transposed implies compact DM
          SparseMatrixOperators::product_SpMatM( int(iN-i0), nextra , vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+nblk*walkerBlock, nW, czero, v.data()+nblk*walkerBlock, nW);
        } else {
          SparseMatrixOperators::product_SpMatTM( int(iN-i0), nextra , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+GFi0*nW+nblk*walkerBlock, nW, zero, v.data()+nblk*walkerBlock, nW);
          if(addBetaBeta && !closed_shell)
            SparseMatrixOperators::product_SpMatTM( int(iN-i0), nextra , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+GFi02*nW+nblk*walkerBlock, nW, one, v.data()+nblk*walkerBlock, nW);
        }
      }           
    }

  }

  // buf is a matrix of walker information which includes the green's function. This matrix
  // was created above by "". The number of columns is walkerBlock and the number of rows is the size of G (since walkers go along columns, the rest of the buffer is irrelevant)
  void PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int i0, int iN, SPValueSMVector& vn, SPComplexSMVector& vnT, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n)
  {
    if(transposed)
      assert(v.size() >= nW*(iN-i0));
    else
      assert(v.size() >= nW*vn.cols());
    const SPComplexType *GF = buff+2*nW;
    int GFi0 = i0;
    int GFi02 = NMO*NMO+i0;
    SPValueType one = SPValueType(1.0);
    const SPValueType zero = SPValueType(0.0);
    if(closed_shell) one = SPValueType(2.0);      
    const SPComplexType czero = SPComplexType(0.0,0.0);
    SPComplexType cone = SPComplexType(1.0);
    if(closed_shell) cone = SPComplexType(2.0);      
    Timer.start("PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer::setup");
    if(!rotated_hamiltonian && !transposed) {  
      // check if there is any work to be done
      // for spinRestricted (addBetaBeta==true)
      //   1: i0/NMO <= occup_alpha[NAEA-1]   
      //   2: iN/NMO >= occup_alpha[0]
      // for addBetaBeta==false
      //   1: i0/NMO <= occup_beta[NAEB-1]   
      //   2: iN/NMO >= occup_beta[0]
      // Notice that if addBetaBeta==true, then i0/iN are within 0:NMO*NMO, otherwise they are from 0:2*NMO*NMO 
      //
      int j0 = i0/NMO;
      int l0 = i0%NMO;
      int jN = iN/NMO;
      int lN = iN%NMO;

      // I have 3 scenarios: 1: addBetaBeta==false, 2: addBetaBeta==true && closed_shell==true, 3: addBetaBeta==true && closed_shell==false
      int mult=1; 
      bool ovlpA = (j0 <= occup_alpha[NAEA-1]) && (jN >= occup_alpha[0]);  
      bool ovlpB; // assume no overlap 
      if(addBetaBeta) {   // {j0,l0,jN,lN} within 0:NMO*NMO
        if(closed_shell)  // no need to build beta sector
          ovlpB=false;  
        else {  // build beta sector and define GFi02
          mult=2;  
          // in this case, must compare shifted indexes 
          ovlpB = (j0 <= occup_beta[NAEB-1]-NMO) && (jN >= occup_beta[0]-NMO);
        }
      } else  // {j0,l0,jN,lN} within 0:2*NMO*NMO
        ovlpB = (j0 <= occup_beta[NAEB-1]) && (jN >= occup_beta[0]);

      if( (!ovlpA && !ovlpB) ) { //nothing to do
        std::fill(v.begin(),v.begin()+nW*vn.cols(),czero);  // v might be bigger, zero out only expected region 
        return;
      }

      // slightly larger resize to make routine simpler
      cGF.resize(nW*mult*(jN-j0+1)*NMO);
      // no need to generate GF matrix outside {i0:iN}
      std::fill(cGF.begin(),cGF.end(),czero); 
      // copy to cGF 
      buff+=2*nW; 
      // copy entire orbital that maps to [j0,jN], shifted back by j0*NMO
      // GFi0 is the location of i0 in the reduced matrix 
      GFi0=l0;
      // GFi02 is the location of i0+NMO*NMO in the reduced matrix
      GFi02=nW*(jN-j0+1)*NMO+l0;
      for(int i=0; i<NAEA; i++, buff+=nW*NMO) { 
        if(occup_alpha[i] >= j0 && occup_alpha[i] <= jN)  
          std::copy(buff,buff+nW*NMO,cGF.begin()+(occup_alpha[i]-j0)*NMO*nW); 
      }
      if(!closed_shell) {  
        if(addBetaBeta) {
          for(int i=0; i<NAEB; i++,buff+=nW*NMO)
            if(occup_beta[i]-NMO >= j0 && occup_beta[i]-NMO <= jN)
              std::copy(buff,buff+nW*NMO,cGF.begin()+nW*(jN-j0+1)*NMO+(occup_beta[i]-NMO-j0)*NMO*nW);
        } else {
          for(int i=0; i<NAEB; i++,buff+=nW*NMO)
            if(occup_beta[i] >= j0 && occup_beta[i] <= jN)
              std::copy(buff,buff+nW*NMO,cGF.begin()+(occup_beta[i]-j0)*NMO*nW);
        }
      }
      GF = cGF.data();
    }
    // if transposed matrix is used, distribution is done over cholesky vectors.
    Timer.stop("PureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer::setup");

    // walkerBlock implies MatV
    if(walkerBlock==1) {
      if(transposed) {
        // transposed implies compact DM
        DenseMatrixOperators::product_Ax( int(iN-i0), vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(), GF, czero, v.data());
      } else {
        DenseMatrixOperators::product_Atx( int(iN-i0), vn.cols(), one, vn.values() + i0*vn.cols(), vn.cols(), GF+GFi0, zero, v.data());
        if(addBetaBeta && !closed_shell)
          DenseMatrixOperators::product_Atx( int(iN-i0), vn.cols(), one, vn.values() + i0*vn.cols(), vn.cols(), GF+GFi02, one, v.data());
      }
    } else {
    // walkerBlock>1, so use MatM  
    // almost certainly, it is highly favourable  to set walkerBlock to nW
   
      int nblk = nW/walkerBlock;
      int nextra = nW%walkerBlock;      

      for(int i=0; i<nblk; i++) {
        if(transposed) {
          // transposed implies compact DM
          DenseMatrixOperators::product( int(iN-i0), walkerBlock, vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(), 
               GF+i*walkerBlock,nW,czero,v.data()+i*walkerBlock,nW);
        } else {
          DenseMatrixOperators::product_AtB( vn.cols(), walkerBlock, int(iN-i0), one, vn.values() + i0*vn.cols(), vn.cols(), 
               GF+GFi0*nW+i*walkerBlock,nW,zero,v.data()+i*walkerBlock,nW);
          if(addBetaBeta && !closed_shell)
            DenseMatrixOperators::product_AtB( vn.cols(), walkerBlock, int(iN-i0), one, vn.values() + i0*vn.cols(), vn.cols(), 
               GF+GFi02*nW+i*walkerBlock,nW,one,v.data()+i*walkerBlock,nW);
        } 
 
      }
      // perform remaining walkers. Should I pad this to the closest power of 2???
      if(nextra > 0) {
        if(transposed) {
          // transposed implies compact DM
          DenseMatrixOperators::product( int(iN-i0), nextra, vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(), 
               GF+nblk*walkerBlock,nW,czero,v.data()+nblk*walkerBlock,nW);
        } else {
          DenseMatrixOperators::product_AtB( vn.cols(), nextra, int(iN-i0), one, vn.values() + i0*vn.cols(), vn.cols(), 
               GF+GFi0*nW+nblk*walkerBlock,nW,zero,v.data()+nblk*walkerBlock,nW);
          if(addBetaBeta && !closed_shell)
            DenseMatrixOperators::product_AtB( vn.cols(), nextra, int(iN-i0), one, vn.values() + i0*vn.cols(), vn.cols(),
               GF+GFi02*nW+nblk*walkerBlock,nW,one,v.data()+nblk*walkerBlock,nW);
        }
      }           
    }
  }

  void PureSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<SPComplexType> >& v2n, const std::vector<IndexType>& v2n_indx, SPValueSMSpMat& vn, std::vector<SPComplexType>& v, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("PureSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators");
#endif
    ComplexType o1,o2;
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2,true);


    SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
    for(int i=0; i<v2n_indx.size()-1; i++) {
      v[i] = static_cast<SPComplexType>(0.0);
      for(int n = v2n_indx[i]; n<v2n_indx[i+1]; n++) {
        IndexType ik = Index2Mat(std::get<0>(v2n[n]),std::get<2>(v2n[n]));
        IndexType jl = Index2Mat(std::get<1>(v2n[n]),std::get<3>(v2n[n]));
        IndexType il = Index2Mat(std::get<0>(v2n[n]),std::get<3>(v2n[n]));
        IndexType jk = Index2Mat(std::get<1>(v2n[n]),std::get<2>(v2n[n]));
        v[i] += (*(itG + ik)) * (*(itG + jl) - *(itG + il)) * (*(itG + jk)) * std::get<4>(v2n[n]);
      }
    }

#ifdef AFQMC_TIMER
    Timer.stop("PureSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators");
#endif

  }

  void PureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMVector& vn, std::vector<SPComplexType>& v, const int n)
  {

    if(trialDensityMatrix_needsupdate) {
      trialDensityMatrix_needsupdate = false;
      ComplexType o1,o2;
      local_evaluateOneBodyTrialDensityMatrix(true);
      SPtrial_density_matrix = trial_density_matrix; 
    }


    SPValueType one = SPValueType(1.0);
    const SPValueType zero = SPValueType(0.0);
    if(closed_shell) one = SPValueType(2.0);
    DenseMatrixOperators::product_Atx(vn.rows(), vn.cols(), one, vn.values(), vn.cols(),SPtrial_density_matrix.data(), zero, v.data());
    if(addBetaBeta && !closed_shell)
    DenseMatrixOperators::product_Atx(vn.rows(), vn.cols(), one, vn.values(), vn.cols(),SPtrial_density_matrix.data()+NMO*NMO, one, v.data());
  }

  void PureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMSpMat& vn, std::vector<SPComplexType>& v, const int n)
  {

    if(trialDensityMatrix_needsupdate) {
      trialDensityMatrix_needsupdate = false;
      ComplexType o1,o2;
      local_evaluateOneBodyTrialDensityMatrix(true);
      SPtrial_density_matrix = trial_density_matrix; 
    }


    SPValueType one = SPValueType(1.0);
    const SPValueType zero = SPValueType(0.0);
    if(closed_shell) one = SPValueType(2.0);
    SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),SPtrial_density_matrix.data(),zero,v.data());
    if(addBetaBeta && !closed_shell)
      SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),SPtrial_density_matrix.data()+NMO*NMO,one,v.data());

  }

}

      
