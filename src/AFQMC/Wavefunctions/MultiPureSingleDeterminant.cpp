#include<cstdlib>
#include<complex>
#include<iostream>
#include<fstream>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "Message/CommOperators.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/SimpleParser.h"
#include "Configuration.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Hamiltonians/SparseGeneralHamiltonian.h"
#include "AFQMC/Wavefunctions/WavefunctionBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHelper.h"
#include "AFQMC/Wavefunctions/MultiPureSingleDeterminant.h"

#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

namespace qmcplusplus
{

bool MultiPureSingleDeterminant::parse(xmlNodePtr cur)
{
    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur;

    std::string type("");
    init_type=std::string("");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.add(type,"type");
    oAttrib.add(init_type,"init");
    oAttrib.put(cur);

    std::transform(type.begin(),type.end(),type.begin(),(int (*)(int)) tolower);
    std::transform(init_type.begin(),init_type.end(),init_type.begin(),(int (*)(int)) tolower);

    if(type != "multipuresd") {
      app_error()<<" ERROR: Problems in MultiPureSingleDeterminant::parse: type should be MultiPureSD. \n"  <<std::endl;
      return false;
    }

    IterCI_maxit=0; 
    IterCI_cut=0; 

    std::string str("no");
    std::string str1("yes");
    std::string str2("no");
    filename = std::string("none");
    filetype = std::string("ascii");
    ParameterSet m_param;
    m_param.add(filename,"filename","std::string");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(str,"save_memory","std::string");
    m_param.add(cutoff,"cutoff","double");
    m_param.add(str1,"diagHam","std::string");
    m_param.add(str2,"iterCI","std::string");
    m_param.add(IterCI_maxit,"iterCI_it","int");
    m_param.add(IterCI_maxit,"iterCI_maxit","int");
    m_param.add(IterCI_cut,"iterCI_cut","double");
    m_param.add(diag_in_steps,"diag_steps","int");
    m_param.put(cur);
  
    std::transform(str.begin(),str.end(),str.begin(),(int (*)(int)) tolower);
    if(str == "yes" || str == "true") runtype = 0; 
    else runtype=1;

    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int)) tolower);
    if(str1 == "no" || str1 == "false") diagHam = false; 

    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int)) tolower);
    if(str2 == "yes" || str2 == "true") iterCI = true; 

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }

    return true;
}

bool MultiPureSingleDeterminant::initFromAscii(std::string fileName)
{

  std::ifstream in;
  in.open(std::string(fileName.c_str()).c_str());
  if(in.fail()) {
     app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
     return false;
  }

  NCA = NCB = 0;

  // fixed for now!!!
  ref = 0;
  std::string format = "standard";

  int na,nb,nm,num;
  std::string type; 
  int nci;

/*
  // use readHeader later 
  in>>nci >>na >>nb >>nm >>num;
  if(in.fail()) {
     app_error()<<"Error with wfn file. Header!!!" <<std::endl; 
     return false;
  }    

  // check for consistency
  if( NAEA>0 )  
    if(NAEA != na) {
      app_error()<<" Error: NAEA incompatible between input file and MultiPureSingleDeterminant init file. \n";
      return false;
    }
  else 
    NAEA=na;

  if( NAEB>0 )
    if(NAEB != nb) {
      app_error()<<" Error: NAEB incompatible between input file and MultiPureSingleDeterminant init file. \n";
      return false;
    }
  else   
    NAEB=nb;

  if( NMO>0 )
    if(NMO != nm) {
      app_error()<<" Error: NMO incompatible between input file and MultiPureSingleDeterminant init file. \n";
      return false;
    }
  else   
    NMO=nm;

  // read format type 
  in>>type;
  if(type == "occ") {
    rotated_hamiltonian = false;
  } else if(type == "rotate") {
    rotated_hamiltonian = true;
  } else {
    app_error()<<" Error in MultiPureSingleDeterminant init file. Unknown file type: " <<type <<std::endl; 
    return false;
  }
*/

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
      } else if(*it == "NCI" || *it == "nci") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NETOT \n";
          return false;
        }
        nci = atoi((++it)->c_str());
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
            app_log()<<"Using a UHF-type trial wave-function in MultiPureSingleDeterminant. \n";
            break;
          }
          case 2:
          {
            app_log()<<"Using a GHF-type trial wave-function in MultiPureSingleDeterminant. \n";
            app_error()<<" GHF type not implemented. \n";
            return false;
            break;
          }
          default:
          {
            app_error()<<"Unknown wave-function type in MultiPureSingleDeterminant: " <<wfntype <<std::endl;
            return false;
          }
        }
      } else if(*it == "FullMO" || *it == "FULLMO") {
        fullMOMat = true;
        app_log()<<" Expecting full MO matrix in MultiPureSingleDeterminant.\n";
      } else if(*it == "FORMAT" || *it == "Format" || *it == "format") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. FORMAT \n";
          return false;
        }
        if( *(it+1) == "phf" || *(it+1) == "PHF" ) format = "phf"; 
        else if( *(it+1) == "standard" ) format = "standard";  
        else {
          app_error()<<" Unknown file format: " <<*(it+1) <<std::endl;
          return false;
        }
        it++;
      } else if(*it == "Type" || *it == "TYPE" || *it == "type") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. TYPE \n";
          return false;
        }
        if( *(it+1) == "rotated" || *(it+1) == "rotate" ) rotated_hamiltonian=true; 
        it++;
      } else if(*it == "CMajor") {
        Cstyle = false;
        app_log()<<" Expecting MO matrix in Column-major format in MultiPureSingleDeterminant.\n";
      }
    }
    getwords(words,in);
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
  } while((words[0].find(std::string("/"))==std::string::npos && words[0].find(std::string("&END"))==std::string::npos));

  ci.resize(nci);
  if(!rotated_hamiltonian) {
    app_log()<<" Requested Pure Slater Determinant version of MultiPureSlaterDeterminant. \n"; 
    if(runtype==0) app_log()<<" Using low memory algorithm (cpu cost ~ O[n*M^3]) " <<std::endl;
    else app_log()<<" Using high memory algorithm (cpu cost ~ O[n^2*M^2]) " <<std::endl;
    occ_orbs.resize(nci*(NAEA+NAEB));
    std::vector<IndexType>::iterator it_occ = occ_orbs.begin();
    for(int i=0; i<nci; i++) {

      if(filetype == "ascii") {
        // later on, apply cutoff
        in>>ci[i];
        if(in.fail()) {
           app_error()<<"Error with wfn file. CI #: " <<i <<std::endl;
           return false;
        }
      }

      for(int k=0; k<NAEA; k++,it_occ++) {
        in>>*it_occ;    
        (*it_occ)--;
        if(in.fail()) {
           app_error()<<"Error with wfn file. CI #: " <<i <<std::endl; 
           return false;
        }     
        if(*it_occ < 0 || *it_occ > NMO) {
           app_error()<<"Error with wfn file. Det definition (NAEA) # : "<<i <<" " <<*it_occ <<std::endl; 
           return false;
        } 
      }
      for(int k=0; k<NAEB; k++,it_occ++) {
        in>>*it_occ;
        if(filetype == "sqc_ascii") *it_occ += NMO;
        (*it_occ)--;
        if(in.fail()) {
           app_error()<<"Error with wfn file. CI #: " <<i <<std::endl;
           return false;
        }
        if(*it_occ < NMO || *it_occ >= 2*NMO) {
           app_error()<<"Error with wfn file. Det definition (NAEB) # : "<<i <<"  " <<*it_occ <<std::endl;
           return false;
        }
      }
      if(filetype == "sqc_ascii") {
        // later on, apply cutoff
        in>>ci[i];
        if(in.fail()) {
          app_error()<<"Error with wfn file. CI #: " <<i <<std::endl;
          return false;
        }
      }
    }
  } else {
    // read OrbMats
    app_log()<<" Requested Rotated Determinant version of MultiPureSlaterDeterminant. \n"; 
    runtype = 1;
   
// debug
//occ_orbs.resize(nci*(NAEA+NAEB));
//for(int i=0; i<NAEA; i++) occ_orbs[i]=i;
//for(int i=0; i<NAEB; i++) occ_orbs[i+NAEA]=i+NMO;

    if(format=="phf") {
      for(int i=0; i<nci; i++) in>>ci[i]; 
    }   

    std::string dum;
    int dumi;
    orbsize = NAEA*NMO;
    if(wfntype==1) orbsize*=2;
    if(wfntype==2) orbsize = NMO*(NAEA+NAEB); 
    OrbMat.resize(nci*orbsize); 
    ComplexType dummy;
    int nread;
    for(int ki=0; ki<nci; ki++) {

      if(format=="standard") in>>ci[ki];
      if(format=="phf") in>>dum >>dumi;    
      if(dum != "Determinant:") {
        app_error()<<" Format error in MultiPureSlaterDeterminant input file. Expecting Determinant: at  det# " <<ki <<", found: " <<dum <<std::endl; 
        return false;
      }

      if (wfntype == 0 ) {

       if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat[ki*orbsize+i*NAEA+j] = dummy;
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
            if(j<NAEA) OrbMat[ki*orbsize+i*NAEA+j] = dummy;
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
            if(j<NAEA) OrbMat[ki*orbsize+i*NAEA+j] = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
        nread = fullMOMat?NMO:NAEB;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEB) OrbMat[ki*orbsize+(i+NMO)*NAEA+j] = dummy;
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
            if(j<NAEA) OrbMat[ki*orbsize+i*NAEA+j] = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
        nread = fullMOMat?NMO:NAEB;
        for(int j=0; j<nread; j++) 
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEB) OrbMat[ki*orbsize+(i+NMO)*NAEA+j] = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant. \n";
              in.close();
              return false;
            }
          }
       }

      } else if(wfntype==2) {
        APP_ABORT("Error: wfntype==2 not uimplemented yet in MultiPureSD. \n\n\n");
      }
    } 
  }
  in.close();

  //trial_density_matrix.resize(NAEA+NAEB,NMO);
  trial_density_matrix.resize(2*NMO,NMO);
  rank_updated_trial_density_matrix.resize(2*NMO,NMO);

  //mixed_density_matrix.resize(NAEA+NAEB,NMO);
  full_mixed_density_matrix.resize(2*NMO,NMO);
  mixed_density_matrix.resize(2*NMO,NMO);
  rank_updated_mixed_density_matrix.resize(2*NMO,NMO);

  StoreSM.resize(2*NMO,NAEA);
  
  overlaps.resize(2*ci.size());
  pe.resize(ci.size());
  ke.resize(ci.size());
  trialDensityMatrix_needsupdate = true;

  // temporary storage
  S0.resize(NAEA,NAEA); 
  S1.resize(NAEB,NAEB); 
  SS0.resize(2*NMO,NAEA); 
  V0.resize(2*NMO*NMO);

  // larger than needed, but just in case
  Cwork.resize(2*NMO);
  pivot.resize(2*NMO);

  Iwork.resize(NAEA+NAEB);

  // calculate relations between determinants
  // setup necessary data structures for low-rank updates

  return true;
}

bool MultiPureSingleDeterminant::getHamiltonian(HamPtr h)
{
  ham0 = h;
  sHam = dynamic_cast<SparseGeneralHamiltonian*>(h);
  if(!sHam) { 
    app_error()<<" Error in MultiPureSingleDeterminant::getHamiltonian. \n"
               <<" Hamiltonian associated with MultiPureSingleDeterminant must of the \n"
               <<" type SparseGeneralHamiltonian. \n";
    APP_ABORT("");
  }
 
  NuclearCoulombEnergy = static_cast<ValueType>(sHam->NuclearCoulombEnergy);

  spinRestricted = sHam->RHF();
  closed_shell = false;
  if(NCA==NCB && NAEA == NAEB && rotated_hamiltonian && spinRestricted && wfntype==0 ) 
    closed_shell = true;

  if(closed_shell) {
    app_log()<<"Found closed shell system. " <<std::endl;
  }

  long m1=0, m2=0, m2c=0, m1c=0;

  if(runtype==0) { // only for pure dets


    APP_ABORT("Error: runtype=0 not yet implemented in MultiPureSD. \n\n");

    std::map<IndexType,bool> all_alpha,all_beta; 
    all_alpha.clear();
    all_beta.clear();
    for(IndexType i=0; i<2*NMO; i++) all_alpha[i]=false;
    for(IndexType i=0; i<2*NMO; i++) all_beta[i]=false;  

    for(std::vector<IndexType>::iterator it=occ_orbs.begin(); it!=occ_orbs.end(); it++)
      if(*it < NMO)
        all_alpha[*it] = true;
      else 
        all_beta[*it] = true;
  
//    if(!sHam->createHamiltonianForPureDeterminant(all_alpha,all_beta,hij,SMSpHijkl,cutoff)) {
//      app_error()<<"Error in createHamiltonianForPureDeterminant. \n";
//      return false;
//    }

  } else {

    app_log()<<" MultiPureSingleDeterminant - Creating Hamiltonians. \n";
    SMSpHijkl.resize(ci.size());  
    if(rotated_hamiltonian) haj.resize(ci.size()); 
    else hij.resize(ci.size()); 
    int wtype = closed_shell?0:1;
    ComplexMatrix Am(2*NMO,NAEA);
    int nr = (wfntype==0)?NMO:2*NMO;
    std::map<IndexType,bool> isOcc_alpha; 
    std::map<IndexType,bool> isOcc_beta;
    for(int ki=0; ki<ci.size(); ki++) {
 
      if(rotated_hamiltonian) {
        for(int i=0; i<nr; i++)  
         for(int j=0; j<NAEA; j++)  
          Am(i,j) = OrbMat[ki*orbsize+i*NAEA+j]; 
        SMSpHijkl[ki].setDims(2*NMO*NMO,2*NMO*NMO);
        SMSpHijkl[ki].setup(head_of_nodes,name+std::string("SMSpHijkl_")+std::to_string(ki),TG.getNodeCommLocal());
        if(!sHam->createHamiltonianForGeneralDeterminant(wtype,Am,haj[ki],SMSpHijkl[ki],cutoff)) {
          app_error()<<"Error in createHamiltonianForGeneralDeterminant. \n";
          return false;
        }
        m1+=haj[ki].size();
        m2+=SMSpHijkl[ki].size();
        m1c+=haj[ki].size()*sizeof(s1D<ComplexType>);
        m2c+=SMSpHijkl[ki].size()*(sizeof(ComplexType)+2*sizeof(int));
      } else {
        isOcc_alpha.clear();
        isOcc_beta.clear();
        for(IndexType i=0; i<2*NMO; i++) isOcc_alpha[i]=false;
        for(IndexType i=0; i<NAEA; i++) isOcc_alpha[ occ_orbs[ki*(NAEA+NAEB)+i] ]=true;
        for(IndexType i=0; i<2*NMO; i++) isOcc_beta[i]=false;
        for(IndexType i=NAEA; i<NAEA+NAEB; i++) isOcc_beta[ occ_orbs[ki*(NAEA+NAEB)+i] ]=true;
        SMSpHijkl[ki].setDims(2*NMO*NMO,2*NMO*NMO);
        SMSpHijkl[ki].setup(head_of_nodes,name+std::string("SMSpHijkl_")+std::to_string(ki),TG.getNodeCommLocal());
        bool closed=false;
        if(NAEA==NAEB) {
          bool diff = false;
          for(int i=0; i<NAEA; i++) 
            if(occ_orbs[ki*(NAEA+NAEB)+i] != occ_orbs[ki*(NAEA+NAEB)+NAEA+i]-NMO) {
              diff=true;
              break;
            }
          if(!diff) closed=true; 
        }
        if(!sHam->createHamiltonianForPureDeterminant(isOcc_alpha,isOcc_beta,hij[ki],SMSpHijkl[ki],cutoff, closed)) {
          app_error()<<"Error in createHamiltonianForPureDeterminant. \n";
          return false;
        }
        m1+=hij[ki].size();
        m2+=SMSpHijkl[ki].size();
        m1c+=hij[ki].size()*sizeof(s1D<ValueType>);
        m2c+=SMSpHijkl[ki].size()*(sizeof(ComplexType)+2*sizeof(int));
      }
    
    }

  }

  if(diagHam) {
    app_log()<<" Diagonalizing trial wave function in MultiPureSingleDeterminant. \n";
    std::vector<RealType> eigVal(ci.size());
    ComplexMatrix eigVec(1,ci.size());
    bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,occ_orbs,ci.size());
    if(sucess) {
      app_log()<<" New trial energy and ci coefficients: " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;
      for(int i=0; i<ci.size(); i++) app_log()<<i <<" old: " <<ci[i] <<" new: " <<eigVec(0,i) <<std::endl; 
      for(int i=0; i<ci.size(); i++) ci[i] = eigVec(0,i); 
      app_log()<<std::endl; 
    } else {
      app_error()<<"Error diagonalizing trial wavefunction. \n";
      return false;
    }
  } 

  if(iterCI) {
    app_log()<<" Iterative CI: " <<std::endl;
    iterativeCI(IterCI_cut,1000,1000,IterCI_maxit);
    app_error()<<" Aborting after IterCI \n";
    return false;
  }

  ComplexType epot,ekin,o1,o2;
  HF.resize(2*NMO,NAEA);
  for(int i=0; i<2*NMO; i++)
  for(int j=0; j<NAEA; j++) HF(i,j)=ComplexType(0.0,0.0);
  if(rotated_hamiltonian) {
    std::copy(OrbMat.data(),OrbMat.data()+NAEA*NMO,HF.data());
    if(wfntype==0)
      std::copy(OrbMat.data(),OrbMat.data()+NAEA*NMO,HF.data()+NAEA*NMO);
    else
      std::copy(OrbMat.data()+NAEA*NMO,OrbMat.data()+2*NAEA*NMO,HF.data()+NAEA*NMO);
//    for(int i=0; i<NAEA; i++) HF(i,i)=ComplexType(1.0,0.0);
//    for(int i=0; i<NAEB; i++) HF(i+NMO,i)=ComplexType(1.0,0.0);
  } else {
    for(int i=0; i<NAEA; i++) HF(occ_orbs[i],i)=ComplexType(1.0,0.0);
    for(int i=0; i<NAEB; i++) HF(occ_orbs[i+NAEA],i)=ComplexType(1.0,0.0);
  }

  app_log()<<std::endl <<"*********************************************************************: \n"
           <<" MultiPureSingleDeterminant: \n"
           <<"     Number of terms in ci expansion: " <<ci.size() <<"\n"
           <<"     Number of terms and memory usage of hij:    " <<m1 <<"  " <<m1c/1.0e6 <<"  MB. " <<std::endl
           <<"     Number of terms and memory usage of Vijkl:  " <<m2 <<"  " <<m2c/1.0e6 <<"  MB. " <<std::endl;
  evaluateLocalEnergy(HF.data(),ekin,epot,o1,o2,-1);
  app_log()<<" <PsiT|H|Det[1]>:   " <<epot+ekin <<std::endl
           <<" Ekin:          " <<ekin <<std::endl
           <<" Epot:          " <<epot <<std::endl;
  evaluateTrialEnergy(ekin,epot);
  app_log()<<" <PsiT|H|PsiT>:   " <<epot+ekin <<std::endl
           <<" Ekin:          " <<ekin <<std::endl
           <<" Epot:          " <<epot <<std::endl;
  app_log()<<"*********************************************************************: \n" <<std::endl <<std::endl;

  return true;
}

bool MultiPureSingleDeterminant::hdf_write(hdf_archive& read, const std::string& tag, bool include_tensors) 
{return true;}

bool MultiPureSingleDeterminant::hdf_write()
{return true;}

void MultiPureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrixFull(const ComplexType* SlaterMat, ComplexType& ovl, ComplexMatrix& dm, bool full)
{

  ComplexType o1;
  ovl = ComplexType(0.0);
  dm = ComplexType(0.0);

  local_evaluateOneBodyMixedDensityMatrix(ref,SlaterMat,overlaps[2*ref],overlaps[2*ref+1],mixed_density_matrix,full);
  o1 = myconj(ci[ref])*overlaps[2*ref]*overlaps[2*ref+1];
  dm = o1*mixed_density_matrix;
  ovl = o1; 

  for(int i=0; i<ci.size(); i++)
  { 
    if(i==ref) 
      continue;

    local_rankUpdateOneBodyMixedDensityMatrix(i,SlaterMat,overlaps[2*i],overlaps[2*i+1],full);
    o1 = myconj(ci[i])*overlaps[2*i]*overlaps[2*i+1];
    dm += o1*rank_updated_mixed_density_matrix;
    ovl += o1; 
  }
  dm *= (ComplexType(1.0,0.0)/ovl); 
}

// might need to save inverse(S0) for fast update formula, don't know yet
void MultiPureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrix(int detn, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, ComplexMatrix& dm, bool full)
{

  // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
  // lots of simplification because A is an identity matrix with possibly exchanged columns, 
  // look at derivations for information

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  IndexType* it_orbs = occ_orbs.data()+(NAEA+NAEB)*detn;
  ComplexType *Orbs = OrbMat.data()+detn*orbsize;
  ComplexType *Orbs_beta = Orbs+((wfntype==0)?0:NAEA*NMO);   

  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,Orbs,NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);
    
  } else {  

    // copy rows corresponding to occupied orbitals to S0
    //   S0(i,:) = SlaterMat(ik,:), where ik is the ith occupied orbital    
    ComplexMatrix::iterator itS0 = S0.begin(); 
    for(IndexType* it=it_orbs; it!=it_orbs+NAEA; it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEA;
      for(; itSM!=itSMend; ++itSM, ++itS0)
        *itS0 = *itSM;
    }
  }

  ovl_alpha = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

  // SS0 = SlaterMat * S0
  DenseMatrixOperators::product(NMO,NAEA,NAEA,SlaterMat,NAEA,S0.data(),NAEA,SS0.data(),NAEA);

  for(ComplexMatrix::iterator it=dm.begin(); it!=dm.end(); it++) *it=zero;

  // copy to dm 
  if(rotated_hamiltonian) {

    if(full) {
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,Orbs,NAEA,zero,dm.data(),NMO);
      DenseMatrixOperators::transpose(NMO,dm.data(),NMO);
    } else {
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        dm(i,j) = SS0(j,i);
    }     
  } else { 
    for(int i=0; i<NAEA; i++)
     for(int j=0; j<NMO; j++)
      dm(*(it_orbs+i),j) = SS0(j,i);
  }

  if(closed_shell) {
    std::copy(dm.data(),dm.data()+NMO*NMO,dm.data()+NMO*NMO);
    ovl_beta=ovl_alpha;
    return;
  }

  // repeat for beta
  if(rotated_hamiltonian) {

    // S0 = transpose(conjg(A))*B     
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,Orbs_beta,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);

  } else {
 
    ComplexMatrix::iterator itS1 = S1.begin(); 
    for(IndexType* it=it_orbs+NAEA; it!=it_orbs+NAEA+NAEB; it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEB;
      for(; itSM!=itSMend; ++itSM, ++itS1)
        *itS1 = *itSM;
    }
  }

  ovl_beta = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());
  
  // SS0(beta) = SlaterMat(beta) * S1
  DenseMatrixOperators::product(NMO,NAEB,NAEB,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,SS0.data()+NAEA*NMO,NAEA);

  // copy to dm
  if(rotated_hamiltonian) {
    if(full) {
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,Orbs_beta,NAEA,zero,dm.data()+NMO*NMO,NMO);
      DenseMatrixOperators::transpose(NMO,dm.data()+NMO*NMO,NMO);
    } else {
      for(int i=0; i<NAEB; i++)
       for(int j=0; j<NMO; j++)
        dm(i+NMO,j) = SS0(j+NMO,i);
    }
  } else {
    for(int i=0; i<NAEB; i++)
     for(int j=0; j<NMO; j++)
      dm(*(it_orbs+NAEA+i),j) = SS0(j+NMO,i);
  }  

} 

// might need to save inverse(S0) for fast update formula, don't know yet
ComplexType MultiPureSingleDeterminant::local_evaluateOverlapSlaterDet(int detn, const ComplexType* SlaterMat)
{

  // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
  // lots of simplification because A is an identity matrix with possibly exchanged columns, 
  // look at derivations for information

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  ComplexType ovl = zero;
  if(rotated_hamiltonian) {

    ComplexType *Orbs = OrbMat.data()+detn*orbsize;
    ComplexType *Orbs_beta = Orbs+((wfntype==0)?0:NAEA*NMO);
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one, Orbs ,NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);
    ovl = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());
    if(closed_shell) return ovl*ovl;
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one, Orbs_beta ,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);
    ovl *= Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());   

  } else {

    IndexType* it_orbs = occ_orbs.data()+(NAEA+NAEB)*detn;

    // copy rows corresponding to occupied orbitals to S0
    //   S0(i,:) = SlaterMat(ik,:), where ik is the ith occupied orbital    
    ComplexMatrix::iterator itS0 = S0.begin(); 
    for(IndexType* it=it_orbs; it!=it_orbs+NAEA; it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEA;
      for(; itSM!=itSMend; ++itSM, ++itS0)
        *itS0 = *itSM;
    }
    ovl = Invert(S0.data(), NAEA, NAEA,Cwork.data(),pivot.data());
  
    // repeat for beta
    ComplexMatrix::iterator itS1 = S1.begin(); 
    for(IndexType* it=it_orbs+NAEA; it!=it_orbs+NAEA+NAEB; it++) {
      const ComplexType* itSM = SlaterMat+ (*it) *NAEA;
      const ComplexType* itSMend = itSM+NAEB;
      for(; itSM!=itSMend; ++itSM, ++itS1)
        *itS1 = *itSM;
    }
    ovl *= Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());
  }
 
  return ovl; 
} 

void MultiPureSingleDeterminant::local_rankUpdateOneBodyMixedDensityMatrix(int n, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool full)
{

  local_evaluateOneBodyMixedDensityMatrix(n,SlaterMat,ovl_alpha,ovl_beta,rank_updated_mixed_density_matrix,full); 

  // implement fast update routine
} 

void MultiPureSingleDeterminant::local_evaluateOneBodyTrialDensityMatrix(bool full)
{
  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  for(int i=0; i<2*NMO; i++)
   for(int j=0; j<NMO; j++)
    trial_density_matrix(i,j) = zero; 
  
  IndexType* it = occ_orbs.data()+(NAEA+NAEB)*ref;
  for(int i=0; i<NAEA; i++,it++)
    trial_density_matrix(*it,*it) = one;

  for(int i=0; i<NAEB; i++,it++)
    trial_density_matrix(*it,*it-NMO) = one;
} 

// no need for low rank updates, this is trivial
void MultiPureSingleDeterminant::local_rankUpdateOneBodyTrialDensityMatrix(int n, bool full)
{
  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  for(int i=0; i<2*NMO; i++)
   for(int j=0; j<NMO; j++)
    rank_updated_trial_density_matrix(i,j) = zero; 

  if(n > ci.size()) {
    app_log()<<"ERROR ERROR ERROR: Not aborting in call to local_rankUpdateevaluateOneBodyTrialDensityMatrix(n) with n>ndets. \n" <<std::endl;
    return;
  }
  IndexType* it = occ_orbs.data()+(NAEA+NAEB)*n;
  for(int i=0; i<NAEA; i++,it++)
    rank_updated_trial_density_matrix(*it,*it) = one;

  for(int i=0; i<NAEB; i++,it++) {
    rank_updated_trial_density_matrix(*it,*it-NMO) = one;
  }
} 

  void MultiPureSingleDeterminant::evaluateOverlap(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
{
  ovl_beta = ComplexType(1.0,0.0);
  ovl_alpha = ComplexType(0.0,0.0);

  for(int i=0; i<ci.size(); i++)
  {
    ovl_alpha += myconj(ci[i]) * local_evaluateOverlapSlaterDet(i,SlaterMat); 
  }

}

  void MultiPureSingleDeterminant::evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

    register ComplexType one = ComplexType(1.0,0.0);
    register ComplexType zero = ComplexType(0.0,0.0);
    register ComplexType lpe=ComplexType(0.0,0.0),lke=ComplexType(0.0,0.0);
    register ComplexType deno=ComplexType(0.0,0.0);
    register ComplexType ke_nume=ComplexType(0.0,0.0), pe_nume=ComplexType(0.0,0.0);
    IndexType ef = 0, NAE = NAEA+NAEB, NAE2 = (NAEA+NAEB)*(NAEA+NAEB);
    int nr1=1, nc1=2*NMO*NMO;

// FIX FIX FIX: This routine is inefficient for orthogonal CI expansions, since it calculates everything from scratch ignoring the fact that different configurations share either alpha or beta determinants. Modify to calculate ONLY the list of different alpha/beta determinants and    

#ifdef AFQMC_TIMER
    Timer.start("MultiPureSingleDeterminant:evaluateLocalEnergy");
#endif

    // clean arrays 
    for(std::vector<ComplexType>::iterator it=pe.begin(); it!=pe.end(); ++it) *it=0.0;
    for(std::vector<ComplexType>::iterator it=ke.begin(); it!=ke.end(); ++it) *it=0.0;


// change this later to do P determinants simultaneously using SpMatMat operations, for higher flop rate possibly. 
    for(int ki=0; ki<ci.size(); ki++)
    {
      lke=lpe=ComplexType(0.0,0.0);

      local_evaluateOneBodyMixedDensityMatrix(ki,SlaterMat,overlaps[2*ki],overlaps[2*ki+1],mixed_density_matrix,false);
      ComplexType scl = myconj(ci[ki])*overlaps[2*ki]*overlaps[2*ki+1];
      deno += scl; 

      if(runtype == 0) {
   
        APP_ABORT("Error: runtype = 0 not yet implemented in MultiPureSlaterDeterminant. \n\n\n");

      } else {

        ComplexMatrix::iterator itG = mixed_density_matrix.begin();
        if(rotated_hamiltonian) {
          std::vector<s1D<ComplexType> >::iterator  end1 = haj[ki].end();
          for(std::vector<s1D<ComplexType> >::iterator it = haj[ki].begin(); it != end1; it++)
            lke += *(itG + std::get<0>(*it)) * std::get<1>(*it);
        } else {
          s1Dit end1 = hij[ki].end();
          for(s1Dit it = hij[ki].begin(); it != end1; it++)
            lke += *(itG + std::get<0>(*it)) * std::get<1>(*it);
        }
        ke[ki]=lke;
        ke_nume += lke*scl;

        SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(nc1,nc1,one,SMSpHijkl[ki],mixed_density_matrix.data(),zero,V0.data());
        itG = mixed_density_matrix.begin();
        ComplexVector::iterator itV = V0.begin();
        for(int i=0; i<nc1; i++,++itG,++itV) lpe += (*itV) * (*itG);
        pe[ki] = 0.5*lpe+NuclearCoulombEnergy;
        pe_nume += (0.5*lpe)*scl;

      } 

    }

    ovl_alpha = deno;
    ovl_beta = ComplexType(1.0,0.0);

    if( std::abs(deno) > 1e-8 ) { 
      epot = pe_nume/deno+NuclearCoulombEnergy;
      ekin = ke_nume/deno;
    } else {
      epot = 0.0;
      ekin = 0.0;
    }
#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:evaluateLocalEnergy");
#endif
   
  }

/*
  void MultiPureSingleDeterminant::evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

    register ComplexType lpe=ComplexType(0.0,0.0),lke=ComplexType(0.0,0.0);
    register ComplexType deno=ComplexType(0.0,0.0);
    register ComplexType ke_nume=ComplexType(0.0,0.0), pe_nume=ComplexType(0.0,0.0);
    IndexType ef = 0, NAE = NAEA+NAEB, NAE2 = (NAEA+NAEB)*(NAEA+NAEB);

#ifdef AFQMC_TIMER
    Timer.start("MultiPureSingleDeterminant:evaluateLocalEnergy");
#endif

    // clean arrays 
    for(std::vector<ComplexType>::iterator it=pe.begin(); it!=pe.end(); ++it) *it=0.0;
    for(std::vector<ComplexType>::iterator it=ke.begin(); it!=ke.end(); ++it) *it=0.0;

    // will access this vector sequentially 
    std::vector<IndexType>::iterator it_orbs;
    std::vector<IndexType>::iterator it_pairs;

    // **********************************
    //         !!!REFERENCE!!!
    // **********************************
    local_evaluateOneBodyMixedDensityMatrix(ref,SlaterMat,overlaps[2*ref],overlaps[2*ref+1],mixed_density_matrix);
    deno += myconj(ci[ref])*overlaps[2*ref]*overlaps[2*ref+1];

    // move it_orbs to the beginning of ref's data
    it_orbs = occ_orbs.begin()+NAE*ref; 
    ComplexMatrix::iterator itG = mixed_density_matrix.begin();
    for(int i=0; i<NAE; i++, it_orbs++)
     for(s1Dit it1 = hij_indx[*it_orbs]; it1 != hij_indx[*it_orbs+1]; it1++) 
      lke += *(itG + std::get<0>(*it1)) * std::get<1>(*it1);
    ke[ref] = lke;
    ke_nume += lke*myconj(ci[ref])*overlaps[2*ref]*overlaps[2*ref+1];

    // move it_pairs to the beginning of ref's data
    int cnt=0;
    for(int i=0; i<ref; i++) cnt+=Vijkl_nterms_per_det[i];
    std::vector<s2Dit>::iterator it = Vijkl_indx.begin()+cnt;
    for(int i=0; i<Vijkl_nterms_per_det[ref]; i++,it++)
     lpe += (*(itG + std::get<0>(**it))) * (*(itG + std::get<1>(**it))) * std::get<2>(**it);
    pe[ref] = 0.5*lpe+NuclearCoulombEnergy;
    pe_nume += (0.5*lpe+NuclearCoulombEnergy)*myconj(ci[ref])*overlaps[2*ref]*overlaps[2*ref+1];

    // process other determinants using low-rank updates  
    
    it_orbs = occ_orbs.begin();
    cnt=0;
    it=Vijkl_indx.begin();
    for(int i=0; i<ci.size(); i++)
    {
      if(i==ref) {
        it+=Vijkl_nterms_per_det[i];
        it_orbs += NAE;
        continue;
      }
      lke=lpe=ComplexType(0.0,0.0);

      local_rankUpdateOneBodyMixedDensityMatrix(i,SlaterMat,overlaps[2*i],overlaps[2*i+1]);
      ComplexType scl = myconj(ci[i])*overlaps[2*i]*overlaps[2*i+1];
      deno += scl; 

      itG = rank_updated_mixed_density_matrix.begin();
      for(int k=0; k<NAE; k++, it_orbs++)
       for(s1Dit it1 = hij_indx[*it_orbs]; it1 != hij_indx[*it_orbs+1]; it1++)
        lke += *(itG + std::get<0>(*it1)) * std::get<1>(*it1);
      ke[i] = lke;
      ke_nume += lke*scl;
      
      for(int k=0; k<Vijkl_nterms_per_det[i]; k++,it++) 
       lpe += (*(itG + std::get<0>(**it))) * (*(itG + std::get<1>(**it))) * std::get<2>(**it);
      pe[i] = 0.5*lpe+NuclearCoulombEnergy;
      pe_nume += (0.5*lpe+NuclearCoulombEnergy)*scl;

    } 
    ovl_alpha = deno;
    ovl_beta = ComplexType(1.0,0.0);
  
    epot = pe_nume/deno;
    ekin = ke_nume/deno;

#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:evaluateLocalEnergy");
#endif
   
  }
*/

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("MultiPureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif
    ComplexType ovlp;
    local_evaluateOneBodyMixedDensityMatrixFull(SlaterMat,ovlp,full_mixed_density_matrix,true);

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(transposed) {
      SparseMatrixOperators::product_SpMatV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta)
        SparseMatrixOperators::product_SpMatV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    } else {
      SparseMatrixOperators::product_SpMatTV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta)
        SparseMatrixOperators::product_SpMatTV<ComplexSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("MultiPureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
    Timer.start("MultiPureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrixFull");
#endif
    ComplexType ovlp;
    local_evaluateOneBodyMixedDensityMatrixFull(SlaterMat,ovlp,full_mixed_density_matrix,true);

#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrixFull");
    Timer.start("MultiPureSingleDeterminant:product_Spvn_x_DMfull");
#endif

    ComplexType one = ComplexType(1.0,0.0);
    ComplexType zero = ComplexType(0.0,0.0);
    if(transposed) {
      SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta)
        SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    } else {
      SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta)
        SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:product_Spvn_x_DMfull");
    Timer.stop("MultiPureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {
    APP_ABORT(" Error: Routine not implemented: MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators. \n\n\n");
  } 
    
  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n)
  {
    APP_ABORT(" Error: Routine not implemented: MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators. \n\n\n");
  } 

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSpMat& Spvn, std::vector<ComplexType>& v, const int n)
  {
  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, const int n)
  {
  }

void MultiPureSingleDeterminant::evaluateTrialEnergy(ComplexType& ke, ComplexType& pe)
{
  ComplexType one = ComplexType(1.0,0.0);
  ComplexType zero = ComplexType(0.0,0.0);
  
  ke=pe=zero;
  if(myComm->rank() == 0) {  // serial for now

   if(rotated_hamiltonian) {  // since I do not have a way to do this right now  

    int nr1=1, nc1=2*NMO*NMO;
    ComplexType oa,ob,pe_nume=zero,deno=zero;
    ComplexType lke,lpe,ke_nume=zero; 
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfntype==0)?0:NAEA*NMO; 

    for(int i=0; i<ci.size(); i++) {

      lke=lpe=ComplexType(0.0,0.0);
      
      for(int ii=0; ii<NMO; ii++)
       for(int j=0; j<NAEA; j++)
        Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
      for(int ii=0; ii<NMO; ii++)
       for(int j=0; j<NAEB; j++)
        Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];
      local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,rank_updated_trial_density_matrix);
      
      ComplexMatrix::iterator itG = rank_updated_trial_density_matrix.begin();
      std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
        lke += *(itG + std::get<0>(*it)) * std::get<1>(*it);
      
      SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(nc1,nc1,one,SMSpHijkl[i],rank_updated_trial_density_matrix.data(),zero,V0.data());
      ComplexVector::iterator itV = V0.begin(); 
      for(int ii=0; ii<nc1; ii++,++itG,++itV) lpe += (*itV) * (*itG);
      
      ComplexType scl = myconj(ci[i])*ci[i]*oa*ob;
      ke_nume += lke*scl;
      pe_nume += 0.5*lpe*scl;
      deno += scl;

      for(int j=i+1; j<ci.size(); j++) {

        lke=lpe=ComplexType(0.0,0.0);

        for(int ii=0; ii<NMO; ii++)
         for(int jj=0; jj<NAEA; jj++)
          Am(ii,jj) = OrbMat[j*orbsize+ii*NAEA+jj];
        for(int ii=0; ii<NMO; ii++)
         for(int jj=0; jj<NAEB; jj++)
          Am(ii+NMO,jj) = OrbMat[j*orbsize+ii*NAEA+jj+del];

        // < Di | c!c | Dj> / <Di | Dj> ,  <Di | Dj> 
        local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,rank_updated_trial_density_matrix,false);

        ComplexMatrix::iterator itG = rank_updated_trial_density_matrix.begin();
        std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
        for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
          lke += *(itG + std::get<0>(*it)) * std::get<1>(*it);

        SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(nc1,nc1,one,SMSpHijkl[i],rank_updated_trial_density_matrix.data(),zero,V0.data());
        ComplexVector::iterator itV = V0.begin();
        for(int ii=0; ii<nc1; ii++,++itG,++itV) lpe += (*itV) * (*itG);

        ComplexType scl = myconj(ci[i])*ci[j]* oa*ob;
        ke_nume += lke*scl + myconj(lke*scl);
        pe_nume += 0.5*lpe*scl + myconj(0.5*lpe*scl);
        deno += scl + myconj(scl);

      }
    }
    pe = pe_nume/deno+NuclearCoulombEnergy;
    ke = ke_nume/deno;
   } else {

      for(int i=0; i<ci.size(); i++) 
        std::sort(occ_orbs.begin()+i*(NAEA+NAEB),occ_orbs.begin()+(i+1)*(NAEA+NAEB));

      ComplexType zero=ComplexType(0);
      pe=ke=zero;
      RealType sg;
      std::vector<IndexType> occ(NAEA+NAEB); 
      IndexType n0,n1,n2,n3;
      std::vector<IndexType> DL(NAEA+NAEB);
      std::vector<IndexType> DR(NAEA+NAEB);
      ComplexType deno=zero;

      for(int ki=0; ki<ci.size(); ki++) {
        std::vector<IndexType>::iterator it = occ_orbs.begin()+ki*(NAEA+NAEB);
        for(int i=0; i<NAEA+NAEB; i++)
        {
          ke += std::conj(ci[ki])*ci[ki]*sHam->H(*(it+i),*(it+i));
          for(int j=i+1; j<NAEA+NAEB; j++) 
            pe += std::conj(ci[ki])*ci[ki]*(sHam->H(*(it+i),*(it+j),*(it+i),*(it+j)) - sHam->H(*(it+i),*(it+j),*(it+j),*(it+i)));
        }
        deno += std::conj(ci[ki])*ci[ki];

        for(int kj=ki+1; kj<ci.size(); kj++) {

          std::vector<IndexType>::iterator jt = occ_orbs.begin()+kj*(NAEA+NAEB);
          std::copy(it,it+NAEA+NAEB,DL.begin());
          std::copy(jt,jt+NAEA+NAEB,DR.begin());
          int cnt = cntExcitations(DL,DR,n0,n1,n2,n3,occ,sg);

          if(cnt==0) {
            app_error()<<" Error: Found repeated determinant in trial wave function in MultiPureSingleDeterminant \n";
            APP_ABORT(" Error: \n\n\n");
          } else if(cnt==2) {
            int nterms = NAEA+NAEB-1;
            ValueType hij = sHam->H(n0,n1);
            ke+=std::conj(ci[ki])*ci[kj]*sg*hij;
            ke+=std::conj(ci[kj])*ci[ki]*sg*myconj(hij);
            ComplexType et=zero; 
            for(int i=0; i<nterms; i++) 
              et += sHam->H(n0,occ[i],n1,occ[i]) - sHam->H(n0,occ[i],occ[i],n1); 
            pe+=std::conj(ci[ki])*ci[kj]*sg*et;
            pe+=std::conj(ci[kj])*ci[ki]*sg*myconj(et);
          } else if(cnt==4) {
            ValueType hij = sHam->H(n0,n1,n2,n3) - sHam->H(n0,n1,n3,n2); 
            pe += std::conj(ci[ki])*ci[kj]*sg*hij;
            pe += std::conj(ci[kj])*ci[ki]*sg*myconj(hij);
          }
        }
      }
      pe = pe/deno+NuclearCoulombEnergy;
      ke /= deno;
   }
  }

  // lazy for now
  myComm->bcast(pe);
  myComm->bcast(ke);

}



void MultiPureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSMSpMat& Spvn, std::vector<ComplexType>& v, const int nstep)
{

  ComplexType one = ComplexType(1.0,0.0);
  ComplexType zero = ComplexType(0.0,0.0);
  int n[4];
  RealType sg;
  ComplexType deno=zero;
 
  trial_density_matrix = zero;  

  if(rotated_hamiltonian) {

    ComplexType o1,oa,ob;
    ComplexType ovl = ComplexType(0.0);
    trial_density_matrix = ComplexType(0.0);
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfntype==0)?0:NAEA*NMO;

    for(int i=0; i<ci.size(); i++) {
     // i==j
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEA; j++)
       Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEB; j++)
       Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];

     local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, rank_updated_trial_density_matrix,true);
     o1 = myconj(ci[i])*ci[i]*oa*ob;
     trial_density_matrix += o1*rank_updated_trial_density_matrix; 
     ovl += o1 ;
     for(int j=i+1; j<ci.size(); j++)
     {
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEA; jj++)
         Am(ii,jj) = OrbMat[j*orbsize+ii*NAEA+jj];
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEB; jj++)
         Am(ii+NMO,jj) = OrbMat[j*orbsize+ii*NAEA+jj+del];
       local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, rank_updated_trial_density_matrix,true);
       o1 = myconj(ci[i])*ci[j]*oa*ob; 
       trial_density_matrix += o1*rank_updated_trial_density_matrix; 
       ovl += o1 ;
       o1 = ci[i]*myconj(ci[j])*myconj(oa*ob);
       // is there an operation for std::complex conjugation???
       for(int ii=0; ii<NMO; ii++)  
        for(int jj=0; jj<NMO; jj++) { 
         trial_density_matrix(ii,jj) += o1*myconj(rank_updated_trial_density_matrix(jj,ii));
         trial_density_matrix(ii+NMO,jj) += o1*myconj(rank_updated_trial_density_matrix(jj+NMO,ii));
        }
       ovl += o1 ;
     }
    }
    trial_density_matrix *= (ComplexType(1.0,0.0)/ovl);

  } else {
    for(int i=0; i<ci.size(); i++) {
     // add diagonal term 
     local_rankUpdateOneBodyTrialDensityMatrix(i); 
     deno += std::norm(ci[i]);
     trial_density_matrix += std::norm(ci[i])*rank_updated_trial_density_matrix;
     for(int j=i+1; j<ci.size(); j++)
     {
       int nex = cmpDets(NAEA,NAEB,n,sg, occ_orbs.begin()+(NAEA+NAEB)*i , occ_orbs.begin()+(NAEA+NAEB)*j, Iwork);
       switch(nex) {
         case 0:  //
         {
           deno += std::conj(ci[i])*ci[j] + std::conj(ci[j])*ci[i];
           trial_density_matrix += (std::conj(ci[i])*ci[j] + std::conj(ci[j])*ci[i])*rank_updated_trial_density_matrix;
           // should I allow this? Shouldn't happen for a well constructed det list 
           break;
         }
         case 2:  // single excitation
         {
           // only non-zero term is G(i,a)
           trial_density_matrix(n[0],n[1]) += std::conj(ci[i])*ci[j];
           trial_density_matrix(n[1],n[0]) += std::conj(ci[j])*ci[i];
           break;
         }
         case 4:   // double excitation
         {
           // do nothing, doesn't contribute to 1-Body operator
           break;
         }
 
       }
     }
    }
    trial_density_matrix *= (one/deno);
  }


  SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,trial_density_matrix.data(),zero,v.data());
  if(addBetaBeta)
    SparseMatrixOperators::product_SpMatTV<ComplexSMSpMat>(Spvn.rows(),Spvn.cols(),one,Spvn,trial_density_matrix.data()+NMO*NMO,one,v.data());

}

void MultiPureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSpMat& Spvn, std::vector<ComplexType>& v, const int nstep)
{

    APP_ABORT("Error in MultiPureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators \n\n\n");

}

bool MultiPureSingleDeterminant::diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ComplexMatrix& eigVec, std::vector<IndexType>& occv, int nci, bool eigV )
{
  ComplexType one = ComplexType(1.0,0.0);
  ComplexType zero = ComplexType(0.0,0.0);
  bool sucess;
  //int nci = ci.size(); 
  
  if(rotated_hamiltonian) { 

   // serial for now
   if(myComm->rank() == 0) { 
  
    int nr1=1, nc1=2*NMO*NMO;
    ComplexType oa,ob,pe_nume=zero,deno=zero;
    ComplexType lke,lpe,ke_nume=zero; 
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfntype==0)?0:NAEA*NMO; 

    ComplexMatrix ov(nci,nci);
    ComplexMatrix hm(nci,nci);
    ComplexType let = zero;

    for(int i=0; i<nci; i++) {
  
      let = zero;
      
      for(int ii=0; ii<NMO; ii++)
       for(int j=0; j<NAEA; j++)
        Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
      for(int ii=0; ii<NMO; ii++)
       for(int j=0; j<NAEB; j++)
        Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];
      local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,rank_updated_trial_density_matrix);
      
      ComplexMatrix::iterator itG = rank_updated_trial_density_matrix.begin();
      SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(nc1,nc1,one,SMSpHijkl[i],rank_updated_trial_density_matrix.data(),zero,V0.data());
      ComplexVector::iterator itV = V0.begin(); 
      for(int ii=0; ii<nc1; ii++,++itG,++itV) let += (*itV) * (*itG);
      let *= 0.5;

      itG = rank_updated_trial_density_matrix.begin();
      std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
        let += *(itG + std::get<0>(*it)) * std::get<1>(*it);
      
      hm(i,i) = oa*ob*let; 
      ov(i,i) = oa*ob;

      for(int j=i+1; j<nci; j++) {

        let=zero;

        for(int ii=0; ii<NMO; ii++)
         for(int jj=0; jj<NAEA; jj++)
          Am(ii,jj) = OrbMat[j*orbsize+ii*NAEA+jj];
        for(int ii=0; ii<NMO; ii++)
         for(int jj=0; jj<NAEB; jj++)
         Am(ii+NMO,jj) = OrbMat[j*orbsize+ii*NAEA+jj+del];

        // < Di | c!c | Dj> / <Di | Dj> ,  <Di | Dj> 
        local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,rank_updated_trial_density_matrix,false);

        ComplexMatrix::iterator itG = rank_updated_trial_density_matrix.begin();
        SparseMatrixOperators::product_SpMatV<ComplexSMSpMat>(nc1,nc1,one,SMSpHijkl[i],rank_updated_trial_density_matrix.data(),zero,V0.data());
        ComplexVector::iterator itV = V0.begin();
        for(int ii=0; ii<nc1; ii++,++itG,++itV) let += (*itV) * (*itG);
        let *= 0.5;

        itG = rank_updated_trial_density_matrix.begin();
        std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
        for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
          let += *(itG + std::get<0>(*it)) * std::get<1>(*it);

        hm(i,j) = oa*ob*let; 
        ov(i,j) = oa*ob;
        hm(j,i) = myconj(hm(i,j)); 
        ov(j,i) = myconj(ov(i,j));

      }
    }

    eigVal.resize(1);
    eigVec.resize(1,nci);
    std::vector<int> ifail(nci); 
    sucess = DenseMatrixOperators::genHermitianEigenSysSelect(nci,hm.data(),nci,ov.data(),nci,1,eigVal.data(),eigV,eigVec.data(),eigVec.size2(),ifail.data()); 
   } else {
    eigVal.resize(1);
    eigVec.resize(1,nci);
   }
    
   myComm->bcast(sucess);
   myComm->bcast(eigVal.data(),eigVal.size(),0,myComm->getMPI());
   myComm->bcast(eigVec.data(),eigVec.size1()*eigVec.size2(),0,myComm->getMPI());

   return sucess;

  } else {

    // since I'm regenerating ci's ,sort occup std::strings
    
    for(int i=0; i<nci; i++) 
      std::sort(occv.begin()+i*(NAEA+NAEB),occv.begin()+(i+1)*(NAEA+NAEB));

    if(myComm->rank()==0) {

      Timer.reset("Generic3");
      Timer.start("Generic3");
      ComplexMatrix hm(nci,nci);
      ComplexMatrix ov(nci,nci);
      ComplexType let;
      RealType sg;
      std::vector<IndexType> occ(NAEA+NAEB); 
      IndexType n0,n1,n2,n3;
      std::vector<IndexType> DL(NAEA+NAEB);
      std::vector<IndexType> DR(NAEA+NAEB);

      // don't rely on H2_2bar in case it is removed
      for(int ki=0; ki<nci; ki++) {
        // i==j
        let=zero;
        std::vector<IndexType>::iterator it = occv.begin()+ki*(NAEA+NAEB); 
        for(int i=0; i<NAEA+NAEB; i++)
        { 
          let += sHam->H(*(it+i),*(it+i));
          for(int j=i+1; j<NAEA+NAEB; j++) { 
            let += sHam->H(*(it+i),*(it+j),*(it+i),*(it+j)) - sHam->H(*(it+i),*(it+j),*(it+j),*(it+i));
          }
        }
 
        ov(ki,ki) = one;
        hm(ki,ki) = let; 
        for(int kj=ki+1; kj<nci; kj++) {

          std::vector<IndexType>::iterator jt = occv.begin()+kj*(NAEA+NAEB); 
          std::copy(it,it+NAEA+NAEB,DL.begin());
          std::copy(jt,jt+NAEA+NAEB,DR.begin());
          int cnt = cntExcitations(DL,DR,n0,n1,n2,n3,occ,sg);     
          
          if(cnt==0) {
            app_error()<<" Error: Found repeated determinant in trial wave function in MultiPureSingleDeterminant \n"; 
            return false;
          } else if(cnt==2) {
            int nterms = NAEA+NAEB-1;
            let=sHam->H(n0,n1);
            for(int i=0; i<nterms; i++)
              let+=sHam->H(n0,occ[i],n1,occ[i]) - sHam->H(n0,occ[i],occ[i],n1);
            hm(ki,kj)=let*sg; 
          } else if(cnt==4) {
            hm(ki,kj) = sg*(sHam->H(n0,n1,n2,n3) - sHam->H(n0,n1,n3,n2)); 
          } else {
            hm(ki,kj) =  zero;
          }

          ov(ki,kj) = ov(kj,ki) = zero;
          hm(kj,ki) = myconj(hm(ki,kj)); 
        }
      }
      Timer.stop("Generic3");
      //app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic2") <<std::endl;

      Timer.reset("Generic4");
      Timer.start("Generic4");

      eigVal.resize(1);
      eigVec.resize(1,nci);
      std::vector<int> ifail(nci);
      sucess = DenseMatrixOperators::genHermitianEigenSysSelect(nci,hm.data(),nci,ov.data(),nci,1,eigVal.data(),eigV,eigVec.data(),eigVec.size2(),ifail.data());

      Timer.stop("Generic4");
      //app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic2") <<std::endl;

    } else {
      eigVal.resize(1);
      eigVec.resize(1,nci);
    }

    myComm->bcast(sucess);
    myComm->bcast(eigVal.data(),eigVal.size(),0,myComm->getMPI());
    myComm->bcast(eigVec.data(),eigVec.size1()*eigVec.size2(),0,myComm->getMPI());

    return sucess;
  }

}

//
void MultiPureSingleDeterminant::iterativeCI(double cutoff, int nmax, int nmax_int, int maxit)
{

// algorithm (simple and dumb for now)

    int ne = NAEA+NAEB; 
    std::vector<IndexType> intm;
    intm.reserve(ne*100000); 
    std::vector<RealType> eigVal(1);
    std::vector<IndexType> vira,virb;
    std::vector<IndexType> indx(ne); 
    vira.reserve(NMO-NAEA);
    virb.reserve(NMO-NAEA);
    ComplexMatrix eigVec;
    
    int nn = occ_orbs.size()/ne;  
    for(int i=0; i<nn; i++)
      std::sort(occ_orbs.begin()+i*ne,occ_orbs.begin()+(i+1)*ne);

    for(int i=0; i<maxit; i++) {
      intm.clear();
      int nci = occ_orbs.size()/ne;  
      int nterms=0;
      Timer.reset("Generic");
      Timer.start("Generic");
      std::vector<IndexType>::iterator it=occ_orbs.begin();   
      for(int nc = 0; nc<nci; nc++, it+=ne) { 
        vira.clear();
        virb.clear();
        for(int iv=0; iv<NMO; iv++)  
          if(!binary_search (it, it+NAEA, iv)) vira.push_back(iv); 
        for(int iv=NMO; iv<2*NMO; iv++)  
          if(!binary_search (it+NAEA, it+ne, iv)) virb.push_back(iv); 
        if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
        if(std::search(intm.begin(),intm.end(),it,it+ne)==intm.end()) {
          for(int ii=0; ii<ne; ii++) intm.push_back(*(it+ii));
          nterms++;
        } 
        // double excitations 
        for(int ia=0; ia<NAEA; ia++) {
          int a = *(it+ia); 
          // aa
          for(int ib=ia+1; ib<NAEA; ib++) {
            int b = *(it+ib); 
            for(int ic=0; ic<NMO-NAEA; ic++) {
              int c = vira[ic]; 
              for(int id=ic+1; id<NMO-NAEA; id++) {
                int d = vira[id]; 
                ValueType Hij = sHam->H(a,b,c,d) - sHam->H(a,b,d,c); 
                if(std::abs(Hij*ci[nc]) > cutoff) {
                  if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                  int sz = intm.size();
                  indx.clear();
                  for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                  indx[ia] = c;
                  indx[ib] = d;
                  std::sort(indx.begin(),indx.end()); 
                  if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                    for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                    nterms++;
                  } 
                  if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                  sz = intm.size();
                  indx.clear();
                  for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                  indx[ia] = d;
                  indx[ib] = c;
                  std::sort(indx.begin(),indx.end());                            
                  if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                    for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                    nterms++;
                  }
                }
              } 
            } 
          }
          // ab
          for(int ib=NAEA; ib<ne; ib++) {
            int b = *(it+ib); 
            for(int ic=0; ic<NMO-NAEA; ic++) {
              int c = vira[ic]; 
              for(int id=0; id<NMO-NAEB; id++) {
                int d = virb[id]; 
                ValueType Hij = sHam->H(a,b,c,d);
                if(std::abs(Hij*ci[nc]) > cutoff) {
                  if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                  int sz = intm.size();
                  indx.clear();
                  for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                  indx[ia] = c; 
                  indx[ib] = d; 
                  std::sort(indx.begin(),indx.end());
                  if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                    for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                    nterms++;
                  } 
                } 
              } 
            } 
          } 
        }  
        for(int ia=NAEA; ia<ne; ia++) {
          int a = *(it+ia);
          for(int ib=ia+1; ib<ne; ib++) {
            int b = *(it+ib);
            for(int ic=0; ic<NMO-NAEB; ic++) {
              int c = virb[ic];
              for(int id=ic+1; id<NMO-NAEB; id++) {
                int d = virb[id];
                ValueType Hij = sHam->H(a,b,c,d) - sHam->H(a,b,d,c);
                if(std::abs(Hij*ci[nc]) > cutoff) {
                  if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                  int sz = intm.size();
                  indx.clear();
                  for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                  indx[ia] = c;
                  indx[ib] = d;
                  std::sort(indx.begin(),indx.end());
                  if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                    for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                    nterms++;
                  }
                  if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                  sz = intm.size();
                  indx.clear();
                  for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                  indx[ia] = d;
                  indx[ib] = c; 
                  std::sort(indx.begin(),indx.end());
                  if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                    for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                    nterms++;
                  } 
                }   
              }   
            }     
          } 
        } 
      }  // states in occ_orbs
      Timer.stop("Generic");

      app_log()<<" Iteration: " <<i <<std::endl; 
      app_log()<<" Intermediate list has " <<nterms <<" terms" <<std::endl; 

      Timer.reset("Generic1");
      Timer.start("Generic1");
      
      bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,intm,nterms);
      if(!sucess) {
        app_error()<<" Error: Problems with diagonalizeTrialWavefunction. \n";
        return;
      }
      app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic3") <<std::endl;
      app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic4") <<std::endl;
      occ_orbs.reserve(intm.size());
      occ_orbs.clear();
      ci.clear();
      ci.reserve(nterms);
      for(int ii=0,nt=0; ii<nterms; ii++) {
        //app_log()<<"ci " <<ii <<" " <<eigVec(0,ii) <<std::endl; 
        if(std::abs(eigVec(0,ii)) > cutoff/10.0) {
          ci.push_back(eigVec(0,ii)); 
          for(int j=0; j<ne; j++)
            occ_orbs.push_back(intm[ii*ne+j]);
          nt++; 
        }    
      }
      Timer.stop("Generic1");

      std::ofstream out("iterativeCI.dat"); 
      if(out.fail()) {
        app_error()<<" Problems opening iterativeCI.dat \n";
        return;
      }
      {

        std::vector<std::tuple<double,int> > dets(ci.size());  
        for(int i=0; i<ci.size(); i++) dets[i] = std::make_tuple( std::abs(ci[i]),i); 
        std::sort( dets.begin(), dets.end(),  
        [] (const std::tuple<double,int>& a, const std::tuple<double,int>& b)
                 {return (std::get<0>(a)>std::get<0>(b));} );
        out<<" &FCI \n NCI = " <<ci.size() <<" \n /\n"; 
        for(int i=0; i<ci.size(); i++) {
          out<<ci[std::get<1>(dets[i])] <<" ";
          for(int j=0; j<ne; j++) out<<occ_orbs[ std::get<1>(dets[i])*ne+j]+1 <<" "; 
          out<<"\n";
        }
        out.close();
      }
 
      app_log()<<" Energy: " <<eigVal[0]+NuclearCoulombEnergy <<std::endl; 
      app_log()<<" Number of determinants after truncation: " <<occ_orbs.size()/ne <<std::endl;
      app_log()<<" Timings: " <<Timer.total("Generic") <<" " <<Timer.total("Generic1") <<std::endl;  

    } // iteration 

    if(diag_in_steps>0) {
      app_log()<<"\n***********************************************\n  #Determinants        Energy: " <<"\n";
      std::vector<std::tuple<double,int> > dets(ci.size());  
      for(int i=0; i<ci.size(); i++) dets[i] = std::make_tuple( std::abs(ci[i]),i); 
      std::sort( dets.begin(), dets.end(),  
        [] (const std::tuple<double,int>& a, const std::tuple<double,int>& b)
                 {return (std::get<0>(a)>std::get<0>(b));} );
      for(int i=1; i<ci.size(); i+=diag_in_steps) {
        intm.clear();
        for(int ki=0; ki<i; ki++) {
         int kk = std::get<1>(dets[ki]);  
         for(int kj=0; kj<ne; kj++) intm.push_back(occ_orbs[ kk*ne+kj]); 
        }
        bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,intm,i,false);
        if(!sucess) {
          app_error()<<" Error: Problems with diagonalizeTrialWavefunction. \n";
          return;
        }
        app_log()<<i <<" " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;  
      }
      app_log()<<"***********************************************" <<std::endl <<std::endl;
    }
 


    // print restart list 

}




// careful 
// FIX FIX FIX
// this routine returns interchanged (i,j)/(k,l), so it is wrong due to a std::complex conjugation for std::complex matrix elements 
int MultiPureSingleDeterminant::cntExcitations(std::vector<IndexType>& DL, std::vector<IndexType>& DR, IndexType& n0, IndexType& n1, IndexType& n2, IndexType& n3, std::vector<IndexType>& occ, RealType& sg)
{
  std::vector<IndexType>::iterator itR = DR.begin();
  std::vector<IndexType>::iterator itL = DL.begin();
  sg = 0.0;
  int cnt=0,pos=0,ind[20],cnt2=0,nq=0,cnt3=0;
  bool found;
  int dummy = 1000000;
  n0=n1=n2=n3=dummy;

  for(int i=0; i<NAEA; i++) {
   found=false;
   for(int j=0; j<NAEA; j++)
     if(*(itL+i) == *(itR+j)) {
       found = true;
       occ[cnt2++] = *(itL+i);
       *(itL+i) = dummy; 
       *(itR+j) = dummy; 
       break;
     }
   if(!found) {
     if(cnt<2) ind[cnt]=i;
     cnt++;
     if(cnt > 2) {
       sg=0.0;
       return 2*cnt;
     }
   }
  }
  for(int i=NAEA; i<NAEA+NAEB; i++) {
   found=false;
   for(int j=NAEA; j<NAEA+NAEB; j++)
     if(*(itL+i) == *(itR+j)) {
       found = true;
       occ[cnt2++] = *(itL+i);
       *(itL+i) = dummy;
       *(itR+j) = dummy;
       break;
     }
   if(!found) {
     if(cnt<2) ind[cnt]=i;
     cnt++;
     if(cnt > 2) {
       sg=0.0;
       return 2*cnt;
     }
   }
  }
  if(cnt == 1) {
    n1=static_cast<IndexType>(*(itL+ind[0]));
    for(int i=0; i<NAEA+NAEB; i++)
     if(*(itR+i) != dummy) {   // there should be only one
       nq = ind[0]-i;
       n0=static_cast<IndexType>(*(itR+i));
       break;
     }
    sg = nq%2==0?1.0:-1.0;
  } else if(cnt == 2)  {
    int iq1=-1,iq2=-1;
    n2=static_cast<IndexType>(*(itL+ind[0]));
    n3=static_cast<IndexType>(*(itL+ind[1]));
    for(int i=0; i<NAEA+NAEB; i++)
     if(*(itR+i) != dummy) {   
       n0=static_cast<IndexType>(*(itR+i));
       iq1=i;
       break;
     }
    for(int i=iq1+1; i<NAEA+NAEB; i++)
     if(*(itR+i) != dummy) {   // there should be only one
       n1=static_cast<IndexType>(*(itR+i)); 
       iq2=i;
       break;
     }
    if(iq1<0 || iq2<0) 
      APP_ABORT("Error in: cntExcitations.\n");
    nq = ind[0]-iq1+ind[1]-iq2;
    sg = nq%2==0?1.0:-1.0;
  } else
    sg=0.0;
  return 2*cnt;
}

}

      
