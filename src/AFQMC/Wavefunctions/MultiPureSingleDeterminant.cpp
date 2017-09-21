#include<cstdlib>
#include<complex>
#include<iostream>
#include<fstream>
#include<bitset>
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

    app_log()<<"\n\n --------------- Parsing MultiPureSD input ------------------ \n\n";

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
    std::string str3("no");
    filename = std::string("none");
    filetype = std::string("ascii");
    ParameterSet m_param;
    m_param.add(filename,"filename","std::string");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(runtype,"runtype","int");
    m_param.add(cutoff,"cutoff","double");
    m_param.add(str1,"diagHam","std::string");
    m_param.add(str2,"iterCI","std::string");
    m_param.add(str3,"fast","std::string");
    m_param.add(IterCI_maxit,"iterCI_it","int");
    m_param.add(IterCI_maxit,"iterCI_maxit","int");
    m_param.add(IterCI_cut,"iterCI_cut","double");
    m_param.add(diag_in_steps,"diag_steps","int");
    m_param.add(write_trial_density_matrix,"trial_density_matrix","std::string");
    m_param.put(cur);
  
    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int)) tolower);
    if(str1 == "no" || str1 == "false") diagHam = false; 

    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int)) tolower);
    if(str2 == "yes" || str2 == "true") iterCI = true; 

    std::transform(str3.begin(),str3.end(),str3.begin(),(int (*)(int)) tolower);
    if(str3 == "yes" || str3 == "true") fast_alg = true; 

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
  int nci;

  NCA = NCB = 0;
  bool fullMOMat = false;
  bool Cstyle = true;
  std::string type("occ");

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
        wfn_type = atoi((++it)->c_str());
        switch(wfn_type) {
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
            app_error()<<"Unknown wave-function type in MultiPureSingleDeterminant: " <<wfn_type <<std::endl;
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
        type = *(it+1);
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
    if(runtype==0) app_log()<<" Using low memory algorithm " <<std::endl;
    else app_log()<<" Using high memory algorithm (cpu cost ~ O[n^2*M^2]) " <<std::endl;
    occ_orbs.resize(nci*(NAEA+NAEB));
    std::vector<IndexType>::iterator it_occ = occ_orbs.begin();
    int nints = std::ceil(NMO/64.0); 
    for(int i=0; i<nci; i++) {

      // later on, apply cutoff
      in>>ci[i];
      if(in.fail()) {
         app_error()<<"Error with wfn file. CI #: " <<i <<std::endl;
         return false;
      }
      if(type == "occ") {
        for(int k=0; k<NAEA; k++,it_occ++) {
          in>>*it_occ;    
          (*it_occ)--;
          if(in.fail()) {
             app_error()<<"Error with wfn file. CI #: " <<i <<std::endl; 
             return false;
          }     
          if(*it_occ < 0 || *it_occ > NMO) {
             app_error()<<"Error with wfn file. Det definition (NAEA) # : "<<i <<" " <<*it_occ <<" " <<NAEA <<" " <<NMO <<std::endl; 
             return false;
          } 
        }
        for(int k=0; k<NAEB; k++,it_occ++) {
          in>>*it_occ;
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
      } else if(type == "binary") {
        // first alpha
        uint64_t a;
        int cnt=0; 
        for(int j=0; j<nints; j++) { 
          in>>a;
          if(in.fail()) {
             app_error()<<"Error with wfn file. CI #: " <<i <<std::endl;
             return false;
          } 
          std::bitset<64> bina(a);
          for(int ii=0, k=j*64; ii<64; ii++,k++) {
            if(k==NMO) break;
            if(bina[ii]) {
              if( ++cnt > NAEA ) {
                app_error()<<" Error with wfn file: CI #:" <<i <<std::endl;
                return false;  
              }
              *(it_occ++) = k;
            } 
          } 
        }    
        if( cnt != NAEA ) {
          app_error()<<" Error with wfn file (cnt!=NAEA): CI #:" <<i <<std::endl;
          return false;  
        }
        cnt=0;
        for(int j=0; j<nints; j++) {
          in>>a;
          if(in.fail()) {
             app_error()<<"Error with wfn file. CI #: " <<i <<std::endl;
             return false;
          }
          std::bitset<64> bina(a);
          for(int ii=0, k=j*64; ii<64; ii++,k++) {
            if(k==NMO) break;
            if(bina[ii]) {
              if( ++cnt > NAEB ) {
                app_error()<<" Error with wfn file: CI #:" <<i <<std::endl;
                return false;
              }
              *(it_occ++) = k+NMO;
            } 
          } 
        } 
        if( cnt != NAEB ) {
          app_error()<<" Error with wfn file (cnt!=NAEB): CI #:" <<i <<std::endl;
          return false;  
        }
      } else {
        APP_ABORT("ERROR: unknown wfn file type. \n\n\n ");
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
    if(wfn_type==1) orbsize*=2;
    if(wfn_type==2) orbsize = 2*NMO*(NAEA+NAEB); 
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

      if (wfn_type == 0 ) {

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

      } else if(wfn_type == 1) {

       if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat[ki*orbsize+i*NAEA+j] = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant, Alpha Block: " <<ki <<" " <<i <<" " <<j <<" \n";
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
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant, Beta Block: " <<ki <<" " <<i <<" " <<j <<" \n";
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
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant, Alpha Block: " <<ki <<" " <<i <<" " <<j <<" \n";
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
              app_error()<<"Problems reading ASCII file in MultiPureSingleDeterminant, Beta Block: " <<ki <<" " <<i <<" " <<j <<" \n";
              in.close();
              return false;
            }
          }
       }

      } else if(wfn_type==2) {
        APP_ABORT("Error: wfn_type==2 not uimplemented yet in MultiPureSD. \n\n\n");
      }
    } 
  }
  in.close();

  //trial_density_matrix.resize(NAEA+NAEB,NMO);
  temp_density_matrix.resize(2*NMO,NMO);
  temp_density_matrix_full.resize(2*NMO,NMO);
  SPtrial_density_matrix.resize(2*NMO,NMO);
  trial_density_matrix.resize(2*NMO,NMO);
  rank_updated_trial_density_matrix.resize(2*NMO,NMO);
  SPrank_updated_trial_density_matrix.resize(2*NMO,NMO);

  //mixed_density_matrix.resize(NAEA+NAEB,NMO);
  full_mixed_density_matrix.resize(2*NMO,NMO);
  mixed_density_matrix.resize(2*NMO,NMO);
  rank_updated_mixed_density_matrix.resize(2*NMO,NMO);

  // DM = QFull * BB0inv
  QFull.resize(NAEA+NAEB,NMO);

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

  // auxiliary matrices
  BB0inv_alpha.resize(NMO,NAEA);
  BB0inv_beta.resize(NMO,NAEB);
  tBB0inv_alpha.resize(NAEA,NMO);
  tBB0inv_beta.resize(NAEB,NMO);

  // larger than needed, but just in case
  Cwork.resize(2*NMO);
  pivot.resize(2*NMO);

  Iwork.resize(NAEA+NAEB);

  // calculate relations between determinants
  // setup necessary data structures for low-rank updates
  //
  
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
  if(NCA==NCB && NAEA == NAEB && rotated_hamiltonian && spinRestricted && wfn_type==0 ) 
    closed_shell = true;

  if(closed_shell) {
    app_log()<<"Found closed shell system. " <<std::endl;
  }

  // walker_type will come from WalkerHandler, but for now only ROHF/UHF is implemented correctly.
  walker_type=1;
  dm_type = std::max(wfn_type,1);
  if(closed_shell) {
    app_log()<<"Found closed shell system. " <<std::endl;
    dm_type=0;
  }

  // if orthogonal expansion, order the determinants by excitation number wrt reference
  if(!rotated_hamiltonian) {
    prepare_excitations();

    Kn.resize(maxEx*maxEx);
    SPKn.resize(maxEx*maxEx);
    VB0.resize(maxEx*NAEA); 
  }

  std::size_t m1=0, m2=0, m2c=0, m1c=0;

  if(runtype==-1) { // no local energy evaluation 

    app_log()<<" MultiPureSingleDeterminant::runtype == -1. No local energy evaluation. \n"; 
 
  } else if(runtype==0) { // only for pure dets

    if(!spinRestricted) 
      APP_ABORT("  Error: MultiPureSingleDeterminant save_memory=yes only works with spin restricted integrals.\n");

    assert(!rotated_hamiltonian);

    // force useFacHam in this case
    useFacHam=true;

    // idea: 1. construct Hamiltonian including all occupied orbitals.
    //       2. for each determinant, construct his own rowIndex vectors assuming non-continuous data
    //         2a. you can build the rowIndex vectors to ignore terms in the beginning and the end of the row that are not going to be included (e.g. ignore all jl when j is smaller than smallest occupied index.)
    //       3. perform sparse matrix vector using global data and column vectors and local rowIndex 

    std::map<IndexType,bool> all_orbs;
    all_orbs.clear();
    for(IndexType i=0; i<2*NMO; i++) all_orbs[i]=false;

    for(std::vector<IndexType>::iterator it=occ_orbs.begin(); it!=occ_orbs.end(); it++)
      all_orbs[*it] = true;
 
    app_log()<<" MultiPureSingleDeterminant - Creating Hamiltonian. \n";
    SMSpHijkl.resize(1);
    hij.resize(1);

    SMSpHijkl[0].setDims(2*NMO*NMO,2*NMO*NMO);
    SMSpHijkl[0].setup(head_of_nodes,name+std::string("SMSpHijkl_0"),TG.getNodeCommLocal());
    // since useFacHam==true, SMSpHijkl only has exchange component
    if(!sHam->createHamiltonianForPureDeterminant(dm_type,useFacHam,all_orbs,all_orbs,hij[0],SMSpHijkl[0],cutoff)) {
      app_error()<<"Error in createHamiltonianForPureDeterminant. \n";
      return false;
    }
    m1=hij[0].size();
    m2=SMSpHijkl[0].size();
    m1c=hij[0].size()*sizeof(s1D<ValueType>);
    m2c=SMSpHijkl[0].size()*(sizeof(SPComplexType)+2*sizeof(int));

    // setup bounds for SMSpHijkl
    prepare_hamiltonian_evaluation();

  } else {

    useFacHam=false;

    app_log()<<" MultiPureSingleDeterminant - Creating Hamiltonians. \n";
    
    if(rotated_hamiltonian) {
      haj.resize(ci.size()); 
      SMSpHabkl.resize(ci.size());  
    } else {
      hij.resize(ci.size()); 
      SMSpHijkl.resize(ci.size());  
    }
    ComplexMatrix Am(2*NMO,NAEA);
    int nr = (wfn_type==0)?NMO:2*NMO;
    std::map<IndexType,bool> isOcc_alpha; 
    std::map<IndexType,bool> isOcc_beta;
    for(int ki=0; ki<ci.size(); ki++) {
 
      if(rotated_hamiltonian) {
        for(int i=0; i<nr; i++)  
         for(int j=0; j<NAEA; j++)  
          Am(i,j) = OrbMat[ki*orbsize+i*NAEA+j]; 
        SMSpHabkl[ki].setDims(2*NMO*NMO,2*NMO*NMO);
        SMSpHabkl[ki].setup(head_of_nodes,name+std::string("SMSpHabkl_")+std::to_string(ki),TG.getNodeCommLocal());
        if(!sHam->createHamiltonianForGeneralDeterminant(dm_type,Am,haj[ki],SMSpHabkl[ki],cutoff)) {
          app_error()<<"Error in createHamiltonianForGeneralDeterminant. \n";
          return false;
        }
        m1+=haj[ki].size();
        m2+=SMSpHabkl[ki].size();
        m1c+=haj[ki].size()*sizeof(s1D<ComplexType>);
        m2c+=SMSpHabkl[ki].size()*(sizeof(SPComplexType)+2*sizeof(int));
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
        if(!sHam->createHamiltonianForPureDeterminant(closed?0:dm_type,useFacHam,isOcc_alpha,isOcc_beta,hij[ki],SMSpHijkl[ki],cutoff)) {
          app_error()<<"Error in createHamiltonianForPureDeterminant. \n";
          return false;
        }
        m1+=hij[ki].size();
        m2+=SMSpHijkl[ki].size();
        m1c+=hij[ki].size()*sizeof(s1D<ValueType>);
        m2c+=SMSpHijkl[ki].size()*(sizeof(SPComplexType)+2*sizeof(int));
      }
    
    }

  }

/*
 * This needs to be done for each determinant
  if(rotated_hamiltonian) {
    split_Ham_rows(SMSpHabkl.rows(),SMSpHabkl.rowIndex_begin(),ik0,ikN);
    pik0 = *(SMSpHabkl.row_index()+ik0);
  } else {
    split_Ham_rows(SMSpHijkl.rows(),SMSpHijkl.rowIndex_begin(),ik0,ikN);
    pik0 = *(SMSpHijkl.row_index()+ik0);
  }
*/

  if(closed_shell)
    local_buff.resize(NMO*NMO+1);
  else
    local_buff.resize(2*NMO*NMO+1);

  if(diagHam && runtype >= 0) {
    app_log()<<" Diagonalizing trial wave function in MultiPureSingleDeterminant. \n";
    std::vector<RealType> eigVal(ci.size());
    ComplexMatrix eigVec(1,ci.size());
    bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,occ_orbs,ci.size());
    if(sucess) {
      app_log()<<" New trial energy and ci coefficients: " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;
      for(int i=0; i<ci.size(); i++) app_log()<<i <<" old: " <<ci[i] <<" new: " <<eigVec(0,i) <<std::endl; 
      for(int i=0; i<ci.size(); i++) ci[i] = eigVec(0,i); 
      if(!rotated_hamiltonian) 
        for(int i=0; i<ci.size(); i++) ci_with_psign[i] = ci[i]*det_sign[i];
      app_log()<<std::endl; 
    } else {
      app_error()<<"Error diagonalizing trial wavefunction. \n";
      return false;
    }
  } 

  if(iterCI && runtype >= 0) {
    app_log()<<" Iterative CI: " <<std::endl;
    iterativeCI(IterCI_cut,1000,1000,IterCI_maxit);
    app_error()<<" Aborting after IterCI \n";
    return false;
  }

  ComplexType epot,ekin,o1,o2;
  HF.resize(2*NMO,NAEA);
  HF = ComplexType(0.0,0.0);  
  if(rotated_hamiltonian) {
    std::copy(OrbMat.data(),OrbMat.data()+NAEA*NMO,HF.data());
    if(wfn_type==0)
      std::copy(OrbMat.data(),OrbMat.data()+NAEA*NMO,HF.data()+NAEA*NMO);
    else
      std::copy(OrbMat.data()+NAEA*NMO,OrbMat.data()+2*NAEA*NMO,HF.data()+NAEA*NMO);
  } else {
    for(int i=0; i<NAEA; i++) HF(occ_orbs[i],i)=ComplexType(1.0,0.0);
    for(int i=0; i<NAEB; i++) HF(occ_orbs[i+NAEA],i)=ComplexType(1.0,0.0);
  }

  if(runtype > 0) {  // can't be called with runtype==0, since Spvn/Dvn has not been set yet
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
  }
  if(write_trial_density_matrix != "" && rank() == 0) {
    local_evaluateOneBodyTrialDensityMatrix();
    std::ofstream out(write_trial_density_matrix.c_str());
    out<<"# trial density matrix: NMO, NAEA, NAEB:" <<NMO <<" " <<NAEA <<" " <<NAEB <<"\n";
    for(int i=0; i<2*NMO; i++) {
      for(int j=0; j<NMO; j++)
        out<<SPtrial_density_matrix(i,j) <<" ";
      out<<"\n";
    }
    out.flush();
    out.close();
  }

  return true;
}

void MultiPureSingleDeterminant::prepare_excitations()
{
    // 0. on the 1st pass, only count excitations, maximum excitation and cnt of terms
    // On second pass
    // 1. count excitations and find {i,j,k,...}/{a,b,c,...}
    // 3. find permutation sign
    // 4. add iajb..., sign, ci_iakb, etc to lists 
    int cnt=0;
    maxEx=maxExa=maxExb=0;
    det_sign.reserve(ci.size());
    ci_with_psign.reserve(ci.size());
    map2unique.reserve(ci.size());
    std::vector<IndexType> loc(NAEA+NAEB);
    std::vector<IndexType> ik(NAEA+NAEB);
    std::vector<IndexType> ak(NAEA+NAEB);
    RealType psign;
    // find unique excitations
    std::vector<IndexType> uniq_a;
    std::vector<std::tuple<int,int,int>> cnter_uniq_a;// <0>=# excitations, <1>=# connections, <2>=pos
    std::vector<IndexType> uniq_b;
    std::vector<std::tuple<int,int,int>> cnter_uniq_b;
    IndexType* it_occ = occ_orbs.data();  
    IndexType* it_ref = occ_orbs.data()+ref*(NAEA+NAEB);  
    uniq_a.reserve(ci.size()*NAEA);
    uniq_b.reserve(ci.size()*NAEB);
    cnter_uniq_a.reserve(ci.size());
    cnter_uniq_b.reserve(ci.size());
    Iwork.resize(NAEA+NAEB);
    // count excitations and connections 
    for(int i=0; i<ci.size(); i++, it_occ+=NAEA+NAEB) {
      RealType psign=1.0;
      int dummy = countExct(NAEA+NAEB,it_ref,it_occ,true,loc.data(),ik.data(),ak.data(),psign); 
      ci_with_psign.push_back(psign*ci[i]);
      det_sign.push_back(psign);

      int pos = -1; // not using search on purpose  
      int nuq = uniq_a.size()/NAEA;
      IndexType* it = uniq_a.data();  
      for(int k=0; k<nuq; k++,it+=NAEA) {
        if( std::equal(it,it+NAEA,it_occ) ) {
          pos = k; 
          break;
        }
      }
      if(pos<0) {
        for(int k=0; k<NAEA; k++) 
          uniq_a.push_back(*(it_occ+k)); 
        // count excitations
        int nex=0;
        for(int i=0; i<NAEA; i++)
          if(!std::binary_search(it_ref,it_ref+NAEA,*(it_occ+i))) 
            nex++;
        if(nex > maxExa) maxExa=nex;
        cnter_uniq_a.push_back(std::make_tuple(nex,1,0));
      } else {
        std::get<1>(cnter_uniq_a[pos])++;
      }

      pos = -1; // not using search on purpose  
      nuq = uniq_b.size()/NAEB;
      it = uniq_b.data();  
      for(int k=0; k<nuq; k++,it+=NAEB) {
        if( std::equal(it,it+NAEB,it_occ+NAEA) ) {
          pos = k; 
          break;
        }
      }
      if(pos<0) {
        for(int k=0; k<NAEB; k++) 
          uniq_b.push_back(*(it_occ+NAEA+k)); 
        // count excitations
        int nex=0;
        for(int k=0; k<NAEB; k++)
          if(!std::binary_search(it_ref+NAEA,it_ref+NAEA+NAEB,*(it_occ+NAEA+k))) 
            nex++;
        if(nex > maxExb) maxExb=nex;
        cnter_uniq_b.push_back(std::make_tuple(nex,1,0));
      } else {
        std::get<1>(cnter_uniq_b[pos])++;
      }
    }
    maxEx = std::max(maxExa,maxExb);
    // data being stored
    // occ_string
    // # excitations
    // loc,i,a for all excitations
    // # connections
    // (index_ci, index_unique_beta/alpha) for each connection
    // size_per_unique: NAEA+3*iex+1+1+2*nconnections
    nunique_alpha = cnter_uniq_a.size(); 
    app_log()<<" Found " <<nunique_alpha <<" unique alpha determinants. \n";
    int nterms=0;
    for(int i=0; i<nunique_alpha; i++) {
      std::get<2>(cnter_uniq_a[i])=nterms;  // where this det begins
      nterms += NAEA + 2 + 3*std::get<0>(cnter_uniq_a[i]) + 2*std::get<1>(cnter_uniq_a[i]);
    }
    iajb_unique_alpha.resize(nterms);
    nunique_beta = cnter_uniq_b.size(); 
    app_log()<<" Found " <<nunique_beta <<" unique beta determinants. \n";
    nterms=0;
    for(int i=0; i<nunique_beta; i++) {
      std::get<2>(cnter_uniq_b[i])=nterms;  // where this det begins
      nterms += NAEB + 2 + 3*std::get<0>(cnter_uniq_b[i]) + 2*std::get<1>(cnter_uniq_b[i]);
    }
    iajb_unique_beta.resize(nterms);

    ovlp_unique_alpha.resize(nunique_alpha);
    ovlp_unique_beta.resize(nunique_beta);

    // now redo and store information
    it_occ = occ_orbs.data();  
    for(int i=0; i<ci.size(); i++, it_occ+=NAEA+NAEB) {
      int apos = -1; // not using search on purpose  
      IndexType* it = uniq_a.data();
      for(int k=0; k<nunique_alpha; k++,it+=NAEA) {
        if( std::equal(it,it+NAEA,it_occ) ) {
          apos = k;
          break;
        }
      }
      assert(apos>=0);
      int bpos = -1; // not using search on purpose  
      it = uniq_b.data();
      for(int k=0; k<nunique_beta; k++,it+=NAEB) {
        if( std::equal(it,it+NAEB,it_occ+NAEA) ) {
          bpos = k;
          break;
        }
      }
      assert(bpos>=0);

      map2unique.push_back(std::make_pair(apos,bpos));

      it = iajb_unique_alpha.data() + std::get<2>(cnter_uniq_a[apos]);
      int nc = *( it + NAEA + 1 + 3*std::get<0>(cnter_uniq_a[apos]) );
      if(nc == 0) { // first time
        RealType psign=1.0;
        int nex = countExct(NAEA,it_ref,it_occ,true,loc.data(),ik.data(),ak.data(),psign); 
        assert(nex==std::get<0>(cnter_uniq_a[apos]));
        for(int k=0; k<NAEA; k++)
          *(it+k) = *(it_ref+k);  
        for(int k=0; k<nex; k++) 
          *(it + loc[k] ) = ak[k];
        it+=NAEA;
        *(it++) = nex;
        for(int k=0; k<nex; k++) {
          *(it++) = loc[k];
          *(it++) = ik[k];
          *(it++) = ak[k];
        }
        *(it++) = 1; 
        *(it++) = i; // loc of ci coeff
        *(it++) = bpos; // loc of unique beta
      } else {
        it+=NAEA;
        it += 3*(*it)+1;
        (*it)++;
        it += 2*((*it) - 1) + 1;
        *(it++) = i; // loc of ci coeff
        *(it++) = bpos; // loc of unique beta 
      }

      it = iajb_unique_beta.data() + std::get<2>(cnter_uniq_b[bpos]);
      nc = *( it + NAEB + 1 + 3*std::get<0>(cnter_uniq_b[bpos]) );
      if(nc == 0) { // first time
        RealType psign=1.0;
        int nex = countExct(NAEB,it_ref+NAEA,it_occ+NAEA,true,loc.data(),ik.data(),ak.data(),psign); 
        assert(nex==std::get<0>(cnter_uniq_b[bpos]));
        for(int l=0; l<NAEB; l++)
          *(it+i) = *(it_ref+NAEA+l);  
        for(int k=0; k<nex; k++) 
          *(it + loc[k] ) = ak[k];
        it+=NAEB;
        *(it++) = nex;
        for(int k=0; k<nex; k++) {
          *(it++) = loc[k];   // [0:NAEB)
          *(it++) = ik[k];    // [NMO:2*NMO)
          *(it++) = ak[k];    // [NMO:2*NMO)
        }
        *(it++) = 1; 
        *(it++) = i; // loc of ci coeff
        *(it++) = apos; // loc of unique alpha
      } else {
        it+=NAEB;
        it += 3*(*it)+1;
        (*it)++;
        it += 2*((*it) - 1) + 1;
        *(it++) = i; // loc of ci coeff
        *(it++) = apos; // loc of unique alpha 
      }
    }
}


// define arrays and structures needed to evaluate the hamiltonian with runtype==0
void MultiPureSingleDeterminant::prepare_hamiltonian_evaluation()
{
  // matrix of boundaries: it is transposed with respect to order in SMSpHijkl
  // storing bounds for orbitals 0,NAEA,NAEA+1,...,NMO,END
  Hijkl_bounds.resize(NMO+2-NAEA,NMO*NMO);
  local_bounds.resize(2, (NAEA+std::max(maxExa,maxExb))*NMO);

  int nend = NMO+2-NAEA-1;
  SPComplexSMSpMat::int_iterator cols = SMSpHijkl[0].cols_begin(); 
  SPComplexSMSpMat::int_iterator indx = SMSpHijkl[0].rowIndex_begin(); 
  for(int i=0; i<NMO*NMO; i++,indx++) {
    int i0 = *indx; 
    int iM = *(indx+1); 
    Hijkl_bounds(0,i) = i0;  
    Hijkl_bounds(nend,i) = iM; 
    for(int n=NAEA; n<NMO; n++) {
      SPComplexSMSpMat::int_iterator p = std::lower_bound(cols+i0,cols+iM,n*NMO); 
      Hijkl_bounds(n-NAEA+1,i) = i0 + std::distance(cols+i0,p);
    }
  }
  // to use runtype==0 with packed storage format for Gia/Gib, you need to redefine the columns
  // of SMSpHijkl
  // CAREFUL when using in any other place!!!
  for(int i=0; i<NMO*NMO; i++) {
    int i0 = Hijkl_bounds(1,i);  // first term outside reference sector
    int iN = Hijkl_bounds(nend,i);  // last term for this row 
    for(int k=i0; k<iN; k++) {
      int c = *(cols+k);
      *(cols+k) = c%NMO;
    }
  }
}

// I'm having issues with setup order. Propg needs to be set up to define Spvn/Dvn here,
// which is not yet done when the wvefunction is set up.
// Calling this at every energy evaluation
void MultiPureSingleDeterminant::allocate_hamiltonian_evaluation()
{

  // reference sector of DMs of unique determinants. refG(ik,ndet)
  // this is mainly used for the sparse hamiltonian
  refG.resize(NAEA*NMO,nunique_alpha+nunique_beta); 

  // packed storage for full DM: Gia(ndet, ik) 
  Gia.resize(nunique_alpha,NMO*(NAEA+maxExa)); 
  Gib.resize(nunique_beta,NMO*(NAEB+maxExb)); 

  // storage for vn*G for all unique determinants
  vb.resize((nunique_alpha+nunique_beta)*nCholVecs);
  vb_helper.resize((nunique_alpha+nunique_beta)*nCholVecs);

  // storage for Muv * Gv = Piu
  refP.resize(NAEA*NMO,(nunique_alpha+nunique_beta));

  // packed storage for Pia 
  Pia.resize(nunique_alpha,NMO*(NAEA+maxExa));   
  Pib.resize(nunique_beta,NMO*(NAEB+maxExb));   

  G0.resize(maxEx,NMO);
}

bool MultiPureSingleDeterminant::hdf_write(hdf_archive& read, const std::string& tag, bool include_tensors) 
{
    return true;
}

bool MultiPureSingleDeterminant::hdf_write()
{
    return true;
}

// on exit, commBuff contains the green function information for all walkers in appropriate place   
// // buffer: (this routine assumes the buffer is large enough) 
// //    wlksz: size of the data of a walker
// //    gfoffset: location of the gf block within a walker
// //    transposed: if true: walker index along columns
// //                if false: walker index along rows (continuous walker data)
void MultiPureSingleDeterminant::evaluateOneBodyMixedDensityMatrix(WalkerHandlerBase* wset, SPComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full)
{

  ComplexType oa,ob;
  int nw = wset->numWalkers(true), cnt=0;
  int nw0 = wset->numWalkers(false);

  // right now there is a problem here, 
  // is called from Propagator, full refers to whether save_memory is true or not.
  // but this routine can be called for both energy or vbias, I still want the option to do either
  // true/false in full
  assert(full);

  int sz = 1 + NMO*NMO;
  if(!closed_shell) sz += NMO*NMO;
  int wstride = transposed?nw0:1;
  int ax = transposed?1:wlksz;

  for(int i=0; i<nw; i++) {
    if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
    if(cnt%ncores_per_TG == core_rank ) {
      local_evaluateOneBodyMixedDensityMatrixFull(wset->getSM(i),oa,full_mixed_density_matrix,full);

// Do I really need a lock?
      {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buf->getMutex()));
        //BLAS::copy(sz,local_buff.data(),1,buf->values() + ax*cnt + wstride*gfoffset ,wstride);
        *(buf->values() + ax*cnt + wstride*gfoffset) = oa;
        BLAS::copy(sz-1,full_mixed_density_matrix.data(),1,buf->values()+ax*cnt+wstride*(gfoffset+1),wstride);
      }
    }
    ++cnt;
  }

}

/*
 * This routine calculates the mixed Density Matrix of the multideterminant trial wavefunction 
 *
 *   DM(i,k) = <PsiT| c+_i ck |SM> / <PsiT|SM>
 *           = sum_n ci*(n) <Dn| c+_i ck |SM> / sum_n ci*(n) <Dn|SM>
 *   For rotated_hamiltonian (or fast_alg=false), this sum is implemented directly.
 *   For orthogonal expansions and fast_alg=true, the following is used:
 *
 *   Woodbury identity:
 *   (A + UV)^-1 = A^-1 + A^-1 * U * (1 + V * A^-1 * U)^-1 * V * A^-1         
 *   Since for pure determinants we can write DM_i^T = SM * (SM_i)^-1 = SM * (SM_0)^-1 * Qi  
 *   Where Qi = 1 + U_i * (1 + V_i * (SM_0)^-1 * U_i)^-1 * V_i * (SM_0)^-1
 *   Then:
 *
 *   DM * <PsiT|SM> = sum_n ci*(n) <Dn|SM> DM_n = sum_n ci*(n) <Dn|SM> Q_n^T * DM_0 
 *                  = QFull * DM_0
 *
 *   QFull = sum_n ci*(n) <Dn|SM> Q_n^T               
 *
 */
void MultiPureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrixFull(const ComplexType* SlaterMat, ComplexType& ovl, SPComplexMatrix& dm, bool full)
{
  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  ComplexType o1, ovl_alpha,ovl_beta;
  ovl = ComplexType(0.0);

//  fast_alg = false;
//  test_cnter++;
//  if(test_cnter==20) fast_alg=true;

  if(fast_alg && !rotated_hamiltonian) { // algorithm using low-rank update formulas

    if(!full) {
      APP_ABORT(" Error: Calling local_evaluateOneBodyMixedDensityMatrixFull with fast_alg and full=false \n");
    }

    // calculate all overlaps first
    calculate_unique_overlaps(o1,ovl,SlaterMat);

    calculate_QFull(o1,true,false,false);
    
    // now DM = QFull^T * BB0inv^T 
    ComplexType invo = ComplexType(1.0,0.0)/ovl;
    DenseMatrixOperators::product_AtBt(NMO,NMO,NAEA,invo,QFull.data(),NMO,BB0inv_alpha.data(),NAEA,zero,dm.data(),NMO);
    DenseMatrixOperators::product_AtBt(NMO,NMO,NAEB,invo,QFull.data()+NAEA*NMO,NMO,BB0inv_beta.data(),NAEB,zero,dm.data()+NMO*NMO,NMO);

    // comment to test fast algorithm
    return;

    // testing 
    ComplexType ovl_ = ovl;
    ComplexMatrix DM_(dm);

    Timer.reset("Generic1");
    Timer.start("Generic1");

    local_evaluateOneBodyMixedDensityMatrix(0,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,full);
    o1 = myconj(ci[0])*overlaps[0]*overlaps[1];
    dm = static_cast<SPComplexType>(o1)*mixed_density_matrix;
    ovl = o1;

    for(int i=1; i<ci.size(); i++)
    {
      local_evaluateOneBodyMixedDensityMatrix(i,SlaterMat,overlaps[2*i],overlaps[2*i+1],mixed_density_matrix,full);
      o1 = myconj(ci[i])*overlaps[2*i]*overlaps[2*i+1];
      dm += static_cast<SPComplexType>(o1)*mixed_density_matrix;
      ovl += o1;
    }
    dm *= (SPComplexType(1.0,0.0)/static_cast<SPComplexType>(ovl)); 

    Timer.stop("Generic1");
    app_log()<<" Time for current alg: " <<Timer.average("Generic1") <<std::endl;
  
    app_log()<<"\n\n\n" <<" ovlps: " <<ovl_ <<"  " <<ovl <<std::endl;
    ComplexMatrix DM2_(DM_);
    DM2_ = DM2_ - dm;  
    app_log()<<"DM difference: (alpha) \n";
    for(int i=0; i<NMO; i++) 
     for(int j=0; j<NMO; j++) 
       if( std::abs(DM2_(i,j)) > 1e-8 ) app_log()<<i <<" " <<j <<" " <<std::abs(DM2_(i,j))
                                                 <<" " <<DM_(i,j) <<" " <<dm(i,j) <<std::endl;
    app_log()<<"DM difference: (beta) \n";
    for(int i=NMO; i<2*NMO; i++) 
     for(int j=0; j<NMO; j++) 
       if( std::abs(DM2_(i,j)) > 1e-8 ) app_log()<<i-NMO <<" " <<j <<" " <<std::abs(DM2_(i,j)) 
                                                 <<" " <<DM_(i,j) <<" " <<dm(i,j) <<std::endl;

    //APP_ABORT(" TESTING \n\n\n");

  } else {

    // accumulating the DM in single precision, is this good enough???
    local_evaluateOneBodyMixedDensityMatrix(0,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,full);
    o1 = myconj(ci[0])*overlaps[0]*overlaps[1];
    dm = static_cast<SPComplexType>(o1)*mixed_density_matrix;
    ovl = o1; 

    for(int i=1; i<ci.size(); i++)
    { 

      local_evaluateOneBodyMixedDensityMatrix(i,SlaterMat,overlaps[2*i],overlaps[2*i+1],mixed_density_matrix,full);
      o1 = myconj(ci[i])*overlaps[2*i]*overlaps[2*i+1];
      dm += static_cast<SPComplexType>(o1)*mixed_density_matrix;
      ovl += o1; 
    }
    dm *= (SPComplexType(1.0,0.0)/static_cast<SPComplexType>(ovl)); 

  }
}

// on input, ovl should be the overlap of the reference determinant.
// on output, it is the overlap of the full wfn
void MultiPureSingleDeterminant::calculate_unique_overlaps(ComplexType& ovl_ref, ComplexType& ovl, const ComplexType* SlaterMat)
{

  // compute SM * (SM_0)^-1. Reference determinant is always the first
  // copy rows corresponding to occupied orbitals to S0
  // S0(i,:) = SlaterMat(ik,:), where ik is the ith occupied orbital
  ComplexType* itS0 = S0.data();
  for(IndexType* it=occ_orbs.data()+ref*(NAEA+NAEB); it!=occ_orbs.data()+NAEA+ref*(NAEA+NAEB); it++, itS0+=NAEA)
    std::copy(SlaterMat+(*it)*NAEA, SlaterMat+((*it)+1)*NAEA ,itS0);

  ovl_ref = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

  // BB0inv_alpha = SlaterMat * S0
  DenseMatrixOperators::product(NMO,NAEA,NAEA,ComplexType(1),SlaterMat,NAEA,S0.data(),NAEA,ComplexType(0),BB0inv_alpha.data(),NAEA);

  // copy rows corresponding to occupied orbitals to S0    
  // S0(i,:) = SlaterMat(ik,:), where ik is the ith occupied orbital
  itS0 = S0.data();
  for(IndexType* it=occ_orbs.data()+ref*(NAEA+NAEB)+NAEA; 
                 it!=occ_orbs.data()+(ref+1)*(NAEB+NAEA); it++, itS0+=NAEB)
    std::copy(SlaterMat+(*it)*NAEA, SlaterMat+(*it)*NAEA+NAEB ,itS0);

  ovl_ref *= Invert(S0.data(), NAEB, NAEB, Cwork.data(),pivot.data());

  // BB0inv_beta = SlaterMat * S0
  DenseMatrixOperators::product(NMO,NAEB,NAEB,ComplexType(1),SlaterMat+NAEA*NMO,NAEA,S0.data(),NAEB,ComplexType(0),BB0inv_beta.data(),NAEB);

  /* iajb_unique_alpha:
   * 0:NAEA-1: occ string
   * NAEA: # excitations = iext
   * NAEA+1:NAEA+1+3*iext: for each excitation: loc,i,a
   * NAEA+1+3*iext+1: # connected dets
   * NAEA+1+3*iext+2:...: (index of det in ci_with_sign,index of unique beta) for all connected dets  
   */ 
  std::vector<IndexType>::iterator it = iajb_unique_alpha.begin();
  std::vector<ComplexType>::iterator ovlp = ovlp_unique_alpha.begin(); 
  for(int n=0; n<nunique_alpha; n++) {
    ComplexType* ptr = Kn.data();
    it += NAEA;
    int iext = *(it++);
    for(int i=0; i<iext; i++) {
      IndexType ik = *(it+3*i+2);  // new orbital 
      for(int k=0; k<iext; k++, ptr++) 
        *ptr = BB0inv_alpha(ik,*(it+3*k));
    }
    if(iext==0)
      *(ovlp++) = ComplexType(1.0,0.0);
    else
      *(ovlp++) = DenseMatrixOperators::DeterminantSmall(Kn.data(), iext, pivot.data());
    it += (3*iext);
    it += 2*(*it) + 1;  
  }

  ovl = ComplexType(0);
  it = iajb_unique_beta.begin();
  ovlp = ovlp_unique_beta.begin();
  for(int n=0; n<nunique_beta; n++, ovlp++) {
    ComplexType* ptr = Kn.data();
    it += NAEB;
    int iext = *(it++);
    for(int i=0; i<iext; i++) {
      IndexType ik = *(it+3*i+2)-NMO;  // new orbital 
      for(int k=0; k<iext; k++, ptr++)  
        *ptr = BB0inv_beta(ik,*(it+3*k));  
    }
    if(iext==0)
      *(ovlp) = ComplexType(1.0,0.0);
    else
      *(ovlp) = DenseMatrixOperators::DeterminantSmall(Kn.data(), iext, pivot.data());
    it += (3*iext); 
    ComplexType ovln = ComplexType(0);
    int nc = *(it++);
    for(int i=0; i<nc; i++, it+=2)
      ovln += myconj(ci_with_psign[*it]) * ovl_ref * (*ovlp) * ovlp_unique_alpha[*(it+1)];
    ovl += ovln;
  }

}

// The transpose of QFull is actually assembled
void MultiPureSingleDeterminant::calculate_QFull(ComplexType ovl_ref, bool getQFull, bool getGa, bool getGb)
{
  if(getQFull) 
    QFull = ComplexType(0);
  ComplexType dummy; 

  std::vector<IndexType>::iterator it = iajb_unique_alpha.begin();
  std::vector<ComplexType>::iterator ovlp = ovlp_unique_alpha.begin();
  for(int n=0; n<nunique_alpha; n++, ovlp++) {
    if(std::abs(*ovlp)<1e-12) {
      it += NAEA;
      it += 3*(*it)+1;
      it += 2*(*it) + 1;
      if(getGa) {
        std::fill(Gia[n],Gia[n+1],ComplexType(0));
        for(int k=0; k<NAEA*NMO; k++)
          refG(k,n) = ComplexType(0);
      }
      continue;
    }
    ComplexType* ptr = Kn.data();
    ComplexType* ptr_VB0 = VB0.data();  
    std::vector<IndexType>::iterator itn = it+NAEA;
    int iext = *(itn++);
    for(int i=0; i<iext; i++, ptr_VB0+=NAEA) {
      IndexType ik = *(itn+3*i+2);  // new orbital 
      std::copy(BB0inv_alpha.data()+NAEA*ik, BB0inv_alpha.data()+NAEA*(ik+1), ptr_VB0); 
      *(ptr_VB0 + *(itn+3*i)) -= ComplexType(1.0,0.0); 
      for(int k=0; k<iext; k++, ptr++)
        *ptr = BB0inv_alpha(ik,*(itn+3*k));
    }
    if(iext!=0) {
      dummy = DenseMatrixOperators::InvertSmall(Kn.data(), iext, Cwork.data(), pivot.data());
      DenseMatrixOperators::product(iext,NAEA,iext,ComplexType(1,0),Kn.data(),iext,VB0.data(),NAEA,ComplexType(0,0),SS0.data(),NAEA);
    }
    itn += (3*iext);
    int nc = *(itn++);
    ComplexType ovln = ComplexType(0);
    for(int i=0; i<nc; i++, itn+=2)  
      ovln += myconj(ci_with_psign[*itn]) * ovl_ref * (*ovlp) * ovlp_unique_beta[*(itn+1)]; 

    // QFull(i,k) = sum_n cn *ovla(n)*ovlb(n)*An^*(i,a)*[1-SS0^T*Un^T](a,k)
    //            --> QFull( lok[k], occ[0:NE] ) -= ovln*SS0(k,0:NE)   // since this is actually QFull^T  
    //            --> QFull( 0:NE , occ[0:NE] ) += ovln                           // since this is actually QFull^T  
    if(getQFull) {
      if(iext!=0) {
        std::vector<IndexType>::iterator ite = it+NAEA+1;               
        for(int k=0; k<iext; k++, ite+=3) {
          IndexType loc = *ite;
          for(int i=0; i<NAEA; i++)  {           
            QFull(loc,*(it+i)) -= ovln*SS0(k,i); 
          }
        }
      }
      for(int i=0; i<NAEA; i++)             
        QFull(i,*(it+i)) += ovln; 
    }
    // Ga does not include factor of conjg(An) and is ordered in special format to be used later
    // Ga = (1 - Un*SS0)^T * BB0^T = BB0^T - SS0^T * (Un^T*BB0^T), 
    // (Un^T*BB0^T)(k,:) = (BB0(:,loc[k]))^T 
    if(getGa) {
      ComplexType* G = Gia[n];  // G = &Gia(n,0)
      // Ga: Using G0 for temporary storage  
      // G0 = (Un^T*BB0^T)
      std::vector<IndexType>::iterator ite = it+NAEA+1;
      for(int i=0; i<iext; i++) {
        int loc = (*(ite+3*i))*NMO;
        std::copy(tBB0inv_alpha.data()+loc,tBB0inv_alpha.data()+loc+NMO,G0.data()+i*NMO);
      }
      // Ga = BB0^T - SS0^T * G0
      std::copy(tBB0inv_alpha.data(),tBB0inv_alpha.data()+NAEA*NMO,G);   
      DenseMatrixOperators::product_AtB(NAEA,NMO,iext,ComplexType(-1.0,0),SS0.data(),NAEA,G0.data(),NMO,ComplexType(1,0),G,NMO);      
      // Ga in packed storage pushed excited orbitals immediately after reference sector
      // Ga( loc[k], :) => Ga( NAEA+k, :) && Ga( loc[k], :)=0 
      for(int i=0; i<iext; i++) {
        int loc = (*(ite+3*i))*NMO;
        std::copy(G+loc,G+loc+NMO,G+(NAEA+i)*NMO);
        std::fill(G+loc,G+loc+NMO,ComplexType(0));
      }
      // right now, also generating transposed copy of reference sector
      G = Gia[n];   
      for(int k=0; k<NAEA*NMO; k++, G++)  
        refG(k,n) = *G;
    }
    // now push it forward
    it = itn;
  }

  it = iajb_unique_beta.begin();
  ovlp = ovlp_unique_beta.begin();
  for(int n=0; n<nunique_beta; n++, ovlp++) {
    if(std::abs(*ovlp)<1e-12) {
      it += NAEB;
      it += 3*(*it)+1;
      it += 2*(*it) + 1;
      if(getGb) {
        std::fill(Gib[n],Gib[n+1],ComplexType(0));
        for(int k=0; k<NAEB*NMO; k++)
          refG(k,n+nunique_alpha) = ComplexType(0);
      }
      continue;
    }
    ComplexType* ptr = Kn.data();
    ComplexType* ptr_VB0 = VB0.data();  
    std::vector<IndexType>::iterator itn = it+NAEB;
    int iext = *(itn++);
    for(int i=0; i<iext; i++, ptr_VB0+=NAEB) {
      IndexType ik = *(itn+3*i+2)-NMO;  // new orbital 
      std::copy(BB0inv_beta.data()+NAEB*ik, BB0inv_beta.data()+NAEB*(ik+1), ptr_VB0); 
      *(ptr_VB0 + *(itn+3*i)) -= ComplexType(1.0,0.0); 
      for(int k=0; k<iext; k++, ptr++)
        *ptr = BB0inv_beta(ik,*(itn+3*k));
    }
    if(iext!=0) {
      dummy = DenseMatrixOperators::InvertSmall(Kn.data(), iext, Cwork.data(), pivot.data());
      DenseMatrixOperators::product(iext,NAEB,iext,ComplexType(1,0),Kn.data(),iext,VB0.data(),NAEB,ComplexType(0,0),SS0.data(),NAEA);
    }
    itn += (3*iext);
    ComplexType ovln = ComplexType(0);
    int nc = *(itn++);
    for(int i=0; i<nc; i++, itn+=2) 
      ovln += myconj(ci_with_psign[*itn]) * ovl_ref * (*ovlp) * ovlp_unique_alpha[*(itn+1)]; 

    // beta section goes in QFull(NAEA+:,:)
    // QFull(i,k) = sum_n cn *ovla(n)*ovlb(n)*An^*(i,a)*[1-SS0^T*Un^T](a,k)
    //            --> QFull( lok[k], occ[0:NE] ) -= ovln*SS0(k,0:NE)   // since this is actually QFull^T  
    //            --> QFull( 0:NE , occ[0:NE] ) += ovln                           // since this is actually QFull^T  
    if(getQFull) {
      if(iext!=0) {
        std::vector<IndexType>::iterator ite = it+NAEB+1;               
        for(int k=0; k<iext; k++, ite+=3) {
          IndexType loc = *ite;
          for(int i=0; i<NAEB; i++)             
            QFull(NAEA+loc,*(it+i)-NMO) -= ovln*SS0(k,i); 
        }
      }
      for(int i=0; i<NAEB; i++)             
        QFull(NAEA+i,*(it+i)-NMO) += ovln; 
    }
    // Gb does not include factor of conjg(An) and is ordered in special format to be used later
    // Gb = (1 - Un*SS0)^T * BB0^T = BB0^T - SS0^T * (Un^T*BB0^T), 
    // (Un^T*BB0^T)(k,:) = (BB0(:,loc[k]))^T 
    if(getGb) {
      ComplexType* G = Gib[n];  // G = &Gib(n,0)
      // Gb: Using G0 for temporary storage  
      // G0 = (Un^T*BB0^T)
      std::vector<IndexType>::iterator ite = it+NAEB+1;
      for(int i=0; i<iext; i++) {
        int loc = (*(ite+3*i))*NMO;
        std::copy(tBB0inv_beta.data()+loc,tBB0inv_beta.data()+loc+NMO,G0.data()+i*NMO);
      }
      // Gb = BB0^T - SS0^T * G0
      std::copy(tBB0inv_beta.data(),tBB0inv_beta.data()+NAEB*NMO,G);   
      DenseMatrixOperators::product_AtB(NAEB,NMO,iext,ComplexType(-1.0,0),SS0.data(),NAEB,G0.data(),NMO,ComplexType(1,0),G,NMO);      
      // Gb in packed storage pushed excited orbitals immediately after reference sector
      // Gb( loc[k], :) => Gb( NAEB+k, :) && Gb( loc[k], :)=0 
      for(int i=0; i<iext; i++) {
        int loc = (*(ite+3*i))*NMO;
        std::copy(G+loc,G+loc+NMO,G+(NAEB+i)*NMO);
        std::fill(G+loc,G+loc+NMO,ComplexType(0));
      }
      // right now, also generating transposed copy of reference sector
      G = Gib[n];   
      for(int k=0; k<NAEB*NMO; k++, G++)  
        refG(k,n+nunique_alpha) = *G;
    }
    // now push it forward
    it = itn;
  }
}

// might need to save inverse(S0) for fast update formula, don't know yet
void MultiPureSingleDeterminant::local_evaluateOneBodyMixedDensityMatrix(int detn, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, SPComplexMatrix& dm, bool full)
{

  // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
  // lots of simplification because A is an identity matrix with possibly exchanged columns, 
  // look at derivations for information

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  IndexType* it_orbs = occ_orbs.data()+(NAEA+NAEB)*detn;
  ComplexType *Orbs = OrbMat.data()+detn*orbsize;
  ComplexType *Orbs_beta = Orbs+((wfn_type==0)?0:NAEA*NMO);   

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
  DenseMatrixOperators::product(NMO,NAEA,NAEA,one,SlaterMat,NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);


  // copy to dm 
  if(rotated_hamiltonian) {

    if(full) {
#if defined(AFQMC_SP)
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,Orbs,NAEA,zero,temp_density_matrix.data(),NMO);
      std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, dm.data());
#else
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,Orbs,NAEA,zero,dm.data(),NMO);
#endif
      DenseMatrixOperators::transpose(NMO,dm.data(),NMO);
    } else {
      std::fill(dm.begin(), dm.begin()+NMO*NMO, 0);
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        dm(i,j) = SS0(j,i);
    }     

  } else { 
    std::fill(dm.begin(), dm.begin()+NMO*NMO, 0);
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
  DenseMatrixOperators::product(NMO,NAEB,NAEB,one,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

  // copy to dm
  if(rotated_hamiltonian) {
    if(full) {
#if defined(AFQMC_SP)
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,Orbs_beta,NAEA,zero,temp_density_matrix.data(),NMO);
      std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, dm.data()+NMO*NMO);
#else
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,Orbs_beta,NAEA,zero,dm.data()+NMO*NMO,NMO);
#endif
      DenseMatrixOperators::transpose(NMO,dm.data()+NMO*NMO,NMO);
    } else {
      std::fill(dm.begin()+NMO*NMO, dm.begin()+2*NMO*NMO, 0);
      for(int i=0; i<NAEB; i++)
       for(int j=0; j<NMO; j++)
        dm(i+NMO,j) = SS0(j+NMO,i);
    }
  } else {
    std::fill(dm.begin()+NMO*NMO, dm.begin()+2*NMO*NMO, 0);
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
    ComplexType *Orbs_beta = Orbs+((wfn_type==0)?0:NAEA*NMO);
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

// MAM: Changed the meaning of this function to be the computation of the full trial density matrix
void MultiPureSingleDeterminant::local_evaluateOneBodyTrialDensityMatrix()
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
    SPtrial_density_matrix = ComplexType(0.0);
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfn_type==0)?0:NAEA*NMO;

    for(int i=0; i<ci.size(); i++) {
     // i==j
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEA; j++)
       Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEB; j++)
       Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];

     local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, SPrank_updated_trial_density_matrix,true);
     o1 = myconj(ci[i])*ci[i]*oa*ob;
     SPtrial_density_matrix += static_cast<SPComplexType>(o1)*SPrank_updated_trial_density_matrix; 
     ovl += o1 ;
     for(int j=i+1; j<ci.size(); j++)
     {
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEA; jj++)
         Am(ii,jj) = OrbMat[j*orbsize+ii*NAEA+jj];
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEB; jj++)
         Am(ii+NMO,jj) = OrbMat[j*orbsize+ii*NAEA+jj+del];
       local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, SPrank_updated_trial_density_matrix,true);
       o1 = myconj(ci[i])*ci[j]*oa*ob; 
       SPtrial_density_matrix += static_cast<SPComplexType>(o1)*SPrank_updated_trial_density_matrix; 
       ovl += o1 ;
       o1 = ci[i]*myconj(ci[j])*myconj(oa*ob);
       // is there an operation for std::complex conjugation???
       for(int ii=0; ii<NMO; ii++)  
        for(int jj=0; jj<NMO; jj++) { 
         SPtrial_density_matrix(ii,jj) += static_cast<SPComplexType>(o1)*myconj(SPrank_updated_trial_density_matrix(jj,ii));
         SPtrial_density_matrix(ii+NMO,jj) += static_cast<SPComplexType>(o1)*myconj(SPrank_updated_trial_density_matrix(jj+NMO,ii));
        }
       ovl += o1 ;
     }
    }
    SPtrial_density_matrix *= (SPComplexType(1.0,0.0)/static_cast<SPComplexType>(ovl));

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
    SPtrial_density_matrix = trial_density_matrix;
  }

/*
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
*/
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

  if(fast_alg && !rotated_hamiltonian) { // algorithm using low-rank update formulas

    // calculate unique overlaps and  total overlap
    calculate_unique_overlaps(ovl_beta,ovl_alpha,SlaterMat);
    ovl_beta = ComplexType(1.0,0.0);

    // comment to test fast algorithm
    return;

    ComplexType tmp=ovl_alpha;

    Timer.reset("Generic1");
    Timer.start("Generic1");

    ovl_beta = ComplexType(1.0,0.0);
    ovl_alpha = ComplexType(0.0,0.0);

    for(int i=0; i<ci.size(); i++)
    {
      ovl_alpha += myconj(ci[i]) * local_evaluateOverlapSlaterDet(i,SlaterMat);
    }

    Timer.stop("Generic1");
    app_log()<<" Time for ovlp current alg: " <<Timer.average("Generic1") <<std::endl;

    app_log()<<" Ovlp: (new/old) " <<tmp <<" " <<ovl_alpha <<std::endl;
    
    //APP_ABORT(" TESTING \n\n\n");

  } else {

    ovl_beta = ComplexType(1.0,0.0);
    ovl_alpha = ComplexType(0.0,0.0);

    for(int i=0; i<ci.size(); i++)
    {
      ovl_alpha += myconj(ci[i]) * local_evaluateOverlapSlaterDet(i,SlaterMat); 
    }

  }

}

  void MultiPureSingleDeterminant::evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

    const SPComplexType one = SPComplexType(1.0,0.0);
    const SPComplexType zero = SPComplexType(0.0,0.0);
    ComplexType lpe=ComplexType(0.0,0.0),lke=ComplexType(0.0,0.0);
    ComplexType deno=ComplexType(0.0,0.0), ovl_ref;
    ComplexType ke_nume=ComplexType(0.0,0.0), pe_nume=ComplexType(0.0,0.0);
    IndexType ef = 0, NAE = NAEA+NAEB, NAE2 = (NAEA+NAEB)*(NAEA+NAEB);
    int nr1=1, nc1=2*NMO*NMO;

    if(runtype==-1) {
      ekin=epot=ComplexType(0);
      ovl_alpha=ovl_beta=ComplexType(1);
      return;
    } 

#ifdef AFQMC_TIMER
    Timer.start("MultiPureSingleDeterminant:evaluateLocalEnergy");
#endif

    // clean arrays 
    for(std::vector<ComplexType>::iterator it=pe.begin(); it!=pe.end(); ++it) *it=0.0;
    for(std::vector<ComplexType>::iterator it=ke.begin(); it!=ke.end(); ++it) *it=0.0;

    if(runtype == 0) {

#if defined(QMC_COMPLEX)
APP_ABORT(" Error: runtype=0 not yet implemented for complex integrals. Need to fix L^+ issue. \n\n\n");
#endif

      // only needed once, hopefuly doesn't waste time
      allocate_hamiltonian_evaluation();
     
      // calculate all overlaps first
      calculate_unique_overlaps(ovl_ref,deno,SlaterMat);

      // generate transpose of BB0inv_alpha, BB0inv_beta 
      for(int i=0; i<NAEA; i++)
        BLAS::copy(NMO,BB0inv_alpha.data()+i,NAEA,tBB0inv_alpha.data()+i*NMO,1);       
      for(int i=0; i<NAEB; i++)
        BLAS::copy(NMO,BB0inv_beta.data()+i,NAEB,tBB0inv_beta.data()+i*NMO,1);       

      // calculate QFull and green functions for unique determinants
      calculate_QFull(ovl_ref,true,true,true);

      // kinetic energy using full DM
      // NOTE: temp_density_matrix_full is missing factor of 1/deno on purpose
      DenseMatrixOperators::product_AtBt(NMO,NMO,NAEA,ComplexType(1.0,0.0),QFull.data(),NMO,BB0inv_alpha.data(),NAEA,zero,temp_density_matrix_full.data(),NMO);
      DenseMatrixOperators::product_AtBt(NMO,NMO,NAEB,ComplexType(1.0,0.0),QFull.data()+NAEA*NMO,NMO,BB0inv_beta.data(),NAEB,zero,temp_density_matrix_full.data()+NMO*NMO,NMO);
      ComplexMatrix::iterator itG = temp_density_matrix_full.begin();
      s1Dit end1 = hij[0].end();
      for(s1Dit it = hij[0].begin(); it != end1; it++) 
        ke_nume += *(itG + std::get<0>(*it)) * std::get<1>(*it);

      // right now, useFacHam is forced with runtype == 0
      // calculate vb for all unique determinants
      
      // reference sector first: only correct if reference is formed by first NAEA states
      if(sparse_vn) // Dvn^T * refG 
        SparseMatrixOperators::product_SpMatTM( int(NAEA*NMO), refG.cols() , Spvn->cols(), SPValueType(1.0), Spvn->values(), Spvn->column_data(), Spvn->row_index(), refG.data() , refG.cols(), SPValueType(0.0), vb_helper.data(), refG.cols());
      else 
        DenseMatrixOperators::product_AtB(Dvn->cols(),refG.cols(),int(NAEA*NMO), SPValueType(1.0), Dvn->values(), Dvn->cols(), refG.data(), refG.cols(), SPValueType(0.0),vb_helper.data(), refG.cols());
      for(int i=0; i<refG.cols(); i++) // vb_helper(nCholVecs,nunique_alpha+nunique_beta)
        BLAS::copy(nCholVecs,vb_helper.data()+i,refG.cols(),vb.data()+i*nCholVecs,1);   

      // now excited states: this can be a lot faster if you pack states by virtual orbital and call gemm  
      std::vector<IndexType>::iterator it = iajb_unique_alpha.begin();
      for(int n=0; n<nunique_alpha; n++) {
        it += NAEA;
        int iext = *(it++);
        for(int i=0; i<iext; i++, it+=3) { 
          int ak = *(it+2);
          if(sparse_vn) { 
            int pi0 = *(Spvn->row_index()+ak*NMO);
            SparseMatrixOperators::product_SpMatTV( NMO, Spvn->cols(), SPValueType(1.0), Spvn->values() + pi0, Spvn->column_data() + pi0, Spvn->row_index()+ak*NMO, Gia[n]+(NAEA+i)*NMO , SPValueType(1.0), vb.data()+n*nCholVecs);
          } else
            DenseMatrixOperators::product_Atx( NMO, Dvn->cols(), SPValueType(1.0), Dvn->values() + ak*NMO*Dvn->cols(), Dvn->cols(), Gia[n]+(NAEA+i)*NMO , SPValueType(1.0), vb.data()+n*nCholVecs);
        } 
        it += 2*(*it)+1;   
      } 
      it = iajb_unique_beta.begin();
      for(int n=0; n<nunique_beta; n++) {
        it += NAEB;
        int iext = *(it++);
        for(int i=0; i<iext; i++, it+=3) {
          int ak = *(it+2)-NMO;
          if(sparse_vn) {
            int pi0 = *(Spvn->row_index()+ak*NMO);
            SparseMatrixOperators::product_SpMatTV( NMO, Spvn->cols(), SPValueType(1.0), Spvn->values() + pi0, Spvn->column_data() + pi0, Spvn->row_index()+ak*NMO, Gib[n]+(NAEB+i)*NMO , SPValueType(1.0), vb.data()+(n+nunique_alpha)*nCholVecs);
          } else
            DenseMatrixOperators::product_Atx( NMO, Dvn->cols(), SPValueType(1.0), Dvn->values() + ak*NMO*Dvn->cols(), Dvn->cols(), Gib[n]+(NAEB+i)*NMO , SPValueType(1.0), vb.data()+(n+nunique_alpha)*nCholVecs);
        }
        it += 2*(*it)+1;
      }

/*
      app_log()<<" Checking DM: \n";
      it = iajb_unique_alpha.begin();
      for(int n=0; n<nunique_alpha; n++) {
        int ne = *(it+NAEA);
        int nd = *(it+NAEA+3*ne+2);
app_log()<<"info: " <<ne <<" " <<nd <<std::endl;
        local_evaluateOneBodyMixedDensityMatrix(nd,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,false);
        for(int i=0,ii=0; i<NAEA; i++) {
         int i2 = *(it+i); 
         int i1 = (i2>=NAEA)?(NAEA+(ii++)):i2; 
         for(int k=0; k<NMO; k++) {
           if(std::abs(Gia(n,i1*NMO+k)-mixed_density_matrix(i2,k)) > 1e-8)
            app_log()<<n <<"  " <<i <<" " <<k <<" " <<std::abs(Gia(n,i1*NMO+k)-mixed_density_matrix(i2,k)) <<" " <<Gia(n,i1*NMO+k) <<" " <<mixed_density_matrix(i2,k) <<std::endl;
         }
        }
        it += NAEA;
        it += 3*(*it)+1;
        it += 2*(*it)+1;
      }

      it = iajb_unique_beta.begin();
      for(int n=0; n<nunique_beta; n++) {
        int ne = *(it+NAEB);
        int nd = *(it+NAEB+3*ne+2);
app_log()<<"info: " <<ne <<" " <<nd <<std::endl;
        local_evaluateOneBodyMixedDensityMatrix(nd,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,false);
        for(int i=0,ii=0; i<NAEB; i++) {
         int i2 = *(it+i)-NMO;
         int i1 = (i2>=NAEB)?(NAEB+(ii++)):i2;
         i2+=NMO;
         for(int k=0; k<NMO; k++) {
           if(std::abs(Gib(n,i1*NMO+k)-mixed_density_matrix(i2,k)) > 1e-8)
            app_log()<<n <<"  " <<i <<" " <<k <<" " <<std::abs(Gia(n,i1*NMO+k)-mixed_density_matrix(i2,k)) <<" " <<Gia(n,i1*NMO+k) <<" " <<mixed_density_matrix(i2,k) <<std::endl;
         }
        }
        it += NAEA;
        it += 3*(*it)+1;
        it += 2*(*it)+1;
      }

      app_log()<<" Checking vb: \n"; 
      it = iajb_unique_alpha.begin();
      for(int n=0; n<nunique_alpha; n++) {
        it += NAEA;
        it += 3*(*it)+1;
        int nd = *(it+1);
        it += 2*(*it)+1;
        local_evaluateOneBodyMixedDensityMatrix(nd,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,false);
        DenseMatrixOperators::product_Atx( NMO*NMO, Dvn->cols(), SPValueType(1.0), Dvn->values(), Dvn->cols(), mixed_density_matrix.data() , SPValueType(0.0), vb_helper.data());
        for(int i=0; i<nCholVecs; i++)
          if(std::abs(vb[n*nCholVecs+i]-vb_helper[i]) > 1e-8)
            app_log()<<n <<"  " <<i <<" " <<std::abs(vb[n*nCholVecs+i]-vb_helper[i]) <<" " <<vb[n*nCholVecs+i] <<" " <<vb_helper[i] <<std::endl;
      }
      it = iajb_unique_beta.begin();
      for(int n=0; n<nunique_beta; n++) {
        it += NAEB;
        it += 3*(*it)+1;
        int nd = *(it+1);
        it += 2*(*it)+1;
        local_evaluateOneBodyMixedDensityMatrix(nd,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,false);
        DenseMatrixOperators::product_Atx( NMO*NMO, Dvn->cols(), SPValueType(1.0), Dvn->values(), Dvn->cols(), mixed_density_matrix.data()+NMO*NMO , SPValueType(0.0), vb_helper.data());
        for(int i=0; i<nCholVecs; i++)
          if(std::abs(vb[(n+nunique_alpha)*nCholVecs+i]-vb_helper[i]) > 1e-8)
            app_log()<<n <<"  " <<i <<" " <<std::abs(vb[(n+nunique_alpha)*nCholVecs+i]-vb_helper[i]) <<" " <<vb[(n+nunique_alpha)*nCholVecs+i] <<" " <<vb_helper[i] <<std::endl;
      }
*/

      // generate Eab
      SPComplexType Eab = zero;
      it = iajb_unique_alpha.begin();
      std::vector<ComplexType>::iterator ovlp = ovlp_unique_alpha.begin();
      // to mitigate unfavourable memory access, try to condense this in a SpMatM oparation
      for(int n=0; n<nunique_alpha; n++,ovlp++) {
        it += NAEA;
        it += 3*(*it)+1;
        if(std::abs(*ovlp)<1e-12) {
          it += 2*(*it)+1;
          continue;
        }
        std::fill(V0.begin(),V0.begin()+nCholVecs,zero);
        int nc = *(it++);
        for(int i=0; i<nc; i++, it+=2) {
          SPComplexType fct = myconj(ci_with_psign[*it]) * ovl_ref * (*ovlp) * ovlp_unique_beta[*(it+1)];    
          BLAS::axpy(nCholVecs,fct,vb.data()+nCholVecs*(nunique_alpha+(*(it+1))),1,V0.data(),1);
        } 
        Eab += BLAS::dot(nCholVecs,vb.data()+n*nCholVecs,V0.data())/dt;
      } 

      // calculate Pia/Pib for all unique determinants

      // reference sector first: SMSpHijkl * refG = refP
      // bounds for orbital ak [Hijkl_bounds(ak-NAEA+1,:), Hijkl_bounds(ak-NAEA+2,:) )
      SparseMatrixOperators::product_SpMatM( int(NAEA*NMO), refG.cols(), refP.rows(), SPValueType(1.0), SMSpHijkl[0].values()+ (*Hijkl_bounds[0]), SMSpHijkl[0].column_data()+(*Hijkl_bounds[0]), Hijkl_bounds[0], Hijkl_bounds[1], refG.data() , refG.cols(), SPValueType(0.0), refP.data(), refG.cols());

      for(int i=0; i<nunique_alpha; i++) //  refP(NAEA*NMO,nunique_alpha+nunique_beta)
        BLAS::copy(NMO*NAEA,refP.data()+i,refP.cols(),Pia[i],1);   
      for(int i=0; i<nunique_beta; i++) //  refP(NAEA*NMO,nunique_alpha+nunique_beta)
        BLAS::copy(NMO*NAEB,refP.data()+(nunique_alpha+i),refP.cols(),Pib[i],1);   

      // add additional contributions
      it = iajb_unique_alpha.begin();
      for(int n=0; n<nunique_alpha; n++) {
        it += NAEA;
        int iext = *(it++);
        if(iext==0) {
          it += 3*iext;
          it += 2*(*it)+1;
          continue;
        }

        for(int i=0; i<iext; i++) {
          int ak = *(it+3*i+2);
          BLAS::copy(NMO,Hijkl_bounds[0]+ak*NMO,local_bounds[0]+i*NMO);
          BLAS::copy(NMO,Hijkl_bounds[1]+ak*NMO,local_bounds[1]+i*NMO);
        }
        SparseMatrixOperators::product_SpMatV( int(iext*NMO), int(NAEA*NMO), SPValueType(1.0), SMSpHijkl[0].values()+ local_bounds(0,0), SMSpHijkl[0].column_data()+ local_bounds(0,0), local_bounds[0], local_bounds[1], Gia[n], SPValueType(0.0), Pia[n]+NAEA*NMO);

        for(int i=0; i<iext; i++) {
          
          int ak = *(it+3*i+2);
          for(int k=0; k<NAEA; k++) {
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEA+1]+k*NMO,local_bounds[0]+k*NMO);
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEA+2]+k*NMO,local_bounds[1]+k*NMO);
          }

          for(int k=0; k<iext; k++) {
            int bk = *(it+3*k+2);
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEA+1]+bk*NMO,local_bounds[0]+(NAEA+k)*NMO);
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEA+2]+bk*NMO,local_bounds[1]+(NAEA+k)*NMO);
          }

          SparseMatrixOperators::product_SpMatV( int((NAEA+iext)*NMO), NMO, SPValueType(1.0), SMSpHijkl[0].values()+ local_bounds(0,0), SMSpHijkl[0].column_data()+ local_bounds(0,0), local_bounds[0], local_bounds[1], Gia[n]+(NAEA+i)*NMO, SPValueType(1.0), Pia[n]);
  
        }
        
        it += 3*iext;
        it += 2*(*it)+1;
      }

      it = iajb_unique_beta.begin();
      for(int n=0; n<nunique_beta; n++) {
        it += NAEB;
        int iext = *(it++);
        if(iext==0) {
          it += 3*iext;
          it += 2*(*it)+1;
          continue;
        }

        for(int i=0; i<iext; i++) {
          int ak = *(it+3*i+2)-NMO;
          BLAS::copy(NMO,Hijkl_bounds[0]+ak*NMO,local_bounds[0]+i*NMO);
          BLAS::copy(NMO,Hijkl_bounds[1]+ak*NMO,local_bounds[1]+i*NMO);
        }
        SparseMatrixOperators::product_SpMatV( int(iext*NMO), int(NAEB*NMO), SPValueType(1.0), SMSpHijkl[0].values()+ local_bounds(0,0), SMSpHijkl[0].column_data()+ local_bounds(0,0), local_bounds[0], local_bounds[1], Gib[n], SPValueType(0.0), Pib[n]+NAEB*NMO);

        for(int i=0; i<iext; i++) {

          int ak = *(it+3*i+2)-NMO;
          for(int k=0; k<NAEB; k++) {
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEB+1]+k*NMO,local_bounds[0]+k*NMO);
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEB+2]+k*NMO,local_bounds[1]+k*NMO);
          }

          for(int k=0; k<iext; k++) {
            int bk = *(it+3*k+2)-NMO;
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEB+1]+bk*NMO,local_bounds[0]+(NAEB+k)*NMO);
            BLAS::copy(NMO,Hijkl_bounds[ak-NAEB+2]+bk*NMO,local_bounds[1]+(NAEB+k)*NMO);
          }

          SparseMatrixOperators::product_SpMatV( int((NAEB+iext)*NMO), NMO, SPValueType(1.0), SMSpHijkl[0].values()+ local_bounds(0,0), SMSpHijkl[0].column_data()+ local_bounds(0,0), local_bounds[0], local_bounds[1], Gib[n]+(NAEB+i)*NMO, SPValueType(1.0), Pib[n]);

        }

        it += 3*iext;
        it += 2*(*it)+1;
      }

/*
      if(first_pass) {
        first_pass=false;
        std::map<IndexType,bool> all_orbs;
        all_orbs.clear();
        for(IndexType i=0; i<2*NMO; i++) all_orbs[i]=false;

        for(std::vector<IndexType>::iterator it=occ_orbs.begin(); it!=occ_orbs.end(); it++)
          all_orbs[*it] = true;

        H2.setDims(2*NMO*NMO,2*NMO*NMO);
        H2.setup(head_of_nodes,name+std::string("Hijkl_H2"),TG.getNodeCommLocal());
        if(!sHam->createHamiltonianForPureDeterminant(dm_type,false,all_orbs,all_orbs,thij,H2,cutoff)) {
          APP_ABORT(" Error during testing. \n\n\n"); 
        }

        thij.clear();
        H1.setDims(2*NMO*NMO,2*NMO*NMO);
        H1.setup(head_of_nodes,name+std::string("Hijkl_H1"),TG.getNodeCommLocal());
        if(!sHam->createHamiltonianForPureDeterminant(dm_type,true,all_orbs,all_orbs,thij,H1,cutoff)) {
          APP_ABORT(" Error during testing. \n\n\n"); 
        }

        all_orbs.clear();
        for(IndexType i=0; i<2*NMO; i++) all_orbs[i]=false;

        for(std::vector<IndexType>::iterator it=occ_orbs.begin(); it!=occ_orbs.begin()+NAEA+NAEB; it++)
          all_orbs[*it] = true;

        H0.setDims(2*NMO*NMO,2*NMO*NMO);
        H0.setup(head_of_nodes,name+std::string("Hijkl_H0"),TG.getNodeCommLocal());
        if(!sHam->createHamiltonianForPureDeterminant(dm_type,false,all_orbs,all_orbs,tthij,H0,cutoff)) {
          APP_ABORT(" Error during testing. \n\n\n");
        }
      }
      app_log()<<" Checking Pia: \n";
      it = iajb_unique_alpha.begin();
      ComplexMatrix testP(NMO,NMO);
      for(int n=0; n<nunique_alpha; n++) {
        int ne = *(it+NAEA);
        int nd = *(it+NAEA+3*ne+2);
        local_evaluateOneBodyMixedDensityMatrix(nd,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,false);
        SparseMatrixOperators::product_SpMatV( NMO*NMO, NMO*NMO, SPValueType(1.0), H1.values(), H1.column_data(), H1.row_index(), mixed_density_matrix.data(), SPValueType(0.0), testP.data());
        for(int i=0,ii=0; i<NAEA; i++) {
         int i2 = *(it+i); 
         int i1 = (i2>=NAEA)?(NAEA+(ii++)):i2; 
         for(int k=0; k<NMO; k++) {
           if(std::abs(Pia(n,i1*NMO+k)-testP(i2,k)) > 1e-8)
            app_log()<<n <<" " <<i1 <<" " <<i2 <<"  " <<i <<" " <<k <<" " <<std::abs(Pia(n,i1*NMO+k)-testP(i2,k)) <<" " <<Pia(n,i1*NMO+k) <<" " <<testP(i2,k) <<std::endl;
         }
        }
        it += NAEA;
        it += 3*(*it)+1;
        it += 2*(*it)+1;
      }
*/

      // generate Eaa+Ebb
      SPComplexType Eaa = zero;
      SPComplexType Ebb = zero;
      it = iajb_unique_alpha.begin();
      ovlp = ovlp_unique_alpha.begin();
      for(int n=0; n<nunique_alpha; n++,ovlp++) {
        it += NAEA;
        it += 3*(*it)+1;
        if(std::abs(*ovlp)<1e-12) {
          it += 2*(*it)+1;
          continue;
        }
        SPComplexType ovln=zero;
        int nc = *(it++);
        for(int i=0; i<nc; i++, it+=2) 
          ovln += myconj(ci_with_psign[*it]) * ovl_ref * (*ovlp) * ovlp_unique_beta[*(it+1)];
        Eaa += ovln*BLAS::dot(Pia.cols(),Gia[n],Pia[n]);
      }
      it = iajb_unique_beta.begin();
      ovlp = ovlp_unique_beta.begin();
      for(int n=0; n<nunique_beta; n++,ovlp++) {
        it += NAEB;
        it += 3*(*it)+1;
        if(std::abs(*ovlp)<1e-12) {
          it += 2*(*it)+1;
          continue;
        }
        SPComplexType ovln=zero;
        int nc = *(it++);
        for(int i=0; i<nc; i++, it+=2)
          ovln += myconj(ci_with_psign[*it]) * ovl_ref * (*ovlp) * ovlp_unique_alpha[*(it+1)];
        Ebb += ovln*BLAS::dot(Pib.cols(),Gib[n],Pib[n]);
      }
      Eaa*=0.5;   
      Ebb*=0.5;   
      pe_nume = (Eab+Eaa+Ebb);

/*
      // testing
      ComplexType ke_nume_t = ke_nume;
      ComplexType pe_nume_t = pe_nume;
      ComplexType deno_ = deno; 
      ke_nume=pe_nume=deno=zero;


      ComplexType eaa_nume = ComplexType(0);
      ComplexType ebb_nume = ComplexType(0);
      for(int ki=0; ki<ci.size(); ki++)
      {
        lke=lpe=ComplexType(0.0,0.0);

        local_evaluateOneBodyMixedDensityMatrix(ki,SlaterMat,overlaps[2*ki],overlaps[2*ki+1],mixed_density_matrix,false);
        ComplexType scl = myconj(ci[ki])*overlaps[2*ki]*overlaps[2*ki+1];
        deno += scl;

        SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
        s1Dit end1 = thij.end();
        for(s1Dit it = thij.begin(); it != end1; it++)
          lke += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
        ke_nume += lke*scl;

        SparseMatrixOperators::product_SpMatV(nc1,nc1,SPValueType(1),H2.values(),H2.column_data(),H2.row_index(),mixed_density_matrix.data(),SPValueType(0),V0.data());
        itG = mixed_density_matrix.begin();
        SPComplexVector::iterator itV = V0.begin();
        for(int i=0; i<nc1; i++,++itG,++itV) lpe += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
        pe_nume += (0.5*lpe)*scl;

        lpe=zero;
        SparseMatrixOperators::product_SpMatV(NMO*NMO,NMO*NMO,SPValueType(1),H1.values(),H1.column_data(),H1.row_index(),mixed_density_matrix.data(),SPValueType(0),V0.data());
        itG = mixed_density_matrix.begin();
        itV = V0.begin();
        for(int i=0; i<NMO*NMO; i++,++itG,++itV) lpe += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
        eaa_nume += (0.5*lpe)*scl;

        lpe=zero;
        SparseMatrixOperators::product_SpMatV(NMO*NMO,NMO*NMO,SPValueType(1),H1.values(),H1.column_data(),H1.row_index(),mixed_density_matrix.data()+NMO*NMO,SPValueType(0),V0.data());
        itG = mixed_density_matrix.begin()+NMO*NMO;
        itV = V0.begin();
        for(int i=0; i<NMO*NMO; i++,++itG,++itV) lpe += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
        ebb_nume += (0.5*lpe)*scl;

//        app_log()<<ki <<" " <<lke <<" "  <<0.5*lpe <<" " <<pe[ki] <<" " <<0.5*lpe <<"\n";

      }

      if(std::abs(ke_nume_t/deno_-ke_nume/deno) > 1e-5) 
        app_log()<<" Differences in ke: " <<std::abs(ke_nume_t/deno_-ke_nume/deno) <<" " <<ke_nume_t/deno_ <<" " <<ke_nume/deno <<"\n";
      if(std::abs(pe_nume_t/deno_-pe_nume/deno) > 1e-5) 
        app_log()<<" Differences in pe: " <<std::abs(pe_nume_t/deno_-pe_nume/deno) <<" " <<pe_nume_t/deno_ <<" " <<pe_nume/deno <<"\n";

      app_log()<<"\n\n\n" <<" Overlap: " <<deno <<" " <<deno_ <<"\n";
      app_log()<<" Total kinetic energy: " <<ke_nume/deno <<" " <<ke_nume_t/deno_ <<"\n";
      app_log()<<" Total potential energy: " <<pe_nume/deno <<" " <<pe_nume_t/deno_ <<"\n";
      app_log()<<" Eaa: " <<Eaa/deno_ <<" Ebb:" <<Ebb/deno_ <<"  Eab:" <<Eab/deno_ <<std::endl; 
      app_log()<<" Eaa: " <<eaa_nume/deno <<" Ebb:" <<ebb_nume/deno <<"  Eab:" <<(pe_nume-eaa_nume-ebb_nume)/deno <<std::endl; 

      local_evaluateOneBodyMixedDensityMatrix(0,SlaterMat,overlaps[0],overlaps[1],mixed_density_matrix,false);

      Timer.reset("Generic1");
      Timer.start("Generic1");

      SparseMatrixOperators::product_SpMatV(nc1,nc1,SPValueType(1),H0.values(),H0.column_data(),H0.row_index(),mixed_density_matrix.data(),SPValueType(0),V0.data());

      Timer.stop("Generic1");
      app_log()<<" Time for 1 SpMatV: " <<Timer.total("Generic1") <<"\n\n\n\n";


      //APP_ABORT(" TESTING TESTING TESTING \n\n\n"); 
*/

    } else {

      for(int ki=0; ki<ci.size(); ki++)
      {
        lke=lpe=ComplexType(0.0,0.0);
 
        local_evaluateOneBodyMixedDensityMatrix(ki,SlaterMat,overlaps[2*ki],overlaps[2*ki+1],mixed_density_matrix,false);
        ComplexType scl = myconj(ci[ki])*overlaps[2*ki]*overlaps[2*ki+1];
        deno += scl; 

        SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
        if(rotated_hamiltonian) {
          std::vector<s1D<ComplexType> >::iterator  end1 = haj[ki].end();
          for(std::vector<s1D<ComplexType> >::iterator it = haj[ki].begin(); it != end1; it++)
            lke += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
        } else {
          s1Dit end1 = hij[ki].end();
          for(s1Dit it = hij[ki].begin(); it != end1; it++)
            lke += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
        }
        ke[ki]=lke;
        ke_nume += lke*scl;

        if(rotated_hamiltonian) {
          SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHabkl[ki].values(),SMSpHabkl[ki].column_data(),SMSpHabkl[ki].row_index(),mixed_density_matrix.data(),zero,V0.data());
        } else {
          SparseMatrixOperators::product_SpMatV(nc1,nc1,SPValueType(1),SMSpHijkl[ki].values(),SMSpHijkl[ki].column_data(),SMSpHijkl[ki].row_index(),mixed_density_matrix.data(),SPValueType(0),V0.data());
        }
        itG = mixed_density_matrix.begin();
        SPComplexVector::iterator itV = V0.begin();
        for(int i=0; i<nc1; i++,++itG,++itV) lpe += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
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
      if( !std::isfinite( std::abs(ovl_alpha*ovl_beta) ) ) 
        ovl_alpha=ovl_beta=1e-12;
      epot = 0.0;
      ekin = 0.0;
    }
#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:evaluateLocalEnergy");
#endif
   
  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMSpMat& vn, SPComplexSMSpMat& vnT, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n)
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

    if(transposed) {
      SPComplexType one = SPComplexType(1.0);
      const SPComplexType zero = SPComplexType(0.0);
      if(closed_shell) one = SPComplexType(2.0);
      SparseMatrixOperators::product_SpMatV(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.column_data(),vnT.row_index(),full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta && ! closed_shell)
        SparseMatrixOperators::product_SpMatV(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.column_data(),vnT.row_index(),full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    } else {
      SPValueType one = SPValueType(1.0);
      SPValueType zero = SPValueType(0.0);
      if(closed_shell) one = SPValueType(2.0);
      SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:product_Spvn_x_DMfull");
    Timer.stop("MultiPureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMVector& vn, SPComplexSMVector& vnT, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n)
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

    if(transposed) {
      SPComplexType one = SPComplexType(1.0);
      const SPComplexType zero = SPComplexType(0.0);
      if(closed_shell) one = SPComplexType(2.0);
      DenseMatrixOperators::product_Ax(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.cols(),full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta && !closed_shell)
        DenseMatrixOperators::product_Ax(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.cols(),full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    } else {
      SPValueType one = SPValueType(1.0);
      SPValueType zero = SPValueType(0.0);
      if(closed_shell) one = SPValueType(2.0);
      DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one,vn.values(),vn.cols(),full_mixed_density_matrix.data(),zero,v.data());
      if(addBetaBeta && !closed_shell)
        DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one,vn.values(),vn.cols(),full_mixed_density_matrix.data()+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("MultiPureSingleDeterminant:product_Spvn_x_DMfull");
    Timer.stop("MultiPureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
#endif

  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int i0, int iN, SPValueSMSpMat& vn, SPComplexSMSpMat& vnT, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n)
  {
    if(transposed)
      assert(v.size() >= nW*(iN-i0));
    else
      assert(v.size() >= nW*vn.cols());
    const SPComplexType *GF = buff+1*nW;
    int GFi0 = i0;
    int GFi02 = NMO*NMO+i0;
    SPValueType one = SPValueType(1.0);
    const SPValueType zero = SPValueType(0.0);
    if(closed_shell) one = SPValueType(2.0);
    const SPComplexType czero = SPComplexType(0.0,0.0);
    SPComplexType cone = SPComplexType(1.0);
    if(closed_shell) cone = SPComplexType(2.0);
    std::size_t p_;
    if(transposed)
      p_ = static_cast<std::size_t>(*(vnT.row_index()+i0));
    else
      p_ = static_cast<std::size_t>(*(vn.row_index()+i0));
    // walkerBlock implies MatV
    if(walkerBlock==1) {

      if(transposed) {
        SparseMatrixOperators::product_SpMatV( int(iN-i0), vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF, czero, v.data());
        if(addBetaBeta && !closed_shell)
          SparseMatrixOperators::product_SpMatV( int(iN-i0), vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+NMO*NMO, cone, v.data());
      } else {
        SparseMatrixOperators::product_SpMatTV( int(iN-i0), vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+i0, zero, v.data());
        if(addBetaBeta && !closed_shell)
          SparseMatrixOperators::product_SpMatTV( int(iN-i0), vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+NMO*NMO+i0, one, v.data());
      }

    } else {    
    // walkerBlock>1, so use MatM  
    // almost certainly, it is highly favourable  to set walkerBlock to nW
      int nblk = nW/walkerBlock;
      int nextra = nW%walkerBlock;

      for(int i=0; i<nblk; i++) {

        if(transposed) {
          SparseMatrixOperators::product_SpMatM( int(iN-i0), walkerBlock , vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+i*walkerBlock, nW, czero, v.data()+i*walkerBlock, nW);
          if(addBetaBeta && !closed_shell)
            SparseMatrixOperators::product_SpMatM( int(iN-i0), walkerBlock , vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+(NMO*NMO)*nW+i*walkerBlock, nW, cone, v.data()+i*walkerBlock, nW);
        } else {
          SparseMatrixOperators::product_SpMatTM( int(iN-i0), walkerBlock , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+i0*nW+i*walkerBlock, nW, zero, v.data()+i*walkerBlock, nW);
          if(addBetaBeta && !closed_shell)
            SparseMatrixOperators::product_SpMatTM( int(iN-i0), walkerBlock , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+(i0+NMO*NMO)*nW+i*walkerBlock, nW, one, v.data()+i*walkerBlock, nW);
        }

      }
      // perform remaining walkers. Should I pad this to the closest power of 2???
      if(nextra > 0) {
        if(transposed) {
          SparseMatrixOperators::product_SpMatM( int(iN-i0), nextra , vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+nblk*walkerBlock, nW, czero, v.data()+nblk*walkerBlock, nW);
          if(addBetaBeta && !closed_shell)
            SparseMatrixOperators::product_SpMatM( int(iN-i0), nextra , vnT.cols(), cone, vnT.values() + p_, vnT.column_data() + p_, vnT.row_index()+i0, GF+NMO*NMO*nW+nblk*walkerBlock, nW, cone, v.data()+nblk*walkerBlock, nW);
        } else {
          SparseMatrixOperators::product_SpMatTM( int(iN-i0), nextra , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+i0*nW+nblk*walkerBlock, nW, zero, v.data()+nblk*walkerBlock, nW);
          if(addBetaBeta && !closed_shell)
            SparseMatrixOperators::product_SpMatTM( int(iN-i0), nextra , vn.cols(), one, vn.values() + p_, vn.column_data() + p_, vn.row_index()+i0, GF+(i0+NMO*NMO)*nW+nblk*walkerBlock, nW, one, v.data()+nblk*walkerBlock, nW);
        }
      }
    }
  }      

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int i0, int iN, SPValueSMVector& vn, SPComplexSMVector& vnT, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n)
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

    // walkerBlock implies MatV
    if(walkerBlock==1) {
      if(transposed) {
        DenseMatrixOperators::product_Ax( int(iN-i0), vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(), GF, czero, v.data());
        if(addBetaBeta && !closed_shell)
          DenseMatrixOperators::product_Ax( int(iN-i0), vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(), GF+NMO*NMO, cone, v.data());
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
          DenseMatrixOperators::product( int(iN-i0), walkerBlock, vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(),
               GF+i*walkerBlock,nW,czero,v.data()+i*walkerBlock,nW);
          if(addBetaBeta && !closed_shell)
            DenseMatrixOperators::product( int(iN-i0), walkerBlock, vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(),
               GF+NMO*NMO*nW+i*walkerBlock,nW,cone,v.data()+i*walkerBlock,nW);
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
          DenseMatrixOperators::product(int(iN-i0), nextra, vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(),
               GF+nblk*walkerBlock,nW,czero,v.data()+nblk*walkerBlock,nW);
          if(addBetaBeta && !closed_shell)
            DenseMatrixOperators::product(int(iN-i0), nextra, vnT.cols(), cone, vnT.values() + i0*vnT.cols(), vnT.cols(),
               GF+NMO*NMO*nW+nblk*walkerBlock,nW,cone,v.data()+nblk*walkerBlock,nW);
        } else {
          DenseMatrixOperators::product_AtB(vn.cols(), nextra, int(iN-i0), one, vn.values() + i0*vn.cols(), vn.cols(),
               GF+GFi0*nW+nblk*walkerBlock,nW,zero,v.data()+nblk*walkerBlock,nW);
          if(addBetaBeta && !closed_shell)
            DenseMatrixOperators::product_AtB(vn.cols(), nextra, int(iN-i0), one, vn.values() + i0*vn.cols(), vn.cols(),
               GF+GFi02*nW+nblk*walkerBlock,nW,one,v.data()+nblk*walkerBlock,nW);
        }
      }
    }

  }

  void MultiPureSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<SPComplexType> >& v2n, const std::vector<IndexType>& vn_indx, SPValueSMSpMat& vn, std::vector<SPComplexType>& v, const int n)
  {
  }

void MultiPureSingleDeterminant::evaluateTrialEnergy(ComplexType& ke, ComplexType& pe)
{
  SPComplexType one = SPComplexType(1.0,0.0);
  SPComplexType zero = SPComplexType(0.0,0.0);
  
  ke=pe=zero;
  if(myComm->rank() == 0) {  // serial for now

   if(rotated_hamiltonian) {  // since I do not have a way to do this right now  

    int nr1=1, nc1=2*NMO*NMO;
    ComplexType oa,ob,pe_nume=zero,deno=zero;
    ComplexType lke,lpe,ke_nume=zero; 
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfn_type==0)?0:NAEA*NMO; 

    for(int i=0; i<ci.size(); i++) {

      lke=lpe=ComplexType(0.0,0.0);
      
      for(int ii=0; ii<NMO; ii++)
       for(int j=0; j<NAEA; j++)
        Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
      for(int ii=0; ii<NMO; ii++)
       for(int j=0; j<NAEB; j++)
        Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];
      local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,SPrank_updated_trial_density_matrix);
      
      SPComplexMatrix::iterator itG = SPrank_updated_trial_density_matrix.begin();
      std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
        lke += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
      
      SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHabkl[i].values(),SMSpHabkl[i].column_data(),SMSpHabkl[i].row_index(),SPrank_updated_trial_density_matrix.data(),zero,V0.data());
      SPComplexVector::iterator itV = V0.begin(); 
      for(int ii=0; ii<nc1; ii++,++itG,++itV) lpe += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
      
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
        local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,SPrank_updated_trial_density_matrix,false);

        SPComplexMatrix::iterator itG = SPrank_updated_trial_density_matrix.begin();
        std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
        for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
          lke += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);

        SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHabkl[i].values(),SMSpHabkl[i].column_data(),SMSpHabkl[i].row_index(),SPrank_updated_trial_density_matrix.data(),zero,V0.data());
        SPComplexVector::iterator itV = V0.begin();
        for(int ii=0; ii<nc1; ii++,++itG,++itV) lpe += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);

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

void MultiPureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMVector& vn, std::vector<SPComplexType>& v, const int nstep)
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
    SPtrial_density_matrix = ComplexType(0.0);
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfn_type==0)?0:NAEA*NMO;

    for(int i=0; i<ci.size(); i++) {
     // i==j
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEA; j++)
       Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEB; j++)
       Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];

     local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, SPrank_updated_trial_density_matrix,true);
     o1 = myconj(ci[i])*ci[i]*oa*ob;
     SPtrial_density_matrix += static_cast<SPComplexType>(o1)*SPrank_updated_trial_density_matrix; 
     ovl += o1 ;
     for(int j=i+1; j<ci.size(); j++)
     {
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEA; jj++)
         Am(ii,jj) = OrbMat[j*orbsize+ii*NAEA+jj];
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEB; jj++)
         Am(ii+NMO,jj) = OrbMat[j*orbsize+ii*NAEA+jj+del];
       local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, SPrank_updated_trial_density_matrix,true);
       o1 = myconj(ci[i])*ci[j]*oa*ob; 
       SPtrial_density_matrix += static_cast<SPComplexType>(o1)*SPrank_updated_trial_density_matrix; 
       ovl += o1 ;
       o1 = ci[i]*myconj(ci[j])*myconj(oa*ob);
       // is there an operation for std::complex conjugation???
       for(int ii=0; ii<NMO; ii++)  
        for(int jj=0; jj<NMO; jj++) { 
         SPtrial_density_matrix(ii,jj) += static_cast<SPComplexType>(o1)*myconj(SPrank_updated_trial_density_matrix(jj,ii));
         SPtrial_density_matrix(ii+NMO,jj) += static_cast<SPComplexType>(o1)*myconj(SPrank_updated_trial_density_matrix(jj+NMO,ii));
        }
       ovl += o1 ;
     }
    }
    SPtrial_density_matrix *= (SPComplexType(1.0,0.0)/static_cast<SPComplexType>(ovl));

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
    SPtrial_density_matrix = trial_density_matrix;
  }

  SPValueType one_=SPValueType(1);
  const SPValueType zero_=SPValueType(0);
  if(closed_shell) one_=SPValueType(2.0);
  DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one_,vn.values(),vn.cols(),SPtrial_density_matrix.data(),zero_,v.data());
  if(addBetaBeta && !closed_shell)
    DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one_,vn.values(),vn.cols(),SPtrial_density_matrix.data()+NMO*NMO,one_,v.data());

}

void MultiPureSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMSpMat& vn, std::vector<SPComplexType>& v, const int nstep)
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
    SPtrial_density_matrix = ComplexType(0.0);
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfn_type==0)?0:NAEA*NMO;

    for(int i=0; i<ci.size(); i++) {
     // i==j
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEA; j++)
       Am(ii,j) = OrbMat[i*orbsize+ii*NAEA+j];
     for(int ii=0; ii<NMO; ii++)
      for(int j=0; j<NAEB; j++)
       Am(ii+NMO,j) = OrbMat[i*orbsize+ii*NAEA+j+del];

     local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, SPrank_updated_trial_density_matrix,true);
     o1 = myconj(ci[i])*ci[i]*oa*ob;
     SPtrial_density_matrix += static_cast<SPComplexType>(o1)*SPrank_updated_trial_density_matrix; 
     ovl += o1 ;
     for(int j=i+1; j<ci.size(); j++)
     {
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEA; jj++)
         Am(ii,jj) = OrbMat[j*orbsize+ii*NAEA+jj];
       for(int ii=0; ii<NMO; ii++)
        for(int jj=0; jj<NAEB; jj++)
         Am(ii+NMO,jj) = OrbMat[j*orbsize+ii*NAEA+jj+del];
       local_evaluateOneBodyMixedDensityMatrix(i, Am.data(), oa, ob, SPrank_updated_trial_density_matrix,true);
       o1 = myconj(ci[i])*ci[j]*oa*ob; 
       SPtrial_density_matrix += static_cast<SPComplexType>(o1)*SPrank_updated_trial_density_matrix; 
       ovl += o1 ;
       o1 = ci[i]*myconj(ci[j])*myconj(oa*ob);
       // is there an operation for std::complex conjugation???
       for(int ii=0; ii<NMO; ii++)  
        for(int jj=0; jj<NMO; jj++) { 
         SPtrial_density_matrix(ii,jj) += static_cast<SPComplexType>(o1)*myconj(SPrank_updated_trial_density_matrix(jj,ii));
         SPtrial_density_matrix(ii+NMO,jj) += static_cast<SPComplexType>(o1)*myconj(SPrank_updated_trial_density_matrix(jj+NMO,ii));
        }
       ovl += o1 ;
     }
    }
    SPtrial_density_matrix *= (SPComplexType(1.0,0.0)/static_cast<SPComplexType>(ovl));

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
    SPtrial_density_matrix = trial_density_matrix;
  }

  SPValueType one_=SPValueType(1);
  const SPValueType zero_=SPValueType(0);
  if(closed_shell) one_=SPValueType(2.0);
  SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one_,vn.values(),vn.column_data(),vn.row_index(),SPtrial_density_matrix.data(),zero_,v.data());
  if(addBetaBeta && !closed_shell)
    SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one_,vn.values(),vn.column_data(),vn.row_index(),SPtrial_density_matrix.data()+NMO*NMO,one_,v.data());

}

bool MultiPureSingleDeterminant::diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ComplexMatrix& eigVec, std::vector<IndexType>& occv, int nci, bool eigV )
{
  SPComplexType one = SPComplexType(1.0,0.0);
  SPComplexType zero = SPComplexType(0.0,0.0);
  bool sucess;
  //int nci = ci.size(); 
  
  if(rotated_hamiltonian) { 

   // serial for now
   if(myComm->rank() == 0) { 
  
    int nr1=1, nc1=2*NMO*NMO;
    ComplexType oa,ob,pe_nume=zero,deno=zero;
    ComplexType lke,lpe,ke_nume=zero; 
    ComplexMatrix Am(2*NMO,NAEA);
    int del = (wfn_type==0)?0:NAEA*NMO; 

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
      local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,SPrank_updated_trial_density_matrix);
      
      SPComplexMatrix::iterator itG = SPrank_updated_trial_density_matrix.begin();
      SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHabkl[i].values(),SMSpHabkl[i].column_data(),SMSpHabkl[i].row_index(),SPrank_updated_trial_density_matrix.data(),zero,V0.data());
      SPComplexVector::iterator itV = V0.begin(); 
      for(int ii=0; ii<nc1; ii++,++itG,++itV) let += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
      let *= 0.5;

      itG = SPrank_updated_trial_density_matrix.begin();
      std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
      for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
        let += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);
      
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
        local_evaluateOneBodyMixedDensityMatrix(i,Am.data(),oa,ob,SPrank_updated_trial_density_matrix,false);

        SPComplexMatrix::iterator itG = SPrank_updated_trial_density_matrix.begin();
        SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHabkl[i].values(),SMSpHabkl[i].column_data(),SMSpHabkl[i].row_index(),SPrank_updated_trial_density_matrix.data(),zero,V0.data());
        SPComplexVector::iterator itV = V0.begin();
        for(int ii=0; ii<nc1; ii++,++itG,++itV) let += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
        let *= 0.5;

        itG = SPrank_updated_trial_density_matrix.begin();
        std::vector<s1D<ComplexType> >::iterator  end1 = haj[i].end();
        for(std::vector<s1D<ComplexType> >::iterator it = haj[i].begin(); it != end1; it++)
          let += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);

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
      app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic3") <<std::endl;

      Timer.reset("Generic4");
      Timer.start("Generic4");

      eigVal.resize(1);
      eigVec.resize(1,nci);
      std::vector<int> ifail(nci);
      sucess = DenseMatrixOperators::genHermitianEigenSysSelect(nci,hm.data(),nci,ov.data(),nci,1,eigVal.data(),eigV,eigVec.data(),eigVec.size2(),ifail.data());

      Timer.stop("Generic4");
      app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic4") <<std::endl;

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

// not very efficient, but simple
int MultiPureSingleDeterminant::countExct(int NE, IndexType* ref, IndexType* D, bool getIndx, IndexType* loc, IndexType* ik, IndexType* ak, RealType& psign) {
  int cnt=0;
  if(getIndx) std::copy(ref,ref+NE,Iwork.data()); 
  IndexType* it = ref;
  for(int i=0; i<NE; i++,it++)
    if(!std::binary_search(D,D+NE,*it)) {
      if(getIndx) { 
        ik[cnt]=*it;
        loc[cnt]=i; 
      }
      cnt++;
    }
  if(!getIndx)
    return cnt;
  it = D;
  int cnt2=0;
  for(int i=0; i<NE; i++,it++)
    if(!std::binary_search(ref,ref+NE,*it)) {
      if(getIndx) {
        ak[cnt2]=*it;
        Iwork[loc[cnt2]]=*it; 
      }
      cnt2++;
    }
  assert(cnt==cnt2);
  psign=RealType(1.0);
  // sort Iwork and count number of exchanges to determine permutation sign
  // sooo slow but sooo simple too
  for(int i=0; i<NE; i++)
    for(int j=i+1; j<NE; j++)
    {
      if(Iwork[j] < Iwork[i])
      {
        psign*=RealType(-1.0);
        std::swap(Iwork[i],Iwork[j]); 
      }
    }  
  return cnt;
}


// careful 
// FIX FIX FIX
// this routine returns interchanged (i,j)/(k,l), so it is wrong due to a std::complex conjugation for std::complex matrix elements 
int MultiPureSingleDeterminant::cntExcitations(std::vector<IndexType>& DL, std::vector<IndexType>& DR, IndexType& n0, IndexType& n1, IndexType& n2, IndexType& n3, std::vector<IndexType>& occ, RealType& sg)
{
  std::vector<IndexType>::iterator itR = DR.begin();
  std::vector<IndexType>::iterator itL = DL.begin();
  sg = 0.0;
  int cnt=0,pos=0,ind[1000],cnt2=0,nq=0,cnt3=0;
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

      
