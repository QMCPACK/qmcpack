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
#include "AFQMC/Wavefunctions/GeneralSingleDeterminant.h"

#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

//#include "mkl.h"
//#include "mkl_service.h"

namespace qmcplusplus
{

bool GeneralSingleDeterminant::parse(xmlNodePtr cur)
{
    if(cur == NULL)
      return false;

    app_log()<<"\n\n --------------- Parsing GenSD input ------------------ \n\n";

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

    if(type != "gensd") {
      app_error()<<" ERROR: Problems in GeneralSingleDeterminant::parse: type should be gensd. \n"  <<std::endl; 

      return false;
    }

    filename = std::string("none");
    filetype = std::string("ascii");
    ParameterSet m_param;
    m_param.add(filename,"filename","std::string");
    m_param.add(filetype,"filetype","std::string");
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

bool GeneralSingleDeterminant::initFromAscii(std::string fileName)
{

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
      } else if(*it == "FullMO" || *it == "FULLMO") {
        fullMOMat = true;
        app_log()<<" Expecting full MO matrix in GeneralSlaterDeterminant.\n";
      } else if(*it == "CMajor") {
        Cstyle = false; 
        app_log()<<" Expecting MO matrix in Column-major format in GeneralSlaterDeterminant.\n";
      } else if(*it == "UHF" || *it == "GHF") {
        if( it+1 == words.end() ) {
          app_error()<<"Format error in ASCII integral file. UHF/GHF \n";
          return false;
        }
        wfn_type = atoi((++it)->c_str());
        switch(wfn_type) {
          case 0:
          {
            app_log()<<"Using a RHF-type trial wave-function in GeneralSlaterDeterminant. \n";  
            break;
          }
          case 1:
          {
            app_log()<<"Using a UHF-type trial wave-function in GeneralSlaterDeterminant. \n";  
            break;
          }
          case 2:
          {
            app_log()<<"Using a GHF-type trial wave-function in GeneralSlaterDeterminant. \n";  
            break;
          }
          default: 
          {
            app_error()<<"Unknown wave-function type in GeneralSlaterDeterminant: " <<wfn_type <<std::endl;
            return false; 
          } 
        }
      } else {
      }
    }
    getwords(words,in);
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
  } while((words[0].find(std::string("/"))==std::string::npos && words[0].find(std::string("&END"))==std::string::npos));

  int ncols = NAEA;
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
           app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
           app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
           app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
           app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
           app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
           app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
         in.close();
         return false;
       }
     }
     for(int j=0; j<nread2; j++) {
       in>>dummy;
       if(j<NAEB) OrbMat(i,j+NAEA) = dummy;
       if(in.fail()) {
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
         in.close();
         return false;
       }
     }
     for(int j=0; j<nread2; j++) {
       in>>dummy; 
       if(j<NAEB) OrbMat(i+NMO,j+NAEA) = dummy;
       if(in.fail()) {
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
         in.close();
         return false;
       }
     }
    }

   } else { 

    // mo matrix in column ordering 
/*    
    nread = fullMOMat?NMO:NAEA;  
    for(int j=0; j<nread; j++) {
      for(int i=0; i<2*NMO; i++) {
       in>>dummy; 
       if(j<NAEA) OrbMat(i,j) = dummy;
       if(in.fail()) {
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
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
         app_error()<<"Problems reading ASCII file in GeneralSingleDeterminant. \n";
         in.close();
         return false;
       }
     }
    }

   }

  }

  in.close();

  return setup_local();

}

// this should be moved to further down in the initialization line, since at this point I
// don't think I know walker_type.
bool GeneralSingleDeterminant::setup_local()
{

  NCA=NCB=0;
  if(NMO <= 0 || NAEA <= 0 || NAEB <= 0 ) {
    app_error()<<"Problems with initial values of NMO,NAEA,NAEB in GeneralSingleDeterminant.\n";
    return false;
  } 

  int ncols = NAEA;
  if(wfn_type == 2)
    ncols = NAEA+NAEB;  

  // NOTE: this is wrong if walker_type > 1
  int dm_ncols = NMO;
  if(wfn_type == 2)
    dm_ncols *= 2; 

  if(init_type == "diagh1" && !readHamFromFile) {
    OrbMat.resize(2*NMO,ncols);
    // initialize later
    // this is mainly for model hamiltonians and projected initial states
  } else if(init_type == "ground" && !readHamFromFile) {
    OrbMat.resize(2*NMO,ncols);
    for(int i=0; i<NAEA; i++) OrbMat(i,i)=ComplexType(1.0); 
    int shft = (wfn_type == 2)?NAEA:0; 
    for(int i=0; i<NAEB; i++) OrbMat(i+NMO,i+shft)=ComplexType(1.0); 
  } 

  //trial_density_matrix.resize(NAEA+NAEB,NMO);
  trial_density_matrix.resize(2*NMO,dm_ncols);
  SPtrial_density_matrix.resize(2*NMO,dm_ncols);

  //mixed_density_matrix.resize(NAEA+NAEB,NMO);
  mixed_density_matrix.resize(2*NMO,dm_ncols);
  temp_density_matrix.resize(2*NMO,dm_ncols);

  // not used right now 
  //overlap_inv.resize(NAEA+NAEB,NAEA); 

  // temporary storage
  S0.resize(NAEA,NAEA); 
  if(wfn_type == 2) {
    T0.resize(NAEA+NAEB,NAEA+NAEB); 
    T1.resize(NAEA+NAEB,2*NMO); 
    SM.resize(2*NMO,NAEA+NAEB); 
    SS0.resize(NAEA,NMO); 
  } else {
    SS0.resize(2*NMO,NAEA); 
  }
  S1.resize(NAEB,NAEB); 
  V0.resize(2*NMO*dm_ncols); 

  Cwork.resize(2*NMO);
  pivot.resize(2*NMO);

  //SpHij.setDims(1,2*NMO*NMO);
  if(!readHamFromFile) SMSpHijkl.setDims(2*NMO*dm_ncols,2*NMO*dm_ncols);
  if(!readHamFromFile) SMSpHijkl.setup(head_of_nodes,"SMSpHijkl",TG.getNodeCommLocal());

  return true;
}

bool GeneralSingleDeterminant::setup(HamPtr cur) {
  return getHamiltonian(cur);
}

bool GeneralSingleDeterminant::getHamiltonian(HamPtr h)
{
  ham0 = h;
  sHam = dynamic_cast<SparseGeneralHamiltonian*>(h);
  if(!sHam) { 
    app_error()<<" Error in GeneralSingleDeterminant::getHamiltonian. \n"
               <<" Hamiltonian associated with GeneralSingleDeterminant must of the \n"
               <<" type SparseGeneralHamiltonian. \n";
    return false;
  }

  spinRestricted = sHam->RHF();
  closed_shell = false;
  if(wfn_type==0 && NCA==NCB && NAEA == NAEB && sHam->RHF() && init_type != "diagh1") {
    closed_shell = std::equal(OrbMat.begin(),OrbMat.begin()+NMO*NAEA,OrbMat.begin()+NMO*NAEA);
  }

  // walker_type will come from WalkerHandler, but for now only ROHF/UHF is implemented correctly.
  walker_type=1;
  dm_type = std::max(wfn_type,1);
  if(closed_shell) {
    app_log()<<"Found closed shell system. " <<std::endl;
    dm_type=0;
  }

  // only allow GHF walkers with GHF wavefunctions. The opposite doesn't make sense to me.
  if(walker_type == 2 && wfn_type < 2)
    APP_ABORT(" Error in GeneralSingleDeterminant::getHamiltonian: GHF walker is only allowed with GHF wavefunction. \n\n\n");

  if(useFacHam) {
    APP_ABORT(" Error: Use of factorized hamiltonian is not implemented in PureSD. \n\n\n");
  }

  NuclearCoulombEnergy = static_cast<ValueType>(sHam->NuclearCoulombEnergy);
  if(!readHamFromFile) {
    std::map<IndexType,bool> isOcc; 
    isOcc.clear();
    for(IndexType i=0; i<2*NMO; i++) isOcc[i]=true;
    if(!sHam->createHamiltonianForPureDeterminant(dm_type,useFacHam,isOcc,isOcc,hij,SMSpHijkl,cutoff)) {
      app_error()<<"Error in createHamiltonianForGeneralDeterminant. \n";
      return false;
    }
  }

  if(init_type == "diagh1" && !readHamFromFile) {
    if(wfn_type==2) {
      app_error()<<" Error in GeneralSlaterDeterminant::getHamiltonian: init_type=diagh1 not implemented for GHF wavefunctions. \n" <<std::endl;
      return false;
    }
    ValueMatrix tH1(NMO), eigVec(NMO);
    RealVector eigVal(NMO);
    s1Dit it; 
    for(it=hij.begin(); it!=hij.end(); it++) {
      if(std::get<0>(*it) >= NMO*NMO) break; 
      int i = std::get<0>(*it)/NMO; 
      int j = std::get<0>(*it)%NMO;
      tH1(i,j) = std::get<1>(*it);
    }
    if(!DenseMatrixOperators::symEigenSysAll(NMO,tH1.data(),NMO,eigVal.data(),eigVec.data(),NMO) ) {
      app_error()<<"Problems with eigenvalue/eigenvector calculation during initialization of wavefunction with init_type=diagH1.\n";
      return false; 
    }    
    std::vector<std::tuple<RealType,int> > v(NMO);
    for(int i=0; i<NMO; i++) v[i] = std::forward_as_tuple(eigVal[i],i); 
    std::sort (v.begin(), v.end(), [] (const std::tuple<RealType,int>& a, const std::tuple<RealType,int>& b) {return std::get<0>(a) < std::get<0>(b);});
    for(int i=0; i<NMO; i++)
     for(int j=0; j<NAEA; j++)
      OrbMat(i,j) = eigVec(i, std::get<1>(v[j]) );   

    if(!spinRestricted) {
      tH1=0;
      eigVal=0;
      eigVec=0; 
      for(; it!=hij.end(); it++) {
        int i = (std::get<0>(*it)-NMO*NMO)/NMO;
        int j = (std::get<0>(*it)-NMO*NMO)%NMO;
        tH1(i,j) = std::get<1>(*it);
      }
      if(!DenseMatrixOperators::symEigenSysAll(NMO,tH1.data(),NMO,eigVal.data(),eigVec.data(),NMO) ) {
        app_error()<<"Problems with eigenvalue/eigenvector calculation during initialization of wavefunction with init_type=diagH1.\n";
      return false;
      }
      std::vector<std::tuple<RealType,int> > v(NMO);
      for(int i=0; i<NMO; i++) v[i] = std::forward_as_tuple(eigVal[i],i);
      std::sort (v.begin(), v.end(), [] (const std::tuple<RealType,int>& a, const std::tuple<RealType,int>& b) {return std::get<0>(a) < std::get<0>(b);});
      for(int i=0; i<NMO; i++)
       for(int j=0; j<NAEA; j++)
        OrbMat(i+NMO,j) = eigVec(i, std::get<1>(v[j]) );
    } else {
      std::copy(OrbMat.begin(),OrbMat.begin()+NAEA*NMO,OrbMat.begin()+NAEA*NMO);
    }

    if(sHam->RHF() && NCA==NCB && NAEA == NAEB && spinRestricted) {
      closed_shell = std::equal(OrbMat.begin(),OrbMat.begin()+NMO*NAEA,OrbMat.begin()+NMO*NAEA);
    }  

  }
  if(closed_shell) {
    app_log()<<"Found closed shell system in GeneralSingleDeterminant. " <<std::endl;
  }

  app_log()<<std::endl <<"*********************************************************************: \n"
           <<" GeneralSingleDeterminant: \n"
           <<"     Number of terms and memory usage of hij:    " <<hij.size() <<"  " <<hij.size()*sizeof(s1D<ValueType>)/1.0e6 <<"  MB. " <<std::endl
           <<"     Number of terms and memory usage of Vijkl:  " <<SMSpHijkl.size() <<"  " <<SMSpHijkl.size()*sizeof(s2D<SPValueType>)/1.0e6 <<"  MB. " <<std::endl; 

  ComplexType e1,e2,o1,o2;
  // Keeping HF and Dalpha/Dbeta in separate places, since HF is used for initialization  

  if(wfn_type==2 && walker_type != 2) {
    HF.resize(2*NMO,NAEA+NAEB);
    HF = OrbMat;
    int wt0=walker_type;
    walker_type=2;
    evaluateLocalEnergy(HF.data(),e1,e2,o1,o2);
    app_log()<<" With <Psi| = <Psi_T(GHF)| \n";
    app_log()<<" Ehf:     " <<std::setprecision(12) <<e1+e2  <<"  \n " //<<ea+eb <<std::endl
             <<" Ekin:     " <<std::setprecision(12) <<e1    <<"  \n " //<<ea <<std::endl
             <<" Epot:     " <<std::setprecision(12) <<e2    <<"  \n "; // <<eb <<std::endl
    app_log()<<" With <Psi| = <Psi_T(UHF)| \n";
    walker_type=wt0;
  } 
  
  // eventually will allow GHF walkers, for now only UHF 
  HF.resize(2*NMO,NAEA);
  // force HF state to be diagonal for now
  for(int i=0; i<NMO; i++)
    for(int j=0; j<NAEA; j++)
      HF(i,j) = OrbMat(i,j); 
  int shft= (wfn_type==2)?NAEA:0;
  for(int i=NMO; i<2*NMO; i++)
    for(int j=0; j<NAEB; j++)
      HF(i,j) = OrbMat(i,j+shft); 
  evaluateLocalEnergy(HF.data(),e1,e2,o1,o2);
  app_log()<<" Ehf:     " <<std::setprecision(12) <<e1+e2  <<"  \n " //<<ea+eb <<std::endl
           <<" Ekin:     " <<std::setprecision(12) <<e1    <<"  \n " //<<ea <<std::endl
           <<" Epot:     " <<std::setprecision(12) <<e2    <<"  \n " // <<eb <<std::endl
           <<"*********************************************************************: \n" <<std::endl <<std::endl;

#ifdef AFQMC_TIMER
    Timer.reset("GeneralSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix");
    Timer.reset("GeneralSingleDeterminant:evaluateLocalEnergy");
#endif

  if(rank() == 0) bool wrote = hdf_write();

  return true;

}

bool GeneralSingleDeterminant::hdf_write()
{
  if(wfn_type==2) {
    app_error()<<" Error in GeneralSlaterDeterminant::hdf_write() not implemented for GHF wavefunctions yet. \n NO RESTART FILE WRITTEN. \n" <<std::endl;
    return true;
  }

  if(hdf_write_file == std::string("")) return true;
  if(readHamFromFile) return true;

  hdf_archive dump(myComm);
  if(!dump.create(hdf_write_file)) {
    app_error()<<" Error opening restart file in GeneralSingleDeterminant. \n"; 
    return false;
  }

  std::string path = "/Wavefunctions/GeneralSingleDeterminant";
  if(dump.is_group( path )) {
    app_error()<<" ERROR: H5Group /Wavefunctions/GeneralSingleDeterminant already exists in restart file. Not over-writing data in file. \n";
    return false;
  }

  dump.push("Wavefunctions");
  dump.push("GeneralSingleDeterminant");

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

// need to write OrbMat and related info (wfn_type, etc)

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
  dump.write(OrbMat,"OrbMat");

  dump.pop();
  dump.pop();

  dump.flush();
  dump.close();

  return true;
}


bool GeneralSingleDeterminant::initFromHDF5(hdf_archive& read,const std::string& tag)
{
  SMSpHijkl.setup(head_of_nodes,"SMSpHijkl",TG.getNodeCommLocal());
  if(head_of_nodes) {

    std::string path = "/Wavefunctions/GeneralSingleDeterminant";
    if(!read.is_group( path )) {
      app_error()<<" ERROR: H5Group /Wavefunctions/GeneralSingleDeterminant does not exist. \n"; 
      return false;
    }

    if(!read.push("Wavefunctions")) return false;
    if(!read.push("GeneralSingleDeterminant")) return false;

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

    OrbMat.resize(2*NMO,NAEA);

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

    if(!read.read(OrbMat,"OrbMat")) return false;
    myComm->bcast(OrbMat);

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

    OrbMat.resize(2*NMO,NAEA);

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

    myComm->bcast(OrbMat);

    myComm->barrier();
    if(!SMSpHijkl.initializeChildren()) return false;
    myComm->barrier();

  }
  readHamFromFile=true;
  return setup_local();
}
    
void GeneralSingleDeterminant::local_evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool diagOnly)
{

  // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
  // lots of simplification because A is an identity matrix with possibly exchanged columns, 
  // look at derivations for information

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  // right now only uhf_walker
  if(wfn_type == 0 || wfn_type==1) {
    // RHF or UHF

    // S0 = transpose(conjg(A))*B
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);

    // S0 = S0^-1       
    ovl_alpha = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

    // SS0 = SlaterMat * S0
    DenseMatrixOperators::product(NMO,NAEA,NAEA,one,SlaterMat,NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);

    // G(alpha) = SS0*transpose(conjg(Dalpha))  
#if defined(AFQMC_SP)
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,temp_density_matrix.data(),NMO);
    std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, mixed_density_matrix.begin());
#else    
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,mixed_density_matrix.data(),NMO);
#endif
    DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data(),NMO);

    if(closed_shell) {
      ovl_beta=ovl_alpha;
      // once this fully works, you will not need the beta sector at all
      std::copy(mixed_density_matrix.begin(),mixed_density_matrix.begin()+NMO*NMO,mixed_density_matrix.begin()+NMO*NMO);   
      return;
    }

    // S1 = transpose(conjg(A))*B  
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NMO*NAEA,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);

    // S0 = S0^-1
    ovl_beta = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

    // SS0(beta) = SlaterMat(beta) * S1
    DenseMatrixOperators::product(NMO,NAEB,NAEB,one,SlaterMat+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

    // G(beta) = SS0*transpose(conjg(Dbeta))  
#if defined(AFQMC_SP)
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,temp_density_matrix.data(),NMO);
    std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, mixed_density_matrix.begin()+NMO*NMO);
#else
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,mixed_density_matrix.data()+NMO*NMO,NMO);
#endif
    DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data()+NMO*NMO,NMO);

  } else {

    // GHF state
    const int ncols = NAEA+NAEB; 

    const ComplexType *ptr=SlaterMat;
    if(walker_type==0) { 
      const ComplexType* ptr1 = SlaterMat;
      ptr = SM.data();
      ComplexType* ptr2 = SM.data();
      for(int i=0; i<NMO; i++,ptr1+=NAEA,ptr2+=ncols) {
        std::copy(ptr1,ptr1+NAEA,ptr2);
        std::copy(ptr1,ptr1+NAEA,ptr2+NAEA+ncols*NMO);
      }
    } else if(walker_type==1) {
      const ComplexType* ptr1 = SlaterMat;
      ptr = SM.data();
      ComplexType* ptr2 = SM.data();
      for(int i=0; i<NMO; i++,ptr1+=NAEA,ptr2+=ncols)
        std::copy(ptr1,ptr1+NAEA,ptr2);
      ptr2 += NAEA;
      for(int i=0; i<NMO; i++,ptr1+=NAEA,ptr2+=ncols)
        std::copy(ptr1,ptr1+NAEB,ptr2);
    } 

    // T0 = transpose(conjg(A))*SM
    DenseMatrixOperators::product_AhB(ncols,ncols,2*NMO,one,OrbMat.data(),ncols,ptr,ncols,zero,T0.data(),ncols);

    // T0 = T0^-1       
    ovl_alpha = Invert(T0.data(), ncols, ncols, Cwork.data(),pivot.data());
    ovl_beta = one;

    // T1 = T*transpose(conjg(OrbMat)) 
    DenseMatrixOperators::product_ABh(ncols,2*NMO,ncols,one,T0.data(),ncols,OrbMat.data(),ncols,zero,T1.data(),2*NMO);
   
    // G = SlaterMat*T1
    if(diagOnly) {
      // return only diagonal block in UHF format/indexing
      // this is used for operations involving 1-body operators in RHF/ROHF/UHF storage format for walker_type<=1 
      // Gaa 
#if defined(AFQMC_SP)
      DenseMatrixOperators::product(NMO,NMO,ncols,one,ptr,ncols,T1.data(),2*NMO,zero,temp_density_matrix.data(),NMO);
      std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, mixed_density_matrix.begin());
#else
      DenseMatrixOperators::product(NMO,NMO,ncols,one,ptr,ncols,T1.data(),2*NMO,zero,mixed_density_matrix.data(),NMO);
#endif
      DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data(),NMO); 

      // Gbb 
#if defined(AFQMC_SP)
      DenseMatrixOperators::product(NMO,NMO,ncols,one,ptr+ncols*NMO,ncols,T1.data()+NMO,2*NMO,zero,temp_density_matrix.data(),NMO);
      std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+NMO*NMO, mixed_density_matrix.begin()+NMO*NMO);
#else
      DenseMatrixOperators::product(NMO,NMO,ncols,one,ptr+ncols*NMO,ncols,T1.data()+NMO,2*NMO,zero,mixed_density_matrix.data()+NMO*NMO,NMO);
#endif
      DenseMatrixOperators::transpose(NMO,mixed_density_matrix.data()+NMO*NMO,NMO); 
    } else {
      // return full matrix [2*NMO,2*NMO]
#if defined(AFQMC_SP)
      DenseMatrixOperators::product(2*NMO,2*NMO,ncols,one,ptr,ncols,T1.data(),2*NMO,zero,temp_density_matrix.data(),2*NMO);
      std::copy(temp_density_matrix.begin(), temp_density_matrix.begin()+4*NMO*NMO, mixed_density_matrix.begin());
#else
      DenseMatrixOperators::product(2*NMO,2*NMO,ncols,one,ptr,ncols,T1.data(),2*NMO,zero,mixed_density_matrix.data(),2*NMO);
#endif
      DenseMatrixOperators::transpose(2*NMO,mixed_density_matrix.data(),2*NMO);
    }

  }
} 

void GeneralSingleDeterminant::local_evaluateOneBodyTrialDensityMatrix()
{

  // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) ) 
  // lots of simplification because A is an identity matrix with possibly exchanged columns, 
  // look at derivations for information

  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0); 

  if(wfn_type == 0 || wfn_type == 1) {
    // S0 = transpose(conjg(A))*B
    DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,OrbMat.data(),NAEA,zero,S0.data(),NAEA);

    // S0 = S0^-1       
    ComplexType ovl = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

    // SS0 = OrbMat * S0
    DenseMatrixOperators::product(NMO,NAEA,NAEA,one,OrbMat.data(),NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);
    // G(alpha) = SS0*transpose(conjg(Dalpha))  
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,OrbMat.data(),NAEA,zero,trial_density_matrix.data(),NMO);

    DenseMatrixOperators::transpose(NMO,trial_density_matrix.data(),NMO);
    if(closed_shell) {
      // once this fully works, you will not need the beta sector at all
      std::copy(trial_density_matrix.begin(),trial_density_matrix.begin()+NMO*NMO,trial_density_matrix.begin()+NMO*NMO);   
      return;
    }

    // S1 = transpose(conjg(A))*B  
    DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NMO*NAEA,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,S1.data(),NAEB);

    // S0 = S0^-1
    ovl *= Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data()); 

    // SS0(beta) = OrbMat(beta) * S1
    DenseMatrixOperators::product(NMO,NAEB,NAEB,one,OrbMat.data()+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

    // G(beta) = SS0*transpose(conjg(Dbeta))  
    DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,OrbMat.data()+NAEA*NMO,NAEA,zero,trial_density_matrix.data()+NMO*NMO,NMO);
    DenseMatrixOperators::transpose(NMO,trial_density_matrix.data()+NMO*NMO,NMO);

  } else {

    // GHF state
    const int ncols = NAEA+NAEB;
    // T0 = transpose(conjg(A))*SM
    DenseMatrixOperators::product_AhB(ncols,ncols,2*NMO,one,OrbMat.data(),ncols,OrbMat.data(),ncols,zero,T0.data(),ncols);

    // T0 = T0^-1       
    ComplexType ovl = Invert(T0.data(), ncols, ncols, Cwork.data(),pivot.data());

    // T1 = T*transpose(conjg(OrbMat)) 
    DenseMatrixOperators::product_ABh(ncols,2*NMO,ncols,one,T0.data(),ncols,OrbMat.data(),ncols,zero,T1.data(),2*NMO);
    // G = SlaterMat*T1
    DenseMatrixOperators::product(2*NMO,2*NMO,ncols,one,OrbMat.data(),ncols,T1.data(),2*NMO,zero,trial_density_matrix.data(),2*NMO);
    DenseMatrixOperators::transpose(2*NMO,trial_density_matrix.data(),2*NMO);

  }

} 
  

  void GeneralSingleDeterminant::evaluateOverlap(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

    const ComplexType one = ComplexType(1.0);
    const ComplexType zero = ComplexType(0.0); 

    if(wfn_type==0 || wfn_type==1) {

      // S0 = transpose(conjg(A))*B
      DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,OrbMat.data(),NAEA,SlaterMat,NAEA,zero,S0.data(),NAEA);

      // S0 = S0^-1       
      ovl_alpha = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

      if(closed_shell) {
        ovl_beta=ovl_alpha;
        return;
      }

      // S1 = transpose(conjg(A))*B  
      DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,OrbMat.data()+NMO*NAEA,NAEA,SlaterMat+NAEA*NMO,NAEA,zero,S1.data(),NAEB);

      // S0 = S0^-1
      ovl_beta = Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

    } else {

      int ncols = NAEA+NAEB; 
      const ComplexType *ptr=SlaterMat;
      if(walker_type==0) {
        // SM = GHF version of SlaterMat 
        const ComplexType* ptr1 = SlaterMat;
        ptr = SM.data();
        ComplexType* ptr2 = SM.data();
        for(int i=0; i<NMO; i++,ptr1+=NAEA,ptr2+=ncols) {
          std::copy(ptr1,ptr1+NAEA,ptr2);
          std::copy(ptr1,ptr1+NAEA,ptr2+NAEA+ncols*NMO);
        }
      } else if(walker_type==1) {
        // SM = GHF version of SlaterMat 
        const ComplexType* ptr1 = SlaterMat;
        ptr = SM.data();
        ComplexType* ptr2 = SM.data();
        for(int i=0; i<NMO; i++,ptr1+=NAEA,ptr2+=ncols)
          std::copy(ptr1,ptr1+NAEA,ptr2);
        ptr2 += NAEA;
        for(int i=0; i<NMO; i++,ptr1+=NAEA,ptr2+=ncols)
          std::copy(ptr1,ptr1+NAEB,ptr2);
      }
    
      // T0 = transpose(conjg(A))*SM
      DenseMatrixOperators::product_AhB(ncols,ncols,2*NMO,one,OrbMat.data(),ncols,ptr,ncols,zero,T0.data(),ncols);

      // T0 = T0^-1       
      ovl_alpha = Invert(T0.data(), ncols, ncols, Cwork.data(),pivot.data());
      ovl_beta = one;

    }

  }

  void GeneralSingleDeterminant::evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n)
  {

#if _DEBUG_AFQMC_
   // assert(SlaterMat.rows() != NMO && SlaterMat.cols() != NMO )
 #endif
#ifdef AFQMC_TIMER
    Timer.start("GeneralSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix"); 
#endif
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,ovl_alpha,ovl_beta,false);
#ifdef AFQMC_TIMER
    Timer.stop("GeneralSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix"); 
#endif

/* For the minimum memory version, store an index that defines the type of conversion. 
 * For example:
 *   n  a1 b1  a2 b2  a3  b3 ..   V
 *   where
 *   n=0 ->   *1   *1   *1     *1
 *   n=1 ->   *1   *-1  conj()    *conj(-1)
 *   n=2 ->   *1    conj()
 *   etc etc     
 */


#ifdef AFQMC_TIMER
    Timer.start("GeneralSingleDeterminant:evaluateLocalEnergy"); 
#endif

    ekin = 0;
    s1Dit end1 = hij.end();
    SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
    for(s1Dit it = hij.begin(); it != end1; it++) 
      ekin += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);

// for closed_shell, you can rewrite this as  sum_ik_jl Ga(ik) * ( 2*<ij|kl> - <ij|lk>  ) * Ga(jl)
// which reduces the cost from ~4M^4 to ~M^4
// but all of the changes happen on the creation of SMSpHijkl, not here 

    epot = 0; 
    const int nr1=1, nc1=SMSpHijkl.cols();
// if(closed_shell) nc1 = NMO*NMO;
    const SPValueType one = SPValueType(1.0);
    const SPValueType zero = SPValueType(0.0);
    SparseMatrixOperators::product_SpMatV(nc1,nc1,one,SMSpHijkl.values(),SMSpHijkl.column_data(),SMSpHijkl.row_index(),mixed_density_matrix.data(),zero,V0.data()); 
    itG = mixed_density_matrix.begin();
    SPComplexVector::iterator itV = V0.begin();
    for(int i=0; i<nc1; i++,++itG,++itV) epot += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG); 
    epot = 0.5*epot+NuclearCoulombEnergy;   

#ifdef AFQMC_TIMER
    Timer.stop("GeneralSingleDeterminant:evaluateLocalEnergy"); 
#endif
    
  }

  void GeneralSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMSpMat& vn, SPComplexSMSpMat& vnT, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("GeneralSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators"); 
#endif
    ComplexType o1,o2;
    const SPComplexType *GF=GG;
    if(needsG) {
     local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2);
     GF = mixed_density_matrix.data();
    }

    if(transposed) {
      SPComplexType one = SPComplexType(1.0);
      const SPComplexType zero = SPComplexType(0.0);
      if(closed_shell) one = SPComplexType(2.0);
      SparseMatrixOperators::product_SpMatV(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.column_data(),vnT.row_index(),GF,zero,v.data());
      if(addBetaBeta && ! closed_shell)
        SparseMatrixOperators::product_SpMatV(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.column_data(),vnT.row_index(),GF+NMO*NMO,one,v.data());
    } else {
      SPValueType one = SPValueType(1.0);
      SPValueType zero = SPValueType(0.0);
      if(closed_shell) one = SPValueType(2.0);
      SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),GF+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("GeneralSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators"); 
#endif
 
  }

  void GeneralSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMVector& vn, SPComplexSMVector& vnT, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("GeneralSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators"); 
#endif
    ComplexType o1,o2;
    const SPComplexType *GF=GG;
    if(needsG) {
     local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2);
     GF = mixed_density_matrix.data();
    }

    if(transposed) {
      SPComplexType one = SPComplexType(1.0);
      const SPComplexType zero = SPComplexType(0.0);
      if(closed_shell) one = SPComplexType(2.0);
      DenseMatrixOperators::product_Ax(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.cols(),GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        DenseMatrixOperators::product_Ax(vnT.rows(),vnT.cols(),one,vnT.values(),vnT.cols(),GF+NMO*NMO,one,v.data());
    } else {
      SPValueType one = SPValueType(1.0);
      SPValueType zero = SPValueType(0.0);
      if(closed_shell) one = SPValueType(2.0);
      DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one,vn.values(),vn.cols(),GF,zero,v.data());
      if(addBetaBeta && !closed_shell)
        DenseMatrixOperators::product_Atx(vn.rows(),vn.cols(),one,vn.values(),vn.cols(),GF+NMO*NMO,one,v.data());
    }

#ifdef AFQMC_TIMER
    Timer.stop("GeneralSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators"); 
#endif

  }

  void GeneralSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int ik0, int ikN, SPValueSMSpMat& vn, SPComplexSMSpMat& vnT, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n)
  {
    APP_ABORT(" Error: Routine not implemented: GeneralSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators. \n\n\n");
  }

  void GeneralSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int ik0, int ikN, SPValueSMVector& vn, SPComplexSMVector& vnT, std::vector<SPComplexType>& v, int walkerBlock, int nW , bool transposed, bool needsG, const int n)
  {
    APP_ABORT(" Error: Routine not implemented: GeneralSingleDeterminant::calculateMixedMatrixElementOfOneBodyOperators. \n\n\n");
  }

  void GeneralSingleDeterminant::calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<SPComplexType> >& v2n, const std::vector<IndexType>& vn_indx, SPValueSMSpMat& vn, std::vector<SPComplexType>& v, const int n)
  {

#ifdef AFQMC_TIMER
    Timer.start("GeneralSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators");
#endif
    ComplexType o1,o2;
    local_evaluateOneBodyMixedDensityMatrix(SlaterMat,o1,o2);


    SPComplexMatrix::iterator itG = mixed_density_matrix.begin();
    for(int i=0; i<vn_indx.size()-1; i++) {
      v[i] = static_cast<SPComplexType>(0.0);
      for(int n = vn_indx[i]; n<vn_indx[i+1]; n++) {
        IndexType ik = Index2Mat(std::get<0>(v2n[n]),std::get<2>(v2n[n]));
        IndexType jl = Index2Mat(std::get<1>(v2n[n]),std::get<3>(v2n[n]));
        IndexType il = Index2Mat(std::get<0>(v2n[n]),std::get<3>(v2n[n]));
        IndexType jk = Index2Mat(std::get<1>(v2n[n]),std::get<2>(v2n[n]));
        v[i] += (*(itG + ik)) * (*(itG + jl) - *(itG + il)) * (*(itG + jk)) * std::get<4>(v2n[n]);
      }
    }

#ifdef AFQMC_TIMER
    Timer.stop("GeneralSingleDeterminant:calculateMixedMatrixElementOfTwoBodyOperators");
#endif

  }

  void GeneralSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMVector& vn, std::vector<SPComplexType>& v, const int n)
  {
    APP_ABORT(" Error: calculateMeanFieldMatrixElementOfOneBodyOperators(Dvn) not implemented yet. \n\n\n");
  }

  void GeneralSingleDeterminant::calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMSpMat& vn, std::vector<SPComplexType>& v, const int n)
  {

    if(trialDensityMatrix_needsupdate) {
      trialDensityMatrix_needsupdate = false;
      ComplexType o1,o2;
      local_evaluateOneBodyTrialDensityMatrix();
      SPtrial_density_matrix = trial_density_matrix;
    }

/*
    app_log()<<"trial_density_matrix: " <<std::endl;
    for(int i=0; i<NMO; i++)
     for(int j=0; j<NMO; j++)
      if(std::abs(trial_density_matrix(i,j)) > 1e-8) app_log()<<i <<" " <<j <<" " <<trial_density_matrix(i,j) <<std::endl;
*/   

    SPValueType one = SPValueType(1.0);
    SPValueType zero = SPValueType(0.0);
    if(closed_shell) one = SPValueType(2.0);
    SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),SPtrial_density_matrix.data(),zero,v.data());
    if(addBetaBeta && !closed_shell)
      SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),SPtrial_density_matrix.data()+NMO*NMO,one,v.data());

  }

}

      
