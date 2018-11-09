#include<cstdlib>
#include<memory>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<map>
#include<utility>
#include<vector>
#include<numeric>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "OhmmsPETE/TinyVector.h"
#include "ParticleBase/ParticleAttrib.h"
#include "type_traits/scalar_traits.h"
#include <Platforms/sysutil.h>
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/SimpleParser.h"
#include "Configuration.h"
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory_Helper.h"
#include "AFQMC/Hamiltonians/createHamiltonian_Helper.hpp"

#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#include "AFQMC/Hamiltonians/SymmetricFactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/KPFactorizedHamiltonian.h"
#include "AFQMC/Hamiltonians/KPTHCHamiltonian.h"
#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D.h"
//#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian_old.h"
//#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D_old.h"

#include "AFQMC/Utilities/readHeader.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"
#include "AFQMC/Utilities/Utils.h"

#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"

#include "AFQMC/Matrix/hdf5_readers.hpp"
#include "AFQMC/Matrix/array_partition.hpp"

namespace qmcplusplus
{

namespace afqmc
{

Hamiltonian HamiltonianFactory::fromASCII(GlobalTaskGroup& gTG, xmlNodePtr cur)
{

  using s2Dit = std::vector<s2D<ValueType> >::iterator;
  using s4Dit = SMDenseVector<s4D<ValueType> >::iterator; 

  // used to sort snD values using only indexes 
  _mySort_snD_ mySort;

  // used to identify equal index sets (value is not compared)  
  _myEqv_snD_ myEqv;

  if(cur == NULL)
    APP_ABORT("Error: NULL xml pointer in HamiltonianFactory::parse(). \n");

  std::string info("info0");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(info,"info");
  oAttrib.put(cur);

  if(InfoMap.find(info) == InfoMap.end()) {
    app_error()<<"ERROR: Undefined info in execute block. \n";
    APP_ABORT("ERROR: Undefined info in execute block. \n");
  }

  AFQMCInfo& AFinfo = InfoMap[info];

  int NMO = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;
  int NCA = AFinfo.NCB;
  int NCB = AFinfo.NCB;

  if(NCA > 0 || NCB > 0)
    APP_ABORT(" Error: Frozen core has been temporarily disabled. \n\n\n");

  // defaults
  RealType cutoff1bar = 1e-8;
  RealType cutoff_cholesky = 1e-6;
  std::string hdf_write_type="";
  std::string fileName = "";
  std::string ascii_write_file = "";
  std::string hdf_write_file = "";
  int number_of_TGs = 1;
  int n_reading_cores = -1;

  std::string str("no");
  ParameterSet m_param;
  m_param.add(cutoff1bar,"cutoff_1bar","double");
  m_param.add(cutoff_cholesky,"cutoff_decomp","double");
  m_param.add(cutoff_cholesky,"cutoff_decomposition","double");
  m_param.add(cutoff_cholesky,"cutoff_factorization","double");
  m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
  m_param.add(fileName,"filename","std::string");
  m_param.add(hdf_write_file,"hdf_write_file","std::string");
  m_param.add(hdf_write_type,"hdf_write_type","std::string");
  m_param.add(ascii_write_file,"ascii_write_file","std::string");
  m_param.add(number_of_TGs,"nblocks","int");
  m_param.add(n_reading_cores,"num_io_cores","int");
  m_param.put(cur);

  std::transform(hdf_write_type.begin(),hdf_write_type.end(),hdf_write_type.begin(),(int (*)(int))tolower);
  std::transform(str.begin(),str.end(),str.begin(),(int (*)(int))tolower);

  if(number_of_TGs>1) {
    app_error()<<" Error: number_of_TGs>1 not allowed with fcidump integral input. \n" <<std::endl;
    APP_ABORT("");
  }

  // make or get TG
  number_of_TGs = std::max(1, std::min(number_of_TGs,gTG.getTotalNodes()));
  TaskGroup_& TG = getTG(gTG,number_of_TGs);

  // processor info
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
  int nread = (n_reading_cores<=0)?(ncores):(std::min(n_reading_cores,ncores));
  int head = TG.getGlobalRank() == 0;

  app_log()<<" Initializing Hamiltonian from file: " <<fileName <<std::endl;

  std::ifstream in;
  std::streampos start;
  if(head) {  
    in.open(fileName.c_str());
    if(in.fail()) {
       app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
       APP_ABORT("Problems opening ASCII integral file.\n"); 
    }

    int nmo_=-1;
    int naea_=-1;
    int naeb_=-1;
    int nca_=-1;
    int ncb_=-1;
    int nmax_,netot_,ms2_;
    bool spinRestricted;
    int ISYM;
    std::vector<IndexType> occup_alpha;
    std::vector<IndexType> occup_beta;
    std::vector<IndexType> orbSymm;
    std::vector<IndexType> occupPerSymm_alpha;
    std::vector<IndexType> occupPerSymm_beta;
    bool orderStates; 
    bool factorizedHamiltonian;

    if(!readHeader(in,nmax_,nmo_,netot_,naea_,naeb_,nca_,ncb_,ms2_,spinRestricted,ISYM,occup_alpha,occup_beta,orbSymm,occupPerSymm_alpha,occupPerSymm_beta,orderStates,factorizedHamiltonian)) {
      app_error()<<" Error: Problem with header section of file. \n";
      APP_ABORT(" Error: Problem with header section of file. \n");
    }

    if(NCA!=nca_ || NCB!=ncb_ || NMO!=nmo_ || NAEA!=naea_ || NAEB!=naeb_)
      APP_ABORT(" Error: Parameters in FCIDUMP file differ from those in the input file. \n\n\n");

    start = in.tellg();

  }

  NMO = NMO-NCA;

  HamiltonianTypes htype;
  if(head) htype = peekHamType(in);
  {
    int htype_ = int(htype);
    TG.Global().broadcast_n(&htype_,1,0);
    htype = HamiltonianTypes(htype_);
  }

  if(NCA != NCB) {
    app_error()<<"Error in readFCIDUMP. NCA!=NCB. Not sure how to implement this! \n";
    APP_ABORT("Error in readFCIDUMP. NCA!=NCB. Not sure how to implement this! \n");
  }

  int nOne_core,nOne,nTwo,nTwo_core,nTwo_mixed;
  int nThree,nThree_core,nThree_mixed;

  std::map<IndexType,IndexType> orbMapA; 
  std::map<IndexType,IndexType> orbMapB; 
  orbMapA[0]=0;
  orbMapB[0]=0;
  for(int i=1; i<=NMO; i++) orbMapA[i]=i;   
  for(int i=1; i<=NMO; i++) orbMapB[i]=i+NMO;   

  int n3Vecs=0;

  if(head) {
    countElementsFromFCIDUMP(in,cutoff1bar,true,NMO,NCA,NAEA,NAEB,nOne,nOne_core,nTwo,nTwo_core,nTwo_mixed,nThree,nThree_core,nThree_mixed,orbMapA,orbMapB,n3Vecs);
    if(nOne_core>0 || nTwo_core>0 || nTwo_mixed>0 || nThree_mixed>0 || nThree_core > 0)
      APP_ABORT(" Error: Frozen core has been temporarily disabled. \n\n\n");
  }

  TG.Global().broadcast_n(&nOne,1,0);
  TG.Global().broadcast_n(&nTwo,1,0);
  TG.Global().broadcast_n(&nThree,1,0);
  TG.Global().broadcast_n(&n3Vecs,1,0);

  std::vector<s2D<ValueType> > H1;

  // no hamiltonian distribution yet
  int min_i = 0;
  int max_i = NMO;

  ValueType FrozenCoreEnergy=0, NuclearCoulombEnergy=0;  

  if(htype == s4DInts) {

    if(n3Vecs > 0) {
      app_error()<<"Found three index terms in FCIDUMP. (Only allowed with factorized hamiltonian. Check!!!" <<std::endl;
      APP_ABORT("Found three index terms in FCIDUMP. (Only allowed with factorized hamiltonian. Check!!!");
    }

    SMDenseVector<s4D<ValueType> >  V2(TG.getCoreID()==0,std::string("SparseGeneralHamiltonian_V2"),&TG.Node());

    if(TG.getNodeID() == 0) {
      V2.reserve(nTwo); 
    } else {
      V2.resize(nTwo); 
    }

    if(head) {

      H1.reserve(nOne);

      Timer.reset("Generic");
      Timer.start("Generic");    
      SPValueSMSpMat* dummy;
      NuclearCoulombEnergy = readElementsFromFCIDUMP(in,cutoff1bar,true,NMO,NCA,NAEA,NAEB,H1,&V2,dummy,orbMapA,orbMapB);
      in.close();
      Timer.stop("Generic");
      app_log()<<" -- Time to read ASCII integral file (2nd time after reorder): " <<Timer.average("Generic") <<"\n";

      Timer.reset("Generic");
      Timer.start("Generic");

      std::sort (H1.begin(), H1.end(),mySort);
      s2Dit ith = std::unique(H1.begin(),H1.end(),myEqv);
      H1.resize( std::distance(H1.begin(),ith) );

      std::sort (V2.begin(), V2.end(),mySort);
      s4Dit itV = std::unique(V2.begin(),V2.end(),myEqv);
      V2.resize_serial( std::distance(V2.begin(),itV));

      Timer.stop("Generic");
      app_log()<<" -- Time to sort sparse integral tables: " <<Timer.average("Generic") <<"\n";

    }
    
    {
      int sz = H1.size();
      TG.Global().broadcast_n(&sz,1,0);
      if(!head)         
        H1.resize(sz);
      //TG.Global().broadcast(H1.begin(),H1.end());
      MPI_Bcast(H1.data(),sizeof(s2D<ValueType>),MPI_CHAR,0,&TG.Global()); 
    }
    if(TG.getCoreID()==0) {
      int sz = V2.size();
      TG.Cores().broadcast_n(&sz,1,0);
      if(!head)
        V2.resize_serial(sz);      
      //TG.Cores().broadcast(V2.begin(),V2.end());    
      MPI_Bcast(V2.values(),V2.size()*sizeof(SPValueType),MPI_CHAR,0,&TG.Cores());
    }

    app_log()<<" Memory used by 2-el integral table (on head node): " <<V2.memoryUsage()/1024.0/1024.0 <<" MB. " <<std::endl;

    TG.global_barrier();

    return Hamiltonian();
    //return Hamiltonian(SparseHamiltonian_s4D_old(AFinfo,cur,std::move(H1),std::move(V2),TG,min_i,max_i,NuclearCoulombEnergy,FrozenCoreEnergy));

  } else if(htype == Factorized) {

    if(nTwo > 0) {
      app_error()<<"Found both 3-Index and 4-Index terms in FCIDUMP. Only one form is allowed. " <<std::endl;
      APP_ABORT("");
    }

    SPValueSMSpMat V2(NMO*NMO,n3Vecs,TG.getCoreID()==0,std::string("SparseGeneralHamiltonian_V2"),&TG.Node());

    if(TG.getNodeID() == 0) {
      V2.reserve(nThree);
    } else {
      V2.resize(nThree);
    }

    if(head) {

      H1.reserve(nOne);

      Timer.reset("Generic");
      Timer.start("Generic");
      SMDenseVector<s4D<ValueType> >* dummy;
      NuclearCoulombEnergy = readElementsFromFCIDUMP(in,cutoff1bar,true,NMO,NCA,NAEA,NAEB,H1,dummy,&V2,orbMapA,orbMapB);
      in.close();
      Timer.stop("Generic");
      app_log()<<" -- Time to read ASCII integral file (2nd time after reorder): " <<Timer.average("Generic") <<"\n";

      std::sort (H1.begin(), H1.end(),mySort);
      s2Dit ith = std::unique(H1.begin(),H1.end(),myEqv);
      H1.resize( std::distance(H1.begin(),ith) );

    }

    if(TG.getNodeID() == 0) {
      Timer.reset("Generic");
      Timer.start("Generic");
      V2.compress(&TG.Node());
      Timer.stop("Generic");
      app_log()<<" -- Time to sort sparse integral tables: " <<Timer.average("Generic") <<"\n";
    }

    {
      int sz = H1.size();
      TG.Global().broadcast_n(&sz,1,0);
      if(!head)
        H1.resize(sz);  
      //TG.Global().broadcast(H1.begin(),H1.end());
      MPI_Bcast(H1.data(),sizeof(s2D<ValueType>),MPI_CHAR,0,&TG.Global()); 
    }
    if(TG.getCoreID()==0) {
      //TG.Cores().broadcast(V2.vals_begin(),V2.vals_end());
      //TG.Cores().broadcast(V2.rows_begin(),V2.rows_end());
      //TG.Cores().broadcast(V2.cols_begin(),V2.cols_end());
      //TG.Cores().broadcast(V2.rowIndex_begin(),V2.rowIndex_end());
      MPI_Bcast(V2.values(),V2.size()*sizeof(SPValueType),MPI_CHAR,0,&TG.Cores());
      MPI_Bcast(V2.row_data(),V2.size()*sizeof(int),MPI_CHAR,0,&TG.Cores());
      MPI_Bcast(V2.column_data(),V2.size()*sizeof(int),MPI_CHAR,0,&TG.Cores());
      MPI_Bcast(V2.row_index(),(V2.rows()+1)*sizeof(int),MPI_CHAR,0,&TG.Cores());
    }

    app_log()<<" Memory used by 2-el integral table (on head node): " <<V2.memoryUsage()/1024.0/1024.0 <<" MB. " <<std::endl;

    TG.global_barrier();

    return Hamiltonian();
    //return Hamiltonian(FactorizedSparseHamiltonian_old(AFinfo,cur,std::move(H1),std::move(V2),TG,min_i,max_i,NuclearCoulombEnergy,FrozenCoreEnergy));

  } else {
    APP_ABORT(" Error: Hamiltonian type not implemented in fromASCII(). \n\n\n");
    return Hamiltonian();
  }
}

Hamiltonian HamiltonianFactory::fromHDF5(GlobalTaskGroup& gTG, xmlNodePtr cur)
{

    if(cur == NULL)
      APP_ABORT("Error: NULL xml pointer in HamiltonianFactory::parse(). \n");

    std::string info("info0");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(info,"info");
    oAttrib.put(cur);

    if(InfoMap.find(info) == InfoMap.end()) {
      app_error()<<"ERROR: Undefined info in execute block. \n";
      APP_ABORT("ERROR: Undefined info in execute block. \n");
    }

    AFQMCInfo& AFinfo = InfoMap[info]; 

    int NMO = AFinfo.NMO; 
    int NAEA = AFinfo.NAEA; 
    int NAEB = AFinfo.NAEB; 

    xmlNodePtr curRoot=cur;

    // defaults
    double cutoff1body = 1e-8;
    double cutoff1bar = 1e-8;
    double cutoff_cholesky = 1e-6;
    std::string hdf_write_type="";
    std::string filetype = "undefined";
    std::string fileName = "";
    std::string ascii_write_file = "";
    std::string hdf_write_file = "";
    int number_of_TGs = 1;
    int n_reading_cores = -1;
    bool orderStates=false; 

    std::string order("no");
    std::string str("no");
    ParameterSet m_param;
    m_param.add(order,"orderStates","std::string");
    m_param.add(cutoff1body,"cutoff_1body","double");
    m_param.add(cutoff1bar,"cutoff_1bar","double");
    m_param.add(cutoff_cholesky,"cutoff_decomp","double");
    m_param.add(cutoff_cholesky,"cutoff_decomposition","double");
    m_param.add(cutoff_cholesky,"cutoff_factorization","double");
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(fileName,"filename","std::string");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");
    m_param.add(hdf_write_type,"hdf_write_type","std::string");
    m_param.add(ascii_write_file,"ascii_write_file","std::string");
    m_param.add(number_of_TGs,"nblocks","int");
    m_param.add(n_reading_cores,"num_io_cores","int");
    m_param.put(cur);

    orderStates=false;
    std::transform(order.begin(),order.end(),order.begin(),(int (*)(int))tolower);
    std::transform(filetype.begin(),filetype.end(),filetype.begin(),(int (*)(int))tolower);
    std::transform(hdf_write_type.begin(),hdf_write_type.end(),hdf_write_type.begin(),(int (*)(int))tolower);
    std::transform(str.begin(),str.end(),str.begin(),(int (*)(int))tolower);

    orderStates = (order == "yes" || order == "true");

    Timer.reset("Generic");
    Timer.start("Generic");

    // make or get TG
    number_of_TGs = std::max(1, std::min(number_of_TGs,gTG.getTotalNodes()));
    TaskGroup_& TG = getTG(gTG,number_of_TGs);

    // processor info
    int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
    int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
    int nread = (n_reading_cores<=0)?(ncores):(std::min(n_reading_cores,ncores));
    int head = TG.getGlobalRank() == 0;

    app_log()<<" Initializing Hamiltonian from file: " <<fileName <<std::endl;

    // FIX FIX FIX 
    hdf_archive dump(TG.Global());
    // these cores will read from hdf file
    if( coreid < nread ) {
      if(!dump.open(fileName,H5F_ACC_RDONLY)) {
        app_error()<<" Error opening integral file in SparseGeneralHamiltonian. \n";
        APP_ABORT("");
      }
      if(!dump.push("Hamiltonian",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group not Hamiltonian found. \n"; 
        APP_ABORT("");
      }
    }

    HamiltonianTypes htype;
    if(head) htype = peekHamType(dump);
    {
      int htype_ = int(htype);
      TG.Global().broadcast_n(&htype_,1,0);
      htype = HamiltonianTypes(htype_);
    }

    int int_blocks,nvecs,nkpts=-1;
    std::vector<int> Idata(8);
    if(head) 
      if(!dump.read(Idata,"dims")) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading dims. \n"; 
        APP_ABORT("");
      } 
    TG.Global().broadcast(Idata.begin(),Idata.end());

    int_blocks = Idata[2];
    if(Idata[3] != NMO) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n"; 
      APP_ABORT(" Error: NMO differs from value in integral file. \n");
    }
    if(Idata[4] != NAEA) {
      app_log()<<" WARNING: NAEA differs from value in integral file. \n"; 
//      APP_ABORT(" ");
    }
    if(Idata[5] != NAEB) {
      app_log()<<" WARNING: NAEB differs from value in integral file. \n"; 
//      APP_ABORT(" ");
    }
    nvecs = Idata[7];
    if(htype == KPFactorized || htype == KPTHC) nkpts=Idata[2];

    // 1 body hamiltonian
    std::vector<s2D<ValueType> > H1;

//    std::vector<IndexType> occup_alpha(NAEA);
//    std::vector<IndexType> occup_beta(NAEB);
    ValueType NuclearCoulombEnergy(0);
    ValueType FrozenCoreEnergy(0); 

    if(head) { 
/*
      std::vector<int> occups(NAEA+NAEB);
      if(!dump.read(occups,"occups")) { 
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading occups dataset. \n"; 
        APP_ABORT(" ");
      }
      for(int i=0; i<NAEA; i++) occup_alpha[i] = occups[i];
      for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) occup_beta[j] = occups[i];
*/
      std::vector<ValueType> Rdata(2);
      if(!dump.read(Rdata,"Energies")) { 
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading  dataset. \n"; 
        APP_ABORT(" ");
      }
      if(Rdata.size()>0)  
        NuclearCoulombEnergy = Rdata[0];
      if(Rdata.size()>1)  
        FrozenCoreEnergy = Rdata[1];
    }
    
//    TG.Global().broadcast(occup_alpha.begin(),occup_alpha.end());
//    TG.Global().broadcast(occup_beta.begin(),occup_beta.end());
    TG.Global().broadcast_n(&NuclearCoulombEnergy,1,0);
    TG.Global().broadcast_n(&FrozenCoreEnergy,1,0);


    if(head) { 

      using std::conj;  
      using std::imag;
      bool foundH1=false;
      boost::multi_array<ValueType,2> hcore(extents[NMO][NMO]);
      if(nkpts > 0) {
        // nothing to do, H1 is read during construction of HamiltonianOperations object.
      } else if(dump.read(hcore,"hcore")) {
        foundH1 = true;

        int nnz = 0;
        for(int i=0; i<hcore.shape()[0]; i++) {
          if(std::abs(imag(hcore[i][i])) > 1e-12) {
            app_error()<<" Error: hcore is not hermitian. \n";
            APP_ABORT("");
          }
          if(std::abs(hcore[i][i]) > cutoff1body) nnz++; 
          for(int j=i+1; j<hcore.shape()[1]; j++) {
            if(std::abs(hcore[i][j]-conj(hcore[j][i])) > 1e-12) {
              app_error()<<" Error: hcore is not hermitian. \n";
              APP_ABORT("");
            }
            if(std::abs(hcore[i][j]) > cutoff1body) nnz++;
          }  
        }
        H1.resize(nnz);
        nnz = 0;
        for(int i=0; i<hcore.shape()[0]; i++) {
          if(std::abs(hcore[i][i]) > cutoff1body) 
            H1[nnz++] = std::make_tuple(OrbitalType(i),OrbitalType(i),ValueType(hcore[i][i]));      
          for(int j=i+1; j<hcore.shape()[1]; j++) {
            if(std::abs(hcore[i][j]) > cutoff1body) 
              H1[nnz++] = std::make_tuple(OrbitalType(i),OrbitalType(j),ValueType(hcore[i][j]));      
          }        
        }
      } else {    

        if(Idata[0] < 1) {
          app_error()<<" Error in HamiltonianFactory::fromHDF5(): Dimensions of H1 < 1.  \n";
          APP_ABORT(" ");
        }

        H1.resize(Idata[0]);

        std::vector<OrbitalType> ivec(2*H1.size());    
        if(!dump.read(ivec,"H1_indx")) {
          app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading H1_indx. \n";
          APP_ABORT(" ");
        } 
        for(int i=0, j=0; i<H1.size(); i++, j+=2)        
          H1[i] = std::make_tuple(ivec[j],ivec[j+1],0);  

        std::vector<ValueType> vvec(H1.size());    
        if(!dump.read(vvec,"H1")) {
          app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading H1.  \n";
          APP_ABORT(" ");
        }

        for(int i=0; i<H1.size(); i++) {
          // keep i<=j by default
          if(std::get<0>(H1[i]) <= std::get<1>(H1[i]))
              std::get<2>(H1[i]) = vvec[i];
          else {
              std::swap(std::get<0>(H1[i]),std::get<1>(H1[i]));
              std::get<2>(H1[i]) = myconj(vvec[i]);
          }
        }
      } 

      std::sort (H1.begin(), H1.end(),
        [] (const s2D<ValueType>& lhs, const s2D<ValueType>& rhs)
          {
            return (bool)(std::get<0>(lhs) < std::get<0>(rhs)) ||
              ( !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) &&
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) );
          }
      );

      int sz = H1.size();
      MPI_Bcast(&sz,1,MPI_INT,0,&TG.Global());
      if(sz>0) MPI_Bcast(H1.data(),H1.size()*sizeof(s2D<ValueType>),MPI_CHAR,0,&TG.Global());
    } else {
      int sz;
      MPI_Bcast(&sz,1,MPI_INT,0,&TG.Global());
      H1.resize(sz);  
      if(sz>0) MPI_Bcast(H1.data(),H1.size()*sizeof(s2D<ValueType>),MPI_CHAR,0,&TG.Global());
    } 

    // now read the integrals
    if(htype == KPTHC) {

      if(coreid < nread && !dump.push("KPTHC",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group not KPTHC found. \n";
        APP_ABORT("");
      }
      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();
      return Hamiltonian(KPTHCHamiltonian(AFinfo,cur,std::move(H1),TG,
                                        NuclearCoulombEnergy,FrozenCoreEnergy));

    } else if(htype == KPFactorized) {

      if(coreid < nread && !dump.push("KPFactorized",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group not KPFactorized found. \n";
        APP_ABORT("");
      }
      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();
      // KPFactorizedHamiltonian matrices are read by THCHamiltonian object when needed, 
      // since their ownership is passed to the HamOps object.
      return Hamiltonian(KPFactorizedHamiltonian(AFinfo,cur,std::move(H1),TG,
                                        NuclearCoulombEnergy,FrozenCoreEnergy));

    } else if(htype == THC) {

      if(coreid < nread && !dump.push("THC",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group not THC found. \n";
        APP_ABORT("");
      }
      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();
      // THC matrices are read by THCHamiltonian object when needed, since their ownership is
      // passed to the HamOps object.
      return Hamiltonian(THCHamiltonian(AFinfo,cur,std::move(H1),TG,
                                        NuclearCoulombEnergy,FrozenCoreEnergy));

    } else if(htype == Factorized) {

      if(coreid < nread && !dump.push("Factorized",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group Factorized not found. \n";
        APP_ABORT("");
      }

      if(TG.getNumberOfTGs() > 1)
        APP_ABORT(" Error: Distributed Factorized hamiltonian not yet implemented. \n\n");

      FactorizedSparseHamiltonian::shm_csr_matrix V2_fact = read_V2fact(dump,TG,nread,NMO,nvecs,cutoff1bar,int_blocks);

      Timer.stop("Generic1");
      app_log()<<" -- Time to read move ucsr into csr matrix: "
               <<Timer.average("Generic1") <<"\n";

      app_log()<<" Memory used by factorized 2-el integral table (on head node): " 
               <<(V2_fact.capacity()*(sizeof(ValueType)+sizeof(IndexType)) + V2_fact.shape()[0]*(2*sizeof(std::size_t)))/1024.0/1024.0 <<" MB. " <<std::endl;

      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();

      Timer.stop("Generic");
      app_log()<<" -- Time to initialize Hamiltonian from h5 file: " <<Timer.average("Generic") <<"\n";

      return Hamiltonian(FactorizedSparseHamiltonian(AFinfo,cur,std::move(H1),std::move(V2_fact),TG,NuclearCoulombEnergy,FrozenCoreEnergy));

    } else if(htype == SymmetricFactorized) {

      if(coreid < nread && !dump.push("SymmetricFactorized",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group SymmetricFactorized not found. \n";
        APP_ABORT("");
      }

      if(TG.getNumberOfTGs() > 1)
        APP_ABORT(" Error: Distributed SymmetricFactorized hamiltonian not yet implemented. \n\n");

      SymmetricFactorizedSparseHamiltonian::shm_csr_matrix V2_fact = read_V2fact(dump,TG,nread,NMO,nvecs,cutoff1bar,int_blocks);

      Timer.stop("Generic1");
      app_log()<<" -- Time to read move ucsr into csr matrix: "
               <<Timer.average("Generic1") <<"\n";

      app_log()<<" Memory used by symmetric factorized 2-el integral table (on head node): "
               <<(V2_fact.capacity()*(sizeof(ValueType)+sizeof(IndexType)) + V2_fact.shape()[0]*(2*sizeof(std::size_t)))/1024.0/1024.0 <<" MB. " <<std::endl;

      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();

      Timer.stop("Generic");
      app_log()<<" -- Time to initialize Hamiltonian from h5 file: " <<Timer.average("Generic") <<"\n";

      return Hamiltonian(SymmetricFactorizedSparseHamiltonian(AFinfo,cur,std::move(H1),std::move(V2_fact),TG,NuclearCoulombEnergy,FrozenCoreEnergy));

    } else if(htype == s4DInts) {

      if(coreid < nread && !dump.push("Integrals",false)) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Group not Integrals found. \n";
        APP_ABORT("");
      }

      int min_i=0;
      int max_i=NMO;
      if(TG.getNumberOfTGs()>1 && number_of_TGs > NMO) 
        APP_ABORT("Error: number_of_TGs > NMO. \n\n\n");

      //SMDenseVector<s4D<ValueType> > V2 = read_V2(dump,TG,nread,NMO,min_i,max_i,cutoff1bar,int_blocks); 
      size_t nr = NMO*(NMO+1)/2;
      size_t nc = NMO*NMO;
      SpVType_shm_csr_matrix V2({nr,nc},{0,0},1,
                                    boost::mpi3::intranode::allocator<SPValueType>(TG.Node())); 

      app_log()<<" Memory used by 2-el integral table (on head node): " <<V2.capacity()*(sizeof(SPValueType)+sizeof(int))/1024.0/1024.0 <<" MB. " <<std::endl;

      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();

      Timer.stop("Generic");
      app_log()<<" -- Time to initialize Hamiltonian from h5 file: " <<Timer.average("Generic") <<"\n";

      return Hamiltonian(SparseHamiltonian_s4D(AFinfo,cur,std::move(H1),std::move(V2),TG,NuclearCoulombEnergy,FrozenCoreEnergy));

    }

}
}  // afqmc
}  // qmcplusplus
