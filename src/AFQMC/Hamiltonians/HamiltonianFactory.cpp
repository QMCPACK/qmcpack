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
#include "io/hdf_multi.h"
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Hamiltonians/HamiltonianFactory_Helper.h"

#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/KPFactorizedHamiltonian.h"
#include "AFQMC/Hamiltonians/KPTHCHamiltonian.h"

#include "AFQMC/Utilities/readHeader.h"
#include "AFQMC/Utilities/Utils.hpp"

#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"

#include "AFQMC/Matrix/hdf5_readers.hpp"
#include "AFQMC/Matrix/array_partition.hpp"

namespace qmcplusplus
{

namespace afqmc
{

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
      if(!dump.readEntry(Idata,"dims")) {
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

    // 1 body hamiltonian: Why isn't this in shared memory!!!
    boost::multi::array<ValueType,2> H1({NMO,NMO});

    ValueType NuclearCoulombEnergy(0);
    ValueType FrozenCoreEnergy(0);

    if(head) {
      std::vector<ValueType> Rdata(2);
      if(!dump.readEntry(Rdata,"Energies")) {
        app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading  dataset. \n";
        APP_ABORT(" ");
      }
      if(Rdata.size()>0)
        NuclearCoulombEnergy = Rdata[0];
      if(Rdata.size()>1)
        FrozenCoreEnergy = Rdata[1];
    }

    TG.Global().broadcast_n(&NuclearCoulombEnergy,1,0);
    TG.Global().broadcast_n(&FrozenCoreEnergy,1,0);

    if(head) {

      using std::conj;
      using std::imag;
      bool foundH1=false;
      if(nkpts > 0) {
        // nothing to do, H1 is read during construction of HamiltonianOperations object.
      } else if(dump.readEntry(H1,"hcore")) {
        foundH1 = true;
      } else {

        H1.reextent({NMO,NMO});
        if(Idata[0] < 1) {
          app_error()<<" Error in HamiltonianFactory::fromHDF5(): Dimensions of H1 < 1.  \n";
          APP_ABORT(" ");
        }

        std::vector<OrbitalType> ivec(2*Idata[0]);
        if(!dump.readEntry(ivec,"H1_indx")) {
          app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading H1_indx. \n";
          APP_ABORT(" ");
        }
        std::vector<ValueType> vvec(Idata[0]);
        if(!dump.readEntry(vvec,"H1")) {
          app_error()<<" Error in HamiltonianFactory::fromHDF5(): Problems reading H1.  \n";
          APP_ABORT(" ");
        }

        for(int i=0; i<Idata[0]; i++) {
          // keep i<=j by default
          if(ivec[i] <= ivec[i]) {
              H1[ivec[2*i]][ivec[2*i+1]] = vvec[i]; 
              H1[ivec[2*i+1]][ivec[2*i]] = conj(vvec[i]); 
          } else {
              H1[ivec[2*i]][ivec[2*i+1]] = conj(vvec[i]); 
              H1[ivec[2*i+1]][ivec[2*i]] = vvec[i]; 
          }
        }
      }
    }
    TG.Global().broadcast_n(to_address(H1.origin()),H1.num_elements(),0);

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

      app_log()<<" Memory used by factorized 2-el integral table (on head node): "
               <<(V2_fact.capacity()*(sizeof(ValueType)+sizeof(IndexType)) + V2_fact.size(0)*(2*sizeof(std::size_t)))/1024.0/1024.0 <<" MB. " <<std::endl;

      if( coreid < nread ) {
        dump.pop();
        dump.pop();
        dump.close();
      }
      TG.global_barrier();

      return Hamiltonian(FactorizedSparseHamiltonian(AFinfo,cur,std::move(H1),std::move(V2_fact),TG,NuclearCoulombEnergy,FrozenCoreEnergy));
    } else {
      app_error()<<" Error in HamiltonianFactory::fromHDF5(): Unknown Hamiltonian Type. \n";
      APP_ABORT("");
    }

}
}  // afqmc
}  // qmcplusplus
