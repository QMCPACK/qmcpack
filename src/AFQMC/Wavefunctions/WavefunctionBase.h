#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTIONBASE_H
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTIONBASE_H

#include "AFQMC/config.h"
#include<Message/MPIObjectBase.h>
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "io/hdf_archive.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"

namespace qmcplusplus
{

// Eventually make this a template to handle walker types
class WavefunctionBase: public MPIObjectBase, public AFQMCInfo
{

  typedef WavefunctionBase* WfnPtr;
  typedef HamiltonianBase* HamPtr;

  public:

    WavefunctionBase(Communicate *c):MPIObjectBase(c),readHamFromFile(false),
       hdf_write_file(""),hdf_read_tag(""),hdf_write_tag(""),wfn_role(""),closed_shell(false),
       filetype(""),TG(c,"WavefunctionTG"),distribute_Ham(false),min_ik(-1),max_ik(-1),ncores_per_TG(1),
       core_rank(0),nnodes_per_TG(1),parallel(false)
    {
    }

    ~WavefunctionBase() {}

    void setHeadComm(bool hd, MPI_Comm comm) {
      head_of_nodes=hd;
      MPI_COMM_HEAD_OF_NODES = comm;
    }

    virtual bool init(std::vector<int>& TGdata, ComplexSMVector *v, hdf_archive& read, const std::string& tag, MPI_Comm tg_comm, MPI_Comm node_comm)
    {

      // setup TG
      ncores_per_TG=TGdata[4]; 
      if(nnodes_per_TG > 1) distribute_Ham=true;
      if(nnodes_per_TG > 1 || ncores_per_TG >1) parallel=true;
      if(!TG.quick_setup(ncores_per_TG,nnodes_per_TG,TGdata[0],TGdata[1],TGdata[2],TGdata[3]))
        return false;
      TG.setBuffer(v);
      core_rank = TG.getCoreRank(); 
      TG.setNodeCommLocal(node_comm);
      TG.setTGCommLocal(tg_comm);   

      // setup WFN
      if(filetype == "none" && init_type=="ground")   
        return setup_local(); 
      else if(filetype == "fcidump" || filetype == "ascii" || filetype == "sqc_ascii")
        return initFromAscii(filename);
      else if(filetype == "xml")
        return initFromXML(filename);
      else if(filetype == "hdf5") {
        hdf_archive readF(myComm);
	if(head_of_nodes) 
          if(!readF.open(filename,H5F_ACC_RDONLY,false)) 
            APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
        if(!initFromHDF5(readF,hdf_read_tag)) {
          app_error()<<" Problems reading hdf5 file in WavefunctionBase::init()";  
          APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
          return false;
        }
        readHamFromFile=true;
        if(head_of_nodes) readF.close();
        return true;
      } else {
        if(!initFromHDF5(read,tag)) {
          app_error()<<" Problems reading restart file in WavefunctionBase::init()"; 
          APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
          return false;
        }
        readHamFromFile=true;
        return true;
      }
      app_error()<<" Could not find a wavefunction initialization type. \n";
      return false;
    }

    bool isClosedShell() {return closed_shell;}

    //virtual bool hdf_write(hdf_archive& read, const std::string& tag, bool include_tensors=true)=0;
    virtual bool hdf_write()=0;
   
    virtual bool setup(HamPtr)=0;

    virtual bool parse(xmlNodePtr)=0;

    ComplexMatrix& getHF() { return HF; }

    virtual int sizeOfInfoForDistributedPropagation() 
    {  
      APP_ABORT("WavefunctionBase::sizeOfInfoForDistributedPropagation() not implemented for this wavefunction type. \n");
    }

    virtual void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSpMat&, std::vector<ComplexType>& v, const int n=-1 )=0;
    virtual void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSMSpMat&, std::vector<ComplexType>& v, const int n=-1 )=0;

    // no need to reimplement this in derived class
    void evaluateLocalEnergy(WalkerHandlerBase* wset, bool first , const int n=-1)
    {
      if(parallel) 
        dist_evaluateLocalEnergy(wset,first,n);
      else
        serial_evaluateLocalEnergy(wset,first,n);
    }

    // no need to reimplement this in derived class
    void serial_evaluateLocalEnergy(WalkerHandlerBase* wset, bool first, const int n=-1)
    { 
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ekin,epot,ovlp_a,ovlp_b;   
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        evaluateLocalEnergy(wset->getSM(i),ekin,epot,ovlp_a,ovlp_b,n);
        if(first) 
          wset->setWalker(i,ekin+epot,ovlp_a,ovlp_b);
        else {
          wset->setEloc2(i,ekin+epot);
          wset->setOvlp2(i,ovlp_a,ovlp_b);
        }
      }
    }

    virtual void dist_evaluateLocalEnergy(WalkerHandlerBase* wset, bool first, const int n=-1)
    {  APP_ABORT("WavefunctionBase::dist_evaluateLocalEnergy not implemented for this wavefunction type. \n"); }

    virtual void evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1)=0;

    virtual void evaluateLocalEnergy(bool addBetaBeta, RealType dt, const ComplexType* SlaterMat, const ComplexSMSpMat& Spvn, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool transposed, const int n=-1) 
    {  APP_ABORT("WavefunctionBase::evaluateLocalEnergy with factorized H not implemented for this wavefunction type. \n"); } 

    // no need to reimplement this in derived class
    void evaluateOverlap(WalkerHandlerBase* wset, bool first, const int n=-1)
    {
      if(parallel) 
        dist_evaluateOverlap(wset,first,n);
      else
        serial_evaluateOverlap(wset,first,n);
    }

    virtual void dist_evaluateOverlap(WalkerHandlerBase* wset, bool first, const int n=-1)
    {  APP_ABORT("WavefunctionBase::dist_evaluateOverlap not implemented for this wavefunction type. \n"); }

    // no need to reimplement this in derived class
    void serial_evaluateOverlap(WalkerHandlerBase* wset, bool first, const int n=-1)
    { 
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ovlp_a,ovlp_b;
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        evaluateOverlap(wset->getSM(i),ovlp_a,ovlp_b,n);
        if(first) 
          wset->setOvlp(i,ovlp_a,ovlp_b);
        else 
          wset->setOvlp2(i,ovlp_a,ovlp_b);
      }
    }

    virtual void evaluateOverlap(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1)=0;

    virtual void evaluateOneBodyMixedDensityMatrix(WalkerHandlerBase* wset, ComplexSMVector* buf, int wlksz, int gfoffset, bool full=true) {
      APP_ABORT(" Error: evaluateOneBodyMixedDensityMatrix not implemented for this wave-function. \n\n\n");
    }

    virtual void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexMatrix& G)=0;

    virtual void evaluateTwoBodyMixedDensityMatrix()=0;

    virtual void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1)=0;
    virtual void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSMSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1)=0;

    virtual void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1)=0;
    virtual void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSMSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1)=0;

    virtual void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSpMat&, std::vector<ComplexType>& v, const int n=-1 )=0;
    virtual void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSMSpMat&, std::vector<ComplexType>& v, const int n=-1 )=0;

    std::string name;
    std::string wfn_role;

    virtual bool check_occ_orbs() {return true; } 

    virtual bool isOccupAlpha( int i) { return true; }
    virtual bool isOccupBeta(int i) { return true; }

    void setCommBuffer(ComplexSMVector& bf)
    {
        //commBuff = bf;
    }

  protected:

    virtual bool setup_local() {
      APP_ABORT(" Error: type=none not allowed for this wavefunction type. \n");
      return false;
    }

    TaskGroup TG; 
    bool distribute_Ham;  
    bool parallel;
    int min_ik, max_ik;    
    int core_rank,ncores_per_TG;
    int nnodes_per_TG;

    std::string filename;
    std::string filetype;
    std::string init_type;

    std::string hdf_write_file;
    std::string hdf_read_tag;
    std::string hdf_write_tag;

    bool closed_shell;

    ComplexMatrix HF;

    bool readHamFromFile;

    bool head_of_nodes;
    MPI_Comm MPI_COMM_HEAD_OF_NODES;

    virtual bool initFromAscii(std::string fileName)=0;  

    virtual bool initFromXML(std::string fileName)=0;  

    virtual bool initFromHDF5(hdf_archive&,const std::string&)=0;  

    virtual bool getHamiltonian(HamPtr )=0;  

    // used to identify the current step.
    // The main purpose of this is to tellthe different
    // wavefunction objects whether we are in the same step
    // or not. This will allow us to reuse information already
    // calculated in a previous section of the current step.
    // e.g. Not to recalculate density matrices if we are redoing
    // a local energy calculation on the same step. 
    // Specially useful for multideterminant calculations
    int time_stamp;

    // Hamiltonian object associated with these wavefuncitons  
    // The hamiltonian gives access to matrix elements, 
    // but trial wavefunction related quantities are calculated
    // in this class.
    // Make sure that there is consistency between the type of hamiltonian
    // and the type of wavefunction object, e.g. sparse versus full matrix 
    HamPtr ham0; 
  
};
}

#endif
