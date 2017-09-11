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
       core_rank(0),nnodes_per_TG(1),parallel(false),dm_type(0),walker_type(1),wfn_type(0),
       useFacHam(false),sparse_vn(false),dt(1.0),initialDet(1),write_trial_density_matrix("")
    {
    }

    ~WavefunctionBase() {}

    virtual bool init(std::vector<int>& TGdata, SPComplexSMVector *v, hdf_archive& read, const std::string& tag, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm node_heads_comm)
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
      MPI_COMM_HEAD_OF_NODES = node_heads_comm;
      TG.setHeadOfNodesComm(MPI_COMM_HEAD_OF_NODES);
      head_of_nodes = (TG.getCoreID()==0);

      // setup WFN
      if(filetype == "none" && init_type=="ground")   
        return setup_local(); 
      else if(filetype == "fcidump" || filetype == "ascii")
        return initFromAscii(filename);
      else if(filetype == "xml")
        return initFromXML(filename);
      else if(filetype == "hdf5") {
        hdf_archive readF(myComm);
	if(head_of_nodes) 
          if(!readF.open(filename,H5F_ACC_RDONLY)) 
            APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
        if(!initFromHDF5(readF,hdf_read_tag)) {
          app_error()<<" Problems reading hdf5 file in WavefunctionBase::init()";  
          APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
          return false;
        }
        readHamFromFile=true;
        if(head_of_nodes) readF.close();
        return true;
//      } else {
//        if(!initFromHDF5(read,tag)) {
//          app_error()<<" Problems reading restart file in WavefunctionBase::init()"; 
//          APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
//          return false;
//        }
//        readHamFromFile=true;
//        return true;
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

    void setupFactorizedHamiltonian(bool sp, SPValueSMSpMat* spvn_, SPValueSMVector* dvn_, RealType dt_, afqmc::TaskGroup* tg_)
    {
      sparse_vn=sp;
      Spvn=spvn_;
      Dvn=dvn_;
      dt=dt_;
      TG_vn = tg_;
      local_nCholVecs = ((sparse_vn)?(Spvn->cols()):(Dvn->cols())); 
      nCholVecs = ((sparse_vn)?(Spvn->cols()):(Dvn->cols())); 
    }

    bool useFactorizedHamiltonian() { return useFacHam; }

    virtual int sizeOfInfoForDistributedPropagation() 
    {  
      APP_ABORT("WavefunctionBase::sizeOfInfoForDistributedPropagation() not implemented for this wavefunction type. \n");
      return 0;
    }

    virtual void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMSpMat&, std::vector<SPComplexType>& v, const int n=-1 )=0;
    virtual void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMVector&, std::vector<SPComplexType>& v, const int n=-1 )=0;

    // no need to reimplement this in derived class
    void evaluateLocalEnergy(WalkerHandlerBase* wset, bool first , const int n=-1)
    {
      if(parallel) 
        dist_evaluateLocalEnergy(wset,first,n);
      else
        serial_evaluateLocalEnergy(wset,first,n);
    }

    void verifyWalkerData(WalkerHandlerBase* wset, bool first , const int n=-1)
    {
      if(parallel)
        dist_verifyWalkerData(wset,first,n);
      else
        serial_verifyWalkerData(wset,first,n);
    }

    // right now just dumping messages to app_error
    void serial_verifyWalkerData(WalkerHandlerBase* wset, bool first, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ekin,epot,ovlp_a,ovlp_b;
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        evaluateLocalEnergy(wset->getSM(i),ekin,epot,ovlp_a,ovlp_b,n);
        if(first) {
          if(std::abs(wset->getEloc(i)-ekin-epot)>1e-8) app_error()<<" diff in verifyWalkerData: eloc (old-new): "<<wset->getEloc(i)<<" "<<(ekin+epot)<<"\n";
          if(std::abs(wset->getOvlpAlpha(i)-ovlp_a)>1e-8) app_error()<<" diff in verifyWalkerData: ovlp_a (old-new): "<<wset->getOvlpAlpha(i)<<" "<<ovlp_a<<"\n";
          if(std::abs(wset->getOvlpBeta(i)-ovlp_b)>1e-8) app_error()<<" diff in verifyWalkerData: ovlp_b (old-new): "<<wset->getOvlpBeta(i)<<" "<<ovlp_b<<"\n";
        } else {
          if(std::abs(wset->getEloc2(i)-ekin-epot)>1e-8) app_error()<<" diff in verifyWalkerData: eloc2 (old-new): "<<wset->getEloc2(i)<<" "<<(ekin+epot)<<"\n";
          if(std::abs(wset->getOvlpAlpha2(i)-ovlp_a)>1e-8) app_error()<<" diff in verifyWalkerData: ovlp_a2 (old-new): "<<wset->getOvlpAlpha2(i)<<" "<<ovlp_a<<"\n";
          if(std::abs(wset->getOvlpBeta2(i)-ovlp_b)>1e-8) app_error()<<" diff in verifyWalkerData: ovlp_b2 (old-new): "<<wset->getOvlpBeta2(i)<<" "<<ovlp_b<<"\n";
        }
      }
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

    virtual void dist_verifyWalkerData(WalkerHandlerBase* wset, bool first, const int n=-1)
    {  
      serial_verifyWalkerData(wset,first,n);
    }  

    virtual void dist_evaluateLocalEnergy(WalkerHandlerBase* wset, bool first, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ekin,epot,ovlp_a,ovlp_b;
      for(int i=0,cnt=0; i<nw; i++,cnt++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        if(cnt%ncores_per_TG == core_rank) {
          evaluateLocalEnergy(wset->getSM(i),ekin,epot,ovlp_a,ovlp_b,n);
          if(first)
            wset->setWalker(i,ekin+epot,ovlp_a,ovlp_b);
          else {
            wset->setEloc2(i,ekin+epot);
            wset->setOvlp2(i,ovlp_a,ovlp_b);
          }
        }
        cnt++;
      }
    } 

    virtual void evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1)=0;

    virtual void evaluateLocalEnergy(bool addBetaBeta, RealType dt, const ComplexType* SlaterMat, const SPValueSMSpMat& Spvn, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool transposed, const int n=-1) 
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
    {
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

    virtual void evaluateOneBodyMixedDensityMatrix(WalkerHandlerBase* wset, SPComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full=true) {
      APP_ABORT(" Error: evaluateOneBodyMixedDensityMatrix not implemented for this wave-function. \n\n\n");
    }

    virtual void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexMatrix& G)=0;

    virtual void evaluateTwoBodyMixedDensityMatrix()=0;

    virtual void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMSpMat&, SPComplexSMSpMat&, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n=-1)=0;
    virtual void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMVector&, SPComplexSMVector&, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n=-1)=0;

    virtual void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int ik0, int ikN, SPValueSMSpMat&, SPComplexSMSpMat&, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n=-1)=0;
    virtual void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int ik0, int ikN, SPValueSMVector&, SPComplexSMVector&, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n=-1)=0;


    // Two body operators: Not used yet
    virtual void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<SPComplexType> >& vn, const std::vector<IndexType>& vn_indx, SPValueSMSpMat&, std::vector<SPComplexType>& v, const int n=-1 )=0;

    std::string name;
    std::string wfn_role;

    virtual bool check_occ_orbs() {return true; } 

    void setCommBuffer(SPComplexSMVector& bf)
    {
        //commBuff = bf;
    }

    // generate transposed of Spvn
    // done here since storage of density matrix is dependent on wavefunction type 
    // base class implementation creates an exact transpose without modification 
    virtual void generateTransposedOneBodyOperator(bool addBetaBeta, SPValueSMVector& Dvn, SPComplexSMVector& DvnT ) {

      if(!addBetaBeta)
        APP_ABORT("  Error: generateTransposedOneBodyOperator not implemented with UHF integrals.");      
      DvnT.setup(head_of_nodes,"DvnT",TG.getNodeCommLocal());
      DvnT.setDims(Dvn.cols(), Dvn.rows());
      DvnT.resize(Dvn.cols()*Dvn.rows());

      if(head_of_nodes) {
        int nr = DvnT.rows(); // == Dvn.cols() 
        int nc = DvnT.cols(); // == Dvn.rows()
        for(int i=0, ii=0; i<nr; i++)
          for(int j=0; j<nc; j++, ii++)
            DvnT[ii] = Dvn[j*nr+i]; // DvnT(i,j) = Dvn(j,i)
      }
      MPI_Barrier(TG.getNodeCommLocal());

    }

    // base class implementation creates an (almost) exact transpose without modification 
    // The number of columns is NMO*NMO, instead of the default 2*NMO*NMO
    virtual void generateTransposedOneBodyOperator(bool addBetaBeta, SPValueSMSpMat& Spvn, SPComplexSMSpMat& SpvnT) {

      assert(Spvn.size() > 0);
      if(!addBetaBeta)
        APP_ABORT("  Error: generateTransposedOneBodyOperator not implemented with UHF integrals.");
      SpvnT.setup(head_of_nodes,"SpvnT",TG.getNodeCommLocal());
      int NMO2 = NMO*NMO;
      SpvnT.setDims(Spvn.cols(),NMO2);
      SpvnT.resize(Spvn.size());
      if(head_of_nodes) {
        std::copy(Spvn.vals_begin(),Spvn.vals_end(),SpvnT.vals_begin());
        std::copy(Spvn.rows_begin(),Spvn.rows_end(),SpvnT.cols_begin());
        std::copy(Spvn.cols_begin(),Spvn.cols_end(),SpvnT.rows_begin());
      }
      app_log()<<" Compressing transposed Cholesky matrix. \n";
      SpvnT.compress(TG.getNodeCommLocal());
    }

  protected:

    virtual bool setup_local() {
      APP_ABORT(" Error: type=none not allowed for this wavefunction type. \n");
      return false;
    }

    afqmc::TaskGroup TG; 
    bool distribute_Ham;  
    bool parallel;
    int min_ik, max_ik;    
    int core_rank,ncores_per_TG;
    int nnodes_per_TG;

    int initialDet; // 0: RHF, 1: UHF

    std::string filename;
    std::string filetype;
    std::string init_type;

    std::string hdf_write_file;
    std::string hdf_read_tag;
    std::string hdf_write_tag;

    std::string write_trial_density_matrix;

    // in case the coulomb energy is evaluated using the factorized hamiltonian
    bool useFacHam;
    bool sparse_vn;
    RealType dt;
    int local_nCholVecs;
    int nCholVecs;
    SPValueSMSpMat *Spvn;
    SPValueSMVector *Dvn;
    afqmc::TaskGroup* TG_vn;     // task group of the factorized hamiltonian

    bool closed_shell;
    // in both cases below: closed_shell=0, UHF/ROHF=1, GHF=2
    int walker_type;
    int wfn_type;
    // dm_type is the DM type, which should be the largest of walker and wfn types
    int dm_type;

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
