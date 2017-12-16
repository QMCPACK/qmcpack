#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTIONHANDLER_H
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTIONHANDLER_H

#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include<Message/MPIObjectBase.h>
#include "AFQMC/Wavefunctions/WavefunctionBase.h"
#include "AFQMC/Wavefunctions/PureSingleDeterminant.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
//#include "AFQMC/Walkers/SlaterDetWalker.h"

namespace qmcplusplus
{

// Eventually make this a template to handle walker types
class WavefunctionHandler: public MPIObjectBase, public AFQMCInfo
{

  typedef WavefunctionBase* WfnPtr;
  typedef HamiltonianBase* HamPtr;

  public:

    WavefunctionHandler(Communicate *c):MPIObjectBase(c),name(""),ham0(NULL),phaselessWfn(NULL),ImpSampWfn(NULL),EstimatorWfn(NULL),distribute_Ham(false),core_rank(0),ncores_per_TG(1),new_algo(false)
    {
      wfns.reserve(10);
    }

    ~WavefunctionHandler() {}

    bool parse(xmlNodePtr cur);

    ComplexMatrix& getHF() {
      return ImpSampWfn->getHF();
    }

    bool init(std::vector<int>& TGdata, SPComplexSMVector *v,hdf_archive&,const std::string&, MPI_Comm, MPI_Comm, MPI_Comm); 

    bool setup(HamPtr); 

    WfnPtr addWfn(xmlNodePtr cur);

    void evaluateMeanFields() {}

    void setupFactorizedHamiltonian(bool sp, SPValueSMSpMat* spvn_, SPValueSMVector* dvn_, RealType dt_, afqmc::TaskGroup* tg_)
    {
      for(int i=0; i<wfns.size(); i++) 
        wfns[i]->setupFactorizedHamiltonian(sp,spvn_,dvn_,dt_,tg_);
    }

    inline int sizeOfInfoForDistributedPropagation(const std::string& type) {
      if(type == std::string("ImportanceSampling")) {
        return ImpSampWfn->sizeOfInfoForDistributedPropagation();
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          return EstimatorWfn->sizeOfInfoForDistributedPropagation();
        } 
        APP_ABORT("Undefined wavefunction in sizeOfInfoForDistributedPropagation('Estimator') \n\n\n");
      } else {
        APP_ABORT("Unknown wavefunction type in sizeOfInfoForDistributedPropagation(). \n");
      }
      return 0;
    } 

 
 
    inline void evaluateLocalEnergyAndOverlap(const std::string& type, const int n, WalkerHandlerBase* wset)
    {
        if(type == std::string("ImportanceSampling")) {
          ImpSampWfn->evaluateLocalEnergy(wset,true,n);
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            EstimatorWfn->evaluateLocalEnergy(wset,false,n);
          } else {
            int nw = wset->numWalkers(true);
            for(int i=0; i<nw; i++)
              wset->setEloc2(i,ComplexType(0,0));
          }
        } else {
          APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap(wset). \n");
        }      
    }

    inline void verifyWalkerData(const std::string& type, const int n, WalkerHandlerBase* wset)
    {
        if(type == std::string("ImportanceSampling")) {
          ImpSampWfn->verifyWalkerData(wset,true,n);
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            EstimatorWfn->verifyWalkerData(wset,false,n);
          }
        } else {
          APP_ABORT("Unknown wavefunction type in verifyWalkerData(wset). \n");
        }
    }

    inline void evaluateLocalEnergyAndOverlap(const std::string& type, const int n, ComplexType* SM, ComplexType& eloc, ComplexType& ovlp_a, ComplexType& ovlp_b)
    {
      ComplexType ekin,epot;
      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->evaluateLocalEnergy(SM,ekin,epot,ovlp_a,ovlp_b,n);
        eloc = ekin+epot;
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          EstimatorWfn->evaluateLocalEnergy(SM,ekin,epot,ovlp_a,ovlp_b,n);
          eloc = ekin+epot;
        } else {
          eloc = 0; 
        }  
      } else {
        APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap(SM). \n"); 
      }
    } 

    inline void evaluateOverlap(const std::string& type, const int n, WalkerHandlerBase* wset) 
    {
        if(type == std::string("ImportanceSampling")) {
          ImpSampWfn->evaluateOverlap(wset,true,n);
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            EstimatorWfn->evaluateOverlap(wset,false,n);
          }   
        } else {
          APP_ABORT("Unknown wavefunction type in evaluateOverlap(wset). \n");
        } 
    }

    inline void evaluateOverlap(const std::string& type, const int n, ComplexType* SM, ComplexType& ovlp_a, ComplexType& ovlp_b)
    {
      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->evaluateOverlap(SM,ovlp_a,ovlp_b,n);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->evaluateOverlap(SM,ovlp_a,ovlp_b,n);
      } else {
        APP_ABORT("Unknown wavefunction type in evaluateOverlap(SM). \n");
      }
    }

    inline void evaluateOneBodyMixedDensityMatrix(const std::string& type, WalkerHandlerBase* wset, SPComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full=true) 
    {
      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->evaluateOneBodyMixedDensityMatrix(wset,buf,wlksz,gfoffset,transposed,full);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->evaluateOneBodyMixedDensityMatrix(wset,buf,wlksz,gfoffset,transposed,full);
      } else {
        APP_ABORT("Unknown wavefunction type in evaluateOneBodyMixedDensityMatrix. \n");
      }
    }

    template<class T1, class T2>
    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const std::string& type, const int n, ComplexType* SM, const SPComplexType* GG, T1& Spvn, T2& SpvnT , std::vector<SPComplexType>& v, bool transposed , bool needsG) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMixedMatrixElementOfOneBodyOperators(addBetaBeta,SM,GG,Spvn,SpvnT,v,transposed,needsG,n);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->calculateMixedMatrixElementOfOneBodyOperators(addBetaBeta,SM,GG,Spvn,SpvnT,v,transposed,needsG,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMixedMatrixElementOfOneBodyOperators. \n");
      }

    } 

    template<class T1, class T2>
    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const std::string& type, const int n, SPComplexType* buff, int ik0, int ikN, T1& Spvn, T2& SpvnT, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed , bool needsG) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(addBetaBeta,buff,ik0,ikN,Spvn,SpvnT,v,walkerBlock,nW,transposed,needsG,n);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(addBetaBeta,buff,ik0,ikN,Spvn,SpvnT,v,walkerBlock,nW,transposed,needsG,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMixedMatrixElementOfOneBodyOperators. \n");
      }

    }

    template<class T>
    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, const std::string& type, const int n, T& Spvn, std::vector<SPComplexType>& v ) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMeanFieldMatrixElementOfOneBodyOperators(addBetaBeta,Spvn,v,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMeanFieldMatrixElementOfOneBodyOperators. \n");
      }

    } 

    template<class T1, class T2>
    void generateTransposedOneBodyOperator(bool addBetaBeta, const std::string& type, T1& Spvn, T2& SpvnT ) {
      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->generateTransposedOneBodyOperator(addBetaBeta,Spvn,SpvnT);
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          EstimatorWfn->generateTransposedOneBodyOperator(addBetaBeta,Spvn,SpvnT);
        } else {
          APP_ABORT("Error: Attempting to access uninitialized Estimator wavefunction \n");
        }
      } else {
        APP_ABORT("Unknown wavefunction type in generateTransposedOneBodyOperator. \n");
      }
    }

    bool check_initialized(const std::string& type)
    {
      if(type == std::string("ImportanceSampling")) {
         return ImpSampWfn!=NULL;
      } else if(type == std::string("Estimator")) {
        return EstimatorWfn!=NULL; 
      }
      return false;
    } 

    void setCommBuffer(std::vector<ComplexType>& bf)
    {
      //for(int i=0; i<wfns.size(); i++) 
        //wfns[i]->setCommBuffer(bf);
    }

    bool check_occ_orbs() {
      return ImpSampWfn->check_occ_orbs();
    }

    std::string name;

    int core_rank,ncores_per_TG;

    bool head_of_nodes;
    MPI_Comm MPI_COMM_HEAD_OF_NODES; 

    bool distribute_Ham;  // implement assuming factorized Ham first  

    ComplexVector local_energy;

    bool new_algo;

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
  
    WfnPtr  phaselessWfn;

    WfnPtr  ImpSampWfn;

    WfnPtr  EstimatorWfn;

    // stores pointers to WavefunctionBase objects owned by this object.
    // This is necessary in case wfns are repeated 
    std::vector<WfnPtr> wfns;

};
}

#endif
