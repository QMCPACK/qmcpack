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

    void setHeadComm(bool hd, MPI_Comm comm) {
      head_of_nodes=hd;
      MPI_COMM_HEAD_OF_NODES = comm;
    }

    bool init(std::vector<int>& TGdata, ComplexSMVector *v,hdf_archive&,const std::string&, MPI_Comm, MPI_Comm); 

    bool setup(HamPtr); 

    WfnPtr addWfn(xmlNodePtr cur);

    void evaluateMeanFields() {}

    bool isClosedShell(const std::string& type) {
      if(type == std::string("ImportanceSampling")) {
        return ImpSampWfn->isClosedShell();
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          return EstimatorWfn->isClosedShell();
        }
        APP_ABORT("Undefined wavefunction in isClosedShell('Estimator') \n\n\n");
      } else {
        APP_ABORT("Unknown wavefunction type in isClosedShell(). \n");
      }
      return  false;
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

/*
      if(new_algo) {
        if(type == std::string("ImportanceSampling")) {
          ImpSampWfn->evaluateLocalEnergy(wset,n);
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            EstimatorWfn->evaluateLocalEnergy(wset,n);
          } else {
            for(int i=0; i<nw; i++)
              wset->setEloc2(i,ComplexType(0,0));
          }
        } else {
          APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap. \n");
        }
      } else {  
        if(type == std::string("ImportanceSampling")) {
          for(int i=0; i<nw; i++) {
            if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue; 
            ImpSampWfn->evaluateLocalEnergy(wset->getSM(i),ekin,epot,ovlp_a,ovlp_b,n);
            wset->setWalker(i,ekin+epot,ovlp_a,ovlp_b); 
          } 
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            for(int i=0; i<nw; i++) {
              if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue; 
              EstimatorWfn->evaluateLocalEnergy(wset->getSM(i),ekin,epot,ovlp_a,ovlp_b,n);
              wset->setEloc2(i,ekin+epot); 
              wset->setOvlp2(i,ovlp_a,ovlp_b); 
            } 
          } else {
            for(int i=0; i<nw; i++) 
              wset->setEloc2(i,ComplexType(0,0)); 
          }
        } else {
          APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap. \n");
        }
      }
*/
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

    inline void evaluateLocalEnergyAndOverlap(bool addBetaBeta, const std::string& type, const int n, ComplexType* SM, ComplexType& eloc, ComplexType& ovlp_a, ComplexType& ovlp_b, const ComplexSMSpMat& Spvn, bool transposed, RealType dt)
    {
      ComplexType ekin,epot;
      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->evaluateLocalEnergy(addBetaBeta,dt,SM,Spvn,ekin,epot,ovlp_a,ovlp_b,transposed,n);
        eloc = ekin+epot;
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          EstimatorWfn->evaluateLocalEnergy(addBetaBeta,dt,SM,Spvn,ekin,epot,ovlp_a,ovlp_b,transposed,n);
          eloc = ekin+epot;
        } else {
          eloc = 0;
        }
      } else {
        APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap(Spvn). \n");
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

/*
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ovlp_a,ovlp_b;
      if(new_algo) {
        if(type == std::string("ImportanceSampling")) {
          //ImpSampWfn->dist_evaluateOverlap(wset,n);
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            //EstimatorWfn->dist_evaluateOverlap(wset,n);
          }   
        } else {
          APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap. \n");
        } 
      } else { 
        if(type == std::string("ImportanceSampling")) {
          for(int i=0; i<nw; i++) {
            if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
            ImpSampWfn->evaluateOverlap(wset->getSM(i),ovlp_a,ovlp_b,n);
            wset->setOvlp(i,ovlp_a,ovlp_b);
          }
        } else if(type == std::string("Estimator")) {
          if(EstimatorWfn!=NULL) {
            for(int i=0; i<nw; i++) {
              if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
              EstimatorWfn->evaluateOverlap(wset->getSM(i),ovlp_a,ovlp_b,n);
              wset->setOvlp2(i,ovlp_a,ovlp_b);
            }
          }
        } else {
          APP_ABORT("Unknown wavefunction type in evaluateLocalEnergyAndOverlap. \n");
        }
      }
*/
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

    inline void evaluateOneBodyMixedDensityMatrix(const std::string& type, WalkerHandlerBase* wset, ComplexSMVector* buf, int wlksz, int gfoffset, bool full=true) 
    {
      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->evaluateOneBodyMixedDensityMatrix(wset,buf,wlksz,gfoffset,full);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->evaluateOneBodyMixedDensityMatrix(wset,buf,wlksz,gfoffset,full);
      } else {
        APP_ABORT("Unknown wavefunction type in evaluateOneBodyMixedDensityMatrix. \n");
      }
    }

//    ComplexType evaluateLocalEnergy(const std::string& type, const int n, ComplexMatrix& SD) {}

//    ComplexType evaluateOverlap(const std::string& type, const int n, ComplexMatrix& SD) {}

    template<class T>
    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const std::string& type, const int n, ComplexType* SM, T& Spvn , std::vector<ComplexType>& v, bool transposed , bool needsG) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMixedMatrixElementOfOneBodyOperators(addBetaBeta,SM,Spvn,v,transposed,needsG,n);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->calculateMixedMatrixElementOfOneBodyOperators(addBetaBeta,SM,Spvn,v,transposed,needsG,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMixedMatrixElementOfOneBodyOperators. \n");
      }

    } 

    template<class T>
    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const std::string& type, const int n, ComplexType* buff, int ik0, int ikN, int pik0, T& Spvn , std::vector<ComplexType>& v, bool transposed , bool needsG) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(addBetaBeta,buff,ik0,ikN,pik0,Spvn,v,transposed,needsG,n);
      } else if(type == std::string("Estimator")) {
        EstimatorWfn->calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(addBetaBeta,buff,ik0,ikN,pik0,Spvn,v,transposed,needsG,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMixedMatrixElementOfOneBodyOperators. \n");
      }

    }

    template<class T>
    void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const std::string& type, const int n, ComplexType* SM, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, T&Spvn, std::vector<ComplexType>& v ) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMixedMatrixElementOfTwoBodyOperators(addBetaBeta,SM,vn,vn_indx,Spvn,v,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMixedMatrixElementOfTwoBodyOperators. \n");
      }
    }


    template<class T>
    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, const std::string& type, const int n, T& Spvn, std::vector<ComplexType>& v ) {

      if(type == std::string("ImportanceSampling")) {
        ImpSampWfn->calculateMeanFieldMatrixElementOfOneBodyOperators(addBetaBeta,Spvn,v,n);
      } else {
        APP_ABORT("Unknown wavefunction type in calculateMeanFieldMatrixElementOfOneBodyOperators. \n");
      }

    } 

    bool isOccupAlpha( const std::string& type, int i) {
      if(type == std::string("ImportanceSampling")) {
        return ImpSampWfn->isOccupAlpha(i);
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          return EstimatorWfn->isOccupAlpha(i);
        } else {
          APP_ABORT("Error: Attempting to access uninitialized Estimator wavefunction \n");
        }
      } else {
        APP_ABORT("Unknown wavefunction type in isOccupAlpha. \n");
      }
      return false;
    } 
 
    bool isOccupBeta( const std::string& type, int i) {
      if(type == std::string("ImportanceSampling")) {
        return ImpSampWfn->isOccupBeta(i);
      } else if(type == std::string("Estimator")) {
        if(EstimatorWfn!=NULL) {
          return EstimatorWfn->isOccupBeta(i);
        } else {
          APP_ABORT("Error: Attempting to access uninitialized Estimator wavefunction \n");
        }
      } else {
        APP_ABORT("Unknown wavefunction type in isOccupBeta. \n");
      }
      return false;
    }

    bool check_initialized(const std::string& type)
    {
      if(type == std::string("ImportanceSampling")) {
         ImpSampWfn!=NULL;
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
