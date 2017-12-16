
#ifndef QMCPLUSPLUS_AFQMC_PHASELESS_WITHIMPSAMPLWITHELOC_FORCEBIAS
#define QMCPLUSPLUS_AFQMC_PHASELESS_WITHIMPSAMPLWITHELOC_FORCEBIAS

#include <cmath>
#include <cfloat>

#include "AFQMC/config.h"
#include <Message/MPIObjectBase.h>
#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
//#include "AFQMC/Walkers/SlaterDetWalker.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"

#include "AFQMC/Sandbox/compare_libraries.h"

namespace qmcplusplus
{


class phaseless_ImpSamp_ForceBias: public PropagatorBase
{

  typedef phaseless_ImpSamp_ForceBias  thisClass;
  typedef phaseless_ImpSamp_ForceBias*  thisClassPtr;

  public:
       
  phaseless_ImpSamp_ForceBias(Communicate *c,  RandomGenerator_t* r) : PropagatorBase(c,r), substractMF(true),use_eig(false),first(true),max_weight(100),apply_constrain(true),save_memory(false),vbias_bound(3.0),imp_sampl(true),hybrid_method(false),eloc_from_Spvn(false),sizeOfG(0),walkerBlock(0),test_cnter(0),cutoff(1e-6),fix_bias(0),n_reading_cores(-1),intgs_per_block(2000000)
  {
  } 

  ~phaseless_ImpSamp_ForceBias() {}

  void Propagate(int steps, int& steps_total, WalkerHandlerBase*, RealType& E1);

  bool setup(std::vector<int>&,SPComplexSMVector*,HamiltonianBase*,WavefunctionHandler*,RealType, hdf_archive&, const std::string&,MPI_Comm,MPI_Comm,MPI_Comm);

  bool parse(xmlNodePtr);   

  bool hdf_write(hdf_archive&, const std::string&);

  bool hdf_read(hdf_archive&, const std::string&);

  bool hdf_write_transposed(hdf_archive&, const std::string&);

  bool hdf_read_transposed(hdf_archive&, const std::string&);

  void benchmark(std::string&,int,int,int,WalkerHandlerBase*);

  SPValueSMVector* getDvn() { return &Dvn; }

  SPValueSMSpMat* getSpvn() { return &Spvn; }

  private:

  std::ofstream out_debug;

  bool hybrid_method;

  bool eloc_from_Spvn;

  bool imp_sampl;

  bool substractMF;

  bool use_eig;

  bool first;

  bool save_memory;

  int fix_bias;

  std::ifstream in_rand;

  RealType cutoff;

  int n_reading_cores;

  int intgs_per_block;

  // one-body operator that defines H1
  // stored in sparse form
  std::vector<s2D<ValueType> > H1;  
  IndexType H1_nalpha;  

  // propagator for the one-body hamiltonian plus any other one-body term 
  // from the HS transformation and/or mean-field removal
  std::vector<s2D<ComplexType> > Propg_H1;
  std::vector<IndexType> Propg_H1_indx;

  // HS potential: sum_n (sigma_n - vbias_n) * v2_n
  // This is stored full
  SPComplexMatrix SPvHS;
  ComplexMatrix vHS;
  ComplexMatrix PHS;

  // potentials that represent the decomposition of Vijkl into a quadratic form
  // Vijkl = -0.5 * sum_n (v2_n)^2 + v0, where v0 is stored somewhere else
  // stored in sparse form and sequentially
  // vn_nterms is a vector that contains the number of non-zero terms (above cutoff)
  // in each term in vn. This allows for easier parsing. 
  //std::vector<s2D<ComplexType> > vn;  
  //std::vector<IndexType> vn_indx;
  //Vector< Vector<s2D<ValueType> >::iterator > vn_bounds;

  int GlobalSpvnSize; 
  std::vector<int> nCholVec_per_node;

  ComplexMatrix vn0;
  SPValueSMSpMat Spvn; 
  SPComplexSMSpMat SpvnT; 
  // this is a pointer to either Spvn or SpvnT, to avoid extra logic in code  
  SPComplexSMSpMat *Spvn_for_onebody; 

  // storing cholesky vectors in dense format as a vector,
  // to avoid having to write a shared memory matrix class.
  // I only use this through library routines that take the vector,
  // so it is ok
  SPValueSMVector Dvn;
  SPComplexSMVector DvnT;
  SPComplexSMVector *Dvn_for_onebody;

  // storage for fields
  std::vector<SPRealType> sigma;
  std::vector<SPComplexType> CV0;

  // Force bias potential, typically mixed potential but can also be MF potential
  std::vector<SPComplexType> vbias;   

  // Mean-field subtraction of vHS
  std::vector<SPComplexType> vMF;   

  // used to calculate mean fields, overlaps, local energies, etc  
  WavefunctionHandler* wfn; 

  // local storage
  ComplexMatrix T1; 
  ComplexMatrix T2; 

  ComplexMatrix S1; 
  ComplexMatrix S2; 

  SPComplexSMVector local_buffer;

  bool walkerBlock;
 
  int test_cnter;

  std::vector<ComplexType> MFfactor;
  std::vector<ComplexType> hybrid_weight;

  std::vector<ComplexType> MFs; 
  std::vector<ComplexType> HWs;

  RealType dEloc;
  ComplexType max_weight;
  RealType vbias_bound;

  bool transposed_walker_buffer;
  bool transposed_generated_vhs;

  int sizeOfG;
  int nCholVecs;  // total number of Cholesky Vectors in the calculation  
  int cvec0,cvecN; // index of first and last+1 Cholesky Vector of this TG 
  int cvec0_loc,cvecN_loc; // if save_memory=false index of first and last+1 Cholesky Vector (relative to cvec0) of this core 
  int vHS_size;
  std::vector<int> walker_per_node;

  // ik breakup of Spvn
  IndexType ik0, ikN;   //  minimum and maximum values of ik index in Spvn
  IndexType pik0, pikN;  // locations of bounds of [ik0,ikN] sector in Spvn  

  void serial_propagation_single_step(WalkerHandlerBase*, RealType& E1);

  void dist_propagation_single_step(WalkerHandlerBase*, RealType& E1);

  // in this case, vbias is calculated for step 0 and fixed for all subsequent sub steps
  void serial_propagation_multiple_steps(int steps, WalkerHandlerBase*, RealType& E1);

  void dist_propagation_multiple_steps(int steps, WalkerHandlerBase*, RealType& E1);

  bool apply_constrain;

  void applyHSPropagator(ComplexMatrix&, ComplexMatrix&, ComplexType& factor, int order=-1, bool calculatevHS=true); 

  void addvHS(SPComplexSMVector *buff, int nw, int sz, WalkerHandlerBase* wset); 

  void addvHS_multiple(int n, SPComplexSMVector *buff, int nw, int sz, WalkerHandlerBase* wset); 

  void sampleGaussianFields();

  void sampleGaussianFields(SPComplexType*,int);

  inline IndexType Index2Mat(const IndexType I, const IndexType J) const {
    return (J<NMO)?(I*NMO+J):(I*NMO+J-NMO);
  }

  inline ComplexType apply_bound_eloc(const ComplexType e, const RealType eshift) const
  {
     // Leaving the imag part untouched, since it is not used. 
     // Only stored in case it might be useful in the future. 
     return ComplexType(std::max( std::min( e.real(), eshift+dEloc ), eshift-dEloc  ),e.imag()); 
  }

  inline ComplexType apply_bound_weight(const ComplexType w ) const
  {
     return (std::abs(w)>std::abs(max_weight))?max_weight:w;
  }

  inline void apply_bound_vbias(SPComplexType* vec, int n)
  {
     SPRealType mag=0.0;
     for(int i=0; i<n; i++,vec++) { 
       mag = std::abs(*vec);
       if(mag > vbias_bound) (*vec)/=(mag/vbias_bound);
     }
  }

  inline void apply_bound_vbias()
  {
     SPRealType mag=0.0;
     for(int i=0; i<vbias.size(); i++) { 
       mag = std::abs(vbias[i]);
       if(mag > vbias_bound) vbias[i]/=(mag/vbias_bound);
     }
  }

  inline SPComplexType apply_bound_vbias(SPComplexType v)
  {
    return (std::abs(v)>vbias_bound)?(v/(std::abs(v)/static_cast<SPValueType>(vbias_bound))):(v);
  }

const char* show_classification(double x) {
    switch(std::fpclassify(x)) {
        case FP_INFINITE:  return "Inf";
        case FP_NAN:       return "NaN";
        case FP_NORMAL:    return "normal";
        case FP_SUBNORMAL: return "subnormal";
        case FP_ZERO:      return "zero";
        default:           return "unknown";
    }
}

};


}

#endif

