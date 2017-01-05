
#ifndef QMCPLUSPLUS_AFQMC_PHASELESS_WITHIMPSAMPLWITHELOC_FORCEBIAS
#define QMCPLUSPLUS_AFQMC_PHASELESS_WITHIMPSAMPLWITHELOC_FORCEBIAS

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
       
  phaseless_ImpSamp_ForceBias(Communicate *c,  RandomGenerator_t* r) : PropagatorBase(c,r), substractMF(true),use_eig(false),first(true),max_weight(100),apply_constrain(true),save_memory(false),vbias_bound(3.0),imp_sampl(true),hybrid_method(false),test_library(false),eloc_from_Spvn(false),sizeOfG(0)
  {
  } 

  ~phaseless_ImpSamp_ForceBias() {}

//  void Propagate(int n, SlaterDetWalker&, RealType& E1, const RealType E2=0);

  void Propagate(int n, WalkerHandlerBase*, RealType& E1, const RealType E2=0);

  bool setup(std::vector<int>&,ComplexSMVector*,HamiltonianBase*,WavefunctionHandler*,RealType, hdf_archive&, const std::string&,MPI_Comm,MPI_Comm);

  bool parse(xmlNodePtr);   

  bool hdf_write(hdf_archive&, const std::string&);

  bool hdf_read(hdf_archive&, const std::string&);

  void benchmark(){
    //PureSingleDeterminant* sd = dynamic_cast<PureSingleDeterminant*>(wfn->ImpSampWfn);
    //compare_libraries(NMO,NAEA,NAEB,Propg_H1,Propg_H1_indx,sd->Vijkl,vn,vn_indx);
  }

  private:

  std::ofstream out_debug;

  bool hybrid_method;

  bool eloc_from_Spvn;

  bool imp_sampl;

  bool substractMF;

  bool use_eig;

  bool first;

  bool save_memory;

  int test_library;

  std::ifstream in_rand;

  RealType cutoff;

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
  ComplexMatrix vHS;   
  ComplexMatrix PHS;

  // potentials that represent the decomposition of Vijkl into a quadratic form
  // Vijkl = -0.5 * sum_n (v2_n)^2 + v0, where v0 is stored somewhere else
  // stored in sparse form and sequentially
  // vn_nterms is a vector that contains the number of non-zero terms (above cutoff)
  // in each term in vn. This allows for easier parsing. 
  std::vector<s2D<ComplexType> > vn;  
  std::vector<IndexType> vn_indx;
  Vector< Vector<s2D<ValueType> >::iterator > vn_bounds;

//  // new storage format for HS potentials
//  ComplexSpMat Spvn; 
//  ComplexSpMat SpvnT; 
//  // this is a pointer to either Spvn or SpvnT, to avoid extra logic in code  
//  ComplexSpMat *Spvn_for_onebody; 

  int GlobalSpvnSize; 
  std:: vector<int> nCholVec_per_node;

  ComplexSMSpMat Spvn; 
  ComplexSMSpMat SpvnT; 
  // this is a pointer to either Spvn or SpvnT, to avoid extra logic in code  
  ComplexSMSpMat *Spvn_for_onebody; 

  // storage for fields
  std::vector<RealType> sigma;
  std::vector<ComplexType> CV0;

  // Force bias potential, typically mixed potential but can also be MF potential
  std::vector<ComplexType> vbias;   

  // Mean-field subtraction of vHS
  std::vector<ComplexType> vMF;   

  // used to calculate mean fields, overlaps, local energies, etc  
  WavefunctionHandler* wfn; 

  // local storage
  ComplexMatrix T1; 
  ComplexMatrix T2; 

  ComplexMatrix S1; 
  ComplexMatrix S2; 

  ComplexSMVector local_buffer;

  std::vector<ComplexType> MFfactor;
  std::vector<ComplexType> hybrid_weight;

  RealType dEloc;
  ComplexType max_weight;
  RealType vbias_bound;

  int sizeOfG;
  int nCholVecs;  // total number of Cholesky Vectors in the calculation  
  int cvec0,cvecN; // index of first and last+1 Cholesky Vector of this core 
  std::vector<int> walker_per_node;
  
  // ik breakup of Spvn
  IndexType ik0, ikN;   //  minimum and maximum values of ik index in Spvn
  IndexType pik0, pikN;  // locations of bounds of [ik0,ikN] sector in Spvn  

  void dist_Propagate(WalkerHandlerBase*);

  bool apply_constrain;

  void applyHSPropagator(ComplexMatrix&, ComplexMatrix&, ComplexType& factor, int order=-1, bool calculatevHS=true); 

  void addvHS(ComplexSMVector *buff, int nw, int sz, WalkerHandlerBase* wset); 

  void sampleGaussianFields();

  void sampleGaussianFields(ComplexType*,int);

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

  inline void apply_bound_vbias(ComplexType* vec, int n)
  {
     RealType mag=0.0;
     for(int i=0; i<n; i++,vec++) { 
       mag = std::abs(*vec);
       if(mag > vbias_bound) (*vec)/=(mag/vbias_bound);
     }
  }

  inline void apply_bound_vbias()
  {
     RealType mag=0.0;
     for(int i=0; i<vbias.size(); i++) { 
       mag = std::abs(vbias[i]);
       if(mag > vbias_bound) vbias[i]/=(mag/vbias_bound);
     }
  }

  inline ComplexType apply_bound_vbias(ComplexType v)
  {
    return (std::abs(v)>vbias_bound)?(v/(std::abs(v)/vbias_bound)):(v);
  }

  void print_tuple(std::vector<s2D<ComplexType> >& v) {
    for(int i=0; i<v.size(); i++) 
      std::cout<<"  -  " <<std::get<0>(v[i]) <<" " <<std::get<1>(v[i]) <<" " <<std::get<2>(v[i]) <<std::endl;
  }

  void print_octave(std::ofstream& out, ComplexMatrix& M) {

    int nC = M.cols();
    
    for(int i=0; i<NMO; i++) {
      for(int j=0; j<nC; j++)
        out<<"complex(" <<M(i,j).real() <<"," <<M(i,j).imag() <<")  ";
      out<<std::endl;
    }

  } 

  void test_linear_algebra();
};


}

#endif

