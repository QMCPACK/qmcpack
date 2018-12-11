#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANFACTORY_HELPER_H
#define QMCPLUSPLUS_AFQMC_HAMILTONIANFACTORY_HELPER_H

#include "Configuration.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D.h"
#include "AFQMC/Hamiltonians/createHamiltonian_Helper.hpp"

namespace qmcplusplus
{

namespace afqmc
{

  void countElementsFromFCIDUMP(std::ifstream& in, const RealType cutoff1bar, const bool spinRestricted, const int NMO, const int NC, const int NAEA, const int NAEB, int& n1, int& n1_c, int& n2, int& n2_c, int& n2_m, int& n3, int& n3_c, int& n3_m, std::map<IndexType,IndexType>& orbMapA, std::map<IndexType,IndexType>& orbMapB, int& n3_vec );


 // This routine assumes that states are ordered (after mapping) in the following way:
 // { core, active+virtual, ignored}, so make sure mappings obey this.   
 // An easy way to implement this is to have a list of all core states and 
 // initialize the mapping (from 1:NC) with the index of the core states.
 // Then complete the mapping with every other orbital index not in the core list.
 // This way you always end up with a mapping with the correct format. 
  template<class shm_vector, class shm_SpMat>
  inline ValueType readElementsFromFCIDUMP(std::ifstream& in, const RealType cutoff1bar, 
         const bool spinRestricted, const int NMO, const int NC, const int NAEA, const int NAEB,
         std::vector<s2D<ValueType> >& V1,
//         std::vector<s2D<ValueType> >& V1_c,
         shm_vector* V2,
//         std::vector<s4D<ValueType> >& V2_c,
//         std::vector<s4D<ValueType> >& V2_m,
         shm_SpMat*  V3,
//         shm_SpMat&  V3_c,
//         shm_SpMat&  V3_m,
         std::map<IndexType,IndexType>& orbMapA, std::map<IndexType,IndexType>& orbMapB) {
     ValueType NuclearCoulombEnergy=0;
     IndexType a,b,c,d, cntS=0, cntD=0,q1;
     IndexType ap,bp,cp,dp,ab,cd;

     if(!spinRestricted)
       APP_ABORT(" Error: spinRestricted has been disabled. \n\n\n"); 

     ValueType val;
     while(!in.eof()) {
       in>>val >>ap >>bp >>cp >>dp;
       if(in.fail()) break;

       // to reduce problems from numerical roundoff 
       if(std::abs(imag(val)) < 1e-8) setImag(val,0.0);

       a=orbMapA[ap];
       b=orbMapA[bp];
       c=orbMapA[cp];
       d=orbMapA[dp];

       if(a>0 && b>0 && c > 0 && d>0) {  // 2-electron ME (ab|cd) = <ac|bd>
         if( std::abs(val) > cutoff1bar ) {
           V2->push_back(find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-1,c-1,b-1,d-1,val) ));
         }
       } else if(a > 0 && b > 0 && c>0) { // factorized 2-electron integral 
         if( std::abs(val) > cutoff1bar ) {
           V3->add((a-1)*NMO+Index2Col(NMO,b-1),cp-1,val);
         }
       } else if(a > 0 && b > 0) { // 1-electron ME (a|b) = <a|b>
         if( std::abs(val) > cutoff1bar ) {
           if(a > b) {q1=b;b=a;a=q1;val = myconj(val);}
           V1.emplace_back(std::make_tuple(a-1,b-1,val));
         }
       } else if(a==0 && b==0 && c==0 && d==0) {
         NuclearCoulombEnergy = val;
       } else if(a!=0 && b==0 && c==0 && d==0) {
         // ignore, these are the eigenvalues of the Fock operator printed by VASP  
       } else {
         app_error()<<"Error with ASCII integral file, bad indexing. " <<std::endl;
         APP_ABORT(" Error in readElementsFromFCIDUMP. \n");
       }
     }
     return NuclearCoulombEnergy;
  }

  inline HamiltonianTypes peekHamType(hdf_archive& dump) {
    if(dump.is_group( std::string("/Hamiltonian/THC") )) return THC;    
    if(dump.is_group( std::string("/Hamiltonian/Factorized") )) return Factorized;    
    if(dump.is_group( std::string("/Hamiltonian/SymmetricFactorized") )) return SymmetricFactorized;    
    if(dump.is_group( std::string("/Hamiltonian/Integrals") )) return s4DInts;    
    APP_ABORT("  Error: Invalid hdf file format in peekHamType(hdf_archive). \n");
    return s4DInts; 
  }

  inline HamiltonianTypes peekHamType(std::ifstream& in) {
    std::streampos start = in.tellg();
    IndexType ap,bp,cp,dp;
    ValueType val;    
    in>>val >>ap >>bp >>cp >>dp; 
    if(in.fail()) 
      APP_ABORT("Problems reading input file in peekHamType. \n");
    in.clear();
    in.seekg(start);
    if( ap != 0 && bp != 0 && cp != 0 && dp != 0 ) return s4DInts;
    else if( ap != 0 && bp != 0 && cp != 0 && dp == 0 ) return Factorized;
    APP_ABORT("  Error: Invalid fcidump file format in peekHamType(ifstream). \n");
    return s4DInts;
  }

  // old routines
  SPValueSMSpMat read_V2fact_old(hdf_archive& dump, TaskGroup_& TG, int nread, int NMO, int& mini, int& max, int nvecs, double cutoff1bar, int int_blocks); 

   FactorizedSparseHamiltonian::shm_csr_matrix read_V2fact(hdf_archive& dump, TaskGroup_& TG, int nread, int NMO, int nvecs, double cutoff1bar, int int_blocks); 
  SMDenseVector<s4D<ValueType> > read_V2(hdf_archive& dump, TaskGroup_& TG, int nread, int NMO, int& min_i, int& max_i, double cutoff1bar, int int_blocks);
  void read_THC(hdf_archive&); 

  void ascii_write(std::string); 
  void hdf_write(std::string); 

}

}

#endif
