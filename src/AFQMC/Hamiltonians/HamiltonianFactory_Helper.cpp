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
#include "AFQMC/Hamiltonians/createHamiltonian_Helper.hpp"

#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/SymmetricFactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D.h"
#include "AFQMC/Utilities/readHeader.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Matrix/hdf5_readers_new.hpp"
#include "AFQMC/Matrix/array_partition.hpp"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/coo_matrix.hpp"

namespace qmcplusplus
{

namespace afqmc
{

  // This routine operates on the FULL MO set, including CORE, ACTIVE and IGNORED states. 
  void countElementsFromFCIDUMP(std::ifstream& in, const RealType cutoff1bar, const bool spinRestricted, const int NMO, const int NC, const int NAEA, const int NAEB, int& n1, int& n1_c, int& n2, int& n2_c, int& n2_m, int& n3, int& n3_c, int& n3_m, std::map<IndexType,IndexType>& orbMapA, std::map<IndexType,IndexType>& orbMapB, int& n3_vec )
  {

     IndexType a,b,c,d;
     IndexType ap,bp,cp,dp;
     ValueType val;
     int uhf_block=0;
     n3_vec=n1=n2=n1_c=n2_c=n2_m=n3=n3_c=n3_m=0;

     while(!in.eof()) {
       in>>val >>ap >>bp >>cp >>dp;
       if(in.fail()) break;

       // to reduce problems from numerical roundoff 
       if(std::abs(imag(val)) < 1e-8) setImag(val,0.0);

       // FCIDUMP format encodes UHF with uhf_block. All indexes are smaller than NMAX/NMO 
       // THIS SHOULD BE DONE AFTER MAPPING!!!
       // BUT YOU CAN'T CREATE THE MAPPING WITHOUT READING ALL STATES, CORRECT???
       // HOW TO SOLVE THIS??? 
       // READ ALL THE FIRST TIME????
      
       // trying to fix problem in 3-index terms
       auto it = orbMapA.find(cp);
       if (it == orbMapA.end())
        orbMapA[cp] = cp; 
       it = orbMapB.find(cp);
       if (it == orbMapB.end())
        orbMapB[cp] = cp; 

       if(spinRestricted) {
        a=orbMapA[ap];
        b=orbMapA[bp];
        c=orbMapA[cp];
        d=orbMapA[dp];
       } else {
         if(uhf_block == 0) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapA[cp];
           d=orbMapA[dp];
         } else if(uhf_block==1) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==2) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==3) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=3.\n";
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n";
              APP_ABORT(" Error in countElementsFromFCIDUMP. \n");
           }
           c=d=0;
         } else if(uhf_block==4) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=4.\n";
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n";
              APP_ABORT(" Error in countElementsFromFCIDUMP. \n");
           }
           c=d=0;
         }
       }

       if(a>0 && b>0 && c>0 && d > 0) {  // 2-electron ME (ab|cd) = <ac|bd>
         if( std::abs(val) > cutoff1bar ) {
           if( !goodSpinSector(a-1,c-1,b-1,d-1,NMO) ) {
             app_error()<<" Problems in countElementsFromFCIDUMP. Inconsistent two body term in integral file: " <<a <<" " <<b <<" " <<c <<" " <<d <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<" " <<val <<std::endl;
             APP_ABORT(" Error in countElementsFromFCIDUMP. \n");
           }
           int tmp1,tmp2,tmp3,tmp4;
           switch(uhf_block) {
             case 0: 
             {  // AAAA
               tmp1 = tmp2= NC;
               break;
             }
             case 1: 
             { // BBBB
               tmp1 = tmp2 = NMO+NC;
               tmp3=0;
               break;
             }
             case 2: 
             { // AABB
               tmp1 = NC;
               tmp2 = NMO+NC;
               break;
             }
           };
           if( a<=tmp1 && b<=tmp1 && c<=tmp2 && d<=tmp2 ) n2_c++;
           else if( a>tmp1 && b>tmp1 && c>tmp2 && d>tmp2 ) n2++; 
           else n2_m++;
         }
       } else if(a>0 && b>0 && c>0) { // factorized 2-electron ME (i,k|n), where (i,k|j,l) = sum_n (i,k|n)*(l,j|n)*  
         if( std::abs(val) > cutoff1bar ) {
           if( !goodSpinSector(a-1,a-1,b-1,b-1,NMO) ) {
             app_error()<<" Problems in countElementsFromFCIDUMP. Inconsistent factorized two body term in integral file: " <<a <<" " <<b <<" " <<c <<" " <<d <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<" " <<val <<std::endl;
             APP_ABORT(" Error in countElementsFromFCIDUMP. \n");
           }
           int tmp1,tmp2,tmp3,tmp4;
           switch(uhf_block) {
             case 0:
             {  // AAAA
               tmp1 = tmp2= NC;
               break;
             }
             case 1:
             { // BBBB
               tmp1 = tmp2 = NMO+NC;
               tmp3=0;
               break;
             }
           };
           if( a<=tmp1 && b<=tmp1 ) n3_c++;
           else if( a>tmp1 && b>tmp1 ) n3++;
           else n3_m++;
           if(cp > n3_vec) n3_vec = cp;
         }
       } else if(a>0 && b>0) { // 1-electron ME (a|b) = <a|b>
         if( std::abs(val) > cutoff1bar ) {

           if(spinRestricted) {
             if( a <= NC && b <= NC ) n1_c++;
             else if( a>NC && b>NC  ) n1++; 
           } else {
             if(uhf_block == 3) {
               if( a <= NC && b <= NC ) n1_c++;
               else if( a>NC && b>NC  ) n1++; 
             } else {
               if( a <= NMO+NC && b <= NMO+NC ) n1_c++;
               else if( a>NMO+NC && b>NMO+NC  ) n1++; 
             }
           }
 
         }
       } else if(a==0 && b==0 && c==0 && d==0) {
         if( std::abs(val)==0 && !spinRestricted ) {
           uhf_block++;
         }
       } else if(a!=0 && b==0 && c==0 && d==0) {
       } else {
         app_error()<<"Error with ASCII integral file, bad indexing. " <<std::endl;
         APP_ABORT(" Error in countElementsFromFCIDUMP. \n");
       }
    }
  }

/*
 // This routine assumes that states are ordered (after mapping) in the following way:
 // { core, active+virtual, ignored}, so make sure mappings obey this.   
 // An easy way to implement this is to have a list of all core states and 
 // initialize the mapping (from 1:NC) with the index of the core states.
 // Then complete the mapping with every other orbital index not in the core list.
 // This way you always end up with a mapping with the correct format. 
  template<class shm_vector, class shm_SpMat>
  ValueType readElementsFromFCIDUMP(std::ifstream& in, const RealType cutoff1bar, 
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
     int uhf_block=0;
     while(!in.eof()) {
       in>>val >>ap >>bp >>cp >>dp;
       if(in.fail()) break;

       // to reduce problems from numerical roundoff 
       if(std::abs(imag(val)) < 1e-8) setImag(val,0.0);

//       if(spinRestricted) { 
        a=orbMapA[ap];
        b=orbMapA[bp];
        c=orbMapA[cp];
        d=orbMapA[dp];
// / *
       } else {
         if(uhf_block == 0) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapA[cp];
           d=orbMapA[dp];
         } else if(uhf_block==1) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==2) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==3) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=3.\n"; 
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n"; 
              APP_ABORT(" Error in readElementsFromFCIDUMP. \n");
           }
           c=d=0;
         } else if(uhf_block==4) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=4.\n"; 
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n"; 
              APP_ABORT(" Error in readElementsFromFCIDUMP. \n");
           }
           c=d=0;
         }
       }
// * /

       if(a>0 && b>0 && c > 0 && d>0) {  // 2-electron ME (ab|cd) = <ac|bd>
         if( std::abs(val) > cutoff1bar ) {

           V2->push_back(find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-1,c-1,b-1,d-1,val) )); 
// / *
           int tmp1,tmp2,tmp4=0,tmp3=0;
           switch(uhf_block) {
             case 0:
             {  // AAAA
               tmp1 = tmp2 = NC;
               tmp3 = tmp4 = NC;
               break;
             }
             case 1:
             { // BBBB
               tmp1 = tmp2 = NMO+NC;
               tmp3 = tmp4 = NC+NC; 
               break;
             }
             case 2:
             { // AABB
               tmp1 = NC;
               tmp3 = NC;
               tmp2 = NMO+NC;
               tmp4 = NC+NC;
               break;
             }
             default:
             {
               app_error()<<"Problems reading FCIDUMP: " 
               <<ap <<" " 
               <<bp <<" " 
               <<cp <<" " 
               <<dp <<" " 
               <<val <<std::endl;
               APP_ABORT(" Error in readElementsFromFCIDUMP. \n");
             }
           };
           if( a<=tmp1 && b<=tmp1 && c<=tmp2 && d<=tmp2 ) 
             V2_c.push_back(find_smaller_equivalent_OneBar_for_integral_list(std::make_tuple(a-1,c-1,b-1,d-1,val))); 
           else if( a>tmp1 && b>tmp1 && c>tmp2 && d>tmp2 ) { 
             V2.push_back(find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-tmp3-1,c-tmp4-1,b-tmp3-1,d-tmp4-1,val) )); 
           } else { 
             V2_m.push_back(find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-1,c-1,b-1,d-1,val) )); 
           }
// * /
         }
       } else if(a > 0 && b > 0 && c>0) { // factorized 2-electron integral 
         if( std::abs(val) > cutoff1bar ) {

         V3->add((a-1)*NMO+Index2Col(NMO,b-1),cp-1,val);
// / *
           int tmp1,tmp2,tmp4=0,tmp3=0;
           switch(uhf_block) {
             case 0:
             {  // AAAA
               tmp1 = tmp2 = NC;
               tmp3 = tmp4 = NC;
               break;
             }
             case 1:
             { // BBBB
               tmp1 = tmp2 = NMO+NC;
               tmp3 = tmp4 = NC+NC;
               break;
             }
             default:
             {
               app_error()<<"Problems reading FCIDUMP: "
               <<ap <<" "
               <<bp <<" "
               <<cp <<" "
               <<dp <<" "
               <<val <<std::endl;
               APP_ABORT(" Error in readElementsFromFCIDUMP. \n");
             }
           };

           if( a<=tmp1 && b<=tmp1 )
             V3_c.insert((a-1)*NMO+Index2Col(NMO,b-1),cp-1,val);
           else if( a>tmp1 && b>tmp1 )
             V3.insert((a-tmp3-1)*NMO+Index2Col(NMO,b-tmp3-1),cp-1,val);
           else
             V3_m.insert((a-1)*NMO+Index2Col(NMO,b-1),cp-1,val);
// * /
         }
       } else if(a > 0 && b > 0) { // 1-electron ME (a|b) = <a|b>
         if( std::abs(val) > cutoff1bar ) { 

           if(a > b) {q1=b;b=a;a=q1;val = myconj(val);}
    
            V1.emplace_back(std::make_tuple(a-1,b-1,val)); 
//           if(spinRestricted) {
//             if( a <= NC && b <= NC ) { 
//               V1_c.emplace_back(std::make_tuple(a-1,b-1,val)); 
//             } else if( a>NC && b>NC  ) {
//               V1.emplace_back(std::make_tuple(a-NC-1,b-NC-1,val)); 
//             }
//           } else {
//             if(uhf_block == 3) {
//               if( a <= NC && b <= NC ) { 
//                 V1_c.emplace_back(std::make_tuple(a-1,b-1,val)); 
//               } else if( a>NC && b>NC  ) { 
//                 V1.emplace_back(std::make_tuple(a-NC-1,b-NC-1,val)); 
//               }
//             } else {
//               if( a <= NMO+NC && b <= NMO+NC ) { 
//                 V1_c.emplace_back(std::make_tuple(a-1,b-1,val)); 
//               } else if( a>NMO+NC && b>NMO+NC  ) { 
//                 V1.emplace_back(std::make_tuple(a-2*NC-1,b-2*NC-1,val)); 
//               }
//             }
//           }
         }
       } else if(a==0 && b==0 && c==0 && d==0) {
         if( std::abs(val)==0 && !spinRestricted ) {
           APP_ABORT(" Error: Format error in FCIDUMP file. \n"); 
           uhf_block++;
         } else {
           NuclearCoulombEnergy = val;
         } 
       } else if(a!=0 && b==0 && c==0 && d==0) {
         // ignore, these are the eigenvalues of the Fock operator printed by VASP  
       } else {
         app_error()<<"Error with ASCII integral file, bad indexing. " <<std::endl;
         APP_ABORT(" Error in readElementsFromFCIDUMP. \n");
       }
     } 
     return NuclearCoulombEnergy;   
  }
*/

  SPValueSMSpMat read_V2fact_old(hdf_archive& dump, TaskGroup_& TG, int nread, int NMO, int& min_i, int& max_i, int nvecs, double cutoff1bar, int int_blocks)
  {

      min_i = 0;
      max_i = nvecs;

      int nrows = NMO*NMO;
      bool distribute_Ham = (TG.getNNodesPerTG() < TG.getTotalNodes());

      Timer.reset("Generic1");
      Timer.start("Generic1");

      simple_matrix_partition<TaskGroup_,IndexType,RealType> split(nrows,nvecs,cutoff1bar);
      std::vector<IndexType> counts;
      // count dimensions of sparse matrix
      count_entries_hdf5_SpMat(dump,split,int_blocks,false,counts,TG,true,nread);

      std::vector<IndexType> sets(TG.getNumberOfTGs()+1);
      if(TG.getCoreID() < nread)
        split.partition_over_TGs(TG,false,counts,sets);

      if(distribute_Ham) {
        if(TG.getGlobalRank()==0) {
          app_log()<<" Partitioning of (factorized) Hamiltonian Vectors: ";
          for(int i=0; i<=TG.getNumberOfTGs(); i++)
            app_log()<<sets[i] <<" ";
          app_log()<<std::endl;
          app_log()<<" Number of terms in each partitioning: ";
          for(int i=0; i<TG.getNumberOfTGs(); i++)
            app_log()<<accumulate(counts.begin()+sets[i],counts.begin()+sets[i+1],0) <<" ";
          app_log()<<std::endl;
        }

        TG.Node().broadcast(sets.begin(),sets.end());
        min_i = sets[TG.getTGNumber()];
        max_i = sets[TG.getTGNumber()+1];
      }

      Timer.stop("Generic1");
      app_log()<<" -- Time to count elements in hdf5 read: "
               <<Timer.average("Generic1") <<"\n";

      // resize Spvn
      long sz; 
      if( TG.getCoreID()==0 )
        sz = std::accumulate(counts.begin()+min_i,counts.begin()+max_i,long(0));
      TG.Node().broadcast_value(sz);

      SPValueSMSpMat V2_fact(nrows,nvecs,TG.getCoreID()==0,std::string("SparseGeneralHamiltonian_V2"),&TG.Node());
      V2_fact.allocate(sz);

      Timer.reset("Generic1");
      Timer.start("Generic1");

      // read Spvn    
      read_hdf5_SpMat(V2_fact,split,dump,int_blocks,TG,true,nread);

      Timer.stop("Generic1");
      app_log()<<" -- Time to read matrix from hdf5: " 
               <<Timer.average("Generic1") <<"\n";

      return V2_fact;
  
  }

  // it is possible to subdivide the rows within equivalent nodes and communicate at the end
  FactorizedSparseHamiltonian::shm_csr_matrix read_V2fact(hdf_archive& dump, TaskGroup_& TG, int nread, int NMO, int nvecs, double cutoff1bar, int int_blocks)
  {
      using counter =  qmcplusplus::afqmc::sparse_matrix_element_counter;
      using Alloc = boost::mpi3::intranode::allocator<ValueType>;
      using ucsr_matrix = ma::sparse::ucsr_matrix<ValueType,int,std::size_t,
                                boost::mpi3::intranode::allocator<ValueType>,
                                ma::sparse::is_root>;

      int min_i = 0;
      int max_i = nvecs;

      int nrows = NMO*NMO;
      bool distribute_Ham = (TG.getNNodesPerTG() < TG.getTotalNodes());
      std::vector<IndexType> row_counts(nrows);

      Timer.reset("Generic1");
      Timer.start("Generic1");

      // calculate column range that belong to this node
      if(distribute_Ham) {  
        // count number of non-zero elements along each column (Chol Vec)
        std::vector<IndexType> col_counts(nvecs);
        csr_hdf5::multiple_reader_global_count(dump,
                                   counter(false,nrows,nvecs,0,nrows,0,nvecs,cutoff1bar),
                                   col_counts,TG,nread);

        std::vector<IndexType> sets(TG.getNumberOfTGs()+1);
        simple_matrix_partition<TaskGroup_,IndexType,RealType> split(nrows,nvecs,cutoff1bar);
        if(TG.getCoreID() < nread)
          split.partition_over_TGs(TG,false,col_counts,sets);

        if(TG.getGlobalRank()==0) {
          app_log()<<" Partitioning of (factorized) Hamiltonian Vectors: ";
          for(int i=0; i<=TG.getNumberOfTGs(); i++)
            app_log()<<sets[i] <<" ";
          app_log()<<std::endl;
          app_log()<<" Number of terms in each partitioning: ";
          for(int i=0; i<TG.getNumberOfTGs(); i++)
            app_log()<<accumulate(col_counts.begin()+sets[i],col_counts.begin()+sets[i+1],0) <<" ";
          app_log()<<std::endl;
        }

        TG.Node().broadcast(sets.begin(),sets.end());
        min_i = sets[TG.getTGNumber()];
        max_i = sets[TG.getTGNumber()+1];

        csr_hdf5::multiple_reader_local_count(dump,
                                 counter(true,nrows,nvecs,0,nrows,min_i,max_i,cutoff1bar),
                                 row_counts,TG,nread);

      } else {  

        // should be faster if ham is not distributed
        csr_hdf5::multiple_reader_global_count(dump,
                                 counter(true,nrows,nvecs,0,nrows,0,nvecs,cutoff1bar),
                                 row_counts,TG,nread);

      }
 
      Timer.stop("Generic1");
      app_log()<<" -- Time to count elements in hdf5 read: "  
               <<Timer.average("Generic1") <<"\n";

      ucsr_matrix ucsr({nrows,max_i-min_i},{0,min_i},row_counts,Alloc(TG.Node()));
      csr::matrix_emplace_wrapper<ucsr_matrix> csr_wrapper(ucsr,TG.Node());

      Timer.reset("Generic1");
      Timer.start("Generic1");
    
      using mat_map =  qmcplusplus::afqmc::matrix_map;
      csr_hdf5::multiple_reader_hdf5_csr<ValueType,int>(csr_wrapper,
                                  mat_map(false,true,nrows,nvecs,0,nrows,min_i,max_i,cutoff1bar),  
                                  dump,TG,nread);
      csr_wrapper.push_buffer();
      TG.node_barrier();

      Timer.stop("Generic1");
      app_log()<<" -- Time to read into ucsr matrix: " 
               <<Timer.average("Generic1") <<"\n";

      // careful here!!!
      Timer.reset("Generic1");
      Timer.start("Generic1");

      return FactorizedSparseHamiltonian::shm_csr_matrix(std::move(ucsr));  
  }

  SMDenseVector<s4D<ValueType> > read_V2(hdf_archive& dump, TaskGroup_& TG, int nread, int NMO, int& min_i, int& max_i, double cutoff1bar, int int_blocks ) 
  {

      bool distribute_Ham = (TG.getNumberOfTGs()>1);
      int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
      int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
      std::vector<int> Idata;
      std::vector<ValueType> vvec;

      std::vector<long> ntpo(NMO,0);
      std::vector<IndexType> indxvec;
      int ntmax;
      std::vector<int> pool_dist;
      if(coreid < nread) {

        // Idata[i]: number of terms per block
        Idata.resize(int_blocks);
        if(!dump.read(Idata,"block_sizes")) {
          app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_block_sizes dataset. \n";
          APP_ABORT(" ");
        }

        ntmax = *std::max_element(Idata.begin(), Idata.end());
        indxvec.reserve(2*ntmax);
        vvec.reserve(ntmax);

        // divide blocks 
        FairDivide(int_blocks,nread,pool_dist);

        int first_local_block = pool_dist[coreid];
        int last_local_block = pool_dist[coreid+1];

        for(int ib=first_local_block; ib<last_local_block; ib++) {
          if( ib%nnodes != nodeid ) continue;
          if(Idata[ib]==0) continue;
          indxvec.resize(2*Idata[ib]);
          vvec.resize(Idata[ib]);
          if(!dump.read(indxvec,std::string("index_")+std::to_string(ib))) {
            app_error()<<" Error in read_V2(): Problems reading index_" <<ib <<" dataset. \n";
            APP_ABORT(" ");
          }
          if(!dump.read(vvec,std::string("vals_")+std::to_string(ib))) {
            app_error()<<" Error in read_V2(): Problems reading vals_" <<ib <<" dataset. \n";
            APP_ABORT(" ");
          }

          std::vector<IndexType>::iterator iti = indxvec.begin();
          for(std::vector<ValueType>::iterator itv = vvec.begin(); itv < vvec.end(); itv++, iti+=2) {
            if(std::abs(*itv) < cutoff1bar )
                continue;
            s4D<ValueType> ijkl = std::make_tuple(  static_cast<OrbitalType>((*iti)/NMO),
                                                    static_cast<OrbitalType>((*(iti+1))/NMO),
                                                    static_cast<OrbitalType>((*iti)%NMO),
                                                    static_cast<OrbitalType>((*(iti+1))%NMO), ValueType(0));
            find_smallest_permutation(ijkl);
            ntpo[std::get<0>(ijkl)]++;
          }
        }
      }
      {  
        std::vector<long> dum(ntpo);
        TG.Global().all_reduce(dum.begin(),dum.end(),ntpo.begin());
      }
  
      if(distribute_Ham) {
        std::vector<long> sets(TG.getNumberOfTGs()+1);
        if(TG.getGlobalRank() == 0) {
          std::vector<long> nv(NMO+1);
          nv[0]=0;
          for(int i=0; i<NMO; i++)
            nv[i+1]=nv[i]+ntpo[i];
          balance_partition_ordered_set(NMO,nv.data(),sets);
          app_log()<<" Hamiltonian partitioning: \n   Orbitals:        ";
          for(int i=0; i<=TG.getNumberOfTGs(); i++) app_log()<<sets[i] <<" ";
          app_log()<<std::endl <<"   Terms per block: ";
          for(int i=0; i<TG.getNumberOfTGs(); i++) app_log()<<nv[sets[i+1]]-nv[sets[i]] <<" ";
          app_log()<<std::endl;
        }
        TG.Global().broadcast(sets.begin(),sets.end());
        min_i = static_cast<int>(sets[TG.getTGNumber()]);
        max_i = static_cast<int>(sets[TG.getTGNumber()+1]);
      }

      long nttot = std::accumulate(ntpo.begin()+min_i,ntpo.begin()+max_i,long(0));
      SMDenseVector<s4D<ValueType> > V2(TG.getCoreID()==0,std::string("SparseGeneralHamiltonian_V2"),&TG.Node()); 
      V2.reserve(nttot);

      if( coreid < nread ) {

        std::vector<IndexType> ivec2;
        std::vector<ValueType> vvec2;
        ivec2.reserve(2*ntmax);
        vvec2.reserve(ntmax);
        int maxv = 10000;
        std::vector<s4D<ValueType>> vals;
        vals.reserve(maxv);

        int first_local_block = pool_dist[coreid];
        int last_local_block = pool_dist[coreid+1];
        int nbtot = last_local_block-first_local_block;
        int niter = nbtot/nnodes + std::min(nbtot%nnodes,1);

        for(int iter=0; iter<niter; iter++) {
          int first_block = first_local_block + nnodes*iter;
          int last_block = std::min(first_block+nnodes,last_local_block);
          int myblock_number = first_block + nodeid;
          if(myblock_number < last_local_block && Idata[myblock_number] > 0) {
            indxvec.resize(2*Idata[myblock_number]);
            vvec.resize(Idata[myblock_number]);
            if(!dump.read(indxvec,std::string("index_")+std::to_string(myblock_number))) {
              app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading index_" <<myblock_number <<" dataset. \n";
              APP_ABORT(" ");
            }
            if(!dump.read(vvec,std::string("vals_")+std::to_string(myblock_number))) {
              app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading vals_" <<myblock_number <<" dataset. \n";
              APP_ABORT(" ");
            }
          }


          for(int k=first_block,ipr=0; k<last_block; k++,ipr++) {
            if(Idata[k] == 0) continue;
            ivec2.resize(2*Idata[k]);
            vvec2.resize(Idata[k]);
            if(ipr==nodeid) {
              assert(myblock_number==k);
              std::copy(indxvec.begin(),indxvec.end(),ivec2.begin());
              std::copy(vvec.begin(),vvec.end(),vvec2.begin());
            }

            TG.Cores().broadcast(ivec2.begin(),ivec2.end(),ipr);
            TG.Cores().broadcast(vvec2.begin(),vvec2.end(),ipr);

            std::vector<IndexType>::iterator iti = ivec2.begin();
            for(std::vector<ValueType>::iterator itv = vvec2.begin(); itv < vvec2.end(); itv++, iti+=2) {
              if(std::abs(*itv) < cutoff1bar )
                  continue;
              s4D<ValueType> ijkl = std::make_tuple(  static_cast<OrbitalType>((*iti)/NMO),
                                                      static_cast<OrbitalType>((*(iti+1))/NMO),
                                                      static_cast<OrbitalType>((*iti)%NMO),
                                                      static_cast<OrbitalType>((*(iti+1))%NMO), *itv);
              find_smallest_permutation(ijkl);
              if( std::get<0>(ijkl) >= min_i && std::get<0>(ijkl) < max_i) {
                vals.push_back(ijkl);
                if(vals.size()==maxv) {
                  V2.push_back(vals,nread>1);
                  vals.clear();
                }
              }
            }
          }
        }
        if(vals.size() > 0)
          V2.push_back(vals,nread>1);
      }

      Timer.reset("Generic2");
      Timer.start("Generic2");
      V2.sort (
        [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs)
          {
            return std::forward_as_tuple(std::get<0>(lhs),std::get<1>(lhs),std::get<2>(lhs),std::get<3>(lhs)) < std::forward_as_tuple(std::get<0>(rhs),std::get<1>(rhs),std::get<2>(rhs),std::get<3>(rhs));
          },
          &TG.Node());
      Timer.stop("Generic2");
      app_log()<<" -- Time to compress Hamiltonian from h5 file: " <<Timer.average("Generic2") <<"\n";
    return V2;

  }

  void ascii_write(std::string ascii_write_file) {

/*
    if(ascii_write_file == std::string("")) return;

    if(factorizedHamiltonian) {
      app_log()<<" factorizedHamiltonian output not yet implemented " <<std::endl;
      return;
    }
 
    std::ofstream out(ascii_write_file.c_str());
    if(out.fail()) {
      app_error()<<" Error opening restart file in SparseGeneralHamiltonian::ascii_write() . " <<std::endl; 
      return;
    }

// hack to write sparse matrix in ascii for testing
/ *
   if(!factorizedHamiltonian) 
     APP_ABORT("Error: Matrix dump only for factorized hamiltonian \n\n\n");
   if(rank()!=0) return;
   int nvecs=V2_fact.rows();
   out<<"# nrows: " <<nvecs <<" ncols: " <<V2_fact.cols() <<"\n";
   for(ComplexSMSpMat::intType i=0; i<nvecs; i++) {
     int k1 = *(V2_fact.rowIndex_begin()+i+1); 
     for(int k=*(V2_fact.rowIndex_begin()+i); k<k1; k++) 
       out<<i <<" " <<*(V2_fact.cols_begin()+k) <<" "
          <<std::setprecision(12) <<*(V2_fact.vals_begin()+k) <<"\n";
   }
   out<<std::endl;
   out.close();
   APP_ABORT(" FINISHED DUMPING ASCII FILE. \n\n\n");   
* /
 
    out<<"&FCI NORB=" <<NMO <<",NELEC=" <<NAEA+NAEB <<",MS2=" <<MS2 <<"," <<"\n"
       <<"ISYM=" <<ISYM <<"," <<"\n/\n"; 

    out.setf(std::ios::scientific, std::ios::floatfield);
    if(factorizedHamiltonian) {
      out<<" factorizedHamiltonian output not yet implemented " <<std::endl;
      return;
    } else {

      SMDenseVector<s4D<ValueType> >::iterator it=V2.begin();
      for( ; it!=V2.end(); it++) {
        if(std::get<0>(*it) >= NMO) break;
        out<<std::setprecision(16) <<std::get<4>(*it) <<std::setw(5) <<std::get<0>(*it)+1 <<std::setw(5) <<std::get<2>(*it)+1 
            <<std::setw(5) <<std::get<1>(*it)+1 <<std::setw(5) <<std::get<3>(*it)+1 <<"\n";
      }
      if(!spinRestricted) {
        out<<" UHF output not yet implemented " <<std::endl;
      }

    }

    s2Dit it=H1.begin();
    for( ; it!=H1.end(); it++) {
      if(std::get<0>(*it) >= NMO) break; 
      out<<std::setprecision(16) <<std::get<2>(*it) <<std::setw(5) <<std::get<0>(*it)+1 <<std::setw(5) <<std::get<1>(*it)+1 <<"    0    0\n"; 
    } 
    if(!spinRestricted) {
      out<<" 0.0                    0    0    0    0\n";
      for( ; it!=H1.end(); it++) {
        if(std::get<0>(*it) >= NMO) break;
        out<<std::setprecision(16) <<std::get<2>(*it) <<" " <<std::get<0>(*it)-NMO+1 <<" " <<std::get<1>(*it)-NMO+1 <<"   0   0\n";
      }
    }
    out<<std::setprecision(16) <<NuclearCoulombEnergy <<"    0    0    0    0\n";

    out.close();
*/    
    
  } 

  // this needs to be modified to wread/write in blocks
  // to avoid having to allocate the full array twice
  // Also, you could add support in hdf_archive for tuples 
  void hdf_write(std::string hdf_write_file) {

/*
    if(skip_V2) 
      return; 

    if(hdf_write_file == std::string("")) return;

    bool writeFactorized = factorizedHamiltonian;
    if(hdf_write_type!="default") {
      if(factorizedHamiltonian) writeFactorized = (hdf_write_type == "factorized");
      else writeFactorized = (hdf_write_type != "integrals");
    }

    // pseudo-hack to avoid a lot of indenting 
    if(myComm->rank() > 0) {


      // only parallelized path
      if(!writeFactorized && factorizedHamiltonian) {
        std::vector<OrbitalType> ivec; 
        std::vector<ValueType> vvec;
        for(int k=0; k<NMO; k++) {
          ivec.clear();
          vvec.clear();
          SparseHamiltonianFromFactorization(k,ivec,vvec,cutoff1bar);
        }
      }
  
      return; 
    }

    hdf_archive dump(myComm);
    if(!dump.create(hdf_write_file)) {
      app_error()<<" Error opening restart file. \n";
      return;
    }

    std::string path = "/Hamiltonian/SparseGeneralHamiltonian";
    if(dump.is_group( path )) {
      app_error()<<" ERROR: H5Group /Hamiltonian/SparseGeneralHamiltonian already exists in restart file. Not over-writing data in file. \n"; 
      return;
    }

    dump.push("Hamiltonian");
    dump.push("SparseGeneralHamiltonian");

    int V2tot=0;
    std::vector<int> Idata(8);
/ * writing at the end since I don't have all the information necessarily
    Idata[0]=H1.size();
    if(factorizedHamiltonian)
      Idata[1]=V2_fact.size();
    else
      Idata[1]=V2.size();
    Idata[2]=0;
    Idata[3]=NMO;
    Idata[4]=NAEA;
    Idata[5]=NAEB;
    Idata[6]=spinRestricted?(0):(1);
    Idata[7]=0;
    if(factorizedHamiltonian)
      Idata[7] = V2_fact.cols(); 
    dump.write(Idata,"dims");
* /

    Idata.resize(NAEA+NAEB);
    for(int i=0; i<NAEA; i++) Idata[i] = occup_alpha[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) Idata[i] = occup_beta[j];
    dump.write(Idata,"occups");

    std::vector<ValueType> Rdata(2);
    Rdata[0] = NuclearCoulombEnergy;
    Rdata[1] = FrozenCoreEnergy; 
    dump.write(Rdata,"Energies");

    // write H1 
    std::vector<OrbitalType> ivec; 
    ivec.resize(2*H1.size());
    for(int i=0, j=0; i<H1.size(); i++, j+=2) 
      std::tie (ivec[j],ivec[j+1],std::ignore) = H1[i];   
    dump.write(ivec,"H1_indx");

    std::vector<ValueType> vvec;
    vvec.resize(H1.size());
    for(int i=0; i<H1.size(); i++)
      std::tie (std::ignore,std::ignore,vvec[i]) = H1[i];
    dump.write(vvec,"H1");

    if(writeFactorized) {

      if(!factorizedHamiltonian) {
        app_error()<<" Error: Can only write factorized hamiltonian if input hamiltonian is already factorized. \n\n\n";
        APP_ABORT(" Error: Can only write factorized hamiltonian if input hamiltonian is already factorized. \n\n\n");
      }

      if(cholesky_residuals.size() != V2_fact.cols()) {
        app_error()<<"Error: Incorrect size of cholesky_residuals: " <<cholesky_residuals.size() <<" " <<V2_fact.cols() <<std::endl;  
        APP_ABORT("Error. \n\n\n");
      }
      dump.write(cholesky_residuals,std::string("V2fact_vec_residual"));

      V2_fact.transpose();

      int nvecs=V2_fact.rows();
      std::vector<int> sz(nvecs);
      for(int i=0; i<nvecs; i++) 
        sz[i]=*(V2_fact.rowIndex_begin()+i+1) - *(V2_fact.rowIndex_begin()+i);
      int nmax = *std::max_element(sz.begin(),sz.end());
      dump.write(sz,std::string("V2fact_vec_sizes"));
      std::vector<IndexType> ivec;
      std::vector<ValueType> vvec;
      ivec.reserve(nmax);
      vvec.reserve(nmax);
      for(ComplexSMSpMat::intType i=0; i<nvecs; i++) {

        ivec.resize(sz[i]);
        vvec.resize(sz[i]);
        std::copy(V2_fact.cols_begin()+(*(V2_fact.rowIndex_begin()+i)),V2_fact.cols_begin()+(*(V2_fact.rowIndex_begin()+i+1)),ivec.begin());
        std::copy(V2_fact.vals_begin()+(*(V2_fact.rowIndex_begin()+i)),V2_fact.vals_begin()+(*(V2_fact.rowIndex_begin()+i+1)),vvec.begin());

        dump.write(ivec,std::string("V2fact_index_")+std::to_string(i));
        dump.write(vvec,std::string("V2fact_vals_")+std::to_string(i));

      }

      V2_fact.transpose();

    } else { 
      // write V2
      // terms with i index are found between [Idata[i-1],Idata[i])  
      Idata.resize(NMO);
      int ntmax=0;
      OrbitalType i0,j0,k0,l0,iold=0;

      if(factorizedHamiltonian) {

        for(int k=0; k<NMO; k++) {
          ivec.clear();
          vvec.clear();
          SparseHamiltonianFromFactorization(k,ivec,vvec,cutoff1bar); 
          Idata[k]= (k==0?0:Idata[k-1]+vvec.size());
          V2tot+=vvec.size();
          dump.write(ivec,std::string("V2_indx_")+std::to_string(k));
          dump.write(vvec,std::string("V2_")+std::to_string(k));
        }

        dump.write(Idata,"V2_block_sizes");

      } else {
        for(int i=0; i<V2.size(); i++) {
          std::tie (i0,j0,k0,l0,std::ignore) = V2[i];
          if(iold != i0) {
            for(int k=iold; k<i0; k++) Idata[k]=i; 
            iold=i0;
          }
        } 
        for(int k=iold; k<NMO; k++)  Idata[k]=V2.size();
        ntmax = Idata[0];
        for(int i=1; i<NMO; i++) {
          if(ntmax < (Idata[i]-Idata[i-1]))
            ntmax = Idata[i]-Idata[i-1];
        }

        dump.write(Idata,"V2_block_sizes");

        ivec.reserve(3*ntmax);
        vvec.reserve(ntmax);

        for(int k=0; k<NMO; k++) { 
          int ik0 = (k==0)?0:Idata[k-1]; 
          if(Idata[k]==ik0) continue; 
          ivec.clear();
          // no point in storing i
          ivec.resize(3*(Idata[k]-ik0));
          for (int i=ik0,jk=0; i<Idata[k]; i++,jk+=3) 
            std::tie (i0,ivec[jk],ivec[jk+1],ivec[jk+2],std::ignore) = V2[i];   
          dump.write(ivec,std::string("V2_indx_")+std::to_string(k));

          vvec.clear();
          vvec.resize(Idata[k]-ik0);
          for(int i=ik0; i<Idata[k]; i++) 
            std::tie (std::ignore,std::ignore,std::ignore,std::ignore,vvec[i-ik0]) = V2[i];
          V2tot+=vvec.size();
          dump.write(vvec,std::string("V2_")+std::to_string(k));
        }
 
      }

    }

    Idata.resize(8);
    Idata[0]=H1.size();
    Idata[1]=V2tot;
    Idata[2]=0;
    Idata[3]=NMO;
    Idata[4]=NAEA;
    Idata[5]=NAEB;
    Idata[6]=spinRestricted?(0):(1);
    Idata[7]=0;
    if(writeFactorized && factorizedHamiltonian)
      Idata[7] = V2_fact.cols();
    dump.write(Idata,"dims");

    dump.pop();
    dump.pop();

    dump.flush();
    dump.close();
*/

  }



}

}
