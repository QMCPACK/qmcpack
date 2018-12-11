//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_SPARSETENSOR_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_SPARSETENSOR_HPP

#include <vector>
#include <type_traits>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/HamiltonianOperations/sparse_matrix_energy.hpp"

namespace qmcplusplus
{

namespace afqmc
{

// T1 depends on whether the integrals are half rotated or not in the real integrals case
// T2 depends on whether the transposed CholMatrix is half rotated or not in the real integrals case
template<class T1, class T2=T1>
class SparseTensor
{
  
#if defined(AFQMC_SP) 
  using SpT1 = typename to_single_precision<T1>::value_type; 
  using SpT2 = typename to_single_precision<T2>::value_type; 
#else
  using SpT1 = T1;
  using SpT2 = T2;
#endif
  using T1shm_csr_matrix = ma::sparse::csr_matrix<SpT1,int,std::size_t,
                                boost::mpi3::intranode::allocator<SpT1>,
                                boost::mpi3::intranode::is_root>;
  using T1shm_csr_matrix_view = typename T1shm_csr_matrix::template matrix_view<int>; 
  using T2shm_csr_matrix = ma::sparse::csr_matrix<SpT2,int,std::size_t,
                                boost::mpi3::intranode::allocator<SpT2>,
                                boost::mpi3::intranode::is_root>;
  using T2shm_csr_matrix_view = typename T2shm_csr_matrix::template matrix_view<int>; 
  using Vshm_csr_matrix = ma::sparse::csr_matrix<SPValueType,int,std::size_t,
                                boost::mpi3::intranode::allocator<SPValueType>,
                                boost::mpi3::intranode::is_root>;
  using Vshm_csr_matrix_view = typename Vshm_csr_matrix::template matrix_view<int>; 
  using CVector = boost::multi_array<ComplexType,1>;  
  using CMatrix = boost::multi_array<ComplexType,2>;  
  using SpT1Vector = boost::multi_array<SpT1,1>;  
  using T1Vector = boost::multi_array<T1,1>;  
  using T1Matrix = boost::multi_array<T1,2>;  
  using SpVector = boost::multi_array<SPComplexType,1>;  
  using this_t = SparseTensor<T1,T2>;

  public:

    SparseTensor(WALKER_TYPES type,
                 CMatrix&& hij_,
                 std::vector<SpT1Vector>&& h1, 
                 std::vector<T1shm_csr_matrix>&& v2, 
                 std::vector<T1shm_csr_matrix_view>&& v2view, 
                 Vshm_csr_matrix&& vn, 
                 Vshm_csr_matrix_view&& vnview, 
                 CMatrix&& vn0_, 
                 T2shm_csr_matrix&& vnT, 
                 T2shm_csr_matrix_view&& vnTview, 
                 ValueType e0_,
                 int gncv, 
                 bool separateeejab_=false ):
        walker_type(type),
        global_nCV(gncv),
        separateEJab(separateeejab_),
        E0(e0_),
        hij(std::move(hij_)),
        haj(std::move(h1)),
        Vakbl(std::move(v2)),
        Vakbl_view(std::move(v2view)),
        Spvn(std::move(vn)),
        Spvn_view(std::move(vnview)),
        SpvnT(std::move(vnT)),
        SpvnT_view(std::move(vnTview)),
        vn0(std::move(vn0_))
    {
        if(separateEJab && haj.size()>1)
          APP_ABORT("Error: separateEJab && haj.size()>1 not yet allowed. \n");
    }

    ~SparseTensor() {}
    
    SparseTensor(const SparseTensor& other) = delete;
    SparseTensor& operator=(const SparseTensor& other) = delete;
    SparseTensor(SparseTensor&& other) = default; 
    SparseTensor& operator=(SparseTensor&& other) = default;

    CMatrix getOneBodyPropagatorMatrix(TaskGroup_& TG, CVector const& vMF) {

      int NMO = hij.shape()[0];  
      // in non-collinear case with SO, keep SO matrix here and add it
      // for now, stay collinear
      CMatrix H1(extents[NMO][NMO]);

      // add sum_n vMF*Spvn, vMF has local contribution only!
      boost::multi_array_ref<ComplexType,1> H1D(H1.origin(),extents[NMO*NMO]);
      std::fill_n(H1D.origin(),H1D.num_elements(),ComplexType(0));
      vHS(vMF, H1D);  
      TG.TG().all_reduce_in_place_n(H1D.origin(),H1D.num_elements(),std::plus<>());

      // add hij + vn0 and symmetrize
      using std::conj;

      for(int i=0; i<NMO; i++) {
        H1[i][i] += hij[i][i] + vn0[i][i];
        for(int j=i+1; j<NMO; j++) {
          H1[i][j] += hij[i][j] + vn0[i][j];
          H1[j][i] += hij[j][i] + vn0[j][i];
          // This is really cutoff dependent!!!  
          if( std::abs( H1[i][j] - conj(H1[j][i]) ) > 1e-6 ) {
            app_error()<<" WARNING in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
            app_error()<<i <<" " <<j <<" " <<H1[i][j] <<" " <<H1[j][i] <<" " 
                       <<hij[i][j] <<" " <<hij[j][i] <<" "
                       <<vn0[i][j] <<" " <<vn0[j][i] <<std::endl;  
            //APP_ABORT("Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n"); 
          }
          H1[i][j] = 0.5*(H1[i][j]+conj(H1[j][i]));
          H1[j][i] = conj(H1[i][j]);
        }
      }

      return H1;  
    }

#if defined(AFQMC_SP) 
 template<class Mat, class MatB,
            typename = typename std::enable_if_t<
                    std::is_same<typename std::decay<MatB>::type::element, SPComplexType>::value> 
    >
    void energy(Mat&& E, const MatB& Gc, int k, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      APP_ABORT(" Error: AFQMC_SP not yet implemented in HamitonianOperations::SparseTensor(). \n"); 
    }
#endif
    
    template<class Mat, class MatB,
            typename = typename std::enable_if_t<
                    std::is_same<typename std::decay<MatB>::type::element, ComplexType>::value> 
    >
    void energy(Mat&& E, const MatB& Gc, int k, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      assert(E.shape()[1]>=3);
      assert(k >= 0 && k < haj.size());  
      assert(k >= 0 && k < Vakbl_view.size());  
      if(Gcloc.num_elements() < Gc.shape()[1] * Vakbl_view[k].shape()[0])
        Gcloc.resize(extents[Vakbl_view[k].shape()[0]*Gc.shape()[1]]);
      boost::multi_array_ref<SPComplexType,2> buff(Gcloc.data(),
                        extents[Vakbl_view[k].shape()[0]][Gc.shape()[1]]);
      shm::calculate_energy(std::forward<Mat>(E),Gc,buff,haj[k],Vakbl_view[k],addH1);
      // testing how to do this right now, make clean design later
      // when you write the FastMSD class 
      if(separateEJab) {
        using ma::T;
        if(haj.size()>1)
          APP_ABORT("Error: separateEJab && haj.size()>1 not yet allowed. \n");
        if(Gcloc.num_elements() < SpvnT.shape()[0] * Gc.shape()[1])
          Gcloc.resize(extents[SpvnT.shape()[0]*Gc.shape()[1]]);
        RealType scl = (walker_type==CLOSED?4.0:1.0); 
        // SpvnT*G
        boost::multi_array_ref<T2,2> v_(Gcloc.origin()+
                                            SpvnT_view.local_origin()[0]*Gc.shape()[1],
                                        extents[SpvnT_view.shape()[0]][Gc.shape()[1]]);
        ma::product(SpvnT_view, Gc, v_); 
        for(int wi=0; wi<Gc.shape()[1]; wi++)
          E[wi][2] = 0.5*scl*ma::dot(v_[indices[range_t()][wi]],v_[indices[range_t()][wi]]); 
      }
      if(addH1) 
        for(int i=0; i<E.shape()[0]; i++)
          E[i][0] += E0;  
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vHS(const MatA& X, MatB&& v, double a=1., double c=0.) {
      assert( Spvn.shape()[1] == X.shape()[0] );
      assert( Spvn.shape()[0] == v.shape()[0] );
      using Type = typename std::decay<MatB>::type::element;

      // Spvn*X 
      boost::multi_array_ref<Type,1> v_(v.origin() + Spvn_view.local_origin()[0], 
                                        extents[Spvn_view.shape()[0]]);
      ma::product(SPValueType(a),Spvn_view,X,SPValueType(c),v_);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(const MatA& X, MatB&& v, double a=1., double c=0.) {
      assert( Spvn.shape()[1] == X.shape()[0] );
      assert( Spvn.shape()[0] == v.shape()[0] );
      assert( X.shape()[1] == v.shape()[1] );
      using Type = typename std::decay<MatB>::type::element;

      // Spvn*X 
      boost::multi_array_ref<Type,2> v_(v[Spvn_view.local_origin()[0]].origin(), 
                                        extents[Spvn_view.shape()[0]][v.shape()[1]]);
      ma::product(SPValueType(a),Spvn_view,X,SPValueType(c),v_);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0.) {
      assert( SpvnT.shape()[1] == G.shape()[0] );
      assert( SpvnT.shape()[0] == v.shape()[0] );
      using Type = typename std::decay<MatB>::type::element ;

      // SpvnT*G
      boost::multi_array_ref<Type,1> v_(v.origin() + SpvnT_view.local_origin()[0], 
                                        extents[SpvnT_view.shape()[0]]);
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view, G, SpT2(c), v_);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0.) {
      assert( SpvnT.shape()[1] == G.shape()[0] );
      assert( SpvnT.shape()[0] == v.shape()[0] );   
      assert( G.shape()[1] == v.shape()[1] );
      using Type = typename std::decay<MatB>::type::element ;

      // SpvnT*G
      boost::multi_array_ref<Type,2> v_(v[SpvnT_view.local_origin()[0]].origin(), 
                                        extents[SpvnT_view.shape()[0]][v.shape()[1]]);
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view, G, SpT2(c), v_);
    }

    bool distribution_over_cholesky_vectors() const{ return true; }
    int local_number_of_cholesky_vectors() const{ return Spvn.shape()[1]; }
    int global_number_of_cholesky_vectors() const{return global_nCV; }

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const{return false;}
    bool transposed_G_for_E() const{return false;} 
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const{return false;} 

  private:

    WALKER_TYPES walker_type;

    int global_nCV;

    bool separateEJab;

    ValueType E0;

    // bare one body hamiltonian
    CMatrix hij;

    // (potentially half rotated) one body hamiltonian
    std::vector<SpT1Vector> haj;

    // sparse 2-body 2-electron integrals in matrix form  
    std::vector<T1shm_csr_matrix> Vakbl; 

    // sparse sub-matrix view
    std::vector<T1shm_csr_matrix_view> Vakbl_view;

    // Cholesky factorization of 2-electron integrals in sparse matrix form 
    Vshm_csr_matrix Spvn;

    // sparse sub-matrix view 
    Vshm_csr_matrix_view Spvn_view;

    // Cholesky factorization of 2-electron integrals in sparse matrix form 
    T2shm_csr_matrix SpvnT;

    // sparse sub-matrix view 
    T2shm_csr_matrix_view SpvnT_view;

    // one-body piece of Hamiltonian factorization
    CMatrix vn0;

    // local storage 
    SpVector Gcloc;

};

}

}

#endif
