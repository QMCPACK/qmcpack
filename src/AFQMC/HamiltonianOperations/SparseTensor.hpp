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
#include "mpi3/shared_communicator.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/HamiltonianOperations/sparse_matrix_energy.hpp"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{

namespace afqmc
{

// T1 depends on whether the integrals are half rotated or not in the real integrals case
// T2 depends on whether the transposed CholMatrix is half rotated or not in the real integrals case
template<class T1, class T2=T1>
class SparseTensor
{

#if defined(MIXED_PRECISION)
  using SpT1 = typename to_single_precision<T1>::value_type;
  using SpT2 = typename to_single_precision<T2>::value_type;
#else
  using SpT1 = T1;
  using SpT2 = T2;
#endif
  using T1shm_csr_matrix = ma::sparse::csr_matrix<SpT1,int,std::size_t,
                                shared_allocator<SpT1>,
                                ma::sparse::is_root>;
  using T1shm_csr_matrix_view = typename T1shm_csr_matrix::template matrix_view<int>;
  using T2shm_csr_matrix = ma::sparse::csr_matrix<SpT2,int,std::size_t,
                                shared_allocator<SpT2>,
                                ma::sparse::is_root>;
  using T2shm_csr_matrix_view = typename T2shm_csr_matrix::template matrix_view<int>;
  using Vshm_csr_matrix = ma::sparse::csr_matrix<SPValueType,int,std::size_t,
                                shared_allocator<SPValueType>,
                                ma::sparse::is_root>;
  using Vshm_csr_matrix_view = typename Vshm_csr_matrix::template matrix_view<int>;
  using CVector = boost::multi::array<ComplexType,1>;
  using CMatrix = boost::multi::array<ComplexType,2>;
  using T1Vector = boost::multi::array<T1,1>;
  using T1Matrix = boost::multi::array<T1,2>;
  using SpVector = boost::multi::array<SPComplexType,1>;
  using shmSpVector = boost::multi::array<SPComplexType,1,shared_allocator<SPComplexType>>;
  using this_t = SparseTensor<T1,T2>;
  using communicator = boost::mpi3::shared_communicator;

  public:

    SparseTensor(communicator& c_, 
                 WALKER_TYPES type,
                 CMatrix&& hij_,
                 std::vector<CVector>&& h1,
                 std::vector<T1shm_csr_matrix>&& v2,
                 std::vector<T1shm_csr_matrix_view>&& v2view,
                 Vshm_csr_matrix&& vn,
                 Vshm_csr_matrix_view&& vnview,
                 CMatrix&& vn0_,
                 std::vector<T2shm_csr_matrix>&& vnT,
                 std::vector<T2shm_csr_matrix_view>&& vnTview,
                 ValueType e0_,
                 int gncv):
/*  2 defined behaviors:
 *  1. NOMSD expected behavior where a single vnT/vnTview is given and it must be consistent
 *     with a full G: NMO*NMO. In this case, the k index in vbias is ignored.
 *     In this case, only EXX is calculated and assumed to also contain EJ.
 *  2. PHMSD expected behavior, where h1.size() == v2.size() == vnT.size(),
 *     vbias must be calculated for each k independently.
 *     v2 is assumed to only contain EXX and EJ is calculated separately.
 *  In summary, vnT.size() determines the behavior of this routine.
 *  NOMSD wavefunctions provide half rotated vnT is single determinant and just transposed vn
 *  if multi-determinant.
 *  PHMSD provides different references (alpha/beta or multi-reference PH) in separate locations
 *  in the std:vector's.
 */
        comm(std::addressof(c_)),
        walker_type(type),
        global_nCV(gncv),
        E0(e0_),
        hij(std::move(hij_)),
        haj(std::move(h1)),
        Vakbl(std::move(v2)),
        Vakbl_view(std::move(v2view)),
        Spvn(std::move(vn)),
        Spvn_view(std::move(vnview)),
        SpvnT(std::move(vnT)),
        SpvnT_view(std::move(vnTview)),
        vn0(std::move(vn0_)),
        SM_TMats(iextensions<1u>{0},shared_allocator<SPComplexType>{c_}),
	separateEJ(true)
    {
	assert(haj.size() == Vakbl.size());
	assert(haj.size() == Vakbl_view.size());
	assert(SpvnT.size() == SpvnT_view.size());
	assert((haj.size() == SpvnT.size()) || (SpvnT.size()==1));
	assert((haj.size() == SpvnT_view.size()) || (SpvnT_view.size()==1));
	if((haj.size() > 1) && (SpvnT.size()==1)) // NOMSD with more than 1 determinant
          separateEJ = false;
    }

    ~SparseTensor() {}

    SparseTensor(const SparseTensor& other) = delete;
    SparseTensor& operator=(const SparseTensor& other) = delete;
    SparseTensor(SparseTensor&& other) = default;
    SparseTensor& operator=(SparseTensor&& other) = default;

    CMatrix getOneBodyPropagatorMatrix(TaskGroup_& TG, CVector const& vMF) {

      int NMO = hij.size(0);
      // in non-collinear case with SO, keep SO matrix here and add it
      // for now, stay collinear
      CMatrix H1({NMO,NMO});

      // add sum_n vMF*Spvn, vMF has local contribution only!
      boost::multi::array_ref<ComplexType,1> H1D(H1.origin(),{NMO*NMO});
      std::fill_n(H1D.origin(),H1D.num_elements(),ComplexType(0));
      vHS(vMF, H1D);
      TG.TG().all_reduce_in_place_n(H1D.origin(),H1D.num_elements(),std::plus<>());

      // add hij + vn0 and symmetrize
      using ma::conj;

      for(int i=0; i<NMO; i++) {
        H1[i][i] += hij[i][i] + vn0[i][i];
        for(int j=i+1; j<NMO; j++) {
          H1[i][j] += hij[i][j] + vn0[i][j];
          H1[j][i] += hij[j][i] + vn0[j][i];
          // This is really cutoff dependent!!!
          if( std::abs( H1[i][j] - ma::conj(H1[j][i]) ) > 1e-6 ) {
            app_error()<<" WARNING in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
            app_error()<<i <<" " <<j <<" " <<H1[i][j] <<" " <<H1[j][i] <<" "
                       <<hij[i][j] <<" " <<hij[j][i] <<" "
                       <<vn0[i][j] <<" " <<vn0[j][i] <<std::endl;
            //APP_ABORT("Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n");
          }
          H1[i][j] = 0.5*(H1[i][j]+ma::conj(H1[j][i]));
          H1[j][i] = ma::conj(H1[i][j]);
        }
      }

      return H1;
    }

    template<class Mat, class MatB>
    void energy(Mat&& E, MatB const& G, int k, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      MatB* Kr(nullptr);
      MatB* Kl(nullptr);
      energy(E,G,k,Kl,Kr,addH1,addEJ,addEXX);
    }

    // Kl and Kr must be in shared memory for this to work correctly
    template<class Mat, class MatB, class MatC, class MatD>
    void energy(Mat&& E, MatB const& Gc, int k, MatC* Kl, MatD* Kr, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      assert(E.size(1)>=3);
      assert(k >= 0 && k < haj.size());
      assert(k >= 0 && k < Vakbl_view.size());
      if(Gcloc.num_elements() < Gc.size(1) * Vakbl_view[k].size(0))
        Gcloc.reextent(iextensions<1u>{Vakbl_view[k].size(0)*Gc.size(1)});
      boost::multi::array_ref<SPComplexType,2> buff(Gcloc.data(),
                        {long(Vakbl_view[k].size(0)),long(Gc.size(1))});

      int nwalk = Gc.size(1);
      int getKr = Kr!=nullptr;
      int getKl = Kl!=nullptr;
      if(E.size(0) != nwalk || E.size(1) < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix dimensions \n");
      for(int n=0; n<nwalk; n++)
        std::fill_n(E[n].origin(),3,ComplexType(0.));
      if(addEJ and getKl)
        assert(Kl->size(0) == nwalk && Kl->size(1) == SpvnT[k].size(0));
      if(addEJ and getKr)
        assert(Kr->size(0) == nwalk && Kr->size(1) == SpvnT[k].size(0));

#if MIXED_PRECISION
      size_t mem_needs = Gc.num_elements();
      set_buffer(mem_needs);
      boost::multi::array_ref<SPComplexType,2> Gsp(to_address(SM_TMats.origin()), Gc.extensions());
      size_t i0, iN;
      std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(Gc.num_elements()),size_t(comm->size()));
      copy_n_cast(to_address(Gc.origin())+i0,iN-i0,to_address(Gsp.origin())+i0);
      comm->barrier();
#else
      auto& Gsp(Gc);  
#endif

      // one-body contribution
      if(addH1) {
        boost::multi::array_cref<ComplexType,1> haj_ref(to_address(haj[k].origin()), iextensions<1u>{haj[k].num_elements()});
        ma::product(ComplexType(1.),ma::T(Gc),haj_ref,ComplexType(1.),E(E.extension(0),0));
        for(int i=0; i<nwalk; i++)
          E[i][0] += E0;
      }

      // move calculation of H1 here
      if(addEXX) {
        shm::calculate_energy(std::forward<Mat>(E),Gsp,buff,Vakbl_view[k]);
      }

      if(separateEJ && addEJ) {
        using ma::T;
        if(Gcloc.num_elements() < SpvnT[k].size(0) * Gc.size(1))
          Gcloc.reextent(iextensions<1u>{SpvnT[k].size(0)*Gc.size(1)});
        assert(SpvnT_view[k].size(1) == Gc.size(0));
        RealType scl = (walker_type==CLOSED?4.0:1.0);
        // SpvnT*G
        boost::multi::array_ref<SPComplexType,2> v_(Gcloc.origin()+
                                            SpvnT_view[k].local_origin()[0]*Gc.size(1),
                                        {long(SpvnT_view[k].size(0)),long(Gc.size(1))});
        ma::product(SpvnT_view[k], Gsp, v_);
        if(getKl || getKr) {
          for(int wi=0; wi<Gc.size(1); wi++) {
            auto _v_ = v_(v_.extension(0),wi); 
            if(getKl) {
              auto Kli = (*Kl)[wi];
              for(int ki=0, qi = SpvnT_view[k].local_origin()[0]; ki<_v_.size(); ki++, qi++)
                Kli[qi] = static_cast<ComplexType>(_v_[ki]);
            }
            if(getKr) {
              auto Kri = (*Kr)[wi];
              for(int ki=0, qi = SpvnT_view[k].local_origin()[0]; ki<_v_.size(); ki++, qi++)
                Kri[qi] = static_cast<ComplexType>(_v_[ki]);
            }
          }
        }
        for(int wi=0; wi<Gc.size(1); wi++)
          E[wi][2] = 0.5*scl*static_cast<ComplexType>(ma::dot(v_(v_.extension(0),wi),v_(v_.extension(0),wi)));
      }
#if MIXED_PRECISION
#endif

    }

    template<class... Args>
    void fast_energy(Args&&... args)
    {
      APP_ABORT(" Error: fast_energy not implemented in SparseTensor. \n");
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      assert( Spvn.size(1) == X.size(0) );
      assert( Spvn.size(0) == v.size(0) );

#if MIXED_PRECISION
      size_t mem_needs = X.num_elements()+v.num_elements();  
      set_buffer(mem_needs);
      boost::multi::array_ref<SPComplexType,1> vsp(to_address(SM_TMats.origin()), v.extensions());  
      boost::multi::array_ref<SPComplexType,1> Xsp(vsp.origin()+vsp.num_elements(), X.extensions());  
      size_t i0, iN;
      std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(X.num_elements()),size_t(comm->size()));
      copy_n_cast(to_address(X.origin())+i0,iN-i0,to_address(Xsp.origin())+i0);
      boost::multi::array_ref<SPComplexType,1> v_(to_address(vsp.origin()) + Spvn_view.local_origin()[0],
                                        iextensions<1u>{Spvn_view.size(0)});
      comm->barrier();
      ma::product(SPValueType(a),Spvn_view,Xsp,SPValueType(c),v_);
      copy_n_cast(to_address(v_.origin()),v_.num_elements(),to_address(v.origin())+Spvn_view.local_origin()[0]);
      comm->barrier();
#else
      using Type = typename std::decay<MatB>::type::element;
      // Spvn*X
      boost::multi::array_ref<Type,1> v_(to_address(v.origin()) + Spvn_view.local_origin()[0],
                                        iextensions<1u>{Spvn_view.size(0)});
      ma::product(SPValueType(a),Spvn_view,X,SPValueType(c),v_);
#endif
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      assert( Spvn.size(1) == X.size(0) );
      assert( Spvn.size(0) == v.size(0) );
      assert( X.size(1) == v.size(1) );
#if MIXED_PRECISION
      size_t mem_needs = X.num_elements()+v.num_elements();
      set_buffer(mem_needs);
      boost::multi::array_ref<SPComplexType,2> vsp(to_address(SM_TMats.origin()), v.extensions());
      boost::multi::array_ref<SPComplexType,2> Xsp(vsp.origin()+vsp.num_elements(), X.extensions());
      size_t i0, iN;
      std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(X.num_elements()),size_t(comm->size()));
      copy_n_cast(to_address(X.origin())+i0,iN-i0,to_address(Xsp.origin())+i0);
      boost::multi::array_ref<SPComplexType,2> v_(to_address(vsp[Spvn_view.local_origin()[0]].origin()),
                                        {long(Spvn_view.size(0)),long(vsp.size(1))});
      comm->barrier();
      ma::product(SPValueType(a),Spvn_view,Xsp,SPValueType(c),v_);
      copy_n_cast(to_address(v_.origin()),v_.num_elements(),
                  to_address(v[Spvn_view.local_origin()[0]].origin()));
      comm->barrier();
#else
      using Type = typename std::decay<MatB>::type::element;
      // Spvn*X
      boost::multi::array_ref<Type,2> v_(to_address(v[Spvn_view.local_origin()[0]].origin()),
                                        {long(Spvn_view.size(0)),long(v.size(1))});
      ma::product(SPValueType(a),Spvn_view,X,SPValueType(c),v_);
#endif
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
      if(not separateEJ) k=0;
      assert( SpvnT[k].size(1) == G.size(0) );
      assert( SpvnT[k].size(0) == v.size(0) );

#if MIXED_PRECISION
      size_t mem_needs = G.num_elements()+v.num_elements();
      set_buffer(mem_needs);
      boost::multi::array_ref<SPComplexType,1> vsp(to_address(SM_TMats.origin()), v.extensions());
      boost::multi::array_ref<SPComplexType,1> Gsp(vsp.origin()+vsp.num_elements(), G.extensions());
      size_t i0, iN;
      std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(G.num_elements()),size_t(comm->size()));
      using std::copy_n;
      copy_n(to_address(G.origin())+i0,iN-i0,to_address(Gsp.origin())+i0);
      boost::multi::array_ref<SPComplexType,1> v_(to_address(vsp.origin()) + SpvnT_view[k].local_origin()[0],
                                        iextensions<1u>{SpvnT_view[k].size(0)});
      comm->barrier();
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view[k], Gsp, SpT2(c), v_);
      copy_n(to_address(v_.origin()),v_.num_elements(),
                  to_address(v.origin()) + SpvnT_view[k].local_origin()[0]);
      comm->barrier();
#else
      using Type = typename std::decay<MatB>::type::element ;
      // SpvnT*G
      boost::multi::array_ref<Type,1> v_(to_address(v.origin()) + SpvnT_view[k].local_origin()[0],
                                        iextensions<1u>{SpvnT_view[k].size(0)});
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view[k], G, SpT2(c), v_);
#endif
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
      if(not separateEJ) k=0;
      assert( SpvnT[k].size(1) == G.size(0) );
      assert( SpvnT[k].size(0) == v.size(0) );
      assert( G.size(1) == v.size(1) );

#if MIXED_PRECISION
      size_t mem_needs = G.num_elements()+v.num_elements();
      set_buffer(mem_needs);
      boost::multi::array_ref<SPComplexType,2> vsp(to_address(SM_TMats.origin()), v.extensions());
      boost::multi::array_ref<SPComplexType,2> Gsp(vsp.origin()+vsp.num_elements(), G.extensions());
      size_t i0, iN;
      std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(G.num_elements()),size_t(comm->size()));
      copy_n(to_address(G.origin())+i0,iN-i0,to_address(Gsp.origin())+i0);
      boost::multi::array_ref<SPComplexType,2> v_(to_address(vsp[SpvnT_view[k].local_origin()[0]].origin()),
                                        {long(SpvnT_view[k].size(0)),long(vsp.size(1))});
      comm->barrier();
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view[k], Gsp, SpT2(c), v_);
      copy_n(to_address(v_.origin()),v_.num_elements(),
                  to_address(to_address(v[SpvnT_view[k].local_origin()[0]].origin()))); 
      comm->barrier();
#else
      using Type = typename std::decay<MatB>::type::element ;
      // SpvnT*G
      boost::multi::array_ref<Type,2> v_(to_address(v[SpvnT_view[k].local_origin()[0]].origin()),
                                        {long(SpvnT_view[k].size(0)),long(v.size(1))});
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view[k], G, SpT2(c), v_);
#endif
    }

    bool distribution_over_cholesky_vectors() const{ return true; }
    int number_of_ke_vectors() const{ return Spvn.size(1); }
    int local_number_of_cholesky_vectors() const{ return Spvn.size(1); }
    int global_number_of_cholesky_vectors() const{ return global_nCV; }
    int global_origin_cholesky_vector() const{ return Spvn.global_origin()[1];}

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const{return false;}
    bool transposed_G_for_E() const{return false;}
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const{return false;}

    bool fast_ph_energy() const { return false; }

    boost::multi::array<ComplexType,2> getHSPotentials()
    {
      int nchol = global_number_of_cholesky_vectors();
      boost::multi::array<ComplexType,2> HSPot({nchol,nchol});
      csr::CSR2MA('T',Spvn,HSPot);
      return HSPot;
    }

  private:

    communicator* comm;

    WALKER_TYPES walker_type;

    int global_nCV;

    bool separateEJ;

    ValueType E0;

    // bare one body hamiltonian
    CMatrix hij;

    // (potentially half rotated) one body hamiltonian
    std::vector<CVector> haj;

    // sparse 2-body 2-electron integrals in matrix form
    std::vector<T1shm_csr_matrix> Vakbl;

    // sparse sub-matrix view
    std::vector<T1shm_csr_matrix_view> Vakbl_view;

    // Cholesky factorization of 2-electron integrals in sparse matrix form
    Vshm_csr_matrix Spvn;

    // sparse sub-matrix view
    Vshm_csr_matrix_view Spvn_view;

    // Cholesky factorization of 2-electron integrals in sparse matrix form
    std::vector<T2shm_csr_matrix> SpvnT;

    // sparse sub-matrix view
    std::vector<T2shm_csr_matrix_view> SpvnT_view;

    // one-body piece of Hamiltonian factorization
    CMatrix vn0;

    // local storage
    SpVector Gcloc;

    shmSpVector SM_TMats;    

    void set_buffer(size_t N) {
      if(SM_TMats.num_elements() < N) {
        SM_TMats.reextent(iextensions<1u>{N});
        using std::fill_n;
        if(comm->root()) 
          fill_n(to_address(SM_TMats.origin()),N,SPComplexType(0.0));
        comm->barrier();
      }
    }

};

}

}

#endif
