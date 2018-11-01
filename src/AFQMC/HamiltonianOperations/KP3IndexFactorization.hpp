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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KP3INDEXFACTORIZATION_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KP3INDEXFACTORIZATION_HPP

#include <vector>
#include <type_traits>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "alf/boost/mpi3/shared_communicator.hpp"
#include "AFQMC/multi/array.hpp"
#include "AFQMC/multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/taskgroup.h"

namespace qmcplusplus
{

namespace afqmc
{

class KP3IndexFactorization
{
  
  using CVector = boost::multi_array<ComplexType,1>;  
  using SpVector = boost::multi_array<SPComplexType,1>;  
  using CMatrix = boost::multi::array<ComplexType,2>;  
  using CMatrix_cref = boost::multi::const_array_ref<ComplexType,2>;
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2>;
  using SpMatrix_cref = boost::multi::const_array_ref<SPComplexType,2>;
  using SpMatrix_ref = boost::multi::array_ref<SPComplexType,2>;
  using C3Tensor = boost::multi::array<ComplexType,3>;
  using SpMatrix = boost::multi::array<SPComplexType,2>;  
  using Sp3Tensor = boost::multi::array<SPComplexType,3>;
  using shmCVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;  
  using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;  
  using shmIMatrix = boost::multi::array<int,2,shared_allocator<int>>;  
  using shmC3Tensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;  
  using shmSpVector = boost::multi::array<SPComplexType,1,shared_allocator<SPComplexType>>;  
  using shmSpMatrix = boost::multi::array<SPComplexType,2,shared_allocator<SPComplexType>>;  
  using shmSp3Tensor = boost::multi::array<SPComplexType,3,shared_allocator<SPComplexType>>;  
  using communicator = boost::mpi3::shared_communicator;
  using this_t = KP3IndexFactorization;

  public:

    KP3IndexFactorization(communicator& c_,
                 WALKER_TYPES type,
                 std::vector<int>&& nopk_,
                 std::vector<int>&& ncholpQ_,
                 shmIMatrix&& nelpk_,
                 shmIMatrix&& QKToK2_,   
                 shmC3Tensor&& hij_,       
                 shmCMatrix&& h1, 
                 //std::vector<shmCVector>&& h1, 
                 std::vector<shmSpMatrix>&& vik, 
                 std::vector<shmSpMatrix>&& vak, 
                 std::vector<shmSpMatrix>&& vlb, 
                 shmC3Tensor&& vn0_, 
                 ValueType e0_,
                 int gncv): 
        comm(std::addressof(c_)),
        walker_type(type),
        global_nCV(gncv),
        E0(e0_),
        H1_Qij(std::move(hij_)),
        haj(std::move(h1)),
        nopk(std::move(nopk_)),
        ncholpQ(std::move(ncholpQ_)),
        nelpk(std::move(nelpk_)),
        QKToK2(std::move(QKToK2_)),
        LQKikn(std::move(vik)),
        LQKank(std::move(vak)),
        LQKlnb(std::move(vlb)),
        vn0_Qij(std::move(vn0_)),
        SM_TMats({1,1},shared_allocator<SPComplexType>{c_}),
        TMats({1,1})
    {
        local_nCV = local_number_of_cholesky_vectors();
    }

    ~KP3IndexFactorization() {}
    
    KP3IndexFactorization(const KP3IndexFactorization& other) = delete;
    KP3IndexFactorization& operator=(const KP3IndexFactorization& other) = delete;
    KP3IndexFactorization(KP3IndexFactorization&& other) = default; 
    KP3IndexFactorization& operator=(KP3IndexFactorization&& other) = default;

    boost::multi_array<ComplexType,2> getOneBodyPropagatorMatrix(TaskGroup_& TG, boost::multi_array<ComplexType,1> const& vMF) {

/*
    // switch to H1_Qij
      int NMO = hij.shape()[0];  
      // in non-collinear case with SO, keep SO matrix here and add it
      // for now, stay collinear
      boost::multi_array<ComplexType,2> H1(extents[NMO][NMO]);

      // add sum_n vMF*Spvn, vMF has local contribution only!
      boost::multi::array_ref<ComplexType,1> H1D(std::addressof(*H1.origin()),{NMO*NMO});
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
*/
      return boost::multi_array<ComplexType,2>(extents[1][1]);
    }

    template<class Mat, class MatB>
    void energy(Mat&& E, MatB const& G, int k, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      MatB* Kr(nullptr);
      MatB* Kl(nullptr);
      energy(E,G,k,Kl,Kr,addH1,addEJ,addEXX);
    }

    // Kl and Kr must be in shared memory for this to work correctly  
    template<class Mat, class MatB, class MatC, class MatD>
    void energy(Mat&& E, MatB const& Gc, int nd, MatC* Kl, MatD* Kr, bool addH1=true, bool addEJ=true, bool addEXX=true) {

      int nkpts = nopk.size(); 
      assert(E.shape()[1]>=3);
      assert(nd >= 0 && nd < nelpk.size());  

      int nwalk = Gc.shape()[0];
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nocca_tot = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+nkpts,0);
      int nocca_max = *std::max_element(nelpk[nd].begin(),nelpk[nd].begin()+nkpts);
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      int noccb_tot = 0;
      if(walker_type==COLLINEAR) noccb_tot = std::accumulate(nelpk[nd].begin()+nkpts,
                                      nelpk[nd].begin()+2*nkpts,0);
      int getKr = Kr!=nullptr;
      int getKl = Kl!=nullptr;
      if(E.shape()[0] != nwalk || E.shape()[1] < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix dimensions \n");

      if(addEJ and getKl)
        assert(Kl->shape()[0] == nwalk && Kl->shape()[1] == local_nCV); 
      if(addEJ and getKr)
        assert(Kr->shape()[0] == nwalk && Kr->shape()[1] == local_nCV); 

      for(int n=0; n<nwalk; n++) 
        std::fill_n(E[n].origin(),3,ComplexType(0.));

      assert(Gc.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      boost::multi::const_array_ref<ComplexType,3> G3Da(std::addressof(*Gc.origin()),
                                                        {nwalk,nocca_tot,nmo_tot} ); 
      boost::multi::const_array_ref<ComplexType,3> G3Db(std::addressof(*Gc.origin())+
                                                        G3Da.num_elements()*(nspin-1),
                                                        {nwalk,noccb_tot,nmo_tot} ); 

      // one-body contribution
      // haj[ndet*nkpts][nocc*nmo]
      // not parallelized for now, since it would require customization of Wfn 
      if(addH1) {
        int na=0, nj=0, nb=0;
        for(int n=0; n<nwalk; n++)
          E[n][0] = E0;  
        for(int K=0; K<nkpts; ++K) {
          boost::multi::array_ref<ComplexType,2> haj_K(std::addressof(*haj[nd*nkpts+K].origin()),
                                                      {nelpk[nd][K],nopk[K]}); 
          for(int n=0; n<nwalk; ++n) { 
            ComplexType E_(0.0);
            for(int a=0; a<nelpk[nd][K]; ++a) 
              //E_ += ma::dot(G3Da[n][na+a]({nj,nj+nopk[K]}),haj_K[a]);
              E_ += ma::dot(G3Da[n][na+a].sliced(nj,nj+nopk[K]),haj_K[a]);
            E[n][0] += E_;
          }    
          na+=nelpk[nd][K];
          if(walker_type==COLLINEAR) {
            boost::multi::array_ref<ComplexType,2> haj_Kb(haj_K.origin()+haj_K.num_elements(),
                                                          {nelpk[nd][nkpts+K],nopk[K]});
            for(int n=0; n<nwalk; ++n) {
              ComplexType E_(0.0);
              for(int a=0; a<nelpk[nd][nkpts+K]; ++a)
                E_ += ma::dot(G3Db[n][nb+a]({nj,nj+nopk[K]}),haj_Kb[a]);
              E[n][0] += E_;
            }
            nb+=nelpk[nd][nkpts+K];
          }  
          nj+=nopk[K];  
        }
      }

      // move calculation of H1 here	
      // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
      //       Not sure how to do it for COLLINEAR.
      if(addEXX) {  

        size_t local_memory_needs = nocca_tot*nchol_max*nmo_max + nmo_max*nocca_max + 10;
        if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});

        RealType scl = (walker_type==CLOSED?2.0:1.0);
        size_t nqk=1;  // start count at 1 to "offset" the calcuation of E1 done at root
        for(int Q=0; Q<nkpts; ++Q) {             // momentum conservation index   
          for(int Kl=0; Kl<nkpts; ++Kl) {           // K is the index of the kpoint pair of (l,k)
            for(int n=0; n<nwalk; ++n) {
              if((nqk++)%comm->size() == comm->rank()) { 
                int nchol = ncholpQ[Q];
                int nl = nopk[Kl];
                int nl0 = std::accumulate(nopk.begin(),nopk.begin()+Kl,0);
                int nb = nelpk[nd][QKToK2[Q][Kl]];
                int nb0 = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+QKToK2[Q][Kl],0);

                SpMatrix_ref Llnb(std::addressof(*LQKlnb[Q][Kl].origin()), 
                                                 {nl,nchol*nb}); 
                SpMatrix_ref TAnb(TMats.origin(),{nocca_tot,nchol*nb}); 

                // T(A,n,b) = sum_l_in_K G(A,l) LQKlnb[Q][K](l,n,b)
                ma::product(G3Da[n]({0,nocca_tot},{nl0,nl0+nl}),Llnb,TAnb);  

                int na0 = 0;
                for(int Ka=0; Ka<nkpts; ++Ka) {  
                  int na = nelpk[nd][Ka]; 
                  int nk = nopk[QKToK2[Q][Ka]];
                  int nk0 = std::accumulate(nopk.begin(),nopk.begin()+QKToK2[Q][Ka],0);

                  SpMatrix_ref Fbk(TMats.origin()+TAnb.num_elements(),{nb,nk}); 
                  SpMatrix_ref Lank(std::addressof(*LQKank[Q][Ka].origin()), 
                                                 {na*nchol,nk}); 
                  SpMatrix_ref Tanb(TAnb[na0].origin(),{na*nchol,nb}); 
              
                  // F(K1,b,k) = sum_a_in_K1 sum_n  T(a,n,b) * LQKank(a,n,k) 
                  ma::product(ma::T(Tanb),Lank,Fbk);

                  // EXX += sum_K1 sum_b_in_Q(K) sum_k_in_Q(K1) F(K1,b,k) * G(b,k)
                  ComplexType E_(0.0);
                  for(int b=0; b<nb; ++b) 
                    E_ += ma::dot(Fbk[b],G3Da[n][nb0+b]({nk0,nk0+nk}));  
                  E[n][1] += scl*E_;
                  na0 += na;  
                }
              }
              if(walker_type!=COLLINEAR) continue; 
              if((nqk++)%comm->size() == comm->rank()) {
              }
            }
          }
        }
      }  

/*
      if(addEJ && not addEXX) {
        // just call routine to calculate vbias!
        using ma::T;
        if(Gcloc.num_elements() < SpvnT[k].shape()[0] * Gc.shape()[1])
          Gcloc.resize(extents[SpvnT[k].shape()[0]*Gc.shape()[1]]);
        assert(SpvnT_view[k].shape()[1] == Gc.shape()[0]);
        RealType scl = (walker_type==CLOSED?4.0:1.0); 
        // SpvnT*G
        boost::multi_array_ref<T2,2> v_(Gcloc.origin()+
                                            SpvnT_view[k].local_origin()[0]*Gc.shape()[1],
                                        extents[SpvnT_view[k].shape()[0]][Gc.shape()[1]]);
        ma::product(SpvnT_view[k], Gc, v_); 
        if(getKl || getKr) { 
          for(int wi=0; wi<Gc.shape()[1]; wi++) {
            auto _v_ = v_[indices[range_t()][wi]];
            if(getKl) {
              auto Kli = (*Kl)[wi];
              for(int ki=0, qi = SpvnT_view[k].local_origin()[0]; ki<_v_.size(); ki++, qi++)
                Kli[qi] = _v_[ki];
            }
            if(getKr) {
              auto Kri = (*Kr)[wi];
              for(int ki=0, qi = SpvnT_view[k].local_origin()[0]; ki<_v_.size(); ki++, qi++)
                Kri[qi] = _v_[ki];
            }
          }
        }
        for(int wi=0; wi<Gc.shape()[1]; wi++) 
          E[wi][2] = 0.5*scl*ma::dot(v_[indices[range_t()][wi]],v_[indices[range_t()][wi]]); 
      }
*/
    }

    template<class... Args>
    void fast_energy(Args&&... args)
    {
      APP_ABORT(" Error: fast_energy not implemented in KP3IndexFactorization. \n"); 
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vHS(const MatA& X, MatB&& v, double a=1., double c=0.) {
/*
      assert( Spvn.shape()[1] == X.shape()[0] );
      assert( Spvn.shape()[0] == v.shape()[0] );
      using Type = typename std::decay<MatB>::type::element;

      // Spvn*X 
      boost::multi_array_ref<Type,1> v_(v.origin() + Spvn_view.local_origin()[0], 
                                        extents[Spvn_view.shape()[0]]);
      ma::product(SPValueType(a),Spvn_view,X,SPValueType(c),v_);
*/
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(const MatA& X, MatB&& v, double a=1., double c=0.) {
/*
      assert( Spvn.shape()[1] == X.shape()[0] );
      assert( Spvn.shape()[0] == v.shape()[0] );
      assert( X.shape()[1] == v.shape()[1] );
      using Type = typename std::decay<MatB>::type::element;

      // Spvn*X 
      boost::multi_array_ref<Type,2> v_(v[Spvn_view.local_origin()[0]].origin(), 
                                        extents[Spvn_view.shape()[0]][v.shape()[1]]);
      ma::product(SPValueType(a),Spvn_view,X,SPValueType(c),v_);
*/
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
/*
      if(not separateEJ) k=0;
      assert( SpvnT[k].shape()[1] == G.shape()[0] );
      assert( SpvnT[k].shape()[0] == v.shape()[0] );
      using Type = typename std::decay<MatB>::type::element ;

      // SpvnT*G
      boost::multi_array_ref<Type,1> v_(v.origin() + SpvnT_view[k].local_origin()[0], 
                                        extents[SpvnT_view[k].shape()[0]]);
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view[k], G, SpT2(c), v_);
*/
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
/*
      if(not separateEJ) k=0;
      assert( SpvnT[k].shape()[1] == G.shape()[0] );
      assert( SpvnT[k].shape()[0] == v.shape()[0] );   
      assert( G.shape()[1] == v.shape()[1] );
      using Type = typename std::decay<MatB>::type::element ;

      // SpvnT*G
      boost::multi_array_ref<Type,2> v_(v[SpvnT_view[k].local_origin()[0]].origin(), 
                                        extents[SpvnT_view[k].shape()[0]][v.shape()[1]]);
      if(walker_type==CLOSED) a*=2.0;
      ma::product(SpT2(a), SpvnT_view[k], G, SpT2(c), v_);
*/
    }

    bool distribution_over_cholesky_vectors() const{ return true; }
    int number_of_ke_vectors() const{ return 0; }
    int local_number_of_cholesky_vectors() const{
        int res=0;
        for(auto& LQ: LQKikn) res += LQ.shape()[1];
        return res;
    } 
    int global_number_of_cholesky_vectors() const{ return global_nCV; }

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const{return false;}
    bool transposed_G_for_E() const{return true;} 
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const{return false;} 

    bool fast_ph_energy() const { return false; }

  private:

    communicator* comm;

    WALKER_TYPES walker_type;

    int global_nCV;
    int local_nCV;

    ValueType E0;

    // bare one body hamiltonian
    shmC3Tensor H1_Qij;

    // (potentially half rotated) one body hamiltonian
    shmCMatrix haj;
    //std::vector<shmCVector> haj;

    // number of orbitals per k-point
    std::vector<int> nopk;

    // number of cholesky vectors per Q-point
    std::vector<int> ncholpQ;

    // number of electrons per k-point
    // nelpk[ndet][nspin*nkpts]
    shmIMatrix nelpk;

    // maps (Q,K) --> k2
    shmIMatrix QKToK2; 

    //Cholesky Tensor Lik[Q][nk][i][k][n]
    std::vector<shmSpMatrix> LQKikn;

    // half-tranformed Cholesky tensor
    std::vector<shmSpMatrix> LQKank;

    // half-tranformed Cholesky tensor
    std::vector<shmSpMatrix> LQKlnb;

    // one-body piece of Hamiltonian factorization
    shmC3Tensor vn0_Qij;

    // shared buffer space
    // using matrix since there are issues with vectors
    shmSpMatrix SM_TMats;
    SpMatrix TMats;

    void set_shm_buffer(size_t N) {
      if(SM_TMats.num_elements() < N) 
        SM_TMats.reextent({N,1});
    }

};

}

}

#endif
