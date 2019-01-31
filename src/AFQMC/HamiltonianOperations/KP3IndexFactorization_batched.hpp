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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KP3INDEXFACTORIZATION_BATCHED_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KP3INDEXFACTORIZATION_BATCHED_HPP

#include <Utilities/NewTimer.h>
#include <vector>
#include <type_traits>
#include <random>
#include <algorithm>

#include "Configuration.h"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/helpers/batched_operations.hpp"
#include "AFQMC/Numerics/helpers/tensor_transpositions.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class KP3IndexFactorization_batched
{


  // allocators
  using Allocator = device_allocator<ComplexType>;
  using SpAllocator = device_allocator<SPComplexType>; 
  using BAllocator = device_allocator<bool>; 
  using IAllocator = device_allocator<int>; 
  using Allocator_shared = node_allocator<ComplexType>; 
  using SpAllocator_shared = node_allocator<SPComplexType>; 
  using IAllocator_shared = node_allocator<int>; 

  // type defs
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;
  using sp_pointer = typename SpAllocator::pointer;
  using const_sp_pointer = typename SpAllocator::const_pointer;
  using pointer_shared = typename Allocator_shared::pointer;
  using const_pointer_shared = typename Allocator_shared::const_pointer;
  using sp_pointer_shared = typename SpAllocator_shared::pointer;
  using const_sp_pointer_shared = typename SpAllocator_shared::const_pointer;

  using stdIVector =  boost::multi::array<int,1>;

  using IVector =  boost::multi::array<int,1,IAllocator>;
  using BoolMatrix =  boost::multi::array<bool,2,BAllocator>;
  using CVector = ComplexVector<Allocator>;  
  using IMatrix = IntegerMatrix<IAllocator>;  
  using CMatrix = ComplexMatrix<Allocator>;  
  using C3Tensor = boost::multi::array<ComplexType,3,Allocator>;

  using SpVector = SPComplexVector<SpAllocator>;  
  using SpMatrix = SPComplexMatrix<SpAllocator>;  
  using Sp3Tensor = boost::multi::array<SPComplexType,3,SpAllocator>;

  using CMatrix_cref = boost::multi::array_ref<ComplexType const,2,const_pointer>;
  using CVector_ref = ComplexVector_ref<pointer>;
  using CMatrix_ref = ComplexMatrix_ref<pointer>;
  using C3Tensor_cref = boost::multi::array_ref<ComplexType const,3,const_pointer>;

  using SpMatrix_cref = boost::multi::array_ref<SPComplexType const,2,sp_pointer>;
  using SpVector_ref = SPComplexVector_ref<sp_pointer>;
  using SpMatrix_ref = SPComplexMatrix_ref<sp_pointer>;
  using Sp3Tensor_ref = SPComplex3Tensor_ref<sp_pointer>;
  using Sp4Tensor_ref = SPComplexArray_ref<4,sp_pointer>;
  using Sp5Tensor_ref = SPComplexArray_ref<5,sp_pointer>;

  using shmCVector = ComplexVector<Allocator_shared>;  
  using shmCMatrix = ComplexMatrix<Allocator_shared>;  
  using shmIMatrix = IntegerMatrix<IAllocator_shared>;  
  using shmC3Tensor = Complex3Tensor<Allocator_shared>;  

  using mpi3C3Tensor = Complex3Tensor<shared_allocator<ComplexType>>;

  using shmSpVector = SPComplexVector<SpAllocator_shared>;  
  using shmSpMatrix = SPComplexMatrix<SpAllocator_shared>;  
  using shmSp3Tensor = SPComplex3Tensor<SpAllocator_shared>;  

  public:

// NOTE: careful with nocc_max, not consistently defined!!!

    // since arrays can be in host, can't assume that types are consistent
    template<class shmCMatrix_, class shmSpMatrix_>
    KP3IndexFactorization_batched(
                 WALKER_TYPES type,
                 stdIVector&& nopk_,
                 stdIVector&& ncholpQ_,
                 stdIVector&& kminus_,
                 boost::multi::array<int,2>&& nelpk_,
                 boost::multi::array<int,2>&& QKToK2_,   
                 mpi3C3Tensor&& hij_,       
                 shmCMatrix_&& h1, 
                 std::vector<shmSpMatrix_>&& vik, 
                 std::vector<shmSpMatrix_>&& vak, 
                 std::vector<shmSpMatrix_>&& vakn, 
                 mpi3C3Tensor&& vn0_, 
                 std::vector<RealType>&& gQ_,
                 int nsampleQ_,
                 ValueType e0_,
                 Allocator const& alloc_,
                 Allocator_shared const& alloc_shared_,
                 int gncv): 
        allocator_(alloc_),
        sp_allocator_(alloc_),
        allocator_shared_(alloc_shared_),
        sp_allocator_shared_(alloc_shared_),
        iallocator_shared_(alloc_shared_),
        walker_type(type),
        global_nCV(gncv),
        E0(e0_),
        H1(std::move(hij_)),
        haj(std::move(h1)),
        nopk(std::move(nopk_)),
        ncholpQ(std::move(ncholpQ_)),
        kminus(std::move(kminus_)),
        nelpk(std::move(nelpk_)),
        QKToK2(std::move(QKToK2_)),
        LQKikn(std::move(move_vector<shmSpMatrix>(std::move(vik)))),
        LQKank(std::move(move_vector<shmSpMatrix>(std::move(vak)))),
        LQKakn(std::move(move_vector<shmSpMatrix>(std::move(vakn)))),
        vn0(std::move(vn0_)),
        gQ(std::move(gQ_)),
        Qwn({1,1}),
        generator(),
        distribution(gQ.begin(),gQ.end()), 
        nsampleQ(nsampleQ_),
        SM_TMats({1,1},sp_allocator_shared_),
        TMats({1,1},sp_allocator_),
        IMats({1,1},IAllocator{allocator_}),
        KKTransID( {nopk.size()+1,nopk.size()}, BAllocator{allocator_}),
        dev_nopk(nopk),
        dev_i0pk( typename IVector::extensions_type{nopk.size()}, IAllocator{allocator_}),
        dev_kminus(kminus),
        dev_ncholpQ(ncholpQ),
        dev_ncholpQ0( typename IVector::extensions_type{nopk.size()}, IAllocator{allocator_}),
        dev_nelpk(nelpk),
        dev_a0pk( typename IMatrix::extensions_type{nelpk.size(0),nelpk.size(1)}, 
                                                                    IAllocator{allocator_}),
        dev_QKToK2(QKToK2), 
        EQ(nopk.size()+2)
    {
      using std::fill_n;
      using std::copy_n;
      nocc_max = *std::max_element(nelpk.origin(),nelpk.origin()+nelpk.num_elements());
      setup_timers(Timers, THCTimerNames, timer_level_coarse);
      local_nCV = std::accumulate(ncholpQ.begin(),ncholpQ.end(),0); 
      fill_n(EQ.data(),EQ.size(),0);
      int nkpts = nopk.size(); 
      Q0=-1;  // if K=(0,0,0) exists, store index here
      for(int Q=0; Q<nkpts; Q++) {
        if(kminus[Q]==Q) {
          bool found=true;
          for(int KI=0; KI<nkpts; KI++)
            if(KI != QKToK2[Q][KI]) {
              found = false;
              break;
            }
          if(found) {
            Q0=Q;
            break;
          }
        }
      }
      // setup dev integer arrays
      std::vector<int> i0(nkpts);
      i0[0]=0;
      for(int i=1; i<nkpts; i++)
        i0[i] = i0[i-1]+2*ncholpQ[i-1];
      copy_n(i0.data(),nkpts,dev_ncholpQ0.origin());
      // dev_nopk  
      i0[0]=0;
      for(int i=1; i<nkpts; i++) 
        i0[i] = i0[i-1]+nopk[i-1];
      copy_n(i0.data(),nkpts,dev_i0pk.origin());  
      // dev_nelpk  
      for(int n=0; n<nelpk.size(0); n++) {
        i0[0]=0;
        for(int i=1; i<nkpts; i++)
          i0[i] = i0[i-1]+nelpk[n][i-1];
        copy_n(i0.data(),nkpts,dev_a0pk[n].origin());
        if(walker_type==COLLINEAR) {
          i0[0]=0;
          for(int i=1; i<nkpts; i++)
            i0[i] = i0[i-1]+nelpk[n][nkpts+i-1];
          copy_n(i0.data(),nkpts,dev_a0pk[n].origin()+nkpts);
        }
      }
      // setup copy/transpose tags
      boost::multi::array<bool,2> KKid({nkpts+1,nkpts});  
      for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
        if(Q < kminus[Q] || Q==Q0) {
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            int QK = QKToK2[Q][K];
            KKid[K][QK] = false;
          }  
        } else if(Q > kminus[Q]) { // use L(-Q)(ki)*
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            int QK = QKToK2[Q][K];
            KKid[K][QK] = true;
          }
        } else { // Q==(-Q) and Q!=Q0
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            int QK = QKToK2[Q][K];
            if(K < QK)  KKid[K][QK] = false;
            else KKid[K][QK] = true;
          }
        }
        KKid[nkpts][Q] = true; // last term stores the second term of Q0
      }
      copy_n(KKid.origin(),KKid.num_elements(),KKTransID.origin());  
    }

    ~KP3IndexFactorization_batched() {}
    
    KP3IndexFactorization_batched(const KP3IndexFactorization_batched& other) = delete;
    KP3IndexFactorization_batched& operator=(const KP3IndexFactorization_batched& other) = delete;
    KP3IndexFactorization_batched(KP3IndexFactorization_batched&& other) = default; 
    KP3IndexFactorization_batched& operator=(KP3IndexFactorization_batched&& other) = default;

    // must have the same signature as shared classes, so keeping it with std::allocator
    boost::multi::array<ComplexType,2> getOneBodyPropagatorMatrix(TaskGroup_& TG, boost::multi::array<ComplexType,1> const& vMF) {

      int nkpts = nopk.size();
      int NMO = std::accumulate(nopk.begin(),nopk.end(),0);
      // in non-collinear case with SO, keep SO matrix here and add it
      // for now, stay collinear

      CVector vMF_(vMF);
      CVector P1D(extensions<1u>{NMO*NMO});
      fill_n(P1D.origin(),P1D.num_elements(),ComplexType(0));
      vHS(vMF_, P1D);
      TG.TG().all_reduce_in_place_n(to_address(P1D.origin()),P1D.num_elements(),std::plus<>());

      boost::multi::array<ComplexType,2> P1({NMO,NMO});
      copy_n(P1D.origin(),NMO*NMO,P1.origin());

      // add H1 + vn0 and symmetrize
      using std::conj;

      for(int K=0, nk0=0; K<nkpts; ++K) {
        for(int i=0, I=nk0; i<nopk[K]; i++, I++) {
          P1[I][I] += H1[K][i][i] + vn0[K][i][i];
          for(int j=i+1, J=I+1; j<nopk[K]; j++, J++) {
            P1[I][J] += H1[K][i][j] + vn0[K][i][j];
            P1[J][I] += H1[K][j][i] + vn0[K][j][i];
            // This is really cutoff dependent!!!
            if( std::abs( P1[I][J] - conj(P1[J][I]) ) > 1e-6 ) {
              app_error()<<" WARNING in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
              app_error()<<I <<" " <<J <<" " <<P1[I][J] <<" " <<P1[j][i] <<" "
                         <<H1[K][i][j] <<" " <<H1[K][j][i] <<" "
                         <<vn0[K][i][j] <<" " <<vn0[K][j][i] <<std::endl;
              //APP_ABORT("Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n");
            }
            P1[I][J] = 0.5*(P1[I][J]+conj(P1[J][I]));
            P1[J][I] = conj(P1[I][J]);
          }
        }
        nk0 += nopk[K];
      }
      return P1;
    }

    template<class Mat, class MatB>
    void energy(Mat&& E, MatB const& G, int k=0, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      MatB* Kr(nullptr);
      MatB* Kl(nullptr);
      energy(E,G,k,Kl,Kr,addH1,addEJ,addEXX);
    }

    // KEleft and KEright must be in shared memory for this to work correctly  
    template<class Mat, class MatB, class MatC, class MatD
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<Mat>())),
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>())),
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatC>())),
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatD>()))
            >
    void energy(Mat&& E, MatB const& Gc, int nd, MatC* KEleft, MatD* KEright, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      if(nsampleQ > 0)
        energy_sampleQ(E,Gc,nd,KEleft,KEright,addH1,addEJ,addEXX);
      else
        energy_exact(E,Gc,nd,KEleft,KEright,addH1,addEJ,addEXX);
    }

    // Catch-all routines for incompatible pointer types!
    void energy(...) {
      print_stacktrace  
      throw std::runtime_error("Error: Calling KP3IndexFactorization_batched::energy catch all.");
    }
    void vHS(...) {
      print_stacktrace  
      throw std::runtime_error("Error: Calling KP3IndexFactorization_batched::energy catch all.");
    }

    // KEleft and KEright must be in shared memory for this to work correctly  
    template<class Mat, class MatB, class MatC, class MatD>
    void energy_exact(Mat&& E, MatB const& Gc, int nd, MatC* KEleft, MatD* KEright, bool addH1=true, bool addEJ=true, bool addEXX=true) {

      using std::fill_n;
      using std::copy_n;
      int nkpts = nopk.size(); 
      assert(E.size(1)>=3);
      assert(nd >= 0 && nd < nelpk.size());  

      int nwalk = Gc.size(1);
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nocca_tot = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+nkpts,0);
      int nocca_max = *std::max_element(nelpk[nd].begin(),nelpk[nd].begin()+nkpts);
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      int nchol_tot = std::accumulate(ncholpQ.begin(),ncholpQ.end(),0);
      int noccb_tot = 0;
      if(walker_type==COLLINEAR) noccb_tot = std::accumulate(nelpk[nd].begin()+nkpts,
                                      nelpk[nd].begin()+2*nkpts,0);
      int getKr = KEright!=nullptr;
      int getKl = KEleft!=nullptr;
      if(E.size(0) != nwalk || E.size(1) < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix dimensions \n");

      size_t mem_needs(nwalk*nkpts*nkpts*nspin*nocca_max*nmo_max);
      size_t cnt(0);  
      if(addEJ) { 
        if(not getKr) mem_needs += nwalk*local_nCV;
        if(not getKl) mem_needs += nwalk*local_nCV;
      }
      set_shm_buffer(mem_needs);

      // messy
      sp_pointer_shared Krptr, Klptr; 
      size_t Knr=0, Knc=0;
      if(addEJ) {
        Knr=nwalk;
        Knc=local_nCV;
        cnt=0;
        if(getKr) {
          assert(KEright->size(0) == nwalk && KEright->size(1) == local_nCV); 
          assert(KEright->stride() == KEright->size(1));
          Krptr = make_device_ptr(KEright->origin()); 
        } else {
          Krptr = SM_TMats.origin(); 
          cnt += nwalk*local_nCV;
        }
        if(getKl) {
          assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV); 
          assert(KEleft->stride(0) == KEleft->size(1));
          Klptr = make_device_ptr(KEleft->origin());
        } else {
          Klptr = SM_TMats.origin()+cnt; 
          cnt += nwalk*local_nCV;
        }
        fill_n(Krptr,Knr*Knc,SPComplexType(0.0));
        fill_n(Klptr,Knr*Knc,SPComplexType(0.0));
      } else if(getKr or getKl) {
        APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
      }   
      SpMatrix_ref Kl(Klptr,{Knr,Knc});
      SpMatrix_ref Kr(Krptr,{Knr,Knc});

      for(int n=0; n<nwalk; n++) 
        fill_n(E[n].origin(),3,ComplexType(0.));

      assert(Gc.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      C3Tensor_cref G3Da(make_device_ptr(Gc.origin()),{nocca_tot,nmo_tot,nwalk} );
      C3Tensor_cref G3Db(make_device_ptr(Gc.origin())+G3Da.num_elements()*(nspin-1),
                                    {noccb_tot,nmo_tot,nwalk} );

      // later on, rewrite routine to loop over spins, to avoid storage of both spin
      // components simultaneously
      Timers[Timer_E1]->start();
      Sp4Tensor_ref GKK(SM_TMats.origin()+cnt,
                        {nspin,nkpts,nkpts,nwalk*nmo_max*nocc_max});
      GKaKjw_to_GKKwaj(G3Da,GKK[0],nelpk[nd].sliced(0,nkpts),dev_nelpk[nd],dev_a0pk[nd]);
      if(walker_type==COLLINEAR)  
        GKaKjw_to_GKKwaj(G3Db,GKK[1],nelpk[nd].sliced(nkpts,2*nkpts),
                                     dev_nelpk[nd].sliced(nkpts,2*nkpts),
                                     dev_a0pk[nd].sliced(nkpts,2*nkpts));
      Timers[Timer_E1]->stop();

      // one-body contribution
      // haj[ndet*nkpts][nocc*nmo]
      // not parallelized for now, since it would require customization of Wfn 
      if(addH1) {
        Timers[Timer_E2]->start();
        for(int n=0; n<nwalk; n++) 
          fill_n(E[n].origin(),1,ComplexType(E0));
        // must use Gc since GKK is is SP
        int na=0, nk=0, nb=0;
        for(int K=0; K<nkpts; ++K) {
#if defined(MIXED_PRECISION) 
          CMatrix_ref haj_K(haj[nd*nkpts+K].origin(),{nocc_max,nmo_max}); 
          for(int a=0; a<nelpk[nd][K]; ++a) 
            ma::product(ComplexType(1.),ma::T(G3Da[na+a].sliced(nk,nk+nopk[K])),
                                        haj_K[a].sliced(0,nopk[K]),
                        ComplexType(1.),E({0,nwalk},0));
          na+=nelpk[nd][K];
          if(walker_type==COLLINEAR) {
            boost::multi::array_ref<ComplexType,2> haj_Kb(haj_K.origin()+haj_K.num_elements(),
                                                      {nocc_max,nmo_max}); 
            for(int b=0; b<nelpk[nd][nkpts+K]; ++b) 
              ma::product(ComplexType(1.),ma::T(G3Db[nb+b].sliced(nk,nk+nopk[K])),
                                        haj_Kb[b].sliced(0,nopk[K]),
                        ComplexType(1.),E({0,nwalk},0));
            nb+=nelpk[nd][nkpts+K];
          }  
          nk+=nopk[K];  
#else
          nk = nopk[K];
          {
            na = nelpk[nd][K];
            CVector_ref haj_K(haj[nd*nkpts+K].origin(),{nocc_max*nmo_max});
            SpMatrix_ref Gaj(GKK[0][K][K].origin(),{nwalk,nocc_max*nmo_max});
            ma::product(ComplexType(1.),Gaj,haj_K,ComplexType(1.),E({0,nwalk},0));
          }
          if(walker_type==COLLINEAR) {
            na = nelpk[nd][nkpts+K];
            CVector_ref haj_K(haj[nd*nkpts+K].origin()+nocc_max*nmo_max,{nocc_max*nmo_max});
            SpMatrix_ref Gaj(GKK[1][K][K].origin(),{nwalk,nocc_max*nmo_max});
            ma::product(ComplexType(1.),Gaj,haj_K,ComplexType(1.),E({0,nwalk},0));
          }  
#endif
        }
        Timers[Timer_E2]->stop();
      }

      // move calculation of H1 here	
      // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
      //       Not sure how to do it for COLLINEAR.
      if(addEXX) {  

        // create some convention for batch_size 
        // simple implementation for now
        // taking 4Gbs for now, is this reasonable???
        size_t Bytes = size_t(4*1024)*size_t(1024*1024);
        Bytes /= size_t(2*nwalk*nocc_max*nocc_max*nchol_max); 
        size_t bz0 = std::max(size_t(1), size_t(std::floor(Bytes)));
        batch_size = std::min(bz0,size_t(nkpts*nkpts));
        int batch_cnt(0);
        using ma::gemmBatched;
        std::vector<sp_pointer> Aarray;
        std::vector<sp_pointer> Barray;
        std::vector<sp_pointer> Carray;
        Aarray.reserve(batch_size);
        Barray.reserve(batch_size);
        Carray.reserve(batch_size);
        std::vector<SPComplexType> scl_factors;
        scl_factors.reserve(batch_size);
        std::vector<int> kdiag;
        kdiag.reserve(batch_size);

        size_t local_memory_needs = batch_size*size_t(2*nwalk*nocc_max*nocc_max*nchol_max) + 
                                    batch_size; 
        if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});
        if(IMats.num_elements() < batch_size) IMats.reextent({batch_size,1});
        cnt=0; 
        SpVector_ref dev_scl_factors(TMats.origin()+cnt,{batch_size}); 
        cnt+=dev_scl_factors.num_elements();
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        size_t nqk=1;  
        for(int Q=0; Q<nkpts; ++Q) {
          bool haveKE=false;

          // simple implementation for now
          Aarray.clear();
          Barray.clear();
          Carray.clear();
          scl_factors.clear();  
          kdiag.clear();  
          batch_cnt=0;
          Sp3Tensor_ref T1(TMats.origin()+cnt,{2*batch_size,nwalk*nocc_max,nocc_max*nchol_max});

          for(int Ka=0; Ka<nkpts; ++Ka) {
            int K0 = ((Q==Q0)?0:Ka);
            for(int Kb=K0; Kb<nkpts; ++Kb) {
              int Qm = kminus[Q];
              int Qm_ = (Q==Q0?nkpts:Qm);
              int Kl_ = QKToK2[Qm][Kb];
              int Kk = QKToK2[Q][Ka];

              if(addEJ && Ka==Kb) 
                kdiag.push_back(batch_cnt);  

              Aarray.push_back(LQKank[nd*nspin*(nkpts+1)+Qm_][Kb].origin());
              Barray.push_back(GKK[0][Ka][Kl_].origin());
              Carray.push_back(T1[batch_cnt++].origin());
              Aarray.push_back(LQKank[nd*nspin*(nkpts+1)+Q][Ka].origin());
              Barray.push_back(GKK[0][Kb][Kk].origin());
              Carray.push_back(T1[batch_cnt++].origin());

              if(Q==Q0 || Ka==Kb)
                scl_factors.push_back(SPComplexType(-scl*0.5));
              else
                scl_factors.push_back(SPComplexType(-scl));

              if( 2*batch_cnt >= batch_size ) {
                Timers[Timer_E3]->start();
                gemmBatched('T','N',nocc_max*nchol_max,nwalk*nocc_max,nmo_max,
                                            SPComplexType(1.0),Aarray.data(),nmo_max,
                                                                          Barray.data(),nmo_max,
                                            SPComplexType(0.0),Carray.data(),nocc_max*nchol_max,
                                                                          Aarray.size());
                Timers[Timer_E3]->stop();

                Timers[Timer_E4]->start();
                copy_n(scl_factors.data(),scl_factors.size(),dev_scl_factors.origin());
                using ma::batched_dot_wabn_wban; 
                batched_dot_wabn_wban(scl_factors.size(),nwalk,nocc_max,nchol_max,
                                         dev_scl_factors.origin(),   
                                         T1.origin(),   
                                         to_address(E[0].origin())+1,E.stride(0));
                Timers[Timer_E4]->stop();

                Timers[Timer_E5]->start();
                if(addEJ) {
                    int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
                    copy_n(kdiag.data(),kdiag.size(),IMats.origin());
                    using ma::batched_Tab_to_Klr;
                    batched_Tab_to_Klr(kdiag.size(),nwalk,nocc_max,nchol_max,nchol_tot,
                                        ncholpQ[Q],nc0,
                                        IMats.origin(),
                                        T1.origin(),
                                        Kl.origin(),
                                        Kr.origin());
                }    
                Timers[Timer_E5]->stop();

                // reset
                Aarray.clear();
                Barray.clear();
                Carray.clear();
                scl_factors.clear();
                kdiag.clear();
                batch_cnt=0;
              }
            }
          }  

          if( batch_cnt > 0 ) {
            Timers[Timer_E3]->start();
            gemmBatched('T','N',nocc_max*nchol_max,nwalk*nocc_max,nmo_max,
                                            SPComplexType(1.0),Aarray.data(),nmo_max,
                                                                          Barray.data(),nmo_max,
                                            SPComplexType(0.0),Carray.data(),nocc_max*nchol_max,
                                                                          Aarray.size());
            Timers[Timer_E3]->stop();

            Timers[Timer_E4]->start();
            copy_n(scl_factors.data(),scl_factors.size(),dev_scl_factors.origin());
            using ma::batched_dot_wabn_wban;
            batched_dot_wabn_wban(scl_factors.size(),nwalk,nocc_max,nchol_max,
                                         dev_scl_factors.origin(),
                                         T1.origin(),
                                         to_address(E[0].origin())+1,E.stride(0));
            Timers[Timer_E4]->stop();

            Timers[Timer_E5]->start();
            if(addEJ) {
                int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
                copy_n(kdiag.data(),kdiag.size(),IMats.origin());
                using ma::batched_Tab_to_Klr;
                batched_Tab_to_Klr(kdiag.size(),nwalk,nocc_max,nchol_max,nchol_tot,
                                        ncholpQ[Q],nc0,
                                        IMats.origin(),
                                        T1.origin(),
                                        Kl.origin(),
                                        Kr.origin());
            }
            Timers[Timer_E5]->stop();
          }
        } // Q
        if(walker_type==COLLINEAR) 
          APP_ABORT("Error: Not implemented.\n");
      }  

      if(addEJ) {
        if(not addEXX) {
          // calculate Kr
          APP_ABORT(" Error: Finish addEJ and not addEXX");
        }
        Timers[Timer_E6]->start();
        size_t nqk=0;  
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        using ma::adotpby;
        for(int n=0; n<nwalk; ++n) {
          adotpby(SPComplexType(0.5*scl*scl),Kl[n],Kr[n],
                      ComplexType(0.0),E[n].origin()+2);
        }
        Timers[Timer_E6]->stop();
      }
    }

    // KEleft and KEright must be in shared memory for this to work correctly  
    template<class Mat, class MatB, class MatC, class MatD>
    void energy_sampleQ(Mat&& E, MatB const& Gc, int nd, MatC* KEleft, MatD* KEright, bool addH1=true, bool addEJ=true, bool addEXX=true) {

      APP_ABORT("Error: energy_sampleQ not yet implemented in batched routine.\n");
      using std::fill_n;
      int nkpts = nopk.size(); 
      assert(E.size(1)>=3);
      assert(nd >= 0 && nd < nelpk.size());  

      int nwalk = Gc.size(1);
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nocca_tot = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+nkpts,0);
      int nocca_max = *std::max_element(nelpk[nd].begin(),nelpk[nd].begin()+nkpts);
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      int noccb_tot = 0;
      if(walker_type==COLLINEAR) noccb_tot = std::accumulate(nelpk[nd].begin()+nkpts,
                                      nelpk[nd].begin()+2*nkpts,0);
      int getKr = KEright!=nullptr;
      int getKl = KEleft!=nullptr;
      if(E.size(0) != nwalk || E.size(1) < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix dimensions \n");

      size_t mem_needs(nwalk*nkpts*nkpts*nspin*nocca_max*nmo_max);
      size_t cnt(0);  
      if(addEJ) { 
        if(not getKr) mem_needs += nwalk*local_nCV;
        if(not getKl) mem_needs += nwalk*local_nCV;
      }
      set_shm_buffer(mem_needs);

      // messy
      sp_pointer_shared Krptr, Klptr;
      size_t Knr=0, Knc=0;
      if(addEJ) {
        Knr=nwalk;
        Knc=local_nCV;
        cnt=0;
        if(getKr) {
          assert(KEright->size(0) == nwalk && KEright->size(1) == local_nCV);
          assert(KEright->stride() == KEright->size(1));
          Krptr = make_device_ptr(KEright->origin());
        } else {
          Krptr = SM_TMats.origin();
          cnt += nwalk*local_nCV;
        }
        if(getKl) {
          assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
          assert(KEleft->stride(0) == KEleft->size(1));
          Klptr = make_device_ptr(KEleft->origin());
        } else {
          Klptr = SM_TMats.origin()+cnt;
          cnt += nwalk*local_nCV;
        }
        fill_n(Krptr,Knr*Knc,SPComplexType(0.0));
        fill_n(Klptr,Knr*Knc,SPComplexType(0.0));
      } else if(getKr or getKl) {
        APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
      }
      SpMatrix_ref Kl(Klptr,{Knr,Knc});
      SpMatrix_ref Kr(Krptr,{Knr,Knc});

      for(int n=0; n<nwalk; n++) 
        fill_n(E[n].origin(),3,ComplexType(0.));

      assert(Gc.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      C3Tensor_cref G3Da(make_device_ptr(Gc.origin()),{nocca_tot,nmo_tot,nwalk} );
      C3Tensor_cref G3Db(make_device_ptr(Gc.origin())+G3Da.num_elements()*(nspin-1),
                            {noccb_tot,nmo_tot,nwalk} );

      Sp4Tensor_ref GKK(SM_TMats.origin()+cnt,
                        {nspin,nkpts,nkpts,nwalk*nmo_max*nocca_max});
      cnt+=GKK.num_elements();
      GKaKjw_to_GKKwaj(G3Da,GKK[0],nelpk[nd].sliced(0,nkpts),dev_nelpk[nd],dev_a0pk[nd]);
      if(walker_type==COLLINEAR)  
        GKaKjw_to_GKKwaj(G3Db,GKK[1],nelpk[nd].sliced(nkpts,2*nkpts),
                                     dev_nelpk[nd].sliced(nkpts,2*nkpts),
                                     dev_a0pk[nd].sliced(nkpts,2*nkpts));

      // one-body contribution
      // haj[ndet*nkpts][nocc*nmo]
      // not parallelized for now, since it would require customization of Wfn 
      if(addH1) {
        // must use Gc since GKK is is SP
        int na=0, nk=0, nb=0;
        for(int n=0; n<nwalk; n++)
          E[n][0] = E0;  
        for(int K=0; K<nkpts; ++K) {
#if defined(MIXED_PRECISION) 
          CMatrix_ref haj_K(haj[nd*nkpts+K].origin(),{nocc_max,nmo_max});
          for(int a=0; a<nelpk[nd][K]; ++a)
            ma::product(ComplexType(1.),ma::T(G3Da[na+a].sliced(nk,nk+nopk[K])),
                                        haj_K[a].sliced(0,nopk[K]),
                        ComplexType(1.),E({0,nwalk},0));
          na+=nelpk[nd][K];
          if(walker_type==COLLINEAR) {
            boost::multi::array_ref<ComplexType,2> haj_Kb(haj_K.origin()+haj_K.num_elements(),
                                                      {nocc_max,nmo_max});
            for(int b=0; b<nelpk[nd][nkpts+K]; ++b)
              ma::product(ComplexType(1.),ma::T(G3Db[nb+b].sliced(nk,nk+nopk[K])),
                                        haj_Kb[b].sliced(0,nopk[K]),
                        ComplexType(1.),E({0,nwalk},0));
            nb+=nelpk[nd][nkpts+K];
          }
          nk+=nopk[K];
#else
          nk = nopk[K];
          {
            na = nelpk[nd][K];
            CVector_ref haj_K(haj[nd*nkpts+K].origin(),{nocc_max*nmo_max});
            SpMatrix_ref Gaj(GKK[0][K][K].origin(),{nwalk,nocc_max*nmo_max});
            ma::product(ComplexType(1.),Gaj,haj_K,ComplexType(1.),E({0,nwalk},0));
          }
          if(walker_type==COLLINEAR) {
            na = nelpk[nd][nkpts+K];
            CVector_ref haj_K(haj[nd*nkpts+K].origin()+nocc_max*nmo_max,{nocc_max*nmo_max});
            SpMatrix_ref Gaj(GKK[1][K][K].origin(),{nwalk,nocc_max*nmo_max});
            ma::product(ComplexType(1.),Gaj,haj_K,ComplexType(1.),E({0,nwalk},0));
          }
#endif
        }
      }

      // move calculation of H1 here	
      // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
      //       Not sure how to do it for COLLINEAR.
      if(addEXX) {  

        if(Qwn.size(0) != nwalk || Qwn.size(1) != nsampleQ)
          Qwn.reextent({nwalk,nsampleQ});
        {
          for(int n=0; n<nwalk; ++n) 
            for(int nQ=0; nQ<nsampleQ; ++nQ) {
              Qwn[n][nQ] = distribution(generator);
/*
              RealType drand = distribution(generator);
              RealType s(0.0);
              bool found=false;
              for(int Q=0; Q<nkpts; Q++) {
                s += gQ[Q];
                if( drand < s ) {
                  Qwn[n][nQ] = Q;
                  found=true;
                  break;
                }
              } 
              if(not found) 
                APP_ABORT(" Error: sampleQ Qwn. \n");  
*/
            }
        }

        size_t local_memory_needs = 2*nocca_max*nocca_max*nchol_max; 
        if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});
        size_t local_cnt=0; 
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        size_t nqk=1;  
        for(int n=0; n<nwalk; ++n) {
          for(int nQ=0; nQ<nsampleQ; ++nQ) {
            int Q = Qwn[n][nQ];
            for(int Ka=0; Ka<nkpts; ++Ka) {
              for(int Kb=0; Kb<nkpts; ++Kb) {
                { 
                  int nchol = ncholpQ[Q];
                  int Qm = kminus[Q];
                  int Qm_ = (Q==Q0?nkpts:Qm);
                  int Kl = QKToK2[Qm][Kb];
                  int Kk = QKToK2[Q][Ka];
                  int nl = nopk[Kl];
                  int nb = nelpk[nd][Kb];
                  int na = nelpk[nd][Ka];
                  int nk = nopk[Kk];

                  SpMatrix_ref Gal(GKK[0][Ka][Kl].origin()+n*na*nl,{na,nl});
                  SpMatrix_ref Gbk(GKK[0][Kb][Kk].origin()+n*nb*nk,{nb,nk});
                  SpMatrix_ref Lank(LQKank[nd*nspin*(nkpts+1)+Q][Ka].origin(),
                                                 {na*nchol,nk});
                  SpMatrix_ref Lbnl(LQKank[nd*nspin*(nkpts+1)+Qm_][Kb].origin(),
                                                 {nb*nchol,nl});

                  SpMatrix_ref Tban(TMats.origin()+local_cnt,{nb,na*nchol});
                  Sp3Tensor_ref T3Dban(TMats.origin()+local_cnt,{nb,na,nchol});
                  SpMatrix_ref Tabn(Tban.origin()+Tban.num_elements(),{na,nb*nchol});
                  Sp3Tensor_ref T3Dabn(Tban.origin()+Tban.num_elements(),{na,nb,nchol});

                  ma::product(Gal,ma::T(Lbnl),Tabn);
                  ma::product(Gbk,ma::T(Lank),Tban);

                  ComplexType E_(0.0);
                  for(int a=0; a<na; ++a)
                    for(int b=0; b<nb; ++b)
                      E_ += ma::dot(T3Dabn[a][b],T3Dban[b][a]);
                  E[n][1] -= scl*0.5*E_/gQ[Q]/double(nsampleQ);

                } // if

                if(walker_type==COLLINEAR) {

                  { 
                    int nchol = ncholpQ[Q];
                    int Qm = kminus[Q];
                    int Qm_ = (Q==Q0?nkpts:Qm);
                    int Kl = QKToK2[Qm][Kb];
                    int Kk = QKToK2[Q][Ka];
                    int nl = nopk[Kl];
                    int nb = nelpk[nd][nkpts+Kb];
                    int na = nelpk[nd][nkpts+Ka];
                    int nk = nopk[Kk];

                    SpMatrix_ref Gal(GKK[1][Ka][Kl].origin()+n*na*nl,{na,nl});
                    SpMatrix_ref Gbk(GKK[1][Kb][Kk].origin()+n*nb*nk,{nb,nk});
                    SpMatrix_ref Lank(LQKank[nd*nspin*(nkpts+1)+nkpts+1+Q][Ka].origin(),
                                                 {na*nchol,nk});
                    SpMatrix_ref Lbnl(LQKank[nd*nspin*(nkpts+1)+nkpts+1+Qm_][Kb].origin(),
                                                 {nb*nchol,nl});

                    SpMatrix_ref Tban(TMats.origin()+local_cnt,{nb,na*nchol});
                    Sp3Tensor_ref T3Dban(TMats.origin()+local_cnt,{nb,na,nchol});
                    SpMatrix_ref Tabn(Tban.origin()+Tban.num_elements(),{na,nb*nchol});
                    Sp3Tensor_ref T3Dabn(Tban.origin()+Tban.num_elements(),{na,nb,nchol});
  
                    ma::product(Gal,ma::T(Lbnl),Tabn);
                    ma::product(Gbk,ma::T(Lank),Tban);
  
                    ComplexType E_(0.0);
                    for(int a=0; a<na; ++a)
                      for(int b=0; b<nb; ++b)
                        E_ += ma::dot(T3Dabn[a][b],T3Dban[b][a]);
                    E[n][1] -= scl*0.5*E_/gQ[Q]/double(nsampleQ);

                  } // if
                } // COLLINEAR 
              } // Kb 
            } // Ka
          } // nQ
        } // n 
      }  

      if(addEJ) {
        size_t local_memory_needs = 2*nchol_max*nwalk; 
        if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});
        cnt=0; 
        SpMatrix_ref Kr_local(TMats.origin(),{nwalk,nchol_max}); 
        cnt+=Kr_local.num_elements();
        SpMatrix_ref Kl_local(TMats.origin()+cnt,{nwalk,nchol_max}); 
        cnt+=Kl_local.num_elements();
        fill_n(Kr_local.origin(),Kr_local.num_elements(),SPComplexType(0.0));
        fill_n(Kl_local.origin(),Kl_local.num_elements(),SPComplexType(0.0));
        size_t nqk=1;  
        for(int Q=0; Q<nkpts; ++Q) {
          bool haveKE=false;
          for(int Ka=0; Ka<nkpts; ++Ka) {
            { 
              haveKE=true;
              int nchol = ncholpQ[Q];
              int Qm = kminus[Q];
              int Qm_ = (Q==Q0?nkpts:Qm);
              int Kl = QKToK2[Qm][Ka];
              int Kk = QKToK2[Q][Ka];
              int nl = nopk[Kl];
              int na = nelpk[nd][Ka];
              int nk = nopk[Kk];

              Sp3Tensor_ref Gwal(GKK[0][Ka][Kl].origin(),{nwalk,na,nl});
              Sp3Tensor_ref Gwbk(GKK[0][Ka][Kk].origin(),{nwalk,na,nk});
              Sp3Tensor_ref Lank(LQKank[nd*nspin*(nkpts+1)+Q][Ka].origin(),
                                                 {na,nchol,nk});
              Sp3Tensor_ref Lbnl(LQKank[nd*nspin*(nkpts+1)+Qm_][Ka].origin(),
                                                 {na,nchol,nl});

              // Twan = sum_l G[w][a][l] L[a][n][l]
              for(int n=0; n<nwalk; ++n) 
                for(int a=0; a<na; ++a)  
                  ma::product(SPComplexType(1.0),Lbnl[a],Gwal[n][a],
                              SPComplexType(1.0),Kl_local[n]);
              for(int n=0; n<nwalk; ++n) 
                for(int a=0; a<na; ++a)  
                  ma::product(SPComplexType(1.0),Lank[a],Gwbk[n][a],
                              SPComplexType(1.0),Kr_local[n]);
            } // if

            if(walker_type==COLLINEAR) {

              { 
                haveKE=true;
                int nchol = ncholpQ[Q];
                int Qm = kminus[Q];
                int Qm_ = (Q==Q0?nkpts:Qm);
                int Kl = QKToK2[Qm][Ka];
                int Kk = QKToK2[Q][Ka];
                int nl = nopk[Kl];
                int na = nelpk[nd][nkpts+Ka];
                int nk = nopk[Kk];

                Sp3Tensor_ref Gwal(GKK[1][Ka][Kl].origin(),{nwalk,na,nl});
                Sp3Tensor_ref Gwbk(GKK[1][Ka][Kk].origin(),{nwalk,na,nk});
                Sp3Tensor_ref Lank(LQKank[nd*nspin*(nkpts+1)+nkpts+1+Q][Ka].origin(),
                                                 {na,nchol,nk});
                Sp3Tensor_ref Lbnl(LQKank[nd*nspin*(nkpts+1)+nkpts+1+Qm_][Ka].origin(),
                                                 {na,nchol,nl});

                // Twan = sum_l G[w][a][l] L[a][n][l]
                for(int n=0; n<nwalk; ++n)
                  for(int a=0; a<na; ++a)  
                    ma::product(SPComplexType(1.0),Lbnl[a],Gwal[n][a],
                                SPComplexType(1.0),Kl_local[n]);
                for(int n=0; n<nwalk; ++n)
                  for(int a=0; a<na; ++a)  
                    ma::product(SPComplexType(1.0),Lank[a],Gwbk[n][a],
                                SPComplexType(1.0),Kr_local[n]);

              } // if
            } // COLLINEAR
          } // Ka
          if(haveKE) {
            int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);  
            using ma::axpy;
            for(int n=0; n<nwalk; n++) {
              axpy(SPComplexType(1.0),Kr_local[n].sliced(0,ncholpQ[Q]),
                                        Kr[n].sliced(nc0,nc0+ncholpQ[Q])); 
              axpy(SPComplexType(1.0),Kl_local[n].sliced(0,ncholpQ[Q]),
                                        Kl[n].sliced(nc0,nc0+ncholpQ[Q])); 
            }
          } // to release the lock
          if(haveKE) { 
            fill_n(Kr_local.origin(),Kr_local.num_elements(),SPComplexType(0.0));
            fill_n(Kl_local.origin(),Kl_local.num_elements(),SPComplexType(0.0));
          }  
        } // Q
        nqk=0;  
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        for(int n=0; n<nwalk; ++n) {
          for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
            {
              int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
              E[n][2] += 0.5*scl*scl*ma::dot(Kl[n]({nc0,nc0+ncholpQ[Q]}),
                                            Kr[n]({nc0,nc0+ncholpQ[Q]}));  
            }
          }
        }
      }
    }

    template<class... Args>
    void fast_energy(Args&&... args)
    {
      APP_ABORT(" Error: fast_energy not implemented in KP3IndexFactorization_batched. \n"); 
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatA>())),
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>())),
             typename = void
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      using BType = typename std::decay<MatB>::type::element ;
      using AType = typename std::decay<MatA>::type::element ;
      boost::multi::array_ref<BType,2,decltype(v.origin())> v_(v.origin(),{v.size(0),1});
      boost::multi::array_ref<AType,2,decltype(X.origin())> X_(X.origin(),{X.size(0),1});
      return vHS(X_,v_,a,c);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatA>())),
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>()))
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {

      int nkpts = nopk.size();
      int nwalk = X.size(1);
      assert(v.size(0)==nwalk);
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      assert(X.num_elements() == nwalk*2*local_nCV);
      assert(v.num_elements() == nwalk*nmo_tot*nmo_tot);
      SPComplexType one(1.0,0.0);
      SPComplexType im(0.0,1.0);
      size_t mem_needs((nkpts+1)*nkpts*nwalk*nmo_max*nmo_max + nwalk*2*nkpts*nchol_max);
      if(SM_TMats.num_elements() < mem_needs) SM_TMats.reextent({mem_needs,1});

      Sp3Tensor_ref vKK(SM_TMats.origin(),{nkpts+1,nkpts,nwalk*nmo_max*nmo_max} );
      Sp4Tensor_ref XQnw(SM_TMats.origin()+vKK.num_elements(),{nkpts,2,nchol_max,nwalk} );
      fill_n(XQnw.origin(),XQnw.num_elements(),SPComplexType(0.0));

      Timers[Timer_vHS1]->start();
      // "rotate" X  
      //  XIJ = 0.5*a*(Xn+ -i*Xn-), XJI = 0.5*a*(Xn+ +i*Xn-)  
// transform to single precision
      SpMatrix_ref Xdev(make_device_ptr(X.origin()),X.extensions());
      for(int Q=0, nq=0; Q<nkpts; ++Q) { 
        if(Q != kminus[Q] || Q==Q0) { 
          auto&& Xp(Xdev.sliced(nq,nq+ncholpQ[Q]));
          auto&& Xm(Xdev.sliced(nq+ncholpQ[Q],nq+2*ncholpQ[Q]));
          ma::add(ComplexType(0.5*a),Xp,ComplexType(-0.5*a*im),Xm,XQnw[Q][0].sliced(0,ncholpQ[Q]));
          ma::add(ComplexType(0.5*a),Xp,ComplexType(0.5*a*im),Xm,XQnw[Q][1].sliced(0,ncholpQ[Q]));
        } else {
          ma::axpy(ComplexType(a),Xdev.sliced(nq,nq+ncholpQ[Q]),XQnw[Q][0].sliced(0,ncholpQ[Q]));
        }  
        nq+=2*ncholpQ[Q];
      } 
      Timers[Timer_vHS1]->stop();
      Timers[Timer_vHS2]->start();
      //  then combine Q/(-Q) pieces
      //  X(Q)np = (X(Q)np + X(-Q)nm)
      for(int Q=0; Q<nkpts; ++Q) {
        if(Q != kminus[Q]) {
          int Qm = kminus[Q];
          ma::axpy(ComplexType(1.0),XQnw[Qm][1],XQnw[Q][0]);
        }
      }
      {
        // assuming contiguous
        ma::scal(c,v);
      }
      Timers[Timer_vHS2]->stop();
        
      Timers[Timer_vHS3]->start();
      using ma::gemmBatched;
      std::vector<sp_pointer> Aarray;
      std::vector<sp_pointer> Barray;
      std::vector<sp_pointer> Carray;
      Aarray.reserve(nkpts*nkpts);
      Barray.reserve(nkpts*nkpts);
      Carray.reserve(nkpts*nkpts);
      for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
        // v[nw][i(in K)][k(in Q(K))] += sum_n LQK[i][k][n] X[Q][0][n][nw]
        if(Q < kminus[Q] || Q==Q0) {
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            int QK = QKToK2[Q][K];
            Aarray.push_back(LQKikn[Q][K].origin());
            Barray.push_back(XQnw[Q][0].origin());
            Carray.push_back(vKK[K][QK].origin());
          }
        } else if(Q > kminus[Q]) { // use L(-Q)(ki)*
        } else { // Q==(-Q) and Q!=Q0
          int kpos(0);
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            if(K < QKToK2[Q][K]) { 
              int QK = QKToK2[Q][K];
              Aarray.push_back(LQKikn[Q][kpos].origin());
              Barray.push_back(XQnw[Q][0].origin());
              Carray.push_back(vKK[K][QK].origin());
              kpos++;
            }
          }
        }
      }
      int nmo_max2 = nmo_max*nmo_max;  
      // C: v = T(X) * T(Lik) --> F: T(Lik) * T(X) = v   
      gemmBatched('T','T',nmo_max2,nwalk,nchol_max,SPComplexType(1.0),Aarray.data(),nchol_max,Barray.data(),nwalk,
                                                   SPComplexType(0.0),Carray.data(),nmo_max2,Aarray.size());

      Aarray.clear();  
      Barray.clear();  
      Carray.clear();  
      for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
        // v[nw][i(in K)][k(in Q(K))] += sum_n LQK[i][k][n] X[Q][0][n][nw]
        if(Q < kminus[Q] || Q==Q0) {
        } else if(Q > kminus[Q]) { // use L(-Q)(ki)*
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            int QK = QKToK2[Q][K];
            Aarray.push_back(LQKikn[kminus[Q]][QK].origin());
            Barray.push_back(XQnw[Q][0].origin());
            Carray.push_back(vKK[K][QK].origin());
          }
        } else { // Q==(-Q) and Q!=Q0
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
            if(K > QKToK2[Q][K]) {
              int kpos(0);
              // find (K,QK)
              for(int q=0, qmin=QKToK2[Q][K]; q<qmin; q++)
                if(q < QKToK2[Q][q]) kpos++;
              int QK = QKToK2[Q][K];
              Aarray.push_back(LQKikn[Q][kpos].origin());
              Barray.push_back(XQnw[Q][0].origin());
              Carray.push_back(vKK[K][QK].origin());
            }
          }
        }
      }
      if(Q0>=0) {
        for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
          Aarray.push_back(LQKikn[Q0][K].origin());
          Barray.push_back(XQnw[Q0][1].origin());
          Carray.push_back(vKK[nkpts][K].origin());
        }
      }
      // C: v = T(X) * T(Lik) --> F: T(Lik) * T(X) = v   
      gemmBatched('C','T',nmo_max2,nwalk,nchol_max,SPComplexType(1.0),Aarray.data(),nchol_max,Barray.data(),nwalk,
                                                   SPComplexType(0.0),Carray.data(),nmo_max2,Aarray.size());

      Timers[Timer_vHS3]->stop();
      Timers[Timer_vHS4]->start();
      Sp3Tensor_ref v3D(make_device_ptr(v.origin()),{nwalk,nmo_tot,nmo_tot});
      vKKwij_to_vwKiKj(vKK,v3D);
      Timers[Timer_vHS4]->stop();
      // do I need to "rotate" back, can be done if necessary
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatA&>())),
//             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>())),
             typename = void
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
      using BType = typename std::decay<MatB>::type::element ;
      using AType = typename std::decay<MatA>::type::element ;
      boost::multi::array_ref<BType,2,decltype(v.origin())> v_(v.origin(),
                                        {v.size(0),1});
      boost::multi::array_ref<AType const,2,decltype(G.origin())> G_(G.origin(),
                                        {G.size(0),1});
      return vbias(G_,v_,a,c,k);  
    }

/*
    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
//             typename = typename std::enable_if_t<(std::is_convertible<typename std::decay<MatA>::type::pointer,pointer>::value)>,  
//             typename = typename std::enable_if_t<(not std::is_convertible<typename std::decay<MatB>::type::pointer,pointer>::value)>,  
              typename = void,  
              typename = void,  
              typename = void,  
              typename = void,  
              typename = void  
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int nd=0) {
    }
*/

    template<class MatA, class MatB,
             typename = std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
//             typename = std::enable_if_t<(std::is_convertible<typename std::decay<MatA>::type::element_ptr,pointer>::value)>, 
//             typename = std::enable_if_t<(std::is_convertible<typename std::decay<MatB>::type::element_ptr,pointer>::value)> 
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int nd=0) {

      using ma::gemmBatched;

      int nkpts = nopk.size();
      assert(nd >= 0 && nd < nelpk.size());
      int nwalk = G.size(1);
      assert(v.size(0)==2*local_nCV);  
      assert(v.size(1)==nwalk);  
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nocca_tot = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+nkpts,0);
      int nocca_max = *std::max_element(nelpk[nd].begin(),nelpk[nd].begin()+nkpts);
      int noccb_max = nocca_max; 
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      int noccb_tot = 0;
      if(walker_type==COLLINEAR) {
        noccb_tot = std::accumulate(nelpk[nd].begin()+nkpts,
                                    nelpk[nd].begin()+2*nkpts,0);
        noccb_max = *std::max_element(nelpk[nd].begin()+nkpts,
                                      nelpk[nd].begin()+2*nkpts);
      }
      RealType scl = (walker_type==CLOSED?2.0:1.0);
      SPComplexType one(1.0,0.0);
      SPComplexType halfa(0.5*a*scl,0.0);
      SPComplexType minusimhalfa(0.0,-0.5*a*scl);
      SPComplexType imhalfa(0.0,0.5*a*scl);

      // space for GQK and for v1
      size_t mem_needs(nkpts*nkpts*nwalk*nocca_max*nmo_max + (nkpts+1)*nchol_max*nwalk); 
      if(SM_TMats.num_elements() < mem_needs) SM_TMats.reextent({mem_needs,1});

      assert(G.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      C3Tensor_cref G3Da(make_device_ptr(G.origin()),{nocca_tot,nmo_tot,nwalk} );
      C3Tensor_cref G3Db(make_device_ptr(G.origin())+G3Da.num_elements()*(nspin-1),
                                                     {noccb_tot,nmo_tot,nwalk} );
      {  
        // assuming contiguous
        ma::scal(c,v);
      }

      for(int spin=0; spin<nspin; spin++) {
        Sp3Tensor_ref v1(SM_TMats.origin(),{nkpts+1,nchol_max,nwalk});
        Sp3Tensor_ref GQ(v1.origin()+v1.num_elements(),{nkpts,nkpts*nocc_max*nmo_max,nwalk} );
        fill_n(GQ.origin(),GQ.num_elements(),SPComplexType(0.0));

        Timers[Timer_vbias1]->start();
        if(spin==0) 
          GKaKjw_to_GQKajw(G3Da,GQ,nelpk[nd],dev_nelpk[nd],dev_a0pk[nd]);
        else
          GKaKjw_to_GQKajw(G3Db,GQ,nelpk[nd].sliced(nkpts,2*nkpts),
                                   dev_nelpk[nd].sliced(nkpts,2*nkpts),
                                   dev_a0pk[nd].sliced(nkpts,2*nkpts));
        Timers[Timer_vbias1]->stop();

        Timers[Timer_vbias2]->start();
        // can use productStridedBatched if LQKakn is changed to a 3Tensor array
        int Kak = nkpts*nocc_max*nmo_max;
        std::vector<sp_pointer> Aarray;
        std::vector<sp_pointer> Barray;
        std::vector<sp_pointer> Carray;
        Aarray.reserve(nkpts+1);
        Barray.reserve(nkpts+1);
        Carray.reserve(nkpts+1);
        for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
          // v_[Q][n][w] = sum_Kak LQ[Kak][n]*G[Q][Kak][w]
          //             F: -->   G[Kak][w] * LQ[Kak][n]
          Aarray.push_back(GQ[Q].origin());
          Barray.push_back(LQKakn[nd*nspin*(nkpts+1)+spin*(nkpts+1)+Q].origin());
          Carray.push_back(v1[Q].origin());
        }
        if(Q0 >= 0) {
          Aarray.push_back(GQ[Q0].origin());
          Barray.push_back(LQKakn[nd*nspin*(nkpts+1)+spin*(nkpts+1)+nkpts].origin());
          Carray.push_back(v1[nkpts].origin());
        }
        gemmBatched('N','T',nwalk,nchol_max,Kak,SPComplexType(1.0),Aarray.data(),nwalk,
                                                                   Barray.data(),nchol_max,
                                                SPComplexType(0.0),Carray.data(),nwalk,
                                                Aarray.size());
        Timers[Timer_vbias2]->stop();


        Timers[Timer_vbias3]->start();
        vbias_from_v1(halfa,v1,v);
        Timers[Timer_vbias3]->stop();
      }
    }

    bool distribution_over_cholesky_vectors() const{ return true; }
    int number_of_ke_vectors() const{ return std::accumulate(ncholpQ.begin(),ncholpQ.end(),0); }
    int local_number_of_cholesky_vectors() const{ return 2*std::accumulate(ncholpQ.begin(),ncholpQ.end(),0); } 
    int global_number_of_cholesky_vectors() const{ return global_nCV; }

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const{return false;}
    bool transposed_G_for_E() const{return false;} 
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const{return true;} 

    bool fast_ph_energy() const { return false; }

  private:

    int Q0;

    int nocc_max;
    size_t batch_size;

    Allocator allocator_;
    SpAllocator sp_allocator_;
    Allocator_shared allocator_shared_;
    SpAllocator_shared sp_allocator_shared_;
    IAllocator_shared iallocator_shared_;

    WALKER_TYPES walker_type;

    int global_nCV;
    int local_nCV;

    ValueType E0;

    // bare one body hamiltonian
    mpi3C3Tensor H1;

    // (potentially half rotated) one body hamiltonian
    shmCMatrix haj;
    //std::vector<shmCVector> haj;

    // number of orbitals per k-point
    boost::multi::array<int,1> nopk;

    // number of cholesky vectors per Q-point
    boost::multi::array<int,1> ncholpQ;

    // position of (-K) in kp-list for every K 
    boost::multi::array<int,1> kminus;

    // number of electrons per k-point
    // nelpk[ndet][nspin*nkpts]
    //shmIMatrix nelpk;
    boost::multi::array<int,2> nelpk;

    // maps (Q,K) --> k2
    //shmIMatrix QKToK2; 
    boost::multi::array<int,2> QKToK2; 

    //Cholesky Tensor Lik[Q][nk][i][k][n]
    std::vector<shmSpMatrix> LQKikn;

    // half-tranformed Cholesky tensor
    std::vector<shmSpMatrix> LQKank;

    // half-tranformed Cholesky tensor
    std::vector<shmSpMatrix> LQKakn;

    // one-body piece of Hamiltonian factorization
    mpi3C3Tensor vn0;

    int nsampleQ;
    std::vector<RealType> gQ;
    boost::multi::array<int,2> Qwn;
    std::default_random_engine generator;
    std::discrete_distribution<int> distribution;

    // shared buffer space
    // using matrix since there are issues with vectors
    shmSpMatrix SM_TMats;
    SpMatrix TMats;
    IMatrix IMats;

    BoolMatrix KKTransID;
    IVector dev_nopk;
    IVector dev_i0pk;
    IVector dev_kminus;
    IVector dev_ncholpQ;
    IVector dev_ncholpQ0;
    IMatrix dev_nelpk;
    IMatrix dev_a0pk;
    IMatrix dev_QKToK2;

//    std::vector<std::unique_ptr<shared_mutex>> mutex;

//    boost::multi::array<ComplexType,3> Qave;
//    int cntQave=0;
    std::vector<ComplexType> EQ;
//    std::default_random_engine generator;
//    std::uniform_real_distribution<RealType> distribution(RealType(0.0),Realtype(1.0));

    void set_shm_buffer(size_t N) {
      if(SM_TMats.num_elements() < N) 
        SM_TMats.reextent({N,1});
    }

    template<class MatA, class MatB, class IVec, class IVec2>
    void GKaKjw_to_GKKwaj(MatA const& GKaKj, MatB && GKKaj,IVec && nocc, IVec2 && dev_no, IVec2 && dev_a0)
    {
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
//      int nocc_max = *std::max_element(nocc.begin(),nocc.end());
      int nmo_tot = GKaKj.size(1);
      int nwalk = GKaKj.size(2);
      int nkpts = nopk.size();
      assert(GKKaj.num_elements() >= nkpts*nkpts*nwalk*nocc_max*nmo_max);

      using ma::KaKjw_to_KKwaj;  
      KaKjw_to_KKwaj(nwalk,nkpts,nmo_max,nmo_tot,nocc_max,
                                dev_nopk.origin(),dev_i0pk.origin(),
                                dev_no.origin(),dev_a0.origin(),
                                GKaKj.origin(),GKKaj.origin());
    }

    template<class MatA, class MatB, class IVec, class IVec2>
    void GKaKjw_to_GQKajw(MatA const& GKaKj, MatB && GQKaj, IVec && nocc, IVec2 && dev_no, IVec2 && dev_a0)
    {
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
//      int nocc_max = *std::max_element(nocc.begin(),nocc.end());
      int nmo_tot = GKaKj.size(1);
      int nwalk = GKaKj.size(2);
      int nkpts = nopk.size();
      assert(GQKaj.num_elements() >= nkpts*nkpts*nwalk*nocc_max*nmo_max);

      using ma::KaKjw_to_QKajw;
      KaKjw_to_QKajw(nwalk,nkpts,nmo_max,nmo_tot,nocc_max,
                                dev_nopk.origin(),dev_i0pk.origin(),
                                dev_no.origin(),dev_a0.origin(),
                                dev_QKToK2.origin(),
                                GKaKj.origin(),GQKaj.origin());  
    }


    /*
     *   vKiKj({nwalk,nmo_tot,nmo_tot});
     *   vKK({nkpts,nkpts,nwalk*nmo_max*nmo_max} );
     */   
    template<class MatA, class MatB>
    void vKKwij_to_vwKiKj(MatA const& vKK, MatB && vKiKj)
    {
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nwalk = vKiKj.size(0);
      int nmo_tot = vKiKj.size(1);
      int nkpts = nopk.size();

      using ma::vKKwij_to_vwKiKj;
      vKKwij_to_vwKiKj(nwalk,nkpts,nmo_max,nmo_tot,KKTransID.origin(),
                                dev_nopk.origin(),dev_i0pk.origin(),
                                vKK.origin(),vKiKj.origin());  
    }

    template<class MatA, class MatB>
    void vbias_from_v1(ComplexType a, MatA const& v1, MatB && vbias)
    {
      int nwalk = vbias.size(1);
      int nkpts = nopk.size();
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());

      using ma::vbias_from_v1;
// using make_device_ptr(vbias.origin()) to catch errors here
      vbias_from_v1(nwalk,nkpts,nchol_max,Q0,dev_kminus.origin(),
                             dev_ncholpQ.origin(),dev_ncholpQ0.origin(),
                             a,v1.origin(),to_address(make_device_ptr(vbias.origin())));
    }

    enum THCTimers
    {
      Timer_vbias1,
      Timer_vbias2,
      Timer_vbias3,
      Timer_vHS1,
      Timer_vHS2,
      Timer_vHS3,
      Timer_vHS4,
      Timer_E1,
      Timer_E2,
      Timer_E3,
      Timer_E4,
      Timer_E5,
      Timer_E6
    };

    TimerNameList_t<THCTimers> THCTimerNames = {
      {Timer_vbias1, "vbias_1"},
      {Timer_vbias2, "vbias_2"},
      {Timer_vbias3, "vbias_3"},
      {Timer_vHS1, "vHS_rotateX"},
      {Timer_vHS2, "vHS_combineX"},
      {Timer_vHS3, "vHS_LX"},
      {Timer_vHS4, "vHS_trans"},
      {Timer_E1, "E_trans"},
      {Timer_E2, "E_H1"},
      {Timer_E3, "E_LG"},
      {Timer_E4, "E_dot"},
      {Timer_E5, "E_Klr"},
      {Timer_E6, "E_EJ"},
    };
    TimerList_t Timers;

};

}

}

#endif
