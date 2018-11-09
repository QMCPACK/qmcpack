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
#include<mutex>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "alf/boost/mpi3/shared_communicator.hpp"
#include "AFQMC/multi/array.hpp"
#include "AFQMC/multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "alf/boost/mpi3/shm/mutex.hpp"

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
  using Sp3Tensor_ref = boost::multi::array_ref<SPComplexType,3>;
  using shmCVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;  
  using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;  
  using shmIMatrix = boost::multi::array<int,2,shared_allocator<int>>;  
  using shmC3Tensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;  
  using shmSpVector = boost::multi::array<SPComplexType,1,shared_allocator<SPComplexType>>;  
  using shmSpMatrix = boost::multi::array<SPComplexType,2,shared_allocator<SPComplexType>>;  
  using shmSp3Tensor = boost::multi::array<SPComplexType,3,shared_allocator<SPComplexType>>;  
  using communicator = boost::mpi3::shared_communicator;
  using shared_mutex = boost::mpi3::shm::mutex;
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
        H1(std::move(hij_)),
        haj(std::move(h1)),
        nopk(std::move(nopk_)),
        ncholpQ(std::move(ncholpQ_)),
        nelpk(std::move(nelpk_)),
        QKToK2(std::move(QKToK2_)),
        LQKikn(std::move(vik)),
        LQKank(std::move(vak)),
        LQKlnb(std::move(vlb)),
        vn0(std::move(vn0_)),
        SM_TMats({1,1},shared_allocator<SPComplexType>{c_}),
        TMats({1,1}),
        mutex(0)
    {
        local_nCV = std::accumulate(ncholpQ.begin(),ncholpQ.end(),0); 
        mutex.reserve(ncholpQ.size());
        for(int nQ=0; nQ<ncholpQ.size(); nQ++)
            mutex.emplace_back(std::make_unique<shared_mutex>(*comm));
    }

    ~KP3IndexFactorization() {}
    
    KP3IndexFactorization(const KP3IndexFactorization& other) = delete;
    KP3IndexFactorization& operator=(const KP3IndexFactorization& other) = delete;
    KP3IndexFactorization(KP3IndexFactorization&& other) = default; 
    KP3IndexFactorization& operator=(KP3IndexFactorization&& other) = default;

    boost::multi_array<ComplexType,2> getOneBodyPropagatorMatrix(TaskGroup_& TG, boost::multi_array<ComplexType,1> const& vMF) {

      int nkpts = nopk.size(); 
      int NMO = std::accumulate(nopk.begin(),nopk.end(),0);
      // in non-collinear case with SO, keep SO matrix here and add it
      // for now, stay collinear
      boost::multi_array<ComplexType,2> P1(extents[NMO][NMO]);

      // making a copy of vMF since it will be modified
      set_shm_buffer(vMF.shape()[0]);
      boost::multi_array_ref<ComplexType,1> vMF_(std::addressof(*SM_TMats.origin()),extents[vMF.shape()[0]]);  

      boost::multi::array_ref<ComplexType,1> P1D(std::addressof(*P1.origin()),{NMO*NMO});
      std::fill_n(P1D.origin(),P1D.num_elements(),ComplexType(0));
      vHS(vMF_, P1D);  
      TG.TG().all_reduce_in_place_n(P1D.origin(),P1D.num_elements(),std::plus<>());

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
    void energy(Mat&& E, MatB const& G, int k, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      MatB* Kr(nullptr);
      MatB* Kl(nullptr);
      energy(E,G,k,Kl,Kr,addH1,addEJ,addEXX);
    }

    // KEleft and KEright must be in shared memory for this to work correctly  
    template<class Mat, class MatB, class MatC, class MatD>
    void energy(Mat&& E, MatB const& Gc, int nd, MatC* KEleft, MatD* KEright, bool addH1=true, bool addEJ=true, bool addEXX=true) {

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
      int getKr = KEright!=nullptr;
      int getKl = KEleft!=nullptr;
      if(E.shape()[0] != nwalk || E.shape()[1] < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix dimensions \n");

double t1, t2, t3, t4=0, t5=0, t6=0;
Timer.reset("T0");
Timer.reset("T1");
Timer.reset("T2");
Timer.reset("T3");

      // messy
      SPComplexType *Krptr, *Klptr; 
      size_t Knr=0, Knc=0;
      if(addEJ) {
        Knr=nwalk;
        Knc=local_nCV;
        size_t mem_needs(0);
        if(not getKr) mem_needs += nwalk*local_nCV;
        if(not getKl) mem_needs += nwalk*local_nCV;
        set_shm_buffer(mem_needs);
        mem_needs=0;
        if(getKr) {
          assert(KEright->shape()[0] == nwalk && KEright->shape()[1] == local_nCV); 
          assert(KEright->strides()[0] == KEright->shape()[1]);
          Krptr = std::addressof(*KEright->origin()); 
        } else {
          Krptr = std::addressof(*SM_TMats.origin()); 
          mem_needs += nwalk*local_nCV;
        }
        if(getKl) {
          assert(KEleft->shape()[0] == nwalk && KEleft->shape()[1] == local_nCV); 
          assert(KEleft->strides()[0] == KEleft->shape()[1]);
          Klptr = std::addressof(*KEleft->origin());
        } else {
          Klptr = std::addressof(*SM_TMats.origin())+mem_needs; 
        }
        if(comm->root()) std::fill_n(Krptr,Knr*Knc,SPComplexType(0.0));
        if(comm->root()) std::fill_n(Klptr,Knr*Knc,SPComplexType(0.0));
      } else if(getKr or getKl) {
        APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
      }   
      SpMatrix_ref Kl(Klptr,{Knr,Knc});
      SpMatrix_ref Kr(Krptr,{Knr,Knc});
      comm->barrier();

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
Timer.start("T0");  
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
Timer.stop("T0"); 
t1 = Timer.total("T0");          

      // move calculation of H1 here	
      // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
      //       Not sure how to do it for COLLINEAR.
Timer.reset("T0");
Timer.start("T0");  
      if(addEXX) {  
        size_t local_memory_needs = nocca_tot*nchol_max*nmo_max + nmo_max*nocca_max + nchol_max;
        if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});

        RealType scl = (walker_type==CLOSED?2.0:1.0);
        size_t nqk=1;  // start count at 1 to "offset" the calcuation of E1 done at root
        // avoiding vectors for now
        SpMatrix_ref Kr_local(TMats.origin(),{1,nchol_max}); 
        std::fill_n(Kr_local.origin(),Kr_local.num_elements(),SPComplexType(0.0));
        for(int n=0; n<nwalk; ++n) {
          for(int Kl=0; Kl<nkpts; ++Kl) {       // K is the index of the kpoint pair of (l,k)
            for(int Q=0; Q<nkpts; ++Q) {            // momentum conservation index   
              bool haveKE=false;
              if((nqk++)%comm->size() == comm->rank()) { 
                int nchol = ncholpQ[Q];
                int nl = nopk[Kl];
                int nl0 = std::accumulate(nopk.begin(),nopk.begin()+Kl,0);
                int nb = nelpk[nd][QKToK2[Q][Kl]];
                int nb0 = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+QKToK2[Q][Kl],0);

                SpMatrix_ref Llnb(std::addressof(*LQKlnb[nd*nkpts+Q][Kl].origin()), 
                                                 {nl,nchol*nb}); 
                SpMatrix_ref TAnb(TMats.origin()+nchol_max,{nocca_tot,nchol*nb}); 

Timer.start("T1");  
                // T(A,n,b) = sum_l_in_K G(A,l) LQKlnb[Q][K](l,n,b)
                ma::product(G3Da[n]({0,nocca_tot},{nl0,nl0+nl}),Llnb,TAnb);  
Timer.stop("T1"); 

Timer.start("T2");  
                if(addEJ) {
                  haveKE=true;
                  int nb0 = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+QKToK2[Q][Kl],0);
                  int nb = nelpk[nd][QKToK2[Q][Kl]]; 
                  boost::multi::array_ref<SPComplexType,3> Tanb(TAnb[nb0].origin(),{nb,nchol,nb});
                  for(int nc=0; nc<nchol; nc++)  
                    for(int b=0; b<nb; b++) 
                      Kr_local[0][nc] += Tanb[b][nc][b]; 
                }
Timer.stop("T2"); 

Timer.start("T3");  
                int na0 = 0;
                for(int Ka=0; Ka<nkpts; ++Ka) {  
                  int na = nelpk[nd][Ka]; 
                  int nk = nopk[QKToK2[Q][Ka]];
                  int nk0 = std::accumulate(nopk.begin(),nopk.begin()+QKToK2[Q][Ka],0);

                  SpMatrix_ref Fbk(TMats.origin()+nchol_max+TAnb.num_elements(),{nb,nk}); 
                  SpMatrix_ref Lank(std::addressof(*LQKank[nd*nkpts+Q][Ka].origin()), 
                                                 {na*nchol,nk}); 
                  SpMatrix_ref Tanb(TAnb[na0].origin(),{na*nchol,nb}); 
              
                  // F(K1,b,k) = sum_a_in_K1 sum_n  T(a,n,b) * LQKank(a,n,k) 
                  ma::product(ma::T(Tanb),Lank,Fbk);

                  // EXX += sum_K1 sum_b_in_Q(K) sum_k_in_Q(K1) F(K1,b,k) * G(b,k)
                  ComplexType E_(0.0);
                  for(int b=0; b<nb; ++b) 
                    E_ += ma::dot(Fbk[b],G3Da[n][nb0+b]({nk0,nk0+nk}));  
                  E[n][1] -= scl*0.5*E_;
                  na0 += na;  
                }
Timer.stop("T3"); 
              }
              if(walker_type==COLLINEAR) {
                if((nqk++)%comm->size() == comm->rank()) {
                }
              }
              if(addEJ && haveKE) {
                std::lock_guard<shared_mutex> guard(*mutex[Q]);
                int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);  
                ma::axpy(SPComplexType(1.0),Kr_local[0].sliced(0,ncholpQ[Q]),Kr[n]({nc0,nc0+ncholpQ[Q]})); 
              } // to release the lock
              if(addEJ && haveKE) 
                std::fill_n(Kr_local.origin(),Kr_local.num_elements(),SPComplexType(0.0));
            }
          }
        }
      }  
Timer.stop("T0"); 
t2 = Timer.total("T0");          

Timer.reset("T0");
Timer.start("T0");  
      if(addEJ) {
        if(not addEXX) {
          // calculate Kr
          APP_ABORT(" Error: Finish addEJ and not addEXX");
        }
        SPComplexType one(1.0);
        size_t local_memory_needs = nchol_max;
        if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        size_t nqk=0;  // start count at 1 to "offset" the calcuation of E1 done at root
        // avoiding vectors for now
        SpMatrix_ref Kl_local(TMats.origin(),{1,nchol_max});
        std::fill_n(Kl_local.origin(),Kl_local.num_elements(),SPComplexType(0.0));
        for(int n=0; n<nwalk; ++n) {
          for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (a,k)
            for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
              bool haveKE=false;
              if((nqk++)%comm->size() == comm->rank()) {
                haveKE=true;
                int nchol = ncholpQ[Q];
                int na = nelpk[nd][K];
                int na0 = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+K,0);
                int nk = nopk[QKToK2[Q][K]];
                int nk0 = std::accumulate(nopk.begin(),nopk.begin()+QKToK2[Q][K],0);

                Sp3Tensor_ref Lank(std::addressof(*LQKank[nd*nkpts+Q][K].origin()),
                                                 {na,nchol,nk});

                for(int a=0; a<na; ++a) 
                  ma::product(one,Lank[a],G3Da[n][na0+a]({nk0,nk0+nk}),one,Kl_local[0]({0,nchol}));
              }
              if(walker_type==COLLINEAR) {
                if((nqk++)%comm->size() == comm->rank()) {
                }
              }
              if(addEJ && haveKE) {
                std::lock_guard<shared_mutex> guard(*mutex[Q]);
                int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
                ma::axpy(one,Kl_local[0].sliced(0,ncholpQ[Q]),Kl[n]({nc0,nc0+ncholpQ[Q]}));
              } // to release the lock
              if(addEJ && haveKE)
                std::fill_n(Kl_local.origin(),Kl_local.num_elements(),SPComplexType(0.0));
            }  
          }  
        }  
        comm->barrier();
        nqk=0;
        for(int n=0; n<nwalk; ++n) {
          for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
            if((nqk++)%comm->size() == comm->rank()) {
              int nc0 = std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
              E[n][2] += 0.5*scl*scl*ma::dot(Kl[n]({nc0,nc0+ncholpQ[Q]}),
                                            Kr[n]({nc0,nc0+ncholpQ[Q]}));  
            }
          }
        }
      }
Timer.stop("T0"); 
t3 = Timer.total("T0");          
/*
std::cout<<" Time in energy - " 
<<"E1: " <<t1 <<"  " 
<<"EX: " <<t2 <<"  (" <<Timer.total("T1") <<"," <<Timer.total("T2") <<"," <<Timer.total("T3") <<")  "
<<"EJ: " <<t3 <<"\n"; 
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
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      using BType = typename std::decay<MatB>::type::element ;
      using AType = typename std::decay<MatA>::type::element ;
      boost::multi_array_ref<BType,2> v_(std::addressof(*v.origin()),
                                        extents[v.shape()[0]][1]);
      boost::multi_array_ref<AType,2> X_(std::addressof(*X.origin()),
                                        extents[X.shape()[0]][1]);
      return vHS(X_,v_,a,c);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      int nkpts = nopk.size();
      int nwalk = X.shape()[1];
      assert(v.shape()[0]==nwalk);
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      assert(X.num_elements() == nwalk*2*local_nCV);
      assert(v.num_elements() == nwalk*nmo_tot*nmo_tot);
      SPComplexType one(1.0,0.0);
      SPComplexType im(0.0,1.0);
      size_t local_memory_needs = nmo_max*nmo_max*nwalk; 
      if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});

      Sp3Tensor_ref v3D(std::addressof(*v.origin()),{nwalk,nmo_tot,nmo_tot});

      // "rotate" X  
      //  XIJ = 0.5*a*(Xn+ -i*Xn-), XJI = 0.5*a*(Xn+ +i*Xn-)  
      for(int Q=0, nq=0; Q<nkpts; ++Q) { 
        int nc0, ncN;  
        std::tie(nc0,ncN) = FairDivideBoundary(comm->rank(),ncholpQ[Q],comm->size()); 
        auto Xnp = std::addressof(*X[nq+nc0].origin());
        auto Xnm = std::addressof(*X[nq+ncholpQ[Q]+nc0].origin());
        for(int n=nc0; n<ncN; ++n) { 
          for(int nw=0; nw<nwalk; ++nw, ++Xnp, ++Xnm) {
            ComplexType Xnp_ = 0.5*a*((*Xnp) -im*(*Xnm)); 
            *Xnm =  0.5*a*((*Xnp) + im*(*Xnm));
            *Xnp = Xnp_;
          }
        }  
        nq+=2*ncholpQ[Q];
      } 
      // scale v by 'c': assuming contiguous data 
      {
        size_t i0, iN;
        std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(v.num_elements()),size_t(comm->size()));
        auto v_ = std::addressof(*v.origin())+i0;
        for(size_t i=i0; i<iN; ++i, ++v_)
          *v_ *= c;
      }
      comm->barrier();     
        
      using ma::T;  
      using ma::H;  
      size_t nqk=0;  
      for(int Q=0, nc0=0; Q<nkpts; ++Q) {      // momentum conservation index   
        for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
          if((nqk++)%comm->size() == comm->rank()) {
            int nchol = ncholpQ[Q];
            int ni = nopk[K];
            int ni0 = std::accumulate(nopk.begin(),nopk.begin()+K,0);
            int nk = nopk[QKToK2[Q][K]];
            int nk0 = std::accumulate(nopk.begin(),nopk.begin()+QKToK2[Q][K],0);

            SpMatrix_ref Likn(std::addressof(*LQKikn[Q][K].origin()),
                              {ni*nk,nchol});
            SpMatrix_ref vik(TMats.origin(),{nwalk,ni*nk});
            Sp3Tensor_ref vik3D(TMats.origin(),{nwalk,ni,nk});

            // v[nw][i(in K)][k(in Q(K))] += sum_n LQK[i][k][n] X[Q][n+][nw]
            ma::product(T(X[indices[range_t(nc0,nc0+nchol)][range_t()]]),
                        T(Likn),vik);

            // it is possible to add the second half here by calculating the (Q*,K*) pair that maps
            // to the JI term corresponding to this (Q,K) pair. Not doing it for now

            for(int nw=0; nw<nwalk; nw++)
              for(int i=0; i<ni; i++)
                ma::axpy(one,vik3D[nw][i],v3D[nw][ni0+i].sliced(nk0,nk0+nk));

          }
        }
        nc0+=2*ncholpQ[Q];
      }
      // adding second half. sync here to avoid need for locks  
      comm->barrier();
      nqk=0;
      for(int Q=0, nc0=0; Q<nkpts; ++Q) {      // momentum conservation index   
        for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (i,k)
          if((nqk++)%comm->size() == comm->rank()) {
            int nchol = ncholpQ[Q];
            int ni = nopk[K];
            int ni0 = std::accumulate(nopk.begin(),nopk.begin()+K,0);
            int nk = nopk[QKToK2[Q][K]];
            int nk0 = std::accumulate(nopk.begin(),nopk.begin()+QKToK2[Q][K],0);

            SpMatrix_ref Likn(std::addressof(*LQKikn[Q][K].origin()),
                              {ni*nk,nchol});
            SpMatrix_ref vik(TMats.origin(),{nwalk,ni*nk});
            Sp3Tensor_ref vik3D(TMats.origin(),{nwalk,ni,nk});

            // v[nw][k(in Q(K))][i(in K)] += sum_n conj(LQK[i][k][n]) X[Q][n-][nw]
            ma::product(T(X[indices[range_t(nc0+nchol,nc0+2*nchol)][range_t()]]),
                        H(Likn),vik);

            for(int nw=0; nw<nwalk; nw++) {
              auto&& vik3D_n = vik3D[nw];   
              for(int k=0; k<nk; k++) {
                ComplexType* v3D_nk = std::addressof(*v3D[nw][nk0+k].origin()) + ni0;    
                for(int i=0; i<ni; i++, ++v3D_nk) 
                  *v3D_nk += vik3D_n[i][k];  
              }
            }
          }
        }
        nc0+=2*ncholpQ[Q];
      }    
      comm->barrier();
      // do I need to "rotate" back, can be done if necessary
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
      using BType = typename std::decay<MatB>::type::element ;
      using AType = typename std::decay<MatA>::type::element ;
      boost::multi_array_ref<BType,2> v_(std::addressof(*v.origin()),
                                        extents[v.shape()[0]][1]);
      boost::const_multi_array_ref<AType,2> G_(std::addressof(*G.origin()),
                                        extents[G.shape()[0]][1]);
      return vbias(G_,v_,a,c,k);  
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int nd=0) {
      int nkpts = nopk.size();
      assert(nd >= 0 && nd < nelpk.size());
      int nwalk = G.shape()[1];
      assert(v.shape()[0]==2*local_nCV);  
      assert(v.shape()[1]==nwalk);  
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
      size_t local_memory_needs = 2*nchol_max*nwalk + std::max(nocca_max,noccb_max)*nwalk;
      if(TMats.num_elements() < local_memory_needs) TMats.reextent({local_memory_needs,1});
      SpMatrix_ref vlocal(TMats.origin(),{2*nchol_max,nwalk});
      std::fill_n(vlocal.origin(),vlocal.num_elements(),SPComplexType(0.0));
      SpMatrix_ref Gl(TMats.origin()+vlocal.num_elements(),{std::max(nocca_max,noccb_max),nwalk});

      assert(G.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      boost::multi::const_array_ref<ComplexType,3> G3Da(std::addressof(*G.origin()),
                                                        {nocca_tot,nmo_tot,nwalk} );
      boost::multi::const_array_ref<ComplexType,3> G3Db(std::addressof(*G.origin())+
                                                        G3Da.num_elements()*(nspin-1),
                                                        {noccb_tot,nmo_tot,nwalk} );


      {  
        size_t i0, iN;
        std::tie(i0,iN) = FairDivideBoundary(size_t(comm->rank()),size_t(v.shape()[0]),size_t(comm->size()));
        for(size_t i=i0; i<iN; ++i) 
          ma::scal(c,v[i]);
      }
      comm->barrier();

      size_t nqk=0;  // start count at 1 to "offset" the calcuation of E1 done at root
      for(int K=0; K<nkpts; ++K) {        // K is the index of the kpoint pair of (a,k)
        for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
          bool haveV=false;
          if((nqk++)%comm->size() == comm->rank()) {
            haveV=true;
            int nchol = ncholpQ[Q];
            int na = nelpk[nd][K];
            int na0 = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+K,0);
            int nk = nopk[QKToK2[Q][K]];
            int nk0 = std::accumulate(nopk.begin(),nopk.begin()+QKToK2[Q][K],0);
            int nb = nelpk[nd][QKToK2[Q][K]];
            int nb0 = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+QKToK2[Q][K],0);
            int nl = nopk[K];
            int nl0 = std::accumulate(nopk.begin(),nopk.begin()+K,0);
            auto&& v1 = vlocal({0,nchol},{0,nwalk});
            auto&& v2 = vlocal({nchol,2*nchol},{0,nwalk});

            Sp3Tensor_ref Lank(std::addressof(*LQKank[nd*nkpts+Q][K].origin()),
                                                 {na,nchol,nk});
            Sp3Tensor_ref Llnb(std::addressof(*LQKlnb[nd*nkpts+Q][K].origin()),
                                                 {nl,nchol,nb});

            // v1[Q][n][nw] += sum_K sum_a_k LQK[a][n][k] G[a][k][nw]
            for(int a=0; a<na; ++a) 
              ma::product(one,Lank[a],G3Da[na0+a]({nk0,nk0+nk},{0,nwalk}),one,v1);

            // v2[Q][n][nw] += sum_K sum_l_b LQK[l][n][b] G[b][l][nw]
            for(int l=0; l<nl; ++l) {
              // might be faster to transpose G -> G[l][nw][b] and then use ma::T(G[l])  
              for(int b=0; b<nb; b++)
                std::copy_n(G3Da[nb0+b][nl0+l].origin(),nwalk,Gl[b].origin());  
              ma::product(one,Llnb[l],Gl({0,nb},{0,nwalk}),one,v2);
            }
          }
          if(walker_type==COLLINEAR) {
            if((nqk++)%comm->size() == comm->rank()) {
            }
          }
          if(haveV) {
            std::lock_guard<shared_mutex> guard(*mutex[Q]);
            int nc0 = 2*std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
            // v+ = 0.5*a*(v1+v2) 
            ma::axpy(halfa,vlocal({0,ncholpQ[Q]},{0,nwalk}),v[indices[range_t(nc0,nc0+ncholpQ[Q])][range_t()]]);
            ma::axpy(halfa,vlocal({ncholpQ[Q],2*ncholpQ[Q]},{0,nwalk}),v[indices[range_t(nc0,nc0+ncholpQ[Q])][range_t()]]);
            // v- = -0.5*a*i*(v1-v2) 
            ma::axpy(minusimhalfa,vlocal({0,ncholpQ[Q]},{0,nwalk}),v[indices[range_t(nc0+ncholpQ[Q],nc0+2*ncholpQ[Q])][range_t()]]);
            ma::axpy(imhalfa,vlocal({ncholpQ[Q],2*ncholpQ[Q]},{0,nwalk}),v[indices[range_t(nc0+ncholpQ[Q],nc0+2*ncholpQ[Q])][range_t()]]);
          } // to release the lock
          if(haveV)
            std::fill_n(vlocal.origin(),vlocal.num_elements(),SPComplexType(0.0));
        }
      }
      comm->barrier();
    }

    bool distribution_over_cholesky_vectors() const{ return true; }
    int number_of_ke_vectors() const{ return std::accumulate(ncholpQ.begin(),ncholpQ.end(),0); }
    int local_number_of_cholesky_vectors() const{ return 2*std::accumulate(ncholpQ.begin(),ncholpQ.end(),0); } 
    int global_number_of_cholesky_vectors() const{ return global_nCV; }

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const{return false;}
    bool transposed_G_for_E() const{return true;} 
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const{return true;} 

    bool fast_ph_energy() const { return false; }

  private:

    communicator* comm;

    WALKER_TYPES walker_type;

    int global_nCV;
    int local_nCV;

    ValueType E0;

    // bare one body hamiltonian
    shmC3Tensor H1;

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
    shmC3Tensor vn0;

    // shared buffer space
    // using matrix since there are issues with vectors
    shmSpMatrix SM_TMats;
    SpMatrix TMats;

    std::vector<std::unique_ptr<shared_mutex>> mutex;

    myTimer Timer;

    void set_shm_buffer(size_t N) {
      if(SM_TMats.num_elements() < N) 
        SM_TMats.reextent({N,1});
    }

};

}

}

#endif
