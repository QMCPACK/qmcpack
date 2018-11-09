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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KPTHCOPS_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KPTHCOPS_HPP

#include<fstream>
#include<mutex>

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "alf/boost/mpi3/shm/mutex.hpp"
#include "alf/boost/mpi3/shared_communicator.hpp"
#include "AFQMC/multi/array.hpp"
#include "AFQMC/multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "type_traits/scalar_traits.h"
#include "AFQMC/Wavefunctions/Excitations.hpp"
#include "AFQMC/Wavefunctions/phmsd_helpers.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class KPTHCOps
{
#if defined(AFQMC_SP) 
  using SpC = typename to_single_precision<ComplexType>::value_type;
#else
  using SpC = ComplexType;  
#endif

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
  using this_t = KPTHCOps;

  public:

    /*
     * NAOA/NAOB stands for number of active orbitals alpha/beta (instead of active electrons)
     */
    KPTHCOps(communicator& c_,
                 WALKER_TYPES type,
                 std::vector<int>&& nopk_,
                 std::vector<int>&& ncholpQ_,
                 shmIMatrix&& nelpk_,
                 shmIMatrix&& QKToK2_,
                 shmIMatrix&& QKToG_,
                 shmC3Tensor&& h1_,
                 shmCMatrix&& haj_,
                 std::vector<shmSpMatrix>&& rotlun_,
                 shmCMatrix&& rotpiu_,
                 std::vector<shmSpMatrix>&& rotpau_,
                 std::vector<shmSpMatrix>&& lun_,
                 shmSpMatrix&& piu_,
                 std::vector<shmCMatrix>&& pau_,
                 shmC3Tensor&& vn0_,
                 ValueType e0_,
                 int gncv,
                 bool low_mem = true,
                 bool verbose=false ):
                comm(std::addressof(c_)),
                low_memory(low_mem),
                walker_type(type),
                global_nCV(gncv), 
                H1(std::move(h1_)),
                haj(std::move(haj_)),
                nopk(std::move(nopk_)),
                ncholpQ(std::move(ncholpQ_)),
                nelpk(std::move(nelpk_)),
                QKToK2(std::move(QKToK2_)),
                QKToG(std::move(QKToG_)),
                rotLQGun(std::move(rotlun_)),
                rotPiu(std::move(rotpiu_)),
                rotcPua(std::move(rotpau_)),
                LQGun(std::move(lun_)),
                Piu(std::move(piu_)),
                cPua(std::move(pau_)),
                //rotMuv(),
                vn0(std::move(vn0_)),
                E0(e0_),
                SM_TMats({1,1},shared_allocator<SPComplexType>{c_}),
                TMats({1,1}),
                mutex(0)
    {
      nGpk.resize(nopk.size());
      for(int Q=0; Q<nopk.size(); Q++) 
        nGpk[Q] = *std::max_element(QKToG[Q].begin(),QKToG[Q].end())+1;
      int nu = Piu.shape()[1]; 
      local_nCV = std::accumulate(ncholpQ.begin(),ncholpQ.end(),0);
      mutex.reserve(ncholpQ.size());
      for(int nQ=0; nQ<ncholpQ.size(); nQ++)
          mutex.emplace_back(std::make_unique<shared_mutex>(*comm));
      if(not low_memory) {  
        APP_ABORT(" Error: Finish implementation of low_memory in KPTHCOps.\n");    
        //rotMuv.emplace_back(shmSpMatrix({nu,nu},shared_allocator<SPComplexType>{c_}));
      }  
/*
      if(haj.size() > 1)
	APP_ABORT(" Error: THC not yet implemented for multiple references.\n");	
      assert(comm);
      // current partition over 'u' for L/Piu
      assert(Luv.shape()[0] == Piu.shape()[1]);
      for(int i=0; i<rotcPua.size(); i++) {
        // rot Ps are not yet distributed
        assert(rotcPua[i].shape()[0] == rotPiu.shape()[1]);
        if(walker_type==CLOSED)
          assert(rotcPua[i].shape()[1]==NAOA);
        else if(walker_type==COLLINEAR)
          assert(rotcPua[i].shape()[1]==NAOA+NAOB);
        else if(walker_type==NONCOLLINEAR)
          assert(rotcPua[i].shape()[1]==NAOA+NAOB);
      }
      for(int i=0; i<cPua.size(); i++) {
        assert(cPua[i].shape()[0]==Luv.shape()[0]);
        if(walker_type==CLOSED)
          assert(cPua[i].shape()[1]==NAOA);
        else if(walker_type==COLLINEAR)
          assert(cPua[i].shape()[1]==NAOA+NAOB);
        else if(walker_type==NONCOLLINEAR)
          assert(cPua[i].shape()[1]==NAOA+NAOB);
      }
      if(walker_type==NONCOLLINEAR) {
        assert(Piu.shape()[0]==2*NMO);
        assert(rotPiu.shape()[0]==2*NMO);
      } else {
        assert(Piu.shape()[0]==NMO);
        assert(rotPiu.shape()[0]==NMO);
      }
*/
    }

    ~KPTHCOps() {}
    
    KPTHCOps(KPTHCOps const& other) = delete;
    KPTHCOps& operator=(KPTHCOps const& other) = delete;

    KPTHCOps(KPTHCOps&& other) = default;
    KPTHCOps& operator=(KPTHCOps&& other) = default; 

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

    // Kl and Kr must be in shared memory for this to work correctly  
    template<class Mat, class MatB, class MatC, class MatD>
    void energy(Mat&& E, MatB const& Gc, int nd, MatC* KEleft, MatD* KEright, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      using ma::T;
      using ma::H;
      using std::conj;
      if(nd>0)
	APP_ABORT(" Error: KPTHC not yet implemented for multiple references.\n");	
      static_assert(E.dimensionality==2);  
      static_assert(Gc.dimensionality==2);  
      assert(E.shape()[0] == Gc.shape()[0]);        
      assert(E.shape()[1] == 3);        
      assert(nd >= 0 && nd < nelpk.size());

      int nkpts = nopk.size();
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

      for(int n=0; n<nwalk; n++)
        std::fill_n(E[n].origin(),3,ComplexType(0.));

      assert(Gc.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      boost::multi::const_array_ref<ComplexType,3> G3Da(std::addressof(*Gc.origin()),
                                                        {nwalk,nocca_tot,nmo_tot} );
      boost::multi::const_array_ref<ComplexType,3> G3Db(std::addressof(*Gc.origin())+
                                                        G3Da.num_elements()*(nspin-1),
                                                        {nwalk,noccb_tot,nmo_tot} );

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

      if(not addEXX and not addEJ) return;

      SPComplexType *Krptr, *Klptr;
      int getKr = KEright!=nullptr;
      int getKl = KEleft!=nullptr;
      int nu = Piu.shape()[1]; 
      size_t memory_needs = nkpts*nkpts*nu*nu; 
      if(addEJ) {
        if(not getKr) memory_needs += nwalk*local_nCV;
        if(not getKl) memory_needs += nwalk*local_nCV;
      }
      set_shm_buffer(memory_needs);
      size_t cnt=0;  
      // Fuv[k1][k2][u][v] = sum_a_l cPua[u][k1][a] * G[k1][a][k2][l] Piu[k2][l][v]
      boost::multi_array_ref<ComplexType,4> Fuv(std::addressof(*SM_TMats.origin()),extents[nkpts][nkpts][nu][nu]);
      cnt+=Fuv.num_elements();

      // messy
      size_t Knr=0, Knc=0;
      if(addEJ) {
        Knr=nwalk;
        Knc=local_nCV;
        if(getKr) {
          assert(KEright->shape()[0] == nwalk && KEright->shape()[1] == local_nCV);
          assert(KEright->strides()[0] == KEright->shape()[1]);
          Krptr = std::addressof(*KEright->origin());
        } else {
          Krptr = std::addressof(*SM_TMats.origin())+cnt;
          cnt += nwalk*local_nCV;
        }
        if(getKl) {
          assert(KEleft->shape()[0] == nwalk && KEleft->shape()[1] == local_nCV);
          assert(KEleft->strides()[0] == KEleft->shape()[1]);
          Klptr = std::addressof(*KEleft->origin());
        } else {
          Klptr = std::addressof(*SM_TMats.origin())+cnt;
          cnt += nwalk*local_nCV;
        }
        if(comm->root()) std::fill_n(Krptr,Knr*Knc,SPComplexType(0.0));
        if(comm->root()) std::fill_n(Klptr,Knr*Knc,SPComplexType(0.0));
      } else if(getKr or getKl) {
        APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
      }
      SpMatrix_ref Kl(Klptr,{Knr,Knc});
      SpMatrix_ref Kr(Krptr,{Knr,Knc});
      comm->barrier();

/*

      if(walker_type==CLOSED || walker_type==NONCOLLINEAR) {
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        for(int wi=0; wi<nwalk; wi++) {
          boost::const_multi_array_ref<ComplexType,2> Gw(G[wi].origin(),extents[nel_][nmo_]);
          boost::const_multi_array_ref<ComplexType,1> G1D(G[wi].origin(),extents[nel_*nmo_]);
          // need a new routine if addEXX is false, 
          // otherwise it is quite inefficient to get Ej only
          Guv_Guu(Gw,Guv,Guu,T1,k);
          if(addEJ) {
            ma::product(rotMuv.get()[indices[range_t(u0,uN)][range_t()]],Guu,
                        Tuu[indices[range_t(u0,uN)]]);
            if(getKl)
              std::copy_n(std::addressof(*Guu.origin())+nu0+u0,uN-u0,std::addressof(*(*Kl)[wi].origin())+u0);
            if(getKr)
              std::copy_n(std::addressof(*Tuu.origin())+u0,uN-u0,std::addressof(*(*Kr)[wi].origin())+u0);
            E[wi][2] = 0.5*scl*scl*ma::dot(Guu[indices[range_t(nu0+u0,nu0+uN)]],Tuu[indices[range_t(u0,uN)]]); 
          }
          if(addEXX) {
            auto Mptr = rotMuv.get()[u0].origin();  
            auto Gptr = Guv[0][u0].origin();  
            for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
              (*Gptr) *= (*Mptr); 
            ma::product(Guv[0][indices[range_t(u0,uN)][range_t()]],rotcPua[k].get(),
                        Qub[indices[range_t(u0,uN)][range_t()]]);
            // using this for now, which should not be much worse
            ma::product(T(Qub[indices[range_t(u0,uN)][range_t()]]),
                        T(rotPiu.get()[indices[range_t()][range_t(nu0+u0,nu0+uN)]]),  
                        Rbk);
            E[wi][1] = -0.5*scl*ma::dot(R1D,G1D);
          }
        }
      } else {
        for(int wi=0; wi<nwalk; wi++) {
          boost::const_multi_array_ref<ComplexType,2> Gw(G[wi].origin(),extents[nel_][nmo_]);
          boost::const_multi_array_ref<ComplexType,1> G1DA(G[wi].origin(),extents[NAOA*nmo_]);
          boost::const_multi_array_ref<ComplexType,1> G1DB(G[wi].origin()+NAOA*nmo_,extents[NAOB*nmo_]);
          Guv_Guu(Gw,Guv,Guu,T1,k);
          // move calculation of Guv/Guu here to avoid storing 2 copies of Guv for alpha/beta
          if(addEJ) {
            ma::product(rotMuv.get()[indices[range_t(u0,uN)][range_t()]],Guu,
                      Tuu[indices[range_t(u0,uN)]]);
            if(getKl)
              std::copy_n(std::addressof(*Guu.origin())+nu0+u0,uN,std::addressof(*(*Kl)[wi].origin())+u0);
            if(getKr)
              std::copy_n(std::addressof(*Tuu.origin())+u0,uN,std::addressof(*(*Kr)[wi].origin())+u0);
            E[wi][2] = 0.5*ma::dot(Guu[indices[range_t(nu0+u0,nu0+uN)]],Tuu[indices[range_t(u0,uN)]]);
          }
          if(addEXX) {
            // alpha
            auto Mptr = rotMuv.get()[u0].origin();
            auto Gptr = Guv[0][u0].origin();
            for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
              (*Gptr) *= (*Mptr);
            ma::product(Guv[indices[0][range_t(u0,uN)][range_t()]],
                      (rotcPua[k].get())[indices[range_t()][range_t(0,NAOA)]],
                      Qub[indices[range_t(u0,uN)][range_t(0,NAOA)]]);
            // using this for now, which should not be much worse
            ma::product(T(Qub[indices[range_t(u0,uN)][range_t(0,NAOA)]]),
                      T(rotPiu.get()[indices[range_t()][range_t(nu0+u0,nu0+uN)]]),
                      Rbk[indices[range_t(0,NAOA)][range_t()]]);
            E[wi][1] = -0.5*ma::dot(R1D[indices[range_t(0,NAOA*nmo_)]],G1DA);
            // beta
            Mptr = rotMuv.get()[u0].origin();
            Gptr = Guv[1][u0].origin();
            for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
              (*Gptr) *= (*Mptr);
            ma::product(Guv[indices[0][range_t(u0,uN)][range_t()]],
                      (rotcPua[k].get())[indices[range_t()][range_t(NAOA,NAOA+NAOB)]],
                      Qub[indices[range_t(u0,uN)][range_t(0,NAOB)]]);
            // using this for now, which should not be much worse
            ma::product(T(Qub[indices[range_t(u0,uN)][range_t(0,NAOB)]]),
                      T(rotPiu.get()[indices[range_t()][range_t(nu0+u0,nu0+uN)]]),
                      Rbk[indices[range_t(0,NAOB)][range_t()]]);
            E[wi][1] -= 0.5*ma::dot(R1D[indices[range_t(0,NAOB*nmo_)]],G1DB);
          }
        }
      }    
      comm->barrier();
*/
    }

    template<class MatE, class MatO, class MatG, class MatQ, class MatB, 
             class index_aos>
    void fast_energy(MatE&& E, MatO&& Ov, MatG const& GrefA, MatG const& GrefB, 
                     MatQ const& QQ0A, MatQ const& QQ0B, MatB&& Qwork,  
                     ph_excitations<int,ComplexType> const& abij,
                     std::array<index_aos,2> const& det_couplings)
    {
/*
      if(haj.size() != 1)
        APP_ABORT(" Error: Single reference implementation currently in KPTHCOps::fast_energy.\n");
      if(walker_type!=CLOSED)
        APP_ABORT(" Error: KPTHCOps::fast_energy requires walker_type==CLOSED.\n");
      static_assert(E.dimensionality==4);
      static_assert(Ov.dimensionality==3);
      static_assert(GrefA.dimensionality==3);
      static_assert(GrefB.dimensionality==3);
      static_assert(QQ0A.dimensionality==3);
      static_assert(QQ0B.dimensionality==3);
      int nspin = E.shape()[0];
      int nrefs = haj.size();
      int nwalk = GrefA.shape()[0];
      int naoa_ = QQ0A.shape()[1];
      int naob_ = QQ0B.shape()[1];
      int naea_ = QQ0A.shape()[2];
      int naeb_ = QQ0B.shape()[2];
      int nmo_ = rotPiu.shape()[0];
      int nu = rotMuv.shape()[0];
      int nu0 = rotMuv.offset()[0];
      int nv = rotMuv.shape()[1];
      int nel_ = rotcPua[0].shape()[1];
      // checking
      assert(E.shape()[2] == nwalk);
      assert(E.shape()[3] == 3);
      assert(Ov.shape()[0] == nspin);
      assert(Ov.shape()[1] == E.shape()[1]);
      assert(Ov.shape()[2] == nwalk);
      assert(GrefA.shape()[1] == naoa_);
      assert(GrefA.shape()[2] == nmo_);
      assert(GrefB.shape()[0] == nwalk);
      assert(GrefB.shape()[1] == naob_);
      assert(GrefB.shape()[2] == nmo_);
      // limited to single reference now
      assert(rotcPua.size() == nrefs);
      assert(nel_ == naoa_);
      assert(nel_ == naob_);

      using ma::T;
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int k0,kN;
      std::tie(k0,kN) = FairDivideBoundary(comm->rank(),nel_,comm->size());
      // right now the algorithm uses 2 copies of matrices of size nuxnv in COLLINEAR case, 
      // consider moving loop over spin to avoid storing the second copy which is not used  
      // simultaneously
      size_t memory_needs = nu*nv + nv + nu  + nel_*(nv+2*nu+2*nel_);
      set_shm_buffer(memory_needs);
      size_t cnt=0;
      // if Alpha/Beta have different references, allocate the largest and 
      // have distinct references for each
      // Guv[nu][nv]
      boost::multi_array_ref<ComplexType,2> Guv(SM_TMats->data(),extents[nu][nv]);
      cnt+=Guv.num_elements();
      // Gvv[v]: summed over spin
      boost::multi_array_ref<ComplexType,1> Gvv(SM_TMats->data()+cnt,extents[nv]);
      cnt+=Gvv.num_elements();
      // S[nel_][nv]
      boost::multi_array_ref<ComplexType,2> Scu(SM_TMats->data()+cnt,extents[nel_][nv]);
      cnt+=Scu.num_elements();
      // Qub[nu][nel_]: 
      boost::multi_array_ref<ComplexType,2> Qub(SM_TMats->data()+cnt,extents[nu][nel_]);
      cnt+=Qub.num_elements();
      boost::multi_array_ref<ComplexType,1> Tuu(SM_TMats->data()+cnt,extents[nu]);
      cnt+=Tuu.num_elements();
      boost::multi_array_ref<ComplexType,2> Jcb(SM_TMats->data()+cnt,extents[nel_][nel_]);
      cnt+=Jcb.num_elements();
      boost::multi_array_ref<ComplexType,2> Xcb(SM_TMats->data()+cnt,extents[nel_][nel_]);
      cnt+=Xcb.num_elements();
      boost::multi_array_ref<ComplexType,2> Tub(SM_TMats->data()+cnt,extents[nu][nel_]);
      cnt+=Tub.num_elements();
      assert(cnt <= memory_needs);
      if(eloc.shape()[0] != 2 || eloc.shape()[1] != nwalk || eloc.shape()[2] != 3) 
        eloc.resize(extents[2][nwalk][3]);

      std::fill_n(eloc.origin(),eloc.num_elements(),ComplexType(0.0));

      RealType scl = (walker_type==CLOSED?2.0:1.0);
      if(comm->root()) {
        std::fill_n(std::addressof(*E.origin()),E.num_elements(),ComplexType(0.0));
        std::fill_n(std::addressof(*Ov[0][1].origin()),nwalk*(Ov.shape()[1]-1),ComplexType(0.0));
        std::fill_n(std::addressof(*Ov[1][1].origin()),nwalk*(Ov.shape()[1]-1),ComplexType(0.0));
        auto Ea = E[0][0];
        auto Eb = E[1][0];
        boost::const_multi_array_ref<ComplexType,2> G2DA(std::addressof(*GrefA.origin()),
                                          extents[nwalk][GrefA[0].num_elements()]);
        ma::product(ComplexType(1.0),G2DA,haj[0],ComplexType(0.0),Ea(Ea.extension(0),0));
        boost::const_multi_array_ref<ComplexType,2> G2DB(std::addressof(*GrefA.origin()),
                                          extents[nwalk][GrefA[0].num_elements()]);
        ma::product(ComplexType(1.0),G2DB,haj[0],ComplexType(0.0),Eb(Eb.extension(0),0));
        for(int i=0; i<nwalk; i++) {
            Ea[i][0] += E0;
            Eb[i][0] += E0;
        }
      }

      for(int wi=0; wi<nwalk; wi++) {

        { // Alpha
          auto Gw = GrefA[wi];
          boost::const_multi_array_ref<ComplexType,1> G1D(std::addressof(*Gw.origin()),
                                                        extents[Gw.num_elements()]);
          Guv_Guu2(Gw,Guv,Gvv,Scu,0);
          if(u0!=uN)
            ma::product(rotMuv.get()[indices[range_t(u0,uN)][range_t()]],Gvv,
                      Tuu[indices[range_t(u0,uN)]]);
          auto Mptr = rotMuv.get()[u0].origin();
          auto Gptr = Guv[u0].origin();
          for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
            (*Gptr) *= (*Mptr);
          if(u0!=uN)
            ma::product(Guv[indices[range_t(u0,uN)][range_t()]],rotcPua[0].get(),
                      Qub[indices[range_t(u0,uN)][range_t()]]);
          comm->barrier();
          if(k0!=kN)
            ma::product(Scu[indices[range_t(k0,kN)][range_t()]],Qub,
                      Xcb[indices[range_t(k0,kN)][range_t()]]);
          // Tub = rotcPua.*Tu
          auto rPptr = rotcPua[0].get()[nu0+u0].origin();
          auto Tuuptr = Tuu.origin()+u0;
          auto Tubptr = Tub[u0].origin();
          for(size_t u_=u0; u_<uN; ++u_, ++Tuuptr)
            for(size_t k=0; k<nel_; ++k, ++rPptr, ++Tubptr)
              (*Tubptr) = (*Tuuptr)*(*rPptr);
          comm->barrier();
          // Jcb = Scu*Tub
          if(k0!=kN)
            ma::product(Scu[indices[range_t(k0,kN)][range_t()]],Tub,
                      Jcb[indices[range_t(k0,kN)][range_t()]]);
          for(int c=k0; c<kN; ++c) 
            eloc[0][wi][1] += -0.5*scl*Xcb[c][c];
          for(int c=k0; c<kN; ++c) 
            eloc[0][wi][2] += 0.5*scl*scl*Jcb[c][c];
        }

        { // Beta: Unnecessary in CLOSED walker type (on Walker)
          auto Gw = GrefB[wi];
          boost::const_multi_array_ref<ComplexType,1> G1D(std::addressof(*Gw.origin()),
                                                        extents[Gw.num_elements()]);
          Guv_Guu2(Gw,Guv,Gvv,Scu,0);
          if(u0!=uN)
            ma::product(rotMuv.get()[indices[range_t(u0,uN)][range_t()]],Gvv,
                      Tuu[indices[range_t(u0,uN)]]);
          auto Mptr = rotMuv.get()[u0].origin();
          auto Gptr = Guv[u0].origin();
          for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
            (*Gptr) *= (*Mptr);
          if(u0!=uN)
            ma::product(Guv[indices[range_t(u0,uN)][range_t()]],rotcPua[0].get(),
                      Qub[indices[range_t(u0,uN)][range_t()]]);
          comm->barrier();
          if(k0!=kN)
            ma::product(Scu[indices[range_t(k0,kN)][range_t()]],Qub,
                      Xcb[indices[range_t(k0,kN)][range_t()]]);
          // Tub = rotcPua.*Tu
          auto rPptr = rotcPua[0].get()[nu0+u0].origin();
          auto Tuuptr = Tuu.origin()+u0;
          auto Tubptr = Tub[u0].origin();
          for(size_t u_=u0; u_<uN; ++u_, ++Tuuptr)
            for(size_t k=0; k<nel_; ++k, ++rPptr, ++Tubptr)
              (*Tubptr) = (*Tuuptr)*(*rPptr);
          comm->barrier();
          // Jcb = Scu*Tub
          if(k0!=kN)
            ma::product(Scu[indices[range_t(k0,kN)][range_t()]],Tub,
                      Jcb[indices[range_t(k0,kN)][range_t()]]);
          for(int c=k0; c<kN; ++c)
            eloc[1][wi][1] += -0.5*scl*Xcb[c][c];
          for(int c=k0; c<kN; ++c)
            eloc[1][wi][2] += 0.5*scl*scl*Jcb[c][c];
        }

      }
      comm->reduce_in_place_n(eloc.origin(),eloc.num_elements(),std::plus<>(),0);
      if(comm->root()) {
        // add Eref contributions to all configurations
        for(int nd=0; nd<E.shape()[1]; ++nd) {
          auto Ea = E[0][nd];
          auto Eb = E[1][nd];
          for(int wi=0; wi<nwalk; wi++) {
            Ea[wi][1] += eloc[0][wi][1];
            Ea[wi][2] += eloc[0][wi][2];
            Eb[wi][1] += eloc[1][wi][1];
            Eb[wi][2] += eloc[1][wi][2];
          }
        }
      }  
      comm->barrier();
*/
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vHS(MatA & X, MatB&& v, double a=1., double c=0.) {
        boost::multi_array_ref<ComplexType,2> X_(X.origin(),extents[X.shape()[0]][1]);
        boost::multi_array_ref<ComplexType,2> v_(v.origin(),extents[1][v.shape()[0]]);
        vHS(X_,v_,a,c);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(MatA & X, MatB&& v, double a=1., double c=0.) {
/*
      int nwalk = X.shape()[1];
#if defined(QMC_COMPLEX)
      int nchol = 2*Luv.shape()[1];
#else
      int nchol = Luv.shape()[1];
#endif
      int nmo_ = Piu.shape()[0];
      int nu = Piu.shape()[1];
      assert(Luv.shape()[0]==nu);
      assert(X.shape()[0]==nchol);
      assert(v.shape()[0]==nwalk);
      assert(v.shape()[1]==nmo_*nmo_);
      using ma::T;
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int k0,kN;
      std::tie(k0,kN) = FairDivideBoundary(comm->rank(),nmo_,comm->size());
      int wk0,wkN;
      std::tie(wk0,wkN) = FairDivideBoundary(comm->rank(),nwalk*nmo_,comm->size());
#define LOW_MEMORY
#if defined(LOW_MEMORY)
      size_t memory_needs = nu*nwalk + nu*nmo_;
#else
      size_t memory_needs = nu*nwalk + nwalk*nu*nmo_;
#endif
      set_shm_buffer(memory_needs);
      boost::multi_array_ref<ComplexType,2> Tuw(SM_TMats->data(),extents[nu][nwalk]);
      // O[nwalk * nmu * nmu]
//Timer.start("T0");
#if defined(QMC_COMPLEX)
      // reinterpret as RealType matrices with 2x the columns
      boost::multi_array_ref<RealType,2> Luv_R(reinterpret_cast<RealType*>(Luv.origin()),
                                                 extents[Luv.shape()[0]][2*Luv.shape()[1]]);
      boost::const_multi_array_ref<RealType,2> X_R(reinterpret_cast<RealType const*>(X.origin()),
                                                 extents[X.shape()[0]][2*X.shape()[1]]);
      boost::multi_array_ref<RealType,2> Tuw_R(reinterpret_cast<RealType*>(Tuw.origin()),
                                                 extents[nu][2*nwalk]);
      ma::product(Luv_R[indices[range_t(u0,uN)][range_t()]],X_R,
                  Tuw_R[indices[range_t(u0,uN)][range_t()]]);  
#else
      ma::product(Luv.get()[indices[range_t(u0,uN)][range_t()]],X,
                  Tuw[indices[range_t(u0,uN)][range_t()]]);  
#endif
      comm->barrier();
//Timer.stop("T0");
#if defined(LOW_MEMORY)
      boost::multi_array_ref<ComplexType,2> Qiu(SM_TMats->data()+nwalk*nu,extents[nmo_][nu]);
      for(int wi=0; wi<nwalk; wi++) {
        // Qiu[i][u] = T[u][wi] * conj(Piu[i][u])
        // v[wi][ik] = sum_u Qiu[i][u] * Piu[k][u]
        // O[nmo * nmu]
//Timer.start("T1");
        for(int i=k0; i<kN; i++) {
          auto p_ = Piu.get()[i].origin();  
          for(int u=0; u<nu; u++, ++p_)
            Qiu[i][u] = Tuw[u][wi]*conj(*p_);
        }
//Timer.stop("T1");
//Timer.start("T2");
        boost::multi_array_ref<ComplexType,2> v_(v[wi].origin(),extents[nmo_][nmo_]);
        // this can benefit significantly from 2-D partition of work
        // O[nmo * nmo * nmu]
        ma::product(a,Qiu[indices[range_t(k0,kN)][range_t()]],T(Piu.get()),
                    c,v_[indices[range_t(k0,kN)][range_t()]]);
//Timer.stop("T2");
      }
#else
      boost::multi_array_ref<ComplexType,2> Qiu(SM_TMats->data()+nwalk*nu,extents[nwalk*nmo_][nu]);
      boost::multi_array_ref<ComplexType,3> Qwiu(SM_TMats->data()+nwalk*nu,extents[nwalk][nmo_][nu]);
      // Qiu[i][u] = T[u][wi] * conj(Piu[i][u])
      // v[wi][ik] = sum_u Qiu[i][u] * Piu[k][u]
      // O[nmo * nmu]
//Timer.start("T1");
      for(int wi=0; wi<nwalk; wi++) 
        for(int i=k0; i<kN; i++) { 
          auto p_ = Piu.get()[i].origin();
          for(int u=0; u<nu; u++, ++p_)
            Qwiu[wi][i][u] = Tuw[u][wi]*conj(*p_);
        }
//Timer.stop("T1");
//Timer.start("T2");
      boost::multi_array_ref<ComplexType,2> v_(v.origin(),extents[nwalk*nmo_][nmo_]);
      // this can benefit significantly from 2-D partition of work
      // O[nmo * nmo * nmu]
      ma::product(a,Qiu[indices[range_t(wk0,wkN)][range_t()]],T(Piu.get()),
                  c,v_[indices[range_t(wk0,wkN)][range_t()]]);
//Timer.stop("T2");
#endif
      comm->barrier();
*/
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(MatA const& G, MatB&& v, double a=1., double c=0., int k=0) {
        boost::const_multi_array_ref<ComplexType,2> G_(G.origin(),extents[1][G.shape()[0]]);
        boost::multi_array_ref<ComplexType,2> v_(v.origin(),extents[v.shape()[0]][1]);
        vbias(G_,v_,a,c,k);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(MatA const& G, MatB&& v, double a=1., double c=0., int k=0) {
/*
      if(k>0)
	APP_ABORT(" Error: THC not yet implemented for multiple references.\n");	
      int nwalk = G.shape()[0];
      int nmo_ = Piu.shape()[0];
      int nu = Piu.shape()[1];  
      int nel_ = cPua[0].shape()[1];
#if defined(QMC_COMPLEX)
      int nchol = 2*Luv.shape()[1];  
#else
      int nchol = Luv.shape()[1];  
#endif
      assert(v.shape()[1]==nwalk);
      assert(v.shape()[0]==nchol);  
      using ma::T;
      int c0,cN;
      std::tie(c0,cN) = FairDivideBoundary(comm->rank(),nchol,comm->size());  
      if(haj.size()==1) {
        size_t memory_needs = nwalk*nu + nel_*nu;
        set_shm_buffer(memory_needs);
        boost::multi_array_ref<ComplexType,2> Guu(SM_TMats->data(),extents[nu][nwalk]);
        boost::multi_array_ref<ComplexType,2> T1(SM_TMats->data()+nwalk*nu,extents[nu][nel_]);
        Guu_from_compact(G,Guu,T1);
#if defined(QMC_COMPLEX)
        // reinterpret as RealType matrices with 2x the columns
        boost::multi_array_ref<RealType,2> Luv_R(reinterpret_cast<RealType*>(Luv.origin()),
                                                 extents[Luv.shape()[0]][2*Luv.shape()[1]]);
        boost::multi_array_ref<RealType,2> Guu_R(reinterpret_cast<RealType*>(Guu.origin()),
                                                 extents[nu][2*nwalk]);
        boost::multi_array_ref<RealType,2> v_R(reinterpret_cast<RealType*>(v.origin()),
                                                 extents[v.shape()[0]][2*v.shape()[1]]);
        ma::product(a,T(Luv_R[indices[range_t()][range_t(c0,cN)]]),Guu_R,
                    c,v_R[indices[range_t(c0,cN)][range_t()]]);
#else
        ma::product(a,T(Luv.get()[indices[range_t()][range_t(c0,cN)]]),Guu,
                    c,v[indices[range_t(c0,cN)][range_t()]]);
#endif
      } else {
        size_t memory_needs = nwalk*nu + nmo_*nu;
        set_shm_buffer(memory_needs);
        boost::multi_array_ref<ComplexType,2> Guu(SM_TMats->data(),extents[nu][nwalk]);
        boost::multi_array_ref<ComplexType,2> T1(SM_TMats->data()+nwalk*nu,extents[nmo_][nu]);
        Guu_from_full(G,Guu,T1);
#if defined(QMC_COMPLEX)
        // reinterpret as RealType matrices with 2x the columns
        boost::multi_array_ref<RealType,2> Luv_R(reinterpret_cast<RealType*>(Luv.origin()),
                                                 extents[Luv.shape()[0]][2*Luv.shape()[1]]);
        boost::multi_array_ref<RealType,2> Guu_R(reinterpret_cast<RealType*>(Guu.origin()),
                                                 extents[nu][2*nwalk]);
        boost::multi_array_ref<RealType,2> v_R(reinterpret_cast<RealType*>(v.origin()),
                                                 extents[v.shape()[0]][2*v.shape()[1]]);
        ma::product(a,T(Luv_R[indices[range_t()][range_t(c0,cN)]]),Guu_R,
                    c,v_R[indices[range_t(c0,cN)][range_t()]]);
#else
        ma::product(a,T(Luv.get()[indices[range_t()][range_t(c0,cN)]]),Guu,
                    c,v[indices[range_t(c0,cN)][range_t()]]);
#endif
      }  
      comm->barrier();
*/
    }

    bool distribution_over_cholesky_vectors() const { return false; }
    int number_of_ke_vectors() const{ 
        return 0; 
    }
#if defined(QMC_COMPLEX)
    int local_number_of_cholesky_vectors() const{ 
        return 0; 
    }
    int global_number_of_cholesky_vectors() const{
        return 0; 
    }
#else
    int local_number_of_cholesky_vectors() const{ 
        return 0; 
    }
    int global_number_of_cholesky_vectors() const{
        return 0; 
    }
#endif

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const {return true;} 
    bool transposed_G_for_E() const {return true;} 
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const {return true;} 

    bool fast_ph_energy() const { return true; }

  protected:

/*
    // Guu[nu][nwalk]
    template<class MatA, class MatB, class MatC>
    void Guu_from_compact(MatA const& G, MatB&& Guu, MatC&& T1) {
      int nmo_ = int(Piu.shape()[0]);
      int nu = int(Piu.shape()[1]);
      int nel_ = cPua[0].shape()[1];
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());  
      int nw=G.shape()[0];  

      assert(G.shape()[0] == Guu.shape()[1]);
      assert(G.shape()[1] == nel_*nmo_);
      assert(Guu.shape()[0] == nu); 
      assert(T1.shape()[0] == nu);
      assert(T1.shape()[1] == nel_);

      using ma::transposed;
      comm->barrier();
      ComplexType a = (walker_type==CLOSED)?ComplexType(2.0):ComplexType(1.0);
      for(int iw=0; iw<nw; ++iw) {
        boost::const_multi_array_ref<ComplexType,2> Giw(G[iw].origin(),extents[nel_][nmo_]);
        // transposing inetermediary to make dot products faster in the next step
        ma::product(transposed(Piu.get()[indices[range_t()][range_t(u0,uN)]]),
                  transposed(Giw),
                  T1[indices[range_t(u0,uN)][range_t()]]);
        for(int u=u0; u<uN; ++u)
          Guu[u][iw] = a*ma::dot(cPua[0].get()[u],T1[u]);               
      }
      comm->barrier();
    }  

    // Guu[nu][nwalk]
    template<class MatA, class MatB, class MatC>
    void Guu_from_full(MatA const& G, MatB&& Guu, MatC&& T1) {
      int nmo_ = int(Piu.shape()[0]);
      int nu = int(Piu.shape()[1]);
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int nw=G.shape()[0];

      assert(G.shape()[0] == Guu.shape()[1]);
      assert(Guu.shape()[0] == nu);
      assert(T1.shape()[1] == nu);
      assert(G.shape()[1] == nmo_*nmo_); 
      assert(T1.shape()[0] == nmo_);

      comm->barrier();
      std::fill_n(Guu[u0].origin(),nw*(uN-u0),ComplexType(0.0));  
      ComplexType a = (walker_type==CLOSED)?ComplexType(2.0):ComplexType(1.0);
      for(int iw=0; iw<nw; ++iw) {
        boost::const_multi_array_ref<ComplexType,2> Giw(G[iw].origin(),extents[nmo_][nmo_]);
        ma::product(Giw,Piu.get()[indices[range_t()][range_t(u0,uN)]],
                  T1[indices[range_t()][range_t(u0,uN)]]);
        for(int i=0; i<nmo_; ++i) { 
          auto Ti = T1[i].origin();
          auto Pi = Piu.get()[i].origin();
          for(int u=u0; u<uN; ++u,++Ti,++Pi)
            Guu[u][iw] += a*(*Pi)*(*Ti);
        }
      }
      comm->barrier();
    }  

    // since this is for energy, only compact is accepted
    // Computes Guv and Guu for a single walker
    // As opposed to the other Guu routines, 
    //  this routine expects G for the walker in matrix form
    // rotMuv is partitioned along 'u'
    // G[nel][nmo]
    // Guv[nspin][nu][nu]
    // Guu[u]: summed over spin
    // T1[nel_][nu]
    template<class MatA, class MatB, class MatC, class MatD>
    void Guv_Guu(MatA const& G, MatB&& Guv, MatC&& Guu, MatD&& T1, int k) {

      static_assert(G.dimensionality == 2);
      static_assert(T1.dimensionality == 2);
      static_assert(Guu.dimensionality == 1);
      static_assert(Guv.dimensionality == 3);
      int nspin = (walker_type==COLLINEAR)?2:1;
      int nmo_ = int(rotPiu.shape()[0]);
      int nu = int(rotMuv.shape()[0]);  // potentially distributed over nodes
      int nv = int(rotMuv.shape()[1]);  // not distributed over nodes
      assert(rotPiu.shape()[1] = nv);
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int nu0 = rotMuv.offset()[0];
      ComplexType zero(0.0,0.0);

      assert(Guu.shape()[0] == nv);
      assert(Guv.shape()[1] == nu);
      assert(Guv.shape()[2] == nv);

      // sync first
      comm->barrier();
      if(walker_type==CLOSED || walker_type==NONCOLLINEAR) {
        int nel_ = (walker_type==CLOSED)?NAOA:(NAOA+NAOB);
        assert(Guv.shape()[0] == 1);
        assert(G.shape()[0] == size_t(nel_));
        assert(G.shape()[1] == size_t(nmo_));
        assert(T1.shape()[0] == size_t(nel_));
        assert(T1.shape()[1] == size_t(nv));

        using ma::transposed;
        ma::product(G,rotPiu.get()[indices[range_t()][range_t(v0,vN)]],
                    T1[indices[range_t()][range_t(v0,vN)]]);
        // This operation might benefit from a 2-D work distribution 
        ma::product(rotcPua[k].get()[indices[range_t(nu0,nu0+nu)][range_t()]],
                    T1[indices[range_t()][range_t(v0,vN)]],
                    Guv[0][indices[range_t()][range_t(v0,vN)]]);
        for(int v=v0; v<vN; ++v) 
          if( v < nu0 || v >= nu0+nu ) {
            Guu[v] = ma::dot(rotcPua[k].get()[v],T1[indices[range_t()][v]]);
          } else
            Guu[v] = Guv[0][v-nu0][v];  
      } else {
        int nel_ = NAOA+NAOB;
        assert(Guv.shape()[0] == 2);
        assert(G.shape()[0] == nel_);
        assert(G.shape()[1] == nmo_);
        assert(T1.shape()[0] == nel_);
        assert(T1.shape()[1] == nv);

        using ma::transposed;
        ma::product(G,rotPiu.get()[indices[range_t()][range_t(v0,vN)]],
                    T1[indices[range_t()][range_t(v0,vN)]]);
        // This operation might benefit from a 2-D work distribution 
        // Alpha
        ma::product(rotcPua[k].get()[indices[range_t(nu0,nu0+nu)][range_t(0,NAOA)]],
                    T1[indices[range_t(0,NAOA)][range_t(v0,vN)]],
                    Guv[0][indices[range_t()][range_t(v0,vN)]]);
        ma::product(rotcPua[k].get()[indices[range_t(nu0,nu0+nu)][range_t(NAOA,nel_)]],
                    T1[indices[range_t(NAOA,nel_)][range_t(v0,vN)]],
                    Guv[1][indices[range_t()][range_t(v0,vN)]]);
        for(int v=v0; v<vN; ++v) 
          if( v < nu0 || v >= nu0+nu ) {
            Guu[v] = ma::dot(rotcPua[k].get()[v],T1[indices[range_t()][v]]);
          } else
            Guu[v] = Guv[0][v-nu0][v]+Guv[1][v-nu0][v];
      } 
      comm->barrier();
    }

    // since this is for energy, only compact is accepted
    // Computes Guv and Guu for a single walker
    // As opposed to the other Guu routines, 
    //  this routine expects G for the walker in matrix form
    // rotMuv is partitioned along 'u'
    // G[nel][nmo]
    // Guv[nu][nu]
    // Guu[u]: summed over spin
    // T1[nel_][nu]
    template<class MatA, class MatB, class MatC, class MatD>
    void Guv_Guu2(MatA const& G, MatB&& Guv, MatC&& Guu, MatD&& T1, int k) {

      static_assert(G.dimensionality == 2);
      static_assert(T1.dimensionality == 2);
      static_assert(Guu.dimensionality == 1);
      static_assert(Guv.dimensionality == 2);
      int nmo_ = int(rotPiu.shape()[0]);
      int nu = int(rotMuv.shape()[0]);  // potentially distributed over nodes
      int nv = int(rotMuv.shape()[1]);  // not distributed over nodes
      assert(rotPiu.shape()[1] = nv);
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int nu0 = rotMuv.offset()[0];
      ComplexType zero(0.0,0.0);

      assert(Guu.shape()[0] == nv);
      assert(Guv.shape()[0] == nu);
      assert(Guv.shape()[1] == nv);

      // sync first
      comm->barrier();
      int nel_ = (walker_type==CLOSED)?NAOA:(NAOA+NAOB);
      assert(G.shape()[0] == size_t(nel_));
      assert(G.shape()[1] == size_t(nmo_));
      assert(T1.shape()[0] == size_t(nel_));
      assert(T1.shape()[1] == size_t(nv));

      using ma::transposed;
      ma::product(G,rotPiu.get()[indices[range_t()][range_t(v0,vN)]],
                  T1[indices[range_t()][range_t(v0,vN)]]);
      // This operation might benefit from a 2-D work distribution 
      ma::product(rotcPua[k].get()[indices[range_t(nu0,nu0+nu)][range_t()]],
                  T1[indices[range_t()][range_t(v0,vN)]],
                  Guv[indices[range_t()][range_t(v0,vN)]]);
      for(int v=v0; v<vN; ++v) 
        if( v < nu0 || v >= nu0+nu ) {
          Guu[v] = ma::dot(rotcPua[k].get()[v],T1[indices[range_t()][v]]);
        } else
         Guu[v] = Guv[v-nu0][v];  
      comm->barrier();
    }
*/
  protected:

    communicator* comm;

    bool low_memory;

    WALKER_TYPES walker_type;

    int global_nCV;
    int local_nCV;

    // bare one body hamiltonian
    shmC3Tensor H1;

    // (potentially half rotated) one body hamiltonian
    shmCMatrix haj;

    // number of orbitals per k-point
    std::vector<int> nopk;

    // number of cholesky vectors per Q-point
    std::vector<int> ncholpQ;

    // number of G per Q-point
    std::vector<int> nGpk;

    // number of electrons per k-point
    // nelpk[ndet][nspin*nkpts]
    shmIMatrix nelpk;

    // maps (Q,K) --> k2
    shmIMatrix QKToK2;

    // maps (Q,K) --> G 
    shmIMatrix QKToG;

    /************************************************/
    // Used in the calculation of the energy
    // Coulomb matrix elements of interpolating vectors
    std::vector<shmSpMatrix> rotLQGun;

    // Orbitals at interpolating points
    shmSpMatrix rotPiu;

    // Half-rotated Orbitals at interpolating points
    std::vector<shmSpMatrix> rotcPua;
    /************************************************/

    /************************************************/
    // Following 3 used in calculation of vbias and vHS
    // Cholesky factorization of Muv 
    std::vector<shmSpMatrix> LQGun;

    // Orbitals at interpolating points
    shmSpMatrix Piu;
 
    // Half-rotated Orbitals at interpolating points
    std::vector<shmSpMatrix> cPua;
    /************************************************/

    // Muv for energy
    std::vector<shmSpMatrix> rotMuv;
     
    // one-body piece of Hamiltonian factorization
    shmC3Tensor vn0;

    ValueType E0;

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

    boost::multi::array<ComplexType,3> eloc;

};

}

}

#endif
