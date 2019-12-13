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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_REAL3INDEXFACTORIZATION_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_REAL3INDEXFACTORIZATION_HPP

#include <vector>
#include <type_traits>
#include <random>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "mpi3/shared_communicator.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/taskgroup.h"

namespace qmcplusplus
{

namespace afqmc
{

// Custom implementation for real build
class Real3IndexFactorization
{

  using IVector = boost::multi::array<int,1>;
  using CVector = boost::multi::array<ComplexType,1>;
  using SpVector = boost::multi::array<SPComplexType,1>;

  using CMatrix = boost::multi::array<ComplexType,2>;
  using CMatrix_cref = boost::multi::array_cref<ComplexType,2>;
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2>;
  using CVector_ref = boost::multi::array_ref<ComplexType,1>;

  using RMatrix = boost::multi::array<RealType,2>;
  using RMatrix_cref = boost::multi::array_cref<RealType,2>;
  using RMatrix_ref = boost::multi::array_ref<RealType,2>;
  using RVector_ref = boost::multi::array_ref<RealType,1>;

  using SpCMatrix = boost::multi::array<SPComplexType,2>;
  using SpCMatrix_cref = boost::multi::array_cref<SPComplexType,2>;
  using SpCVector_ref = boost::multi::array_ref<SPComplexType,1>;
  using SpCMatrix_ref = boost::multi::array_ref<SPComplexType,2>;

  using SpRMatrix = boost::multi::array<SPRealType,2>;
  using SpRMatrix_cref = boost::multi::array_cref<SPRealType,2>;
  using SpRVector_ref = boost::multi::array_ref<SPRealType,1>;
  using SpRMatrix_ref = boost::multi::array_ref<SPRealType,2>;

  using C3Tensor = boost::multi::array<ComplexType,3>;
  using SpC3Tensor = boost::multi::array<SPComplexType,3>;
  using SpC3Tensor_ref = boost::multi::array_ref<SPComplexType,3>;
  using SpC4Tensor_ref = boost::multi::array_ref<SPComplexType,4>;

  using shmCVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using shmRMatrix = boost::multi::array<RealType,2,shared_allocator<RealType>>;
  using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using shmC3Tensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;

  using shmSpRVector = boost::multi::array<SPRealType,1,shared_allocator<SPRealType>>;
  using shmSpRMatrix = boost::multi::array<SPRealType,2,shared_allocator<SPRealType>>;
  using shmSpCMatrix = boost::multi::array<SPComplexType,2,shared_allocator<SPComplexType>>;
  using shmSpC3Tensor = boost::multi::array<SPComplexType,3,shared_allocator<SPComplexType>>;

  using this_t = Real3IndexFactorization;

  public:

    Real3IndexFactorization(afqmc::TaskGroup_& tg_,
                 WALKER_TYPES type,
                 shmRMatrix&& hij_,
                 shmCMatrix&& haj_,
                 shmSpRMatrix&& vik,
                 shmSpCMatrix&& vak,
                 std::vector<shmSpC3Tensor>&& vank,
                 shmCMatrix&& vn0_,
                 ValueType e0_,
                 int cv0,
                 int gncv):
        TG(tg_),
        walker_type(type),
        global_origin(cv0),
        global_nCV(gncv),
        local_nCV(0),
        E0(e0_),
        hij(std::move(hij_)),
        haj(std::move(haj_)),
        Likn(std::move(vik)),
        Lank(std::move(vank)),
        Lakn(std::move(vak)),
        vn0(std::move(vn0_)),
        SM_TMats({1,1},shared_allocator<SPComplexType>{TG.TG_local()})
    {
      local_nCV=Likn.size(1);
      TG.Node().barrier();
    }

    ~Real3IndexFactorization() {}

    Real3IndexFactorization(const Real3IndexFactorization& other) = delete;
    Real3IndexFactorization& operator=(const Real3IndexFactorization& other) = delete;
    Real3IndexFactorization(Real3IndexFactorization&& other) = default;
    Real3IndexFactorization& operator=(Real3IndexFactorization&& other) = default;

    CMatrix getOneBodyPropagatorMatrix(TaskGroup_& TG, boost::multi::array<ComplexType,1> const& vMF) 
    {
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

    // KEleft and KEright must be in shared memory for this to work correctly
    template<class Mat, class MatB, class MatC, class MatD>
    void energy(Mat&& E, MatB const& Gc, int nd, MatC* KEleft, MatD* KEright, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      assert(E.size(1)>=3);
      assert(nd >= 0);
      assert(nd < haj.size());
      if(walker_type==COLLINEAR)
        assert(2*nd+1 < Lank.size());
      else
        assert(nd < Lank.size());

      int nwalk = Gc.size(0);
      int nspin = (walker_type==COLLINEAR?2:1);
      int NMO = hij.size(0); 
      int nel[2];
      nel[0] = Lank[nspin*nd].size(0);
      nel[1] = ((nspin==2)?Lank[nspin*nd+1].size(0):0);
      assert(Lank[nspin*nd].size(1) == local_nCV);
      assert(Lank[nspin*nd].size(2) == NMO);
      if(nspin==2) {
        assert(Lank[nspin*nd+1].size(1) == local_nCV);
        assert(Lank[nspin*nd+1].size(2) == NMO);
      }
      assert(Gc.num_elements() == nwalk*(nel[0]+nel[1])*NMO);

      int getKr = KEright!=nullptr;
      int getKl = KEleft!=nullptr;
      if(E.size(0) != nwalk || E.size(1) < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/Real3IndexFactorization::energy(...). Incorrect matrix dimensions \n");

      // T[nwalk][nup][nup][local_nCV] + D[nwalk][nwalk][local_nCV]
      size_t mem_needs(0);
      size_t cnt(0);
      if(addEJ) {
#if MIXED_PRECISION
        mem_needs += nwalk*local_nCV; 
#else
        if(not getKl) mem_needs += nwalk*local_nCV;
#endif
      }
      if(addEXX) {
        mem_needs += nwalk*nel[0]*nel[0]*local_nCV;        
#if MIXED_PRECISION
        mem_needs += nwalk*nel[0]*NMO;
#else
        if(nspin == 2) mem_needs += nwalk*nel[0]*NMO;
#endif
      }
      set_shm_buffer(mem_needs);

      // messy
      SPComplexType *Klptr(nullptr);
      long Knr=0, Knc=0;
      if(addEJ) {
        Knr=nwalk;
        Knc=local_nCV;
        if(getKr) {
          assert(KEright->size(0) == nwalk && KEright->size(1) == local_nCV);
          assert(KEright->stride(0) == KEright->size(1));
        }
#if MIXED_PRECISION
        if(getKl) {
          assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
          assert(KEleft->stride(0) == KEleft->size(1));
        }
#else
        if(getKl) {
          assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
          assert(KEleft->stride(0) == KEleft->size(1));
          Klptr = to_address(KEleft->origin());
        } else 
#endif
        {
          Klptr = to_address(SM_TMats.origin())+cnt;
          cnt += Knr*Knc; 
        }
        if(TG.TG_local().root()) std::fill_n(Klptr,Knr*Knc,SPComplexType(0.0));
      } else if(getKr or getKl) {
        APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
      }
      SpCMatrix_ref Kl(Klptr,{long(Knr),long(Knc)});

      for(int n=0; n<nwalk; n++)
        std::fill_n(E[n].origin(),3,ComplexType(0.));


      // one-body contribution
      // haj[ndet][nocc*nmo]
      // not parallelized for now, since it would require customization of Wfn
      if(addH1) {
        boost::multi::array_cref<ComplexType,1> haj_ref(to_address(haj[nd].origin()), iextensions<1u>{haj[nd].num_elements()});
        ma::product(ComplexType(1.),Gc,haj_ref,ComplexType(1.),E(E.extension(0),0));
        for(int i=0; i<nwalk; i++)
          E[i][0] += E0;
      }

      // move calculation of H1 here
      // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
      //       Not sure how to do it for COLLINEAR.
      if(addEXX) {
        SPRealType scl = (walker_type==CLOSED?2.0:1.0);

        for(int ispin=0, is0=0; ispin<nspin; ispin++) {

          size_t cnt_(cnt);
          SPComplexType *ptr(nullptr);
#if MIXED_PRECISION
          ptr = to_address(SM_TMats.origin())+cnt_;
          cnt_ += nwalk*nel[ispin]*NMO;
          for(int n=0; n<nwalk; ++n) {
            if( n%TG.TG_local().size() != TG.TG_local().rank() ) continue;
            copy_n_cast(to_address(Gc[n].origin())+is0,nel[ispin]*NMO,ptr+n*nel[ispin]*NMO);  
          }
          TG.TG_local().barrier();
#else
          if(nspin==1) {
            ptr = to_address(Gc.origin());
          } else {
            ptr = to_address(SM_TMats.origin())+cnt_;
            cnt_ += nwalk*nel[ispin]*NMO;
            for(int n=0; n<nwalk; ++n) {
              if( n%TG.TG_local().size() != TG.TG_local().rank() ) continue;  
              std::copy_n(to_address(Gc[n].origin())+is0,nel[ispin]*NMO,ptr+n*nel[ispin]*NMO);  
            }
            TG.TG_local().barrier();
          }
#endif
 
          SpCMatrix_ref GF(ptr,{nwalk*nel[ispin],NMO});  
          SpCMatrix_ref Lan(to_address(Lank[nd*nspin + ispin].origin()),
                                                 {nel[ispin]*local_nCV,NMO});
          SpCMatrix_ref Twban(to_address(SM_TMats.origin())+cnt_,{nwalk*nel[ispin],nel[ispin]*local_nCV});
          SpC4Tensor_ref T4Dwban(Twban.origin(),{nwalk,nel[ispin],nel[ispin],local_nCV});

          long i0, iN;
          std::tie(i0,iN) = FairDivideBoundary(long(TG.TG_local().rank()),long(nel[ispin]*local_nCV),
                                               long(TG.TG_local().size()));
          ma::product(GF,ma::T(Lan.sliced(i0,iN)),Twban(Twban.extension(0),{i0,iN}));
          TG.TG_local().barrier();

          for(int n=0, an=0; n<nwalk; ++n) {
            ComplexType E_(0.0);
            for(int a=0; a<nel[ispin]; ++a, an++) {
              if( an%TG.TG_local().size() != TG.TG_local().rank() ) continue;  
              for(int b=0; b<nel[ispin]; ++b)
                E_ += static_cast<ComplexType>(ma::dot(T4Dwban[n][a][b],T4Dwban[n][b][a]));
            }    
            E[n][1] -= 0.5*scl*E_;
          }

          if(addEJ) {
            for(int n=0; n<nwalk; ++n) {
              if( n%TG.TG_local().size() != TG.TG_local().rank() ) continue;  
              for(int a=0; a<nel[ispin]; ++a) 
                ma::axpy(SPComplexType(1.0),T4Dwban[n][a][a],Kl[n]);
            }
          }
          is0 += nel[ispin]*NMO; 

        } // if
      }   
      TG.TG_local().barrier();

      if(addEJ) {
        if(not addEXX) {
          // calculate Kr
          APP_ABORT(" Error: Finish addEJ and not addEXX");
        }
        TG.TG_local().barrier();
        SPRealType scl = (walker_type==CLOSED?2.0:1.0);
        for(int n=0; n<nwalk; ++n) {
          if(n%TG.TG_local().size() == TG.TG_local().rank()) 
            E[n][2] += 0.5*static_cast<ComplexType>(scl*scl*ma::dot(Kl[n],Kl[n]));
        }
#if MIXED_PRECISION
        if(getKl) {
          long i0, iN;
          std::tie(i0,iN) = FairDivideBoundary(long(TG.TG_local().rank()),long(KEleft->num_elements()),
                                               long(TG.TG_local().size()));
          copy_n_cast(Klptr+i0,iN-i0,to_address(KEleft->origin())+i0);
        } 
#endif
        if(getKr) {
          long i0, iN;
          std::tie(i0,iN) = FairDivideBoundary(long(TG.TG_local().rank()),long(KEright->num_elements()),
                                               long(TG.TG_local().size()));
          copy_n_cast(Klptr+i0,iN-i0,to_address(KEright->origin())+i0);
        }
        TG.TG_local().barrier();
      }
    }

    template<class... Args>
    void fast_energy(Args&&... args)
    {
      APP_ABORT(" Error: fast_energy not implemented in Real3IndexFactorization. \n");
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      using BType = typename std::decay<MatB>::type::element ;
      using AType = typename std::decay<MatA>::type::element ;
      boost::multi::array_ref<BType,2> v_(to_address(v.origin()),
                                        {v.size(0),1});
      boost::multi::array_ref<const AType,2> X_(to_address(X.origin()),
                                        {X.size(0),1});
      return vHS(X_,v_,a,c);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(MatA& X, MatB&& v, double a=1., double c=0.) {
      assert( Likn.size(1) == X.size(0) );
      assert( Likn.size(0) == v.size(0) );
      assert( X.size(1) == v.size(1) );
      long ik0, ikN;
      std::tie(ik0,ikN) = FairDivideBoundary(long(TG.TG_local().rank()),long(Likn.size(0)),long(TG.TG_local().size()));
#if MIXED_PRECISION
      size_t mem_needs = X.num_elements()+v.num_elements();
      set_shm_buffer(mem_needs);
      boost::multi::array_ref<SPComplexType,2> vsp(to_address(SM_TMats.origin()), v.extensions());
      boost::multi::array_ref<SPComplexType,2> Xsp(vsp.origin()+vsp.num_elements(), X.extensions());
      long i0, iN;
      std::tie(i0,iN) = FairDivideBoundary(long(TG.TG_local().rank()),long(X.num_elements()),long(TG.TG_local().size()));
      copy_n_cast(to_address(X.origin())+i0,iN-i0,to_address(Xsp.origin())+i0);
      TG.TG_local().barrier();
      ma::product(SPValueType(a),Likn.sliced(ik0,ikN),Xsp,
                  SPValueType(c),vsp.sliced(ik0,ikN));
      copy_n_cast(to_address(vsp[ik0].origin()),vsp.size(1)*(ikN-ik0),
                  to_address(v[ik0].origin()));
#else
      ma::product(SPValueType(a),Likn.sliced(ik0,ikN),X,
                  SPValueType(c),v.sliced(ik0,ikN));
#endif
      TG.TG_local().barrier();
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
      using BType = typename std::decay<MatB>::type::element ;
      using AType = typename std::decay<MatA>::type::element ;
      boost::multi::array_ref<BType,2> v_(to_address(v.origin()),
                                        {v.size(0),1});
      if(haj.size(0) == 1) {

        boost::multi::array_cref<AType,2> G_(to_address(G.origin()),
                                        {1,G.size(0)});
        return vbias(G_,v_,a,c,k);

      } else {  

        boost::multi::array_cref<AType,2> G_(to_address(G.origin()),
                                        {G.size(0),1});
        return vbias(G_,v_,a,c,k);
      }
    }

    // v(n,w) = sum_ak L(ak,n) G(w,ak)
    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int k=0) {
      if(haj.size(0) == 1) {
        assert( Lakn.size(0) == G.size(1) );
        assert( Lakn.size(1) == v.size(0) );
        assert( G.size(0) == v.size(1) );
        long ic0, icN;
        std::tie(ic0,icN) = FairDivideBoundary(long(TG.TG_local().rank()),long(Lakn.size(1)),long(TG.TG_local().size()));

#if MIXED_PRECISION
        size_t mem_needs = G.num_elements()+v.num_elements();
        set_shm_buffer(mem_needs);
        boost::multi::array_ref<SPComplexType,2> vsp(to_address(SM_TMats.origin()), v.extensions());
        boost::multi::array_ref<SPComplexType,2> Gsp(vsp.origin()+vsp.num_elements(), G.extensions());
        long i0, iN;
        std::tie(i0,iN) = FairDivideBoundary(long(TG.TG_local().rank()),long(G.num_elements()),long(TG.TG_local().size()));
        copy_n_cast(to_address(G.origin())+i0,iN-i0,to_address(Gsp.origin())+i0);
        TG.TG_local().barrier();
        if(walker_type==CLOSED) a*=2.0;
        ma::product(SPValueType(a),ma::T(Lakn(Lakn.extension(0),{ic0,icN})),ma::T(Gsp),
                  SPValueType(c),vsp.sliced(ic0,icN));
        copy_n_cast(to_address(vsp[ic0].origin()),vsp.size(1)*(icN-ic0),
                  to_address(v[ic0].origin()));
        TG.TG_local().barrier();
#else
        if(walker_type==CLOSED) a*=2.0;
        ma::product(SPValueType(a),ma::T(Lakn(Lakn.extension(0),{ic0,icN})),ma::T(G),
                  SPValueType(c),v.sliced(ic0,icN));
#endif
      } else {
        // multideterminant is not half-rotated, so use Likn
        assert( Likn.size(0) == G.size(0) );
        assert( Likn.size(1) == v.size(0) );
        assert( G.size(1) == v.size(1) );
        long ic0, icN;
        std::tie(ic0,icN) = FairDivideBoundary(long(TG.TG_local().rank()),long(Likn.size(1)),long(TG.TG_local().size()));

#if MIXED_PRECISION
        size_t mem_needs = G.num_elements()+v.num_elements();
        set_shm_buffer(mem_needs);
        boost::multi::array_ref<SPComplexType,2> vsp(to_address(SM_TMats.origin()), v.extensions());
        boost::multi::array_ref<SPComplexType,2> Gsp(vsp.origin()+vsp.num_elements(), G.extensions());
        long i0, iN;
        std::tie(i0,iN) = FairDivideBoundary(long(TG.TG_local().rank()),long(G.num_elements()),long(TG.TG_local().size()));
        copy_n_cast(to_address(G.origin())+i0,iN-i0,to_address(Gsp.origin())+i0);
        TG.TG_local().barrier();
        if(walker_type==CLOSED) a*=2.0;
        ma::product(SPValueType(a),ma::T(Likn(Likn.extension(0),{ic0,icN})),Gsp,
                  SPValueType(c),vsp.sliced(ic0,icN));
        copy_n_cast(to_address(vsp[ic0].origin()),vsp.size(1)*(icN-ic0),
                  to_address(v[ic0].origin()));
        TG.TG_local().barrier();
#else
        if(walker_type==CLOSED) a*=2.0;
        ma::product(SPValueType(a),ma::T(Likn(Likn.extension(0),{ic0,icN})),G,
                  SPValueType(c),v.sliced(ic0,icN));
#endif
      }
      TG.TG_local().barrier();
    }

    bool distribution_over_cholesky_vectors() const{ return true; }
    int number_of_ke_vectors() const{ return local_nCV; }
    int local_number_of_cholesky_vectors() const{ return local_nCV; }
    int global_number_of_cholesky_vectors() const{ return global_nCV; }
    int global_origin_cholesky_vector() const{ return global_origin; }

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const{ return (haj.size(0) == 1); } 
    bool transposed_G_for_E() const{return true;}
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const{return false;}

    bool fast_ph_energy() const { return false; }

    boost::multi::array<ComplexType,2> getHSPotentials()
    {
      return boost::multi::array<ComplexType,2>{};
    }

  private:

    afqmc::TaskGroup_& TG;

    WALKER_TYPES walker_type;

    int global_origin;
    int global_nCV;
    int local_nCV;

    ValueType E0;

    // bare one body hamiltonian
    shmRMatrix hij;

    // (potentially half rotated) one body hamiltonian
    shmCMatrix haj;

    //Cholesky Tensor Lik[i][k][n]
    shmSpRMatrix Likn;

    // permuted half-tranformed Cholesky tensor
    // Lank[ 2*idet + ispin ]
    std::vector<shmSpC3Tensor> Lank;

    // half-tranformed Cholesky tensor
    // only used in single determinant case, haj.size(0)==1.
    shmSpCMatrix Lakn;

    // one-body piece of Hamiltonian factorization
    shmCMatrix vn0;

    // shared buffer space
    // using matrix since there are issues with vectors
    shmSpCMatrix SM_TMats;

    myTimer Timer;

    void set_shm_buffer(size_t N) {
      if(SM_TMats.num_elements() < N)
        SM_TMats.reextent({N,1});
    }

};

}

}

#endif
