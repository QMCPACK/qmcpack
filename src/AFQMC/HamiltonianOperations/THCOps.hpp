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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_THCOPS_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_THCOPS_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "Utilities/FairDivide.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "mpi3/shared_communicator.hpp"
#include "AFQMC/Matrix/mpi3_shared_ma_proxy.hpp"
#include "type_traits/scalar_traits.h"
#include "AFQMC/Wavefunctions/Excitations.hpp"
#include "AFQMC/Wavefunctions/phmsd_helpers.hpp"

namespace qmcplusplus
{

namespace afqmc
{

template<class T>
class THCOps
{
#if defined(MIXED_PRECISION)
  using SpT = typename to_single_precision<T>::value_type;
  using SpC = typename to_single_precision<ComplexType>::value_type;
#else
  using SpT = T;
  using SpC = ComplexType;
#endif

  using TVector = boost::multi::array<T,1>;
  using SpTVector = boost::multi::array<SpT,1>;
  using CVector = boost::multi::array<ComplexType,1>;
  using CMatrix = boost::multi::array<ComplexType,2>;
  using TMatrix = boost::multi::array<T,2>;
  using shmCMatrix = mpi3_shared_ma_proxy<ComplexType>;
  using shmVMatrix = mpi3_shared_ma_proxy<T>;
  using communicator = boost::mpi3::shared_communicator;
  using shmSpMatrix = boost::multi::array<SPComplexType,2,shared_allocator<SPComplexType>>;

  public:

    static const HamiltonianTypes HamOpType = THC;
    HamiltonianTypes getHamType() const { return HamOpType; }

    /*
     * NAOA/NAOB stands for number of active orbitals alpha/beta (instead of active electrons)
     */
    THCOps(communicator& c_,
           int nmo_, int naoa_, int naob_,
           WALKER_TYPES type,
           CMatrix&& hij_,
           std::vector<CVector>&& h1,
           shmVMatrix&& rotmuv_,
           shmCMatrix&& rotpiu_,
           std::vector<shmCMatrix>&& rotpau_,
           shmVMatrix&& luv_,
           shmCMatrix&& piu_,
           std::vector<shmCMatrix>&& pau_,
           CMatrix&& v0_,
           ValueType e0_,
           bool verbose=false ):
                comm(std::addressof(c_)),
                NMO(nmo_),NAOA(naoa_),NAOB(naob_),
                walker_type(type),
                hij(std::move(hij_)),
                haj(std::move(h1)),
                rotMuv(std::move(rotmuv_)),
                rotPiu(std::move(rotpiu_)),
                rotcPua(std::move(rotpau_)),
                Luv(std::move(luv_)),
                Piu(std::move(piu_)),
                cPua(std::move(pau_)),
                v0(std::move(v0_)),
                E0(e0_),
                SM_TMats({1,1},shared_allocator<SPComplexType>{c_})
    {
/*
for(int i=0; i<haj[0].size(0); i++)
  std::cout<<i <<" " <<haj[0][i]  <<"\n";
std::cout<<"H1:\n";
for(int i=0; i<hij.size(0); i++)
  for(int j=0; j<hij.size(1); j++)
    std::cout<<i <<" " <<j <<" " <<hij[i][j]  <<"\n";
std::cout<<"\n";
*/
      if(haj.size() > 1)
	APP_ABORT(" Error: THC not yet implemented for multiple references.\n");
      assert(comm);
      // current partition over 'u' for L/Piu
      assert(Luv.size(0) == Piu.size(1));
      for(int i=0; i<rotcPua.size(); i++) {
        // rot Ps are not yet distributed
        assert(rotcPua[i].size(0) == rotPiu.size(1));
        if(walker_type==CLOSED)
          assert(rotcPua[i].size(1)==NAOA);
        else if(walker_type==COLLINEAR)
          assert(rotcPua[i].size(1)==NAOA+NAOB);
        else if(walker_type==NONCOLLINEAR)
          assert(rotcPua[i].size(1)==NAOA+NAOB);
      }
      for(int i=0; i<cPua.size(); i++) {
        assert(cPua[i].size(0)==Luv.size(0));
        if(walker_type==CLOSED)
          assert(cPua[i].size(1)==NAOA);
        else if(walker_type==COLLINEAR)
          assert(cPua[i].size(1)==NAOA+NAOB);
        else if(walker_type==NONCOLLINEAR)
          assert(cPua[i].size(1)==NAOA+NAOB);
      }
      if(walker_type==NONCOLLINEAR) {
        assert(Piu.size(0)==2*NMO);
        assert(rotPiu.size(0)==2*NMO);
      } else {
        assert(Piu.size(0)==NMO);
        assert(rotPiu.size(0)==NMO);
      }
    }

    ~THCOps() {}

    THCOps(THCOps const& other) = delete;
    THCOps& operator=(THCOps const& other) = delete;

    THCOps(THCOps&& other) = default;
    THCOps& operator=(THCOps&& other) = default;

    CMatrix getOneBodyPropagatorMatrix(TaskGroup_& TG, CVector const& vMF) {
      int NMO = hij.size(0);
      // in non-collinear case with SO, keep SO matrix here and add it
      // for now, stay collinear
      CMatrix H1({NMO,NMO});

      // add sum_n vMF*Spvn, vMF has local contribution only!
      boost::multi::array_ref<ComplexType,1> H1D(H1.origin(),iextensions<1u>{NMO*NMO});
      std::fill_n(H1D.origin(),H1D.num_elements(),ComplexType(0));
      vHS(vMF, H1D);
      TG.TG().all_reduce_in_place_n(H1D.origin(),H1D.num_elements(),std::plus<>());

      // add hij + v0 and symmetrize
      using ma::conj;
      for(int i=0; i<NMO; i++) {
        H1[i][i] += hij[i][i] + v0[i][i];
        for(int j=i+1; j<NMO; j++) {
          H1[i][j] += hij[i][j] + v0[i][j];
          H1[j][i] += hij[j][i] + v0[j][i];
          if( std::abs( H1[i][j] - ma::conj(H1[j][i]) ) > 1e-8 ) {
            app_error()<<" Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
            app_error()<<i <<" " <<j <<" " <<H1[i][j] <<" " <<H1[j][i] <<" "
                       <<H1[i][j]-(hij[i][j] + v0[i][j]) <<" " <<H1[j][i]-(hij[j][i] + v0[j][i]) <<" "
                       <<hij[i][j] <<" " <<hij[j][i] <<" "
                       <<v0[i][j] <<" " <<v0[j][i] <<std::endl;
            APP_ABORT("Error in getOneBodyPropagatorMatrix. H1 is not hermitian. \n");
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
    void energy(Mat&& E, MatB const& G, int k, MatC* Kl, MatD* Kr, bool addH1=true, bool addEJ=true, bool addEXX=true) {
      if(k>0)
	APP_ABORT(" Error: THC not yet implemented for multiple references.\n");
      // G[nel][nmo]
      //static_assert(E.dimensionality==2);
      //static_assert(G.dimensionality==2);
      assert(E.size(0) == G.size(0));
      assert(E.size(1) == 3);
      int nwalk = G.size(0);
      int getKr = Kr!=nullptr;
      int getKl = Kl!=nullptr;

      // addH1
      std::fill_n(E.origin(),E.num_elements(),ComplexType(0.0));
      if(addH1) {
        ma::product(ComplexType(1.0),G,haj[k],ComplexType(0.0),E(E.extension(0),0));
        for(int i=0; i<nwalk; i++) E[i][0] += E0;
      }
      if(not (addEJ || addEXX)) return;

      int nmo_ = rotPiu.size(0);
      int nu = rotMuv.size(0);
      int nu0 = rotMuv.global_offset()[0];
      int nv = rotMuv.size(1);
      int nel_ = rotcPua[0].size(1);
      int nspin = (walker_type==COLLINEAR)?2:1;
      assert(G.size(1) == nel_*nmo_);
      if(addEJ and getKl)
        assert(Kl->size(0) == nwalk && Kl->size(1) == nu);
      if(addEJ and getKr)
        assert(Kr->size(0) == nwalk && Kr->size(1) == nu);
      using ma::T;
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      // right now the algorithm uses 2 copies of matrices of size nuxnv in COLLINEAR case,
      // consider moving loop over spin to avoid storing the second copy which is not used
      // simultaneously
      size_t memory_needs = nspin*nu*nv + nv + nu  + nel_*(nv+nu);
      set_shm_buffer(memory_needs);
      size_t cnt=0;
      // Guv[nspin][nu][nv]
      boost::multi::array_ref<ComplexType,3> Guv(to_address(SM_TMats.origin()),{nspin,nu,nv});
      cnt+=Guv.num_elements();
      // Guu[u]: summed over spin
      boost::multi::array_ref<ComplexType,1> Guu(to_address(SM_TMats.origin())+cnt,iextensions<1u>{nv});
      cnt+=Guu.num_elements();
      // T1[nel_][nv]
      boost::multi::array_ref<ComplexType,2> T1(to_address(SM_TMats.origin())+cnt,{nel_,nv});
      cnt+=T1.num_elements();
      boost::multi::array_ref<ComplexType,1> Tuu(to_address(SM_TMats.origin())+cnt,iextensions<1u>{nu});
      cnt+=Tuu.num_elements();

      auto&& M_(rotMuv.get());
      int bsz = 256;
      int nbu = ((uN-u0) + bsz - 1) / bsz;
      int nbv = (nv + bsz - 1) / bsz;

      if(walker_type==CLOSED || walker_type==NONCOLLINEAR) {
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        for(int wi=0; wi<nwalk; wi++) {
          boost::multi::array_cref<ComplexType,2> Gw(to_address(G[wi].origin()),{nel_,nmo_});
          boost::multi::array_cref<ComplexType,1> G1D(to_address(G[wi].origin()),iextensions<1u>{nel_*nmo_});
          // need a new routine if addEXX is false,
          // otherwise it is quite inefficient to get Ej only
          Guv_Guu(Gw,Guv,Guu,T1,k);
          if(addEJ) {
            ma::product(rotMuv.get().sliced(u0,uN),Guu,
                        Tuu.sliced(u0,uN));
            if(getKl)
              std::copy_n(to_address(Guu.origin())+nu0+u0,uN-u0,to_address((*Kl)[wi].origin())+u0);
            if(getKr)
              std::copy_n(to_address(Tuu.origin())+u0,uN-u0,to_address((*Kr)[wi].origin())+u0);
            E[wi][2] = 0.5*scl*scl*ma::dot(Guu.sliced(nu0+u0,nu0+uN),Tuu.sliced(u0,uN));
          }
          if(addEXX) {
//Timer.reset("T0");
//Timer.start("T0");
            ComplexType E_(0.0);
            for(int bu=0; bu<nbu; ++bu) {
              int i0 = u0+bu*bsz;
              int iN = std::min(u0+(bu+1)*bsz,uN);
              for(int bv=0; bv<nbv; ++bv) {
                int j0 = bv*bsz;
                int jN = std::min((bv+1)*bsz,nv);
                for(int i=i0; i<iN; ++i) {
                  for(int j=j0; j<jN; ++j)
                    E_ += Guv[0][i][j] * M_[i][j] * Guv[0][j][i];
                }
              }
            }
            E[wi][1] = -0.5*scl*E_;
          }
        }
      } else {
        for(int wi=0; wi<nwalk; wi++) {
          boost::multi::array_cref<ComplexType,2> Gw(to_address(G[wi].origin()),{nel_,nmo_});
          boost::multi::array_cref<ComplexType,1> G1DA(to_address(G[wi].origin()),iextensions<1u>{NAOA*nmo_});
          boost::multi::array_cref<ComplexType,1> G1DB(to_address(G[wi].origin())+NAOA*nmo_,iextensions<1u>{NAOB*nmo_});
          Guv_Guu(Gw,Guv,Guu,T1,k);
          // move calculation of Guv/Guu here to avoid storing 2 copies of Guv for alpha/beta
          if(addEJ) {
            ma::product(rotMuv.get().sliced(u0,uN),Guu,
                      Tuu.sliced(u0,uN));
            if(getKl)
              std::copy_n(to_address(Guu.origin())+nu0+u0,uN,to_address((*Kl)[wi].origin())+u0);
            if(getKr)
              std::copy_n(to_address(Tuu.origin())+u0,uN,to_address((*Kr)[wi].origin())+u0);
            E[wi][2] = 0.5*ma::dot(Guu.sliced(nu0+u0,nu0+uN),Tuu.sliced(u0,uN));
          }
          if(addEXX) {
            // alpha
            ComplexType E_(0.0);
            for(int bu=0; bu<nbu; ++bu) {
              int i0 = u0+bu*bsz;
              int iN = std::min(u0+(bu+1)*bsz,uN);
              for(int bv=0; bv<nbv; ++bv) {
                int j0 = bv*bsz;
                int jN = std::min((bv+1)*bsz,nv);
                for(int i=i0; i<iN; ++i) {
                  for(int j=j0; j<jN; ++j)
                    E_ += Guv[0][i][j] * M_[i][j] * Guv[0][j][i];
                }
              }
            }
            for(int bu=0; bu<nbu; ++bu) {
              int i0 = u0+bu*bsz;
              int iN = std::min(u0+(bu+1)*bsz,uN);
              for(int bv=0; bv<nbv; ++bv) {
                int j0 = bv*bsz;
                int jN = std::min((bv+1)*bsz,nv);
                for(int i=i0; i<iN; ++i) {
                  for(int j=j0; j<jN; ++j)
                    E_ += Guv[1][i][j] * M_[i][j] * Guv[1][j][i];
                }
              }
            }
            E[wi][1] = -0.5*E_;
          }
        }
      }
      comm->barrier();
    }

    template<class MatE, class MatO, class MatG, class MatQ, class MatB,
             class index_aos>
    void fast_energy(MatE&& E, MatO&& Ov, MatG const& GrefA, MatG const& GrefB,
                     MatQ const& QQ0A, MatQ const& QQ0B, MatB&& Qwork,
                     ph_excitations<int,ComplexType> const& abij,
                     std::array<index_aos,2> const& det_couplings)
    {
      if(haj.size() != 1)
        APP_ABORT(" Error: Single reference implementation currently in THCOps::fast_energy.\n");
      if(walker_type!=CLOSED)
        APP_ABORT(" Error: THCOps::fast_energy requires walker_type==CLOSED.\n");
      /*
       * E[nspins][maxn_unique_confg][nwalk][3]
       * Ov[nspins][maxn_unique_confg][nwalk]
       * GrefA[nwalk][NAOA][NMO]
       * GrefB[nwalk][NAOB][NMO]
       * QQ0A[nwalk][NAOA][NAEA]
       * QQ0B[nwalk][NAOA][NAEA]
       */
      static_assert(std::decay<MatE>::type::dimensionality==4, "Wrong dimensionality");
      static_assert(std::decay<MatO>::type::dimensionality==3, "Wrong dimensionality");
      static_assert(std::decay<MatG>::type::dimensionality==3, "Wrong dimensionality");
      static_assert(std::decay<MatQ>::type::dimensionality==3, "Wrong dimensionality");
      //static_assert(std::decay<MatB>::type::dimensionality==3, "Wrong dimensionality");
      int nspin = E.size(0);
      int nrefs = haj.size();
      int nwalk = GrefA.size(0);
      int naoa_ = QQ0A.size(1);
      int naob_ = QQ0B.size(1);
      int nmo_ = rotPiu.size(0);
      int nu = rotMuv.size(0);
      int nu0 = rotMuv.global_offset()[0];
      int nv = rotMuv.size(1);
      int nel_ = rotcPua[0].size(1);
      // checking
      assert(E.size(2) == nwalk);
      assert(E.size(3) == 3);
      assert(Ov.size(0) == nspin);
      assert(Ov.size(1) == E.size(1));
      assert(Ov.size(2) == nwalk);
      assert(GrefA.size(1) == naoa_);
      assert(GrefA.size(2) == nmo_);
      assert(GrefB.size(0) == nwalk);
      assert(GrefB.size(1) == naob_);
      assert(GrefB.size(2) == nmo_);
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
      boost::multi::array_ref<ComplexType,2> Guv(to_address(SM_TMats.origin()),{nu,nv});
      cnt+=Guv.num_elements();
      // Gvv[v]: summed over spin
      boost::multi::array_ref<ComplexType,1> Gvv(to_address(SM_TMats.origin())+cnt,iextensions<1u>{nv});
      cnt+=Gvv.num_elements();
      // S[nel_][nv]
      boost::multi::array_ref<ComplexType,2> Scu(to_address(SM_TMats.origin())+cnt,{nel_,nv});
      cnt+=Scu.num_elements();
      // Qub[nu][nel_]:
      boost::multi::array_ref<ComplexType,2> Qub(to_address(SM_TMats.origin())+cnt,{nu,nel_});
      cnt+=Qub.num_elements();
      boost::multi::array_ref<ComplexType,1> Tuu(to_address(SM_TMats.origin())+cnt,iextensions<1u>{nu});
      cnt+=Tuu.num_elements();
      boost::multi::array_ref<ComplexType,2> Jcb(to_address(SM_TMats.origin())+cnt,{nel_,nel_});
      cnt+=Jcb.num_elements();
      boost::multi::array_ref<ComplexType,2> Xcb(to_address(SM_TMats.origin())+cnt,{nel_,nel_});
      cnt+=Xcb.num_elements();
      boost::multi::array_ref<ComplexType,2> Tub(to_address(SM_TMats.origin())+cnt,{nu,nel_});
      cnt+=Tub.num_elements();
      assert(cnt <= memory_needs);
      if(eloc.size(0) != 2 || eloc.size(1) != nwalk || eloc.size(2) != 3)
        eloc.reextent({2,nwalk,3});

      std::fill_n(eloc.origin(),eloc.num_elements(),ComplexType(0.0));

      RealType scl = (walker_type==CLOSED?2.0:1.0);
      if(comm->root()) {
        std::fill_n(to_address(E.origin()),E.num_elements(),ComplexType(0.0));
        std::fill_n(to_address(Ov[0][1].origin()),nwalk*(Ov.size(1)-1),ComplexType(0.0));
        std::fill_n(to_address(Ov[1][1].origin()),nwalk*(Ov.size(1)-1),ComplexType(0.0));
        auto Ea = E[0][0];
        auto Eb = E[1][0];
        boost::multi::array_cref<ComplexType,2> G2DA(to_address(GrefA.origin()),
                                          {nwalk,GrefA[0].num_elements()});
        ma::product(ComplexType(1.0),G2DA,haj[0],ComplexType(0.0),Ea(Ea.extension(0),0));
        boost::multi::array_cref<ComplexType,2> G2DB(to_address(GrefA.origin()),
                                          {nwalk,GrefA[0].num_elements()});
        ma::product(ComplexType(1.0),G2DB,haj[0],ComplexType(0.0),Eb(Eb.extension(0),0));
        for(int i=0; i<nwalk; i++) {
            Ea[i][0] += E0;
            Eb[i][0] += E0;
        }
      }

      for(int wi=0; wi<nwalk; wi++) {

        { // Alpha
          auto Gw = GrefA[wi];
          boost::multi::array_cref<ComplexType,1> G1D(to_address(Gw.origin()),
                                                        iextensions<1u>{Gw.num_elements()});
          Guv_Guu2(Gw,Guv,Gvv,Scu,0);
          if(u0!=uN)
            ma::product(rotMuv.get().sliced(u0,uN),Gvv,
                      Tuu.sliced(u0,uN));
          auto Mptr = rotMuv.get()[u0].origin();
          auto Gptr = to_address(Guv[u0].origin());
          for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
            (*Gptr) *= (*Mptr);
          if(u0!=uN)
            ma::product(Guv.sliced(u0,uN),rotcPua[0].get(),
                      Qub.sliced(u0,uN));
          comm->barrier();
          if(k0!=kN)
            ma::product(Scu.sliced(k0,kN),Qub,
                      Xcb.sliced(k0,kN));
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
            ma::product(Scu.sliced(k0,kN),Tub,
                      Jcb.sliced(k0,kN));
          for(int c=k0; c<kN; ++c)
            eloc[0][wi][1] += -0.5*scl*Xcb[c][c];
          for(int c=k0; c<kN; ++c)
            eloc[0][wi][2] += 0.5*scl*scl*Jcb[c][c];
/*
          calculate_ph_energies(0,comm->rank(),comm->size(),
                                E[0],Ov[0],QQ0A,Qwork,
                                rotMuv,
                                abij,det_couplings);
*/
        }

        { // Beta: Unnecessary in CLOSED walker type (on Walker)
          auto Gw = GrefB[wi];
          boost::multi::array_cref<ComplexType,1> G1D(to_address(Gw.origin()),
                                                        iextensions<1u>{Gw.num_elements()});
          Guv_Guu2(Gw,Guv,Gvv,Scu,0);
          if(u0!=uN)
            ma::product(rotMuv.get().sliced(u0,uN),Gvv,
                      Tuu.sliced(u0,uN));
          auto Mptr = rotMuv.get()[u0].origin();
          auto Gptr = to_address(Guv[u0].origin());
          for(size_t k=0, kend=(uN-u0)*nv; k<kend; ++k, ++Gptr, ++Mptr)
            (*Gptr) *= (*Mptr);
          if(u0!=uN)
            ma::product(Guv.sliced(u0,uN),rotcPua[0].get(),
                      Qub.sliced(u0,uN));
          comm->barrier();
          if(k0!=kN)
            ma::product(Scu.sliced(k0,kN),Qub,
                      Xcb.sliced(k0,kN));
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
            ma::product(Scu.sliced(k0,kN),Tub,
                      Jcb.sliced(k0,kN));
          for(int c=k0; c<kN; ++c)
            eloc[1][wi][1] += -0.5*scl*Xcb[c][c];
          for(int c=k0; c<kN; ++c)
            eloc[1][wi][2] += 0.5*scl*scl*Jcb[c][c];
        }

      }
      comm->reduce_in_place_n(eloc.origin(),eloc.num_elements(),std::plus<>(),0);
      if(comm->root()) {
        // add Eref contributions to all configurations
        for(int nd=0; nd<E.size(1); ++nd) {
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
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vHS(MatA & X, MatB&& v, double a=1., double c=0.) {
        boost::multi::array_cref<ComplexType,2> X_(to_address(X.origin()),{X.size(0),1});
        boost::multi::array_ref<ComplexType,2> v_(to_address(v.origin()),{1,v.size(0)});
        vHS(X_,v_,a,c);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vHS(MatA & X, MatB&& v, double a=1., double c=0.) {
      int nwalk = X.size(1);
#if defined(QMC_COMPLEX)
      int nchol = 2*Luv.size(1);
#else
      int nchol = Luv.size(1);
#endif
      int nmo_ = Piu.size(0);
      int nu = Piu.size(1);
      assert(Luv.size(0)==nu);
      assert(X.size(0)==nchol);
      assert(v.size(0)==nwalk);
      assert(v.size(1)==nmo_*nmo_);
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
      boost::multi::array_ref<ComplexType,2> Tuw(to_address(SM_TMats.origin()),{nu,nwalk});
      // O[nwalk * nmu * nmu]
#if defined(QMC_COMPLEX)
      // reinterpret as RealType matrices with 2x the columns
      boost::multi::array_ref<RealType,2> Luv_R(reinterpret_cast<RealType*>(Luv.origin()),
                                                 {Luv.size(0),2*Luv.size(1)});
      boost::multi::array_cref<RealType,2> X_R(reinterpret_cast<RealType const*>(to_address(X.origin())),
                                                 {X.size(0),2*X.size(1)});
      boost::multi::array_ref<RealType,2> Tuw_R(reinterpret_cast<RealType*>(Tuw.origin()),
                                                 {nu,2*nwalk});
      ma::product(Luv_R.sliced(u0,uN),X_R,
                  Tuw_R.sliced(u0,uN));
#else
      ma::product(Luv.get().sliced(u0,uN),X,
                  Tuw.sliced(u0,uN));
#endif
      comm->barrier();
#if defined(LOW_MEMORY)
      boost::multi::array_ref<ComplexType,2> Qiu(to_address(SM_TMats.origin())+nwalk*nu,{nmo_,nu});
      for(int wi=0; wi<nwalk; wi++) {
        // Qiu[i][u] = T[u][wi] * conj(Piu[i][u])
        // v[wi][ik] = sum_u Qiu[i][u] * Piu[k][u]
        // O[nmo * nmu]
        for(int i=k0; i<kN; i++) {
          auto p_ = Piu.get()[i].origin();
          for(int u=0; u<nu; u++, ++p_)
            Qiu[i][u] = Tuw[u][wi]*ma::conj(*p_);
        }
        boost::multi::array_ref<ComplexType,2> v_(to_address(v[wi].origin()),{nmo_,nmo_});
        // this can benefit significantly from 2-D partition of work
        // O[nmo * nmo * nmu]
        ma::product(a,Qiu.sliced(k0,kN),T(Piu.get()),
                    c,v_.sliced(k0,kN));
      }
#else
      boost::multi::array_ref<ComplexType,2> Qiu(to_address(SM_TMats.origin())+nwalk*nu,{nwalk*nmo_,nu});
      boost::multi::array_ref<ComplexType,3> Qwiu(to_address(SM_TMats.origin())+nwalk*nu,{nwalk,nmo_,nu});
      // Qiu[i][u] = T[u][wi] * conj(Piu[i][u])
      // v[wi][ik] = sum_u Qiu[i][u] * Piu[k][u]
      // O[nmo * nmu]
      for(int wi=0; wi<nwalk; wi++)
        for(int i=k0; i<kN; i++) {
          auto p_ = Piu.get()[i].origin();
          for(int u=0; u<nu; u++, ++p_)
            Qwiu[wi][i][u] = Tuw[u][wi]*ma::conj(*p_);
        }
      boost::multi::array_ref<ComplexType,2> v_(to_address(v.origin()),{nwalk*nmo_,nmo_});
      // this can benefit significantly from 2-D partition of work
      // O[nmo * nmo * nmu]
      ma::product(a,Qiu.sliced(wk0,wkN),T(Piu.get()),
                  c,v_.sliced(wk0,wkN));
#endif
      comm->barrier();
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==1)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
             typename = void
            >
    void vbias(MatA const& G, MatB&& v, double a=1., double c=0., int k=0) {
        boost::multi::array_cref<ComplexType,2> G_(to_address(G.origin()),{1,G.size(0)});
        boost::multi::array_ref<ComplexType,2> v_(to_address(v.origin()),{v.size(0),1});
        vbias(G_,v_,a,c,k);
    }

    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==2)>
            >
    void vbias(MatA const& G, MatB&& v, double a=1., double c=0., int k=0) {
      if(k>0)
	APP_ABORT(" Error: THC not yet implemented for multiple references.\n");
      int nwalk = G.size(0);
      int nmo_ = Piu.size(0);
      int nu = Piu.size(1);
      int nel_ = cPua[0].size(1);
#if defined(QMC_COMPLEX)
      int nchol = 2*Luv.size(1);
#else
      int nchol = Luv.size(1);
#endif
      assert(v.size(1)==nwalk);
      assert(v.size(0)==nchol);
      using ma::T;
      int c0,cN;
      std::tie(c0,cN) = FairDivideBoundary(comm->rank(),nchol,comm->size());
      if(haj.size()==1) {
        size_t memory_needs = nwalk*nu + nel_*nu;
        set_shm_buffer(memory_needs);
        boost::multi::array_ref<ComplexType,2> Guu(to_address(SM_TMats.origin()),{nu,nwalk});
        boost::multi::array_ref<ComplexType,2> T1(to_address(SM_TMats.origin())+nwalk*nu,{nu,nel_});
        Guu_from_compact(G,Guu,T1);
#if defined(QMC_COMPLEX)
        // reinterpret as RealType matrices with 2x the columns
        boost::multi::array_ref<RealType,2> Luv_R(reinterpret_cast<RealType*>(Luv.origin()),
                                                 {Luv.size(0),2*Luv.size(1)});
        boost::multi::array_ref<RealType,2> Guu_R(reinterpret_cast<RealType*>(Guu.origin()),
                                                 {nu,2*nwalk});
        boost::multi::array_ref<RealType,2> v_R(reinterpret_cast<RealType*>(to_address(v.origin())),
                                                 {v.size(0),2*v.size(1)});
        ma::product(a,T(Luv_R(Luv_R.extension(0),{c0,cN})),Guu_R,
                    c,v_R.sliced(c0,cN));
#else
        ma::product(a,T(Luv.get()(Luv.get().extension(0),{c0,cN})),Guu,
                    c,v.sliced(c0,cN));
#endif
      } else {
        size_t memory_needs = nwalk*nu + nmo_*nu;
        set_shm_buffer(memory_needs);
        boost::multi::array_ref<ComplexType,2> Guu(to_address(SM_TMats.origin()),{nu,nwalk});
        boost::multi::array_ref<ComplexType,2> T1(to_address(SM_TMats.origin())+nwalk*nu,{nmo_,nu});
        Guu_from_full(G,Guu,T1);
#if defined(QMC_COMPLEX)
        // reinterpret as RealType matrices with 2x the columns
        boost::multi::array_ref<RealType,2> Luv_R(reinterpret_cast<RealType*>(Luv.origin()),
                                                 {Luv.size(0),2*Luv.size(1)});
        boost::multi::array_ref<RealType,2> Guu_R(reinterpret_cast<RealType*>(Guu.origin()),
                                                 {nu,2*nwalk});
        boost::multi::array_ref<RealType,2> v_R(reinterpret_cast<RealType*>(to_address(v.origin())),
                                                 {v.size(0),2*v.size(1)});
        ma::product(a,T(Luv_R(Luv_R.extension(0),{c0,cN})),Guu_R,
                    c,v_R.sliced(c0,cN));
#else
        ma::product(a,T(Luv.get()(Luv.get().extension(0),{c0,cN})),Guu,
                    c,v.sliced(c0,cN));
#endif
      }
      comm->barrier();
    }

    template<class Mat, class MatB>
    void generalizedFockMatrix(Mat&& G, MatB&& Fp, MatB&& Fm)
    {
      APP_ABORT(" Error: generalizedFockMatrix not implemented for this hamiltonian.\n"); 
    }

    bool distribution_over_cholesky_vectors() const { return false; }
    int number_of_ke_vectors() const{
        return rotMuv.size(0);
    }
#if defined(QMC_COMPLEX)
    int local_number_of_cholesky_vectors() const{
        return 2*Luv.size(1);
    }
    int global_number_of_cholesky_vectors() const{
        return 2*Luv.size(1);
    }
#else
    int local_number_of_cholesky_vectors() const{
        return Luv.size(1);
    }
    int global_number_of_cholesky_vectors() const{
        return Luv.size(1);
    }
#endif
    int global_origin_cholesky_vector() const{
        return 0;
    }

    // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
    bool transposed_G_for_vbias() const {return true;}
    bool transposed_G_for_E() const {return true;}
    // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
    bool transposed_vHS() const {return true;}

    bool fast_ph_energy() const { return true; }

    boost::multi::array<ComplexType,2> getHSPotentials()
    {
      return boost::multi::array<ComplexType,2>{};
    }

  protected:

    // Guu[nu][nwalk]
    template<class MatA, class MatB, class MatC>
    void Guu_from_compact(MatA const& G, MatB&& Guu, MatC&& T1) {
      int nmo_ = int(Piu.size(0));
      int nu = int(Piu.size(1));
      int nel_ = cPua[0].size(1);
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int nw=G.size(0);

      assert(G.size(0) == Guu.size(1));
      assert(G.size(1) == nel_*nmo_);
      assert(Guu.size(0) == nu);
      assert(T1.size(0) == nu);
      assert(T1.size(1) == nel_);

      using ma::transposed;
      comm->barrier();
      ComplexType a = (walker_type==CLOSED)?ComplexType(2.0):ComplexType(1.0);
      for(int iw=0; iw<nw; ++iw) {
        boost::multi::array_cref<ComplexType,2> Giw(to_address(G[iw].origin()),{nel_,nmo_});
        // transposing inetermediary to make dot products faster in the next step
        ma::product(transposed(Piu.get()({0,nmo_},{u0,uN})),
                  transposed(Giw),
                  T1.sliced(u0,uN));
        for(int u=u0; u<uN; ++u)
          Guu[u][iw] = a*ma::dot(cPua[0].get()[u],T1[u]);
      }
      comm->barrier();
    }

    // Guu[nu][nwalk]
    template<class MatA, class MatB, class MatC>
    void Guu_from_full(MatA const& G, MatB&& Guu, MatC&& T1) {
      int nmo_ = int(Piu.size(0));
      int nu = int(Piu.size(1));
      int u0,uN;
      std::tie(u0,uN) = FairDivideBoundary(comm->rank(),nu,comm->size());
      int nw=G.size(0);

      assert(G.size(0) == Guu.size(1));
      assert(Guu.size(0) == nu);
      assert(T1.size(1) == nu);
      assert(G.size(1) == nmo_*nmo_);
      assert(T1.size(0) == nmo_);

      comm->barrier();
      std::fill_n(Guu[u0].origin(),nw*(uN-u0),ComplexType(0.0));
      ComplexType a = (walker_type==CLOSED)?ComplexType(2.0):ComplexType(1.0);
      for(int iw=0; iw<nw; ++iw) {
        boost::multi::array_cref<ComplexType,2> Giw(to_address(G[iw].origin()),{nmo_,nmo_});
        ma::product(Giw,Piu.get()({0,nmo_},{u0,uN}),
                  T1(T1.extension(0),{u0,uN}));
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

      static_assert(std::decay<MatA>::type::dimensionality == 2, "Wrong dimensionality");
      static_assert(std::decay<MatB>::type::dimensionality == 3, "Wrong dimensionality");
      static_assert(std::decay<MatC>::type::dimensionality == 1, "Wrong dimensionality");
      static_assert(std::decay<MatD>::type::dimensionality == 2, "Wrong dimensionality");
      int nmo_ = int(rotPiu.size(0));
      int nu = int(rotMuv.size(0));  // potentially distributed over nodes
      int nv = int(rotMuv.size(1));  // not distributed over nodes
      assert(rotPiu.size(1) == nv);
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int nu0 = rotMuv.global_offset()[0];
      ComplexType zero(0.0,0.0);

      assert(Guu.size(0) == nv);
      assert(Guv.size(1) == nu);
      assert(Guv.size(2) == nv);

      // sync first
      comm->barrier();
      if(walker_type==CLOSED || walker_type==NONCOLLINEAR) {
        int nel_ = (walker_type==CLOSED)?NAOA:(NAOA+NAOB);
        assert(Guv.size(0) == 1);
        assert(G.size(0) == size_t(nel_));
        assert(G.size(1) == size_t(nmo_));
        assert(T1.size(0) == size_t(nel_));
        assert(T1.size(1) == size_t(nv));

        using ma::transposed;
        ma::product(G,rotPiu.get()({0,nmo_},{v0,vN}),
                    T1(T1.extension(0),{v0,vN}));
        // This operation might benefit from a 2-D work distribution
        ma::product(rotcPua[k].get().sliced(nu0,nu0+nu),
                    T1(T1.extension(0),{v0,vN}),
                    Guv[0]({0,nu},{v0,vN}));
        for(int v=v0; v<vN; ++v)
          if( v < nu0 || v >= nu0+nu ) {
            Guu[v] = ma::dot(rotcPua[k].get()[v],T1(T1.extension(0),v)); 
          } else
            Guu[v] = Guv[0][v-nu0][v];
      } else {
        int nel_ = NAOA+NAOB;
        assert(Guv.size(0) == 2);
        assert(G.size(0) == nel_);
        assert(G.size(1) == nmo_);
        assert(T1.size(0) == nel_);
        assert(T1.size(1) == nv);

        using ma::transposed;
        ma::product(G,rotPiu.get()({0,nmo_},{v0,vN}),
                    T1(T1.extension(0),{v0,vN}));
        // This operation might benefit from a 2-D work distribution
        // Alpha
        ma::product(rotcPua[k].get()({nu0,nu0+nu},{0,NAOA}),
                    T1({0,NAOA},{v0,vN}),
                    Guv[0]({0,nu},{v0,vN}));
        ma::product(rotcPua[k].get()({nu0,nu0+nu},{NAOA,nel_}),
                    T1({NAOA,nel_},{v0,vN}),
                    Guv[1]({0,nu},{v0,vN}));
        for(int v=v0; v<vN; ++v)
          if( v < nu0 || v >= nu0+nu ) {
            Guu[v] = ma::dot(rotcPua[k].get()[v],T1(T1.extension(0),v));
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

      static_assert(std::decay<MatA>::type::dimensionality == 2, "Wrong dimensionality");
      static_assert(std::decay<MatB>::type::dimensionality == 2, "Wrong dimensionality");
      static_assert(std::decay<MatC>::type::dimensionality == 1, "Wrong dimensionality");
      static_assert(std::decay<MatD>::type::dimensionality == 2, "Wrong dimensionality");
      int nmo_ = int(rotPiu.size(0));
      int nu = int(rotMuv.size(0));  // potentially distributed over nodes
      int nv = int(rotMuv.size(1));  // not distributed over nodes
      assert(rotPiu.size(1) == nv);
      int v0,vN;
      std::tie(v0,vN) = FairDivideBoundary(comm->rank(),nv,comm->size());
      int nu0 = rotMuv.global_offset()[0];
      ComplexType zero(0.0,0.0);

      assert(Guu.size(0) == nv);
      assert(Guv.size(0) == nu);
      assert(Guv.size(1) == nv);

      // sync first
      comm->barrier();
      int nel_ = (walker_type==CLOSED)?NAOA:(NAOA+NAOB);
      assert(G.size(0) == size_t(nel_));
      assert(G.size(1) == size_t(nmo_));
      assert(T1.size(0) == size_t(nel_));
      assert(T1.size(1) == size_t(nv));

      using ma::transposed;
      ma::product(G,rotPiu.get()({0,nmo_},{v0,vN}),
                  T1(T1.extension(0),{v0,vN}));
      // This operation might benefit from a 2-D work distribution
      ma::product(rotcPua[k].get().sliced(nu0,nu0+nu),
                  T1(T1.extension(0),{v0,vN}),
                  Guv(Guv.extension(0),{v0,vN}));
      for(int v=v0; v<vN; ++v)
        if( v < nu0 || v >= nu0+nu ) {
          Guu[v] = ma::dot(rotcPua[k].get()[v],T1(T1.extension(0),v)); 
        } else
         Guu[v] = Guv[v-nu0][v];
      comm->barrier();
    }

  protected:

    communicator* comm;

    int NMO,NAOA,NAOB;

    WALKER_TYPES walker_type;

    // bare one body hamiltonian
    CMatrix hij;

    // (potentially half rotated) one body hamiltonian
    std::vector<CVector> haj;

    /************************************************/
    // Used in the calculation of the energy
    // Coulomb matrix elements of interpolating vectors
    shmVMatrix rotMuv;

    // Orbitals at interpolating points
    shmCMatrix rotPiu;

    // Half-rotated Orbitals at interpolating points
    std::vector<shmCMatrix> rotcPua;
    /************************************************/

    /************************************************/
    // Following 3 used in calculation of vbias and vHS
    // Cholesky factorization of Muv
    shmVMatrix Luv;

    // Orbitals at interpolating points
    shmCMatrix Piu;

    // Half-rotated Orbitals at interpolating points
    std::vector<shmCMatrix> cPua;
    /************************************************/

    // one-body piece of Hamiltonian factorization
    CMatrix v0;

    ValueType E0;

    // shared memory for intermediates
    //shmSpMatrix SM_TMats;
    boost::multi::array<ComplexType,2,shared_allocator<ComplexType>> SM_TMats;

    boost::multi::array<ComplexType,2> Rbk;

    boost::multi::array<ComplexType,3> eloc;

    myTimer Timer;

    void set_shm_buffer(size_t N) {
      if(SM_TMats.num_elements() < N)
        SM_TMats.reextent({N,1});
    }

};

}

}

#endif
