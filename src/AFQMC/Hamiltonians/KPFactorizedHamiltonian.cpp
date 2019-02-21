#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<map>
#include<utility>
#include<vector>
#include<numeric>
#include <functional>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "Configuration.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Utilities/kp_utilities.hpp"
#include "AFQMC/Hamiltonians/KPFactorizedHamiltonian.h"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"
//#include "AFQMC/HamiltonianOperations/KP3IndexFactorizationIO.hpp"

namespace qmcplusplus
{

namespace afqmc
{

HamiltonianOperations KPFactorizedHamiltonian::getHamiltonianOperations(bool pureSD,
	     bool addCoulomb, WALKER_TYPES type, std::vector<PsiT_Matrix>& PsiT,
	     double cutvn, double cutv2, TaskGroup_& TGprop, TaskGroup_& TGwfn,
	     hdf_archive& hdf_restart) {

  using shmIMatrix = boost::multi::array<int,2,shared_allocator<int>>;
  using shmCVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using shmCMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using shmCTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;
  using shmSpMatrix = boost::multi::array<SPComplexType,2,shared_allocator<SPComplexType>>;
  using shmSpTensor = boost::multi::array<SPComplexType,3,shared_allocator<SPComplexType>>;
  using SpMatrix = boost::multi::array<SPComplexType,2>;
  using SpMatrix_ref = boost::multi::array_ref<SPComplexType,2>;
  using Sp3Tensor_ref = boost::multi::array_ref<SPComplexType,3>;

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if(TGwfn.Global().root()) write_hdf = !hdf_restart.closed();
  //  if(TGwfn.Global().root()) write_hdf = (hdf_restart.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  if(type==COLLINEAR)
    assert(PsiT.size()%2 == 0);
  int nspins = ((type!=COLLINEAR)?1:2);
  int ndet = PsiT.size()/nspins;

  if(ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  long nkpts, Qbeg=0, Qend, nQ;
  hdf_archive dump(TGwfn.Global());
  // right now only Node.root() reads
  if( TG.Node().root() ) {
    if(!dump.open(fileName,H5F_ACC_RDONLY)) {
      app_error()<<" Error opening integral file in THCHamiltonian. \n";
      APP_ABORT("");
    }
    if(!dump.push("Hamiltonian",false)) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Group not Hamiltonian found. \n";
      APP_ABORT("");
    }
  }

  std::vector<int> Idata(8);
  if( TG.Global().root() ) {
    if(!dump.readEntry(Idata,"dims")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading dims. \n";
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(Idata.begin(),8,0);
  nkpts = Qend = Idata[2];
  app_log()<<" nkpts: " <<nkpts <<std::endl;

  // partition Q over nodes if distributed Q

  nQ = Qend-Qbeg;

  std::vector<int> nmo_per_kp(nkpts);
  std::vector<int> nchol_per_kp(nkpts);
  std::vector<int> kminus(nkpts);
  shmIMatrix QKtok2({nkpts,nkpts},shared_allocator<int>{TG.Node()});
  ValueType E0;
  if( TG.Global().root() ) {
    if(!dump.readEntry(nmo_per_kp,"NMOPerKP")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading NMOPerKP. \n";
      APP_ABORT("");
    }
    if(!dump.readEntry(nchol_per_kp,"NCholPerKP")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading NCholPerKP. \n";
      APP_ABORT("");
    }
    if(!dump.readEntry(kminus,"MinusK")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading MinusK. \n";
      APP_ABORT("");
    }
    for(int k=0; k<nkpts; k++) {
      if(kminus[k] < k) nchol_per_kp[k] = nchol_per_kp[kminus[k]];
    }
    if(!dump.readEntry(QKtok2,"QKTok2")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading QKTok2. \n";
      APP_ABORT("");
    }
    std::vector<ValueType> E_(2);
    if(!dump.readEntry(E_,"Energies")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Energies. \n";
      APP_ABORT("");
    }
    E0 = E_[0]+E_[1];
    if(nmo_per_kp.size() != nkpts ||
       nchol_per_kp.size() != nkpts ||
       kminus.size() != nkpts ||
       QKtok2.shape()[0] != nkpts ||
       QKtok2.shape()[1] != nkpts
      ) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Inconsistent dimension (NMOPerKP,NCholPerKP,QKtTok2): "
                 <<nkpts <<" "
                 <<nmo_per_kp.size() <<" "
                 <<nchol_per_kp.size() <<" "
                 <<kminus.size() <<" "
                 <<QKtok2.shape()[0] <<" "
                 <<QKtok2.shape()[1] <<std::endl;
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(&E0,1,0);
  TG.Global().broadcast_n(nmo_per_kp.begin(),nmo_per_kp.size(),0);
  TG.Global().broadcast_n(nchol_per_kp.begin(),nchol_per_kp.size(),0);
  TG.Global().broadcast_n(kminus.begin(),kminus.size(),0);
  if(TG.Node().root())
    TG.Cores().broadcast_n(std::addressof(*QKtok2.origin()),QKtok2.num_elements(),0);
  TG.Node().barrier();

  int Q0=-1;  // stores the index of the Q=(0,0,0) Q-point
              // this must always exist
  for(int Q=0; Q<nkpts; Q++) {
    if(kminus[Q]==Q) {
      bool found=true;
      for(int KI=0; KI<nkpts; KI++)
        if(KI != QKtok2[Q][KI]) {
          found = false;
          break;
        }
      if(found) {
        if( Q0 > 0 )
          APP_ABORT(" Error: @ Q-points satisfy Q=0 condition.\n");
        Q0=Q;
      } else {
        if( nkpts%2 != 0)
          APP_ABORT(" Error: Unexpected situation: Q==(-Q)!=Q0 and odd Nk. \n");
      }
    }
  }
  if(Q0<0)
    APP_ABORT(" Error: Can not find Q=0 Q-point.\n");
  int nmo_max = *std::max_element(nmo_per_kp.begin(),nmo_per_kp.end());
  int nchol_max = *std::max_element(nchol_per_kp.begin(),nchol_per_kp.end());
  shmCTensor H1({nkpts,nmo_max,nmo_max},shared_allocator<ComplexType>{TG.Node()});
  std::vector<shmSpMatrix> LQKikn;
  LQKikn.reserve(nkpts);
  for(int Q=0; Q<nkpts; Q++)
    if( Q==Q0 )
      LQKikn.emplace_back( shmSpMatrix({nkpts,nmo_max*nmo_max*nchol_per_kp[Q]},
                                   shared_allocator<SPComplexType>{TG.Node()}) );
    else if( kminus[Q] == Q )   // only storing half of K points and using symmetry
      LQKikn.emplace_back( shmSpMatrix({nkpts/2,nmo_max*nmo_max*nchol_per_kp[Q]},
                                   shared_allocator<SPComplexType>{TG.Node()}) );
    else if(Q < kminus[Q])
      LQKikn.emplace_back( shmSpMatrix({nkpts,nmo_max*nmo_max*nchol_per_kp[Q]},
                                   shared_allocator<SPComplexType>{TG.Node()}) );
    else // Q > kminus[Q]
      LQKikn.emplace_back( shmSpMatrix({1,1},
                                   shared_allocator<SPComplexType>{TG.Node()}) );

  if( TG.Node().root() ) {
    // now read H1_kpQ
    for(int Q=0; Q<nkpts; Q++) {
      // until double_hyperslabs work!
      boost::multi::array<ComplexType,2> h1({nmo_per_kp[Q],nmo_per_kp[Q]});
      if(!dump.readEntry(h1,std::string("H1_kp")+std::to_string(Q))) {
        app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading /Hamiltonian/H1_kp" <<Q <<". \n";
        APP_ABORT("");
      }
      H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}) = h1;
      //H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}) = h1;
    }
    // read LQ
    if(!dump.push("KPFactorized",false)) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Group KPFactorized not found. \n";
      APP_ABORT("");
    }
    for(int Q=0; Q<nkpts; Q++) {
      using std::conj;
      if(Q==Q0) {
        if(!dump.readEntry(LQKikn[Q],std::string("L")+std::to_string(Q))) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading /Hamiltonian/KPFactorized/L" <<Q <<". \n";
          APP_ABORT("");
        }
      } else if(kminus[Q]==Q) {
        SpMatrix L_({nkpts,nmo_max*nmo_max*nchol_per_kp[Q]});
        auto&& LQ(LQKikn[Q]);
        if(!dump.readEntry(L_,std::string("L")+std::to_string(Q))) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading /Hamiltonian/KPFactorized/L" <<Q <<". \n";
          APP_ABORT("");
        }
        int kpos=0;
        for(int K=0; K<nkpts; K++) {
          assert(K != QKtok2[Q][K]);
          if(K < QKtok2[Q][K]) {
            LQ[kpos] = L_[K];
            kpos++;
          }
        }
      } else if(Q < kminus[Q]) {
        if(!dump.readEntry(LQKikn[Q],std::string("L")+std::to_string(Q))) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading /Hamiltonian/KPFactorized/L" <<Q <<". \n";
          APP_ABORT("");
        }
        if(LQKikn[Q].shape()[0] != nkpts || LQKikn[Q].shape()[1] != nmo_max*nmo_max*nchol_per_kp[Q]) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading /Hamiltonian/KPFactorized/L" <<Q <<". \n"
                 <<" Unexpected dimensins: " <<LQKikn[Q].shape()[0] <<" " <<LQKikn[Q].shape()[1] <<std::endl;
          APP_ABORT("");
        }
      }
    }
    dump.pop();
  }
  TG.Node().barrier();

  // calculate vn0
  shmCTensor vn0({nkpts,nmo_max,nmo_max},shared_allocator<ComplexType>{TG.Node()});

  // generate nocc_per_kp using PsiT and nmo_per_kp
  shmIMatrix nocc_per_kp({ndet,nspins*nkpts},shared_allocator<int>{TG.Node()});
  if(TG.Node().root()) {
    if(type==COLLINEAR) {
      for(int i=0; i<ndet; i++) {
        if(not get_nocc_per_kp(nmo_per_kp,PsiT[2*i],nocc_per_kp[i]({0,nkpts}))) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
        if(not get_nocc_per_kp(nmo_per_kp,PsiT[2*i+1],nocc_per_kp[i]({nkpts,2*nkpts}))) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
      }
    } else {
      for(int i=0; i<ndet; i++)
        if(not get_nocc_per_kp(nmo_per_kp,PsiT[i],nocc_per_kp[i])) {
          app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
    }
  }
  TG.Node().barrier();
  int nocc_max = *std::max_element(std::addressof(*nocc_per_kp.origin()),
                                   std::addressof(*nocc_per_kp.origin())+nocc_per_kp.num_elements());

  /* half-rotate LQ and H1:
   * Given that PsiT = H(SM),
   * h[K][a][k] = sum_i PsiT[K][a][i] * h[K][i][k]
   * L[Q][K][a][k][n] = sum_i PsiT[K][a][i] * L[Q][K][i][k][n]
   * L[Q][K][l][b][n] = sum_i PsiT[K][b][] * L[Q][K][l][k][n]*
   * LQKak has a special transposition to facilitate computations
   * of the energy, and they are stored with padding to max linear dimension
   * LQKank[Q][K][...] = LQKank[Q][K][a][n][k] = LQKakn[Q][K][a][k][n]
   */
  std::vector<shmSpMatrix> LQKank;
  LQKank.reserve(ndet*nspins*(nkpts+1));  // storing 2 components for Q=0, since it is not assumed symmetric
  shmCMatrix haj({ndet*nkpts,(type==COLLINEAR?2:1)*nocc_max*nmo_max},shared_allocator<ComplexType>{TG.Node()});
  if(TG.Node().root()) std::fill_n(haj.origin(),haj.num_elements(),ComplexType(0.0));
  int ank_max = nocc_max*nchol_max*nmo_max;
  for(int nd=0; nd<ndet; nd++) {
    for(int Q=0; Q<(nkpts+1); Q++) {
      LQKank.emplace_back(shmSpMatrix({nkpts,ank_max},shared_allocator<ComplexType>{TG.Node()}));
    }
    if(type==COLLINEAR) {
      for(int Q=0; Q<(nkpts+1); Q++) {
        LQKank.emplace_back(shmSpMatrix({nkpts,ank_max},shared_allocator<ComplexType>{TG.Node()}));
      }
    }
  }

  for(int nd=0, nt=0, nq0=0; nd<ndet; nd++, nq0+=(nkpts+1)*nspins) {
    for(int Q=0; Q<(nkpts+1); Q++) {
      for(int K=0; K<nkpts; K++, nt++) {
        if(nt%TG.Node().size() == TG.Node().rank()) {
          std::fill_n(std::addressof(*LQKank[nq0+Q][K].origin()),LQKank[nq0+Q][K].num_elements(),SPComplexType(0.0));
          if(type==COLLINEAR) {
            std::fill_n(std::addressof(*LQKank[nq0+nkpts+1+Q][K].origin()),LQKank[nq0+nkpts+1+Q][K].num_elements(),SPComplexType(0.0));
          }
        }
      }
    }
  }
  TG.Node().barrier();
  boost::multi::array<SPComplexType,2> buff({nmo_max,nchol_max});
  for(int nd=0, nt=0, nq0=0; nd<ndet; nd++, nq0+=(nkpts+1)*nspins) {
    for(int Q=0; Q<nkpts; Q++) {
      for(int K=0; K<nkpts; K++, nt++) {
        if(nt%TG.Global().size() == TG.Global().rank()) {
          // haj and add half-transformed right-handed rotation for Q=0
          int Qm = kminus[Q];
          int QK = QKtok2[Q][K];
          int na = nocc_per_kp[nd][K];
          int nb = nocc_per_kp[nd][nkpts+K];
          int ni = nmo_per_kp[K];
          int nk = nmo_per_kp[QK];
          int nchol = nchol_per_kp[Q];
          if(Q==0) {
            if(type==COLLINEAR) {
              { // Alpha
                auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd],K);
                assert(Psi.shape()[0] == na);
                boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+K].origin()),
                                                             {na,ni});
                ma::product(Psi,H1[K]({0,ni},{0,ni}),haj_r);
              }
              { // Beta
                auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd+1],K);
                assert(Psi.shape()[0] == nb);
                boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+K].origin())+
                                                                            na*ni,
                                                             {nb,ni});
                ma::product(Psi,H1[K]({0,ni},{0,ni}),haj_r);
              }
            } else {
              auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[nd],K);
              assert(Psi.shape()[0] == na);
              boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+K].origin()),
                                                           {na,ni});
              ma::product(ComplexType(2.0),Psi,H1[K]({0,ni},{0,ni}),
                          ComplexType(0.0),haj_r);
            }
          }
          if(type==COLLINEAR) {
            { // Alpha
// change get_PsiK to cast to the value_type of the result
              auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[2*nd],K);
              assert(Psi.shape()[0] == nocc_per_kp[nd][K]);
              if(Q < Qm || Q==Q0 || ((Q==Qm)&&(K<QK))) {
                int kpos = K;
                if( Q==Qm && Q!=Q0 ) { //find position in symmetric list
                  kpos=0;
                  for(int K_=0; K_<K; K_++)
                    if(K_ < QKtok2[Q][K_]) kpos++;
                }
                Sp3Tensor_ref Likn(std::addressof(*LQKikn[Q][kpos].origin()),
                                              {ni,nk,nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+Q][K].origin()),
                                              {na,nchol,nk});
                ma_rotate::getLank(Psi,Likn,Lank,buff);
                if(Q==Q0) {
                  assert(K==QK);
                  Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+nkpts][K].origin()),
                                              {na,nchol,nk});
                  ma_rotate::getLank_from_Lkin(Psi,Likn,Lank,buff);
                }
              } else {
                int kpos = QK;
                if( Q==Qm ) { //find position in symmetric list
                  kpos=0;
                  for(int K_=0; K_<QK; K_++)
                    if(K_ < QKtok2[Q][K_]) kpos++;
                }
                Sp3Tensor_ref Lkin(std::addressof(*LQKikn[Qm][QK].origin()),
                                              {nk,ni,nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+Q][K].origin()),
                                              {na,nchol,nk});
                ma_rotate::getLank_from_Lkin(Psi,Lkin,Lank,buff);
              }
            }
            { // Beta
// change get_PsiK to cast to the value_type of the result
              auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[2*nd+1],K);
              assert(Psi.shape()[0] == nb);
              if(Q < Qm || Q==Q0 || ((Q==Qm)&&(K<QK))) {
                int kpos = K;
                if( Q==Qm && Q!=Q0 ) { //find position in symmetric list
                  kpos=0;
                  for(int K_=0; K_<K; K_++)
                    if(K_ < QKtok2[Q][K_]) kpos++;
                }
                Sp3Tensor_ref Likn(std::addressof(*LQKikn[Q][kpos].origin()),
                                                {ni,nk,nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+nkpts+1+Q][K].origin()),
                                                {nb,nchol,nk});
                ma_rotate::getLank(Psi,Likn,Lank,buff);
                if(Q==Q0) {
                  assert(K==QK);
                  Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+nkpts+1+nkpts][K].origin()),
                                                {nb,nchol,nk});
                  ma_rotate::getLank_from_Lkin(Psi,Likn,Lank,buff);
                }
              } else {
                int kpos = QK;
                if( Q==Qm ) { //find position in symmetric list
                  kpos=0;
                  for(int K_=0; K_<QK; K_++)
                    if(K_ < QKtok2[Q][K_]) kpos++;
                }
                Sp3Tensor_ref Lkin(std::addressof(*LQKikn[Qm][QK].origin()),
                                                {nk,ni,nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+nkpts+1+Q][K].origin()),
                                                {nb,nchol,nk});
                ma_rotate::getLank_from_Lkin(Psi,Lkin,Lank,buff);
              }
            }
          } else {
// change get_PsiK to cast to the value_type of the result
            auto Psi = get_PsiK<SpMatrix>(nmo_per_kp,PsiT[nd],K);
            assert(Psi.shape()[0] == na);
            if(Q < Qm || Q==Q0 || ((Q==Qm)&&(K<QK))) {
              int kpos = K;
              if( Q==Qm && Q!=Q0 ) { //find position in symmetric list
                kpos=0;
                for(int K_=0; K_<K; K_++)
                  if(K_ < QKtok2[Q][K_]) kpos++;
              }
              Sp3Tensor_ref Likn(std::addressof(*LQKikn[Q][kpos].origin()),
                                              {ni,nk,nchol});
              Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+Q][K].origin()),
                                              {na,nchol,nk});
              ma_rotate::getLank(Psi,Likn,Lank,buff);
              if(Q==Q0) {
                assert(K==QK);
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+nkpts][K].origin()),
                                              {na,nchol,nk});
                ma_rotate::getLank_from_Lkin(Psi,Likn,Lank,buff);
              }
            } else {
              int kpos = QK;
              if( Q==Qm ) { //find position in symmetric list
                kpos=0;
                for(int K_=0; K_<QK; K_++)
                  if(K_ < QKtok2[Q][K_]) kpos++;
              }
              Sp3Tensor_ref Lkin(std::addressof(*LQKikn[Qm][kpos].origin()),
                                              {nk,ni,nchol});
              Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0+Q][K].origin()),
                                              {na,nchol,nk});
              ma_rotate::getLank_from_Lkin(Psi,Lkin,Lank,buff);
            }
          }
        }
      }
    }
  }
  TG.Global().barrier();
  if(TG.Node().root()) {
    TG.Cores().all_reduce_in_place_n(std::addressof(*haj.origin()),
                                     haj.num_elements(),std::plus<>());
    for(int Q=0; Q<LQKank.size(); Q++)
      TG.Cores().all_reduce_in_place_n(std::addressof(*LQKank[Q].origin()),
                                       LQKank[Q].num_elements(),std::plus<>());
    std::fill_n(std::addressof(*vn0.origin()),vn0.num_elements(),ComplexType(0.0));
  }
  TG.Node().barrier();

  // calculate (only Q=0) vn0(I,L) = -0.5 sum_K sum_j sum_n L[0][K][i][j][n] conj(L[0][K][l][j][n])
  for(int K=0; K<nkpts; K++) {
    if(K%TG.Node().size() == TG.Node().rank()) {
      boost::multi::array_ref<SPComplexType,2> Likn(std::addressof(*LQKikn[Q0][K].origin()),
                                                   {nmo_per_kp[K],nmo_per_kp[K]*nchol_per_kp[0]});
      using ma::H;
      ma::product(-0.5,Likn,H(Likn),0.0,vn0[K]({0,nmo_per_kp[K]},{0,nmo_per_kp[K]}));
    }
  }
  TG.Node().barrier();
  // in parallel, whoever has Q=0 calculates and bcasts

//  TG.Node().barrier();
//  if(TG.Node().root())
//    TG.Cores().all_reduce_in_place_n(std::addressof(*vn0.origin()),vn0.num_elements(),std::plus<>());

  if( TG.Node().root() ) {
    dump.pop();
    dump.close();
  }

  int global_ncvecs = std::accumulate(nchol_per_kp.begin(),nchol_per_kp.end(),0);
/*
  ComplexType E_(0.0);
  boost::multi::array<ComplexType,3> G({nkpts,nmo_per_kp[0],nmo_per_kp[0]});
  for(int K=0; K<nkpts; K++) {
    auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[0],K);
    ma::product(ma::H(Psi),Psi,G[K]);
    ma::transpose(G[K]);
  }
  boost::multi::array<ComplexType,3> Gc({nkpts,nocc_per_kp[0][0],nmo_per_kp[0]});
  for(int K=0; K<nkpts; K++) {
    auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[0],K);
    for(int a=0; a<nocc_per_kp[0][K]; a++)
      for(int j=0; j<nmo_per_kp[K]; j++)
        Gc[K][a][j] = std::conj(Psi[a][j]);
  }
  for(int K=0; K<nkpts; K++) {
    auto&& G_=G[K];
    for(int i=0; i<nmo_per_kp[K]; i++)
      for(int j=0; j<nmo_per_kp[K]; j++)
        E_ += H1[K][i][j]*G_[i][j];
  }
  std::cout<<" E1+E0: " <<std::setprecision(12) <<E0+2*E_ <<"  ";
  E_=0.0;
  for(int K=0; K<nkpts; K++) {
    auto&& G_=Gc[K];
    for(int a=0, aj=0; a<nocc_per_kp[0][K]; a++)
      for(int j=0; j<nmo_per_kp[K]; j++, ++aj)
        E_ += haj[K][aj]*G_[a][j];
  }
  std::cout<<E0+E_ <<std::endl;

  boost::multi::array<int,2> KK2Q({nkpts,nkpts});
  for(int KI=0; KI<nkpts; KI++)
  for(int KK=0; KK<nkpts; KK++) {
    KK2Q[KI][KK]=-1;
    for(int Q=0; Q<nkpts; Q++)
      if(QKtok2[Q][KI]==KK) {
        KK2Q[KI][KK]=Q;
        break;
      }
    assert(KK2Q[KI][KK]>=0);
  }

  std::ofstream out("fact_ints.dat");
  boost::multi::array<ComplexType,2> IJKL({nmo_max*nmo_max,nmo_max*nmo_max});
  boost::multi::array_ref<ComplexType,4> IJKL4D(IJKL.origin(),{nmo_max,nmo_max,nmo_max,nmo_max});
  ComplexType EX(0.0);
  ComplexType EJ(0.0);
  for(int KI=0; KI<nkpts; KI++)
  for(int KL=0; KL<nkpts; KL++) {
    for(int KK=0; KK<nkpts; KK++)
    {
      int Q = KK2Q[KI][KK];
      int KJ = QKtok2[Q][KL];
      boost::multi::array_ref<ComplexType,2> LKI(std::addressof(*LQKikn[Q][KI].origin()),{nmo_max*nmo_max,nchol_per_kp[Q]});
      boost::multi::array_ref<ComplexType,2> LKL(std::addressof(*LQKikn[Q][KL].origin()),{nmo_max*nmo_max,nchol_per_kp[Q]});
      ma::product(LKI,ma::H(LKL),IJKL);
      for(int i=0; i<nmo_per_kp[0]; i++)
      for(int k=0; k<nmo_per_kp[0]; k++)
      for(int j=0; j<nmo_per_kp[0]; j++)
      for(int l=0; l<nmo_per_kp[0]; l++)
        out<<Q <<" " <<KI <<" " <<KK <<" " <<KJ <<" " <<KL <<" " <<i <<" " <<k <<" " <<j <<" " <<l <<" " <<IJKL4D[i][k][l][j] <<"\n";
      if(KI==KK && KL == KJ) { // EJ
        for(int i=0; i<nmo_per_kp[0]; i++)
        for(int k=0; k<nmo_per_kp[0]; k++)
        for(int l=0; l<nmo_per_kp[0]; l++)
        for(int j=0; j<nmo_per_kp[0]; j++) {
          EJ += 0.5*IJKL4D[i][k][l][j] * G[KI][i][k] * G[KJ][j][l];
        }
      }
      if(KI==KL && KJ == KK) { // EX
        for(int i=0; i<nmo_per_kp[0]; i++)
        for(int k=0; k<nmo_per_kp[0]; k++)
        for(int l=0; l<nmo_per_kp[0]; l++)
        for(int j=0; j<nmo_per_kp[0]; j++) {
          EX += 0.5*IJKL4D[i][k][l][j] * G[KI][i][l] * G[KJ][j][k];
        }
      }
    }
  }
  out.close();
  std::cout<<" EX: " <<std::setprecision(12) <<EX*2 <<std::endl;
  std::cout<<" EJ: " <<std::setprecision(12) <<EJ*4 <<std::endl;
// */
/*
  boost::multi::array<ComplexType,2> ABKL({nocc_max*nmo_max,nocc_max*nmo_max});
  boost::multi::array_ref<ComplexType,4> ABKL4D(ABKL.origin(),{nocc_max,nmo_max,nmo_max,nocc_max});
  EX = EJ = ComplexType(0.0);
  for(int KI=0; KI<nkpts; KI++)
  for(int KL=0; KL<nkpts; KL++) {
    for(int KK=0; KK<nkpts; KK++)
    {
      int Q = KK2Q[KI][KK];
      int Qm = kminus[Q];
      assert(QKtok2[Q][KI]==KK);
      int KJ = QKtok2[Q][KL];
      boost::multi::array_ref<ComplexType,3> LKI(std::addressof(*LQKank[Q][KI].origin()),{nocc_max,nchol_per_kp[Q],nmo_max});
      //boost::multi::array_ref<ComplexType,3> LKL(std::addressof(*LQKlnb[Q][KL].origin()),{nmo_max,nchol_per_kp[Q],nocc_max});
      boost::multi::array_ref<ComplexType,3> LKL(std::addressof(*LQKank[Qm][KJ].origin()),{nmo_max,nchol_per_kp[Q],nocc_max});
      if((KI==KK && KL == KJ) || (KI==KL && KJ == KK)) {
        for(int a=0; a<nocc_per_kp[0][0]; a++)
        for(int k=0; k<nmo_per_kp[0]; k++)
        for(int l=0; l<nmo_per_kp[0]; l++)
        for(int b=0; b<nocc_per_kp[0][0]; b++) {
          ABKL4D[a][k][l][b] = ComplexType(0.0);
          for(int n=0; n<nchol_per_kp[Q]; n++)
            ABKL4D[a][k][l][b] += LKI[a][n][k] * conj(LKL[b][n][l]);
        }
      }
      if(KI==KK && KL == KJ) { // EJ
        for(int a=0; a<nocc_per_kp[0][0]; a++)
        for(int k=0; k<nmo_per_kp[0]; k++)
        for(int l=0; l<nmo_per_kp[0]; l++)
        for(int b=0; b<nocc_per_kp[0][0]; b++) {
          EJ += 0.5*ABKL4D[a][k][l][b] * Gc[KI][a][k] * Gc[KJ][b][l];
        }
      }
      if(KI==KL && KJ == KK) { // EX
        for(int a=0; a<nocc_per_kp[0][0]; a++)
        for(int k=0; k<nmo_per_kp[0]; k++)
        for(int l=0; l<nmo_per_kp[0]; l++)
        for(int b=0; b<nocc_per_kp[0][0]; b++) {
          EX += 0.5*ABKL4D[a][k][l][b] * Gc[KI][a][l] * Gc[KJ][b][k];
        }
      }
    }
  }
  std::cout<<" EX: " <<std::setprecision(12) <<EX*2 <<std::endl;
  std::cout<<" EJ: " <<std::setprecision(12) <<EJ*4 <<std::endl;
// */

  std::vector<RealType> gQ(nkpts);
  if(nsampleQ>0) {
    app_log()<<" Sampling EXX energy using distribution over Q vector obtained from "
             <<" trial energy. \n";

    RealType scl = (type==CLOSED?2.0:1.0);
    size_t nqk=0;
    for(int Q=0; Q<nkpts; ++Q) {            // momentum conservation index
      int Qm = kminus[Q];
      int Qm_ = (Q==Q0?nkpts:Qm);
      for(int Ka=0; Ka<nkpts; ++Ka) {
        int Kk = QKtok2[Q][Ka];
        int Kb = Kk;
        int Kl = QKtok2[Qm][Kb];
        if( (Ka != Kl) || (Kb != Kk) )
          APP_ABORT(" Error: Problems with EXX.\n");
        if((nqk++)%TG.Global().size() == TG.Global().rank()) {
          int nchol = nchol_per_kp[Q];
          int nl = nmo_per_kp[Kl];
          int nb = nocc_per_kp[0][Kb];
          int nk = nmo_per_kp[Kk];
          int na = nocc_per_kp[0][Ka];

          SpMatrix_ref Lank(std::addressof(*LQKank[Q][Ka].origin()),
                                           {na*nchol,nk});
          SpMatrix_ref Lbnl(std::addressof(*LQKank[Qm_][Kb].origin()),
                                           {nb*nchol,nl});
          SpMatrix Tban({nb,na*nchol});
          Sp3Tensor_ref T3ban(Tban.origin(),{nb,na,nchol});
          SpMatrix Tabn({na,nb*nchol});
          Sp3Tensor_ref T3abn(Tabn.origin(),{na,nb,nchol});

          auto Gal = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[0],Ka);
          auto Gbk = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[0],Kb);
          for(int a=0; a<na; ++a)
            for(int l=0; l<nl; ++l)
              Gal[a][l] = conj(Gal[a][l]);
          for(int b=0; b<nb; ++b)
            for(int k=0; k<nk; ++k)
              Gbk[b][k] = conj(Gbk[b][k]);

          ma::product(Gal,ma::T(Lbnl),Tabn);
          ma::product(Gbk,ma::T(Lank),Tban);

          ComplexType E_(0.0);
          for(int a=0; a<na; ++a)
            for(int b=0; b<nb; ++b)
              E_ += ma::dot(T3abn[a][b],T3ban[b][a]);
          gQ[Q] -= scl*0.5*real(E_);
        }
        if(type==COLLINEAR) {
          APP_ABORT(" Finish UHF.\n ");
        }
      }
    }
    TG.Global().all_reduce_in_place_n(gQ.begin(),nkpts,std::plus<>());
    RealType E_ = std::accumulate(gQ.begin(),gQ.end(),RealType(0.0));
    for( auto& v: gQ ) v /= E_;
    app_log()<<" EXX: " <<E_ <<std::endl;
    for( auto v: gQ ) {
      if(v < 0.0)
        APP_ABORT(" Error: g(Q) < 0.0, implement shift to g(Q). \n")
    }
  }

  return HamiltonianOperations(KP3IndexFactorization(TGwfn.TG_local(), type,std::move(nmo_per_kp),
            std::move(nchol_per_kp),std::move(kminus),std::move(nocc_per_kp),
            std::move(QKtok2),std::move(H1),std::move(haj),std::move(LQKikn),
            std::move(LQKank),std::move(vn0),std::move(gQ),nsampleQ,E0,global_ncvecs));

}

}
}
