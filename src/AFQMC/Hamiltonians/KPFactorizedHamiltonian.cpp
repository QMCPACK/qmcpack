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

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if(TGwfn.Global().root()) write_hdf = (hdf_restart.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  if(type==COLLINEAR)
    assert(PsiT.size()%2 == 0);
  int nspins = ((type!=COLLINEAR)?1:2);
  int ndet = PsiT.size()/nspins;

  if(ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  long nkpts, Q0=0, QN, nQ; 
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
    if(!dump.read(Idata,"dims")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading dims. \n";
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(Idata.begin(),8,0);
  nkpts = QN = Idata[2];
  app_log()<<" nkpts: " <<nkpts <<std::endl;
    
  // partition Q over nodes if distributed Q
   
  nQ = QN-Q0; 

  std::vector<int> nmo_per_kp(nkpts);
  std::vector<int> nchol_per_kp(nkpts);
  std::vector<int> kminus(nkpts);
  shmIMatrix QKtok2({nkpts,nkpts},shared_allocator<int>{TG.Node()});
  ValueType E0; 
  if( TG.Global().root() ) {
    if(!dump.read(nmo_per_kp,"NMOPerKP")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading NMOPerKP. \n";
      APP_ABORT("");
    }
    if(!dump.read(nchol_per_kp,"NCholPerKP")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading NCholPerKP. \n";
      APP_ABORT("");
    }
    if(!dump.read(kminus,"MinusK")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading MinusK. \n";
      APP_ABORT("");
    }
    if(!dump.read(QKtok2,"QKTok2")) {
      app_error()<<" Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading QKTok2. \n";
      APP_ABORT("");
    }
    std::vector<ValueType> E_(2);
    if(!dump.read(E_,"Energies")) {
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

  int nmo_max = *std::max_element(nmo_per_kp.begin(),nmo_per_kp.end());
  int nchol_max = *std::max_element(nchol_per_kp.begin(),nchol_per_kp.end());
  shmCTensor H1({nkpts,nmo_max,nmo_max},shared_allocator<ComplexType>{TG.Node()}); 
  std::vector<shmSpMatrix> LQKikn;
  LQKikn.reserve(nkpts);   
  for(int Q=0; Q<nkpts; Q++) 
    LQKikn.emplace_back( shmSpMatrix({nkpts,nmo_max*nmo_max*nchol_per_kp[Q]},
                                   shared_allocator<SPComplexType>{TG.Node()}) );

  if( TG.Node().root() ) {
    // now read H1_kpQ 
    for(int Q=0; Q<nkpts; Q++) {
      // until double_hyperslabs work!
      boost::multi::array<ComplexType,2> h1({nmo_per_kp[Q],nmo_per_kp[Q]}); 
      if(!dump.read(h1,std::string("H1_kp")+std::to_string(Q))) {
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
      if(!dump.read(LQKikn[Q],std::string("L")+std::to_string(Q))) {
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
   * LQKak and LQKlb have a special transposition to facilitate computations 
   * of the energy, and they are stored with padding to max linear dimension
   * LQKank[Q][K][...] = LQKank[Q][K][a][n][k] = LQKakn[Q][K][a][k][n] 
   * LQKlnb[Q][K][...] = LQKlnb[Q][K][l][n][b] = LQKakn[Q][K][l][b][n] 
   */
  std::vector<shmSpMatrix> LQKank;
  LQKank.reserve(ndet*nspins*nkpts);   
  std::vector<shmSpMatrix> LQKlnb;
  LQKlnb.reserve(ndet*nspins*nkpts);   
/*
  std::vector<shmCVector> haj;
  haj.reserve(ndet*nkpts);
  // allocation only
  for(int nd=0; nd<ndet; nd++) {
    for(int Q=0; Q<nkpts; Q++) {
      int nak = nocc_per_kp[nd][Q]*nmo_per_kp[Q];
      if(type==COLLINEAR) nak += nocc_per_kp[nd][Q+nkpts]*nmo_per_kp[Q];
      haj.emplace_back(shmCVector({boost::multi::index_extension{nak}},shared_allocator<ComplexType>{TG.Node()}));
    }
  }
*/
  shmCMatrix haj({ndet*nkpts,(type==COLLINEAR?2:1)*nocc_max*nmo_max},shared_allocator<ComplexType>{TG.Node()});
  if(TG.Node().root()) std::fill_n(haj.origin(),haj.num_elements(),ComplexType(0.0));
  int ank_max = nocc_max*nchol_max*nmo_max;
  for(int nd=0; nd<ndet; nd++) {
    for(int Q=0; Q<nkpts; Q++) {
      LQKank.emplace_back(shmSpMatrix({nkpts,ank_max},shared_allocator<ComplexType>{TG.Node()}));
      LQKlnb.emplace_back(shmSpMatrix({nkpts,ank_max},shared_allocator<ComplexType>{TG.Node()}));
    }
    if(type==COLLINEAR) {
      for(int Q=0; Q<nkpts; Q++) {
        LQKank.emplace_back(shmSpMatrix({nkpts,ank_max},shared_allocator<ComplexType>{TG.Node()}));
        LQKlnb.emplace_back(shmSpMatrix({nkpts,ank_max},shared_allocator<ComplexType>{TG.Node()}));
      }
    }
  }

  for(int nd=0, nt=0; nd<ndet; nd++) {
    for(int Q=0; Q<nkpts; Q++) {
      for(int K=0; K<nkpts; K++, nt++) {
        if(nt%TG.Node().size() == TG.Node().rank()) {
          std::fill_n(std::addressof(*LQKank[Q][K].origin()),LQKank[Q][K].num_elements(),SPComplexType(0.0));
          std::fill_n(std::addressof(*LQKlnb[Q][K].origin()),LQKlnb[Q][K].num_elements(),SPComplexType(0.0));
          if(type==COLLINEAR) {
            std::fill_n(std::addressof(*LQKank[nkpts+Q][K].origin()),LQKank[nkpts+Q][K].num_elements(),SPComplexType(0.0));
            std::fill_n(std::addressof(*LQKlnb[nkpts+Q][K].origin()),LQKlnb[nkpts+Q][K].num_elements(),SPComplexType(0.0));
          }
        }
      }
    }
  }
  TG.Node().barrier();
  boost::multi::array<SPComplexType,2> buff({nmo_max,nchol_max});  
  for(int nd=0, nt=0; nd<ndet; nd++) {
    for(int Q=0; Q<nkpts; Q++) {
      for(int K=0; K<nkpts; K++, nt++) {
        if(nt%TG.Global().size() == TG.Global().rank()) {
          // haj  
          if(Q==0) {
            if(type==COLLINEAR) {
              { // Alpha 
                auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd],K); 
                assert(Psi.shape()[0] == nocc_per_kp[nd][K]);
                boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+K].origin()),
                                                             {nocc_per_kp[nd][K],nmo_per_kp[K]});
                ma::product(Psi,H1[K]({0,nmo_per_kp[K]},{0,nmo_per_kp[K]}),haj_r);
              }
              { // Beta 
                auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd+1],K);
                assert(Psi.shape()[0] == nocc_per_kp[nd][nkpts+K]);
                boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+K].origin())+
                                                                            nocc_per_kp[nd][K]*nmo_per_kp[K],
                                                             {nocc_per_kp[nd][nkpts+K],nmo_per_kp[K]});
                ma::product(Psi,H1[K]({0,nmo_per_kp[K]},{0,nmo_per_kp[K]}),haj_r);
              }
            } else {
              auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[nd],K);
              assert(Psi.shape()[0] == nocc_per_kp[nd][K]);
              boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+K].origin()),
                                                           {nocc_per_kp[nd][K],nmo_per_kp[K]});
              ma::product(ComplexType(2.0),Psi,H1[K]({0,nmo_per_kp[K]},{0,nmo_per_kp[K]}),
                          ComplexType(0.0),haj_r);
            }
          }  
          if(type==COLLINEAR) {
            { // Alpha 
// change get_PsiK to cast to the value_type of the result 
              auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[2*nd],K);
              assert(Psi.shape()[0] == nocc_per_kp[nd][K]);
              boost::multi::array_ref<SPComplexType,2> Likn(std::addressof(*LQKikn[Q][K].origin()),
                                                           {nmo_per_kp[K],nmo_per_kp[QKtok2[Q][K]]*nchol_per_kp[Q]});
              boost::multi::array_ref<SPComplexType,2> Lakn(std::addressof(*LQKank[Q][K].origin()),
                                                           {nocc_per_kp[nd][K],nmo_per_kp[QKtok2[Q][K]]*nchol_per_kp[Q]});
              ma::product(Psi,Likn,Lakn);
              // transpose to form expected by KP3IndexFactorization  
              for(int a=0; a<nocc_per_kp[nd][K]; a++) {
                boost::multi::array_ref<SPComplexType,2> Lkn(Lakn[a].origin(),
                                                           {nmo_per_kp[QKtok2[Q][K]],nchol_per_kp[Q]});
                boost::multi::array_ref<SPComplexType,2> Lnk(Lakn[a].origin(),
                                                           {nchol_per_kp[Q],nmo_per_kp[QKtok2[Q][K]]});
                buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}) = Lkn;
                ma::transpose(buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}),Lnk);      
              }  
              // doing this "by-hand" now
              Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[nd],QKtok2[Q][K]);
              boost::multi::array_ref<SPComplexType,3> Llbn(std::addressof(*LQKlnb[Q][K].origin()),
                                                           {nmo_per_kp[K],nocc_per_kp[nd][QKtok2[Q][K]],nchol_per_kp[Q]});
              for(int l=0; l<nmo_per_kp[K]; ++l) {
                auto psi_bj = Psi.origin();
                for(int b=0; b<nocc_per_kp[nd][QKtok2[Q][K]]; ++b) {
                  auto Likn_jn = Likn[l].origin();
                  for(int j=0, jn=0; j<nmo_per_kp[QKtok2[Q][K]]; ++j, ++psi_bj) {
                    auto Llbn_lbn = Llbn[l][b].origin();
                    for(int n=0; n<nchol_per_kp[Q]; ++n, ++Likn_jn, ++Llbn_lbn) 
                      (*Llbn_lbn) += (*psi_bj) * conj(*Likn_jn);
                  }  
                }  
              }  
              for(int l=0; l<nmo_per_kp[K]; ++l) {
                boost::multi::array_ref<SPComplexType,2> Lbn(Llbn[l].origin(),
                                                           {nocc_per_kp[nd][QKtok2[Q][K]],nchol_per_kp[Q]});
                boost::multi::array_ref<SPComplexType,2> Lnb(Llbn[l].origin(),
                                                           {nchol_per_kp[Q],nocc_per_kp[nd][QKtok2[Q][K]]});
                buff({0,Lbn.shape()[0]},{0,Lbn.shape()[1]}) = Lbn;
                ma::transpose(buff({0,Lbn.shape()[0]},{0,Lbn.shape()[1]}),Lnb);
              }
            }
            { // Beta 
// change get_PsiK to cast to the value_type of the result 
              auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[2*nd+1],K);
              assert(Psi.shape()[0] == nocc_per_kp[nd][nkpts+K]);
              boost::multi::array_ref<SPComplexType,2> Likn(std::addressof(*LQKikn[Q][K].origin()),
                                                           {nmo_per_kp[K],nmo_per_kp[QKtok2[Q][K]]*nchol_per_kp[Q]});
              boost::multi::array_ref<SPComplexType,2> Lakn(std::addressof(*LQKank[Q][nkpts+K].origin()),
                                                           {nocc_per_kp[nd][nkpts+K],nmo_per_kp[QKtok2[Q][K]]*nchol_per_kp[Q]});
              ma::product(Psi,Likn,Lakn);
              // transpose to form expected by KP3IndexFactorization  
              for(int a=0; a<nocc_per_kp[nd][nkpts+K]; a++) {
                boost::multi::array_ref<SPComplexType,2> Lkn(Lakn[a].origin(),
                                                           {nmo_per_kp[QKtok2[Q][K]],nchol_per_kp[Q]});
                boost::multi::array_ref<SPComplexType,2> Lnk(Lakn[a].origin(),
                                                           {nchol_per_kp[Q],nmo_per_kp[QKtok2[Q][K]]});
                buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}) = Lkn;
                ma::transpose(buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}),Lnk);
              }
              // doing this "by-hand" now
              Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[nd],QKtok2[Q][nkpts+K]);
              boost::multi::array_ref<SPComplexType,3> Llbn(std::addressof(*LQKlnb[Q][K].origin()),
                                                           {nmo_per_kp[K],nocc_per_kp[nd][QKtok2[Q][nkpts+K]],nchol_per_kp[Q]});
              for(int l=0; l<nmo_per_kp[K]; ++l) {
                auto psi_bj = Psi.origin();
                for(int b=0; b<nocc_per_kp[nd][QKtok2[Q][nkpts+K]]; ++b) {
                  auto Likn_jn = Likn[l].origin();
                  for(int j=0, jn=0; j<nmo_per_kp[QKtok2[Q][nkpts+K]]; ++j, ++psi_bj) {
                    auto Llbn_lbn = Llbn[l][b].origin();
                    for(int n=0; n<nchol_per_kp[Q]; ++n, ++Likn_jn, ++Llbn_lbn) 
                      (*Llbn_lbn) += (*psi_bj) * conj(*Likn_jn);
                  }
                }
              }
              for(int l=0; l<nmo_per_kp[K]; ++l) {
                boost::multi::array_ref<SPComplexType,2> Lbn(Llbn[l].origin(),
                                                           {nocc_per_kp[nd][QKtok2[Q][nkpts+K]],nchol_per_kp[Q]});
                boost::multi::array_ref<SPComplexType,2> Lnb(Llbn[l].origin(),
                                                           {nchol_per_kp[Q],nocc_per_kp[nd][QKtok2[Q][nkpts+K]]});
                buff({0,Lbn.shape()[0]},{0,Lbn.shape()[1]}) = Lbn;
                ma::transpose(buff({0,Lbn.shape()[0]},{0,Lbn.shape()[1]}),Lnb);
              }
            }
          } else {
// change get_PsiK to cast to the value_type of the result 
            auto Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[nd],K);
            assert(Psi.shape()[0] == nocc_per_kp[nd][K]);
            boost::multi::array_ref<SPComplexType,2> Likn(std::addressof(*LQKikn[Q][K].origin()),
                                                         {nmo_per_kp[K],nmo_per_kp[QKtok2[Q][K]]*nchol_per_kp[Q]});
            boost::multi::array_ref<SPComplexType,2> Lakn(std::addressof(*LQKank[Q][K].origin()),
                                                         {nocc_per_kp[nd][K],nmo_per_kp[QKtok2[Q][K]]*nchol_per_kp[Q]});
            ma::product(Psi,Likn,Lakn);
            // transpose to form expected by KP3IndexFactorization  
            for(int a=0; a<nocc_per_kp[nd][K]; a++) {
              boost::multi::array_ref<SPComplexType,2> Lkn(Lakn[a].origin(),
                                                         {nmo_per_kp[QKtok2[Q][K]],nchol_per_kp[Q]});
              boost::multi::array_ref<SPComplexType,2> Lnk(Lakn[a].origin(),
                                                         {nchol_per_kp[Q],nmo_per_kp[QKtok2[Q][K]]});
/*
              //buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}) = Lkn;
              auto b_(buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}));
              b_ = Lkn;
              ma::transpose(buff({0,Lkn.shape()[0]},{0,Lkn.shape()[1]}),Lnk);
*/
              for(int k=0; k<nmo_per_kp[QKtok2[Q][K]]; k++)
                for(int n=0; n<nchol_per_kp[Q]; n++)
                  buff[k][n] = Lkn[k][n];
              for(int k=0; k<nmo_per_kp[QKtok2[Q][K]]; k++)
                for(int n=0; n<nchol_per_kp[Q]; n++)
                  Lnk[n][k] = buff[k][n];
            }
            Psi = get_PsiK<boost::multi::array<SPComplexType,2>>(nmo_per_kp,PsiT[nd],QKtok2[Q][K]);
            // doing this "by-hand" now
            boost::multi::array_ref<SPComplexType,3> Llbn(std::addressof(*LQKlnb[Q][K].origin()),
                           {nmo_per_kp[K],nocc_per_kp[nd][QKtok2[Q][K]],nchol_per_kp[Q]});
            for(int l=0; l<nmo_per_kp[K]; ++l) {
              auto psi_bj = Psi.origin();
              for(int b=0; b<nocc_per_kp[nd][QKtok2[Q][K]]; ++b) {
                auto Likn_jn = Likn[l].origin();
                for(int j=0, jn=0; j<nmo_per_kp[QKtok2[Q][K]]; ++j, ++psi_bj) {
                  auto Llbn_lbn = Llbn[l][b].origin();
                  for(int n=0; n<nchol_per_kp[Q]; ++n, ++Likn_jn, ++Llbn_lbn) 
                    (*Llbn_lbn) += (*psi_bj) * conj(*Likn_jn);
                }
              }
            }
            for(int l=0; l<nmo_per_kp[K]; ++l) {
              boost::multi::array_ref<SPComplexType,2> Lbn(Llbn[l].origin(),
                        {nocc_per_kp[nd][QKtok2[Q][K]],nchol_per_kp[Q]});
              boost::multi::array_ref<SPComplexType,2> Lnb(Llbn[l].origin(),
                        {nchol_per_kp[Q],nocc_per_kp[nd][QKtok2[Q][K]]});
/*
              buff({0,Lbn.shape()[0]},{0,Lbn.shape()[1]}) = Lbn;
              ma::transpose(buff({0,Lbn.shape()[0]},{0,Lbn.shape()[1]}),Lnb);
*/
              for(int b=0; b<nocc_per_kp[nd][QKtok2[Q][K]]; b++)
                for(int n=0; n<nchol_per_kp[Q]; n++)
                  buff[b][n] = Lbn[b][n];
              for(int b=0; b<nocc_per_kp[nd][QKtok2[Q][K]]; b++)
                for(int n=0; n<nchol_per_kp[Q]; n++)
                  Lnb[n][b] = buff[b][n];
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
    for(int Q=0; Q<LQKank.size(); Q++) {
      TG.Cores().all_reduce_in_place_n(std::addressof(*LQKank[Q].origin()),
                                       LQKank[Q].num_elements(),std::plus<>());
      TG.Cores().all_reduce_in_place_n(std::addressof(*LQKlnb[Q].origin()),
                                       LQKlnb[Q].num_elements(),std::plus<>());
    }
    std::fill_n(std::addressof(*vn0.origin()),vn0.num_elements(),ComplexType(0.0));
  }
  TG.Node().barrier();

  // calculate (only Q=0) vn0(I,L) = -0.5 sum_K sum_j sum_n L[0][K][i][j][n] conj(L[0][K][l][j][n])
  for(int K=0; K<nkpts; K++) {
    if(K%TG.Node().size() == TG.Node().rank()) {
      boost::multi::array_ref<SPComplexType,2> Likn(std::addressof(*LQKikn[0][K].origin()),
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
  std::cout<<" EX: " <<std::setprecision(12) <<EX*2 <<std::endl;
  std::cout<<" EJ: " <<std::setprecision(12) <<EJ*4 <<std::endl;

  boost::multi::array<ComplexType,2> ABKL({nocc_max*nmo_max,nocc_max*nmo_max});
  boost::multi::array_ref<ComplexType,4> ABKL4D(ABKL.origin(),{nocc_max,nmo_max,nmo_max,nocc_max});
  EX = EJ = ComplexType(0.0);
  for(int KI=0; KI<nkpts; KI++)
  for(int KL=0; KL<nkpts; KL++) {
    for(int KK=0; KK<nkpts; KK++)
    {
      int Q = KK2Q[KI][KK];
      assert(QKtok2[Q][KI]==KK);
      int KJ = QKtok2[Q][KL];
      boost::multi::array_ref<ComplexType,3> LKI(std::addressof(*LQKank[Q][KI].origin()),{nocc_max,nchol_per_kp[Q],nmo_max});
      boost::multi::array_ref<ComplexType,3> LKL(std::addressof(*LQKlnb[Q][KL].origin()),{nmo_max,nchol_per_kp[Q],nocc_max});
      if((KI==KK && KL == KJ) || (KI==KL && KJ == KK)) { 
        for(int a=0; a<nocc_per_kp[0][0]; a++)
        for(int k=0; k<nmo_per_kp[0]; k++)
        for(int l=0; l<nmo_per_kp[0]; l++)
        for(int b=0; b<nocc_per_kp[0][0]; b++) {
          ABKL4D[a][k][l][b] = ComplexType(0.0); 
          for(int n=0; n<nchol_per_kp[Q]; n++) 
            ABKL4D[a][k][l][b] += LKI[a][n][k] * LKL[l][n][b];
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
*/

  return HamiltonianOperations(KP3IndexFactorization(TGwfn.TG_local(), type,std::move(nmo_per_kp),
            std::move(nchol_per_kp),std::move(kminus),std::move(nocc_per_kp),
            std::move(QKtok2),std::move(H1),std::move(haj),std::move(LQKikn),
            std::move(LQKank),std::move(LQKlnb),
            std::move(vn0),E0,global_ncvecs));

}

}
}
