#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
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
#include "AFQMC/Hamiltonians/KPTHCHamiltonian.h"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"
//#include "AFQMC/HamiltonianOperations/KPTHCOpsIO.hpp"

namespace qmcplusplus
{

namespace afqmc
{

HamiltonianOperations KPTHCHamiltonian::getHamiltonianOperations(bool pureSD, 
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
  long nmu, rotnmu;
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
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
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
  shmIMatrix QKtoG({nkpts,nkpts},shared_allocator<int>{TG.Node()});
  ValueType E0; 
  if( TG.Global().root() ) {
    if(!dump.read(nmo_per_kp,"NMOPerKP")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading NMOPerKP. \n";
      APP_ABORT("");
    }
    if(!dump.read(nchol_per_kp,"NCholPerKP")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading NCholPerKP. \n";
      APP_ABORT("");
    }
    if(!dump.read(kminus,"MinusK")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading MinusK. \n";
      APP_ABORT("");
    }
#define SYMMETRIZE_LQ
#ifdef SYMMETRIZE_LQ
    for(int k=0; k<nkpts; k++) {
      if(kminus[k] < k) nchol_per_kp[k] = nchol_per_kp[kminus[k]];
    }
#endif
    if(!dump.read(QKtok2,"QKTok2")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading QKTok2. \n";
      APP_ABORT("");
    }
    if(!dump.read(QKtoG,"QKToG")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading QKToG. \n";
      APP_ABORT("");
    }
#ifdef SYMMETRIZE_LQ
    // since I will force L[Q][G][u][n] = L[-Q][G][u][n]*, make sure QKToG is consistent
    std::vector<int> Gmap(nkpts);
    for(int Q=0; Q<nkpts; Q++) {
      if(kminus[Q] < Q) {
        int Qm = kminus[Q];
        for(int k=0; k<nkpts; k++) Gmap[k]=-1;
        for(int KI=0; KI<nkpts; KI++){
          int G = QKtoG[Q][KI];
          int KJ = QKtok2[Q][KI];
          int Gm = QKtoG[Qm][KJ];
          Gmap[KI] = Gm;
          // checking for consistency, any (Q,KI) that maps to G must also map (-Q,KJ) to Gm 
          for(int KI_=0; KI_<nkpts; KI_++){
            int G_ = QKtoG[Q][KI_];
            int KJ_ = QKtok2[Q][KI_];
            int Gm_ = QKtoG[Qm][KJ_];
            if( G_==G && Gm!=Gm_ )
              APP_ABORT(" Error: Inconsistent QKtoG map.\n");
          }
        }
        for(int KI=0; KI<nkpts; KI++)
          QKtoG[Q][KI] = Gmap[KI];
      }   
    }
#endif
    std::vector<ValueType> E_(2);
    if(!dump.read(E_,"Energies")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Energies. \n";
      APP_ABORT("");
    }
    E0 = E_[0]+E_[1];
    E0 = nkpts*(E0-1.34789140434);
    if(nmo_per_kp.size() != nkpts ||
       nchol_per_kp.size() != nkpts ||
       kminus.size() != nkpts ||
       QKtok2.shape()[0] != nkpts ||  
       QKtok2.shape()[1] != nkpts || 
       QKtoG.shape()[0] != nkpts ||  
       QKtoG.shape()[1] != nkpts  
      ) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Inconsistent dimension (NMOPerKP,NCholPerKP,QKtTok2,QKToG): " 
                 <<nkpts <<" " 
                 <<nmo_per_kp.size() <<" " 
                 <<nchol_per_kp.size() <<" " 
                 <<kminus.size() <<" "
                 <<QKtok2.shape()[0] <<" " 
                 <<QKtok2.shape()[1] <<" "
                 <<QKtoG.shape()[0] <<" " 
                 <<QKtoG.shape()[1] <<std::endl; 
      APP_ABORT("");
    }
  }
  TG.Global().broadcast_n(&E0,1,0);
  TG.Global().broadcast_n(nmo_per_kp.begin(),nmo_per_kp.size(),0);
  TG.Global().broadcast_n(nchol_per_kp.begin(),nchol_per_kp.size(),0);
  TG.Global().broadcast_n(kminus.begin(),kminus.size(),0);
  if(TG.Node().root()) { 
    TG.Cores().broadcast_n(std::addressof(*QKtok2.origin()),QKtok2.num_elements(),0);
    TG.Cores().broadcast_n(std::addressof(*QKtoG.origin()),QKtoG.num_elements(),0);
  }
  TG.Node().barrier();

  Idata.resize(2);
  int nmo_max = *std::max_element(nmo_per_kp.begin(),nmo_per_kp.end());
  int nmo_tot = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.end(),0);
  int nchol_max = *std::max_element(nchol_per_kp.begin(),nchol_per_kp.end());
  shmCTensor H1({nkpts,nmo_max,nmo_max},shared_allocator<ComplexType>{TG.Node()}); 
  if( TG.Node().root() ) {
    // now read H1_kpQ 
    for(int Q=0; Q<nkpts; Q++) {
      // until double_hyperslabs work!
      boost::multi::array<ComplexType,2> h1({nmo_per_kp[Q],nmo_per_kp[Q]}); 
      if(!dump.read(h1,std::string("H1_kp")+std::to_string(Q))) {
        app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading /Hamiltonian/H1_kp" <<Q <<". \n";
        APP_ABORT("");
      }
      H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}) = h1; 
      //H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}) = h1; 
    }    

    // read LQ
    if(!dump.push("KPTHC",false)) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Group KPTHC not found. \n";
      APP_ABORT("");
    }

    if(!dump.read(Idata,"dims")) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading KPTHC/dims. \n";
      APP_ABORT("");
    }
  }
  TG.Node().broadcast_n(Idata.begin(),2,0);
  nmu = Idata[0];
  rotnmu = 260;//Idata[1];  

  std::vector<shmSpMatrix> LQGun;
  LQGun.reserve(nkpts);
  shmSpMatrix Piu({nmo_tot,nmu},shared_allocator<SPComplexType>{TG.Node()});
  for(int Q=0; Q<nkpts; Q++) {
    int nG = *std::max_element(QKtoG[Q].begin(),QKtoG[Q].end())+1;  
    LQGun.emplace_back( shmSpMatrix({nG*nmu,nchol_per_kp[Q]},
                                   shared_allocator<SPComplexType>{TG.Node()}) );
  }  
  std::vector<shmSpMatrix> rotMuv;
  rotMuv.reserve(nkpts);
  shmSpMatrix rotPiu({nmo_tot,rotnmu},shared_allocator<SPComplexType>{TG.Node()});
  for(int Q=0; Q<nkpts; Q++) {
    int nG = *std::max_element(QKtoG[Q].begin(),QKtoG[Q].end())+1;  
    rotMuv.emplace_back( shmSpMatrix({nG*rotnmu,nG*rotnmu},
                                   shared_allocator<SPComplexType>{TG.Node()}) );
  }  
  if( TG.Node().root() ) {
    for(int Q=0; Q<nkpts; Q++) {
      int nG = *std::max_element(QKtoG[Q].begin(),QKtoG[Q].end())+1;  
#ifdef SYMMETRIZE_LQ
      using std::conj;
      if(kminus[Q] < Q) {
        int Qm = kminus[Q];
        int nGm = *std::max_element(QKtoG[Qm].begin(),QKtoG[Qm].end())+1;  
        if(nG != nGm)
          APP_ABORT(" Error: nG != nGm. \n");
        auto LQ(std::addressof(*LQGun[Q].origin()));
        auto LQm(std::addressof(*LQGun[Qm].origin()));
        for(int Gun=0; Gun<LQGun[Q].num_elements(); Gun++, ++LQ, ++LQm) 
          (*LQ) = conj(*LQm);
      } else
#endif
      {  
        if(!dump.read(LQGun[Q],std::string("L")+std::to_string(Q))) {
          app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading /Hamiltonian/KPTHC/L" <<Q <<". \n";
          APP_ABORT("");
        }
      }
      if(LQGun[Q].shape()[0] != nG*nmu || LQGun[Q].shape()[1] != nchol_per_kp[Q]) {
        app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading /Hamiltonian/KPTHC/L" <<Q <<". \n"
                 <<" Unexpected dimensins: " <<LQGun[Q].shape()[0] <<" " <<LQGun[Q].shape()[1] <<std::endl; 
        APP_ABORT("");
      }    
// adding scaling factor until Fionn does it on his end
/*
      {
        RealType scl = 1.0/std::sqrt(double(nkpts));
        auto LQ(std::addressof(*LQGun[Q].origin()));
        for(int Gun=0; Gun<LQGun[Q].num_elements(); Gun++, ++LQ)
          (*LQ) *= scl; 
      }
*/
      if(!dump.read(rotMuv[Q],std::string("HalfTransformedMuv")+std::to_string(Q))) {
        app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading /Hamiltonian/KPTHC/HalfTransformedMuv0" <<Q <<". \n";
        APP_ABORT("");
      }  
      if(rotMuv[Q].shape()[0] != nG*rotnmu || rotMuv[Q].shape()[1] != nG*rotnmu) {
        app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading /Hamiltonian/KPTHC/L" <<Q <<". \n"
                 <<" Unexpected dimensins: " <<rotMuv[Q].shape()[0] <<" " <<rotMuv[Q].shape()[1] <<std::endl;
        APP_ABORT("");
      }
    }
    if(!dump.read(Piu,std::string("Orbitals"))) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
               <<" Problems reading /Hamiltonian/KPTHC/Orbitals. \n";
      APP_ABORT("");
    }
    if(Piu.shape()[0] != nmo_tot || Piu.shape()[1] != nmu) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
               <<" Problems reading /Hamiltonian/KPTHC/Orbitals. \n"
               <<" Unexpected dimensins: " <<Piu.shape()[0] <<" " <<Piu.shape()[1] <<std::endl; 
      APP_ABORT("");
    }    
    if(!dump.read(rotPiu,std::string("HalfTransformedFullOrbitals"))) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
               <<" Problems reading /Hamiltonian/KPTHC/HalfTransformedFullOrbitals. \n";
      APP_ABORT("");
    }
    if(rotPiu.shape()[0] != nmo_tot || rotPiu.shape()[1] != rotnmu) {
      app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
               <<" Problems reading /Hamiltonian/KPTHC/Orbitals. \n"
               <<" Unexpected dimensins: " <<rotPiu.shape()[0] <<" " <<rotPiu.shape()[1] <<std::endl;
      APP_ABORT("");
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
          app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Only wavefunctions in block-diagonal form are accepted. \n"; 
          APP_ABORT("");
        }
        if(not get_nocc_per_kp(nmo_per_kp,PsiT[2*i+1],nocc_per_kp[i]({nkpts,2*nkpts}))) {
          app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Only wavefunctions in block-diagonal form are accepted. \n"; 
          APP_ABORT("");
        }
      }  
    } else {
      for(int i=0; i<ndet; i++) 
        if(not get_nocc_per_kp(nmo_per_kp,PsiT[i],nocc_per_kp[i])) {
          app_error()<<" Error in KPTHCHamiltonian::getHamiltonianOperations():"
                 <<" Only wavefunctions in block-diagonal form are accepted. \n"; 
          APP_ABORT("");
        }
    }
  }
  TG.Node().barrier();
  int nocc_max = *std::max_element(std::addressof(*nocc_per_kp.origin()),
                                   std::addressof(*nocc_per_kp.origin())+nocc_per_kp.num_elements());

  /* half-rotate Piu and H1:
   * Given that PsiT = H(SM), 
   * h[K][a][k] = sum_i PsiT[K][a][i] * h[K][i][k]
   * cPua[u][K][a] = sum_i PsiT[K](a,i) * Piu[K][i][u]
   */
  std::vector<shmSpMatrix> cPua;
  cPua.reserve(ndet*nspins);   
  std::vector<shmSpMatrix> rotcPua;
  rotcPua.reserve(ndet*nspins);   
  shmCMatrix haj({ndet*nkpts,(type==COLLINEAR?2:1)*nocc_max*nmo_max},shared_allocator<ComplexType>{TG.Node()});
  std::pair<int,int> nel;
  nel.first = std::accumulate(std::addressof(*nocc_per_kp[0].origin()),std::addressof(*nocc_per_kp[0].origin())+nkpts,0);
  if(type==COLLINEAR) 
    nel.second = std::accumulate(std::addressof(*nocc_per_kp[0].origin())+nkpts,
                                 std::addressof(*nocc_per_kp[0].origin())+2*nkpts,0);
  for(int nd=0; nd<ndet; nd++) {
    cPua.emplace_back(shmSpMatrix({nmu,nel.first},shared_allocator<SPComplexType>{TG.Node()}));
    if(type==COLLINEAR) 
      cPua.emplace_back(shmSpMatrix({nmu,nel.second},shared_allocator<SPComplexType>{TG.Node()}));
    rotcPua.emplace_back(shmSpMatrix({rotnmu,nel.first},shared_allocator<SPComplexType>{TG.Node()}));
    if(type==COLLINEAR) 
      rotcPua.emplace_back(shmSpMatrix({rotnmu,nel.second},shared_allocator<SPComplexType>{TG.Node()}));
  }
  if(TG.Node().root()) {
    std::fill_n(haj.origin(),haj.num_elements(),ComplexType(0.0));
    for(auto& v:cPua)
      std::fill_n(v.origin(),v.num_elements(),SPComplexType(0.0));
    for(auto& v:rotcPua)
      std::fill_n(v.origin(),v.num_elements(),SPComplexType(0.0));
  }
  TG.Node().barrier();

  boost::multi::array<SPComplexType,2> buff({nmo_max,nchol_max});  
  for(int nd=0, nt=0; nd<ndet; nd++) {
    for(int Q=0; Q<nkpts; Q++) {
      if((nt++)%TG.Global().size() == TG.Global().rank()) {
        // haj  
        if(type==COLLINEAR) {
          { // Alpha 
            auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd],Q); 
            assert(Psi.shape()[0] == nocc_per_kp[nd][Q]);
            boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+Q].origin()),
                                                         {nocc_per_kp[nd][Q],nmo_per_kp[Q]});
            ma::product(Psi,H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}),haj_r);
          }
          { // Beta 
            auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd+1],Q);
            assert(Psi.shape()[0] == nocc_per_kp[nd][nkpts+Q]);
            boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+Q].origin())+
                                                                        nocc_per_kp[nd][Q]*nmo_per_kp[Q],
                                                         {nocc_per_kp[nd][nkpts+Q],nmo_per_kp[Q]});
            ma::product(Psi,H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}),haj_r);
          }
        } else {
          auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[nd],Q);
          assert(Psi.shape()[0] == nocc_per_kp[nd][Q]);
          boost::multi::array_ref<ComplexType,2> haj_r(std::addressof(*haj[nd*nkpts+Q].origin()),
                                                       {nocc_per_kp[nd][Q],nmo_per_kp[Q]});
          ma::product(ComplexType(2.0),Psi,H1[Q]({0,nmo_per_kp[Q]},{0,nmo_per_kp[Q]}),
                      ComplexType(0.0),haj_r);
        }
      }   
      if((nt++)%TG.Global().size() == TG.Global().rank()) {
        if(type==COLLINEAR) {
          APP_ABORT(" Finish .\n"); 
          auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd],Q);
          int ne0 = std::accumulate(std::addressof(*nocc_per_kp[nd].origin()),std::addressof(*nocc_per_kp[nd].origin())+Q,0);
          int ni0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+Q,0);
          assert(Psi.shape()[0] == nocc_per_kp[nd][Q]);
          ma::product(ma::H(Piu({ni0,ni0+nmo_per_kp[Q]},{0,nmu})),ma::T(Psi),
                      cPua[2*nd]({0,nmu},{ne0,ne0+Psi.shape()[0]}));
          ma::product(ma::H(rotPiu({ni0,ni0+nmo_per_kp[Q]},{0,rotnmu})),ma::T(Psi),
                      rotcPua[2*nd]({0,rotnmu},{ne0,ne0+Psi.shape()[0]}));

          Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[2*nd+1],Q);
          ne0 = std::accumulate(std::addressof(*nocc_per_kp[nd].origin())+nkpts,
                                std::addressof(*nocc_per_kp[nd].origin())+nkpts+Q,0);
          assert(Psi.shape()[0] == nocc_per_kp[nd][nkpts+Q]);
          ma::product(ma::H(Piu({ni0,ni0+nmo_per_kp[Q]},{0,nmu})),ma::T(Psi),
                      cPua[2*nd+1]({0,nmu},{ne0,ne0+Psi.shape()[0]}));
          assert(Psi.shape()[0] == nocc_per_kp[nd][nkpts+Q]);
          ma::product(ma::H(rotPiu({ni0,ni0+nmo_per_kp[Q]},{0,rotnmu})),ma::T(Psi),
                      rotcPua[2*nd+1]({0,rotnmu},{ne0,ne0+Psi.shape()[0]}));

        } else {
          int ne0 = std::accumulate(std::addressof(*nocc_per_kp[nd].origin()),std::addressof(*nocc_per_kp[nd].origin())+Q,0);  
          int ni0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+Q,0);  
          auto Psi = get_PsiK<boost::multi::array<ComplexType,2>>(nmo_per_kp,PsiT[nd],Q);
          assert(Psi.shape()[0] == nocc_per_kp[nd][Q]);
          ma::product(ma::H(Piu({ni0,ni0+nmo_per_kp[Q]},{0,nmu})),ma::T(Psi),
                      cPua[nd]({0,nmu},{ne0,ne0+Psi.shape()[0]}));
          ma::product(ma::H(rotPiu({ni0,ni0+nmo_per_kp[Q]},{0,rotnmu})),ma::T(Psi),
                      rotcPua[nd]({0,rotnmu},{ne0,ne0+Psi.shape()[0]}));
        }
      }
    }
  }
  TG.Global().barrier();
  if(TG.Node().root()) {
    TG.Cores().all_reduce_in_place_n(std::addressof(*haj.origin()),
                                     haj.num_elements(),std::plus<>());
    for(int n=0; n<cPua.size(); n++) 
      TG.Cores().all_reduce_in_place_n(std::addressof(*cPua[n].origin()),
                                       cPua[n].num_elements(),std::plus<>());
    for(int n=0; n<rotcPua.size(); n++) 
      TG.Cores().all_reduce_in_place_n(std::addressof(*rotcPua[n].origin()),
                                       rotcPua[n].num_elements(),std::plus<>());
    std::fill_n(std::addressof(*vn0.origin()),vn0.num_elements(),ComplexType(0.0));
  }
  TG.Node().barrier();

/*
  // calculate (only Q=0) vn0(I,L) = -0.5 sum_K sum_j sum_n L[0][K][i][j][n] conj(L[0][K][l][j][n])
  for(int K=0; K<nkpts; K++) {
    if(K%TG.Node().size() == TG.Node().rank()) {
      boost::multi::array_ref<SPComplexType,2> Likn(std::addressof(*LQKikn[0][K].origin()),
                                                   {nmo_per_kp[K],nmo_per_kp[K]*nchol_per_kp[0]});        
      using ma::H;
      ma::product(-0.5,Likn,H(Likn),0.0,vn0[K]({0,nmo_per_kp[K]},{0,nmo_per_kp[K]}));
    }
  }
*/
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
  app_log()<<" E1+E0: " <<std::setprecision(12) <<E0+2*E_ <<"  ";
  E_=0.0; 
  for(int K=0; K<nkpts; K++) {
    auto&& G_=Gc[K];
    for(int a=0, aj=0; a<nocc_per_kp[0][K]; a++)
      for(int j=0; j<nmo_per_kp[K]; j++, ++aj)
        E_ += haj[K][aj]*G_[a][j];
  }
  //app_log()<<nkpts*(E0-1.34789140434)+E_ <<std::endl;
  app_log()<<E0+E_ <<std::endl;

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
/*
  // LIJ[I][J][Q][n] = sum_u u[KI][i][u] u[KJ][j][u] LQGun[Q][G][u][n]
  boost::multi::array<ComplexType,5> LIJ({nkpts,nkpts,nmo_max,nmo_max,nchol_max});
  boost::multi::array_ref<ComplexType,4> LIJ4D(LIJ.origin(),{nkpts,nkpts,nmo_max*nmo_max,nchol_max});
  std::fill_n(LIJ.origin(),LIJ.num_elements(),0);
  boost::multi::array<ComplexType,2> uij({nmo_max*nmo_max,nmu});
  for(int KI=0, n_=0; KI<nkpts; KI++)
  for(int KJ=0; KJ<nkpts; KJ++) {  
    if((n_++)%TGwfn.Global().size() == TGwfn.Global().rank()) {
//Timer.reset("T0");
//Timer.start("T0");
      int Q = KK2Q[KI][KJ];
      int G1 = QKtoG[Q][KI];
      int ni0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+KI,0);
      int nj0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+KJ,0);
      for(int i=0, I=ni0; i<nmo_per_kp[0]; i++, ++I) {
        auto&& uI(Piu[ni0+i]);
        for(int j=0, J=nj0; j<nmo_per_kp[0]; j++, ++J) {
          auto&& uJ(Piu[nj0+j]);
          auto LIJ_(LIJ[KI][KJ][i][j].origin());  
          for(int u=0; u<nmu; u++) { 
            ComplexType uij = conj(uI[u])*uJ[u];
            auto lq(std::addressof(*LQGun[Q][nmu*G1+u].origin()));
            for(int n=0; n<nchol_per_kp[Q]; n++) 
              LIJ_[n] += uij * lq[n]; 
          }
        }
      }
//Timer.stop("T0");
//app_log()<<KI <<" " <<KJ <<" " <<Timer.total("T0") <<std::endl;
    }
  } 
  TGwfn.Global().all_reduce_n(LIJ.origin(), LIJ.num_elements(), std::plus<>());

  boost::multi::array<ComplexType,2> IJKL({nmo_max*nmo_max,nmo_max*nmo_max});
  boost::multi::array<ComplexType,2> Muv({nmu,nmu});
  boost::multi::array_ref<ComplexType,4> IJKL4D(IJKL.origin(),{nmo_max,nmo_max,nmo_max,nmo_max});
  ComplexType EX(0.0);
  ComplexType EJ(0.0);
  myTimer Timer;
  Timer.reset("T0");
  for(int KI=0, n_=0; KI<nkpts; KI++)
  for(int KL=0; KL<nkpts; KL++) { 
   for(int KK=0; KK<nkpts; KK++) 
   {
    int Q = KK2Q[KI][KK];
    int KJ = QKtok2[Q][KL];
    if( not (KI==KK && KL == KJ) ) continue;
    if((n_++)%TGwfn.Global().size() == TGwfn.Global().rank()) {    
      int ni0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+KI,0);  
      int nj0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+KJ,0);  
      int nk0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+KK,0);  
      int nl0 = std::accumulate(nmo_per_kp.begin(),nmo_per_kp.begin()+KL,0);  
      //boost::multi::array_ref<ComplexType,2> LKI(std::addressof(*LIJ[KI][KK].origin()),{nmo_max*nmo_max,nchol_per_kp[Q]}); 
      //boost::multi::array_ref<ComplexType,2> LKL(std::addressof(*LIJ[KL][KJ].origin()),{nmo_max*nmo_max,nchol_per_kp[Q]}); 
      ma::product(LIJ4D[KI][KK]({0,nmo_max*nmo_max},{0,nchol_per_kp[Q]}),
                  ma::H(LIJ4D[KL][KJ]({0,nmo_max*nmo_max},{0,nchol_per_kp[Q]})),IJKL);
//Timer.stop("T0");
//app_log()<<i <<" " <<Timer.total("T0") <<std::endl;
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
  }
  TGwfn.Global().all_reduce_n(&EX, 1, std::plus<>());
  TGwfn.Global().all_reduce_n(&EJ, 1, std::plus<>());
  app_log()<<" EX: " <<std::setprecision(12) <<EX*2 <<std::endl;
  app_log()<<" EJ: " <<std::setprecision(12) <<EJ*4 <<std::endl;
*/

  return HamiltonianOperations(KPTHCOps(TGwfn.TG_local(),type,
                            std::move(nmo_per_kp),std::move(nchol_per_kp),
                            std::move(kminus),std::move(nocc_per_kp),
                            std::move(QKtok2),std::move(QKtoG),
                            std::move(H1),std::move(haj),
                            std::move(rotMuv),std::move(rotPiu),std::move(rotcPua),
                            std::move(LQGun),std::move(Piu),std::move(cPua),
                            std::move(vn0),E0,global_ncvecs));

}

}
}
