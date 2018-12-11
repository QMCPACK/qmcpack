#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<map>
#include<utility>
#include<vector>
#include<numeric>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "Configuration.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Hamiltonians/SymmetricFactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/HSPotential_Helpers.h"
#include "AFQMC/Matrix/array_partition.hpp"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"

namespace qmcplusplus
{

namespace afqmc
{

SpVType_shm_csr_matrix SymmetricFactorizedSparseHamiltonian::calculateHSPotentials(double cut, 
        TaskGroup_& TGprop, boost::multi_array<ComplexType,2>& vn0) {

  using Alloc = boost::mpi3::intranode::allocator<SPValueType>;
  if(TG.getNumberOfTGs() > 1)
    APP_ABORT("Error: HSPotential not implemented with distributed Hamiltonian. \n");

  vn0.resize(extents[NMO][NMO]);
  std::fill_n(vn0.data(),NMO*NMO,ComplexType(0));
  for(int i=0, cnt=0; i<NMO; i++)
    for(int l=i; l<NMO; l++, cnt++) {
      if(cnt%TG.Global().size() != TG.Global().rank()) continue;
      ValueType vl = ValueType(0);
      for(int j=0; j<NMO; j++)
        vl += H(i,j,j,l);
      vn0[i][l] -= 0.5*vl;
      if(i!=l) vn0[l][i] -= 0.5*myconj(vl);
    }
  TG.Global().all_reduce_in_place_n(vn0.origin(),vn0.num_elements(),std::plus<>());

  if(TG.getNumberOfTGs() > 1) {

    APP_ABORT(" Finish HS. \n");
    return SpVType_shm_csr_matrix({0,0},{0,0},0,Alloc(TG.Node()));

  } else {

    // you always need to count since some vectors might be empty 
    auto nnz_per_cv = HamHelperSymmetric::count_nnz_per_cholvec(cut,TG,V2_fact,NMO);

    // partition and 
    std::size_t cv0, cvN;
    if(TGprop.getNNodesPerTG() == 1 ) { // Spvn is not distributed
      cv0 = 0;
      cvN = nnz_per_cv.size();
      if(TG.Global().root()) {
        app_log()<<std::endl <<"Partition of cholesky vector over nodes in TG: ";
        app_log()<<std::count_if(nnz_per_cv.begin(),
                                 nnz_per_cv.begin()+cvN,
                                 [] (std::size_t i) { return i>0; } );
        app_log()<<std::endl;
        app_log()<<"Number of terms in Cholesky Matrix per node in TG: ";
        app_log()<<std::accumulate(nnz_per_cv.begin(),
                                   nnz_per_cv.begin()+cvN,std::size_t(0));
        app_log()<<std::endl;
      }
    } else {
      std::vector<std::size_t> cv_boundaries(TGprop.getNNodesPerTG()+1);
      simple_matrix_partition<TaskGroup_,std::size_t,double> split(V2_fact.shape()[0],
                                                                   nnz_per_cv.size(),cut);
      // no need for all cores to do this
      if(TG.Global().root())
        split.partition(TGprop,false,nnz_per_cv,cv_boundaries);
      TG.Global().broadcast_n(cv_boundaries.begin(),cv_boundaries.size());  
      cv0 = cv_boundaries[TGprop.getLocalNodeNumber()];  
      cvN = cv_boundaries[TGprop.getLocalNodeNumber()+1];  
      // no need for all cores to do this
      if(TG.Global().root()) {
        app_log()<<std::endl <<"Partition of cholesky vector over nodes in TG: ";
        for(int i=0; i<TGprop.getNNodesPerTG(); i++)
            app_log()<<std::count_if(nnz_per_cv.begin()+cv_boundaries[i],
                       nnz_per_cv.begin()+cv_boundaries[i+1],
                       [] (std::size_t i) { return i>0; } ) <<" ";
        app_log()<<std::endl;    
        app_log()<<"Number of terms in Cholesky Matrix per node in TG: ";
        for(int i=0; i<TGprop.getNNodesPerTG(); i++)
            app_log()<<std::accumulate(nnz_per_cv.begin()+cv_boundaries[i],
                       nnz_per_cv.begin()+cv_boundaries[i+1],std::size_t(0)) <<" ";
        app_log()<<std::endl <<std::endl;    
      }
    }

    auto nnz_per_ik = HamHelperSymmetric::count_nnz_per_ik(cut,TG,V2_fact,NMO,cv0,cvN); 

    std::size_t nvec = std::count_if(nnz_per_cv.begin()+cv0,nnz_per_cv.begin()+cvN,
               [] (std::size_t const& i) { return i>0; } );

    std::size_t cv_origin = std::count_if(nnz_per_cv.begin(),nnz_per_cv.begin()+cv0,
               [] (std::size_t const& i) { return i>0; } );

    // can build csr directly since cores work on non-overlapping rows
    // and can use emplace_back
    SpVType_shm_csr_matrix csr({NMO*NMO,nvec}, {0,cv_origin}, nnz_per_ik, Alloc(TG.Node()));

    HamHelperSymmetric::generateHSPotential(csr,cut,TG,V2_fact,NMO,cv0,cvN); 
    TG.node_barrier();

    return csr;
  }

}

 SpCType_shm_csr_matrix SymmetricFactorizedSparseHamiltonian::halfRotatedHijkl(WALKER_TYPES type, TaskGroup_& TGHam, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, RealType const cut){
    check_wavefunction_consistency(type,Alpha,Beta,NMO,NAEA,NAEB);
    using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
    std::size_t nr=NAEA*NMO;
    if(type==COLLINEAR) nr = (NAEA+NAEB)*NMO;
    else if(type==NONCOLLINEAR) nr = (NAEA+NAEB)*2*NMO;
    if(TGHam.getNNodesPerTG() > 1) {
      using tvec = std::vector<std::tuple<int,int,SPComplexType>>;
      tvec tmat;
      tmat.reserve(100000); // reserve some space  
      rotateHijklSymmetric<tvec>(type,TG,tmat,Alpha,Beta,V2_fact,
            cut,maximum_buffer_size,false,false);
      TG.node_barrier();
      return csr::shm::construct_distributed_csr_matrix_from_distributed_containers<SpCType_shm_csr_matrix>(tmat,nr,nr,TGHam);
    } else {
      using ucsr_mat = SpCType_shm_csr_matrix::base;
      ucsr_mat ucsr({nr,nr},{0,0},0,Alloc(TG.Node()));
      rotateHijklSymmetric<ucsr_mat>(type,TG,ucsr,Alpha,Beta,
            V2_fact,cut,maximum_buffer_size,true,true);
      TG.node_barrier();
      return csr::shm::construct_csr_matrix_from_distributed_ucsr<SpCType_shm_csr_matrix,TaskGroup_>(std::move(ucsr),TG);
    }
  }

  HamiltonianOperations SymmetricFactorizedSparseHamiltonian::getHamiltonianOperations(
            bool pureSD, WALKER_TYPES type, std::vector<PsiT_Matrix>& PsiT, 
            double cutvn, double cutv2, TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& dump) {

    if(type==COLLINEAR)
      assert(PsiT.size()%2 == 0);

    Timer.reset("Generic");
    Timer.start("Generic");   
    boost::multi_array<ComplexType,2> vn0(extents[NMO][NMO]);
    auto Spvn(std::move(calculateHSPotentials(cutvn,TGprop,vn0)));
    auto Spvnview(csr::shm::local_balanced_partition(Spvn,TGprop));
    Timer.stop("Generic");
    app_log()<<" Time for calculateHSPotentials: " <<Timer.total("Generic") <<std::endl;

    // trick: the last node always knows what the total # of chol vecs is
    int global_ncvecs=0;
    if(TG.getNodeID() == TG.getTotalNodes()-1 && TG.getCoreID()==0)
      global_ncvecs = Spvn.global_origin()[1] + Spvn.shape()[1];
    global_ncvecs = TG.Global().all_reduce_value(global_ncvecs,std::plus<>());

    ValueType E0 = OneBodyHamiltonian::NuclearCoulombEnergy + 
                   OneBodyHamiltonian::FrozenCoreEnergy;

    // several posibilities
    int ndet = ((type!=COLLINEAR)?(PsiT.size()):(PsiT.size()/2));
    // SparseTensor<Integrals_Type, SpvnT_Type>
    if(ndet==1) {

      Timer.reset("Generic");
      Timer.start("Generic");   
      std::vector<boost::multi_array<SPComplexType,1>> hij;
      hij.reserve(ndet);
      hij.emplace_back(halfRotatedHij(type,&PsiT[0],((type==COLLINEAR)?(&PsiT[1]):(&PsiT[0])))); 
      auto SpvnT(sparse_rotate::halfRotateCholeskyMatrixForBias(type,TGprop,
                              &PsiT[0],((type==COLLINEAR)?(&PsiT[1]):(&PsiT[0])),
                              Spvn,cutv2));
      auto SpvnTview(csr::shm::local_balanced_partition(SpvnT,TGprop));
      Timer.stop("Generic");
      app_log()<<" Time for halfRotateCholeskyMatrixForBias: " <<Timer.total("Generic") <<std::endl;

      // in single determinant, SpvnT is the half rotated cholesky matrix 
      if(pureSD) {
        // in pureSD: V2 is ValueType
        using sparse_ham = SparseTensor<ValueType,ComplexType>;

        return HamiltonianOperations{};
        //return HamiltonianOperations(sparse_ham(getH1(),std::move(hij),std::move(V2),
        //    std::move(V2view),std::move(Spvn),std::move(Spvnview),
        //    std::move(vn0),std::move(SpvnT),std::move(SpvnTview),E0,global_ncvecs));
      } else {
        // in this case: V2 is ComplexType since it is half rotated
        using sparse_ham = SparseTensor<ComplexType,ComplexType>;

        Timer.reset("Generic");
        Timer.start("Generic");   
        std::vector<SpCType_shm_csr_matrix> V2;
        V2.reserve(ndet);
        V2.emplace_back(halfRotatedHijkl(type,TGwfn,&PsiT[0],
                                           ((type==COLLINEAR)?(&PsiT[1]):(&PsiT[0])),cutv2));
        std::vector<SpCType_shm_csr_matrix::template matrix_view<int>> V2view;
        V2view.reserve(ndet);
        V2view.emplace_back(csr::shm::local_balanced_partition(V2[0],TGwfn));
        Timer.stop("Generic");
        app_log()<<" Time for halfRotateHijkl: " <<Timer.total("Generic") <<std::endl;

        return HamiltonianOperations(sparse_ham(type,getH1(),std::move(hij),std::move(V2),
            std::move(V2view),std::move(Spvn),std::move(Spvnview),
            std::move(vn0),std::move(SpvnT),std::move(SpvnTview),E0,global_ncvecs,true));
      }
    } else {
      // in multi determinant, SpvnT is transposed(Spvn) 

      Timer.reset("Generic");
      Timer.start("Generic");   
      std::vector<boost::multi_array<SPComplexType,1>> hij;
      hij.reserve(ndet);
      int skp=((type==COLLINEAR)?1:0); 
      for(int n=0, nd=0; n<ndet; ++n, nd+=(skp+1)) 
         hij.emplace_back(halfRotatedHij(type,&PsiT[nd],&PsiT[nd+skp])); 
      auto SpvnT(csr::shm::transpose(Spvn));
      auto SpvnTview(csr::shm::local_balanced_partition(SpvnT,TGprop));
      Timer.stop("Generic");
      app_log()<<" Time for halfRotateCholeskyMatrixForBias: " <<Timer.total("Generic") <<std::endl;

      if(pureSD) {
        // in pureSD: V2 is ValueType
        using sparse_ham = SparseTensor<ValueType,ValueType>;

        return HamiltonianOperations{};
        //return HamiltonianOperations(sparse_ham(getH1(),std::move(hij),std::move(V2),
        //    std::move(V2view),std::move(Spvn),std::move(Spvnview),
        //    std::move(vn0),std::move(SpvnT),std::move(SpvnTview),E0,global_ncvecs));
      } else {
        // in this case: V2 is ComplexType since it is half rotated
        using sparse_ham = SparseTensor<ComplexType,ValueType>;

        Timer.reset("Generic");
        Timer.start("Generic");   
        std::vector<SpCType_shm_csr_matrix> V2;
        V2.reserve(ndet);
        for(int n=0, nd=0; n<ndet; ++n, nd+=(skp+1)) 
          V2.emplace_back(halfRotatedHijkl(type,TGwfn,&PsiT[nd],&PsiT[nd+skp],cutv2));
        std::vector<SpCType_shm_csr_matrix::template matrix_view<int>> V2view;
        V2view.reserve(ndet);
        for(auto& v:V2) 
          V2view.emplace_back(csr::shm::local_balanced_partition(v,TGwfn));
        Timer.stop("Generic");
        app_log()<<" Time for halfRotateHijkl (for all dets): " <<Timer.total("Generic") <<std::endl;

        return HamiltonianOperations(sparse_ham(type,getH1(),std::move(hij),std::move(V2),
            std::move(V2view),std::move(Spvn),std::move(Spvnview),
            std::move(vn0),std::move(SpvnT),std::move(SpvnTview),E0,global_ncvecs,true));
      }
    }

  }

/*
  bool FactorizedSparseHamiltonian::createHamiltonianForPureDeterminant(int walker_type, bool aa_only, std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b , std::vector<s1D<ValueType> >& hij, SPValueSMSpMat& Vijkl, const RealType cut)
  {
    if(walker_type==2) {
      APP_ABORT("Error: GHF density matrix only implemented with spinRestricted integrals. \n");
    }

    if(distribute_Ham) {
      APP_ABORT("Error: createHamiltonianForPureDeterminant not implemented with distributed Hamiltonian. \n");
    }

    // teporary until mpi3 is fully integrated
    TaskGroup_ TGham(TG,"DummyHS",1,TG.getNCoresPerTG());

    createHij(walker_type,NMO,occ_a,occ_b,H1,hij,cut);

    Timer.reset("Generic");
    Timer.start("Generic");

    ValueType zero = ValueType(0);

    std::vector<s4D<ValueType> > vs4D;
    vs4D.reserve(24);

    std::size_t cnt2=0;
    long cnter=0, npr = TG.getGlobalSize(), rk = TG.getGlobalRank();
    OrbitalType i,j,k,l,j1,k1,l1,j2,k2,l2;
    ValueType J1,J2,J3,J1a=zero,J2a=zero,J3a=zero,fct;

    DiagHam.resize(NMO,NMO);
    for(IndexType i=0; i<NMO; i++)
    for(IndexType k=i; k<NMO; k++, cnter++) {
      if( cnter%npr != rk ) continue;
      DiagHam(i,k) =  H(i,k,k,i);
      DiagHam(k,i) = DiagHam(i,k);
#if defined(QMC_COMPLEX)
      if(DiagHam(i,k).imag() > 1e-8) {
          app_error()<<" Error: Found complex diagonal on hamiltonian. " <<i <<" " <<k <<" " <<DiagHam(i,k) <<std::endl;
          APP_ABORT("Error: Found complex diagonal on hamiltonian.");
      }
      ComplexType Hik = H(i,k,i,k);
      if(Hik.imag() > cutoff_cholesky) {
          app_error()<<" Error: Problems with factorized Hamiltonian <ik|ik> should be real: "
                     <<i <<" " <<k <<" " <<Hik <<std::endl;
          APP_ABORT("Error: Problems with factorized hamiltonian.");
      }
#endif
    }
    {
      ValueMatrix dum(DiagHam);  
      TG.Global().all_reduce(dum.begin(),dum.end(),DiagHam.begin()); 
    }

    int occi, occj, occk, occl;
    cnter=0;
    // Approximation: (similar to non factorized case)
    //   - if <ij|kl>  <  cutoff1bar, set to zero 
    for(IndexType i=0; i<NMO; i++)  {
    for(IndexType j=i; j<NMO; j++,cnter++)  {
     if( cnter%npr != rk ) continue;
     occi = (occ_a[i]||occ_b[i])?1:0;
     occj = (occ_a[j]||occ_b[j])?1:0;
    for(IndexType k=j; k<NMO; k++)  {
     occk = (occ_a[k]||occ_b[k])?1:0;
    for(IndexType l=k; l<NMO; l++)  {

        occl = (occ_a[l]||occ_b[l])?1:0;
        if( occi + occj + occk + occl < 2 ) continue;

        // NOTE NOTE NOTE:
        // <ik|ik> can get a complex component due to truncation error, eliminate it here

        J1=J2=J3=zero;
        // J1 = <ij|kl>   
        // J1 < sqrt( <ik|ki> * <lj|jl> )  
        if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(l,j)) ) > cutoff1bar ) {
          J1 = H(i,j,k,l);
#if defined(QMC_COMPLEX)
          if(i==k && j==l) J1 = ValueType(J1.real(),0);
#endif
        }

        // J2 = <ij|lk> 
        if(i==j || l==k) {
          J2=J1;
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(k,j)) ) > cutoff1bar ) {
            J2 = H(i,j,l,k);
#if defined(QMC_COMPLEX)
            if(i==l && j==k) J2 = ValueType(J2.real(),0);
#endif
          }
        }

        // J3 = <ik|jl> 
        if(j==k) {
          J3=J1;
        } else if(i==l) {
          J3 = myconj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(l,k)) ) > cutoff1bar ) {
            J3 = H(i,k,j,l);
#if defined(QMC_COMPLEX)
            if(i==j && k==l) J3 = ValueType(J3.real(),0);
#endif
          }
        }

#if defined(QMC_COMPLEX)
        J1a=J2a=J3a=zero;

        //  J2a = <ik|lj>
        if(l==j) {
          J2a=J3;
        } else if(i==k) {
          J2a=J3;
        } else if(k==j) {
          J2a=J2;
        } else if(i==l) {
          J2a=std::conj(J2);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(j,k)) ) > cutoff1bar )
            J2a = H(i,k,l,j);
        }

        //  J3a = <il|jk> 
        if(l==j) {
          J3a=J2;
        } else if(i==k) {
          J3a=std::conj(J2);
        } else if(k==l) {
          J3a=J3;
        } else if(i==j) {
          J3a=std::conj(J3);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(k,l)) ) > cutoff1bar )
            J3a = H(i,l,j,k);
        }

        //  For complex, there are 3 extra non-symmetric terms:
        //  J1a = <il|kj>
        if(k==l) {
          J1a=J2a;
        } else if(i==j) {
          J1a=std::conj(J2a);
        } else if(j==k) {
          J1a=J3a;
        } else if(i==l) {
          J1a=J3a;
        } else if(l==j) {
          J1a=J1;
        } else if(i==k) {
          J1a=std::conj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(j,l)) ) > cutoff1bar )
            J1a = H(i,l,k,j);
        }
#endif



        if( std::abs(J1)<cutoff1bar && std::abs(J2)<cutoff1bar && std::abs(J3)<cutoff1bar && std::abs(J1a)<cutoff1bar && std::abs(J2a)<cutoff1bar && std::abs(J3a)<cutoff1bar) continue;

        vs4D.clear();
        if(walker_type==0) {
          find_all_contributions_to_hamiltonian_closed_shell(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==1) {
          find_all_contributions_to_hamiltonian_collinear(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==2) {
          find_all_contributions_to_hamiltonian_ghf(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else {
          APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
        }
        cnt2+=count_allowed_terms(vs4D,occ_a,occ_b);

    }
    }
    }
    }

    cnt2 = TG.Global().all_reduce_value(cnt2);
    //number_of_terms = cnt2;
    app_log()<<" Number of terms in hamiltonian: " <<cnt2 <<std::endl;
    Vijkl.reserve(cnt2);
    cnt2=0;

    cnter=0;
    // Approximation: (similar to non factorized case)
    //   - if <ij|kl>  <  cutoff1bar, set to zero 
    for(IndexType i=0; i<NMO; i++)  {
    for(IndexType j=i; j<NMO; j++,cnter++)  {
     if( cnter%npr != rk ) continue;
     occi = (occ_a[i]||occ_b[i])?1:0;
     occj = (occ_a[j]||occ_b[j])?1:0;
    for(IndexType k=j; k<NMO; k++)  {
     occk = (occ_a[k]||occ_b[k])?1:0;
    for(IndexType l=k; l<NMO; l++)  {

        occl = (occ_a[l]||occ_b[l])?1:0;
        if( occi + occj + occk + occl < 2 ) continue;

        J1=J2=J3=zero;
        // J1 = <ij|kl>   
        // J1 < sqrt( <ik|ki> * <lj|jl> )  
        if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(l,j)) ) > cutoff1bar ) {
          J1 = H(i,j,k,l);
#if defined(QMC_COMPLEX)
          if(i==k && j==l) J1 = ValueType(J1.real(),0);
#endif
        }

        // J2 = <ij|lk> 
        if(i==j || l==k) {
          J2=J1;
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(k,j)) ) > cutoff1bar ) {
            J2 = H(i,j,l,k);
#if defined(QMC_COMPLEX)
            if(i==l && j==k) J2 = ValueType(J2.real(),0);
#endif
          }
        }

        // J3 = <ik|jl> 
        if(j==k) {
          J3=J1;
        } else if(i==l) {
          J3 = myconj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(l,k)) ) > cutoff1bar ) {
            J3 = H(i,k,j,l);
#if defined(QMC_COMPLEX)
            if(i==j && k==l) J3 = ValueType(J3.real(),0);
#endif
          }
        }

#if defined(QMC_COMPLEX)
        J1a=J2a=J3a=zero;

        //  J2a = <ik|lj>
        if(l==j) {
          J2a=J3;
        } else if(i==k) {
          J2a=J3;
        } else if(k==j) {
          J2a=J2;
        } else if(i==l) {
          J2a=std::conj(J2);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(j,k)) ) > cutoff1bar )
            J2a = H(i,k,l,j);
        }

        //  J3a = <il|jk> 
        if(l==j) {
          J3a=J2;
        } else if(i==k) {
          J3a=std::conj(J2);
        } else if(k==l) {
          J3a=J3;
        } else if(i==j) {
          J3a=std::conj(J3);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(k,l)) ) > cutoff1bar )
            J3a = H(i,l,j,k);
        }

        //  For complex, there are 3 extra non-symmetric terms:
        //  J1a = <il|kj>
        if(k==l) {
          J1a=J2a;
        } else if(i==j) {
          J1a=std::conj(J2a);
        } else if(j==k) {
          J1a=J3a;
        } else if(i==l) {
          J1a=J3a;
        } else if(l==j) {
          J1a=J1;
        } else if(i==k) {
          J1a=std::conj(J1);
        } else {
          if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(j,l)) ) > cutoff1bar )
            J1a = H(i,l,k,j);
        }

#endif

        if( std::abs(J1)<cutoff1bar && std::abs(J2)<cutoff1bar && std::abs(J3)<cutoff1bar && std::abs(J1a)<cutoff1bar && std::abs(J2a)<cutoff1bar && std::abs(J3a)<cutoff1bar) continue;

        vs4D.clear();
        if(walker_type==0) {
          find_all_contributions_to_hamiltonian_closed_shell(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==1) {
          find_all_contributions_to_hamiltonian_collinear(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else if(walker_type==2) {
          find_all_contributions_to_hamiltonian_ghf(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
        } else {
          APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
        }
        cnt2+=add_allowed_terms(NMO,vs4D,occ_a,occ_b, Vijkl, true, walker_type==2);                   
    }
    }
    }
    }
    TG.global_barrier();
    Timer.stop("Generic");
    app_log()<<"Time to generate full Hamiltonian from factorized form: " <<Timer.total("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");
    redistribute_sparse_matrix(TGham,Vijkl);
    Timer.stop("Generic");
    app_log()<<"Time to communicate Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");

    if(!Vijkl.remove_repeated_and_compress(TG.Node().impl_)) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }

    Timer.stop("Generic");
    app_log()<<"Time to remove_repeated_and_compress Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    return true;

  }
*/
}
}
