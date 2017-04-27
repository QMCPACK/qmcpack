#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_H
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_H

#include<fstream>

#include "AFQMC/config.h"
#include <Message/MPIObjectBase.h>
#include "Numerics/DeterminantOperators.h"
#include "Numerics/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include "Message/Communicate.h"
//#include "Message/CommOperators.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Walkers/SlaterDetWalker.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

namespace qmcplusplus
{

//class SlaterDetOperations: public EstimatorBase 
class SlaterDetOperations: public MPIObjectBase, public AFQMCInfo 
{
  public:

    typedef HamiltonianBase* HamPtr;

    SlaterDetOperations(Communicate *c):MPIObjectBase(c), ham(NULL) {}

    void setup(HamPtr h, myTimer* timer_) {
      ham=h;
      timer=timer_;
      GF.resize(2*NMO,NMO);
      tGF.resize(2*NMO,NMO);
      V0.resize(2*NMO*NMO);
      Cwork.resize(2*NMO);
      pivot.resize(2*NMO);
    }

    void green_function(ComplexType* A, ComplexType* B, ComplexType& ovlp, SPComplexMatrix& G, bool getG=true) {
      const ComplexType one = ComplexType(1.0);
      const ComplexType zero = ComplexType(0.0); 

      // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) )
      ComplexMatrix S0(NAEA,NAEA);
      ComplexMatrix S1(NAEB,NAEB);
      ComplexMatrix SS0(2*NMO,NAEA);
      

      if(getG) tGF = ComplexType(0.0); 
      S0 = ComplexType(0.0); 
      S1 = ComplexType(0.0); 

      // S0 = transpose(conjg(A))*B  
      DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,A,NAEA,B,NAEA,zero,S0.data(),NAEA);

      // S0 = S0^-1
      ovlp = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

      // SS0 = B * S0
      if(getG) DenseMatrixOperators::product(NMO,NAEA,NAEA,one,B,NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);   
      // G(beta) = SS0*transpose(conjg(A))  
      if(getG) DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,A,NAEA,zero,tGF.data(),NMO);

      // S1 = transpose(conjg(A))*B  
      DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,A+NMO*NAEA,NAEA,B+NAEA*NMO,NAEA,zero,S1.data(),NAEB);

      // S0 = S0^-1
      ovlp *= Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

      if(!getG) return;

      if(std::abs(ovlp) < 1e-6) {
        G = SPComplexType(0.0);
        return;
      }

      // SS0(beta) = B(beta) * S1
      DenseMatrixOperators::product(NMO,NAEB,NAEB,one,B+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

      // G(beta) = SS0*transpose(conjg(A))  
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,A+NAEA*NMO,NAEA,zero,tGF.data()+NMO*NMO,NMO);

      for(int i=0; i<NMO; i++)
       for(int j=0; j<i; j++) {
         std::swap(tGF(i,j),tGF(j,i));
         std::swap(tGF(i+NMO,j),tGF(j+NMO,i));
       }
      std::copy(tGF.begin(),tGF.end(),G.begin());

    }
    
    void green_function(ComplexMatrix& A, ComplexMatrix& B, ComplexType& ovlp, ComplexMatrix& G, bool getG=true) {
      const ComplexType one = ComplexType(1.0);
      const ComplexType zero = ComplexType(0.0); 

      // G = transpose( B * ( transpose(conjg(A)) * B )^-1 * transpose(conjg(A)) )
      ComplexMatrix S0(NAEA,NAEA);
      ComplexMatrix S1(NAEB,NAEB);
      ComplexMatrix SS0(2*NMO,NAEA);

      if(getG) G = ComplexType(0.0); 
      S0 = ComplexType(0.0); 
      S1 = ComplexType(0.0); 

      // S0 = transpose(conjg(A))*B  
      DenseMatrixOperators::product_AhB(NAEA,NAEA,NMO,one,A.data(),NAEA,B.data(),NAEA,zero,S0.data(),NAEA);

      // S0 = S0^-1
      ovlp = Invert(S0.data(), NAEA, NAEA, Cwork.data(),pivot.data());

      // SS0 = B * S0
      if(getG) DenseMatrixOperators::product(NMO,NAEA,NAEA,one,B.data(),NAEA,S0.data(),NAEA,zero,SS0.data(),NAEA);   
      // G(beta) = SS0*transpose(conjg(A))  
      if(getG) DenseMatrixOperators::product_ABh(NMO,NMO,NAEA,one,SS0.data(),NAEA,A.data(),NAEA,zero,G.data(),NMO);

      // S1 = transpose(conjg(A))*B  
      DenseMatrixOperators::product_AhB(NAEB,NAEB,NMO,one,A.data()+NMO*NAEA,NAEA,B.data()+NAEA*NMO,NAEA,zero,S1.data(),NAEB);

      // S0 = S0^-1
      ovlp *= Invert(S1.data(), NAEB, NAEB, Cwork.data(),pivot.data());

      if(!getG) return;

      if(std::abs(ovlp) < 1e-6) {
        G = ComplexType(0.0);
        return;
      }

      // SS0(beta) = B(beta) * S1
      DenseMatrixOperators::product(NMO,NAEB,NAEB,one,B.data()+NAEA*NMO,NAEA,S1.data(),NAEB,zero,SS0.data()+NAEA*NMO,NAEA);

      // G(beta) = SS0*transpose(conjg(A))  
      DenseMatrixOperators::product_ABh(NMO,NMO,NAEB,one,SS0.data()+NAEA*NMO,NAEA,A.data()+NAEA*NMO,NAEA,zero,G.data()+NMO*NMO,NMO);

      for(int i=0; i<NMO; i++)
       for(int j=0; j<i; j++) {
         std::swap(G(i,j),G(j,i));
         std::swap(G(i+NMO,j),G(j+NMO,i));
       }

    }

    void matrix_element_and_overlap(ComplexType* A, ComplexType* B, ComplexType& ovlp, ComplexType& hamME) {
      
      green_function(A,B,ovlp,GF);

      if( std::abs(ovlp) < 1e-6 ) {
        ovlp = ComplexType(0.0);
        hamME = ComplexType(0.0);
        return;
      }

      SPValueSMSpMat *V;
      std::vector<s1D<ValueType> > *h;
      int nr1=1, nc1=2*NMO*NMO;
      SPValueType one = SPValueType(1.0);
      SPValueType zero = SPValueType(0.0);

      ham->getFullHam(h,V);   

      hamME = ComplexType(0.0); 

      SparseMatrixOperators::product_SpMatV(nc1,nc1,one,V->values(), V->column_data(),  V->row_index() ,GF.data(),zero,V0.data());

      SPComplexMatrix::iterator itG = GF.begin();
      SPComplexVector::iterator itV = V0.begin();
      for(int i=0; i<nc1; i++,++itG,++itV) hamME += static_cast<ComplexType>(*itV) * static_cast<ComplexType>(*itG);
      hamME = 0.5*hamME+ham->NuclearCoulombEnergy;    

      std::vector<s1D<ValueType> >::iterator end1 = h->end();
      itG = GF.begin();
      for(std::vector<s1D<ValueType> >::iterator it = h->begin(); it != end1; it++)
        hamME += static_cast<ComplexType>(*(itG + std::get<0>(*it))) * std::get<1>(*it);

      hamME *= ovlp;

    }

    void diag( std::vector<SlaterDetWalker>::iterator itbegin, std::vector<SlaterDetWalker>::iterator itend, int nstates, std::vector<RealType>& eigVal, ComplexMatrix& eigVec, ComplexType& exactEnergy, ComplexMatrix* HF, bool getEigV=false ) {

      if(myComm->size() > 1 ) 
        APP_ABORT(" ERROR: Estimators::SlaterDetOperations::diag(): Only implemented in serial. \n");

      int N = 0;
      for( std::vector<SlaterDetWalker>::iterator it1 = itbegin; it1!=itend; it1++)
        if(it1->alive && std::abs(it1->weight) > 1e-3) N++;
      if(N == 0) return;
      if(HF!=NULL) N++;
      nstates = std::min(nstates,N);
      ComplexMatrix H(N),S(N);
      exactEnergy = ComplexType(0.0);
      ComplexType nume=ComplexType(0.0);
      ComplexType deno=ComplexType(0.0);
#ifdef AFQMC_TIMER
      timer->start("SlaterDetOperations::diag::evaluate_H_S");
#endif
      int i=0;
      if(HF!=NULL){ // always include HF state in list
        int j=i;
        matrix_element_and_overlap(HF->data(),HF->data(),S(i,j),H(i,j));
        j++; 
        for( std::vector<SlaterDetWalker>::iterator it2 = itbegin; it2!=itend; it2++) {
          if( !(it2->alive && std::abs(it2->weight) > 1e-3) ) continue;
          matrix_element_and_overlap(HF->data(),(it2->SlaterMat).data(),S(i,j),H(i,j));
          if(i!=j) {
            H(j,i)=std::conj(H(i,j));
            S(j,i)=std::conj(S(i,j));
          }
          j++;
        }
        i++;
      }
      for( std::vector<SlaterDetWalker>::iterator it1 = itbegin; it1!=itend; it1++) {
        if( !(it1->alive && std::abs(it1->weight) > 1e-3) ) continue;
        int j=i;
        for( std::vector<SlaterDetWalker>::iterator it2 = it1; it2!=itend; it2++) {
          if( !(it2->alive && std::abs(it2->weight) > 1e-3) ) continue;
          matrix_element_and_overlap((it1->SlaterMat).data(),(it2->SlaterMat).data(),S(i,j),H(i,j));
          nume += it1->weight*it2->weight*H(i,j); 
          deno += it1->weight*it2->weight*S(i,j); 
          if(i!=j) {
            H(j,i)=std::conj(H(i,j));
            S(j,i)=std::conj(S(i,j));
            nume += std::conj(it1->weight)*it2->weight*H(j,i)/(std::conj(std::get<0>(it1->overlap_alpha)*std::get<0>(it1->overlap_beta)) * std::get<0>(it2->overlap_alpha)*std::get<0>(it2->overlap_beta) ); 
            deno += std::conj(it1->weight)*it2->weight*S(j,i)/(std::conj(std::get<0>(it1->overlap_alpha)*std::get<0>(it1->overlap_beta)) * std::get<0>(it2->overlap_alpha)*std::get<0>(it2->overlap_beta) ); 
          }        
          j++; 
        }
        i++;
      }
      exactEnergy = nume/deno;
#ifdef AFQMC_TIMER
      timer->stop("SlaterDetOperations::diag::evaluate_H_S");
#endif
#ifdef AFQMC_TIMER
      timer->start("SlaterDetOperations::diag::solve_GEV");
#endif
      if(nstates > 0) {
        std::vector<int> ifail(N);
        eigVal.resize(nstates);
        getEigV=true;
        if(getEigV)
          eigVec.resize(nstates,N);   
/*
std::ofstream out1("Hr.dat");
std::ofstream out2("Hc.dat");
std::ofstream out3("Sr.dat");
std::ofstream out4("Sc.dat");
for(int i=0; i<N; i++) { 
for(int j=0; j<N; j++) out1<<H(i,j).real() <<" ";  
for(int j=0; j<N; j++) out2<<H(i,j).imag() <<" ";  
out1<<std::endl;
out2<<std::endl;
}
for(int i=0; i<N; i++) { 
for(int j=0; j<N; j++) out3<<S(i,j).real() <<" ";  
for(int j=0; j<N; j++) out4<<S(i,j).imag() <<" ";  
out3<<std::endl;
out4<<std::endl;
}
out1.close();
out2.close();
out3.close();
out4.close();
APP_ABORT("Testing. \n");
*/
        bool sucess = DenseMatrixOperators::genHermitianEigenSysSelect(N,H.data(),N,S.data(),N,nstates,eigVal.data(),getEigV,eigVec.data(),eigVec.size2(),ifail.data());
        if(!sucess) for(int i=0; i<nstates; i++) eigVal[i]=0.0;
        else {
          std::ofstream out("diag.dat",std::ios_base::app | std::ios_base::out);
          std::vector<double> coeff(N);
          for(int i=0; i<N; i++) coeff[i] = std::abs(eigVec(1,i));
          std::sort(coeff.begin(),coeff.end()); 
          for(int i=0; i<N; i++) out<<coeff[i] <<" ";
          out<<std::endl;
          out.close();
        } 
      }
#ifdef AFQMC_TIMER
      timer->stop("SlaterDetOperations::diag::solve_GEV");
#endif

    }

    ComplexType overlap( std::vector<SlaterDetWalker>::iterator itbegin, std::vector<SlaterDetWalker>::iterator itend ) {

      std::vector<ComplexType> ovlp(2,ComplexType(0.0));
      ComplexType sum_w = ComplexType(0.0);
      ComplexMatrix A(2*NMO,NAEA);
      ComplexMatrix B(2*NMO,NAEA);
      ComplexMatrix G(1);
      std::vector<char> buffer_in;
      std::vector<char> buffer_out;
      std::vector<int> to(1);
      std::vector<int> from(myComm->size());
      int sz = itbegin->sizeForDump();
      int nWtot=0, nW = 0, nWmax=0; 
      for(std::vector<SlaterDetWalker>::iterator it1 = itbegin; it1!=itend; it1++)
        if(it1->alive) nW++;
  
      to[0]=nW;
      myComm->allgather(to,from,1);
      for(int i=0; i<myComm->size(); i++) {
        nWtot += from[i];
        if(from[i] > nWmax) nWmax=from[i];
      }

      buffer_out.resize(nW*sz);

      int cnt=0;
      for(std::vector<SlaterDetWalker>::iterator it=itbegin; it!=itend; it++)
        if(it->alive) {
          it->dumpToChar( buffer_out.data()+cnt );
          cnt+=sz;
        }
      
      ovlp[0] = ovlp[1] = ComplexType(0.0);
      ComplexType w1, w2, o1, o2, e1, e2, ov;
      // diagonal contribution
      for(int i=0; i<nW; i++) {
        itbegin->unpackFromChar(buffer_out.data()+sz*i,A,w1,e1,o1);
        ovlp[0] += w1;    
        for(int j=i; j<nW; j++) {
          itbegin->unpackFromChar(buffer_out.data()+j*sz,B,w2,e2,o2);
          green_function(A,B,ov,G,false);
          ovlp[1] += std::conj(w1)*w2*ov/(std::conj(o1)*o2)
                  +  std::conj(w2)*w1*std::conj(ov)/(std::conj(o2)*o1);
        }
      }

      if(myComm->size() == 1)
        return ovlp[0]/std::sqrt(std::abs(ovlp[1])); 
      buffer_in.resize(nWmax*sz);
      int rec = (myComm->rank()+1)%(myComm->size()); 
      int send = (myComm->rank()-1)%(myComm->size()); 
      for(int i=0; i<myComm->size()-1; i++) {

//        myComm->isend(send, send*myComm->size()+myComm->rank() ,buffer_out);
//        myComm->irecv(rec, myComm->rank()*myComm->size()+rec ,buffer_in);

        // dump way to avoid double counting, but efficiency depends heavily on load balance
        if( rec < myComm->rank() ) {
          // I only do the top half
          for(int i=0; i<from[rec]/2; i++) {
            itbegin->unpackFromChar(buffer_in.data()+sz*i,A,w1,e1,o1);
            for(int j=0; j<nW; j++) {
              itbegin->unpackFromChar(buffer_out.data()+j*sz,B,w2,e2,o2);
              green_function(A,B,ov,G,false);
              ovlp[1] += std::conj(w1)*w2*ov/(std::conj(o1)*o2)
                      +  std::conj(w2)*w1*std::conj(ov)/(std::conj(o2)*o1);
            }
          }    
        } else {
          // I only do the bottom half
          for(int i=nW/2; i<nW; i++) {
            itbegin->unpackFromChar(buffer_out.data()+sz*i,A,w1,e1,o1);
            for(int j=0; j<from[rec]; j++) {
              itbegin->unpackFromChar(buffer_in.data()+j*sz,B,w2,e2,o2);
              green_function(A,B,ov,G,false);
              ovlp[1] += std::conj(w1)*w2*ov/(std::conj(o1)*o2)
                      +  std::conj(w2)*w1*std::conj(ov)/(std::conj(o2)*o1);
            }
          }
        } 

        rec = (rec+1)%(myComm->size()); 
        send = (send-1)%(myComm->size()); 
      }
  
      std::vector<ComplexType> res(2);
      //myComm->gsum(ovlp); 
      return res[0]/std::sqrt(std::abs(res[1]));   
    }

  private:

    std::vector<ComplexType> Cwork;
    std::vector<int> pivot;

    HamPtr ham; 
    ComplexMatrix tGF;
    SPComplexMatrix GF;
    SPComplexVector V0;
    myTimer* timer;
};

}

#endif
