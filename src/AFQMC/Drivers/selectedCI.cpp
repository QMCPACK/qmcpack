#include<tuple>
#include<map>
#include<string>
#include<iomanip>

#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include<Message/MPIObjectBase.h>
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsData/libxmldefs.h"
#include "Configuration.h"
#include <qmc_common.h>

#include "AFQMC/config.h"
#include "AFQMC/Drivers/selectedCI.h"

#include "AFQMC/Utilities/Utils.h"

namespace qmcplusplus {

bool selectedCI::run()
{

  app_log()<<" Running selectedCI. \n";

  int ne = NAEA+NAEB;
  std::vector<IndexType> intm;
  std::vector<IndexType> new_dets;
  std::vector<ValueType> new_ci;
  intm.reserve(ne*100000);
  std::vector<RealType> eigVal(1);
  std::vector<IndexType> vira,virb;
  std::vector<IndexType> indx(ne);
  std::vector<OrbitalType> KLs;
  vira.reserve(NMO-NAEA);
  virb.reserve(NMO-NAEA);
  ValueMatrix eigVec;

  // generate excitation tables and hamiltonian
  sHam->generate_selCI_Ham(cutoff_list);

  occ_orbs.clear();
  ci.clear();
  occ_orbs.reserve(ne*100000);
  ci.reserve(100000);
  new_dets.reserve(ne*100000);
  new_ci.reserve(100000);
  // occ_orbs: determinants in the permanent list that have already been excited from 
  // new_dets: determinants generated in the previous iteration, from which excitations are generated
  //           in the current iteration
  for(int i=0; i<NAEA; i++)
    new_dets.push_back(i);
  for(int i=NMO; i<NMO+NAEB; i++)
    new_dets.push_back(i);
  new_ci.push_back(1.0);

// ideas for speed-ups
// 1. Use the fact that H is stored in sparse ordered form to
//    only consider terms that have hmat*ci > cut   
//    To do this, make a routine taht given (i,j) it gets all
//    (k,l) such that h(i,j,k,l)*ci > cut     
// 2. Keep a different list for those determinants added in the last
//    step (those new to the list that haven't been "excited" before)
//    This is a reasonable approximation since the ci coeffs do not change
//    much. Then only excite from this list.  
// 3. For parallelization, we need to solve both hamiltonian storage (which will be bigger now)
//    and fare share of work. For now, assume hamiltonian is  
// 4. use (real) parallel? sparse eigenvalue solver that produces only a given # of states 

  for(int ite=0; ite<maxit; ite++) {
    app_log()<<" Iteration: " <<ite <<std::endl;
    intm.clear();
    int nci = new_dets.size()/ne;
    int nterms=0;
    Timer.reset("Generic");
    Timer.start("Generic");
    std::vector<IndexType>::iterator it=new_dets.begin();
    if((occ_orbs.size()+new_dets.size()) > occ_orbs.capacity()) occ_orbs.reserve( occ_orbs.size()+new_dets.size() + 10000*ne );
    if((ci.size()+new_ci.size()) > ci.capacity()) ci.reserve( ci.size()+new_ci.size() + 10000 );

    for(int nc = 0; nc<nci; nc++, it+=ne) {
      vira.clear();
      virb.clear();
      for(int iv=0; iv<NMO; iv++)
        if(!binary_search (it, it+NAEA, iv)) vira.push_back(iv);
      for(int iv=NMO; iv<2*NMO; iv++)
        if(!binary_search (it+NAEA, it+ne, iv)) virb.push_back(iv);
      if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
      // add new_dets[nc] to occ_orbs
      ci.push_back(new_ci[nc]);
      for(int ii=0; ii<ne; ii++) 
        occ_orbs.push_back(*(it+ii));

      double cut = cutoff_list/std::abs(new_ci[nc]);
      // double excitations 
      for(int ia=0; ia<NAEA; ia++) {
        OrbitalType a = *(it+ia);
        // aa
        for(int ib=ia+1; ib<NAEA; ib++) {
          OrbitalType b = *(it+ib);

          sHam->get_selCI_excitations(a,b,0,cut,&(*it),KLs);
          if((nterms+KLs.size())*ne >= intm.capacity()) intm.reserve(ne*(nterms+KLs.size()+10000));
          for(std::vector<OrbitalType>::iterator itkl=KLs.begin(); itkl<KLs.end(); itkl+=2) { 
            for(int ii=0; ii<ne; ii++) indx[ii] = *(it+ii); 
            indx[ia] = *(itkl);
            indx[ib] = *(itkl+1);
            std::sort(indx.begin(),indx.end());
            intm.insert(intm.end(),std::begin(indx),std::end(indx)); 
            nterms++;
          }
        }

        // ab
        for(int ib=NAEA; ib<ne; ib++) {
          OrbitalType b = *(it+ib);
          sHam->get_selCI_excitations(a,b,1,cut,&(*it),KLs);
          if((nterms+KLs.size())*ne >= intm.capacity()) intm.reserve(ne*(nterms+KLs.size()+10000));
          for(std::vector<OrbitalType>::iterator itkl=KLs.begin(); itkl<KLs.end(); itkl+=2) {
            for(int ii=0; ii<ne; ii++) indx[ii] = *(it+ii);
            indx[ia] = *(itkl);
            indx[ib] = *(itkl+1);
            std::sort(indx.begin(),indx.end());
            intm.insert(intm.end(),std::begin(indx),std::end(indx));
            nterms++;
          }
        }

      } 

      for(int ia=NAEA; ia<ne; ia++) {
        int a = *(it+ia);
        for(int ib=ia+1; ib<ne; ib++) {
          int b = *(it+ib);

          sHam->get_selCI_excitations(a,b,3,cut,&(*it),KLs);
          if((nterms+KLs.size())*ne >= intm.capacity()) intm.reserve(ne*(nterms+KLs.size()+10000));
          for(std::vector<OrbitalType>::iterator itkl=KLs.begin(); itkl<KLs.end(); itkl+=2) {
            for(int ii=0; ii<ne; ii++) indx[ii] = *(it+ii);
            indx[ia] = *(itkl);
            indx[ib] = *(itkl+1);
            std::sort(indx.begin(),indx.end());
            intm.insert(intm.end(),std::begin(indx),std::end(indx));
            nterms++;
          }

        }
      }
    }  // states in new_dets 
    Timer.stop("Generic");
    app_log()<<"Time to generate excitations: " <<Timer.total("Generic") <<std::endl;

    if(occ_orbs.size() == 0) {
      app_error()<<" Error in selectedCI::run(): Main determinant list is empty. \n";
      APP_ABORT(" Error in selectedCI::run(): Main determinant list is empty. \n"); 
    }
    if(ci.size() == 0) {
      app_error()<<" Error in selectedCI::run(): Main ci determinant list is empty. \n";
      APP_ABORT(" Error in selectedCI::run(): Main ci determinant list is empty. \n"); 
    }
    if(ci.size() != occ_orbs.size()/ne) {
      app_error()<<" Error in selectedCI::run(): Main determinant list size is inconsistent, ci, dets: " <<ci.size() <<" " <<occ_orbs.size() <<" \n";
      APP_ABORT(" Error in selectedCI::run(): Main determinant list size is inconsistent. \n"); 
    }

    Timer.reset("Generic");
    Timer.start("Generic");
    // sort occ_orbs/ci
    new_dets.clear();
    new_ci.clear();   
    if(ci.size() > 1)
      sort_multiple(occ_orbs.begin(),occ_orbs.end()-ne,ci.begin(),ci.end()-1);

/*
//debug
int ntr = occ_orbs.size()/ne;
for(int i=0; i<ntr; i++)
 for(int k=i+1; k<ntr; k++)
  if(list_order(occ_orbs.begin()+k*ne,occ_orbs.begin()+i*ne)) {
    app_error()<<" Error with order of occ_orbs. \n";
    return false;
  }
*/
    // clean intm list
    app_log()<<" Intermediate list has " <<nterms <<" new terms (before cleanup)" <<std::endl;
    if(intm.size() > 0) { 
      val_at_pivot.resize(NAEA+NAEB);
      sort_list(intm.begin(),intm.end()-ne);
/*
//debug
ntr = intm.size()/ne;
for(int i=0; i<ntr; i++)
 for(int k=i+1; k<ntr; k++)
  if(list_order(intm.begin()+k*ne,intm.begin()+i*ne)) {
    app_error()<<" Error with order of intm. \n";
    return false;
  }
*/
      remove_repeated(intm,occ_orbs);  
/*
//debug
int ntr1 = occ_orbs.size()/ne;
ntr = intm.size()/ne;
for(int i=0; i<ntr1; i++)
 for(int k=0; k<ntr; k++)
  if(list_equal(occ_orbs.begin()+i*ne,intm.begin()+k*ne)) {
    app_error()<<" Error: Repeated elements after call to remove_repeated. \n"; 
    return false;
  }
*/
      if(intm.size()%ne != 0) APP_ABORT("Error: After remove_repeated. \n\n\n");
      nterms = intm.size()/ne;
    } else {
      nterms=0;
    }
    
    Timer.stop("Generic");
    app_log()<<"Time to sort and clean new list of dets: " <<Timer.total("Generic") <<std::endl;
  
    app_log()<<" Intermediate list has " <<nterms <<" new terms" <<std::endl;

    if(nterms == 0) {
      app_log()<<" Intermediate determinant list is empty. Stopping iterations. \n";
      break; 
    }
    bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,occ_orbs,occ_orbs.size()/ne,intm,intm.size()/ne,true);
    if(!sucess) {
      app_error()<<" Error: Problems with diagonalizeTrialWavefunction. \n";
      return false;
    }
    app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic3") <<std::endl;
    app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic4") <<std::endl;
    int cnt=0;
    nci = ci.size();
    RealType normlz = 0;
    for(int ii=0; ii<nci; ii++) 
      normlz += mynorm(eigVec(0,ii));
    for(int ii=0; ii<nterms; ii++) 
      if(std::abs(eigVec(0,ii+nci)) > cutoff_diag) { 
        cnt++;
        normlz += mynorm(eigVec(0,ii+nci));
      }
    normlz = std::sqrt(normlz); 
    app_log()<<" Normalization of ci vector: " <<normlz <<std::endl;
    new_dets.reserve(cnt);
    new_ci.reserve(cnt);
    for(int ii=0; ii<nci; ii++) 
      ci[ii] = eigVec(0,ii)/normlz;
    for(int ii=0; ii<nterms; ii++) {
      if(std::abs(eigVec(0,ii+nci)) > cutoff_diag) {
        new_ci.push_back(eigVec(0,ii+nci)/normlz);
        for(int j=0; j<ne; j++)
          new_dets.push_back(intm[ii*ne+j]);
      }   
    }

    std::ofstream out("iterativeCI.dat");
    if(out.fail()) {
      app_error()<<" Problems opening iterativeCI.dat \n";
      return false;
    }
    {
      std::vector<std::tuple<double,int> > dets;
      dets.reserve(ci.size()+new_ci.size());
      for(int i=0; i<ci.size(); i++) dets.push_back(std::make_tuple( std::abs(ci[i]),i+1));
      for(int i=0; i<new_ci.size(); i++) dets.push_back(std::make_tuple( std::abs(new_ci[i]),-(i+1)));
      std::sort( dets.begin(), dets.end(),
      [] (const std::tuple<double,int>& a, const std::tuple<double,int>& b)
               {return (std::get<0>(a)>std::get<0>(b));} );
      out<<" &FCI \n NCI = " <<dets.size() <<" \n /\n";
      for(int i=0; i<dets.size(); i++) {
        if(std::get<1>(dets[i]) > 0) {
          out<<ci[std::get<1>(dets[i])-1] <<" ";
          int nt = (std::get<1>(dets[i])-1)*ne;
          for(int j=0; j<ne; j++) out<<occ_orbs[nt+j]+1 <<" ";
          out<<"\n";
        } else {
          int nt = -std::get<1>(dets[i])-1; 
          out<<new_ci[nt] <<" ";
          for(int j=0; j<ne; j++) out<<new_dets[nt*ne+j]+1 <<" ";
          out<<"\n";
        }
      }
      out.close();
    }

    app_log()<<" Energy: " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;
    app_log()<<" Number of determinants after truncation: " <<(occ_orbs.size()+new_dets.size())/ne <<std::endl;

  } // iteration 

  if(diag_in_steps>0) {
    if((occ_orbs.size()+new_dets.size()) > occ_orbs.capacity()) occ_orbs.reserve( occ_orbs.size()+new_dets.size() );
    if((ci.size()+new_ci.size()) > ci.capacity()) ci.reserve( ci.size()+new_ci.size() );
    for(int nc=0; nc<new_ci.size(); nc++) {
      ci.push_back(new_ci[nc]);
      for(int ii=0; ii<ne; ii++)
        occ_orbs.push_back(*(new_dets.begin()+nc*ne+ii));
    }
    new_dets.clear();
    new_ci.clear();
    sort_multiple(occ_orbs.begin(),occ_orbs.end()-ne,ci.begin(),ci.end()-1);

    app_log()<<"\n***********************************************\n  #Determinants        Energy: " <<"\n";
    std::vector<std::tuple<double,int> > dets(ci.size());
    for(int i=0; i<ci.size(); i++) dets[i] = std::make_tuple( std::abs(ci[i]),i);
    std::sort( dets.begin(), dets.end(),
      [] (const std::tuple<double,int>& a, const std::tuple<double,int>& b)
               {return (std::get<0>(a)>std::get<0>(b));} );
    for(int i=1; i<ci.size(); i+=diag_in_steps) {
      intm.clear();
      for(int ki=0; ki<i; ki++) {
       int kk = std::get<1>(dets[ki]);
       for(int kj=0; kj<ne; kj++) intm.push_back(occ_orbs[ kk*ne+kj]);
      }
      bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,intm,i,new_dets,0,false);
      if(!sucess) {
        app_error()<<" Error: Problems with diagonalizeTrialWavefunction. \n";
        return false;
      }
      app_log()<<i <<" " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;
    }
    app_log()<<"***********************************************" <<std::endl <<std::endl;
  }

  return true; 
}

void selectedCI::sort_list(std::vector<IndexType>::iterator left, std::vector<IndexType>::iterator right)
{
  std::vector<IndexType>::iterator i = left, j = right;
  std::vector<IndexType>::iterator pivot = left; 
  std::advance(pivot, (std::distance(left,right)/(NAEA+NAEB))/2*(NAEA+NAEB)); 
  std::copy(pivot,pivot+NAEA+NAEB,val_at_pivot.begin()); 
  pivot = val_at_pivot.begin();

  /* partition */
  while (i <= j) {
    while (list_order(i,pivot))
      i+=(NAEA+NAEB);
    while (list_order(pivot,j))
      j-=(NAEA+NAEB);
    if(i <= j) {
      for(int k=0; k<NAEA+NAEB; k++) std::swap( *(i+k), *(j+k) );  
      i+=(NAEA+NAEB);
      j-=(NAEA+NAEB);
    }
  };

  /* recursion */
  if (left < j)
    sort_list(left, j);
  if (i < right)
    sort_list(i, right);
}

void selectedCI::sort_multiple(std::vector<IndexType>::iterator left, std::vector<IndexType>::iterator right, std::vector<ValueType>::iterator vl, std::vector<ValueType>::iterator vr)
{
  if(left==right) return;
  std::vector<IndexType>::iterator i = left, j = right;
  std::vector<ValueType>::iterator vi = vl, vj = vr;
  std::vector<IndexType>::iterator pivot = left; 
  std::advance(pivot, (std::distance(left,right)/(NAEA+NAEB))/2*(NAEA+NAEB)); 
  std::copy(pivot,pivot+NAEA+NAEB,val_at_pivot.begin()); 
  pivot = val_at_pivot.begin();

  /* partition */
  while (i <= j) {
    while (list_order(i,pivot)) {
      i+=(NAEA+NAEB);
      vi++;
    }
    while (list_order(pivot,j)) {
      j-=(NAEA+NAEB);
      vj--;
    }
    if(i <= j) {
      for(int k=0; k<NAEA+NAEB; k++) std::swap( *(i+k), *(j+k) );  
      std::swap(*vi,*vj);
      i+=(NAEA+NAEB);
      j-=(NAEA+NAEB);
      vi++;
      vj--;
    }
  };

  /* recursion */
  if (left < j)
    sort_multiple(left, j, vl, vj);
  if (i < right)
    sort_multiple(i, right, vi, vr);
}

// 1. remove repeated from vnew
// 2. remove from vnew intersection with vold
void selectedCI::remove_repeated(std::vector<IndexType>& vnew, std::vector<IndexType>& vold)
{

  std::vector<IndexType>::iterator first = vnew.begin();
  std::vector<IndexType>::iterator last = vnew.end();
  std::vector<IndexType>::iterator new_last = last; 
  
  if (first!=last) { 
    new_last=first;  
    first+=(NAEA+NAEB);
    while (first < last)
    {
      if (!list_equal(new_last,first)) { 
        new_last+=(NAEA+NAEB); 
        for(int i=0; i<NAEA+NAEB; i++) *(new_last+i) = *(first++);
      } else
        first+=(NAEA+NAEB); 
    }
    new_last+=(NAEA+NAEB); 
  }
  // cut part of the vector with repeated elements
  vnew.resize( std::distance(vnew.begin(),new_last) );

  // now new_last points to the end of the good segment of the list 
  // remove intersection with vold
  if(vnew.size() == 0 || vold.size()==0) return; 

   
  std::vector<IndexType>::iterator vold_first = vold.begin();
  first = vnew.begin();
  last = vnew.end(); 
  new_last = first;

  // loop through vnew, when a good (not found in vold) element is found, copy to new_last 
  while (first < last)
  { 
    if (!mysearch(vold_first,vold.end(),first)) {
      for(int i=0; i<NAEA+NAEB; i++) *(new_last++) = *(first++);
    } else
      first+=(NAEA+NAEB);
  }

  // cut part of the vector with repeated elements
  vnew.resize( std::distance(vnew.begin(),new_last) );

}

// right now I'm building the Hamiltonian from scratch
// Later on store it and and rotate it according to changes in the ordering in occ_orbs
bool selectedCI::diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ValueMatrix& eigVec, std::vector<IndexType>& occ1, int nci1, std::vector<IndexType>& occ2, int nci2, bool eigV )
{
  ValueType one = ValueType(1.0);
  ValueType zero = ValueType(0.0);
  bool sucess;

    for(int i=0; i<nci1; i++)
      std::sort(occ1.begin()+i*(NAEA+NAEB),occ1.begin()+(i+1)*(NAEA+NAEB));
    for(int i=0; i<nci2; i++)
      std::sort(occ2.begin()+i*(NAEA+NAEB),occ2.begin()+(i+1)*(NAEA+NAEB));

    int nci = nci1+nci2;

    if(myComm->rank()==0) {

      Timer.reset("Generic3");
      Timer.start("Generic3");
      ValueMatrix hm(nci,nci);
      ValueType let;
      RealType sg;
      std::vector<IndexType> occ(NAEA+NAEB);
      IndexType n0,n1,n2,n3;
      std::vector<IndexType> DL(NAEA+NAEB);
      std::vector<IndexType> DR(NAEA+NAEB);

      for(int ki=0; ki<nci; ki++) {
        // i==j
        let=zero;
        std::vector<IndexType>::iterator it;
        if(ki < nci1)
          it = occ1.begin()+ki*(NAEA+NAEB);
        else
          it = occ2.begin()+(ki-nci1)*(NAEA+NAEB);
        for(int i=0; i<NAEA+NAEB; i++)
        {
          let += sHam->H(*(it+i),*(it+i));
          for(int j=i+1; j<NAEA+NAEB; j++) {
            let += sHam->H(*(it+i),*(it+j),*(it+i),*(it+j)) - sHam->H(*(it+i),*(it+j),*(it+j),*(it+i));
          }
        }

        hm(ki,ki) = let;
        for(int kj=ki+1; kj<nci; kj++) {

          std::vector<IndexType>::iterator jt;
          if(kj < nci1)
            jt = occ1.begin()+kj*(NAEA+NAEB);
          else
            jt = occ2.begin()+(kj-nci1)*(NAEA+NAEB);
          std::copy(it,it+NAEA+NAEB,DL.begin());
          std::copy(jt,jt+NAEA+NAEB,DR.begin());
          int cnt = cntExcitations(NAEA,NAEB,DL,DR,n0,n1,n2,n3,occ,sg);

          if(cnt==0) {
            app_error()<<" Error: Found repeated determinant in trial wave function in MultiPureSingleDeterminant \n";
            return false;
          } else if(cnt==2) {
            int nterms = NAEA+NAEB-1;
            let=sHam->H(n0,n1);
            for(int i=0; i<nterms; i++)
              let+=sHam->H(n0,occ[i],n1,occ[i]) - sHam->H(n0,occ[i],occ[i],n1);
            hm(ki,kj)=let*sg;
          } else if(cnt==4) {
            hm(ki,kj) = sg*(sHam->H(n0,n1,n2,n3) - sHam->H(n0,n1,n3,n2));
          } else {
            hm(ki,kj) =  zero;
          }

          hm(kj,ki) = myconj(hm(ki,kj));
        }
      }
      Timer.stop("Generic3");
      //app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic2") <<std::endl;

      Timer.reset("Generic4");
      Timer.start("Generic4");

      eigVal.resize(1);
      if(eigV) eigVec.resize(1,nci);
      sucess = DenseMatrixOperators::symEigenSysSelect(nci,hm.data(),nci,1,eigVal.data(),eigV,eigVec.data(),eigVec.size2());

      Timer.stop("Generic4");
      //app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic2") <<std::endl;

    } else {
      eigVal.resize(1);
      if(eigV) eigVec.resize(1,nci);
    }

    myComm->bcast(sucess);
    myComm->bcast(eigVal.data(),eigVal.size(),0,myComm->getMPI());
    if(eigV) myComm->bcast(eigVec.data(),eigVec.size1()*eigVec.size2(),0,myComm->getMPI());

    return sucess;

}

bool selectedCI::parse(xmlNodePtr cur)
{
  if(cur==NULL) return false;

  if(cur == NULL)
    return false;

  xmlNodePtr curRoot=cur;

  maxit=0;
  cutoff_list=cutoff_diag=0;
  build_full_hamiltonian = true;

  std::string str("yes");
  ParameterSet m_param;
  m_param.add(str,"full_hamiltonian","std::string");
  m_param.add(output_filename,"output_filename","std::string");
  m_param.add(output_filename,"output","std::string");
  m_param.add(cutoff_list,"cutoff_list","double");
  m_param.add(cutoff_diag,"cutoff_diag","double");
  m_param.add(maxit,"maxit","int");
  m_param.add(diag_in_steps,"diag_steps","int");
  m_param.put(cur);

  std::transform(str.begin(),str.end(),str.begin(),(int (*)(int)) tolower);
  if(str == "no" || str == "false") build_full_hamiltonian = false; 

  return true;
}


bool selectedCI::setup(HamPtr h0, WSetPtr w0, PropPtr p0, WfnPtr wf0)
{

  if(h0==NULL) {
    app_error()<<" Error: Null Hamiltonian pointer in selectedCI::setup(). \n";
    return false; 
  }

  ham0=h0;
  wlkBucket=NULL;
  prop0=NULL;
  wfn0=NULL;

  sHam = dynamic_cast<SparseGeneralHamiltonian*>(ham0);
  if(!sHam) {
    app_error()<<" Error in MultiPureSingleDeterminant::getHamiltonian. \n"
               <<" Hamiltonian associated with MultiPureSingleDeterminant must of the \n"
               <<" type SparseGeneralHamiltonian. \n";
    APP_ABORT("");
  }

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Beginning Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  app_log()<<" Using " <<ncores_per_TG <<" cores per node in a TaskGroup. \n";
  // right now this TG is not used. It is needed for setup purposes and to  
  // get a unique TG number for every group of cores on a node (used in the WalkerSet)
  TG.setup(ncores_per_TG,1,false);  
  std::vector<int> TGdata(5);
  TG.getSetupInfo(TGdata); 

  // setup local-to-node MPI Comm
  // TGdata[0]: node_number
  myComm->split_comm(TGdata[0],MPI_COMM_NODE_LOCAL);
  TG.setNodeCommLocal(MPI_COMM_NODE_LOCAL);
  int key = TG.getTGNumber();
  myComm->split_comm(key,MPI_COMM_TG_LOCAL);
  TG.setTGCommLocal(MPI_COMM_TG_LOCAL);

  key = TG.getCoreID();
  myComm->split_comm(key,MPI_COMM_HEAD_OF_NODES);
  TG.setHeadOfNodesComm(MPI_COMM_HEAD_OF_NODES);

  CommBuffer.setup(TG.getCoreRank()==0,std::string("COMMBuffer_")+std::to_string(myComm->rank()),MPI_COMM_TG_LOCAL);
  TG.setBuffer(&CommBuffer);

  app_log()<<"\n****************************************************\n"
           <<"               Initializating Hamiltonian \n"
           <<"****************************************************\n"
           <<std::endl;

  // hamiltonian
  if(!ham0->init(TGdata,&CommBuffer,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,MPI_COMM_HEAD_OF_NODES)) {
    app_error()<<"Error initializing Hamiltonian in selectedCI::setup" <<std::endl; 
    return false; 
  }   

  NuclearCoulombEnergy = static_cast<ValueType>(sHam->NuclearCoulombEnergy); 

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Finished Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  return true;
}

// writes checkpoint file
bool selectedCI::checkpoint(int block, int step) 
{
  return true;
}

// sets up restart archive and reads  
bool selectedCI::restart(hdf_archive& read)
{
  return true;
}

bool selectedCI::clear()
{
  return true;
}

}

