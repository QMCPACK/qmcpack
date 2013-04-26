//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "QMCWaveFunctions/Fermion/DiracDeterminantTruncation.h"
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

extern "C" {
#include "ILUGMRESInterface.h"
}

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
DiracDeterminantTruncation::DiracDeterminantTruncation(SPOSetBasePtr const &spos, int first):
  DiracDeterminantBase(spos,first)
{
}

///default destructor
DiracDeterminantTruncation::~DiracDeterminantTruncation() {}

DiracDeterminantTruncation& DiracDeterminantTruncation::operator=(const DiracDeterminantTruncation& s)
{
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}


void
DiracDeterminantTruncation::set_truncation(int first, int nel,double &temp_cutoff,double &temp_radius)
{
  cutoff=temp_cutoff;
  FirstIndex = first;
  resize(nel,nel);
  radius=temp_radius;
}


DiracDeterminantTruncation::DiracDeterminantTruncation(const DiracDeterminantTruncation& s):
  DiracDeterminantBase(s)
{
  this->resize(s.NumPtcls,s.NumOrbitals);
}

void
DiracDeterminantTruncation::ChooseNearbyParticles(int ptcl,list<int> &closePtcls)
{
  closePtcls.clear();
  closePtcls.push_back(ptcl);
  int i=ptcl;
  for (int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
  {
    int j=d_table->J[nn];
    assert(j<NumPtcls);
    PosType disp=d_table->r(nn);
    if (std::sqrt(dot(disp,disp))<radius)
      closePtcls.push_back(j);
  }
  for (int ptcl=0; ptcl<i; ptcl++)
  {
    int nn=d_table->M[ptcl]+(i-ptcl)-1;
    int j=d_table->J[nn];
    assert(j<NumPtcls);
    PosType disp=d_table->r(nn);
    if (std::sqrt(dot(disp,disp))<radius)
      closePtcls.push_back(ptcl);
  }
}


//   void
//   DiracDeterminantTruncation::set_iterative(int first, int nel,double &temp_cutoff) {
//     cutoff=temp_cutoff;
//     FirstIndex = first;
//     resize(nel,nel);

//   }


void DiracDeterminantTruncation::resize(int nel, int morb)
{
  psiM_actual.resize(nel,morb);
  temp_psiM2.resize(nel,morb);
  psiM2.resize(nel,morb);
  psi_diff.resize(nel);
  particleLists.resize(nel);
  DiracDeterminantBase::resize(nel,morb);
}


//   void DiracDeterminantTruncation::SparseToCSR(vector<int> &Arp, vector<int> &Ari,vector<double> &Arx)
//   {
//     int systemSize=LastIndex-FirstIndex;
//     int nnz_index=0;
//     Arp.push_back(0);
//     cerr<<"Particle list size is "<<particleLists.size()<<endl;
//     for (int ptcl=0;ptcl<particleLists.size();ptcl++){
//       for (list<pair<int,double> >::iterator orb=particleLists[ptcl].begin();orb!=particleLists[ptcl].end();orb++){
// 	pair<int,double> myPair=*orb;
// 	int orbitalIndex=myPair.first;
// 	double value=myPair.second;
// 	//	if (abs(myValue)>=cutoff){
// 	  Ari.push_back(orbitalIndex);
// 	  Arx.push_back(value);
// 	  nnz_index++;
// 	}
//       //      }
//       Arp.push_back(nnz_index);
//     }
//     cerr<<"Ari size is "<<Ari.size()<<endl;

//   }

void DiracDeterminantTruncation::UpdatePsiM2(ValueVector_t &vec,int ptcl)
{
  for (int i=0; i<vec.size(); i++)
  {
    psi_diff(i)=vec(i)-psiM_actual(ptcl,i);
  }
  temp_psiM2=psiM2;
  MatrixOperators::half_outerProduct(psiM_actual,psi_diff,ptcl,temp_psiM2);
  MatrixOperators::other_half_outerProduct(psiM_actual,psi_diff,ptcl,temp_psiM2);
  temp_psiM2(ptcl,ptcl)+=MatrixOperators::innerProduct(psi_diff,psi_diff);
}

DiracDeterminantBase::ValueType DiracDeterminantTruncation::ratio(ParticleSet& P, int iat)
{
  //    cerr<<"Using local ratio "<<TestMe()<<endl;
  UpdateMode=ORB_PBYP_RATIO;
  WorkingIndex = iat-FirstIndex;
  assert(FirstIndex==0);
  list<int> nearbyPtcls;
  ChooseNearbyParticles(iat,nearbyPtcls);
  int num_nearbyPtcls=nearbyPtcls.size();
  ValueMatrix_t psiM2_small_old(num_nearbyPtcls,num_nearbyPtcls);
  ValueMatrix_t psiM2_small_new(num_nearbyPtcls,num_nearbyPtcls);
  int i=-1;
  for (list<int>::iterator  iIter=nearbyPtcls.begin(); iIter!=nearbyPtcls.end(); iIter++)
  {
    i++;
    int j=-1;
    for (list<int>::iterator  jIter=nearbyPtcls.begin(); jIter!=nearbyPtcls.end(); jIter++)
    {
      j++;
      //	cerr<<i<<" "<<j<<" "<<*iIter<<" "<<*jIter<<" "<<psiM2_small_old.extent(0)<<endl;
      assert(i<psiM2_small_old.extent(0));
      assert(j<psiM2_small_old.extent(1));
      psiM2_small_old(i,j)=psiM2(*iIter,*jIter);
      //	psiM2_small_old(i,j)=psiM2(i,j);
    }
  }
  Phi->evaluate(P, iat, psiV);
  UpdatePsiM2(psiV,iat);
  i=-1;
  for (list<int>::iterator  iIter=nearbyPtcls.begin(); iIter!=nearbyPtcls.end(); iIter++)
  {
    i++;
    int j=-1;
    for (list<int>::iterator  jIter=nearbyPtcls.begin(); jIter!=nearbyPtcls.end(); jIter++)
    {
      j++;
      psiM2_small_new(i,j)=temp_psiM2(*iIter,*jIter);
      //	psiM2_small_new(i,j)=temp_psiM2(i,j);
    }
  }
  int sign;
  //    cerr<<"Size "<<nearbyPtcls.size()<<endl;
  ValueType new_det=invert_matrix_log(psiM2_small_new,sign,true);
  ValueType old_det=invert_matrix_log(psiM2_small_old,sign,true);
  //    ValueType test_det=invert_matrix_log(psiM_actual,sign,true);
  //    invert_matrix_log(psiM_actual,sign,true);
  //    cerr<<"Preparing stuff "<<nearbyPtcls.size()<<" "<<new_det<<" "<<old_det<<" "<<test_det<<endl;
  //    vector<int> Arp; vector<int> Ari; vector<double> Arx;
  //    SparseToCSR(Arp,Ari,Arx);
  //    vector<int> Arp2(Arp.size()); vector<int> Ari2(Ari.size()); vector<double> Arx2(Arx.size());
  //    Arp2=Arp; Ari2=Ari; Arx2=Arx;
  //    int particleMoved=iat;
  //    int systemSize=LastIndex-FirstIndex;
  //    int nnzUpdatedPassed=Ari.size();
  //    vector<double> uPassed(psiV.size());
  //    double detRatio_ILU=0;
  //    for (int i=0;i<uPassed.size();i++)
  //      uPassed[i]=psiV[i];
  //    cerr<<"Calling stuff"<<systemSize<<" "<<Arp.size()<<endl;
  //    assert(systemSize+1==Arp.size());
  //    assert(Ari.size()<=nnzUpdatedPassed);
  //    cerr<<"Entering"<<endl;
  //    calcDeterminantILUGMRES(&particleMoved, &systemSize, &nnzUpdatedPassed, uPassed.data(), Arp.data(), Ari.data(), Arx.data(), Arp2.data(), Ari2.data(), Arx2.data(), &detRatio_ILU);
  //    int *Arp_ptr; int *Ari_ptr; double *Arx_ptr;
  //    DenseToCSR(psiM_actual,Arp_ptr,Ari_ptr,Arx_ptr);
  //    cerr<<"The size of my particle list is "<<particleLists[iat].size()<<" "<<cutoff<<endl;
  oldPtcl.clear();
  assert(iat<particleLists.size());
  particleLists[iat].swap(oldPtcl);
  for (int i=0; i<psiV.size(); i++)
  {
    if (abs(psiV(i))>=cutoff)
    {
      pair <int,double> temp(i,psiV(i));
      particleLists[iat].push_back(temp);
    }
  }
  curRatio=std::exp((new_det-old_det)/2);
  return curRatio;
#ifdef DIRAC_USE_BLAS
  curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
  //    cerr<<"RATIOS: "<<curRatio<<" "<<std::exp((new_det-old_det)/2)<<endl;
  return curRatio;
#else
  curRatio = DetRatio(psiM, psiV.begin(),iat-FirstIndex);
  //    cerr<<"RATIOS: "<<curRatio<<" "<<std::exp((new_det-old_det)/2)<<endl;
  return curRatio;
#endif
}



DiracDeterminantBase::ValueType DiracDeterminantTruncation::ratio(ParticleSet& P, int iat,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  //    cerr<<"doing large update"<<endl;
  UpdateMode=ORB_PBYP_ALL;
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  WorkingIndex = iat-FirstIndex;
#ifdef DIRAC_USE_BLAS
  curRatio = BLAS::dot(NumOrbitals,psiM_temp[WorkingIndex],&psiV[0]);
#else
  curRatio= DetRatio(psiM_temp, psiV.begin(),WorkingIndex);
#endif
  if(abs(curRatio)<numeric_limits<RealType>::epsilon())
  {
    UpdateMode=ORB_PBYP_RATIO; //singularity! do not update inverse
    return 0.0;
  }
  //update psiM_temp with the row substituted
  DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);
  //update dpsiM_temp and d2psiM_temp
  for(int j=0; j<NumOrbitals; j++)
  {
    dpsiM_temp(WorkingIndex,j)=dpsiV[j];
    d2psiM_temp(WorkingIndex,j)=d2psiV[j];
  }
  int kat=FirstIndex;
  const ValueType* restrict yptr=psiM_temp.data();
  const ValueType* restrict d2yptr=d2psiM_temp.data();
  const GradType* restrict dyptr=dpsiM_temp.data();
  for(int i=0; i<NumPtcls; i++,kat++)
  {
    //This mimics gemm with loop optimization
    GradType rv;
    ValueType lap=0.0;
    for(int j=0; j<NumOrbitals; j++,yptr++)
    {
      rv += *yptr * *dyptr++;
      lap += *yptr * *d2yptr++;
    }
    //using inline dot functions
    //GradType rv=dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
    //ValueType lap=dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);
    //Old index: This is not pretty
    //GradType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
    //ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
    //for(int j=1; j<NumOrbitals; j++) {
    //  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
    //  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
    //}
    lap -= dot(rv,rv);
    dG[kat] += rv - myG[kat];
    myG_temp[kat]=rv;
    dL[kat] += lap -myL[kat];
    myL_temp[kat]=lap;
  }
  return curRatio;
}



/** move was accepted, update the real container
 */
void DiracDeterminantTruncation::acceptMove(ParticleSet& P, int iat)
{
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  switch(UpdateMode)
  {
  case ORB_PBYP_RATIO:
    std::copy(psiV.begin(),psiV.end(),psiM_actual[WorkingIndex]);
    //	DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
    psiM2=temp_psiM2;
    break;
  case ORB_PBYP_PARTIAL:
    assert(2==3);
    std::copy(psiV.begin(),psiV.end(),psiM_actual[WorkingIndex]);
    DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    psiM2=temp_psiM2;
    //////////////////////////////////////
    ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
    //myG = myG_temp;
    //myL = myL_temp;
    ///////////////////////
    break;
  default:
    assert(1==2);
    myG = myG_temp;
    myL = myL_temp;
    psiM = psiM_temp;
    std::copy(psiV.begin(),psiV.end(),psiM_actual[WorkingIndex]);
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    psiM2=temp_psiM2;
    break;
  }
  curRatio=1.0;
}



///SHOULD THIS REALLY NOT RETURN ANYTHING!
DiracDeterminantBase::RealType
DiracDeterminantTruncation::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  return DiracDeterminantBase::evaluateLog(P,buf);
}




DiracDeterminantBase::RealType
DiracDeterminantTruncation::evaluateLog(ParticleSet& P,
                                        ParticleSet::ParticleGradient_t& G,
                                        ParticleSet::ParticleLaplacian_t& L)
{
  d_table=DistanceTable::add(P);
  Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
  ///I think at this point psiM has particles as the first index
  for (int ptcl=0; ptcl<psiM.extent(1); ptcl++)
  {
    particleLists[ptcl].clear();
    for (int orbital=0; orbital<psiM.extent(0); orbital++)
    {
      if (abs(psiM(orbital,ptcl))>=cutoff)
      {
        pair<int,double> temp(orbital,psiM(orbital,ptcl));
        particleLists[ptcl].push_back(temp);
      }
    }
  }
  if(NumPtcls==1)
  {
    //CurrentDet=psiM(0,0);
    ValueType det=psiM(0,0);
    ValueType y=1.0/det;
    psiM(0,0)=y;
    GradType rv = y*dpsiM(0,0);
    G(FirstIndex) += rv;
    L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
    LogValue = evaluateLogAndPhase(det,PhaseValue);
  }
  else
  {
    psiM_actual=psiM;
    MatrixOperators::transpose(psiM_actual);
    MatrixOperators::ABt(psiM_actual,psiM_actual,psiM2);
    LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
    int sign;
    ValueType new_det=invert_matrix_log(psiM2,sign,true);
    invert_matrix_log(psiM2,sign,true);
    ValueType old_det=invert_matrix_log(temp_psiM2,sign,true);
    ValueType test_det=invert_matrix_log(psiM_actual,sign,true);
    invert_matrix_log(psiM_actual,sign,true);
    //	cerr<<"init stuff "<<" "<<new_det<<" "<<old_det<<" "<<test_det<<endl;
    const ValueType* restrict yptr=psiM.data();
    const ValueType* restrict d2yptr=d2psiM.data();
    const GradType* restrict dyptr=dpsiM.data();
    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      GradType rv;
      ValueType lap=0.0;
      for(int j=0; j<NumOrbitals; j++,yptr++)
      {
        rv += *yptr * *dyptr++;
        lap += *yptr * *d2yptr++;
      }
      G(iat) += rv;
      L(iat) += lap - dot(rv,rv);
    }
  }
  return LogValue;
}



}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3265 $   $Date: 2008-10-15 09:20:33 -0500 (Wed, 15 Oct 2008) $
 * $Id: DiracDeterminantTruncation.cpp 3265 2008-10-15 14:20:33Z jnkim $
 ***************************************************************************/
