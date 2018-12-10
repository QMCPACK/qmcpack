//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_LRBREAKUP_H
#define QMCPLUSPLUS_LRBREAKUP_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "LongRange/KContainer.h"
#include "Numerics/OhmmsBlas.h"
#include <cassert>

namespace qmcplusplus
{

template<class BreakupBasis>
struct LRBreakup
{

  DECLARE_COULOMB_TYPES

  //Typedef for the lattice-type. We don't need the full particle-set.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;

  //We use an internal k-list with degeneracies to do the breakup.
  //We do this because the number of vectors is much larger than we'd
  //use elsewhere.
  void AddKToList(mRealType k, mRealType degeneracy=1.0);

  ///The basis to be used for breakup.
  BreakupBasis& Basis;
  /// For each k,  KList[k][0] = |k| and KList[k][1] = degeneracy
  std::vector<TinyVector<mRealType,2> > KList;
  /** setup KList
   * @param kc k-space cutoff for long-range sums
   * @param kcont  k at which approximate (spherical shell) degeneracies are used.
   * @param kmax largest k used for performing the breakup
   * @return the maximum kshell for the given kc
   */
  int SetupKVecs(mRealType kc, mRealType kcont, mRealType kmax);

  //Fk is FT of F_full(r) up to kmax
  //adjust is used for constraining values in the breakup

  /* REPLACED SO WE CAN USE TYPES OTHER THAN STL VECTOR.
      mRealType DoBreakup(const std::vector<mRealType> &Fk, std::vector<mRealType> &t,
  		       const std::vector<bool> &adjust);
      mRealType DoBreakup(const std::vector<mRealType> &Fk, std::vector<mRealType> &t);
  */
  mRealType DoBreakup(mRealType *Fk, mRealType *t, mRealType *adjust);
  mRealType DoGradBreakup(mRealType *Fk, mRealType *t, mRealType *adjust);
  mRealType DoStrainBreakup(mRealType *Fk, mRealType *dFk, mRealType *t, mRealType *adjust);
  void DoAllBreakup(mRealType *chisqr, mRealType *Fk, mRealType *dFk, mRealType *t, 
                    mRealType *gt, mRealType *dt, mRealType *adjust); 
  mRealType DoBreakup(mRealType *Fk, mRealType *t)
  {
    const mRealType tolerance = std::numeric_limits<mRealType>::epsilon();
    //t must be allocated up to Basis.NumBasisElem();
    //Fk must be allocated and filled up to KList.size();
    //  assert(t.size()==Basis.NumBasisElem());
    Matrix<mRealType> A;
    std::vector<mRealType> b;
    Matrix<mRealType> cnk;
    int numElem = Basis.NumBasisElem(); //t.size();
    A.resize(numElem, numElem);
    b.resize(numElem,0.0);
    cnk.resize(numElem,KList.size());
    // Fill in cnk.
    // app_log() << "Check OMP size : numElem, KList.size : " << numElem << " , " << KList.size() << std::endl;
    #pragma omp parallel for shared(cnk)
    for (int n=0; n<numElem; n++)
    {
      for (int ki=0; ki<KList.size(); ki++)
      {
        mRealType k = KList[ki][0];
        cnk(n,ki) = Basis.c(n,k);
      }
    }
    // Now, fill in A and b
    A = 0.0;
    for (int l=0; l<numElem; l++)
    {
      for (int ki=0; ki<KList.size(); ki++)
      {
        b[l] += KList[ki][1]*Fk[ki] * cnk(l, ki);
        for (int n=0; n<numElem; n++)
          A(l,n) += KList[ki][1]*cnk(l,ki)*cnk(n,ki);
      }
    }
    //////////////////////////
    //Do the SVD:
    // Matrix<mRealType> U(numElem, numElem), V(numElem, numElem);
    // std::vector<mRealType> S(numElem), Sinv(numElem);
    //////////////////////////
    //  SVdecomp(A, U, S, V);
    //////////////////////////
    int M = A.rows();
    int N = A.cols();
    Matrix<mRealType> Atrans(N,M);
    Matrix<mRealType> U, V;
    U.resize(std::min(M,N),M);
    V.resize(N,std::min(M,N));
    std::vector<mRealType> S, Sinv;
    S.resize(std::min(N,M));
    //Do the transpose
    for(int i=0; i<M; i++)
    {
      for(int j=0; j<N; j++)
        Atrans(j,i) = A(i,j);
    }
    char JOBU = 'S';
    char JOBVT = 'S';
    int LDA = M;
    int LDU = M;
    int LDVT = std::min(M,N);
    int LWORK = 10*std::max(3*std::min(N,M)+std::max(M,N),5*std::min(M,N));
    std::vector<mRealType> WORK(LWORK);
    int INFO;
    LAPACK::gesvd(&JOBU,&JOBVT,&M,&N,Atrans.data(),&LDA,&S[0],U.data(),
           &LDU,V.data(),&LDVT,&WORK[0],&LWORK,&INFO);
    assert(INFO==0);
    int ur = U.rows();
    int uc = U.cols();
    Matrix<mRealType> Utrans(uc,ur);
    for(int i=0; i<ur; i++)
    {
      for(int j=0; j<uc; j++)
        Utrans(j,i) = U(i,j);
    }
    U.resize(uc,ur);
    U=Utrans;
    ///////////////////////////////////
    // Zero out near-singular values
    mRealType Smax=S[0];
    for (int i=1; i<S.size(); i++)
      Smax = std::max(S[i],Smax);
    Sinv.resize(S.size());
    for (int i=0; i<S.size(); i++)
      Sinv[i] = (S[i] < (tolerance*Smax)) ? 0.0 : (1.0/S[i]);
    int numSingular = 0;
    for (int i=0; i<Sinv.size(); i++)
      if (Sinv[i] == 0.0)
        numSingular++;
    if (numSingular > 0)
      std::cout << "There were " << numSingular << " singular values in breakup.\n";
    for(int i=0; i<numElem; i++)
      t[i] = 0.0;
    // Compute t_n, removing singular values
    for (int i=0; i<numElem; i++)
    {
      mRealType coef = 0.0;
      for (int j=0; j<numElem; j++)
        coef += U(j,i) * b[j];
      coef *= Sinv[i];
      for (int k=0; k<numElem; k++)
        t[k] += coef * V(k,i);
    }
    // Calculate chi-squared
    mRealType Yk, chi2;
    chi2 = 0.0;
    for (int ki=0; ki<KList.size(); ki++)
    {
      Yk = Fk[ki];
      for (int n=0; n<numElem; n++)
      {
        Yk -= cnk(n,ki)*t[n];
      }
      chi2 += KList[ki][1]*Yk*Yk;
    }
    return (chi2);
  }


  //The constructor. Call the constructor of basis...
  //set up the basis parameters too.
  LRBreakup (BreakupBasis& bref) : Basis(bref)
  {
    /*Do Nothing*/
  }


  mRealType DoGradBreakup(mRealType *Fk, mRealType *t)
  {
    const mRealType tolerance = std::numeric_limits<mRealType>::epsilon();
    //t must be allocated up to Basis.NumBasisElem();
    //Fk must be allocated and filled up to KList.size();
    //  assert(t.size()==Basis.NumBasisElem());
    Matrix<mRealType> A;
    std::vector<mRealType> b;
    Matrix<mRealType> cnk;
    int numElem = Basis.NumBasisElem(); //t.size();
    A.resize(numElem, numElem);
    b.resize(numElem,0.0);
    cnk.resize(numElem,KList.size());
    // Fill in cnk.
    for (int n=0; n<numElem; n++)
    {
      for (int ki=0; ki<KList.size(); ki++)
      {
        mRealType k = KList[ki][0];
        cnk(n,ki) = Basis.c(n,k);
      }
    }
    // Now, fill in A and b
    
    A = 0.0;
    for (int l=0; l<numElem; l++)
    {
      for (int ki=0; ki<KList.size(); ki++)
      {
		mRealType k2=KList[ki][0]*KList[ki][0];
        b[l] += k2*KList[ki][1]*Fk[ki] * cnk(l, ki);
        for (int n=0; n<numElem; n++)
          A(l,n) += k2*KList[ki][1]*cnk(l,ki)*cnk(n,ki);
      }
    }
    //////////////////////////
    //Do the SVD:
    // Matrix<mRealType> U(numElem, numElem), V(numElem, numElem);
    // std::vector<mRealType> S(numElem), Sinv(numElem);
    //////////////////////////
    //  SVdecomp(A, U, S, V);
    //////////////////////////
    int M = A.rows();
    int N = A.cols();
    Matrix<mRealType> Atrans(N,M);
    Matrix<mRealType> U, V;
    U.resize(std::min(M,N),M);
    V.resize(N,std::min(M,N));
    std::vector<mRealType> S, Sinv;
    S.resize(std::min(N,M));
    //Do the transpose
    for(int i=0; i<M; i++)
    {
      for(int j=0; j<N; j++)
        Atrans(j,i) = A(i,j);
    }
    char JOBU = 'S';
    char JOBVT = 'S';
    int LDA = M;
    int LDU = M;
    int LDVT = std::min(M,N);
    int LWORK = 10*std::max(3*std::min(N,M)+std::max(M,N),5*std::min(M,N));
    std::vector<mRealType> WORK(LWORK);
    int INFO;
    LAPACK::gesvd(&JOBU,&JOBVT,&M,&N,Atrans.data(),&LDA,&S[0],U.data(),
           &LDU,V.data(),&LDVT,&WORK[0],&LWORK,&INFO);
    assert(INFO==0);
    int ur = U.rows();
    int uc = U.cols();
    Matrix<mRealType> Utrans(uc,ur);
    for(int i=0; i<ur; i++)
    {
      for(int j=0; j<uc; j++)
        Utrans(j,i) = U(i,j);
    }
    U.resize(uc,ur);
    U=Utrans;
    ///////////////////////////////////
    // Zero out near-singular values
    mRealType Smax=S[0];
    for (int i=1; i<S.size(); i++)
      Smax = std::max(S[i],Smax);
    Sinv.resize(S.size());
    for (int i=0; i<S.size(); i++)
      Sinv[i] = (S[i] < (tolerance*Smax)) ? 0.0 : (1.0/S[i]);
    int numSingular = 0;
    for (int i=0; i<Sinv.size(); i++)
      if (Sinv[i] == 0.0)
        numSingular++;
    if (numSingular > 0)
      std::cout << "There were " << numSingular << " singular values in breakup.\n";
    for(int i=0; i<numElem; i++)
      t[i] = 0.0;
    // Compute t_n, removing singular values
    for (int i=0; i<numElem; i++)
    {
      mRealType coef = 0.0;
      for (int j=0; j<numElem; j++)
        coef += U(j,i) * b[j];
      coef *= Sinv[i];
      for (int k=0; k<numElem; k++)
        t[k] += coef * V(k,i);
    }
    // Calculate chi-squared
    mRealType Yk, chi2;
    chi2 = 0.0;
    for (int ki=0; ki<KList.size(); ki++)
    {
	  mRealType k2=KList[ki][0]*KList[ki][0];
      Yk = Fk[ki];
      for (int n=0; n<numElem; n++)
      {
        Yk -= cnk(n,ki)*t[n];
      }
      chi2 += k2*KList[ki][1]*Yk*Yk;
    }
    return (chi2);
  }

};

template<class BreakupBasis>
void
LRBreakup<BreakupBasis>::AddKToList(mRealType k,
                                    mRealType degeneracy /* =1.0 */)
{
  //Search for this k already in list
  int ki=0;
  while((ki < KList.size()) && (std::abs(k-KList[ki][0]) > 1.0e-12))
    ki++;
  if(ki==KList.size())
  {
    TinyVector<mRealType,2> temp(k,degeneracy);
    KList.push_back(temp);
  }
  else
    KList[ki][1] += degeneracy;
}


template<class BreakupBasis>
int
LRBreakup<BreakupBasis>::SetupKVecs(mRealType kc, mRealType kcont, mRealType kmax)
{
  //Add low |k| ( < kcont) k-points with exact degeneracy
  KContainer kexact;
  kexact.UpdateKLists(Basis.get_Lattice(),kcont);
  bool findK=true;
  mRealType kc2=kc*kc;
  //use at least one shell
  size_t ks=0;
  kc2 = std::max(kc2,static_cast<mRealType>(kexact.ksq[kexact.kshell[ks]]));
  while(findK)
  {
    if(kexact.ksq[kexact.kshell[ks]]>kc2)
      findK=false;
    else
      ks++;
  }
  size_t maxkshell=ks;
  size_t numk =kexact.numk-kexact.kshell[ks];
  for(; ks<kexact.kshell.size()-1; ks++)
    AddKToList(std::sqrt(kexact.ksq[kexact.kshell[ks]]), kexact.kshell[ks+1]-kexact.kshell[ks]);
  ////Add these vectors to the internal list
  //int numk=0;
  //mRealType modk2;
  //for(int ki=0; ki<kexact.numk; ki++) {
  //  modk2 = dot(kexact.kpts_cart[ki],kexact.kpts_cart[ki]);
  //  if(modk2 > (kc*kc)) { //Breakup needs kc < k < kcont.
  //    AddKToList(std::sqrt(modk2));
  //    numk++;
  //  }
  //}
  //Add high |k| ( >kcont, <kmax) k-points with approximate degeneracy
  //Volume of 1 K-point is (2pi)^3/(a1.a2^a3)
#if OHMMS_DIM==3
  mRealType kelemvol = 8*M_PI*M_PI*M_PI/Basis.get_CellVolume();
  //Generate 4000 shells:
  const int N=4000;
  mRealType deltak = (kmax-kcont)/N;
  for(int i=0; i<N; i++)
  {
    mRealType k1 = kcont + deltak*i;
    mRealType k2 = k1 + deltak;
    mRealType kmid = 0.5*(k1+k2);
    mRealType shellvol = 4.0*M_PI*(k2*k2*k2-k1*k1*k1)/3.0;
    mRealType degeneracy = shellvol/kelemvol;
    AddKToList(kmid,degeneracy);
    numk += static_cast<int>(degeneracy);
  }
#elif OHMMS_DIM==2
  mRealType kelemvol = 4*M_PI*M_PI/Basis.get_CellVolume();
  //Generate 8000 shells:
  const int N=8000;
  mRealType deltak = (kmax-kcont)/N;
  for(int i=0; i<N; i++)
  {
    mRealType k1 = kcont + deltak*i;
    mRealType k2 = k1 + deltak;
    mRealType kmid = 0.5*(k1+k2);
    mRealType shellvol = M_PI*(k2*k2-k1*k1);
    mRealType degeneracy = shellvol/kelemvol;
    AddKToList(kmid,degeneracy);
    numk += static_cast<int>(degeneracy);
  }
#endif
  app_log()<<"  NUMBER OF OPT_BREAK KVECS = "<<numk<< std::endl;

  return maxkshell;
  //numk now contains the total number of vectors.
  //this->klist.size() contains the number of unique vectors.
}

//Do the constrained breakup
template<class BreakupBasis>
typename LRBreakup<BreakupBasis>::mRealType
LRBreakup<BreakupBasis>::DoBreakup(mRealType *Fk,
                                   mRealType *t,
                                   mRealType *adjust)
{
  const mRealType tolerance = std::numeric_limits<mRealType>::epsilon();
  //t and adjust must be allocated up to Basis.NumBasisElem();
  //Fk must be allocated and filled up to KList.size();
  //  assert(t.size()==adjust.size());
  //  assert(t.size()==Basis.NumBasisElem());
  Matrix<mRealType> A;
  std::vector<mRealType> b;
  Matrix<mRealType> cnk;
  int N = Basis.NumBasisElem(); //t.size();
  A.resize(N,N);
  b.resize(N,0.0);
  cnk.resize(N,KList.size());
  //Fill in cnk.
  for (int n=0; n<N; n++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k = KList[ki][0];
      cnk(n,ki) = Basis.c(n,k);
    }
  }
  //Fill in A and b
  A = 0.0;
  for (int l=0; l<N; l++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      b[l] += KList[ki][1]*Fk[ki] * cnk(l, ki);
      for (int n=0; n<N; n++)
        A(l,n) += KList[ki][1]*cnk(l,ki)*cnk(n,ki);
    }
  }
  //Reduce for constraints
  int M = N;
  for (int i=0; i<N; i++)
    if (!adjust[i])
      M--;
  //The c is for "constrained"
  Matrix<mRealType> Ac;
  Ac.resize(M,M);
  std::vector<mRealType> bc(M,0.0), tc(M,0.0);
  //Build constrained Ac and bc
  int j=0;
  for (int col=0; col<N; col++)
  {
    if (adjust[col])
    {
      // Copy column a A to Ac
      int i=0;
      for (int row=0; row<N; row++)
        if (adjust[row])
        {
          Ac(i,j) = A(row,col);
          i++;
        }
      j++;
    }
    else
    {
      // Otherwise, subtract t(col)*A(:,col) from bc
      for (int row=0; row<N; row++)
        b[row] -= A(row,col)*t[col];
    }
  }
  j=0;
  for (int row=0; row<N; row++)
    if (adjust[row])
    {
      bc[j] = b[row];
      j++;
    }
  // Do SVD:
  // -------
  // Matrix<mRealType> U(M, M), V(M, M);
  // std::vector<mRealType> S(M), Sinv(M);
  // SVdecomp(Ac, U, S, V);
  ////////////////////////////////
  int m = Ac.rows();
  int n = Ac.cols();
  Matrix<mRealType> Atrans(n,m);
  Matrix<mRealType> U, V;
  U.resize(std::min(m,n),m);
  V.resize(n,std::min(m,n));
  std::vector<mRealType> S, Sinv;
  S.resize(std::min(n,m));
  //do the transpose
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
      Atrans(j,i) = Ac(i,j);
  }
  char JOBU = 'S';
  char JOBVT = 'S';
  int LDA = m;
  int LDU = m;
  int LDVT = std::min(m,n);
  int LWORK = 10*std::max(3*std::min(n,m)+std::max(m,n),5*std::min(m,n));
  std::vector<mRealType> WORK(LWORK);
  int INFO;
  LAPACK::gesvd(&JOBU,&JOBVT,&m,&n,Atrans.data(),&LDA,&S[0],U.data(),
         &LDU,V.data(),&LDVT,&WORK[0],&LWORK,&INFO);
  assert(INFO==0);
  int ur = U.rows();
  int uc = U.cols();
  Matrix<mRealType> Utrans(uc,ur);
  for(int i=0; i<ur; i++)
  {
    for(int j=0; j<uc; j++)
      Utrans(j,i) = U(i,j);
  }
  U.resize(uc,ur);
  U=Utrans;
  //////////////////////////////////
  // Zero out near-singular values
  mRealType Smax=S[0];
  for (int i=1; i<M; i++)
    Smax = std::max(S[i],Smax);
  for (int i=0; i<M; i++)
    if (S[i] < 0.0)
      std::cout << "negative singlar value.\n";
  //  perr << "Smax = " << Smax << std::endl;
  Sinv.resize(S.size());
  for (int i=0; i<M; i++)
    Sinv[i] = (S[i] < (tolerance*Smax)) ? 0.0 : (1.0/S[i]);
  int numSingular = 0;
  for (int i=0; i<Sinv.size(); i++)
    if (Sinv[i] == 0.0)
      numSingular++;
  if (numSingular > 0)
    std::cout << "There were " << numSingular << " singular values in breakup.\n";
  // Compute t_n, removing singular values
  for (int i=0; i<M; i++)
  {
    mRealType coef = 0.0;
    for (int j=0; j<M; j++)
      coef += U(j,i) * bc[j];
    coef *= Sinv[i];
    for (int k=0; k<M; k++)
      tc[k] += coef * V(k,i);
  }
  // Now copy tc values into t
  j=0;
  for (int i=0; i<N; i++)
    if (adjust[i])
    {
      t[i] = tc[j];
      j++;
    }
  // Calculate chi-squared
  mRealType Yk, chi2;
  chi2 = 0.0;
  for (int ki=0; ki<KList.size(); ki++)
  {
    Yk = Fk[ki];
    for (int n=0; n<N; n++)
    {
      Yk -= cnk(n,ki)*t[n];
    }
    chi2 += KList[ki][1]*Yk*Yk;
  }
  return (chi2);
}


template<class BreakupBasis>
typename LRBreakup<BreakupBasis>::mRealType
LRBreakup<BreakupBasis>::DoGradBreakup(mRealType *Fk,
                                   mRealType *t,
                                   mRealType *adjust)
{
  const mRealType tolerance = std::numeric_limits<mRealType>::epsilon();
  //t and adjust must be allocated up to Basis.NumBasisElem();
  //Fk must be allocated and filled up to KList.size();
  //  assert(t.size()==adjust.size());
  //  assert(t.size()==Basis.NumBasisElem());
  Matrix<mRealType> A;
  std::vector<mRealType> b;
  Matrix<mRealType> cnk;
  int N = Basis.NumBasisElem(); //t.size();
  A.resize(N,N);
  b.resize(N,0.0);
  cnk.resize(N,KList.size());
  //Fill in cnk.
  for (int n=0; n<N; n++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k = KList[ki][0];
      cnk(n,ki) = Basis.c(n,k);
    }
  }
  //Fill in A and b
  A = 0.0;
  for (int l=0; l<N; l++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k2=KList[ki][0]*KList[ki][0];
      b[l] += k2*KList[ki][1]*Fk[ki] * cnk(l, ki);
      for (int n=0; n<N; n++)
        A(l,n) += k2*KList[ki][1]*cnk(l,ki)*cnk(n,ki);
    }
  }
  //Reduce for constraints
  int M = N;
  for (int i=0; i<N; i++)
    if (!adjust[i])
      M--;
  //The c is for "constrained"
  Matrix<mRealType> Ac;
  Ac.resize(M,M);
  std::vector<mRealType> bc(M,0.0), tc(M,0.0);
  //Build constrained Ac and bc
  int j=0;
  for (int col=0; col<N; col++)
  {
    if (adjust[col])
    {
      // Copy column a A to Ac
      int i=0;
      for (int row=0; row<N; row++)
        if (adjust[row])
        {
          Ac(i,j) = A(row,col);
          i++;
        }
      j++;
    }
    else
    {
      // Otherwise, subtract t(col)*A(:,col) from bc
      for (int row=0; row<N; row++)
        b[row] -= A(row,col)*t[col];
    }
  }
  j=0;
  for (int row=0; row<N; row++)
    if (adjust[row])
    {
      bc[j] = b[row];
      j++;
    }
  // Do SVD:
  // -------
  // Matrix<mRealType> U(M, M), V(M, M);
  // std::vector<mRealType> S(M), Sinv(M);
  // SVdecomp(Ac, U, S, V);
  ////////////////////////////////
  int m = Ac.rows();
  int n = Ac.cols();
  Matrix<mRealType> Atrans(n,m);
  Matrix<mRealType> U, V;
  U.resize(std::min(m,n),m);
  V.resize(n,std::min(m,n));
  std::vector<mRealType> S, Sinv;
  S.resize(std::min(n,m));
  //do the transpose
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
      Atrans(j,i) = Ac(i,j);
  }
  char JOBU = 'S';
  char JOBVT = 'S';
  int LDA = m;
  int LDU = m;
  int LDVT = std::min(m,n);
  int LWORK = 10*std::max(3*std::min(n,m)+std::max(m,n),5*std::min(m,n));
  std::vector<mRealType> WORK(LWORK);
  int INFO;
  LAPACK::gesvd(&JOBU,&JOBVT,&m,&n,Atrans.data(),&LDA,&S[0],U.data(),
         &LDU,V.data(),&LDVT,&WORK[0],&LWORK,&INFO);
  assert(INFO==0);
  int ur = U.rows();
  int uc = U.cols();
  Matrix<mRealType> Utrans(uc,ur);
  for(int i=0; i<ur; i++)
  {
    for(int j=0; j<uc; j++)
      Utrans(j,i) = U(i,j);
  }
  U.resize(uc,ur);
  U=Utrans;
  //////////////////////////////////
  // Zero out near-singular values
  mRealType Smax=S[0];
  for (int i=1; i<M; i++)
    Smax = std::max(S[i],Smax);
  for (int i=0; i<M; i++)
    if (S[i] < 0.0)
      std::cout << "negative singlar value.\n";
  //  perr << "Smax = " << Smax << std::endl;
  Sinv.resize(S.size());
  for (int i=0; i<M; i++)
    Sinv[i] = (S[i] < (tolerance*Smax)) ? 0.0 : (1.0/S[i]);
  int numSingular = 0;
  for (int i=0; i<Sinv.size(); i++)
    if (Sinv[i] == 0.0)
      numSingular++;
  if (numSingular > 0)
    std::cout << "There were " << numSingular << " singular values in breakup.\n";
  // Compute t_n, removing singular values
  for (int i=0; i<M; i++)
  {
    mRealType coef = 0.0;
    for (int j=0; j<M; j++)
      coef += U(j,i) * bc[j];
    coef *= Sinv[i];
    for (int k=0; k<M; k++)
      tc[k] += coef * V(k,i);
  }
  // Now copy tc values into t
  j=0;
  for (int i=0; i<N; i++)
    if (adjust[i])
    {
      t[i] = tc[j];
      j++;
    }
  // Calculate chi-squared
  mRealType Yk, chi2;
  chi2 = 0.0;
  for (int ki=0; ki<KList.size(); ki++)
  {
    Yk = Fk[ki];
    for (int n=0; n<N; n++)
    {
      Yk -= cnk(n,ki)*t[n];
    }
    chi2 += KList[ki][0]*KList[ki][0]*KList[ki][1]*Yk*Yk;
  }
  return (chi2);
}
template<class BreakupBasis>
typename LRBreakup<BreakupBasis>::mRealType
LRBreakup<BreakupBasis>::DoStrainBreakup(mRealType *Fk, mRealType *dFk,
                                   mRealType *t,
                                   mRealType *adjust)
{
  const mRealType tolerance = std::numeric_limits<mRealType>::epsilon();
  //t and adjust must be allocated up to Basis.NumBasisElem();
  //Fk must be allocated and filled up to KList.size();
  //  assert(t.size()==adjust.size());
  //  assert(t.size()==Basis.NumBasisElem());
  Matrix<mRealType> A;
  std::vector<mRealType> b;
  Matrix<mRealType> dcnk;
  int N = Basis.NumBasisElem(); //t.size();
  A.resize(N,N);
  b.resize(N,0.0);
  dcnk.resize(N,KList.size());
  //Fill in cnk.
  for (int n=0; n<N; n++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k = KList[ki][0];
      dcnk(n,ki) = Basis.dc_dk(n,k); //-Basis.c(n,k);
    }
  }
  //Fill in A and b
  A = 0.0;
  for (int l=0; l<N; l++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k2=KList[ki][0]*KList[ki][0];
  //    b[l] += k2*KList[ki][1]*(dFk[ki]-Fk[ki]) * dcnk(l, ki);
       b[l] += k2*KList[ki][1]*(dFk[ki]) * dcnk(l, ki);
      for (int n=0; n<N; n++)
        A(l,n) += k2*KList[ki][1]*dcnk(l,ki)*dcnk(n,ki);
    }
  }
  //Reduce for constraints
  int M = N;
  for (int i=0; i<N; i++)
    if (!adjust[i])
      M--;
  //The c is for "constrained"
  Matrix<mRealType> Ac;
  Ac.resize(M,M);
  std::vector<mRealType> bc(M,0.0), tc(M,0.0);
  //Build constrained Ac and bc
  int j=0;
  for (int col=0; col<N; col++)
  {
    if (adjust[col])
    {
      // Copy column a A to Ac
      int i=0;
      for (int row=0; row<N; row++)
        if (adjust[row])
        {
          Ac(i,j) = A(row,col);
          i++;
        }
      j++;
    }
    else
    {
      // Otherwise, subtract t(col)*A(:,col) from bc
      for (int row=0; row<N; row++)
        b[row] -= A(row,col)*t[col];
    }
  }
  j=0;
  for (int row=0; row<N; row++)
    if (adjust[row])
    {
      bc[j] = b[row];
      j++;
    }
  // Do SVD:
  // -------
  // Matrix<mRealType> U(M, M), V(M, M);
  // std::vector<mRealType> S(M), Sinv(M);
  // SVdecomp(Ac, U, S, V);
  ////////////////////////////////
  int m = Ac.rows();
  int n = Ac.cols();
  Matrix<mRealType> Atrans(n,m);
  Matrix<mRealType> U, V;
  U.resize(std::min(m,n),m);
  V.resize(n,std::min(m,n));
  std::vector<mRealType> S, Sinv;
  S.resize(std::min(n,m));
  //do the transpose
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
      Atrans(j,i) = Ac(i,j);
  }
  char JOBU = 'S';
  char JOBVT = 'S';
  int LDA = m;
  int LDU = m;
  int LDVT = std::min(m,n);
  int LWORK = 10*std::max(3*std::min(n,m)+std::max(m,n),5*std::min(m,n));
  std::vector<mRealType> WORK(LWORK);
  int INFO;
  LAPACK::gesvd(&JOBU,&JOBVT,&m,&n,Atrans.data(),&LDA,&S[0],U.data(),
         &LDU,V.data(),&LDVT,&WORK[0],&LWORK,&INFO);
  assert(INFO==0);
  int ur = U.rows();
  int uc = U.cols();
  Matrix<mRealType> Utrans(uc,ur);
  for(int i=0; i<ur; i++)
  {
    for(int j=0; j<uc; j++)
      Utrans(j,i) = U(i,j);
  }
  U.resize(uc,ur);
  U=Utrans;
  //////////////////////////////////
  // Zero out near-singular values
  mRealType Smax=S[0];
  for (int i=1; i<M; i++)
    Smax = std::max(S[i],Smax);
  for (int i=0; i<M; i++)
    if (S[i] < 0.0)
      std::cout << "negative singlar value.\n";
  //  perr << "Smax = " << Smax << std::endl;
  Sinv.resize(S.size());
  for (int i=0; i<M; i++)
    Sinv[i] = (S[i] < (tolerance*Smax)) ? 0.0 : (1.0/S[i]);
  int numSingular = 0;
  for (int i=0; i<Sinv.size(); i++)
    if (Sinv[i] == 0.0)
      numSingular++;
  if (numSingular > 0)
    std::cout << "There were " << numSingular << " singular values in breakup.\n";
  // Compute t_n, removing singular values
  for (int i=0; i<M; i++)
  {
    mRealType coef = 0.0;
    for (int j=0; j<M; j++)
      coef += U(j,i) * bc[j];
    coef *= Sinv[i];
    for (int k=0; k<M; k++)
      tc[k] += coef * V(k,i);
  }
  // Now copy tc values into t
  j=0;
  for (int i=0; i<N; i++)
    if (adjust[i])
    {
      t[i] = tc[j];
      j++;
    }
  // Calculate chi-squared
  mRealType Yk, chi2;
  chi2 = 0.0;
  for (int ki=0; ki<KList.size(); ki++)
  {
    Yk = dFk[ki]; //-Fk[ki];
    for (int n=0; n<N; n++)
    {
      Yk -= dcnk(n,ki)*t[n];
    }
    chi2 += KList[ki][0]*KList[ki][0]*KList[ki][1]*Yk*Yk;
  }
  return (chi2);
}

template<class BreakupBasis>
void LRBreakup<BreakupBasis>::DoAllBreakup(mRealType *chisqrlist, mRealType *Fk, mRealType *dFk, mRealType *t, mRealType *gt, mRealType *dt, mRealType *adjust)
{
  const mRealType tolerance = std::numeric_limits<mRealType>::epsilon();
  //t and adjust must be allocated up to Basis.NumBasisElem();
  //Fk must be allocated and filled up to KList.size();
  //  assert(t.size()==adjust.size());
  //  assert(t.size()==Basis.NumBasisElem());
  Matrix<mRealType> A;
  Matrix<mRealType> Af;
  Matrix<mRealType> As;
  
  std::vector<mRealType> b;
  std::vector<mRealType> bf;
  std::vector<mRealType> bs;
  
  Matrix<mRealType> cnk;
  Matrix<mRealType> dcnk;
  
  int N = Basis.NumBasisElem(); //t.size();
 
  A.resize(N,N);
  Af.resize(N,N);
  As.resize(N,N);
  
  b.resize(N,0.0);
  bf.resize(N,0.0);
  bs.resize(N,0.0);
  
  cnk.resize(N,KList.size());
  dcnk.resize(N,KList.size());
  //Fill in cnk.
  for (int n=0; n<N; n++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k = KList[ki][0];
      cnk(n,ki) = Basis.c(n,k);
      dcnk(n,ki) = Basis.dc_dk(n,k); //-Basis.c(n,k);
    }
  }
  //Fill in A and b
  A = 0.0;
  Af = 0.0;
  As = 0.0;
  
  for (int l=0; l<N; l++)
  {
    for (int ki=0; ki<KList.size(); ki++)
    {
      mRealType k2=KList[ki][0]*KList[ki][0];
      
      mRealType temp=KList[ki][1]*Fk[ki]*cnk(l,ki);
  //    b[l] += k2*KList[ki][1]*(dFk[ki]-Fk[ki]) * dcnk(l, ki);
      b[l] += temp;
      bf[l] += k2*temp;
      bs[l] += k2*KList[ki][1]*dFk[ki]* dcnk(l, ki);
      
      for (int n=0; n<N; n++)
      {
        temp=KList[ki][1]*cnk(l,ki)*cnk(n,ki);
        A(l,n) += temp;
        Af(l,n) += k2*temp;
        As(l,n) += k2*KList[ki][1]*dcnk(l,ki)*dcnk(n,ki);
      }
    }
  }
  //************************************
  //FOR POTENTIAL AND FORCE
  //************************************
  
  //Reduce for constraints
  int M = N;
  for (int i=0; i<N; i++)
    if (!adjust[i])
      M--;
  //The c is for "constrained"
  Matrix<mRealType> Ac;
  Matrix<mRealType> Afc;
  Matrix<mRealType> Asc;
  Ac.resize(M,M);
  Afc.resize(M,M);
  Asc.resize(M,M);
  std::vector<mRealType> bc(M,0.0), bfc(M,0.0), bsc(M,0.0),
                   tc(M,0.0), tfc(M,0.0), tsc(M,0.0);
  //Build constrained Ac and bc
  int j=0;
  for (int col=0; col<N; col++)
  {
    if (adjust[col])
    {
      // Copy column a A to Ac
      int i=0;
      for (int row=0; row<N; row++)
        if (adjust[row])
        {
          Ac(i,j) = A(row,col);
          Afc(i,j) = Af(row,col);
          Asc(i,j) = As(row,col);
          i++;
        }
      j++;
    }
    else
    {
      // Otherwise, subtract t(col)*A(:,col) from bc
      for (int row=0; row<N; row++)
      {
        b[row] -= A(row,col)*t[col];
        bf[row] -= Af(row,col)*gt[col];
        bs[row] -= As(row,col)*dt[col];
      }
    }
  }
  j=0;
  for (int row=0; row<N; row++)
    if (adjust[row])
    {
      bc[j] = b[row];
      bfc[j] = bf[row];
      bsc[j] = bs[row];
      j++;
    }
  // Do SVD:
  // -------
  // Matrix<mRealType> U(M, M), V(M, M);
  // std::vector<mRealType> S(M), Sinv(M);
  // SVdecomp(Ac, U, S, V);
  ////////////////////////////////
  int m = Ac.rows();
  int n = Ac.cols();
  Matrix<mRealType> A_trans(n,m);
  Matrix<mRealType> Af_trans(n,m);
  Matrix<mRealType> As_trans(n,m);
  Matrix<mRealType> U, V;
  Matrix<mRealType> Uf, Vf;
  Matrix<mRealType> Us, Vs;
  
  U.resize(std::min(m,n),m);
  V.resize(n,std::min(m,n));
  
  Uf.resize(std::min(m,n),m);
  Vf.resize(n,std::min(m,n));
  
  Us.resize(std::min(m,n),m);
  Vs.resize(n,std::min(m,n));
  
  std::vector<mRealType> S, Sinv;
  S.resize(std::min(n,m));
  
  std::vector<mRealType> Sf, Sfinv;
  Sf.resize(std::min(n,m));
  
  std::vector<mRealType> Ss, Ssinv;
  Ss.resize(std::min(n,m));
  //do the transpose
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
    {
      A_trans(j,i) = Ac(i,j);
      Af_trans(j,i) = Afc(i,j);
      As_trans(j,i) = Asc(i,j);
    }
  }
  char JOBU = 'S';
  char JOBVT = 'S';
  int LDA = m;
  int LDU = m;
  int LDVT = std::min(m,n);
  int LWORK = 10*std::max(3*std::min(n,m)+std::max(m,n),5*std::min(m,n));
  std::vector<mRealType> WORK(LWORK);
  int INFO;
  LAPACK::gesvd(&JOBU,&JOBVT,&m,&n,A_trans.data(),&LDA,&S[0],U.data(),
         &LDU,V.data(),&LDVT,&WORK[0],&LWORK,&INFO);
  assert(INFO==0);
  
  LAPACK::gesvd(&JOBU,&JOBVT,&m,&n,Af_trans.data(),&LDA,&Sf[0],Uf.data(),
         &LDU,Vf.data(),&LDVT,&WORK[0],&LWORK,&INFO);
  assert(INFO==0);

  LAPACK::gesvd(&JOBU,&JOBVT,&m,&n,As_trans.data(),&LDA,&Ss[0],Us.data(),
         &LDU,Vs.data(),&LDVT,&WORK[0],&LWORK,&INFO);
  assert(INFO==0);

  int ur = U.rows();
  int uc = U.cols();
  Matrix<mRealType> U_trans(uc,ur);
  Matrix<mRealType> Uf_trans(uc,ur);
  Matrix<mRealType> Us_trans(uc,ur);
  for(int i=0; i<ur; i++)
  {
    for(int j=0; j<uc; j++)
    {
      U_trans(j,i) = U(i,j);
      Uf_trans(j,i) = Uf(i,j);
      Us_trans(j,i) = Us(i,j);
    }
  }
  U.resize(uc,ur);
  U=U_trans;
  
  Uf.resize(uc,ur);
  Uf=Uf_trans;
  
  Us.resize(uc,ur);
  Us=Us_trans;
  //////////////////////////////////
  // Zero out near-singular values
  
  //First, do normal breakup.
  mRealType Smax=S[0];
  for (int i=1; i<M; i++)
    Smax = std::max(S[i],Smax);
  for (int i=0; i<M; i++)
    if (S[i] < 0.0)
      std::cout << "negative singlar value.\n";
  //  perr << "Smax = " << Smax << std::endl;
  Sinv.resize(S.size());
  for (int i=0; i<M; i++)
    Sinv[i] = (S[i] < (tolerance*Smax)) ? 0.0 : (1.0/S[i]);
  int numSingular = 0;
  for (int i=0; i<Sinv.size(); i++)
    if (Sinv[i] == 0.0)
      numSingular++;
  if (numSingular > 0)
    std::cout << "There were " << numSingular << " singular values in energy breakup.\n";
  // Compute t_n, removing singular values
  
  //Second, do force.
  Smax=Sf[0];
  for (int i=1; i<M; i++)
    Smax = std::max(Sf[i],Smax);
  for (int i=0; i<M; i++)
    if (Sf[i] < 0.0)
      std::cout << "negative singlar value.\n";
  //  perr << "Smax = " << Smax << std::endl;
  Sfinv.resize(Sf.size());
  for (int i=0; i<M; i++)
    Sfinv[i] = (Sf[i] < (tolerance*Smax)) ? 0.0 : (1.0/Sf[i]);
  numSingular = 0;
  for (int i=0; i<Sfinv.size(); i++)
    if (Sfinv[i] == 0.0)
      numSingular++;
  if (numSingular > 0)
    std::cout << "There were " << numSingular << " singular values in force breakup.\n";
  // Compute t_n, removing singular values
  
  //First, do strain.
  Smax=Ss[0];
  for (int i=1; i<M; i++)
    Smax = std::max(Ss[i],Smax);
  for (int i=0; i<M; i++)
    if (Ss[i] < 0.0)
      std::cout << "negative singlar value.\n";
  //  perr << "Smax = " << Smax << std::endl;
  Ssinv.resize(Ss.size());
  for (int i=0; i<M; i++)
    Ssinv[i] = (Ss[i] < (tolerance*Smax)) ? 0.0 : (1.0/Ss[i]);
  numSingular = 0;
  for (int i=0; i<Ssinv.size(); i++)
    if (Ssinv[i] == 0.0)
      numSingular++;
  if (numSingular > 0)
    std::cout << "There were " << numSingular << " singular values in strain breakup.\n";
  // Compute t_n, removing singular values
  
  for (int i=0; i<M; i++)
  {
    mRealType coef = 0.0;
    mRealType coef_f=0.0;
    mRealType coef_s=0.0;
    
    for (int j=0; j<M; j++)
    {
      coef += U(j,i) * bc[j];
      coef_f += Uf(j,i) * bfc[j];
      coef_s += Us(j,i) * bsc[j];
    }
    
    coef *= Sinv[i];
    coef_f *= Sfinv[i];
    coef_s *= Ssinv[i];
    
    for (int k=0; k<M; k++)
    {
      tc[k] += coef * V(k,i);
      tfc[k] += coef_f * Vf(k,i);
      tsc[k] += coef_s * Vs(k,i);
    }
  }
  // Now copy tc values into t
  j=0;
  for (int i=0; i<N; i++)
    if (adjust[i])
    {
      t[i] = tc[j];
      gt[i] = tfc[j];
      dt[i] = tsc[j];
      j++;
    }
  // Calculate chi-squared
  mRealType Yk(0.0), chi2(0.0);
  mRealType Yk_f(0.0), chi2_f(0.0);
  mRealType Yk_s(0.0), chi2_s(0.0);

  for (int ki=0; ki<KList.size(); ki++)
  {
    Yk = Fk[ki]; //-Fk[ki];
    Yk_f = Fk[ki];
    Yk_s= dFk[ki];
    
    for (int n=0; n<N; n++)
    {
      Yk -= cnk(n,ki)*t[n];
      Yk_f -= cnk(n,ki)*gt[n];
      Yk_s -= dcnk(n,ki)*dt[n];
    }
    chi2 += KList[ki][1]*Yk*Yk;
    chi2_f += KList[ki][0]*KList[ki][0]*KList[ki][1]*Yk_f*Yk_f;
    chi2_s += KList[ki][0]*KList[ki][0]*KList[ki][1]*Yk_s*Yk_s;
  }
//  std::vector<mRealType> chisqrtmp(3);
  
  chisqrlist[0]=chi2;
  chisqrlist[1]=chi2_f;
  chisqrlist[2]=chi2_s;
  
  //chisqrlist=chisqrtmp;

}

}

#endif
