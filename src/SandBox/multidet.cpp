//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Numerics/MatrixOperators.h"
#include "SandBox/determiant_operators.h"
#include "Utilities/IteratorUtility.h"
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo Welcome(argc,argv,OHMMS::Controller->rank());
  int M=64;
  int nc=16;
  int ic=1;
  int unocc=-1;
  int niters=10000;
  while(ic<argc)
  {
    std::string c(argv[ic]);
    if(c=="-v")//number of valence states
      M=atoi(argv[++ic]);
    else
      if(c== "-c")//number of conduction states
        nc=atoi(argv[++ic]);
      else
        if(c=="-u")//valence state from the top, starting zero
          unocc=M-atoi(argv[++ic])-1;
        else
          if(c=="-i")//number of iterations
            niters=atoi(argv[++ic]);
    ++ic;
  }
  if(unocc<0)
    unocc=M-1;
  int bigM=M+nc;
  Matrix<double>  psi_0(M,M),psi_1(M,M),psi_2(M,M),psi_saved(M,M),Identity(M,M);
  Matrix<double>  psi_big(bigM,M);
  Vector<double> psi_v(M);
  for(int i=0; i<psi_big.size(); ++i)
    psi_big(i)=Random();
  for(int i=0; i<M; ++i)
    psi_v(i)=Random();
  copy(psi_big[0],psi_big[M],psi_0.data());
  psi_saved=psi_0;
  {
    psi_1=psi_0;
    double det_0=invert_matrix(psi_0,true);
    Vector<double> newcols(bigM);
    for(int i=0; i<bigM; ++i)
      newcols(i)=Random();
    //replace a column
    //double rc=DetRatioTranspose(psi_0,newcols.begin(),0);
    //double rc=getRatioByColSubstitution(psi_0.data(),newcols.data(),M,0);
    double rc=getRatioByColSubstitution(psi_0.data(),newcols.data(),M,0);
    //replace this
    for(int i=0; i<M; ++i)
      psi_1(0,i)=newcols(i);
    double det_1=invert_matrix(psi_1,true);
    std::cout << "Checking column substitution= " << rc
         << " " << det_1/det_0 << std::endl;
  }
  //save the original matrices
  psi_1=psi_saved;
  double phase=0;
  double logdet_0=InvertWithLog(psi_0.data(),M,M,phase);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  MatrixOperators::product(psi_0,psi_saved,Identity);
  for(int i=0; i<M; ++i)
    Identity(i,i)-=1.0;
  std::cout << "Checking identity " << BLAS::norm2(Identity.size(),Identity.data())<< std::endl;
  Vector<double> ratios(nc), ratios_0(nc);
  std::vector<int> imap(nc);
  for(int i=0; i<nc; ++i)
    imap[i]=i;
  Timer myclock;
  for(int iter=0; iter<niters; ++iter)
    getRatiosByRowSubstitution(psi_big[M],psi_0[unocc],ratios_0.data(),M,imap);
  double t_naive=myclock.elapsed();
  myclock.restart();
  for(int iter=0; iter<niters; ++iter)
    getRatiosByRowSubstitution(psi_big[M],psi_0[unocc],ratios.data(),M,nc);
  double t_gemv=myclock.elapsed();
  std::cout << "Timing of ratios gemv= " << t_gemv << " naive=" << t_naive << " better=" << t_naive/t_gemv<< std::endl;
  ratios_0 -= ratios;
  std::cout << "Error in the ratios= " <<  BLAS::norm2(ratios_0.size(),ratios_0.data())<< std::endl;
  //Matrix<double> newinv(nc,M*M);
  std::vector<Matrix<double>*> newinv_v(nc);
  for(int i=0; i<nc; ++i)
    newinv_v[i]=new Matrix<double>(M,M);
  double t_multi=0.0, t_dir=0.0;
  for(int iter=0; iter<niters; ++iter)
  {
    for(int i=0; i<psi_big.size(); ++i)
      psi_big(i)=Random();
    myclock.restart();
    getRatiosByRowSubstitution(psi_big[M],psi_0[unocc],ratios.data(),M,imap);
    //multidet_row_update(psi_0.data(),psi_big[M],ratios.data(),newinv.data(),M,unocc,nc);
    //multidet_row_update(psi_0.data(),psi_big[M],ratios.data(),newinv.data(),M,unocc,imap);
    multidet_row_update(psi_0,psi_big,ratios,newinv_v,M,unocc,imap);
    t_multi+=myclock.elapsed();
    for(int k=0; k<imap.size(); ++k)
    {
      int i=imap[k]+M;
      //std::copy(newinv[k],newinv[k+1],psi_2.data());
      copy(newinv_v[k]->begin(),newinv_v[k]->end(),psi_2.data());
      psi_1=psi_saved;
      myclock.restart();
      for(int j=0; j<M; ++j)
        psi_1(j,unocc)=psi_big(i,j);
      double newphase;
      double logdet_1=InvertWithLog(psi_1.data(),M,M,newphase);
      t_dir+=myclock.elapsed();
      if(iter==4)
      {
        psi_2-=psi_1;
        std::cout << "Error in Inverse matrix = " << BLAS::norm2(psi_2.size(),psi_2.data()) << std::endl;
        std::cout << "Inverse matrix norm2 = " <<BLAS::norm2(psi_1.size(),psi_1.data()) << std::endl;
        if(phase==newphase)//too lazy
          std::cout << "DetRatioTranspose = " << ratios[k] << " error=" << ratios[k]-std::exp(logdet_1-logdet_0) << std::endl;
      }
    }
  }
  std::cout << M << " " << nc << " " << t_dir/niters << " " << t_multi/niters << " " << t_dir/t_multi
       << std::endl;
  for(int i=0; i<nc; ++i)
    delete newinv_v[i];
  OHMMS::Controller->finalize();
  return 0;
}

