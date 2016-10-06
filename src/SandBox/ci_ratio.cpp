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
    
    



/**@file excitation.cpp
 * @brief Test code for multideterminant tree
 */
#include <Utilities/OhmmsInfo.h>
#include <Utilities/RandomGenerator.h>
#include <Utilities/Timer.h>
#include <Message/Communicate.h>
#include <Message/OpenMP.h>
#include <Numerics/MatrixOperators.h>
#include <Numerics/DeterminantOperators.h>
#include <QMCWaveFunctions/Fermion/ci_node.h>
#include <QMCWaveFunctions/Fermion/ci_builder.h>
#include <Utilities/IteratorUtility.h>
#include <fstream>
#include <limits>
using namespace qmcplusplus;

template<typename CT, typename T>
inline void check_ratios(const CT& a, const CT& b, T eps)
{
  T del;
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  for(int i=0; i<a.size(); ++i)
    if(std::abs(del=(a[i]/b[i]-1.0))>eps)
      std::cout << i << " recursive=" << a[i] << " error=" << del << std::endl;
}

int main(int argc, char** argv)
{
  const int max_states=8;
  const int vmax=4;
  const int cmax=8;
  typedef ci_node_proxy node_type;
  std::vector<node_type> excitations;
  int multiplet=3;
  int nparent=0;
  int nchildren=0;
  std::vector<int> det_offset(multiplet+2,1);
  //add zero
  excitations.push_back(node_type());
  {
    ci_builder<max_states> ci_maker(vmax,cmax);
    //single is treated differently
    det_offset[1]=ci_maker.singles(excitations);
    for(int level=0,k=1; level<multiplet; ++level,++k)
      det_offset[k+1]=ci_maker.promote(excitations,level);
    for(int i=0; i<excitations.size(); ++i)
    {
      nchildren+=excitations[i].children.size();
      if(excitations[i].closed())
        continue;
      nparent++;
    }
  }
  std::ofstream fout("tree.xml");
  int count=0;
  excitations[0].write_node<max_states>(fout,0,count,excitations);
  count=0;
  const int M=8;
  ci_node<double> CI;
  CI.build_tree(M,excitations);
  CI.set_peers(M,vmax);
  std::ofstream fout1("tree_1.xml");
  CI.write_node(fout1);
  typedef Matrix<double> mat_t;
  typedef Vector<double> vec_t;
  mat_t psi0(M,M), psiv(M,M), psic(cmax,M), Identity(M,M);
  vec_t psiv_big(M+cmax);
  std::vector<mat_t*> inverse(excitations.size());
  for(int i=0; i<inverse.size(); ++i)
    inverse[i]=new mat_t(M,M);
  //build up big  psi
  for(int i=0; i<psiv.size(); ++i)
    psiv(i)=Random();
  for(int i=0; i<psic.size(); ++i)
    psic(i)=Random();
  for(int i=0; i<psiv_big.size(); ++i)
    psiv_big(i)=Random();
  psi0=psiv;//save
  int mdtot=excitations.size();
  std::vector<double> dets(mdtot),detbyratio(mdtot);//store determinants
  std::vector<double> ratios_rec(mdtot), ptcl_ratios(mdtot), ptcl_ratios_rec(mdtot);
  const bool substitute_col=true;
  const double eps=1e-12;
  double det_0,det_inv;
  //get multi determinant values using recursive
  CI.inverse=psiv;
  det_0=invert_matrix(CI.inverse,true);
  CI.getRatios(psic,ratios_rec,substitute_col);
  CI.debugRatios(psiv,psic,dets,inverse,substitute_col);
  std::cout << "Comparing promotion " << std::endl;
  for(int i=0; i<mdtot; ++i)
    detbyratio[i]=det_0*ratios_rec[i];
  check_ratios(dets,detbyratio,eps);
  //consider a particle update
  int iat=0;
  CI.debugRatioByRowSubstitution(psiv_big,iat,ptcl_ratios,inverse);
  CI.inverse=*inverse[0];
  //CI.ratioByRowSubstitution(psic,psiv_big,iat,ptcl_ratios_rec);
  CI.ratioByRowSubstitution(psiv_big,iat,ptcl_ratios_rec);
  std::cout << "Comparing particle-update ratios " << std::endl;
  check_ratios(ptcl_ratios,ptcl_ratios_rec,eps);
  {
    //use direct method
    psiv=psi0;
    for(int j=0; j<M; ++j)
      psiv(j,iat)=psiv_big(j);
    for(int j=0; j<cmax; ++j)
      psic(j,iat)=psiv_big(j+M);
    CI.debugRatios(psiv,psic,dets,inverse,substitute_col);
    for(int i=0; i<mdtot; ++i)
      detbyratio[i]=det_0*ratios_rec[i]*ptcl_ratios[i];
    std::cout << "Comparing particle-update against direct method" << std::endl;
    check_ratios(dets,detbyratio,eps);
  }
  Timer clock;
  for(int i=0; i<1000; ++i)
    CI.debugRatioByRowSubstitution(psiv_big,iat,ptcl_ratios,inverse);
  std::cout << "Time for direct ratio = " << clock.elapsed() << std::endl;
  clock.restart();
  for(int i=0; i<1000; ++i)
    CI.ratioByRowSubstitution(psiv_big,iat,ptcl_ratios);
  std::cout << "Time for recusrive ratio = " << clock.elapsed() << std::endl;
  //for(int i=0; i<ratios.size(); ++i) detbyratio[i]=det_0*ratios[i]*ptcl_ratios[i];
  //check_ratios(dets,detbyratio,eps);
  //psi_big=psi0;
  //for(int i=0; i<psi_big.rows(); ++i) psi_big(i,iat)=psiv_big(i);
  //std::copy(psi_big[0],psi_big[M],psi_test.data());
  //CI.inverse=psi_test;
  //det_0=invert_matrix(CI.inverse,true);
  //dets[0]=det_0;
  //CI.debugRatios(psi_test,psi_big,dets,substitute_col);
  //CI.getRatios(psi_big,ratios,substitute_col);
  //det_inv=1.0/det_0;
  //for(int i=0; i<ratios.size(); ++i) dets[i]*=det_inv;
  //check_ratios(dets,ratios,eps);
  return 0;
}
