//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
using namespace std;

template<typename CT, typename T>
inline void check_ratios(const CT& a, const CT& b, T eps)
{
  T del;
  cout.setf(std::ios::scientific, std::ios::floatfield);
  for(int i=0; i<a.size(); ++i)
    if(abs(del=(a[i]/b[i]-1.0))>eps)
     cout << i << " recursive=" << a[i] << " error=" << del << endl;
}

int main(int argc, char** argv)
{
  const int max_states=8;
  const int vmax=4;
  const int cmax=8;
  typedef ci_node_proxy node_type;
  vector<node_type> excitations;

  int multiplet=3;
  int nparent=0;
  int nchildren=0;

  vector<int> det_offset(multiplet+2,1);
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
      if(excitations[i].closed()) continue;
      nparent++;
    }
  }

  ofstream fout("tree.xml");
  int count=0;
  excitations[0].write_node<max_states>(fout,0,count,excitations);

  count=0;
  const int M=8;
  ci_node<double> CI;
  CI.build_tree(M,excitations);
  CI.set_peers(M,vmax);

  ofstream fout1("tree_1.xml");
  CI.write_node(fout1);

  typedef Matrix<double> mat_t;
  typedef Vector<double> vec_t;

  mat_t psi0(M,M), psiv(M,M), psic(cmax,M), Identity(M,M);
  vec_t psiv_big(M+cmax);

  //build up big  psi
  for(int i=0; i<psiv.size();++i) psiv(i)=Random();
  for(int i=0; i<psic.size();++i) psic(i)=Random();
  psi0=psiv;//save

  vector<double> dets(excitations.size()),ratios(excitations.size()), ptcl_ratios(excitations.size());
  const bool substitute_col=false;
  const double eps=1e-12;
  double det_0,det_inv;

  det_0=CI.getRatios(psiv,psic,ratios,substitute_col);
  CI.debugRatios(psiv,psic,dets,substitute_col);

  for(int i=0; i<ratios.size(); ++i) ratios[i]*=det_0;
  check_ratios(dets,ratios,eps);

  //pick a particle index
  int iat=0;
  CI.ratioByRowSubstitution(psiv_big,iat,ptcl_ratios);

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
