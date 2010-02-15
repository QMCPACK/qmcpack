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
  int npeers=0;
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
  mat_t psi_big(M+cmax,M), psi_test(M,M), Identity(M,M);

  //build up big  psi
  for(int i=0; i<psi_big.size();++i) psi_big(i)=Random();

  std::copy(psi_big[0],psi_big[M],psi_test.data());
  CI.inverse=psi_test;
  double det_0=invert_matrix(CI.inverse,true);

  vector<double> dets(excitations.size()),ratios(excitations.size());

  //calculate the determinants by replacement
  dets[0]=det_0;
  CI.debugRatios(psi_test,psi_big,dets,true);

  //calculate the determinants by replacement
  CI.getRatios(psi_big,ratios,true);

  cout << "Checking row replacement " << endl;
  cout.setf(std::ios::scientific, std::ios::floatfield);
  cerr.setf(std::ios::scientific, std::ios::floatfield);
  double det_inv=1.0/det_0,del;
  double eps=1e-12;//1.e3*numeric_limits<double>::epsilon();
  for(int i=0; i<ratios.size(); ++i) dets[i]*=det_inv;
  for(int i=0; i<ratios.size(); ++i)
    if(abs(del=(ratios[i]/dets[i]-1.0))>eps)
     cerr << i << " recursive=" << ratios[i] << " error=" << del << endl;

  cout << "Checking col replacement " << endl;
  dets[0]=det_0;
  CI.debugRatios(psi_test,psi_big,dets,false);
  CI.getRatios(psi_big,ratios,false);
  for(int i=0; i<ratios.size(); ++i) dets[i]*=det_inv;
  for(int i=0; i<ratios.size(); ++i)
    if(abs(del=(ratios[i]/dets[i]-1.0))>eps)
     cerr << i << " recursive=" << ratios[i] << " error=" << del << endl;
  return 0;
}
