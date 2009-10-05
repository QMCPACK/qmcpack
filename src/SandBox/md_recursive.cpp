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
#include <QMCWaveFunctions/Fermion/excitation_node.h>
#include <QMCWaveFunctions/Fermion/ci_builder.h>
#include <Utilities/IteratorUtility.h>
using namespace qmcplusplus;
using namespace std;

int main(int argc, char** argv)
{
  const int max_states=8;
  const int vmax=4;
  const int cmax=4;
  typedef excitation_node<double> node_type;
  vector<node_type> excitations;

  int multiplet=2;
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

  int M=8;
  typedef Matrix<double> mat_t;
  vector<mat_t*> dets(excitations.size());
  mat_t psi_big(M+cmax,M), psi_test(M,M);

  for(int i=0; i<psi_big.size();++i) psi_big(i)=Random();
  //allocate the inverse matrix
  for(int i=0; i<excitations.size(); ++i) excitations[i].resize(M);
  for(int i=0; i<dets.size(); ++i) dets[i]=new mat_t(M,M);

  mat_t uv(vmax,cmax);
  double phase0,phase1;

  //psi(i,j) for i=orbital index, j=particle index offset by the starting particle
  //psi_big.rows()=occupied states + unoccupied states.
  std::copy(psi_big[0],psi_big[M],dets[0]->data());

  //double logdet_0=InvertWithLog(dets[0]->data(),M,M,phase0);
  double det_0=invert_matrix(*dets[0],true);

  for(int i=0; i<vmax; ++i)
    for(int j=0; j<cmax; ++j)
      uv(i,j)=BLAS::dot(M,(*dets[0])[M-i-1],psi_big[j+M]);

  vector<double> ratios(excitations.size(),1.0), ratios_impl(excitations.size(),1.0);

  //save the inverse
  excitations[0].inverse=(*dets[0]);

  //using implicit method
  excitations[0].get_ratios_root(psi_big,excitations,ratios_impl,M-vmax);

  //using bebug method
  excitations[0].get_ratios_root_debug(psi_big,excitations,ratios);

  for(int i=0; i<ratios.size(); ++i)
    cout << "Ratio " << i << " " << ratios_impl[i] << " " << ratios[i] << endl;


  for(int i=0; i<excitations.size(); ++i) 
  {
    node_type cur(excitations[i]);
    int pid=cur.parent_id;
    node_type parent(excitations[pid]);
    bitset<cmax> valence(cur.ground);
    cout << "<node level=\"" << valence.count()
      << "\" g=\"" <<  valence
      << "\" e=\"" << bitset<max_states>(cur.excited)
      << "\" g_id=\"" << cur.ground
      << "\" e_id=\"" << cur.excited
      << "\" my_id=\"" << cur.my_id
      << "\" p_id=\"" << pid
      << "\" from=\"" << cur.from
      << "\" to=\"" << cur.to
      << "\"/>" 
      << " ratio = " << ratios[i] 
      << endl;
  }

  for(int i=det_offset[0]; i<det_offset[1];++i)
  {
    int v=M-excitations[i].from-1;
    int c=M+excitations[i].to;

    std::copy(psi_big[0],psi_big[M],dets[i]->data());
    for(int j=0; j<M; ++j) (*dets[i])(j,v)=psi_big(c,j);
    //double logdet_1=InvertWithLog(dets_ref[i]->data(),M,M,phase1);
    double det_1=invert_matrix(*dets[i],true);
    psi_test=*dets[i] - excitations[i].inverse;
    cout << setw(4) << i << " from=" <<setw(4) <<  v << " to=" <<setw(4) <<  c << " "
      << setw(20) << ratios[i]<< setw(20) << (1.0-det_1/(det_0*ratios[i]))
      << setw(20) << BLAS::norm2(psi_test.size(),psi_test.data()) << endl;
  }

  cout << "----------------------------------------- " << endl;
  cout << "Double excitation " << det_offset[1] << "-" << det_offset[2] << endl;
  for(int i=det_offset[1]; i<det_offset[2];++i)
  {
    int parent=excitations[i].parent_id;
    if(parent==0) continue;
    int v0=M-excitations[parent].from-1;
    int c0=M+excitations[parent].to;
    int v=M-excitations[i].from-1;
    int c=M+excitations[i].to;
    std::copy(psi_big[0],psi_big[M],dets[i]->data());
    for(int j=0; j<M; ++j) (*dets[i])(j,v0)=psi_big(c0,j);
    for(int j=0; j<M; ++j) (*dets[i])(j,v)=psi_big(c,j);
    double det_1=invert_matrix(*dets[i],true);
    cout <<setw(4) <<  i << " from=" << setw(4) << v0 << " " << setw(4) << v << " to=" <<setw(4) <<  c0 <<setw(4) <<  c  << " "
      << setw(20) << ratios[i]<< setw(20) << (1.0-det_1/(det_0*ratios[i]))
      << setw(20) << BLAS::norm2(psi_test.size(),psi_test.data()) << endl;
  }

  cout << "----------------------------------------- " << endl;
  cout << "Triple excitation " << det_offset[2] << "-" << det_offset[3] << endl;
  for(int i=det_offset[2]; i<det_offset[3];++i)
  {
    int parent=excitations[i].parent_id;
    int grand_parent=excitations[parent].parent_id;
    int v0=M-excitations[grand_parent].from-1;
    int c0=M+excitations[grand_parent].to;
    int v1=M-excitations[parent].from-1;
    int c1=M+excitations[parent].to;
    int v=M-excitations[i].from-1;
    int c=M+excitations[i].to;
    std::copy(psi_big[0],psi_big[M],dets[i]->data());
    for(int j=0; j<M; ++j) (*dets[i])(j,v0)=psi_big(c0,j);
    for(int j=0; j<M; ++j) (*dets[i])(j,v1)=psi_big(c1,j);
    for(int j=0; j<M; ++j) (*dets[i])(j,v)=psi_big(c,j);
    double det_1=invert_matrix(*dets[i],true);
    cout <<setw(4) <<  i << " from=" << setw(4) << v0 << setw(4) << v1 << setw(4) << v 
      << " to=" <<setw(4) << c0 << setw(4) << c1 << setw(4) <<  c  << " "
      << setw(20) << ratios[i]<< setw(20) << (1.0-det_1/(det_0*ratios[i]))
      << setw(20) << BLAS::norm2(psi_test.size(),psi_test.data()) << endl;
  }

  cout << "----------------------------------------- " << endl;
  cout << "Qud excitation " << det_offset[3] << "-" << det_offset[4] << endl;
  for(int i=det_offset[3]; i<det_offset[4];++i)
  {
    int parent=excitations[i].parent_id;
    int parent1=excitations[parent].parent_id;
    int parent2=excitations[parent1].parent_id;
    int v0=M-excitations[parent2].from-1;
    int c0=M+excitations[parent2].to;
    int v1=M-excitations[parent1].from-1;
    int c1=M+excitations[parent1].to;
    int v2=M-excitations[parent].from-1;
    int c2=M+excitations[parent].to;
    int v=M-excitations[i].from-1;
    int c=M+excitations[i].to;
    std::copy(psi_big[0],psi_big[M],dets[i]->data());
    for(int j=0; j<M; ++j) (*dets[i])(j,v0)=psi_big(c0,j);
    for(int j=0; j<M; ++j) (*dets[i])(j,v1)=psi_big(c1,j);
    for(int j=0; j<M; ++j) (*dets[i])(j,v2)=psi_big(c2,j);
    for(int j=0; j<M; ++j) (*dets[i])(j,v)=psi_big(c,j);
    double det_1=invert_matrix(*dets[i],true);
    cout <<setw(4) <<  i << " from=" << setw(4) << v0 << setw(4) << v1 << setw(4) << v2 << setw(4) << v 
      << " to=" <<setw(4) << c0 << setw(4) << c1 << setw(4) << c2 << setw(4) <<  c  << " "
      << setw(20) << ratios[i]<< setw(20) << (1.0-det_1/(det_0*ratios[i]))
      << setw(20) << BLAS::norm2(psi_test.size(),psi_test.data()) << endl;
  }

  delete_iter(dets.begin(),dets.end());
  return 0;
}
