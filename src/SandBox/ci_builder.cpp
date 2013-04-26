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
      if(excitations[i].closed())
        continue;
      nparent++;
    }
  }
  ofstream fout("tree.xml");
  int npeers=0;
  int count=0;
  excitations[0].write_node<max_states>(fout,0,count,excitations);
  count=0;
  ci_node<double> real_nodes;
  real_nodes.build_tree(64,excitations);
  real_nodes.set_peers(64,vmax);
  ofstream fout1("tree_1.xml");
  real_nodes.write_node(fout1);
  //copy constructor
  ci_node<double> copied(real_nodes);
  ofstream fout2("tree_2.xml");
  copied.write_node(fout2);
  return 0;
}
