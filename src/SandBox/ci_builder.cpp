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
using namespace qmcplusplus;

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
  int npeers=0;
  int count=0;
  excitations[0].write_node<max_states>(fout,0,count,excitations);
  count=0;
  ci_node<double> real_nodes;
  real_nodes.build_tree(64,excitations);
  real_nodes.set_peers(64,vmax);
  std::ofstream fout1("tree_1.xml");
  real_nodes.write_node(fout1);
  //copy constructor
  ci_node<double> copied(real_nodes);
  std::ofstream fout2("tree_2.xml");
  copied.write_node(fout2);
  return 0;
}
