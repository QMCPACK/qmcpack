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
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include <QMCWaveFunctions/Fermion/excitation_node.h>
#include <QMCWaveFunctions/Fermion/ci_builder.h>

using namespace qmcplusplus;
using namespace std;

int main(int argc, char** argv)
{
  const int max_states=8;
  const int vmax=4;
  const int cmax=8;
  typedef excitation_node node_type;
  vector<node_type> excitations;

  //add zero
  excitations.push_back(node_type());

  ci_builder<max_states> ci_maker(vmax,cmax);

  int multiplet=3;
  ci_maker.singles(excitations);
  for(int level=0; level<multiplet; ++level) ci_maker.promote(excitations,level);

  int nparent=0;
  int nchildren=0;
  for(int i=0; i<excitations.size(); ++i) 
  {
    nchildren+=excitations[i].children.size();
    if(excitations[i].closed()) continue;
    nparent++;
  }

  int level=0;
  int count=0;
  cout << "<ci_basis max_level=\"" << multiplet+1 << "\" size=\"" << nchildren+1 << "\">" << endl;
  cout << "<!-- " <<endl;
  cout << " Number of child nodes " << nchildren << "== " << excitations.size()-1 << endl;
  cout << " Number of parent nodes " << nparent << endl;
  cout << "-->" <<endl;
  cout << "<ci_vector>" << endl;
  for(int i=0; i<excitations.size(); ++i) 
  {
    node_type cur(excitations[i]);
    int pid=cur.ref_state;
    node_type parent(excitations[pid]);
    bitset<cmax> valence(cur.ground);
    cout << "<node level=\"" << valence.count()
      << "\" g=\"" <<  valence
      << "\" e=\"" << bitset<max_states>(cur.excited)
      << "\" g_id=\"" << cur.ground
      << "\" e_id=\"" << cur.excited
      << "\" loc=\"" << i
      << "\" p_id=\"" << pid
      << "\" from=\"" << cur.from
      << "\" to=\"" << cur.to
      << "\"/>" 
      << endl;
  }
  cout << "</ci_vector>" << endl;

  cout << "<ci_tree>" << endl;
  excitations[0].write_node<max_states>(level,count,excitations);
  cout << "</ci_tree>" << endl;
  cout << "</ci_basis>" << endl;
  return 0;
}
