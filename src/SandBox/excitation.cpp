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
#include <bitset>
#include <set>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"

using namespace qmcplusplus;

/** a node of multideterminant tree
 */
struct excitation_node
{
  unsigned long ground;
  unsigned long excited;
  int ref_state;
  int from;
  int to;
  vector<int> children;

  /** default constructor
   * Invalidate the excitation
   */
  inline excitation_node():ground(0u),excited(0u),ref_state(0),from(0),to(0) {}

  inline explicit excitation_node(unsigned long vid, unsigned long cid)
    : ground(vid),excited(cid),ref_state(0),from(0),to(0) {}

  inline void add_child(int i)
  {
    children.push_back(i);
  }

  inline bool closed() const
  {
    return children.empty();
  }

  template<unsigned CMAX>
  inline void traverse(int & level, int count, vector<excitation_node>& ci)
  {
    //if(children.size()) 
    std::cout << "<node level=\""<<level
      << "\" g=\"" << std::bitset<CMAX>(ground) 
      << "\" e=\"" << std::bitset<CMAX>(excited) 
      << "\" gid=\"" << ground
      << "\" eid=\"" << excited
      << "\" from=\"" << from 
      << "\" to=\""<< to 
      <<"\" loc=\"" << count
      << "\"";
    if(children.size()) std::cout << ">" << std::endl;
    else std::cout << "/>" << std::endl;
    for(int i=0; i<children.size(); ++i) 
    {
      int next=level+1;
      ci[ children[i] ].traverse<CMAX>(next,children[i],ci);
    }
    if(children.size()) std::cout << "</node>"<<endl;
  }
}; 


using namespace std;

template<unsigned CMAX>
struct ci_builder
{
  typedef set<unsigned long> iset_type;
  int max_vbm;
  int max_cbm;
  vector<int> v_offset;
  vector<iset_type*> occupied;
  vector<iset_type*> excited;

  /** default constructor
   *
   * @param vmax maximum valance states 
   * @param cmax maximum conduction states
   */
  inline explicit ci_builder(int vmax=1, int cmax=CMAX)
    :max_vbm(vmax), max_cbm(cmax)
  {
    v_offset.reserve(CMAX); 
  }

  /** create the single excitations */
  template<typename NODETYPE> inline void singles(vector<NODETYPE>& excitations)
  {
    if(occupied.empty()) occupied.push_back(new iset_type);
    if(excited.empty()) excited.push_back(new iset_type);

    int level=0;
    iset_type& v_p(*occupied[level]);
    iset_type& c_p(*excited[level]);

    typedef bitset<CMAX> id_type;
    v_offset.clear();
    v_offset.push_back(1);
    for(int v=0; v<max_vbm; ++v)
    {
      id_type g; g.set(v,1);
      v_p.insert(g.to_ulong());
      for(int c=0; c<max_cbm;++c) 
      {
        id_type e; e.set(c,1);
        c_p.insert(e.to_ulong());
        //add child node to the root
        excitations[0].add_child(excitations.size());
        NODETYPE anode(g.to_ulong(),e.to_ulong());
        anode.from=v; anode.to=c;
        excitations.push_back(anode);
      }
      v_offset.push_back(excitations.size());
    }
  }

  /** promote from level CIs
   *
   * @param excitations container for the CIs upto level
   * @param level the reference level from which promotions are made
   */
  template<typename NODETYPE>
  inline void promote(vector<NODETYPE>& excitations, int level)
  {
    int dn=occupied.size()-level;
    while(dn)//allocate a new level
    {
      occupied.push_back(new iset_type);
      excited.push_back(new iset_type);
      dn--;
    }

    const iset_type& v_base(*occupied[level]);
    const iset_type& c_base(*excited[level]);
    iset_type& v_p(*occupied[level+1]);
    iset_type& c_p(*excited[level+1]);

    typedef bitset<CMAX> id_type;
    for(typename iset_type::const_iterator i=v_base.begin(); i!=v_base.end(); ++i)
    {
      id_type v0(*i);
      for(int i=0; i<max_vbm; ++i) 
      {
        if(v0[i]) continue;
        id_type v1(v0); v1.set(i,1); v_p.insert(v1.to_ulong());
      }
    }

    for(typename iset_type::const_iterator i=c_base.begin(); i!=c_base.end(); ++i)
    {
      id_type v0(*i);
      for(int i=0; i<max_cbm; ++i) 
      {
        if(v0[i]) continue;
        id_type v1(v0); v1.set(i,1); c_p.insert(v1.to_ulong());
      }
    }

    vector<int> v_tmp;
    v_tmp.push_back(excitations.size());
    for(typename iset_type::iterator i=v_p.begin(); i!=v_p.end(); ++i)
    {
      for(typename iset_type::iterator j=c_p.begin(); j!=c_p.end(); ++j)
      {
        NODETYPE anode(*i,*j);
        int pid=find_parent(anode,excitations);
        if(pid>=0)
        {
          excitations[pid].add_child(excitations.size());
          excitations.push_back(anode);//(NODETYPE(*i,*j));
        }
      }
      v_tmp.push_back(excitations.size());
    }
    v_offset=v_tmp;
  }

  inline int changed_bit(std::bitset<CMAX>& x)
  {
    int count = 0;
    while(!x[count]) count++;
    return count;
  }

  template<typename NODETYPE>
  int find_parent(NODETYPE& anode,vector<NODETYPE>& p)
  {
    bool findv=true;
    std::bitset<CMAX> vid(anode.ground);
    std::bitset<CMAX> cid(anode.excited);
    unsigned long vid_p;
    unsigned long cid_p;
    int from;
    int to;

    int i=0;
    int imax=v_offset.size();
    while(findv && i<imax)
    {
      int where=v_offset[i];
      std::bitset<CMAX> parent(p[where].ground);
      parent^=vid;
      if(parent.count()==1) 
      {
        findv=false;
        vid_p=p[where].ground;
        from=changed_bit(parent);
      }
      ++i;
    }

    if(findv) 
    {
      cerr << " Parent is not found. Error " << endl;
      return -1;
    }

    findv=true;
    int j=v_offset[i-1];
    int jmax=v_offset[i];
    while(findv && j<jmax)
    {
      std::bitset<CMAX> parent(p[j].excited);
      parent^=cid;
      if(parent.count() ==1)
      {
        findv=false;
        cid_p=p[j].excited;
        to=changed_bit(parent);
      }
      ++j;
    }

    if(findv)
    {
      cerr << " Parent is not found. Error " << endl;
      return -1;
    }
    anode.from=from;
    anode.to=to;
    return anode.ref_state=j-1;
  }
};


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
  excitations[0].traverse<max_states>(level,count,excitations);
  cout << "</ci_tree>" << endl;
  cout << "</ci_basis>" << endl;
  return 0;
}
