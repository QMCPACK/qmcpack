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
    
    
#ifndef QMCPLUSPLUS_CI_BUILDER_H
#define QMCPLUSPLUS_CI_BUILDER_H
#include <bitset>
#include <vector>
#include <set>
#include <iostream>
#include <Message/Communicate.h>
#include <Utilities/IteratorUtility.h>
namespace qmcplusplus
{
template<unsigned CMAX>
struct ci_builder
{
  typedef std::set<unsigned long> iset_type;
  int max_vbm;
  int max_cbm;
  std::vector<int> v_offset;
  std::vector<iset_type*> occupied;
  std::vector<iset_type*> excited;

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

  ~ci_builder()
  {
    delete_iter(occupied.begin(),occupied.end());
    delete_iter(excited.begin(),excited.end());
  }

  /** create the single excitations
   *
   * @param exciations serialized list of excitations
   * @return the number of single generated max_vmc*max_cbm
   */
  template<typename NODETYPE> inline int singles(std::vector<NODETYPE>& excitations)
  {
    if(occupied.empty())
      occupied.push_back(new iset_type);
    if(excited.empty())
      excited.push_back(new iset_type);
    int level=0;
    iset_type& v_p(*occupied[level]);
    iset_type& c_p(*excited[level]);
    typedef std::bitset<CMAX> id_type;
    v_offset.clear();
    v_offset.push_back(1);
    for(int v=0; v<max_vbm; ++v)
    {
      id_type g;
      g.set(v,1);
      v_p.insert(g.to_ulong());
      for(int c=0; c<max_cbm; ++c)
      {
        id_type e;
        e.set(c,1);
        c_p.insert(e.to_ulong());
        //add child node to the root
        excitations[0].add_child(excitations.size());
        NODETYPE anode(g.to_ulong(),e.to_ulong());
        anode.from=v;
        anode.to=c;
        anode.my_id=excitations.size();
        excitations.push_back(anode);
      }
      v_offset.push_back(excitations.size());
    }
    return excitations.size();
  }

  /** promote from level CIs
   *
   * @param excitations container for the CIs upto level
   * @param level the reference level from which promotions are made
   */
  template<typename NODETYPE>
  inline int promote(std::vector<NODETYPE>& excitations, int level)
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
    typedef std::bitset<CMAX> id_type;
    for(typename iset_type::const_iterator i=v_base.begin(); i!=v_base.end(); ++i)
    {
      id_type v0(*i);
      for(int i=0; i<max_vbm; ++i)
      {
        if(v0[i])
          continue;
        id_type v1(v0);
        v1.set(i,1);
        v_p.insert(v1.to_ulong());
      }
    }
    for(typename iset_type::const_iterator i=c_base.begin(); i!=c_base.end(); ++i)
    {
      id_type v0(*i);
      for(int i=0; i<max_cbm; ++i)
      {
        if(v0[i])
          continue;
        id_type v1(v0);
        v1.set(i,1);
        c_p.insert(v1.to_ulong());
      }
    }
    std::vector<int> v_tmp;
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
          anode.my_id=excitations.size();
          excitations.push_back(anode);//(NODETYPE(*i,*j));
        }
      }
      v_tmp.push_back(excitations.size());
    }
    v_offset=v_tmp;
    //return the number of excitations created
    return excitations.size();
  }

  inline int changed_bit(std::bitset<CMAX>& x)
  {
    int count = 0;
    while(!x[count])
      count++;
    return count;
  }

  template<typename NODETYPE>
  int find_parent(NODETYPE& anode,std::vector<NODETYPE>& p)
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
      APP_ABORT("Parent is not found. Error ");
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
      APP_ABORT("Parent is not found. Error ");
    }
    anode.from=from;
    anode.to=to;
    return anode.parent_id=j-1;
  }
};
}
#endif
