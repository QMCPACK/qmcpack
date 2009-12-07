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
#ifndef QMCPLUSPLUS_CI_NODE_PROXY_H
#define QMCPLUSPLUS_CI_NODE_PROXY_H
#include <bitset>
#include <vector>
#include <iostream>

namespace qmcplusplus
{
  /** helper class to build ci recursively
   */
  struct ci_node_proxy
  {
    int my_id;
    int parent_id;
    int from;
    int to;
    unsigned long ground;
    unsigned long excited;
    /** index of the child nodes
     */
    std::vector<int> children;
    /** default constructor
     * Invalidate the excitation
     */
    inline explicit ci_node_proxy(unsigned long vid=0u, unsigned long cid=0u)
      : my_id(0), parent_id(0),from(0),to(0), ground(vid),excited(cid) 
    { }

    inline void add_child(int i)
    {
      children.push_back(i);
    }

    inline bool closed() const
    {
      return children.empty();
    }

    template<unsigned CMAX>
      inline void write_node(ostream& os, int level, int& count, std::vector<ci_node_proxy>& ci)
      {
        my_id=count;
        //if(children.size()) 
        os  << "<node "
          //<< "<node level=\""<<level
          //<< "\" g=\"" << std::bitset<CMAX>(ground) 
          //<< "\" e=\"" << std::bitset<CMAX>(excited) 
          //<< "\" g_id=\"" << ground
          //<< "\" e_id=\"" << excited
          //<< "\" from=\"" << from 
          << " from=\"" << from 
          << "\" to=\""<< to 
          <<"\" p_id=\"" << parent_id
          <<"\" my_id=\"" << my_id //count
          << "\"";
        if(children.size()) os << ">" << std::endl;
        else os << "/>" << std::endl;
        count++;
        for(int i=0; i<children.size(); ++i) 
        {
          int next=level+1;
          //ci[ children[i] ].write_node<CMAX>(os,next,children[i],ci);
          //reassign the id
          ci[ children[i] ].parent_id=my_id;
          ci[ children[i] ].write_node<CMAX>(os,next,count,ci);
        }
        if(children.size()) os << "</node>"<<std::endl;
      }


    template<typename NT>
    void build_tree(int m, int& count, std::vector<ci_node_proxy>& ci, NT* my)
    {
      my->my_id=my_id;
      my->parent_id=parent_id;
      my->from=m-1-from;
      my->to=m+to;
      count++;
      if(children.size())
      {
        my->children.resize(children.size(),0);
        for(int i=0; i<children.size(); ++i)
        {
          NT* anode=new NT;
          ci[children[i]].build_tree(m,count,ci,anode);
          my->children[i]=anode;
        }
      }
    }
  }; 

  struct ci_node
  {
    int my_id;
    int parent_id;
    int from;
    int to;
    vector<ci_node*> children;

    inline ci_node():my_id(0),parent_id(0),from(0),to(0){}

    ///copy constructor
    inline ci_node(const ci_node& rhs)
      :my_id(rhs.my_id),parent_id(rhs.parent_id),from(rhs.from),to(rhs.to)
    {
      children.resize(rhs.children.size(),0);
      for(int i=0; i<children.size();++i)
        children[i]=new ci_node(*rhs.children[i]);
    }

    ///destructor
    ~ci_node()
    {
      cerr << "Deleting node = " << my_id << endl;
      for(int i=0; i<children.size(); ++i) delete children[i];
    }

    void build_tree(int m, std::vector<ci_node_proxy>& ci)
    {
      int count=0;
      ci[0].build_tree(m,count,ci,this);
    }

    //this is just debug
    void write_node(ostream& os)
    {
      os  << "<node "
        << " from=\"" << from 
        << "\" to=\""<< to 
        <<"\" p_id=\"" << parent_id
        <<"\" my_id=\"" << my_id //count
        << "\"";
      if(children.size()) os << ">" << std::endl;
      else os << "/>" << std::endl;
      for(int i=0; i<children.size(); ++i) children[i]->write_node(os);
      if(children.size()) os << "</node>"<<std::endl;
    }

  };
}
#endif
