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
#ifndef QMCPLUSPLUS_EXCITATION_NODE_H
#define QMCPLUSPLUS_EXCITATION_NODE_H
#include <bitset>
#include <vector>
#include <iostream>
namespace qmcplusplus
{
  /** a node of multideterminant tree
  */
  struct excitation_node
  {
    unsigned long ground;
    unsigned long excited;
    int ref_state;
    int from;
    int to;
    std::vector<int> children;

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
      inline void write_node(int & level, int count, std::vector<excitation_node>& ci)
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
          ci[ children[i] ].write_node<CMAX>(next,children[i],ci);
        }
        if(children.size()) std::cout << "</node>"<<std::endl;
      }
  }; 
}
#endif
