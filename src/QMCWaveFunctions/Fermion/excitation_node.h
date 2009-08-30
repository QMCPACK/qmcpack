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
#include <OhmmsPETE/OhmmsMatrix.h>
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
  /** a node of multideterminant tree
  */
  template<typename T>
  struct excitation_node
  {
    typedef Matrix<T> matrix_type;
    int my_id;
    int parent_id;
    int from;
    int to;
    ///ratio with respect to the parent
    T ratio;
    unsigned long ground;
    unsigned long excited;
    matrix_type inverse;
    std::vector<int> children;

    /** default constructor
     * Invalidate the excitation
     */
    inline explicit excitation_node(unsigned long vid=0u, unsigned long cid=0u)
      : my_id(0), parent_id(0),from(0),to(0),ratio(1.0),ground(vid),excited(cid) {}

    inline void resize(int m)
    {
      inverse.resize(m,m);
      //if(children.size()) inverse.resize(m,m);
    }

    inline void add_child(int i)
    {
      children.push_back(i);
    }

    inline bool closed() const
    {
      return children.empty();
    }

    inline void get_ratios_root(const matrix_type& psi_big
        , std::vector<excitation_node>& ci
        , std::vector<T>& ratios
        , int m)
    {
      ratios[0]=1.0;
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].get_ratios(psi_big,ci,ratios,m,1.0);
    }

    inline void get_ratios(const matrix_type& psi_big
        , std::vector<excitation_node>& ci
        , std::vector<T>& ratios
        , int m
        , T ratio_base)
    {
      int v=m-1-from;
      int c=m+to;
      cout << "my_id =" << my_id << " " << " ratio by promotion from="<< v << " to=" << c << endl;
      ratio=BLAS::dot(m,ci[parent_id].inverse[v],psi_big[c]);
      ratios[my_id]=ratio_base*ratio;
      //T r=ratios[my_id]=BLAS::dot(m,ci[parent_id].inverse.data()+v,m,psi_big[c],1);
      //if(children.size()) 
      //{
        inverse=ci[parent_id].inverse;
        det_row_update(inverse.data(),psi_big[c],m,v,ratio);
      //}

        ratio_base *= ratio;
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].get_ratios(psi_big,ci,ratios,m,ratio_base);
    }

    template<unsigned CMAX>
      inline void write_node(int level, int count, std::vector<excitation_node>& ci)
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
