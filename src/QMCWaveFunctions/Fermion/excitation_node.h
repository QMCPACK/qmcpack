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
#include <Numerics/MatrixOperators.h>
#include <Numerics/determinant_operators.h>

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
    ///product
    T utv;
    unsigned long ground;
    unsigned long excited;
    matrix_type inverse;
    std::vector<int> children;

    /** default constructor
     * Invalidate the excitation
     */
    inline explicit excitation_node(unsigned long vid=0u, unsigned long cid=0u)
      : my_id(0), parent_id(0),from(0),to(0),ratio(1.0),utv(0.0), ground(vid),excited(cid) {}

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

    inline void get_ratios_root(const matrix_type& psi_big, std::vector<excitation_node>& ci
        , std::vector<T>& ratios, int vmin)
    {
      ratios[0]=1.0;
      const int m=inverse.rows();
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].get_ratios_partial(inverse,psi_big,ci,ratios,m,vmin,1.0,1.0);
    }

    inline void get_ratios_partial(const matrix_type& inv0
        , const matrix_type& psi_big
        , std::vector<excitation_node>& ci
        , std::vector<T>& ratios
        , int m
        , int vmin
        , T ratio_base
        , T inv_utv
        )
    {
      const int n=m-vmin;
      T temp[n];

      BLAS::copy((m-vmin)*m,inv0[vmin],inverse[vmin]);
      //update inverse matrix partially
      //do all the blocks to be updated: this can be further optimized by bit masks
      if(parent_id>0)
      {
        int pv=m-1-ci[parent_id].from;
        int pc=m+ci[parent_id].to;
        for(int vv=vmin,iv=0; vv<m; ++vv,++iv)
        {
          T  gamma=BLAS::dot(m,inv0[vv],psi_big[pc])*inv_utv;
          for(int j=0; j<m; ++j) inverse(vv,j) -= gamma*inv0(pv,j);
        }
      }

      int v=m-1-from;
      int c=m+to;
      utv=BLAS::dot(m,inverse[v],psi_big[c]);
      ratios[my_id]=ratio_base*utv;

      for(int i=0; i<children.size(); ++i)
        ci[children[i]].get_ratios_partial(inverse,psi_big,ci,ratios,m,vmin,ratios[my_id],1.0/utv);
    }

    inline void get_ratios_root_debug(const matrix_type& psi_big
        , std::vector<excitation_node>& ci
        , std::vector<T>& ratios
        )
    {
      const int m=inverse.rows();
      ratios[0]=1.0;
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].get_ratios_debug(psi_big,ci,ratios,m,1.0);
    }

    void get_ratios_debug(const matrix_type& psi_big , std::vector<excitation_node>& ci , std::vector<T>& ratios
        , int m, T ratio_base)
    {
      int v=m-1-from;
      int c=m+to;
      ratio=BLAS::dot(m,ci[parent_id].inverse[v],psi_big[c]);
      ratios[my_id]=ratio_base*ratio;

      //if(children.size())
      //{
        inverse=ci[parent_id].inverse;
        det_row_update(inverse.data(),psi_big[c],m,v,ratio);

        ratio_base *= ratio;
        for(int i=0; i<children.size(); ++i)
          ci[ children[i] ].get_ratios_debug(psi_big,ci,ratios,m,ratio_base);
      //}
    }

    template<unsigned CMAX>
      inline void write_node(int level, int count, std::vector<excitation_node>& ci)
      {
        //if(children.size()) 
        std::cout << "<node level=\""<<level
          << "\" g=\"" << std::bitset<CMAX>(ground) 
          << "\" e=\"" << std::bitset<CMAX>(excited) 
          << "\" g_id=\"" << ground
          << "\" e_id=\"" << excited
          << "\" from=\"" << from 
          << "\" to=\""<< to 
          <<"\" my_id=\"" << count
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
