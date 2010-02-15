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
#include <Utilities/IteratorUtility.h>
#include <type_traits/scalar_traits.h>

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
      my->to=to;
      //my->to=m+to;
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

  /** CI node with a recursive tree
   *
   * Each node handles an excitation from an occupied state to an excited state.
   */
  template<typename T>
  struct ci_node
  {
    typedef typename scalar_traits<T>::real_type real_type;
    typedef typename scalar_traits<T>::value_type value_type;
    typedef Matrix<T> matrix_type;
    typedef Vector<T> vector_type;
    /// seiralized ID  of this node
    int my_id;
    /**  seiralized ID of the parent
     *
     * This is not needed in a tree presentation
     */
    int parent_id;
    /** the index of the occupied state to be replaced
     */
    int from;
    /** the index of the unoccupied state to replace from row/column
     */
    int to;
    /** inverse matrix
     */
    matrix_type inverse;
    ///the occupied states that should be updated
    std::vector<int> peers;
    ///child nodes
    std::vector<ci_node*> children;

    ///default constructor
    inline ci_node():my_id(0),parent_id(0),from(0),to(0){}

    ///copy constructor
    inline ci_node(const ci_node& rhs)
      :my_id(rhs.my_id),parent_id(rhs.parent_id)
       ,from(rhs.from),to(rhs.to),peers(rhs.peers)
    {
      if(rhs.children.size())
      {
        children.resize(rhs.children.size(),0);
        for(int i=0; i<children.size();++i)
          children[i]=new ci_node(*rhs.children[i]);
      }
    }

    ///destructor
    ~ci_node()
    {
      delete_iter(children.begin(),children.end());
    }

    /** set the peers for recursive updates 
     * @param m number of occupied states
     * @param nv number of occupied states to be promoted
     *
     * This is called by the root node to build up peers for recursive updates
     */
    void set_peers(int m, int nv)
    {
      inverse.resize(m,m);
      peers.reserve(m);
      for(int i=m-nv;i<m; ++i) peers.push_back(i);
      for(int i=0; i<children.size(); ++i) children[i]->set_peers(m,peers);
    }

    /** set the peers from the peers of its parent
     * @param p_peers list of peers of the parent
     */
    void set_peers(int m, const std::vector<int>& p_peers)
    {
      inverse.resize(m,m);
      if(children.size())
      {
        peers.reserve(p_peers.capacity());
        for(int i=0; i<p_peers.size(); ++i) if(p_peers[i]!=from) peers.push_back(p_peers[i]);
        for(int i=0; i<children.size(); ++i) children[i]->set_peers(m,peers);
      }
    }

    void build_tree(int m, std::vector<ci_node_proxy>& ci)
    {
      int count=0;
      ci[0].build_tree(m,count,ci,this);
    }

    /** get ratios with respect to \f$D_0\f$
     * @param psiv container of the valence states
     * @param psic container of the conduction states
     * @param ratios ratios[my_id] is evaluated
     * @param replace_row true, if rows are replaced
     */
    T getRatios(const matrix_type& psiv, const matrix_type& psic, std::vector<T>& ratios, bool replace_row)
    {
      const int m=inverse.rows();
      inverse=psiv;
      T det_0=invert_matrix(inverse,true);
      T ratio_base=ratios[0]=1.0;
      if(replace_row)
        for(int i=0; i<children.size();++i)
          children[i]->getRatiosByRowPromotion(inverse,psic,m,ratios,ratio_base);
      else
        for(int i=0; i<children.size();++i)
          children[i]->getRatiosByColPromotion(inverse,psic,m,ratios,ratio_base);
      return det_0;
    }

    inline void getRatiosByRowPromotion(const matrix_type& inv0
        , const matrix_type& psic
        , int m
        , std::vector<T>& ratios
        , T ratio_base
        )
    {
      T utv=BLAS::dot(m,inv0.data()+from,m,psic[to],1);
      ratios[my_id]=ratio_base*utv;

      if(children.size())
      {
        ratio_base=ratios[my_id];

        inverse=inv0;
        //fully update inverse and propagate it
        //det_col_update(inverse.data(),psi_big[to],m,from,utv);
        //only partial update
        BLAS::axpy(m,utv,inv0.data()+from,m,inverse.data()+from,m);
        T inv_utv=1.0/utv;
        for(int iv=0; iv<peers.size(); ++iv)
        {
          const int vv=peers[iv];
          T  gamma=-inv_utv*BLAS::dot(m,inv0.data()+vv,m,psic[to],1);
          BLAS::axpy(m,gamma,inv0.data()+from,m,inverse.data()+vv,m);
        }

        for(int i=0; i<children.size(); ++i)
          children[i]->getRatiosByRowPromotion(inverse,psic,m,ratios,ratio_base);
      }
    }

    inline void getRatiosByColPromotion(const matrix_type& inv0
        , const matrix_type& psic
        , int m
        , std::vector<T>& ratios
        , T ratio_base
        )
    {
      T utv=BLAS::dot(m,inv0[from],psic[to]);
      ratios[my_id]=ratio_base*utv;
      if(children.size())
      {
        ratio_base=ratios[my_id];
        inverse=inv0;

        //fully update the inverse and propagate it
        //det_row_update(inverse.data(),psi_big[to],m,from,utv);
        //partial update
        BLAS::axpy(m,utv,inv0[from],inverse[from]);
        T inv_utv=1.0/utv;
        for(int iv=0; iv<peers.size(); ++iv)
        {
          const int vv=peers[iv];
          T  gamma=-inv_utv*BLAS::dot(m,inv0[vv],psic[to]);
          BLAS::axpy(m,gamma,inv0[from],inverse[vv]);
        }

        for(int i=0; i<children.size(); ++i)
          children[i]->getRatiosByColPromotion(inverse,psic,m,ratios,ratio_base);
      }
    }

    void debugRatios(const matrix_type& psi0, const matrix_type& psic,std::vector<T>& dets
        , bool replace_row)
    {
      const int m=inverse.rows();
      inverse=psi0;
      dets[my_id]=invert_matrix(inverse,true);
      if(replace_row)
        for(int i=0; i<children.size();++i) children[i]->debugRatiosByRowPromotion(psi0,psic,m,dets);
      else
        for(int i=0; i<children.size();++i) children[i]->debugRatiosByColPromotion(psi0,psic,m,dets);
    }

    void debugRatiosByRowPromotion(const matrix_type& psi0, const matrix_type& psic,int m, std::vector<T>& dets)
    {
      inverse=psi0;
      for(int j=0; j<m;++j) inverse(from,j)=psic(to,j);
      matrix_type psip(inverse);
      dets[my_id]=invert_matrix(inverse,true);
      for(int i=0; i<children.size();++i)
        children[i]->debugRatiosByRowPromotion(psip,psic,m,dets);
    }

    void debugRatiosByColPromotion(const matrix_type& psi0, const matrix_type& psic,int m, std::vector<T>& dets)
    {
      inverse=psi0;
      for(int j=0; j<m;++j) inverse(j,from)=psic(to,j);
      matrix_type psip(inverse);
      dets[my_id]=invert_matrix(inverse,true);
      for(int i=0; i<children.size();++i)
        children[i]->debugRatiosByColPromotion(psip,psic,m,dets);
    }

    void ratioByRowSubstitution(const vector_type& psiv_big, int irow, std::vector<T>& ratios)
    {
      const int m=inverse.rows();
      ratios[my_id]=BLAS::dot(m,inverse[irow],psiv_big.data());

      //split a big vector into valance and conduction vectors
      vector_type psiv(m);
      vector_type psic(psiv_big.size()-m);
      BLAS::copy(psic.size(),psiv_big.data()+m,psic.data());
      for(int i=0; i<children.size();++i)
      {
        std::copy(psiv_big.data(),psiv_big.data()+m,psiv.data());
        children[i]->ratioByRowSubstitution(psiv, psiv_big, irow, m, ratios);
      }
    }

    void ratioByRowSubstitution(vector_type& psiv, const vector_type& psic, int irow, int m, std::vector<T>& ratios)
    {
      psiv(from)=psic(to);
      ratios[my_id]=BLAS::dot(m,inverse[irow],psiv.data());
      if(children.empty()) return;
      vector_type psiv_temp(m);
      for(int i=0; i<children.size();++i)
      {
        std::copy(psiv.begin(),psiv.end(),psiv_temp.data());
        children[i]->ratioByRowSubstitution(psiv_temp, psic, irow, m, ratios);
      }
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
      if(peers.size())
      {
        os << " peers=\"";
        int i=0;
        for(i=0; i<peers.size()-1;i++) os << peers[i] << " ";
        os << peers[i] << "\"";
      }
      if(children.size()) os << ">" << std::endl;
      else os << "/>" << std::endl;
      for(int i=0; i<children.size(); ++i) children[i]->write_node(os);
      if(children.size()) os << "</node>"<<std::endl;
    }

  };
}
#endif
