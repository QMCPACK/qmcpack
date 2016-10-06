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
 *
 * column and row denote those of the original matrix.
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
  /** inverse matrix
   *
   * The size can be reduced to (vmax,m)
   * where vmax is the occupied states for the excitations.
   */
  matrix_type inverse;
  /** index of the child nodes
   */
  std::vector<int> children;
  /** index of the rows to be updated by this node
   */
  std::vector<int> peers;

  /** default constructor
   * Invalidate the excitation
   */
  inline explicit excitation_node(unsigned long vid=0u, unsigned long cid=0u)
    : my_id(0), parent_id(0),from(0),to(0),ratio(1.0),utv(0.0), ground(vid),excited(cid)
  { }

  //inline void resize(int m)
  //{
  //  inverse.resize(m,m);
  //}

  /** reset from and to with respect to m ground states
   * @param m number of ground states
   * @param nv number of ground states that are excited
   * @param ci CI expansion
   */
  void resize(int m, int nv, std::vector<excitation_node>& ci)
  {
    inverse.resize(m,m);
    peers.reserve(nv);
    for(int i=0, iv=m-1; i<nv; ++i, --iv)
      peers.push_back(iv);
    for(int i=0; i<children.size(); ++i)
      ci[ children[i] ].resetOffset(m,peers,ci);
  }

  void resetOffset(int m, std::vector<int>& p_peers, std::vector<excitation_node>& ci)
  {
    inverse.resize(m,m);
    if(children.empty())
      return;
    //build-up the rows this node has to update during recursion
    peers.reserve(p_peers.capacity());
    int v=m-1-from;
    for(int i=0, iv=0; i<p_peers.size(); ++i)
      if(p_peers[i]!=v)
        peers.push_back(p_peers[i]);
    for(int i=0; i<children.size(); ++i)
      ci[ children[i] ].resetOffset(m,peers,ci);
  }

  inline void add_child(int i)
  {
    children.push_back(i);
  }

  inline bool closed() const
  {
    return children.empty();
  }

  /** evaluate the determinants recursively starting from the ground state
   * @param psi_big SPO set including the ground and excited states
   * @param ci vector of excitation_node
   * @param ratios ratios[i]=det(i)/det(0)
   * @param promote_row if true, the row index of the original matrix denotes the state index
   */
  inline void getRatios(const matrix_type& psi_big, std::vector<excitation_node>& ci
                        , std::vector<T>& ratios, bool promote_row)
  {
    ratios[0]=1.0;
    const int m=inverse.rows();
    if(promote_row)
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].getRatiosByRowPromotion(inverse,psi_big,ci,ratios,m,1.0);
    else
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].getRatiosByColPromotion(inverse,psi_big,ci,ratios,m,1.0);
  }

  inline void getRatiosByColPromotion(const matrix_type& inv0
                                      , const matrix_type& psi_big
                                      , std::vector<excitation_node>& ci
                                      , std::vector<T>& ratios
                                      , int m
                                      , T ratio_base
                                     )
  {
    inverse=inv0;
    const int v=m-1-from;
    const T* restrict u_old=inv0[v];
    const T* restrict u_new=psi_big[m+to];
    utv=BLAS::dot(m,u_old,u_new);
    BLAS::axpy(m,utv,u_old,inverse[v]);
    //utv=BLAS::dot(m,inv0[v],psi_big[c]);
    //BLAS::axpy(m,utv,inv0[v],inverse[v]);
    T inv_utv=1.0/utv;
    for(int iv=0; iv<peers.size(); ++iv)
    {
      const int vv=peers[iv];
      //T  gamma=-inv_utv*BLAS::dot(m,inv0[vv],psi_big[c]);
      //BLAS::axpy(m,gamma,inv0[v],inverse[vv]);
      T  gamma=-inv_utv*BLAS::dot(m,inv0[vv],u_new);
      BLAS::axpy(m,gamma,u_old,inverse[vv]);
    }
    ratios[my_id]=ratio_base*utv;
    for(int i=0; i<children.size(); ++i)
      ci[children[i]].getRatiosByColPromotion(inverse,psi_big,ci,ratios,m,ratios[my_id]);
  }

  inline void getRatiosByRowPromotion(const matrix_type& inv0
                                      , const matrix_type& psi_big
                                      , std::vector<excitation_node>& ci
                                      , std::vector<T>& ratios
                                      , int m
                                      , T ratio_base
                                     )
  {
    inverse=inv0;
    const int v=m-1-from;
    const int c=m+to;
    utv=BLAS::dot(m,inv0.data()+v,m,psi_big[c],1);
    BLAS::axpy(m,utv,inv0.data()+v,m,inverse.data()+v,m);
    T inv_utv=1.0/utv;
    for(int iv=0; iv<peers.size(); ++iv)
    {
      const int vv=peers[iv];
      T  gamma=-inv_utv*BLAS::dot(m,inv0.data()+vv,m,psi_big[c],1);
      BLAS::axpy(m,gamma,inv0.data()+v,m,inverse.data()+vv,m);
    }
    ratios[my_id]=ratio_base*utv;
    for(int i=0; i<children.size(); ++i)
      ci[children[i]].getRatiosByRowPromotion(inverse,psi_big,ci,ratios,m,ratios[my_id]);
  }

  inline void inverseUpdate(const matrix_type& psi_big, std::vector<excitation_node>& ci, std::vector<T>& ratios
                            , bool row_excitation)
  {
    const int m=inverse.rows();
    ratios[0]=1.0;
    if(row_excitation)
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].inverseUpdateByRow(psi_big,ci,ratios,m,1.0);
    else
      for(int i=0; i<children.size(); ++i)
        ci[ children[i] ].inverseUpdateByColumn(psi_big,ci,ratios,m,1.0);
  }

  void inverseUpdateByColumn(const matrix_type& psi_big , std::vector<excitation_node>& ci , std::vector<T>& ratios
                             , int m, T ratio_base)
  {
    const int v=m-1-from;
    const int c=m+to;
    ratio=BLAS::dot(m,ci[parent_id].inverse[v],psi_big[c]);
    ratios[my_id]=ratio_base*ratio;
    //if(children.size())
    //{
    inverse=ci[parent_id].inverse;
    det_row_update(inverse.data(),psi_big[c],m,v,ratio);
    ratio_base *= ratio;
    for(int i=0; i<children.size(); ++i)
      ci[ children[i] ].inverseUpdateByColumn(psi_big,ci,ratios,m,ratio_base);
    //}
  }

  void inverseUpdateByRow(const matrix_type& psi_big , std::vector<excitation_node>& ci , std::vector<T>& ratios
                          , int m, T ratio_base)
  {
    int v=m-1-from;
    int c=m+to;
    ratio=BLAS::dot(m,ci[parent_id].inverse.data()+v,m,psi_big[c],1);
    ratios[my_id]=ratio_base*ratio;
    //if(children.size())
    //{
    inverse=ci[parent_id].inverse;
    det_col_update(inverse.data(),psi_big[c],m,v,ratio);
    ratio_base *= ratio;
    for(int i=0; i<children.size(); ++i)
      ci[ children[i] ].inverseUpdateByRow(psi_big,ci,ratios,m,ratio_base);
    //}
  }

  /** evaluate the ratio with a new column
   *
   * inverse matrix is updated by inverseUpdate(...,true)
   */
  inline void getRatioByColSubstitution(std::vector<T>& u_c
                                        , std::vector<excitation_node>& ci
                                        , std::vector<T>& ratios
                                        , int col_id
                                       )
  {
    const int m=inverse.rows();
    ratios[0]=ratio=BLAS::dot(m,inverse[col_id],u_c.data());
    std::vector<T> temp(u_c);
    for(int i=0; i<children.size(); ++i)
    {
      ci[ children[i] ].getRatioByColSubstitution(temp,ci,ratios,m,col_id);
      //restore the original u_c
      for(int iv=0,vv=m-1; iv<peers.size(); ++iv,--vv)
        temp[vv]=u_c[vv];
    }
  }

  void getRatioByColSubstitution(std::vector<T>& u_c , std::vector<excitation_node>& ci
                                 , std::vector<T>& ratios , int m, int col_id)
  {
    int v=m-1-from;
    int c=m+to;
    u_c[v]=u_c[c];
    ratios[my_id]=ratio=BLAS::dot(m,inverse[col_id],u_c.data());
    std::vector<T> temp(u_c);
    for(int i=0; i<children.size(); ++i)
    {
      ci[ children[i] ].getRatioByColSubstitution(temp,ci,ratios,m,col_id);
      for(int iv=0; iv<peers.size(); ++iv)
      {
        int vv=peers[iv];
        temp[vv]=u_c[vv];
      }
    }
  }

  /** evaluate the ratio with a new row
   *
   * inverse matrix is updated by inverseUpdateByColumn
   */
  inline void getRatioByRowSubstitution(std::vector<T>& u_c
                                        , std::vector<excitation_node>& ci
                                        , std::vector<T>& ratios
                                        , int row_id
                                       )
  {
    const int m=inverse.rows();
    ratios[0]=ratio=BLAS::dot(m,inverse.data()+row_id,m,u_c.data(),1);
    std::vector<T> temp(u_c);
    for(int i=0; i<children.size(); ++i)
    {
      ci[ children[i] ].getRatioByRowSubstitution(temp,ci,ratios,m,row_id);
      //restore the original u_c
      for(int iv=0,vv=m-1; iv<peers.size(); ++iv,--vv)
        temp[vv]=u_c[vv];
    }
  }

  void getRatioByRowSubstitution(std::vector<T>& u_c , std::vector<excitation_node>& ci
                                 , std::vector<T>& ratios , int m, int row_id)
  {
    const int v=m-1-from;
    const int c=m+to;
    u_c[v]=u_c[c];
    ratios[my_id]=ratio=BLAS::dot(m,inverse.data()+row_id,m,u_c.data(),1);
    std::vector<T> temp(u_c);
    for(int i=0; i<children.size(); ++i)
    {
      ci[ children[i] ].getRatioByRowSubstitution(temp,ci,ratios,m,row_id);
      for(int iv=0; iv<peers.size(); ++iv)
      {
        int vv=peers[iv];
        temp[vv]=u_c[vv];
      }
    }
  }


  template<unsigned CMAX>
  inline void write_node(std::ostream& os, int level, int& count, std::vector<excitation_node>& ci)
  {
    my_id=count;
    //if(children.size())
    os << "<node level=\""<<level
       << "\" g=\"" << std::bitset<CMAX>(ground)
       << "\" e=\"" << std::bitset<CMAX>(excited)
       << "\" g_id=\"" << ground
       << "\" e_id=\"" << excited
       << "\" from=\"" << from
       << "\" to=\""<< to
       <<"\" p_id=\"" << parent_id
       <<"\" my_id=\"" << my_id //count
       << "\"";
    if(children.size())
      os << ">" << std::endl;
    else
      os << "/>" << std::endl;
    count++;
    for(int i=0; i<children.size(); ++i)
    {
      int next=level+1;
      //ci[ children[i] ].write_node<CMAX>(os,next,children[i],ci);
      //reassign the id
      ci[ children[i] ].parent_id=my_id;
      ci[ children[i] ].write_node<CMAX>(os,next,count,ci);
    }
    if(children.size())
      os << "</node>"<<std::endl;
  }
};
}
#endif
