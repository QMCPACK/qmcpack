//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_CI_NODE_H
#define QMCPLUSPLUS_CI_NODE_H
#include <bitset>
#include <vector>
#include <algorithm>
#include <iostream>
#include <Utilities/IteratorUtility.h>
#include <type_traits/scalar_traits.h>
#include "Numerics/OhmmsBlas.h"
#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/Tensor.h"
#include "Numerics/DeterminantOperators.h"
#include "QMCWaveFunctions/Fermion/ci_configuration.h"

namespace qmcplusplus
{

/** CI node with a recursive tree
 *
 * Each node handles an excitation from an occupied state to an excited state.
 */
template<typename T>
struct ci_node
{
  typedef typename scalar_traits<T>::real_type real_type;
  typedef typename scalar_traits<T>::value_type value_type;
  typedef TinyVector<value_type,OHMMS_DIM>  GradType;
  typedef TinyVector<real_type,OHMMS_DIM>   PosType;
  typedef Tensor<value_type,OHMMS_DIM>      HessType;
  typedef Vector<value_type>        ValueVector;
  typedef Matrix<value_type>        ValueMatrix;
  typedef Vector<GradType>          GradVector;
  typedef Matrix<GradType>          GradMatrix;
  typedef Vector<HessType>          HessVector;
  typedef Matrix<HessType>          HessMatrix;
  typedef TinyVector<HessType, OHMMS_DIM>   GGGType;
  typedef Vector<GGGType>           GGGVector;
  typedef Matrix<GGGType>           GGGMatrix;

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
  ValueMatrix inverse;
  ///the occupied states that should be updated
  std::vector<int> peers;
  ///child nodes
  std::vector<ci_node*> children;
  /// the configuration of the node
  configuration myC;
  /// temporary storage
  GradMatrix dpsiM;
  GradVector dpsiV;
  ValueMatrix d2psiM;
  ValueVector d2psiV;

  ///default constructor
  inline ci_node(int m, int n):my_id(0),parent_id(0),from(0),to(0)
  {
    inverse.resize(m,m);
    peers.reserve(m);
    dpsiM.resize(m,n);
    d2psiM.resize(m,n);
    dpsiV.resize(n);
    d2psiV.resize(n);
  }

  ///copy constructor
  inline ci_node(const ci_node& rhs)
    :my_id(rhs.my_id),parent_id(rhs.parent_id)
    ,from(rhs.from),to(rhs.to),peers(rhs.peers),
    myC(rhs.myC)
  {
    int m=rhs.inverse.cols();
    int n=rhs.dpsiM.cols();
    inverse.resize(m,m);
    peers.reserve(m);
    dpsiM.resize(m,n);
    d2psiM.resize(m,n);
    dpsiV.resize(n);
    d2psiV.resize(n);
    if(rhs.children.size())
    {
      children.resize(rhs.children.size(),0);
      for(int i=0; i<children.size(); ++i)
        children[i]=new ci_node(*rhs.children[i]);
    }
  }

  ///destructor
  ~ci_node()
  {
    delete_iter(children.begin(),children.end());
  }

  // assumes that tree is fully connected. Generalize later.
  void build_tree(std::vector<configuration>& confg, configuration& baseC, std::vector<int>& map2det)
  {
    int count=1,rem,add;
    int m=inverse.cols();
    int n=dpsiM.cols();
    for(std::vector<configuration>::iterator it=confg.begin()+1; it!=confg.end(); ++it)
      it->taken=false;
    for(int i=0; i<confg.size(); i++)
    {
      if(baseC.isSingle(confg[i],rem,add) && !(confg[i].taken) )
      {
        confg[i].taken=true;
        map2det[i] = count;
        children.push_back(new ci_node(m,n,count++,0,rem,add,confg[i]));
        children.back()->look_for_peers(confg,count,map2det);
        children.back()->set_peers(peers);
        // peers needs to be ordered
        sort (peers.begin(), peers.end());
      }
    }
    for(std::vector<configuration>::iterator it=confg.begin()+1; it!=confg.end(); ++it)
    {
      if(!(it->taken) )
      {
        APP_ABORT("Error in ci_node::build_tree. Either tree is not fully connected or error in algorithm.\n");
      }
    }
  }

  /** get ratios with respect to \f$D_0\f$
   * @param psiv container of the valence states
   * @param psic container of the conduction states
   * @param ratios ratios[my_id] is evaluated
   * @param replace_row true, if rows are replaced
   */
  /*
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
  */

  void getRatios(const ValueMatrix& inv0, const ValueMatrix& psi, const GradMatrix& dpsi, const ValueMatrix& d2psi, ValueVector& ratios, GradMatrix& grads, ValueMatrix& lapls, bool replace_row=false)
  {
    const int m=inverse.rows();
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        inverse(i,j)=inv0(j,i);
    T ratio_base=ratios[0]=1.0;
    if(replace_row)
    {
      APP_ABORT("Need to implement getRatiosByRowPromotion. \n");
      //for(int i=0; i<children.size();++i)
      //  children[i]->getRatiosByRowPromotion(inverse,psic,m,ratios,ratio_base);
    }
    else
    {
      for(int i=0; i<children.size(); ++i)
        children[i]->getRatiosByColPromotion(inverse,psi,dpsi,d2psi,ratios,grads,lapls,ratio_base);
    }
  }

  void getRatios(const int iat, const ValueMatrix& inv0, const ValueMatrix& psi, const GradVector& dpsi, const ValueVector& d2psi, ValueVector& ratios, GradMatrix& grads, ValueMatrix& lapls, bool replace_row=false)
  {
    const int m=inverse.rows();
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        inverse(i,j)=inv0(j,i);
    T ratio_base=ratios[0]=1.0;
    if(replace_row)
    {
      APP_ABORT("Need to implement getRatiosByRowPromotion. \n");
      //for(int i=0; i<children.size();++i)
      //  children[i]->getRatiosByRowPromotion(inverse,psic,m,ratios,ratio_base);
    }
    else
    {
      for(int i=0; i<children.size(); ++i)
        children[i]->getRatiosByColPromotion(iat,inverse,psi,dpsi,d2psi,ratios,grads,lapls,ratio_base);
    }
  }

  void getRatios(const ValueMatrix& inv0, const ValueMatrix& psi, ValueVector& ratios, bool replace_row=false)
  {
    const int m=inverse.rows();
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        inverse(i,j)=inv0(j,i);
    T ratio_base=ratios[0]=1.0;
    if(replace_row)
    {
      APP_ABORT("Need to implement getRatiosByRowPromotion. \n");
      //for(int i=0; i<children.size();++i)
      //  children[i]->getRatiosByRowPromotion(inverse,psic,m,ratios,ratio_base);
    }
    else
    {
      for(int i=0; i<children.size(); ++i)
        children[i]->getRatiosByColPromotion(inverse,psi,ratios,ratio_base);
    }
  }

  /*
      inline void getRatiosByRowPromotion(
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
  */

  inline void getRatiosByColPromotion(const ValueMatrix& inv0
                                      , const ValueMatrix& psi
                                      , ValueVector& ratios
                                      , T ratio_base)
  {
    const int m=inv0.rows();
    const int m1=psi.rows();
// To who may read this:
// is it worth transposing psiM,dpsiM and d2psiM in SPOSetProxy
// to avoid using incremented versions of dot?
    //T utv=BLAS::dot(m,inv0[from],1,psi.data()+to,m1);
    T utv=0.0;
    for(int i=0; i<m; i++)
      utv += inv0(from,i)*psi(i,to);
    ratios[my_id]=ratio_base*utv;
    if(children.size())
    {
      ratio_base=ratios[my_id];
      inverse=inv0;
      //fully update the inverse and propagate it
      //det_row_update(inverse.data(),psi_big[to],m,from,utv);
      //partial update
      T inv_utv=1.0/utv;
      //BLAS::axpy(m,inv_utv-1.0,inv0[from],inverse[from]);
      for(int iv=0; iv<peers.size(); ++iv)
      {
        const int vv=peers[iv];
        //T  gamma=-inv_utv*BLAS::dot(m,inv0[vv],1,psi.data()+to,m1);
        //BLAS::axpy(m,gamma,inv0[from],inverse[vv]);
        T gamma=0.0;
        for(int i=0; i<m; i++)
          gamma-=inv0(vv,i)*psi(i,to);
        gamma*=inv_utv;
        for(int i=0; i<m; i++)
          inverse(vv,i) += gamma*inv0(from,i);
      }
      for(int i=0; i<children.size(); ++i)
        children[i]->getRatiosByColPromotion(inverse,psi,ratios,ratio_base);
    }
  }

  inline void getRatiosByColPromotion(const int iat
                                      , const ValueMatrix& inv0
                                      , const ValueMatrix& psi
                                      , const GradVector& dpsi
                                      , const ValueVector& d2psi
                                      , ValueVector& ratios
                                      , GradMatrix& grads
                                      , ValueMatrix& lapls
                                      , T ratio_base)
  {
    const int m=inv0.rows();
    T utv=0.0;
    for(int i=0; i<m; i++)
      utv += inv0(from,i)*psi(i,to);
    ratios[my_id]=ratio_base*utv;
    T inv_utv = 1.0/utv;
    dpsiV=dpsi;
    d2psiV=d2psi;
    dpsiV[from]=dpsi[to];
    d2psiV[from]=d2psi[to];
    if(children.size())
    {
      ratio_base=ratios[my_id];
      inverse=inv0;
      int iv=0;
      for(int k=0; k<from; k++)
      {
        T gamma=0.0;
        for(int i=0; i<m; i++)
          gamma-=inv0(k,i)*psi(i,to);
        gamma*=inv_utv;
        //    if(iv<m && k==peers[iv]) {
        //      iv++;
        for(int i=0; i<m; i++)
          inverse(k,i) += gamma*inv0(from,i);
        //    } else {
        //      inverse(k,iat)+= gamma*inv0(from,iat);
        //    }
      }
      //  inverse(from,iat)=inv0(from,iat)*inv_utv;
      for(int i=0; i<m; i++)
        inverse(from,i)=inv0(from,i)*inv_utv;
      for(int k=from+1; k<m; k++)
      {
        T gamma=0.0;
        for(int i=0; i<m; i++)
          gamma-=inv0(k,i)*psi(i,to);
        gamma*=inv_utv;
        //    if(iv<m && k==peers[iv]) {
        //      iv++;
        for(int i=0; i<m; i++)
          inverse(k,i) += gamma*inv0(from,i);
        //    } else {
        //      inverse(k,iat)+= gamma*inv0(from,iat);
        //    }
      }
      grads(my_id,iat)=0.0;
      lapls(my_id,iat)=0.0;
      for(int j=0; j<m; j++)
        grads(my_id,iat)+=inverse(j,iat)*dpsiV[j];
      for(int j=0; j<m; j++)
        lapls(my_id,iat)+=inverse(j,iat)*d2psiV[j];
      for(int i=0; i<children.size(); ++i)
        children[i]->getRatiosByColPromotion(iat,inverse,psi,dpsiV,d2psiV,ratios,grads,lapls,ratio_base);
    }
    else
    {
      inverse=inv0;
      for(int k=0; k<from; k++)
      {
        T gamma=0.0;
        for(int i=0; i<m; i++)
          gamma-=inv0(k,i)*psi(i,to);
        gamma*=inv_utv;
        inverse(k,iat)+= gamma*inv0(from,iat);
      }
      inverse(from,iat)=inv0(from,iat)*inv_utv;
      for(int k=from+1; k<m; k++)
      {
        T gamma=0.0;
        for(int i=0; i<m; i++)
          gamma-=inv0(k,i)*psi(i,to);
        gamma*=inv_utv;
        inverse(k,iat)+= gamma*inv0(from,iat);
      }
      grads(my_id,iat)=0.0;
      lapls(my_id,iat)=0.0;
      for(int j=0; j<m; j++)
        grads(my_id,iat)+=inverse(j,iat)*dpsiV[j];
      for(int j=0; j<m; j++)
        lapls(my_id,iat)+=inverse(j,iat)*d2psiV[j];
    }
  }

  inline void getRatiosByColPromotion(const ValueMatrix& inv0
                                      , const ValueMatrix& psi
                                      , const GradMatrix& dpsi
                                      , const ValueMatrix& d2psi
                                      , ValueVector& ratios
                                      , GradMatrix& grads
                                      , ValueMatrix& lapls
                                      , T ratio_base)
  {
    const int m=inv0.rows();
    T utv=0.0;
    for(int i=0; i<m; i++)
      utv += inv0(from,i)*psi(i,to);
    ratios[my_id]=ratio_base*utv;
    T inv_utv = 1.0/utv;
    dpsiM=dpsi;
    d2psiM=d2psi;
    for(int j=0; j<m; ++j)
      dpsiM(j,from)=dpsi(j,to);
    for(int j=0; j<m; ++j)
      d2psiM(j,from)=d2psi(j,to);
    inverse=inv0;
    for(int k=0; k<from; k++)
    {
      T gamma=0.0;
      for(int i=0; i<m; i++)
        gamma-=inv0(k,i)*psi(i,to);
      gamma*=inv_utv;
      for(int i=0; i<m; i++)
        inverse(k,i) += gamma*inv0(from,i);
    }
    for(int i=0; i<m; i++)
      inverse(from,i)=inv0(from,i)*inv_utv;
    for(int k=from+1; k<m; k++)
    {
      T gamma=0.0;
      for(int i=0; i<m; i++)
        gamma-=inv0(k,i)*psi(i,to);
      gamma*=inv_utv;
      for(int i=0; i<m; i++)
        inverse(k,i) += gamma*inv0(from,i);
    }
    for(int i=0; i<m; i++)
    {
      grads(my_id,i)=0.0;
      lapls(my_id,i)=0.0;
      for(int j=0; j<m; j++)
        grads(my_id,i)+=inverse(j,i)*dpsiM(i,j);
      for(int j=0; j<m; j++)
        lapls(my_id,i)+=inverse(j,i)*d2psiM(i,j);
    }
    if(children.size())
    {
      ratio_base=ratios[my_id];
      for(int i=0; i<children.size(); ++i)
        children[i]->getRatiosByColPromotion(inverse,psi,dpsiM,d2psiM,ratios,grads,lapls,ratio_base);
    }
  }

  void debugRatios(const ValueMatrix& psi, ValueVector& dets
                   , bool replace_row=false)
  {
    const int m=inverse.rows();
    ValueMatrix psi0(m,m);
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        psi0(i,j) = psi(i,j);
    inverse=psi0;
    dets[my_id]=invert_matrix<ValueMatrix>(inverse,true);
    if(replace_row)
    {
      APP_ABORT("Need to implement getRatiosByRowPromotion. \n");
      //for(int i=0; i<children.size();++i) children[i]->debugRatiosByRowPromotion(psi0,psic,m,dets);
    }
    else
      for(int i=0; i<children.size(); ++i)
        children[i]->debugRatiosByColPromotion(psi0,psi,dets);
  }

  void debugRatios(const ValueMatrix& psi
                   , const GradMatrix& dpsi
                   , const ValueMatrix& d2psi
                   , ValueVector& ratios
                   , GradMatrix& grads
                   , ValueMatrix& lapls
                   , bool replace_row=false)
  {
    const int m=inverse.rows();
    ValueMatrix psi0(m,m);
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        psi0(i,j) = psi(i,j);
    inverse=psi0;
    ratios[my_id]=invert_matrix<ValueMatrix>(inverse,true);
    if(replace_row)
    {
      APP_ABORT("Need to implement getRatiosByRowPromotion. \n");
      //for(int i=0; i<children.size();++i) children[i]->debugRatiosByRowPromotion(psi0,psic,m,dets);
    }
    else
      for(int i=0; i<children.size(); ++i)
        children[i]->debugRatiosByColPromotion(psi0,psi,dpsi,d2psi,ratios,grads,lapls);
  }

  /*
      void debugRatiosByRowPromotion(const ValueMatrix& psi0, const ValueMatrix& psic,int m, ValueVector& dets)
      {
        inverse=psi0;
        for(int j=0; j<m;++j) inverse(from,j)=psic(to,j);
        matrix_type psip(inverse);
        dets[my_id]=invert_matrix(inverse,true);
        for(int i=0; i<children.size();++i)
          children[i]->debugRatiosByRowPromotion(psip,psic,m,dets);
      }
  */

  void debugRatiosByColPromotion(const ValueMatrix& psi0, const ValueMatrix& psic, ValueVector& dets)
  {
    const int m=inverse.rows();
    //inverse=psi0;
    int nq = 0;
    for(int i=0; i<myC.occup.size(); i++)
    {
      if(myC.occup[i])
      {
        for(int j=0; j<m; ++j)
          inverse(j,nq)=psic(j,i);
        nq++;
      }
    }
    ValueMatrix psip(inverse);
    dets[my_id]=invert_matrix<ValueMatrix>(inverse,true);
    for(int i=0; i<children.size(); ++i)
      children[i]->debugRatiosByColPromotion(psip,psic,dets);
  }

  void debugRatiosByColPromotion(const ValueMatrix& psi0
                                 , const ValueMatrix& psi
                                 , const GradMatrix& dpsi
                                 , const ValueMatrix& d2psi
                                 , ValueVector& ratios
                                 , GradMatrix& grads
                                 , ValueMatrix& lapls)
  {
    const int m=inverse.rows();
    //inverse=psi0;
    dpsiM=dpsi;
    d2psiM=d2psi;
    int nq = 0;
    for(int i=0; i<myC.occup.size(); i++)
    {
      if(myC.occup[i])
      {
        for(int j=0; j<m; ++j)
          inverse(j,nq)=psi(j,i);
        for(int j=0; j<m; ++j)
          dpsiM(j,nq)=dpsi(j,i);
        for(int j=0; j<m; ++j)
          d2psiM(j,nq)=d2psi(j,i);
        nq++;
      }
    }
    //for(int j=0; j<m;++j) inverse(j,from)=psi(j,to);
    //for(int j=0; j<m;++j) dpsiM(j,from)=dpsi(j,to);
    //for(int j=0; j<m;++j) d2psiM(j,from)=d2psi(j,to);
    ValueMatrix psip(inverse);
    ratios[my_id]=invert_matrix<ValueMatrix>(inverse,true);
    for(int i=0; i<m; i++)
    {
      grads(my_id,i)=0.0;
      lapls(my_id,i)=0.0;
      for(int j=0; j<m; j++)
        grads(my_id,i)+=inverse(j,i)*dpsiM(i,j);
      for(int j=0; j<m; j++)
        lapls(my_id,i)+=inverse(j,i)*d2psiM(i,j);
    }
    for(int i=0; i<children.size(); ++i)
      children[i]->debugRatiosByColPromotion(psip,psi,dpsiM,d2psiM,ratios,grads,lapls);
  }

  /*
      void ratioByRowSubstitution(const ValueVector& psiv_big, int irow, ValueVector& ratios)
      {
        const int m=inverse.rows();
        ratios[my_id]=BLAS::dot(m,inverse[irow],psiv_big.data());

        //split a big vector into valance and conduction vectors
        ValueVector psiv(m);
        ValueVector psic(psiv_big.size()-m);
        BLAS::copy(psic.size(),psiv_big.data()+m,psic.data());
        for(int i=0; i<children.size();++i)
        {
          copy(psiv_big.data(),psiv_big.data()+m,psiv.data());
          children[i]->ratioByRowSubstitution(psiv, psiv_big, irow, m, ratios);
        }
      }

      void ratioByRowSubstitution(ValueVector& psiv, const ValueVector& psic, int irow, int m, ValueVector& ratios)
      {
        psiv(from)=psic(to);
        ratios[my_id]=BLAS::dot(m,inverse[irow],psiv.data());
        if(children.empty()) return;
        ValueVector psiv_temp(m);
        for(int i=0; i<children.size();++i)
        {
          copy(psiv.begin(),psiv.end(),psiv_temp.data());
          children[i]->ratioByRowSubstitution(psiv_temp, psic, irow, m, ratios);
        }
      }
  */

  //this is just debug
  void write_node(std::ostream& os)
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
      for(i=0; i<peers.size()-1; i++)
        os << peers[i] << " ";
      os << peers[i] << "\"";
    }
    if(children.size())
      os << ">" << std::endl;
    else
      os << "/>" << std::endl;
    for(int i=0; i<children.size(); ++i)
      children[i]->write_node(os);
    if(children.size())
      os << "</node>"<<std::endl;
  }

private:

  inline ci_node(int m, int n, int id, int parent, int r, int a, configuration c):my_id(id),parent_id(parent),from(r),to(a),myC(c)
  {
    inverse.resize(m,m);
    peers.reserve(m);
    dpsiM.resize(m,n);
    d2psiM.resize(m,n);
    dpsiV.resize(n);
    d2psiV.resize(n);
  }

  void look_for_peers(std::vector<configuration>& confg, int& count, std::vector<int>& map2det)
  {
    int rem,add;
    int m=inverse.cols();
    int n=dpsiM.cols();
    for(int i=0; i<confg.size(); i++)
    {
      if(myC.isSingle(confg[i],rem,add) && !(confg[i].taken) )
      {
        confg[i].taken=true;
        map2det[i]=count;
        children.push_back(new ci_node(m,n,count++,my_id,rem,add,confg[i]));
        peers.push_back(rem);
        children.back()->look_for_peers(confg,count,map2det);
        children.back()->set_peers(peers);
      }
    }
  }

  void set_peers(std::vector<int>& p_peers)
  {
    if(children.size())
    {
      peers.clear();
      for(int i=0; i<children.size(); ++i)
        peers.push_back(children[i]->from);
      for(int i=0; i<children.size(); ++i)
        children[i]->set_peers(peers);
      for(int i=0; i<peers.size(); ++i)
      {
        bool found=false;
        for(int j=0; j<p_peers.size(); ++j)
          if(p_peers[j]==peers[i])
            found=true;
        if(!found)
          p_peers.push_back(peers[i]);
      }
    }
  }


};
}
#endif
