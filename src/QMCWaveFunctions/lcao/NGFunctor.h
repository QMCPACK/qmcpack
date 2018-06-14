//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_NGFUNCTOR_FOR_SOA_H
#define QMCPLUSPLUS_NGFUNCTOR_FOR_SOA_H

#include "Numerics/OneDimQuinticSpline.h"
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
  //this is a temporary solution before switching to BsplineFunctor<T>
  template<typename T>
  struct NGFunctor//: public OptimizableFunctorBase
  {
    typedef T                    real_type;
    typedef T                    value_type;
    typedef T                    point_type;
    typedef OneDimGridBase<real_type>    grid_type;
    typedef OneDimQuinticSpline<real_type> functor_type;
    functor_type myFunc;
    real_type Y, dY, d2Y, d3Y;

    NGFunctor(grid_type* agrid):myFunc(agrid) { }

    template<typename VV>
      NGFunctor(grid_type* agrid, const VV& nv):myFunc(agrid,nv) { }

    NGFunctor(const NGFunctor& in)=default;

    NGFunctor(const NGFunctor& in, grid_type* agrid, bool grid_manager): myFunc(in.myFunc)
    {
      myFunc.m_grid=agrid;
      myFunc.setGridManager(grid_manager);
    }

    void checkInVariables(opt_variables_type& active) {}
    void checkOutVariables(const opt_variables_type& active) {}
    void resetParameters(const opt_variables_type& active) {}
    void reset() {}
    inline real_type f(real_type r)
    {
      return myFunc.f(r);
    }
    inline real_type df(real_type r)
    {
      return myFunc.df(r);
    }
    bool put(xmlNodePtr cur)
    {
      return true;
    }

    //OptimizableFunctorBase* makeClone()
    //NGFunctor<T>* makeClone()
    //{
    //  return new NGFunctor<T>(*this);
    //  //NGFunctor<T> *myclone=new NGFunctor<T>(*this);
    //  //myclone->myFunc.m_grid=myFunc.m_grid->makeClone();
    //  //myclone->setGridManager(true);
    //  //return myclone;
    //}

    inline value_type operator()(int i) const
    {
      return myFunc(i);
    }
    inline value_type& operator()(int i)
    {
      return myFunc(i);
    }

    inline real_type evaluate(real_type r)
    {
      return myFunc.splint(r);
    }

    inline real_type  evaluate(real_type r, real_type& du, real_type& d2u)
    {
      return myFunc.splint(r,du,d2u);
    }

    inline real_type
      evaluate(real_type r, real_type& du, real_type& d2u, real_type& d3u)
      {
        return myFunc.splint(r,du,d2u,d3u);
      }
    inline grid_type& grid()
    {
      return myFunc.grid();
    }
    inline void setGridManager(bool willmanage)
    {
      myFunc.setGridManager(willmanage);
    }

    inline void spline(int imin, value_type yp1, int imax, value_type ypn)
    {
      myFunc.spline(imin,yp1,imax,ypn);
    }

    inline void resize(int n)
    {
      myFunc.resize(n);
    }
  };
}
#endif
