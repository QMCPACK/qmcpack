//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Jeongnim Kim and Kenneth P Esler
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
/**@file einspline_engine.hpp
 */
#ifndef QMCPLUSPLUS_EINSPLINE_ENGINE_HPP
#define QMCPLUSPLUS_EINSPLINE_ENGINE_HPP
#include <spline/bspline_traits.hpp>
#include <spline/einspline_impl.hpp>
namespace qmcplusplus
{

  /** einspline_engine
   *
   * The copy constructor is disabled.
   */
  template<typename ENGT>
    class einspline_engine
    {
      public:
        enum {D=bspline_engine_traits<ENGT>::DIM};
        typedef typename bspline_engine_traits<ENGT>::real_type real_type;
        typedef typename bspline_engine_traits<ENGT>::value_type value_type;
        typedef value_type* pointer;
        ///number of splines
        int num_splines;
        ///spline engine
        ENGT* spliner; 

        einspline_engine(ENGT* s=0):num_splines(0),spliner(s)
      { 
        if(spliner) num_splines=spliner->num_splines;
      }

        void create(TinyVector<real_type,D>& start, TinyVector<real_type,D>& end
            , TinyVector<int,D>& ng, bc_code bc, int nstates)
        {
          num_splines=nstates;
          spliner=einspline::create(spliner,start,end,ng,bc,nstates);
        }

        inline void set(int i, Array<value_type,D>& data)
        {
          einspline::set(spliner,i,data.data());
        }

        template<typename PT, typename VT> 
          inline void evaluate(const PT& r, VT& psi)
          { 
            einspline::evaluate(spliner,r,psi); 
          }

        template<typename PT, typename VT, typename GT> 
          inline void evaluate(const PT& r, VT& psi, GT& grad)
          { 
            einspline::evaluate(spliner,r,psi,grad); 
          }

        template<typename PT, typename VT, typename GT>
          inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, VT& lap)
          { 
            einspline::evaluate_vgl(spliner,r,psi,grad,lap); 
          }

        template<typename PT, typename VT, typename GT, typename HT>
          inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess)
          { 
            einspline::evaluate_vgh(spliner,r,psi,grad,hess); 
          }
      private:
        einspline_engine(const einspline_engine<ENGT>& rhs) {}
    };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/
