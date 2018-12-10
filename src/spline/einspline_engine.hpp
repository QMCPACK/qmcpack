//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
