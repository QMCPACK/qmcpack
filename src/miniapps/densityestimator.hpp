//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Qin, qinkenjaist@icloud.com, JAIST.
//                
//
// File created by: Ken Qin, qinkenjaist@icloud.com, JAIST.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file densityestimator.hpp
 */

#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <spline2/MultiBspline.hpp>
#include <simd/allocator.hpp>
#include <iostream>

namespace qmcplusplus
{
  template<typename T>
    struct density_estimator
    {
      ///number of electrons
      int neles;
      ///bumber of walkers
      int nwalkers;
      ///position of the walkers
      int lastBlock;
      ///number of the grad
      int ngrad;
      
      aligned_vector<gContainer_type*> estimators;
      aligned_vector<vContainer_type*> psi;
      aligned_vector<gContainer_type*> grad;
      aligned_vector<hContainer_type*> hess;
            
      ///destructors
      ~density_estimator()
      {
        if(psi.size()) clean();
        if(Owner)
          for(int i=0; i<nBlocks; ++i) delete einsplines[i];
      }

      void clean()
      {
        const int n=psi.size();
        for(int i=0; i<n; ++i)
        {
          delete psi[i];
          delete grad[i];
          delete hess[i];
        }
      }

      void calculate()
      {
        ///calculate the density
        
      }
  };
}

#endif
