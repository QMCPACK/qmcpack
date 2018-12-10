//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"

//using namespace qmcplusplus;
//RandomGenerator_t
//Random(CommCreate::get()->getNodeID(), CommCreate::get()->getNumNodes());

qmcplusplus::RandomGenerator_t qmcplusplus::Random;

// /**class GaussinRandomGenerator
//  *\brief A filter class that converts random numbers [0,1) -> gaussian
//  */
// class GaussianRandomGenerator {
// public:

//   typedef RandomGenerator_t::Return_t Return_t;

//   GaussianRandomGenerator(RandomGenerator_t& rg):d_engine(rg) { }

//   inline Return_t operator()(){
//     if(newpair) {
//       d_engine.bivariate(gauss0,gauss1);
//       newpair = false;
//       return gauss0;
//     } else {
//       newpair = true;
//       return gauss1;
//     }
//   }
// private:
//   RandomGenerator_t d_engine;
//   bool newpair;
//   Return_t gauss0, gauss1;
// };
// GaussianRandomGenerator GaussianRandom(Random);


//   class GaussianRandom {
//   public:
//     typedef RandomGenerator_t::Return_t Return_t;
//     GaussianRandom(RandomGenerator_t& rg, Return_t sig=1.0, Return_t c0=0.0):
//       d_engine(rg), newpair(true){ Sigma2 = sig*sig; Center = c0;}
//     inline Return_t operator()(){
//       if(newpair) {
// 	d_engine.bivariate(gauss0,gauss1);
// 	newpair = false;
// 	return Sigma2*gauss0+Center;
//       } else {
// 	newpair = true;
// 	return Sigma2*gauss1+Center;
//       }
//     }
//   private:
//     RandomGenerator_t& d_engine;
//     bool newpair;
//     Return_t gauss0, gauss1, Sigma2, Center;

//   };

