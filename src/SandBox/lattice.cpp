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
    
    



#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Utilities/RandomGenerator.h>

namespace qmcplusplus
{
template<typename T, unsigned D>
inline Tensor<T,D> toCart(const Tensor<T,D>& g, const Tensor<T,D>& hess)
{
  Tensor<T,D> res;
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
    {
      T s=0;
      for(int k=0; k<D; ++k)
        for(int l=0; l<D; ++l)
          s+=g(k,i)*hess(k,l)*g(l,j);
      res(i,j)=s;
    }
  return res;
}
template<typename T, unsigned D>
inline Tensor<T,D> toCartT(const Tensor<T,D>& g, const Tensor<T,D>& hess)
{
  Tensor<T,D> res;
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
    {
      T s=0;
      for(int k=0; k<D; ++k)
        for(int l=0; l<D; ++l)
          s+=g(i,k)*hess(k,l)*g(j,l);
      res(i,j)=s;
    }
  return res;
}

template<typename T, unsigned D>
inline Tensor<T,D> outerProductSymm(const TinyVector<T,D>& v, const TinyVector<T,D>& w)
{
  Tensor<T,D> res;
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
      res(i,j)=v(i)*w(j)+v(j)*w(i);
  return res;
}
}

using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->rank());
  Random.init(0,1,11);
  ParticleSet ions;
  typedef ParticleSet::ParticleLayout_t LatticeType;
  typedef ParticleSet::TensorType       TensorType;
  typedef ParticleSet::RealType         RealType;
  RealType scale=5.432;
  TensorType fcc(0.1,0.5,0.5,0.5,-0.1,0.5,0.5,0.5,0.0);
  TensorType lat=scale*fcc;
  //set the BC to be periodic
  ions.Lattice.BoxBConds=1;
  //assign lattice
  ions.Lattice.set(lat);
  //print lattice
  ions.Lattice.print(std::cout);
  TensorType hess;
  TensorType gT=transpose(ions.Lattice.G);
  for(int i=0; i<hess.size(); ++i)
    hess(i)=Random();
  TensorType tmphs = dot(transpose(ions.Lattice.G),hess);
  TensorType hs = dot(tmphs,ions.Lattice.G);
  TensorType newhs=toCart(ions.Lattice.G,hess);
  TensorType newhsT=toCartT(gT,hess);
  std::cout << hs << std::endl<< std::endl;
  std::cout << newhs << std::endl<< std::endl;
  std::cout << newhsT-hs << std::endl;
  OHMMS::Controller->finalize();
  return 0;
}
