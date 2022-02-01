//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RANDOMSEQUENCEGENERATOR_H
#define QMCPLUSPLUS_RANDOMSEQUENCEGENERATOR_H
#include <algorithm>
#include <type_traits>
#include "OhmmsPETE/OhmmsMatrix.h"
#include "ParticleBase/ParticleAttrib.h"
#include "Particle/MCCoords.hpp"
#include "config/stdlib/Constants.h"

/*!\fn template<class T> void assignGaussRand(T* restrict a, unsigned n)
  *\param a the starting pointer
  *\param n the number of type T to be assigned
  *\brief Assign Gaussian distributed random numbers using Box-Mueller algorithm. Called by overloaded funtions makeGaussRandom
  */
namespace qmcplusplus
{
template<class T, class RG>
inline void assignGaussRand(T* restrict a, unsigned n, RG& rng)
{
  OHMMS_PRECISION_FULL slightly_less_than_one = 1.0 - std::numeric_limits<OHMMS_PRECISION_FULL>::epsilon();
  int nm1                                     = n - 1;
  OHMMS_PRECISION_FULL temp1, temp2;
  for (int i = 0; i < nm1; i += 2)
  {
    temp1    = std::sqrt(-2.0 * std::log(1.0 - slightly_less_than_one * rng()));
    temp2    = 2.0 * M_PI * rng();
    a[i]     = temp1 * std::cos(temp2);
    a[i + 1] = temp1 * std::sin(temp2);
  }
  if (n % 2 == 1)
  {
    temp1  = std::sqrt(-2.0 * std::log(1.0 - slightly_less_than_one * rng()));
    temp2  = 2.0 * M_PI * rng();
    a[nm1] = temp1 * std::cos(temp2);
  }
}

/*!\fn template<class T> void assignUniformRand(T* restrict a, unsigned n)
  *\param a the starting pointer
  *\param n the number of type T to be assigned
  *\brief Assign unifor distributed random numbers [0,1)
  */
//template<class T>
//inline void assignUniformRand(T* restrict a, unsigned n) {
//  for(int i=0; i<n;i++) a[i] = Random();
//  //generate(a,a+n,Random);
//}
template<class T, class RG>
inline void assignUniformRand(T* restrict a, unsigned n, RG& rng)
{
  for (int i = 0; i < n; i++)
    a[i] = rng();
}

template<typename T, unsigned D, class RG>
inline void makeGaussRandomWithEngine(ParticleAttrib<TinyVector<T, D>>& a, RG& rng)
{
  assignGaussRand(&(a[0][0]), a.size() * D, rng);
}

template<typename T, unsigned D, class RG>
inline void makeGaussRandomWithEngine(std::vector<TinyVector<T, D>>& a, RG& rng)
{
  assignGaussRand(&(a[0][0]), a.size() * D, rng);
}

template<CoordsType CT, class RG>
inline void makeGaussRandomWithEngine(MCCoords<CT>& a, RG& rng)
{
  makeGaussRandomWithEngine(a.positions, rng);
  if constexpr (std::is_same<MCCoords<CT>, MCCoords<CoordsType::POS_SPIN>>::value)
    makeGaussRandomWithEngine(a.spins, rng);
}

template<typename T, class RG>
inline void makeGaussRandomWithEngine(std::vector<T>& a, RG& rng)
{
  static_assert(std::is_floating_point<T>::value,
                "makeGaussRandomWithEngine(std::vector<T>...) only implemented for floating point T");
  assignGaussRand(&(a[0]), a.size(), rng);
}

template<typename T, class RG>
inline void makeGaussRandomWithEngine(ParticleAttrib<T>& a, RG& rng)
{
  assignGaussRand(&(a[0]), a.size(), rng);
}

} // namespace qmcplusplus
#endif
