//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RANDOMSEQUENCEGENERATORGLOBAL_H
#define QMCPLUSPLUS_RANDOMSEQUENCEGENERATORGLOBAL_H

#include <algorithm>
#include <RandomGenerator.h>
#include "RandomSeqGenerator.h"

namespace qmcplusplus
{
template<typename T, unsigned D>
inline void makeGaussRandom(std::vector<TinyVector<T, D>>& a)
{
  assignGaussRand(&(a[0][0]), a.size() * D, Random);
}

///specialized functions: stick to overloading
template<typename T, unsigned D>
inline void makeGaussRandom(Matrix<TinyVector<T, D>>& a)
{
  assignGaussRand(&(a(0, 0)[0]), a.size() * D, Random);
}

template<typename T, unsigned D>
inline void makeGaussRandom(ParticleAttrib<TinyVector<T, D>>& a)
{
  assignGaussRand(&(a[0][0]), a.size() * D, Random);
}

template<typename T>
inline void makeGaussRandom(ParticleAttrib<T>& a)
{
  assignGaussRand(&(a[0]), a.size(), Random);
}

template<typename T, unsigned D>
inline void makeUniformRandom(ParticleAttrib<TinyVector<T, D>>& a)
{
  assignUniformRand(&(a[0][0]), a.size() * D, Random);
}

template<typename T>
inline void makeUniformRandom(ParticleAttrib<T>& a)
{
  assignUniformRand(&(a[0]), a.size(), Random);
}

template<typename T>
inline void makeSphereRandom(ParticleAttrib<TinyVector<T, 3>>& a)
{
  for (int i = 0; i < a.size(); i++)
  {
    bool failed = true;
    while (failed)
    {
      T x   = 1.0 - 2.0 * Random();
      T y   = 1.0 - 2.0 * Random();
      T z   = 1.0 - 2.0 * Random();
      T sep = std::sqrt(x * x + y * y + z * z);
      if (sep < 1)
      {
        T rinv  = 1.0 / sep;
        a[i][0] = x * rinv;
        a[i][1] = y * rinv;
        a[i][2] = z * rinv;
        failed  = false;
      }
    }
  }
}

template<typename T>
inline void makeSphereRandom(ParticleAttrib<TinyVector<T, 2>>& a)
{
  for (int i = 0; i < a.size(); i++)
  {
    bool failed = true;
    while (failed)
    {
      T x   = 1.0 - 2.0 * Random();
      T y   = 1.0 - 2.0 * Random();
      T sep = std::sqrt(x * x + y * y);
      if (sep < 1)
      {
        T rinv  = 1.0 / sep;
        a[i][0] = x * rinv;
        a[i][1] = y * rinv;
        failed  = false;
      }
    }
  }
}

} // namespace qmcplusplus


#endif
