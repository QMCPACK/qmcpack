//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "FakeRandom.h"

namespace qmcplusplus
{
template<class T>
FakeRandom<T>::FakeRandom() = default;

template<class T>
void FakeRandom<T>::set_value(double val)
{
  m_val = val;
}

template<class T>
T FakeRandom<T>::operator()()
{
  return m_val;
}

template class FakeRandom<float>;
template class FakeRandom<double>;

} // namespace qmcplusplus