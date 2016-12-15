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


#include "Utilities/FakeRandom.h"

FakeRandom::FakeRandom() : m_val(0.5) {}

void FakeRandom::set_value(double val)
{
  m_val = val;
}

double FakeRandom::operator()()
{
  return m_val;
}

