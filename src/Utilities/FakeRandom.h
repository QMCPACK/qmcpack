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



#ifndef FAKERANDOM_H
#define FAKERANDOM_H

// Deterministic generator for testing.

class FakeRandom
{
public:
  FakeRandom();
  typedef unsigned int uint_type;
  double operator()();
  void set_value(double val);
private:
  double m_val;
};

#endif
