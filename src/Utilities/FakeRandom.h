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

#include "RandomBase.h"

// Deterministic generator for testing.

namespace qmcplusplus
{
template<typename T = double>
class FakeRandom : public RandomBase<T>
{
public:
  using result_type = typename RandomBase<T>::result_type;
  using uint_type   = typename RandomBase<T>::uint_type;

  FakeRandom();
  T operator()() override;

  void init(int iseed) override{};
  void seed(uint_fast32_t aseed) override{};
  void write(std::ostream& rout) const override { rout << m_val; };
  void read(std::istream& rin) override { rin >> m_val; };
  void load(const std::vector<uint_type>& newstate) override{};
  void save(std::vector<uint_type>& curstate) const override{};
  size_t state_size() const override { return 0; };
  std::unique_ptr<RandomBase<T>> makeClone() const override { return std::make_unique<FakeRandom<T>>(*this); }

  void set_value(double val);

private:
  T m_val{0.5};
};

} // namespace qmcplusplus
#endif
